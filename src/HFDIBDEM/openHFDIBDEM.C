/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /---| |/ / | |___| |  | |
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/    |___/  |_____|_|  |_|
      | |                     H ybrid F ictitious D omain - I mmersed B oundary
      |_|                                        and D iscrete E lement M ethod
-------------------------------------------------------------------------------
License

    openHFDIB-DEM is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIB. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    Foam

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*)
\*---------------------------------------------------------------------------*/
#include "openHFDIBDEM.H"
#include "polyMesh.H"
#include "fvCFD.H"
#include "fvMatrices.H"
#include "geometricOneField.H"

#include "interpolationCellPoint.H"
#include "interpolationCell.H"

#include "scalarMatrices.H"
#include "OFstream.H"
#include <iostream>
#include "defineExternVars.H"
#include "parameters.H"
#include "finiteWallGeometry.H"

#define ORDER 2

#include <fstream>
#include <iomanip>

using namespace Foam;
using namespace contactModel;

namespace
{

label dominantWallDirection(const vector& normal)
{
    const scalar ax = mag(normal[0]);
    const scalar ay = mag(normal[1]);
    const scalar az = mag(normal[2]);

    if (ax >= ay && ax >= az)
    {
        return 0;
    }
    else if (ay >= ax && ay >= az)
    {
        return 1;
    }

    return 2;
}

label axisAlignedWallDirection(const vector& normal)
{
    if (mag(normal) <= VSMALL)
    {
        return -1;
    }

    const vector normalHat(normal/mag(normal));
    label normalDirection = -1;
    label nonZeroDirections = 0;

    for (label dir = 0; dir < 3; dir++)
    {
        if (mag(normalHat[dir]) > 1e-8)
        {
            normalDirection = dir;
            nonZeroDirections++;
        }
    }

    return nonZeroDirections == 1 ? normalDirection : -1;
}

}

//---------------------------------------------------------------------------//
openHFDIBDEM::openHFDIBDEM(const Foam::fvMesh& mesh)
:
mesh_(mesh),
HFDIBDEMDict_
(
    IOobject
    (
        "HFDIBDEMDict",
        "constant",
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
transportProperties_
(
    IOobject
    (
        "transportProperties",
        "constant",
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
bodyNames_(HFDIBDEMDict_.lookup("bodyNames")),
prtcInfoTable_(0),
stepDEM_(readScalar(HFDIBDEMDict_.lookup("stepDEM"))),
recordSimulation_(readBool(HFDIBDEMDict_.lookup("recordSimulation")))
{
    wallMotionInfo::clear();

    materialProperties::matProps_insert(
        "None",
        materialInfo("None", 1, 1, 1, 1, 1)
    );

    if(HFDIBDEMDict_.found("recordFirstTimeStep"))
    {
        recordFirstTimeStep_ = readBool(HFDIBDEMDict_.lookup("recordFirstTimeStep"));
    }

    if(HFDIBDEMDict_.found("nSolidsInDomain"))
    {
        solverInfo::setNSolidsTreshnold(readLabel(HFDIBDEMDict_.lookup("nSolidsInDomain")));
    }

    if(HFDIBDEMDict_.found("weightDragByLambda"))
    {
        solverInfo::setDragOverLambdaCeil(readBool(HFDIBDEMDict_.lookup("weightDragByLambda")));
        if(solverInfo::getDragWeightByLambda())
        {
            Info << " -- Drag force is integrated over the immersedbody cells and weighted by solid volume fraction"  << endl;
        }
        else
        {
            Info << " -- Drag force is integrated over the immersedbody cells intergration is over the whole cells" << endl;
        }
    }
    else
    {
        Info << " -- Switch not found, assuming intergration over the immersedbody cells and weighted by solid volume fraction" << endl;
    }


    dictionary demDic = HFDIBDEMDict_.subDict("DEM");
    dictionary materialsDic = demDic.subDict("materials");
    List<word> materialsNames = materialsDic.toc();
    if(demDic.found("increasedDamping"))
    {
        contactModelInfo::setIncreasedDamping(readBool(demDic.lookup("increasedDamping")));
    }
    else
    {
        contactModelInfo::setIncreasedDamping(false);
    }

    forAll(materialsNames, matI)
    {
        dictionary matIDic = materialsDic.subDict(materialsNames[matI]);
        scalar eps = readScalar(matIDic.lookup("eps"));
        if(!contactModelInfo::getIncreasedDamping())
        {
            eps = 0.906463027*eps + 0.093538298;//LinearRegression based on LIGGGHTS data testing
            eps = min(eps, 1.0);
        }

        materialProperties::matProps_insert(
            materialsNames[matI],
            materialInfo(
                materialsNames[matI],
                readScalar(matIDic.lookup("Y")),
                readScalar(matIDic.lookup("nu")),
                readScalar(matIDic.lookup("mu")),
                readScalar(matIDic.lookup("adhN")),
                // readScalar(matIDic.lookup("eps"))
                eps
            )
        );
    }

    if(demDic.found("interfaceAdh"))
    {
        dictionary interfAdhDic = demDic.subDict("interfaceAdh");
        List<word> interNames = interfAdhDic.toc();
        forAll(interNames, interI)
        {
            dictionary interDicI = interfAdhDic.subDict(interNames[interI]);
            wordList interMat = interDicI.lookup("materials");
            string interKey;
            if(interMat[0] < interMat[1])
            {
                interKey += interMat[0];
                interKey += "-";
                interKey += interMat[1];
            }
            else
            {
                interKey += interMat[1];
                interKey += "-";
                interKey += interMat[0];
            }

            interAdhesion::interAdhesion_insert(
                interKey,
                readScalar(interDicI.lookup("value"))
            );
        }
    }

    if(demDic.found("LcCoeff"))
    {
        contactModelInfo::setLcCoeff(readScalar(demDic.lookup("LcCoeff")));
    }
    else
    {
        contactModelInfo::setLcCoeff(4.0);
    }

    if(demDic.found("rotationModel"))
    {
        word rotModel = demDic.lookup("rotationModel");
        if(rotModel == "chen2012")
        {
            contactModelInfo::setRotationModel(0);
        }
        else if(rotModel == "mindlin1953")
        {
            contactModelInfo::setRotationModel(1);
        }
        else
        {
            Info << "Rotation Model not recognized, setting to default mindlin1953" << endl;
            contactModelInfo::setRotationModel(1);
        }    
    }
    else
    {
        Info << "Rotation Model not recognized, setting to default mindlin1953" << endl;
        contactModelInfo::setRotationModel(1);
    }
    

    Info <<" -- Coefficient for characteristic Lenght Lc is set to : "<< contactModelInfo::getLcCoeff() << endl;

    dictionary patchDic = demDic.subDict("collisionPatches");
    List<word> patchNames = patchDic.toc();

    forAll(patchNames, patchI)
    {
        dictionary patchIDic = patchDic.subDict(patchNames[patchI]);

        word patchMaterial = patchIDic.lookup("material");
        vector patchNVec = patchIDic.lookup("nVec");
        vector planePoint = patchIDic.lookup("planePoint");

        bool finitePatch = false;
        vector minBound = vector::zero;
        vector maxBound = vector::zero;

        if (patchIDic.found("finitePatch"))
        {
           finitePatch = readBool(patchIDic.lookup("finitePatch"));
        }

        if (finitePatch && patchIDic.found("vertices"))
        {
            pointField patchVertices(patchIDic.lookup("vertices"));

            if (patchVertices.size() < 3)
            {
                FatalIOErrorInFunction(HFDIBDEMDict_)
                    << "finite collisionPatch " << patchNames[patchI]
                    << " needs at least three ordered vertices."
                    << exit(FatalIOError);
            }

            wallPlaneInfo::wallPlaneInfo_insert
            (
                patchNames[patchI],
                patchNVec,
                planePoint,
                patchVertices
            );

            Info<< " -- collisionPatch " << patchNames[patchI]
                << " is treated as a finite polygon with "
                << patchVertices.size() << " vertices"
                << endl;
        }
        else if
        (
            finitePatch
        && patchIDic.found("minBound")
        && patchIDic.found("maxBound")
        )
        {
            minBound = patchIDic.lookup("minBound");
            maxBound = patchIDic.lookup("maxBound");

            wallPlaneInfo::wallPlaneInfo_insert
            (
                patchNames[patchI],
                patchNVec,
                planePoint,
                true,
                minBound,
                maxBound
            );

            Info<< " -- collisionPatch " << patchNames[patchI]
                << " is treated as finite patch"
                << " minBound " << minBound
                << " maxBound " << maxBound
                << endl;
       }
        else
       {
           if
            (
               finitePatch
            && !patchIDic.found("vertices")
            && (!patchIDic.found("minBound") || !patchIDic.found("maxBound"))
            )
            {
                Info<< " -- finitePatch requested for collisionPatch "
                    << patchNames[patchI]
                    << " but neither vertices nor minBound/maxBound "
                    << "is available. "
                    << "Treating it as an infinite wall."
                    << endl;
            }

            wallPlaneInfo::wallPlaneInfo_insert
            (
                patchNames[patchI],
                patchNVec,
                planePoint
            );
        }

        if
        (
            wallPlaneInfo::getWallFiniteInfo().found(patchNames[patchI])
         && wallPlaneInfo::getWallFiniteInfo()[patchNames[patchI]]
        )
        {
            // Construct once during dictionary parsing so malformed,
            // non-coplanar, or non-convex patches fail before body creation.
            (void)finiteWallGeometry::lookup(patchNames[patchI]);

            if
            (
                wallPlaneInfo::usesLegacyAxisAlignedFinitePath
                (
                    patchNames[patchI]
                )
            )
            {
                Info<< " -- finite-wall contact path: Git-main "
                    << "axis-aligned rectangular fast path" << endl;
            }
            else
            {
                Info<< " -- finite-wall contact path: general polygon path"
                    << endl;
            }
        }

        wallMatInfo::wallMatInfo_insert
        (
            patchNames[patchI],
            materialProperties::getMatProps()[patchMaterial]
        );

        if (patchIDic.isDict("motion"))
        {
            const dictionary motionDic(patchIDic.subDict("motion"));
            const word motionType(motionDic.lookup("type"));

            if (motionType == "meshPatchTranslation")
            {
                const word sourcePatch(motionDic.lookup("patch"));

                wallMotionInfo::insert
                (
                    patchNames[patchI],
                    sourcePatch,
                    planePoint,
                    mesh_.time().timeIndex()
                );

                Info<< " -- collisionPatch " << patchNames[patchI]
                    << " follows translating mesh patch " << sourcePatch
                    << endl;
            }
            else
            {
                FatalIOErrorInFunction(HFDIBDEMDict_)
                    << "Unsupported motion type '" << motionType
                    << "' for collisionPatch " << patchNames[patchI] << nl
                    << "Only meshPatchTranslation is implemented."
                    << exit(FatalIOError);
            }
        }

        if (patchIDic.isDict("supportLimit"))
        {
            const dictionary supportDic(patchIDic.subDict("supportLimit"));
            const word supportType(supportDic.lookup("type"));
            const word leaderWall(supportDic.lookup("wall"));
            const word followedBound(supportDic.lookup("bound"));

            if (supportType != "followWall")
            {
                FatalIOErrorInFunction(HFDIBDEMDict_)
                    << "Unsupported supportLimit type '" << supportType
                    << "' for collisionPatch " << patchNames[patchI] << nl
                    << "Only followWall is implemented."
                    << exit(FatalIOError);
            }

            if (followedBound != "min" && followedBound != "max")
            {
                FatalIOErrorInFunction(HFDIBDEMDict_)
                    << "supportLimit bound must be min or max for "
                    << patchNames[patchI]
                    << exit(FatalIOError);
            }

            wallSupportLeader_.insert(patchNames[patchI], leaderWall);
            wallSupportBound_.insert(patchNames[patchI], followedBound);
        }
    }

    const List<string> supportWalls(wallSupportLeader_.sortedToc());

    forAll(supportWalls, wallI)
    {
        const string& followerWall = supportWalls[wallI];
        const string& leaderWall = wallSupportLeader_[followerWall];

        if (!wallMotionInfo::isMoving(leaderWall))
        {
            FatalIOErrorInFunction(HFDIBDEMDict_)
                << "collisionPatch " << followerWall
                << " follows wall " << leaderWall
                << ", but that wall has no meshPatchTranslation motion."
                << exit(FatalIOError);
        }

        if
        (
            !wallPlaneInfo::getWallFiniteInfo().found(followerWall)
         || !wallPlaneInfo::getWallFiniteInfo()[followerWall]
        )
        {
            FatalIOErrorInFunction(HFDIBDEMDict_)
                << "supportLimit requires finitePatch true for "
                << followerWall
                << exit(FatalIOError);
        }

        const vector& leaderNormal =
            wallPlaneInfo::getWallPlaneInfo()[leaderWall][0];

        const label supportDirection =
            axisAlignedWallDirection(leaderNormal);

        if (supportDirection < 0)
        {
            FatalIOErrorInFunction(HFDIBDEMDict_)
                << "supportLimit leader wall " << leaderWall
                << " must translate along one Cartesian axis. nVec is "
                << leaderNormal << "."
                << exit(FatalIOError);
        }

        if
        (
            wallPlaneInfo::getWallPatchVerticesInfo().found(followerWall)
        )
        {
            const vector followerNormal =
                wallPlaneInfo::getWallPlaneInfo()[followerWall][0]
               /mag(wallPlaneInfo::getWallPlaneInfo()[followerWall][0]);

            const vector leaderNormalHat = leaderNormal/mag(leaderNormal);

            if (mag(followerNormal & leaderNormalHat) > 1e-8)
            {
                FatalIOErrorInFunction(HFDIBDEMDict_)
                    << "polygonal supportLimit follower " << followerWall
                    << " surface must be parallel to the leader motion "
                    << "direction. Follower nVec " << followerNormal
                    << " is not perpendicular to leader nVec "
                    << leaderNormalHat << "."
                    << exit(FatalIOError);
            }

            const pointField& followerVertices =
                wallPlaneInfo::getWallPatchVerticesInfo()[followerWall];

            const vector& minBound =
                wallPlaneInfo::getWallMinBoundInfo()[followerWall];

            const vector& maxBound =
                wallPlaneInfo::getWallMaxBoundInfo()[followerWall];

            const scalar supportSpan =
                maxBound[supportDirection] - minBound[supportDirection];

            const scalar supportTolerance = max
            (
                scalar(1e-12),
                scalar(1e-10)*max(mag(maxBound - minBound), scalar(1))
            );

            if (supportSpan <= supportTolerance)
            {
                FatalIOErrorInFunction(HFDIBDEMDict_)
                    << "polygonal supportLimit follower " << followerWall
                    << " has no extent along leader direction "
                    << supportDirection << "."
                    << exit(FatalIOError);
            }

            label minCapVertices = 0;
            label maxCapVertices = 0;

            forAll(followerVertices, vertexI)
            {
                const scalar coordinate =
                    followerVertices[vertexI][supportDirection];

                if
                (
                    mag(coordinate - minBound[supportDirection])
                 <= supportTolerance
                )
                {
                    minCapVertices++;
                }
                else if
                (
                    mag(coordinate - maxBound[supportDirection])
                 <= supportTolerance
                )
                {
                    maxCapVertices++;
                }
                else
                {
                    FatalIOErrorInFunction(HFDIBDEMDict_)
                        << "polygonal supportLimit follower "
                        << followerWall << " must be a two-cap extrusion "
                        << "along direction " << supportDirection << ". "
                        << "Vertex " << followerVertices[vertexI]
                        << " lies on neither cap."
                        << exit(FatalIOError);
                }
            }

            if (minCapVertices < 2 || maxCapVertices < 2)
            {
                FatalIOErrorInFunction(HFDIBDEMDict_)
                    << "polygonal supportLimit follower " << followerWall
                    << " needs at least two vertices on each support cap. "
                    << "Found " << minCapVertices << " and "
                    << maxCapVertices << "."
                    << exit(FatalIOError);
            }

            Info<< " -- polygonal collisionPatch " << followerWall
                << " safely follows " << leaderWall
                << " on its " << wallSupportBound_[followerWall]
                << " support cap" << endl;
        }
    }

    resetMovingWallsFromMesh();

    if(demDic.found("cyclicPatches"))
    {
        Info << "CyclicPatches Found " << endl;
        dictionary cyclicPatchDic = demDic.subDict("cyclicPatches");
        List<word> cyclicPatchNames = cyclicPatchDic.toc();
        forAll(cyclicPatchNames, patchI)
        {
            vector patchNVec = cyclicPatchDic.subDict(cyclicPatchNames[patchI]).lookup("nVec");
            vector planePoint = cyclicPatchDic.subDict(cyclicPatchNames[patchI]).lookup("planePoint");
            word neighbourPatch = cyclicPatchDic.subDict(cyclicPatchNames[patchI]).lookup("neighbourPatch");

            cyclicPlaneInfo::insert(
                cyclicPatchNames[patchI],
                patchNVec,
                planePoint,
                neighbourPatch
            );
        }
        Info << "CyclicPatches  " <<  cyclicPatchNames <<endl;
    }

    if (HFDIBDEMDict_.found("geometricD"))
    {
        geometricD = HFDIBDEMDict_.lookup("geometricD");
    }
    else
    {
        geometricD = mesh_.geometricD();
    }

    forAll (geometricD, direction)
    {
        if (geometricD[direction] == -1)
        {
            case3D = false;
            emptyDir[direction] = 1;
            emptyDim = direction;
            break;
        }
    }

    if (HFDIBDEMDict_.isDict("virtualMesh"))
    {
        dictionary vMDic = HFDIBDEMDict_.subDict("virtualMesh");
        virtualMeshLevel::setVirtualMeshLevel(readScalar(vMDic.lookup("level")),readScalar(vMDic.lookup("charCellSize")));
        Info <<" -- VirtMesh Decomposition Level is set to        : "<< virtualMeshLevel::getVirtualMeshLevel() << endl;
        Info <<" -- VirtMesh charCellSize for boundary is set to  : "<< virtualMeshLevel::getCharCellSize() << endl;

    }
    else
    {
        virtualMeshLevel::setVirtualMeshLevel(1,1);
        Info <<" -- VirtMesh Decomposition Level is set to        : "<< virtualMeshLevel::getVirtualMeshLevel() << endl;
        Info <<" -- VirtMesh charCellSize for boundary is set to  : "<< virtualMeshLevel::getCharCellSize() << endl;

    }

    recordOutDir_ = mesh_.time().rootPath() + "/" + mesh_.time().globalCaseName() + "/bodiesInfo";
}
//---------------------------------------------------------------------------//
openHFDIBDEM::~openHFDIBDEM()
{}
//---------------------------------------------------------------------------//
vector openHFDIBDEM::movingWallPatchPlanePoint
(
    const string& wallName
) const
{
    const wallMotionState& motionState = wallMotionInfo::state(wallName);
    const word& sourcePatch = motionState.sourcePatch();
    const label patchId = mesh_.boundaryMesh().findPatchID(sourcePatch);

    if (patchId < 0)
    {
        FatalErrorInFunction
            << "Cannot find mesh boundary patch '" << sourcePatch
            << "' for moving DEM wall '" << wallName << "'."
            << exit(FatalError);
    }

    const vectorField& patchCf = mesh_.Cf().boundaryField()[patchId];
    const vectorField& patchSf = mesh_.Sf().boundaryField()[patchId];

    scalar patchArea = 0;
    vector areaWeightedCenter(vector::zero);
    vector areaVector(vector::zero);

    forAll(patchSf, faceI)
    {
        const scalar faceArea = mag(patchSf[faceI]);

        patchArea += faceArea;
        areaWeightedCenter += faceArea*patchCf[faceI];
        areaVector += patchSf[faceI];
    }

    reduce(patchArea, sumOp<scalar>());
    reduce(areaWeightedCenter, sumOp<vector>());
    reduce(areaVector, sumOp<vector>());

    if (patchArea <= VSMALL || mag(areaVector) <= VSMALL)
    {
        FatalErrorInFunction
            << "Moving DEM wall '" << wallName << "' is linked to patch '"
            << sourcePatch << "', but its global area is zero."
            << exit(FatalError);
    }

    const vector planePoint(areaWeightedCenter/patchArea);
    const vector patchNormal(areaVector/mag(areaVector));

    const vector configuredNormal
    (
        wallPlaneInfo::getWallPlaneInfo()[wallName][0]
    );

    if (mag(configuredNormal) <= VSMALL)
    {
        FatalErrorInFunction
            << "Moving DEM wall '" << wallName << "' has a zero nVec."
            << exit(FatalError);
    }

    const vector normalHat(configuredNormal/mag(configuredNormal));

    if ((normalHat & patchNormal) < 1.0 - 1e-6)
    {
        FatalErrorInFunction
            << "Normal mismatch for moving DEM wall '" << wallName << "'." << nl
            << "Configured nVec: " << normalHat << nl
            << "Mesh patch normal: " << patchNormal << nl
            << "Source patch: " << sourcePatch
            << exit(FatalError);
    }

    scalar maxPlaneDistance = 0;

    forAll(patchCf, faceI)
    {
        maxPlaneDistance = max
        (
            maxPlaneDistance,
            mag((patchCf[faceI] - planePoint) & normalHat)
        );
    }

    reduce(maxPlaneDistance, maxOp<scalar>());

    const scalar planeTolerance =
        max(scalar(1e-10), scalar(1e-8)*mag(mesh_.bounds().span()));

    if (maxPlaneDistance > planeTolerance)
    {
        FatalErrorInFunction
            << "Mesh patch '" << sourcePatch << "' for moving DEM wall '"
            << wallName << "' is not planar within tolerance "
            << planeTolerance << ". Maximum distance is "
            << maxPlaneDistance << "."
            << exit(FatalError);
    }

    return planePoint;
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::setMovingWallsAtFraction(const scalar fraction)
{
    const scalar alpha = min(max(fraction, scalar(0)), scalar(1));
    const List<string> movingWalls(wallMotionInfo::all().sortedToc());

    forAll(movingWalls, wallI)
    {
        const string& wallName = movingWalls[wallI];
        const wallMotionState& motionState = wallMotionInfo::state(wallName);

        const vector planePoint =
            (1.0 - alpha)*motionState.previousPlanePoint()
          + alpha*motionState.nextPlanePoint();

        wallPlaneInfo::setWallPlanePoint(wallName, planePoint);

        if
        (
            wallPlaneInfo::getWallFiniteInfo().found(wallName)
         && wallPlaneInfo::getWallFiniteInfo()[wallName]
         && !wallPlaneInfo::getWallPatchVerticesInfo().found(wallName)
        )
        {
            vector minBound =
                wallPlaneInfo::getWallMinBoundInfo()[wallName];

            vector maxBound =
                wallPlaneInfo::getWallMaxBoundInfo()[wallName];

            const label normalDirection = dominantWallDirection
            (
                wallPlaneInfo::getWallPlaneInfo()[wallName][0]
            );

            // The finite-wall algorithm uses only tangential bounds.  Keeping
            // the normal components on the moving plane makes the registry
            // geometrically self-consistent without altering that algorithm.
            minBound[normalDirection] = planePoint[normalDirection];
            maxBound[normalDirection] = planePoint[normalDirection];

            wallPlaneInfo::setWallFiniteBounds
            (
                wallName,
                minBound,
                maxBound
            );
        }
    }

    const List<string> followerWalls(wallSupportLeader_.sortedToc());

    forAll(followerWalls, wallI)
    {
        const string& followerWall = followerWalls[wallI];
        const string& leaderWall = wallSupportLeader_[followerWall];

        const vector& leaderPoint =
            wallPlaneInfo::getWallPlaneInfo()[leaderWall][1];

        const label leaderNormalDirection = dominantWallDirection
        (
            wallPlaneInfo::getWallPlaneInfo()[leaderWall][0]
        );

        vector minBound =
            wallPlaneInfo::getWallMinBoundInfo()[followerWall];

        vector maxBound =
            wallPlaneInfo::getWallMaxBoundInfo()[followerWall];

        const scalar supportTolerance = max
        (
            scalar(1e-12),
            scalar(1e-10)*max(mag(maxBound - minBound), scalar(1))
        );

        if (wallSupportBound_[followerWall] == "max")
        {
            if
            (
                leaderPoint[leaderNormalDirection]
             <= minBound[leaderNormalDirection] + supportTolerance
            )
            {
                FatalErrorInFunction
                    << "Moving support wall " << leaderWall
                    << " crossed the fixed cap of follower "
                    << followerWall << "."
                    << exit(FatalError);
            }

            maxBound[leaderNormalDirection] =
                leaderPoint[leaderNormalDirection];
        }
        else
        {
            if
            (
                leaderPoint[leaderNormalDirection]
             >= maxBound[leaderNormalDirection] - supportTolerance
            )
            {
                FatalErrorInFunction
                    << "Moving support wall " << leaderWall
                    << " crossed the fixed cap of follower "
                    << followerWall << "."
                    << exit(FatalError);
            }

            minBound[leaderNormalDirection] =
                leaderPoint[leaderNormalDirection];
        }

        wallPlaneInfo::setWallFiniteBounds
        (
            followerWall,
            minBound,
            maxBound
        );
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::resetMovingWallsFromMesh()
{
    const List<string> movingWalls(wallMotionInfo::all().sortedToc());

    forAll(movingWalls, wallI)
    {
        const string& wallName = movingWalls[wallI];
        const vector planePoint = movingWallPatchPlanePoint(wallName);

        wallMotionState& motionState = wallMotionInfo::state(wallName);
        motionState.setPreviousPlanePoint(planePoint);
        motionState.setNextPlanePoint(planePoint);
        motionState.setVelocity(vector::zero);
        motionState.setTimeIndex(mesh_.time().timeIndex());

        Info<< " -- initialized moving DEM wall " << wallName
            << " from patch " << motionState.sourcePatch()
            << " at " << planePoint << endl;
    }

    setMovingWallsAtFraction(1.0);
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::updateMovingWallsFromMesh()
{
    const List<string> movingWalls(wallMotionInfo::all().sortedToc());
    const label currentTimeIndex = mesh_.time().timeIndex();
    const scalar deltaTime = mesh_.time().deltaT().value();

    forAll(movingWalls, wallI)
    {
        const string& wallName = movingWalls[wallI];
        wallMotionState& motionState = wallMotionInfo::state(wallName);

        if (motionState.timeIndex() != currentTimeIndex)
        {
            motionState.setPreviousPlanePoint
            (
                motionState.nextPlanePoint()
            );
            motionState.setTimeIndex(currentTimeIndex);
        }

        const vector nextPlanePoint = movingWallPatchPlanePoint(wallName);

        const vector displacement =
            nextPlanePoint - motionState.previousPlanePoint();

        const vector normal =
            wallPlaneInfo::getWallPlaneInfo()[wallName][0]
           /mag(wallPlaneInfo::getWallPlaneInfo()[wallName][0]);

        const vector tangentialDisplacement =
            displacement - normal*(displacement & normal);

        const scalar motionTolerance =
            max(scalar(1e-10), scalar(1e-8)*mag(mesh_.bounds().span()));

        if (mag(tangentialDisplacement) > motionTolerance)
        {
            FatalErrorInFunction
                << "Moving DEM wall '" << wallName
                << "' has tangential displacement "
                << tangentialDisplacement << "." << nl
                << "Only translation along the fixed wall normal is supported."
                << exit(FatalError);
        }

        motionState.setNextPlanePoint(nextPlanePoint);

        if (deltaTime > VSMALL)
        {
            motionState.setVelocity
            (
                normal*(displacement & normal)/deltaTime
            );
        }
        else
        {
            motionState.setVelocity(vector::zero);
        }
    }

    // CFD and body reconstruction use the end-of-step mesh.  DEM contact
    // later interpolates back through this same interval.
    setMovingWallsAtFraction(1.0);
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::writeMovingWallDEMForces
(
    const List<string>& wallNames,
    const List<vector>& wallImpulses,
    const scalar deltaTime
) const
{
    if (!Pstream::master() || wallNames.size() == 0 || deltaTime <= VSMALL)
    {
        return;
    }

    const fileName caseDir
    (
        mesh_.time().rootPath()/mesh_.time().globalCaseName()
    );

    const fileName postDir(caseDir/"postProcessing");
    const fileName forceDir(postDir/"movingWallDEMForces");
    const fileName forceFile(forceDir/"force.dat");

    if (!isDir(postDir))
    {
        mkDir(postDir);
    }

    if (!isDir(forceDir))
    {
        mkDir(forceDir);
    }

    const bool writeHeader = !isFile(forceFile);
    std::ofstream output(forceFile.c_str(), std::ios::out | std::ios::app);

    if (!output.good())
    {
        FatalErrorInFunction
            << "Cannot open moving-wall DEM force file " << forceFile
            << exit(FatalError);
    }

    output << std::setprecision(12);

    if (writeHeader)
    {
        output << "# Time wall Px Py Pz Ux Uy Uz Fx Fy Fz" << std::endl;
        output << "# Force is the DEM reaction on the wall, averaged over the CFD step."
               << std::endl;
    }

    forAll(wallNames, wallI)
    {
        const wallMotionState& motionState =
            wallMotionInfo::state(wallNames[wallI]);

        const vector& planePoint = motionState.nextPlanePoint();
        const vector& wallVelocity = motionState.velocity();
        const vector meanForce(wallImpulses[wallI]/deltaTime);

        output << mesh_.time().value() << ' '
               << wallNames[wallI] << ' '
               << planePoint[0] << ' '
               << planePoint[1] << ' '
               << planePoint[2] << ' '
               << wallVelocity[0] << ' '
               << wallVelocity[1] << ' '
               << wallVelocity[2] << ' '
               << meanForce[0] << ' '
               << meanForce[1] << ' '
               << meanForce[2] << std::endl;
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::initialize
(
    volScalarField& body,
    volVectorField& U,
    volScalarField& refineF,
    label recomputeM0,
    word runTime
)
{
    if(HFDIBDEMDict_.found("outputSetup"))
    {
        dictionary outputDic = HFDIBDEMDict_.subDict("outputSetup");
        bool basicOutput = readBool(outputDic.lookup("basic"));
        bool iBoutput = readBool(outputDic.lookup("iB"));
        bool DEMoutput = readBool(outputDic.lookup("DEM"));
        bool addModelOutput = readBool(outputDic.lookup("addModel"));
        bool parallelDEMOutput = readBool(outputDic.lookup("parallelDEM"));
        InfoH.setOutput(
            basicOutput,
            iBoutput,
            DEMoutput,
            addModelOutput,
            parallelDEMOutput
        );
    }

    preCalculateCellPoints();

    if(HFDIBDEMDict_.found("interpolationSchemes"))
    {
        HFDIBinterpDict_ = HFDIBDEMDict_.subDict("interpolationSchemes");

        if(HFDIBinterpDict_.found("method"))
        {
            word intMethod = HFDIBinterpDict_.lookup("method");

            if(intMethod == "leastSquares")
            {
                dictionary lsCoeffsDict
                    = HFDIBinterpDict_.subDict("leastSquaresCoeffs");
                ibInterp_.set(new leastSquaresInt(
                    mesh_,
                    readScalar(lsCoeffsDict.lookup("distFactor")),
                    readScalar(lsCoeffsDict.lookup("radiusFactor")),
                    readScalar(lsCoeffsDict.lookup("angleFactor")),
                    readScalar(lsCoeffsDict.lookup("maxCCRows"))
                ));
            }
            else if(intMethod == "line")
            {
                ibInterp_.set(new lineInt(HFDIBinterpDict_));
            }
        }
    }

    bool startTime0(runTime == "0");

    // initialize addModels
    addModels_.setSize(bodyNames_.size());
    immersedBodies_.setSize(0);                                         //on the fly creation
    refineF *= 0;
    recomputeM0_ = recomputeM0;

    if(!startTime0)
    {
        if(!isDir(recordOutDir_))
            mkDir(recordOutDir_);
        else
        {
            fileNameList entries(readDir(recordOutDir_,fileType::directory)); // OF version 8, For version 6 use fileName::DIRECTORY instead of fileType::directory
            scalar runTimeS(stod(runTime));
            forAll(entries,entry)
            {
                scalar dirTime(stod(entries[entry].name()));
                if(dirTime > runTimeS)
                {
                    word pathI(recordOutDir_ + "/" + entries[entry]);
                    rmDir(pathI);
                }
            }
        }

        restartSimulation(body, refineF, runTime);
    }
    else
    {
        if(!isDir(recordOutDir_))
            mkDir(recordOutDir_);
        else
        {
            rmDir(recordOutDir_);
            mkDir(recordOutDir_);
        }
    }

    #include "initializeAddModels.H"

    forAll (addModels_,modelI)
    {
        word bodyName(bodyNames_[modelI]);
        InfoH << basic_Info << "Creating immersed body based on: " << bodyName << endl;

        label maxAdditions(1000);
        label cAddition(0);

        while (addModels_[modelI].shouldAddBody(body) and cAddition < maxAdditions and immersedBodies_.size() < solverInfo::getNSolidsTreshnold())
        {
            InfoH << addModel_Info << "addModel invoked action, trying to add new body" << endl;
            std::shared_ptr<geomModel> bodyGeomModel(addModels_[modelI].addBody(body, immersedBodies_));
            cAddition++;

            // initialize the immersed bodies
            if (addModels_[modelI].getBodyAdded())
            {
                label newIBSize(immersedBodies_.size()+1);
                label addIBPos(newIBSize - 1);
                immersedBodies_.setSize(newIBSize);

                InfoH << addModel_Info << "Trying to set immersedBodies" << endl;
                immersedBodies_.set
                (
                    addIBPos,
                    new immersedBody
                    (
                        bodyName,
                        mesh_,
                        HFDIBDEMDict_,
                        transportProperties_,
                        addIBPos,
                        recomputeM0_,
                        bodyGeomModel,
                        ibInterp_,
                        cellPoints_
                    )
                );
                immersedBodies_[addIBPos].createImmersedBody(body,refineF);
                immersedBodies_[addIBPos].computeBodyCharPars();
                if (immersedBodies_[addIBPos].getStartSynced())
                {
                    immersedBodies_[addIBPos].initSyncWithFlow(U);
                }
                verletList_.addBodyToVList(immersedBodies_[addIBPos]);
                InfoH << addModel_Info << "Body based on: " << bodyName << " successfully added" << endl;
                cAddition = 0;
            }
            else
            {
                InfoH << addModel_Info << "Body based on: "
                    << bodyName << " should have been added but was not "
                    << "(probably overlap with an already existing body)"
                    << endl;
            }
        }
    }

    verletList_.initialSorting();
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::createBodies(volScalarField& body,volScalarField& refineF,const bool dynamicMesh)
{
forAll (immersedBodies_,bodyId)
{
    if (immersedBodies_[bodyId].getIsActive())
    {
        immersedBodies_[bodyId].postContactUpdateBodyField
        (
            body,
            refineF,
            dynamicMesh
        );
    }
}

    DynamicList<scalar> particleMasses;
    DynamicList<label> particleCells;
    DynamicList<symmTensor> particleInertiaTensors;

    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].syncImmersedBodyParralell1(body,refineF);
            if (immersedBodies_[bodyId].getGeomModel().isCluster())
            {
                clusterBody& cBody = dynamic_cast<clusterBody&>(immersedBodies_[bodyId].getGeomModel());
                std::vector<std::shared_ptr<geomModel>>& cBodies = cBody.getClusterBodies();
                for (auto& cB : cBodies)
                {
                    particleMasses.append(cB->getM());
                    particleCells.append(cB->getNCells());
                    particleInertiaTensors.append(cB->getI());
                }
            }
            else
            {
                particleMasses.append(immersedBodies_[bodyId].getGeomModel().getM());
                particleCells.append(immersedBodies_[bodyId].getGeomModel().getNCells());
                particleInertiaTensors.append(immersedBodies_[bodyId].getGeomModel().getI());
            }
        }
    }
    reduce(particleMasses,sumOp<List<scalar>>());
    reduce(particleCells,sumOp<List<label>>());
    reduce(particleInertiaTensors,sumOp<List<symmTensor>>());

    label bodyIndex(0);
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            if (immersedBodies_[bodyId].getGeomModel().isCluster())
            {
                clusterBody& cBody = dynamic_cast<clusterBody&>(immersedBodies_[bodyId].getGeomModel());
                std::vector<std::shared_ptr<geomModel>>& cBodies = cBody.getClusterBodies();
                for (auto& cB : cBodies)
                {
                    cB->setM(particleMasses[bodyIndex]);
                    cB->setNCells(particleCells[bodyIndex]);
                    cB->setI(particleInertiaTensors[bodyIndex]);
                    bodyIndex++;
                }
                cBody.setMassAndInertia();
            }
            else
            {
                immersedBodies_[bodyId].getGeomModel().setM(particleMasses[bodyIndex]);
                immersedBodies_[bodyId].getGeomModel().setNCells(particleCells[bodyIndex]);
                immersedBodies_[bodyId].getGeomModel().setI(particleInertiaTensors[bodyIndex]);
                bodyIndex++;
            }

            immersedBodies_[bodyId].syncImmersedBodyParralell2(body,refineF);
            immersedBodies_[bodyId].checkIfInDomain(body);
            immersedBodies_[bodyId].updateOldMovementVars();
        }
    }

    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].chceckBodyOp();
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::preUpdateBodies
(
    volScalarField& body,
    volVectorField& f
)
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            // create body or compute body-fluid coupling and estimate
            // potential contacts with walls
            immersedBodies_[bodyId].inContactWithStatic(false);

            immersedBodies_[bodyId].updateOldMovementVars();
            immersedBodies_[bodyId].printStats();
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::postUpdateBodies
(
    volScalarField& body,
    volVectorField& f
)
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].clearIntpInfo();
            immersedBodies_[bodyId].postPimpleUpdateImmersedBody(body,f);
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::recreateBodies
(
    volScalarField& body,
    volScalarField& refineF,
    const bool dynamicMesh
)
{
    refineF *= 0;
    preCalculateCellPoints();
    forAll (addModels_,modelI)
    {
        addModels_[modelI].recreateBoundBox();
    }
forAll (immersedBodies_,bodyId)
{
    if (immersedBodies_[bodyId].getIsActive())
    {
        if (dynamicMesh)
        {
            immersedBodies_[bodyId].recreateBodyFieldDynamic
            (
                body,
                refineF
            );
        }
        else
        {
            immersedBodies_[bodyId].recreateBodyField
            (
                body,
                refineF
            );
        }
    }
}
    DynamicList<scalar> particleMasses;
    DynamicList<label> particleCells;
    DynamicList<symmTensor> particleInertiaTensors;
    
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].syncImmersedBodyParralell1(body,refineF);
            if (immersedBodies_[bodyId].getGeomModel().isCluster())
            {
                clusterBody& cBody = dynamic_cast<clusterBody&>(immersedBodies_[bodyId].getGeomModel());
                std::vector<std::shared_ptr<geomModel>>& cBodies = cBody.getClusterBodies();
                for (auto& cB : cBodies)
                {
                    particleMasses.append(cB->getM());
                    particleCells.append(cB->getNCells());
                    particleInertiaTensors.append(cB->getI());
                }
            }
            else
            {
                particleMasses.append(immersedBodies_[bodyId].getGeomModel().getM());
                particleCells.append(immersedBodies_[bodyId].getGeomModel().getNCells());
                particleInertiaTensors.append(immersedBodies_[bodyId].getGeomModel().getI());
            }
        }
    }
    reduce(particleMasses,sumOp<List<scalar>>());
    reduce(particleCells,sumOp<List<label>>());
    reduce(particleInertiaTensors,sumOp<List<symmTensor>>());

    label bodyIndex(0);
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            if (immersedBodies_[bodyId].getGeomModel().isCluster())
            {
                clusterBody& cBody = dynamic_cast<clusterBody&>(immersedBodies_[bodyId].getGeomModel());
                std::vector<std::shared_ptr<geomModel>>& cBodies = cBody.getClusterBodies();
                for (auto& cB : cBodies)
                {
                    cB->setM(particleMasses[bodyIndex]);
                    cB->setNCells(particleCells[bodyIndex]);
                    cB->setI(particleInertiaTensors[bodyIndex]);
                    bodyIndex++;
                }
                cBody.setMassAndInertia();
            }
            else
            {
                immersedBodies_[bodyId].getGeomModel().setM(particleMasses[bodyIndex]);
                immersedBodies_[bodyId].getGeomModel().setNCells(particleCells[bodyIndex]);
                immersedBodies_[bodyId].getGeomModel().setI(particleInertiaTensors[bodyIndex]);
                bodyIndex++;
            }

            immersedBodies_[bodyId].syncImmersedBodyParralell2(body,refineF);
            immersedBodies_[bodyId].checkIfInDomain(body);
        }
    }
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].checkIfInDomain(body);
            if(immersedBodies_[bodyId].getrecomputeM0() > 0)
            {
                immersedBodies_[bodyId].computeBodyCharPars();
                immersedBodies_[bodyId].recomputedM0();
            }
            InfoH << iB_Info << "-- body "
                << immersedBodies_[bodyId].getBodyId() << " Re-created" << endl;
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::interpolateIB( volVectorField & V
                              ,volVectorField & Vs
                              ,volScalarField & body)
{
    if(ibInterp_.valid())
    {
        ibInterp_->resetInterpolator(V);
    }
    // reset imposed field
    Vs = V;

    // loop over all the immersed bodies
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            // update imposed field according to body
            immersedBodies_[bodyId].updateVectorField(Vs, V.name(),body);

            if(ibInterp_.valid())
            {
                ibInterp_->ibInterpolate
                (
                    immersedBodies_[bodyId].getIntpInfo(),
                    Vs,
                    immersedBodies_[bodyId].getUatIbPoints(),
                    mesh_
                );
            }
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::writeBodiesInfo()
{
    if(!recordSimulation_)
        return;

    word curOutDir(recordOutDir_ + "/" + mesh_.time().timeName());


    mkDir(curOutDir);
    mkDir(curOutDir +"/stlFiles");
    DynamicLabelList activeIB;
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            activeIB.append(bodyId);
        }
    }
    wordList bodyNames;
    scalar listZize(activeIB.size());
    label bodiesPerProc = ceil(listZize/Pstream::nProcs());
    InfoH << basic_Info << "Active IB listZize      : " << listZize<< endl;
    InfoH << basic_Info << "bodiesPerProc : " << bodiesPerProc<< endl;
    // Pout << "Processor "<< Pstream::myProcNo() << endl;

    for(int assignProc = Pstream::myProcNo()*bodiesPerProc; assignProc < min((Pstream::myProcNo()+1)*bodiesPerProc,activeIB.size()); assignProc++)
    {
        const label bodyId(activeIB[assignProc]);
        // Pout <<"Processor "<< Pstream::myProcNo() << " writes Body " << bodyId << endl;
        word path(curOutDir + "/body" + std::to_string(immersedBodies_[bodyId].getBodyId()) +".info");
        OFstream ofStream(path);
        IOobject outClass
            (
                path,
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            );
        IOdictionary outDict(outClass);

        outDict.writeHeader(ofStream);
        immersedBodies_[bodyId].recordBodyInfo(outDict,curOutDir);
        outDict.writeData(ofStream);
    }

}
//---------------------------------------------------------------------------//
void openHFDIBDEM::updateDEM(volScalarField& body,volScalarField& refineF)
{
    if (cyclicPlaneInfo::getCyclicPlaneInfo().size() > 0)
    {
        forAll (immersedBodies_,bodyId)
        {
            if (!immersedBodies_[bodyId].getGeomModel().isCluster())
            {
                vector transVec = vector::zero;

                if (detectCyclicContact(
                    immersedBodies_[bodyId].getWallCntInfo(),
                    transVec
                ))
                {
                    verletList_.removeBodyFromVList(immersedBodies_[bodyId]);

                    scalar thrSurf(readScalar(HFDIBDEMDict_.lookup("surfaceThreshold")));
                    std::shared_ptr<periodicBody> newPeriodicBody
                        = std::make_shared<periodicBody>(mesh_, thrSurf);

                    newPeriodicBody->setRhoS(immersedBodies_[bodyId].getGeomModel().getRhoS());
                    std::shared_ptr<geomModel> iBcopy(immersedBodies_[bodyId].getGeomModel().getCopy());
                    iBcopy->bodyMovePoints(transVec);
                    newPeriodicBody->addBodyToCluster(immersedBodies_[bodyId].getGeomModelPtr());
                    newPeriodicBody->addBodyToCluster(iBcopy);
                    immersedBodies_[bodyId].getGeomModelPtr() = newPeriodicBody;

                    verletList_.addBodyToVList(immersedBodies_[bodyId]);
                    Info << "Periodic body created for body " << bodyId << endl;
                }
            }
            else
            {
                periodicBody& cBody = dynamic_cast<periodicBody&>(immersedBodies_[bodyId].getGeomModel());

                if(cBody.shouldBeUnclustered())
                {
                    verletList_.removeBodyFromVList(immersedBodies_[bodyId]);

                    immersedBodies_[bodyId].getGeomModelPtr() = cBody.getRemGeomModel();

                    verletList_.addBodyToVList(immersedBodies_[bodyId]);
                    Info << "Periodic body unclustered for body " << bodyId << endl;
                }
            }
        }
    }

    scalar deltaTime(mesh_.time().deltaT().value());
    scalar pos(0.0);
    scalar step(stepDEM_);
    // scalar timeStep(step*deltaTime);
    List<DynamicList<pointField>> bodiesPositionList(Pstream::nProcs());
    // Infos <<bodiesPositionList.size() << endl;
    HashTable <label,Tuple2<label, label>,Hash<Tuple2<label, label>>> syncOutForceKeyTable;
    HashTable <label,Tuple2<label, label>,Hash<Tuple2<label, label>>> contactResolvedKeyTable;
    HashTable <label,label,Hash<label>> wallContactIBTable;

    const List<string> movingWallNames
    (
        wallMotionInfo::all().sortedToc()
    );

    forAll(movingWallNames, movingWallI)
    {
        const wallMotionState& motionState =
            wallMotionInfo::state(movingWallNames[movingWallI]);

        if (motionState.timeIndex() != mesh_.time().timeIndex())
        {
            FatalErrorInFunction
                << "Moving DEM wall '" << movingWallNames[movingWallI]
                << "' was not updated for CFD time "
                << mesh_.time().timeName() << "." << nl
                << "Call updateMovingWallsFromMesh() immediately after "
                << "mesh.update() before updateDEM()."
                << exit(FatalError);
        }
    }

    List<vector> movingWallImpulses
    (
        movingWallNames.size(),
        vector::zero
    );

    while( pos < 1)
    {
        bodiesPositionList[Pstream::myProcNo()].clear();

        InfoH << DEM_Info << " Start DEM pos: " << pos
            << " DEM step: " << step << endl;

        InfoH << basic_Info << " DEM - CFD Time: "
            << mesh_.time().value() + deltaTime*pos << endl;

        forAll (immersedBodies_,ib)
        {
            immersedBodies_[ib].updateMovement(deltaTime*step*0.5);

            if(Pstream::myProcNo() == 0 )
            {
                immersedBodies_[ib].moveImmersedBody(deltaTime*step);
                if(immersedBodies_[ib].getGeomModel().getcType() != cluster)
                {
                    bodiesPositionList[Pstream::myProcNo()].append(immersedBodies_[ib].getGeomModel().getBodyPoints());
                }
                else
                {
                    clusterBody& cBody = dynamic_cast<clusterBody&>(immersedBodies_[ib].getGeomModel());
                    std::vector<std::shared_ptr<geomModel>>& cBodies = cBody.getClusterBodies();
                    for (auto& cB : cBodies)
                    {
                        bodiesPositionList[Pstream::myProcNo()].append(cB->getBodyPoints());
                    }
                }
            }
        }

        Pstream::gatherList(bodiesPositionList,0);
        Pstream::scatterList(bodiesPositionList,0);

        label bodyIndex(0);
        forAll (immersedBodies_,ib)
        {
            if(immersedBodies_[ib].getGeomModel().getcType() != cluster)
            {
                immersedBodies_[ib].getGeomModel().setBodyPosition(bodiesPositionList[0][bodyIndex++]);
            }
            else
            {
                clusterBody& cBody = dynamic_cast<clusterBody&>(immersedBodies_[ib].getGeomModel());
                std::vector<std::shared_ptr<geomModel>>& cBodies = cBody.getClusterBodies();
                for (auto& cB : cBodies)
                {
                    cB->setBodyPosition(bodiesPositionList[0][bodyIndex++]);
                }
            }
        }

        bodiesPositionList[Pstream::myProcNo()].clear();

        // Particle positions above are at the end of this DEM substep.
        // Evaluate contact against the piston at that same fractional time.
        setMovingWallsAtFraction(pos + step);

        verletList_.update(immersedBodies_);

        DynamicLabelList wallContactIB;
        wallContactIBTable.clear();
        forAll (immersedBodies_,bodyId)
        {
            immersedBody& cIb(immersedBodies_[bodyId]);
            if (cIb.getIsActive())
            {
                // set F_ and T_ to zero.
                cIb.resetContactForces();

                if(cIb.getbodyOperation() != 0)
                {
                    // detect wall contact
                    if(detectWallContact
                    (
                        mesh_,
                        cIb.getibContactClass(),
                        cIb.getWallCntInfo()
                    ))
                    {
                        cIb.getibContactClass().setWallContact(true);

                        if (cIb.getWallCntInfo().hasStaticContactPatch())
                        {
                            cIb.getibContactClass().inContactWithStatic(true);
                        }

                        wallContactIB.append(bodyId);
                        wallContactIBTable.insert(bodyId,wallContactIB.size()-1);
                        // cIb.getWallCntInfo().registerSubContactList(wallContactList);
                    }
                }
            }
        }
        // possibleWallContacts = wallContactIB.size();
        List<bool> wallContactResolvedList(wallContactIB.size(),false);

        if(wallContactIB.size() > 0)
        {
            label wallContactPerProc(ceil(double(wallContactIB.size())/Pstream::nProcs()));
            // Info <<" wallContactPerProc : "<< wallContactPerProc << endl;
            if( wallContactIB.size() <= Pstream::nProcs())
            {
                wallContactPerProc = 1;
            }
            for(int assignProc = Pstream::myProcNo()*wallContactPerProc; assignProc < min((Pstream::myProcNo()+1)*wallContactPerProc,wallContactIB.size()); assignProc++)
            {
                immersedBody& cIb(immersedBodies_[wallContactIB[assignProc]]);
                if(cIb.getGeomModel().getcType() != sphere && cIb.getGeomModel().getcType() != cluster)
                {
                    cIb.getWallCntInfo().findContactAreas();
                }

                DynamicList<wallSubContactInfo*> wallContactList;
                cIb.getWallCntInfo().registerSubContactList(wallContactList);
    
                List<bool> wallcRList(wallContactList.size(),false);

                forAll(wallContactList,sC)
                {
                    wallSubContactInfo* sCW = wallContactList[sC];
                    bool resolved(solveWallContact(
                        mesh_,
                        cIb.getWallCntInfo(),
                        deltaTime*step,
                        *sCW
                        ));

                    sCW->setResolvedContact(resolved);
                    wallContactResolvedList[assignProc] += resolved;
                    wallcRList[sC] = resolved;

                    if (resolved && sCW->getContactPatches().size() == 1)
                    {
                        const string& wallName =
                            sCW->getContactPatches()[0];

                        if (wallMotionInfo::isMoving(wallName))
                        {
                            forAll(movingWallNames, movingWallI)
                            {
                                if (movingWallNames[movingWallI] == wallName)
                                {
                                    // getOutForce() is the force on the
                                    // particle.  The wall reaction has the
                                    // opposite sign.  Accumulate impulse so
                                    // variable DEM substeps are handled.
                                    movingWallImpulses[movingWallI] -=
                                        sCW->getOutForce().F
                                       *(deltaTime*step);
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            reduce(wallContactResolvedList,sumOp<List<bool>>());

            List<vector> iBodyOutForceList(wallContactIB.size(),vector::zero);
            List<vector> iBodyOutTorqueList(wallContactIB.size(),vector::zero);

            forAll (wallContactIB,iB)
            {
                immersedBody& cIb(immersedBodies_[wallContactIB[iB]]);
                if(wallContactIBTable.found(cIb.getBodyId()))
                {
                    label cKey(wallContactIBTable[cIb.getBodyId()]);
                    if(wallContactResolvedList[cKey])
                    {
                        std::vector<std::shared_ptr<wallSubContactInfo>>& subCList
                            = cIb.getWallCntInfo().getWallSCList();

                        for(auto sCW : subCList)
                        {
                            iBodyOutForceList[cKey] += sCW->getOutForce().F;
                            iBodyOutTorqueList[cKey] += sCW->getOutForce().T;
                        }
                    }
                }
            }
            reduce(iBodyOutForceList,sumOp<List<vector>>());
            reduce(iBodyOutTorqueList,sumOp<List<vector>>());

            forAll (wallContactIB,iB)
            {
                immersedBody& cIb(immersedBodies_[wallContactIB[iB]]);
                forces cF;
                cF.F = iBodyOutForceList[iB];
                cF.T = iBodyOutTorqueList[iB];

                cIb.updateContactForces(cF);
                cIb.getWallCntInfo().clearOldContact();
            }
        }

        wallContactIB.clear();

        DynamicList<prtSubContactInfo*> contactList;
        // check only pairs whose bounding boxes are intersected for the contact
        label vListSize(0);
        for (auto it = verletList_.begin(); it != verletList_.end(); ++it)
        {
            const Tuple2<label, label> cPair = Tuple2<label, label>(it->first, it->second);

            label cInd(cPair.first());
            bool cStatic(immersedBodies_[cInd].getbodyOperation() == 0);

            label tInd(cPair.second());
            bool tStatic(immersedBodies_[tInd].getbodyOperation() == 0);

            if((immersedBodies_[cInd].getIsActive() && immersedBodies_[tInd].getIsActive())
                &&
                !(cStatic && tStatic)
            )
            {
                if(cStatic)
                    immersedBodies_[tInd].inContactWithStatic(true);

                if(tStatic)
                    immersedBodies_[cInd].inContactWithStatic(true);

                prtContactInfo& prtcInfo(getPrtcInfo(
                    cPair)
                );

                prtcInfo.clearData();
                getContacts(
                    mesh_,
                    prtcInfo
                );

                // prtcInfo.syncContactList();

                prtcInfo.registerContactList(contactList);
            }
            vListSize++;
        }

        List<bool> contactResolved(contactList.size(),false);
        List<label> contactResolvedcKey(contactList.size(),0);
        List<label> contactResolvedtKey(contactList.size(),0);
        bool syncedData(true);
        reduce(syncedData, orOp<bool>());

        if(contactList.size() > 0 )
        {
            label contactPerProc(ceil(double(contactList.size())/Pstream::nProcs()));
            if( contactList.size() <= Pstream::nProcs())
            {
                contactPerProc = 1;
            }

            for(int assignProc = Pstream::myProcNo()*contactPerProc; assignProc < min((Pstream::myProcNo()+1)*contactPerProc,contactList.size()); assignProc++)
            {
                prtSubContactInfo* sCI = contactList[assignProc];
                const Tuple2<label, label>& cPair = sCI->getCPair();

                contactResolvedcKey[assignProc] = cPair.first();
                contactResolvedtKey[assignProc] = cPair.second();

                ibContactClass& cClass(immersedBodies_[cPair.first()].getibContactClass());
                ibContactClass& tClass(immersedBodies_[cPair.second()].getibContactClass());

                if(detectPrtPrtContact(mesh_,cClass,tClass,*sCI))
                {
                    // resolvedPrtContacts++;
                    prtContactInfo& prtcInfo(getPrtcInfo(cPair));

                    bool resolved(solvePrtContact(mesh_, prtcInfo, *sCI, deltaTime*step));
                    sCI->setResolvedContact(resolved);

                    contactResolved[assignProc] += resolved;
                }
            }
        }

        reduce(contactResolved,sumOp<List<bool>>());
        reduce(contactResolvedcKey,sumOp<List<label>>());
        reduce(contactResolvedtKey,sumOp<List<label>>());
        contactResolvedKeyTable.clear();
        forAll(contactResolvedcKey,cKey)
        {
            contactResolvedKeyTable.insert(Tuple2<label, label>(contactResolvedcKey[cKey],contactResolvedtKey[cKey]),cKey);
        }
        List<vector> cBodyOutForceList(vListSize,vector::zero);
        List<vector> cBodyOutTorqueList(vListSize,vector::zero);
        List<vector> tBodyOutForceList(vListSize,vector::zero);
        List<vector> tBodyOutTorqueList(vListSize,vector::zero);
        syncOutForceKeyTable.clear();

        label nIter(0);
        for (auto it = verletList_.begin(); it != verletList_.end(); ++it)
        {
            const Tuple2<label, label> cPair = Tuple2<label, label>(it->first, it->second);

            prtContactInfo& prtcInfo(getPrtcInfo(cPair));

            if(contactResolvedKeyTable.found(cPair))
            {
                label nSubContact(0);
                std::vector<std::shared_ptr<prtSubContactInfo>>& subCList
                    = prtcInfo.getPrtSCList();
                for(auto sC : subCList)
                {
                    nSubContact++;
                    cBodyOutForceList[nIter] += sC->getOutForce().first().F;
                    cBodyOutTorqueList[nIter] += sC->getOutForce().first().T;
                    tBodyOutForceList[nIter] += sC->getOutForce().second().F;
                    tBodyOutTorqueList[nIter] += sC->getOutForce().second().T;
                }
            }
            syncOutForceKeyTable.insert(cPair,nIter);
            nIter++;
        }

        reduce(cBodyOutForceList,sumOp<List<vector>>());
        reduce(cBodyOutTorqueList,sumOp<List<vector>>());
        reduce(tBodyOutForceList,sumOp<List<vector>>());
        reduce(tBodyOutTorqueList,sumOp<List<vector>>());

        label nvListIter(0);

        for (auto it = verletList_.begin(); it != verletList_.end(); ++it)
        {
            const Tuple2<label, label> cPair = Tuple2<label, label>(it->first, it->second);
            label cInd(cPair.first());
            label tInd(cPair.second());

            if(!contactResolvedKeyTable.found(cPair))
            {
                if(prtcInfoTable_.found(cPair))
                {
                    prtcInfoTable_.erase(cPair);
                    continue;
                }
            }
            else if(contactResolved[contactResolvedKeyTable[cPair]])
            {
                if(!syncOutForceKeyTable.found(cPair))
                {
                    Pout <<" -- cPair  "<<cInd << " - "<<tInd << " not found in syncOutForceKeyTable" << endl;
                    continue;
                }

                nvListIter = syncOutForceKeyTable[cPair];
                if(nvListIter > cBodyOutForceList.size())
                {
                    Pout <<" -- cPair  "<<cInd << " - "<<tInd << " nvListIter > bodiesOutForceList[Pstream::myProcNo()].size()" << endl;
                    continue;
                }

                vector F1 = vector::zero;
                vector T1 = vector::zero;
                vector F2 = vector::zero;
                vector T2 = vector::zero;

                F1 += cBodyOutForceList[nvListIter];
                T1 += cBodyOutTorqueList[nvListIter];
                F2 += tBodyOutForceList[nvListIter];
                T2 += tBodyOutTorqueList[nvListIter];

                forces cF;
                cF.F = F1;
                cF.T = T1;
                forces tF;
                tF.F = F2;
                tF.T = T2;

                immersedBodies_[cInd].updateContactForces(cF);
                immersedBodies_[tInd].updateContactForces(tF);
            }
            else
            {
                if(prtcInfoTable_.found(cPair))
                {
                    prtcInfoTable_.erase(cPair);
                    continue;
                }
            }
        }

        forAll (immersedBodies_,ib)
        {
            immersedBodies_[ib].updateMovement(deltaTime*step*0.5);
            immersedBodies_[ib].printBodyInfo();
            // immersedBodies_[ib].computeBodyCoNumber();
            // if (maxCoNum < immersedBodies_[ib].getCoNum())
            // {
                // maxCoNum = immersedBodies_[ib].getCoNum();
                // bodyId = ib;
            // }
        }
        // InfoH << basic_Info << "Max CoNum = " << maxCoNum << " at body " << bodyId << endl;

        pos += step;

        if (pos + step + SMALL >= 1)
            step = 1 - pos;
//OS Time effitiency Testing
        // demItegrationTime_ = DEMIntergrationRun.timeIncrement();
//OS Time effitiency Testing
    }

    setMovingWallsAtFraction(1.0);

    if (movingWallNames.size() > 0)
    {
        reduce(movingWallImpulses, sumOp<List<vector>>());
        writeMovingWallDEMForces
        (
            movingWallNames,
            movingWallImpulses,
            deltaTime
        );
    }
}
//---------------------------------------------------------------------------//
prtContactInfo& openHFDIBDEM::getPrtcInfo(Tuple2<label,label> cPair)
{
    if(!prtcInfoTable_.found(cPair))
    {
        prtcInfoTable_.insert(cPair, autoPtr<prtContactInfo>( new prtContactInfo(
            immersedBodies_[cPair.first()].getibContactClass(),
            immersedBodies_[cPair.first()].getContactVars(),
            immersedBodies_[cPair.second()].getibContactClass(),
            immersedBodies_[cPair.second()].getContactVars()
        )));
    }

    return prtcInfoTable_[cPair]();
}
//---------------------------------------------------------------------------//
// function to either add or remove bodies from the simulation
void openHFDIBDEM::addRemoveBodies
(
    volScalarField& body,
    volVectorField& U,
    volScalarField& refineF
)
{
    forAll (addModels_,modelI)
    {
        word bodyName(bodyNames_[modelI]);

        label maxAdditions(50);
        label cAddition(0);

        while (addModels_[modelI].shouldAddBody(body) and cAddition < maxAdditions)
        {
            InfoH << addModel_Info << "addModel invoked action, trying to add new body" << endl;
            std::shared_ptr<geomModel> bodyGeomModel(addModels_[modelI].addBody(body, immersedBodies_));

            cAddition++;

            if (addModels_[modelI].getBodyAdded())
            {
                InfoH << addModel_Info << "STL file correctly generated, registering the new body" << endl;

                // prepare pointer list for IBs (increase its size)
                label newIBSize(immersedBodies_.size()+1);
                label addIBPos(newIBSize - 1);
                immersedBodies_.setSize(newIBSize);

                // create the new body
                immersedBodies_.set
                (
                    addIBPos,
                    new immersedBody
                    (
                        bodyName,
                        mesh_,
                        HFDIBDEMDict_,
                        transportProperties_,
                        addIBPos,
                        recomputeM0_,
                        bodyGeomModel,
                        ibInterp_,
                        cellPoints_
                    )
                );

                // get reference for further processing
                immersedBody& nBody(immersedBodies_[addIBPos]);
                nBody.createImmersedBody(body,refineF);
                nBody.computeBodyCharPars();
                if (nBody.getStartSynced())
                {
                    nBody.initSyncWithFlow(U);
                }
                verletList_.addBodyToVList(nBody);

                InfoH << addModel_Info
                    << "new body included into the simulation" << endl;
                cAddition = 0;
            }
            else
            {
                InfoH << addModel_Info
                    << "new body should have been added but was not "
                    << "(probably overlap with an existing body)"
                    << endl;
            }
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::updateFSCoupling
(
    volScalarField& body,
    volVectorField& f
)
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].pimpleUpdate(body,f);
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::restartSimulation
(
    volScalarField& body,
    volScalarField& refineF,
    word runTime
)
{
    word timePath(recordOutDir_+"/"+runTime);
    fileNameList files(readDir(timePath));
    scalar thrSurf(readScalar(HFDIBDEMDict_.lookup("surfaceThreshold")));

    forAll(files,f)
    {
        IOdictionary bodyDict
        (
            IOobject
            (
                timePath + "/" + files[f],
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        word bodyId(std::to_string(readLabel(bodyDict.lookup("bodyId"))));
        word bodyName(bodyDict.lookup("bodyName"));
        vector Vel(bodyDict.lookup("Vel"));
        scalar omega(readScalar(bodyDict.lookup("omega")));
        vector Axis(bodyDict.lookup("Axis"));
        bool isStatic(readBool(bodyDict.lookup("static")));
        label timeStepsInContWStatic(readLabel(bodyDict.lookup("timeStepsInContWStatic")));

        std::shared_ptr<geomModel> bodyGeomModel;
        word bodyGeom;
        // check if the immersedDict_ contains bodyGeom
        if (HFDIBDEMDict_.subDict(bodyName).found("bodyGeom"))
        {
            word input = HFDIBDEMDict_.subDict(bodyName).lookup("bodyGeom");
            bodyGeom = input;
            InfoH << iB_Info << "Found bodyGeom for "
                << bodyName << ", the body is: " << bodyGeom << endl;
        }
        else
        {
            bodyGeom = "convex";
            InfoH << iB_Info << "Did not find bodyGeom for "
                << bodyName << ", using bodyGeom: " << bodyGeom << endl;
        }

        if(bodyGeom == "convex")
        {
            word stlPath(timePath + "/stlFiles/"+bodyId+".stl");
            bodyGeomModel = std::make_shared<convexBody>(mesh_,stlPath,thrSurf);
        }
        else if(bodyGeom == "nonConvex")
        {
            word stlPath(timePath + "/stlFiles/"+bodyId+".stl");
            bodyGeomModel = std::make_shared<nonConvexBody>(mesh_,stlPath,thrSurf);
        }
        else if(bodyGeom == "sphere")
        {
            vector startPosition = bodyDict.subDict("sphere").lookup("position");
            scalar radius = readScalar(bodyDict.subDict("sphere").lookup("radius"));

            bodyGeomModel = std::make_shared<sphereBody>(mesh_,startPosition,radius,thrSurf);
        }
        else
        {
            word stlPath(timePath + "/stlFiles/"+bodyId+".stl");
            InfoH << iB_Info << "bodyGeom: " << bodyGeom
                << " not supported, using bodyGeom nonConvex" << endl;
            bodyGeom = "nonConvex";
            bodyGeomModel = std::make_shared<nonConvexBody>(mesh_,stlPath,thrSurf);
        }

        label newIBSize(immersedBodies_.size()+1);
        label addIBPos(newIBSize - 1);
        immersedBodies_.setSize(newIBSize);

        InfoH << iB_Info << "Restarting body: " << bodyId << " as "
            << addIBPos << " bodyName: " << bodyName << endl;
        immersedBodies_.set
        (
            addIBPos,
            new immersedBody
            (
                bodyName,
                mesh_,
                HFDIBDEMDict_,
                transportProperties_,
                addIBPos,
                recomputeM0_,
                bodyGeomModel,
                ibInterp_,
                cellPoints_
            )
        );

        immersedBodies_[addIBPos].createImmersedBody(body,refineF);
        immersedBodies_[addIBPos].computeBodyCharPars();
        immersedBodies_[addIBPos].setRestartSim(Vel,omega,Axis,isStatic,timeStepsInContWStatic);
        verletList_.addBodyToVList(immersedBodies_[addIBPos]);
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::preCalculateCellPoints()
{
    cellPoints_.clear();
    cellPoints_.setSize(mesh_.nCells());
    forAll(mesh_.C(), cellI)
    {
        cellPoints_[cellI] = mesh_.cellPoints()[cellI];
    }

    forAll (immersedBodies_,bodyId)
    {
        immersedBodies_[bodyId].getGeomModel().resetHashTable();
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::writeFirtsTimeBodiesInfo()
{
    word curOutDir(recordOutDir_ + "/" + mesh_.time().timeName());
    bool checkExistance(false);
    if(!recordSimulation_ || isDir(curOutDir))
        return;
    reduce(checkExistance,orOp<bool>());
    if(Pstream::myProcNo() == 0)
    {
        mkDir(curOutDir);
        mkDir(curOutDir +"/stlFiles");
    }
    reduce(checkExistance,orOp<bool>());

    DynamicLabelList activeIB;
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            activeIB.append(bodyId);
        }
    }

    wordList bodyNames;
    scalar listZize(activeIB.size());
    label bodiesPerProc = ceil(listZize/Pstream::nProcs());
    InfoH << basic_Info << "Active IB listZize      : " << listZize<< endl;
    InfoH << basic_Info << "bodiesPerProc : " << bodiesPerProc<< endl;

    for(int assignProc = Pstream::myProcNo()*bodiesPerProc; assignProc < min((Pstream::myProcNo()+1)*bodiesPerProc,activeIB.size()); assignProc++)
    {
        const label bodyId(activeIB[assignProc]);
        word path(curOutDir + "/body" + std::to_string(immersedBodies_[bodyId].getBodyId()) +".info");
        OFstream ofStream(path);
        IOobject outClass
            (
                path,
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            );
        IOdictionary outDict(outClass);

        outDict.writeHeader(ofStream);
        immersedBodies_[bodyId].recordBodyInfo(outDict,curOutDir);
        outDict.writeData(ofStream);
    }

}
//---------------------------------------------------------------------------//
void openHFDIBDEM::setSolverInfo()
{
    solverInfo::setOnlyDEM(true);
}
//---------------------------------------------------------------------------//
