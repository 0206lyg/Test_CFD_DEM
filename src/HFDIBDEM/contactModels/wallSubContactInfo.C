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
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/
#include "wallSubContactInfo.H"

#include "interAdhesion.H"
#include "wallMatInfo.H"

#include "virtualMeshLevel.H"
#include "wallPlaneInfo.H"
#include "contactModelInfo.H"
using namespace Foam;

namespace
{

scalar minScalarWSCI(const scalar a, const scalar b)
{
    return a < b ? a : b;
}

scalar maxScalarWSCI(const scalar a, const scalar b)
{
    return a > b ? a : b;
}

bool finiteWallPatchActiveWSCI(const string& wallName)
{
    const HashTable<bool,string,Hash<string>>& finiteInfo =
        wallPlaneInfo::getWallFiniteInfo();

    const HashTable<vector,string,Hash<string>>& minInfo =
        wallPlaneInfo::getWallMinBoundInfo();

    const HashTable<vector,string,Hash<string>>& maxInfo =
        wallPlaneInfo::getWallMaxBoundInfo();

    return
    (
        finiteInfo.found(wallName)
     && finiteInfo[wallName]
     && minInfo.found(wallName)
     && maxInfo.found(wallName)
    );
}

bool singleFiniteWallPatchWSCI
(
    const List<string>& contactPatches,
    string& wallName
)
{
    if (contactPatches.size() != 1)
    {
        return false;
    }

    wallName = contactPatches[0];

    return finiteWallPatchActiveWSCI(wallName);
}

label wallNormalDirectionWSCI(const string& wallName)
{
    const vector& nVec = wallPlaneInfo::getWallPlaneInfo()[wallName][0];

    const scalar ax = mag(nVec[0]);
    const scalar ay = mag(nVec[1]);
    const scalar az = mag(nVec[2]);

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

void clipBBToFinitePatchSupportWSCI
(
    boundBox& bb,
    const string& wallName
)
{
    const vector& minBound =
        wallPlaneInfo::getWallMinBoundInfo()[wallName];

    const vector& maxBound =
        wallPlaneInfo::getWallMaxBoundInfo()[wallName];

    const label nDir = wallNormalDirectionWSCI(wallName);

    for (label dir = 0; dir < 3; dir++)
    {
        if (dir == nDir)
        {
            continue;
        }

        const scalar lo =
            minScalarWSCI(minBound[dir], maxBound[dir]);

        const scalar hi =
            maxScalarWSCI(minBound[dir], maxBound[dir]);

        bb.min()[dir] = maxScalarWSCI(bb.min()[dir], lo);
        bb.max()[dir] = minScalarWSCI(bb.max()[dir], hi);
    }
}

bool finitePatchTangentialSpanValidWSCI
(
    const boundBox& bb,
    const string& wallName
)
{
    const label nDir = wallNormalDirectionWSCI(wallName);

    for (label dir = 0; dir < 3; dir++)
    {
        if (dir == nDir)
        {
            continue;
        }

        if ((bb.max()[dir] - bb.min()[dir]) <= VSMALL)
        {
            return false;
        }
    }

    return true;
}

scalar componentProductWSCI(const vector& nVec)
{
    return
    (
        maxScalarWSCI(nVec[0], scalar(1))
       *maxScalarWSCI(nVec[1], scalar(1))
       *maxScalarWSCI(nVec[2], scalar(1))
    );
}

vector generalFiniteSubVolumeNVectorWSCI
(
    const boundBox& bb,
    const bool allowNormalPadding,
    boundBox& correctedBB,
    bool& valid
)
{
    correctedBB = bb;
    valid = true;

    const scalar h = virtualMeshLevel::getCharCellSize();
    const scalar level = virtualMeshLevel::getLevelOfDivision();
    const scalar subH = h/level;

    vector subVolumeNVector(vector::zero);

    for (label dir = 0; dir < 3; dir++)
    {
        const scalar spanI =
            correctedBB.max()[dir] - correctedBB.min()[dir];

        if (spanI <= VSMALL)
        {
            if (allowNormalPadding)
            {
                subVolumeNVector[dir] = 1;

                correctedBB.min()[dir] -= subH*0.5;
                correctedBB.max()[dir] += subH*0.5;
            }
            else
            {
                valid = false;
                subVolumeNVector[dir] = 0;
            }
        }
        else
        {
            const scalar nRaw = spanI/subH;

            subVolumeNVector[dir] =
                maxScalarWSCI
                (
                    scalar(1),
                    ceil(nRaw - SMALL)
                );
        }
    }

    if (correctedBB.volume() <= VSMALL)
    {
        valid = false;
    }

    return subVolumeNVector;
}

vector legacyFiniteSubVolumeNVectorWSCI
(
    const boundBox& bb,
    const string& wallName,
    const bool allowNormalPadding,
    boundBox& correctedBB,
    bool& valid
)
{
    correctedBB = bb;
    valid = true;

    const scalar h = virtualMeshLevel::getCharCellSize();
    const scalar level = virtualMeshLevel::getLevelOfDivision();
    const scalar subH = h/level;

    vector subVolumeNVector(vector::zero);
    const label nDir = wallNormalDirectionWSCI(wallName);

    for (label dir = 0; dir < 3; dir++)
    {
        const scalar spanI =
            correctedBB.max()[dir] - correctedBB.min()[dir];

        if (spanI <= VSMALL)
        {
            if (allowNormalPadding && dir == nDir)
            {
                subVolumeNVector[dir] = 1;
                correctedBB.min()[dir] -= subH*0.5;
                correctedBB.max()[dir] += subH*0.5;
            }
            else
            {
                valid = false;
                subVolumeNVector[dir] = 0;
            }
        }
        else
        {
            const scalar nRaw = spanI/subH;

            subVolumeNVector[dir] =
                maxScalarWSCI(scalar(1), ceil(nRaw - SMALL));
        }
    }

    clipBBToFinitePatchSupportWSCI(correctedBB, wallName);

    if (!finitePatchTangentialSpanValidWSCI(correctedBB, wallName))
    {
        valid = false;
    }

    if (correctedBB.volume() <= VSMALL)
    {
        valid = false;
    }

    return subVolumeNVector;
}

} // End anonymous namespace


//---------------------------------------------------------------------------//
wallSubContactInfo::wallSubContactInfo
(
    List<Tuple2<point,boundBox>> contactBBData,
    List<Tuple2<point,boundBox>> planeBBData,
    List<string> contactPatches,
    List<Tuple2<point,boundBox>> internalBBData,
    HashTable<physicalProperties,string,Hash<string>> wallMeanPars,
    boundBox BB,
    label bodyId
)
:
contactPatches_(contactPatches),
internalBBData_(internalBBData),
wallMeanPars_(wallMeanPars),
BB_(BB),
bodyId_(bodyId)
{
    string finiteWallName;
    const bool finiteSinglePatch =
        singleFiniteWallPatchWSCI(contactPatches, finiteWallName);

    const bool legacyAxisAlignedFinitePatch =
        finiteSinglePatch
     && wallPlaneInfo::usesLegacyAxisAlignedFinitePath(finiteWallName);

    const scalar h = virtualMeshLevel::getCharCellSize();
    const scalar level = virtualMeshLevel::getLevelOfDivision();
    const scalar subH = h/level;

    forAll(contactBBData,cBD)
    {
        vector subVolumeNVector(vector::zero);
        scalar subVolumeV = pow(subH, 3);

        if (finiteSinglePatch)
        {
            if (legacyAxisAlignedFinitePatch)
            {
                // Preserve the complete Git-main fast-path preprocessing.
                clipBBToFinitePatchSupportWSCI
                (
                    contactBBData[cBD].second(),
                    finiteWallName
                );
            }

            bool validFiniteBB = true;
            boundBox correctedBB;

            if (legacyAxisAlignedFinitePatch)
            {
                subVolumeNVector =
                    legacyFiniteSubVolumeNVectorWSCI
                    (
                        contactBBData[cBD].second(),
                        finiteWallName,
                        false,
                        correctedBB,
                        validFiniteBB
                    );
            }
            else
            {
                subVolumeNVector =
                    generalFiniteSubVolumeNVectorWSCI
                (
                    contactBBData[cBD].second(),
                    false,
                    correctedBB,
                    validFiniteBB
                );
            }

            if (!validFiniteBB)
            {
                continue;
            }

            contactBBData[cBD].second() = correctedBB;
            contactBBData[cBD].first() = correctedBB.midpoint();

            subVolumeV =
                correctedBB.volume()
               /componentProductWSCI(subVolumeNVector);
        }
        else
        {
            subVolumeNVector = vector(
                floor((contactBBData[cBD].second().span()[0]/h)*level),
                floor((contactBBData[cBD].second().span()[1]/h)*level),
                floor((contactBBData[cBD].second().span()[2]/h)*level)
            );

            if(cmptMin(subVolumeNVector)<SMALL)
            {
                // Pout <<" Trubble with subVolumeNVector "<< subVolumeNVector << endl;
                for(int i=0;i<3;i++)
                {
                    if(subVolumeNVector[i] <SMALL)
                    {                
                        subVolumeNVector[i] = 1;
                        contactBBData[cBD].second().min()[i] -= subH*0.5;
                        contactBBData[cBD].second().max()[i] += subH*0.5;
                    }  
                }
                // Pout <<" Corrected subVolumeNVector "<< subVolumeNVector << endl;
            }
        }

        autoPtr<virtualMeshWallInfo> vmWInfo(
            new virtualMeshWallInfo(
                contactBBData[cBD].second(),
                contactBBData[cBD].first(),
                subVolumeNVector,
                h,
                subVolumeV,
                finiteSinglePatch,
                finiteSinglePatch ? finiteWallName : string("")
            )  
        );
        vmWInfoList_.append(vmWInfo);
    }

    forAll(planeBBData,pBD)
    {
        vector subVolumeNVector(vector::zero);
        scalar subVolumeV = pow(subH, 3);

        if (finiteSinglePatch)
        {
            if (legacyAxisAlignedFinitePatch)
            {
                // The original plane box is clipped before one-cell normal
                // padding and fitted-grid construction.
                clipBBToFinitePatchSupportWSCI
                (
                    planeBBData[pBD].second(),
                    finiteWallName
                );
            }

            bool validFinitePlaneBB = true;
            boundBox correctedPlaneBB;

            if (legacyAxisAlignedFinitePatch)
            {
                subVolumeNVector =
                    legacyFiniteSubVolumeNVectorWSCI
                    (
                        planeBBData[pBD].second(),
                        finiteWallName,
                        true,
                        correctedPlaneBB,
                        validFinitePlaneBB
                    );
            }
            else
            {
                subVolumeNVector =
                    generalFiniteSubVolumeNVectorWSCI
                (
                    planeBBData[pBD].second(),
                    true,
                    correctedPlaneBB,
                    validFinitePlaneBB
                );
            }

            if (!validFinitePlaneBB)
            {
                continue;
            }

            planeBBData[pBD].second() = correctedPlaneBB;
            planeBBData[pBD].first() = correctedPlaneBB.midpoint();

            subVolumeV =
                correctedPlaneBB.volume()
               /componentProductWSCI(subVolumeNVector);
        }
        else
        {
            subVolumeNVector = vector(
                ceil((planeBBData[pBD].second().span()[0]/h))*level,
                ceil((planeBBData[pBD].second().span()[1]/h))*level,
                ceil((planeBBData[pBD].second().span()[2]/h))*level
            );

            for(int i=0;i<3;i++)
            {
                if(subVolumeNVector[i] == planeBBData[pBD].second().minDim())
                {                
                    subVolumeNVector[i] = 1;
                    planeBBData[pBD].second().min()[i] -= subH*0.5;
                    planeBBData[pBD].second().max()[i] += subH*0.5;
                }  
            }
        }

        autoPtr<virtualMeshWallInfo> vmWInfo(
            new virtualMeshWallInfo(
                planeBBData[pBD].second(),
                planeBBData[pBD].first(),
                subVolumeNVector,
                h,
                subVolumeV,
                finiteSinglePatch,
                finiteSinglePatch ? finiteWallName : string("")
            )
        );
        vmPlaneInfoList_.append(vmWInfo);
    }
}

wallSubContactInfo::~wallSubContactInfo()
{
}
//---------------------------------------------------------------------------//
vector wallSubContactInfo::getLVec(wallContactVars& wallCntvar, ibContactClass ibCClass)
{
    return ibCClass.getGeomModel().getLVec(wallCntvar.contactCenter_);
}
//---------------------------------------------------------------------------//
vector wallSubContactInfo::getVeli(wallContactVars& wallCntvar, ibContactVars& cVars)
{
    return (-((wallCntvar.lVec_-cVars.Axis_
        *((wallCntvar.lVec_) & cVars.Axis_))
        ^ cVars.Axis_)*cVars.omega_+ cVars.Vel_);
}
//---------------------------------------------------------------------------//
void wallSubContactInfo::evalVariables(
    wallContactVars& wallCntvar,
    ibContactClass& ibCClass,
    ibContactVars& cVars
)
{
    reduceM_ =
    (
        ibCClass.getGeomModel().getM0()
        *ibCClass.getGeomModel().getM0()
        /(ibCClass.getGeomModel().getM0()
        +ibCClass.getGeomModel().getM0())
    );

    wallCntvar.lVec_ = getLVec(wallCntvar,ibCClass);
    // wallCntvar.lVec_ = wallCntvar.contactCenter_ - ibCClass.getGeomModel().getCoM();
    wallCntvar.Veli_ = getVeli(wallCntvar, cVars);

    wallCntvar.relativeVelocity_ =
        wallCntvar.Veli_ - wallCntvar.wallVelocity_;

    wallCntvar.Vn_ =
        -(wallCntvar.relativeVelocity_ & wallCntvar.contactNormal_);
    wallCntvar.Lc_ = (contactModelInfo::getLcCoeff())*mag(wallCntvar.lVec_)*mag(wallCntvar.lVec_)/(mag(wallCntvar.lVec_) + mag(wallCntvar.lVec_));
    

    wallCntvar.curAdhN_ = min
    (
        wallCntvar.getMeanCntPar().maxAdhN_,
        max(wallCntvar.curAdhN_, wallCntvar.getMeanCntPar().aY_
            *wallCntvar.contactVolume_
            /(sqr(wallCntvar.Lc_)*8*Foam::constant::mathematical::pi))
    );
}
//---------------------------------------------------------------------------//
vector wallSubContactInfo::getFNe(wallContactVars& wallCntvar)
{
    return (wallCntvar.getMeanCntPar().aY_*wallCntvar.contactVolume_
        /(wallCntvar.Lc_+SMALL))*wallCntvar.contactNormal_;
}
//---------------------------------------------------------------------------//
vector wallSubContactInfo::getFA(wallContactVars& wallCntvar)
{
    return ((sqrt(8*Foam::constant::mathematical::pi
        *wallCntvar.getMeanCntPar().aY_
        *wallCntvar.curAdhN_*wallCntvar.contactVolume_))
        *wallCntvar.contactNormal_);
}
//---------------------------------------------------------------------------//
vector wallSubContactInfo::getFNd(wallContactVars& wallCntvar)
{
    physicalProperties& meanCntPar(wallCntvar.getMeanCntPar());
    return ((meanCntPar.reduceBeta_*sqrt(meanCntPar.aY_
            *reduceM_*wallCntvar.contactArea_/(wallCntvar.Lc_+SMALL))*
            wallCntvar.Vn_)*wallCntvar.contactNormal_);

}
//---------------------------------------------------------------------------//
vector wallSubContactInfo::getFt(wallContactVars& wallCntvar, scalar deltaT)
{
    physicalProperties& meanCntPar(wallCntvar.getMeanCntPar());
    // project last Ft into a new direction
    vector FtLastP(wallCntvar.FtPrev_
        - (wallCntvar.FtPrev_ & wallCntvar.contactNormal_)
        *wallCntvar.contactNormal_);
    // scale projected Ft to have same magnitude as FtLast
    vector FtLastS(mag(wallCntvar.FtPrev_) * (FtLastP/(mag(FtLastP)+SMALL)));
    // Orthogonal projection of relative velocity onto the wall normal.
    const vector relativeVelocityNorm =
        wallCntvar.contactNormal_
       *(wallCntvar.relativeVelocity_ & wallCntvar.contactNormal_);

    const vector Vt =
        wallCntvar.relativeVelocity_ - relativeVelocityNorm;
    // compute tangential force
    if(contactModelInfo::getUseMindlinRotationalModel())
    {
        
        scalar kT = 200*8*meanCntPar.aG_*(wallCntvar.contactArea_/(wallCntvar.Lc_+SMALL));
        vector deltaFt(kT*Vt*deltaT + 2*meanCntPar.reduceBeta_*sqrt(kT*reduceM_)*Vt);
        wallCntvar.FtPrev_ = - FtLastS - deltaFt;
    }

    if(contactModelInfo::getUseChenRotationalModel())
    {
   
        vector Ftdi(meanCntPar.reduceBeta_*sqrt(meanCntPar.aG_*reduceM_*wallCntvar.Lc_)*Vt);
        Ftdi += meanCntPar.aG_*wallCntvar.Lc_*Vt*deltaT;
        wallCntvar.FtPrev_ = - FtLastS - Ftdi;
    }

    return wallCntvar.FtPrev_;
}
//---------------------------------------------------------------------------//
void wallSubContactInfo::syncData()
{
    reduce(outForce_.F, sumOp<vector>());
    reduce(outForce_.T, sumOp<vector>());
}
//---------------------------------------------------------------------------//
void wallSubContactInfo::syncContactResolve()
{
    reduce(contactResolved_,orOp<bool>());
}
//---------------------------------------------------------------------------//
autoPtr<virtualMeshWallInfo>& wallSubContactInfo::getVMContactInfo
(
    label ID
)
{
    return vmWInfoList_[ID];
}
//---------------------------------------------------------------------------//
autoPtr<virtualMeshWallInfo>& wallSubContactInfo::getVMPlaneInfo
(
    label ID
)
{
    return vmPlaneInfoList_[ID];
}

//---------------------------------------------------------------------------//
