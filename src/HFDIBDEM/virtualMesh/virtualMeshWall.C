/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| "_ \ / _ \ "_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
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
#include "virtualMeshWall.H"

#include "virtualMeshLevel.H"

using namespace Foam;

//---------------------------------------------------------------------------//
virtualMeshWall::virtualMeshWall
(
    virtualMeshWallInfo& vMeshWallInfo,
    geomModel& cGeomModel
)
:
cGeomModel_(cGeomModel),
vMeshWallInfo_(vMeshWallInfo),
bbMatrix_(vMeshWallInfo.subVolumeNVector,
    vMeshWallInfo.bBox,
    vMeshWallInfo.charCellSize,
    vMeshWallInfo.subVolumeV,
    vMeshWallInfo.fitToBBox)
{}

virtualMeshWall::~virtualMeshWall()
{
}
//---------------------------------------------------------------------------//
bool virtualMeshWall::detectFirstContactPoint()
{
    if (vMeshWallInfo_.usesFittedBBoxGrid())
    {
        // A fitted finite-wall VM is a complete bounded integration domain.
        // It needs no seed-to-cell discovery; scan its valid index range.
        return evaluateContact() > VSMALL;
    }

    autoPtr<DynamicVectorList> nextToCheck(
        new DynamicVectorList);

    autoPtr<DynamicVectorList> auxToCheck(
        new DynamicVectorList);

    nextToCheck->append(bbMatrix_.getSVIndexForPoint_Wall(vMeshWallInfo_.getStartingPoint()));
    nextToCheck->append
    (
        bbMatrix_.cornerNeighbourSubVolumes(nextToCheck()[0])
    );
    // InfoH << DEM_Info << " -- VM firstSV : " << nextToCheck()[0] << " point " << bbMatrix_[nextToCheck()[0]].center << endl;
    while (nextToCheck->size() > 0)
    {
        auxToCheck->clear();
        forAll (nextToCheck(),sV)
        {
            subVolumeProperties& cSubVolume = bbMatrix_[nextToCheck()[sV]];
            if (!cSubVolume.toCheck)
            {
                continue;
            }
            checkSubVolume(cSubVolume, nextToCheck()[sV]);

            if (cSubVolume.isCBody)
            {
                vMeshWallInfo_.startingPoint = cSubVolume.center;
                resetSubVolume(cSubVolume);

                return true;

            }
            auxToCheck().append
            (
                bbMatrix_.cornerNeighbourSubVolumes(nextToCheck()[sV])
            );
        }
        const autoPtr<DynamicVectorList> helpPtr(nextToCheck.ptr());
        nextToCheck.set(auxToCheck.ptr());
        auxToCheck = helpPtr;
    }
    return false;
}
//---------------------------------------------------------------------------////---------------------------------------------------------------------------//
bool virtualMeshWall::detectFirstFaceContactPoint()
{
    if (vMeshWallInfo_.usesFittedBBoxGrid())
    {
        // The finite plane VM is likewise a complete one-cell-thick slab.
        // Avoid the legacy starting-point lookup used by infinite walls.
        return evaluateContact() > VSMALL;
    }

    autoPtr<DynamicVectorList> nextToCheck(
        new DynamicVectorList);

    autoPtr<DynamicVectorList> auxToCheck(
        new DynamicVectorList);

    nextToCheck->append(bbMatrix_.getSVIndexForPoint_Wall(vMeshWallInfo_.getStartingPoint()));
    nextToCheck->append(bbMatrix_.faceNeighbourSubVolumes(nextToCheck()[0]));
    // InfoH << DEM_Info << " -- VM firstSV : " << nextToCheck()[0] << " point " << bbMatrix_[nextToCheck()[0]].center << endl;
    while (nextToCheck->size() > 0)
    {
        auxToCheck->clear();
        forAll (nextToCheck(),sV)
        {
            subVolumeProperties& cSubVolume = bbMatrix_[nextToCheck()[sV]];
            if (!cSubVolume.toCheck)
            {
                continue;
            }
            checkSubVolume(cSubVolume, nextToCheck()[sV]);

            if (cSubVolume.isCBody)
            {
                vMeshWallInfo_.startingPoint = cSubVolume.center;
                resetSubVolume(cSubVolume);

                return true;

            }
            auxToCheck().append(bbMatrix_.faceNeighbourSubVolumes(nextToCheck()[sV]));
        }
        const autoPtr<DynamicVectorList> helpPtr(nextToCheck.ptr());
        nextToCheck.set(auxToCheck.ptr());
        auxToCheck = helpPtr;
    }
    return false;
}
//---------------------------------------------------------------------------//
scalar virtualMeshWall::evaluateContact()
{
    contactCenter_ = vector::zero;

    if (vMeshWallInfo_.usesFittedBBoxGrid())
    {
        scalar contactVolume = 0;
        const vector& matrixSize =
            vMeshWallInfo_.subVolumeNVector;

        // A fitted finite VM represents the whole clipped box.  Iterate its
        // guaranteed-valid integer range directly: no seed, neighbour queue,
        // duplicate entries, or out-of-range starting-point index is needed.
        for (label i = 0; i < label(matrixSize[0]); i++)
        {
            for (label j = 0; j < label(matrixSize[1]); j++)
            {
                for (label k = 0; k < label(matrixSize[2]); k++)
                {
                    const vector subVolumeIndex(i, j, k);

                    subVolumeProperties& cSubVolume =
                        bbMatrix_[subVolumeIndex];

                    checkSubVolume(cSubVolume, subVolumeIndex);

                    if (cSubVolume.isCBody)
                    {
                        contactVolume += cSubVolume.cBodyVolume;
                        contactCenter_ +=
                            cSubVolume.cBodyCenter
                           *cSubVolume.cBodyVolume;
                    }
                }
            }
        }

        if (contactVolume > VSMALL)
        {
            contactCenter_ /= contactVolume;
        }

        return contactVolume;
    }

    // Original infinite-wall evaluation path.
    label volumeCount = 0;
    autoPtr<DynamicVectorList> nextToCheck(
        new DynamicVectorList);
    autoPtr<DynamicVectorList> auxToCheck(
        new DynamicVectorList);
    nextToCheck->append(bbMatrix_.getSVIndexForPoint_Wall(vMeshWallInfo_.getStartingPoint()));
    label iterCount(0);
    while (nextToCheck().size() > 0)
    {
        auxToCheck().clear();

        forAll (nextToCheck(),sV)
        {
            subVolumeProperties& cSubVolume = bbMatrix_[nextToCheck()[sV]];
            iterCount++;
            if (!cSubVolume.toCheck)
            {
                continue;
            }

            checkSubVolume(cSubVolume, nextToCheck()[sV]);
            if (cSubVolume.isCBody)
            {
                volumeCount++;
                contactCenter_ += cSubVolume.center;
                auxToCheck->append(bbMatrix_.faceNeighbourSubVolumes(nextToCheck()[sV]));
            }
        }
        const autoPtr<DynamicVectorList> helpPtr(nextToCheck.ptr());
        nextToCheck.set(auxToCheck.ptr());
        auxToCheck = helpPtr;
    }
    if (volumeCount > 0)
    {
        contactCenter_ /= volumeCount;
    }

    return volumeCount*bbMatrix_.getSubVolumeV();
}
//---------------------------------------------------------------------------//
void virtualMeshWall::checkSubVolume
(
    subVolumeProperties& subVolumePropertiesI,
    const vector& subVolumeIndex
)
{
    if (!subVolumePropertiesI.toCheck)
    {
        return;
    }

    if (!vMeshWallInfo_.usesFittedBBoxGrid())
    {
        // Preserve the original infinite-wall centre-sampling behavior.
        subVolumePropertiesI.isCBody =
            cGeomModel_.pointInside(subVolumePropertiesI.center);
        subVolumePropertiesI.toCheck = false;
        return;
    }

    const boundBox cellBB = bbMatrix_.getSubVolumeBB(subVolumeIndex);
    Foam::subVolume particleCell(cellBB);
    const volumeType cellType =
        cGeomModel_.getVolumeType(particleCell, true);

    scalar occupiedVolume = 0;
    point occupiedCenter = cellBB.midpoint();

    if (cellType == volumeType::inside)
    {
        occupiedVolume = cellBB.volume();
    }
    else if (cellType == volumeType::mixed)
    {
        boundBox limitedBB(cellBB);

        if
        (
            cGeomModel_.limitFinalSubVolume
            (
                particleCell,
                true,
                limitedBB
            )
        )
        {
            vector limitedMin = limitedBB.min();
            vector limitedMax = limitedBB.max();
            bool validLimitedBB = true;

            for (label dir = 0; dir < 3; dir++)
            {
                limitedMin[dir] = max
                (
                    limitedMin[dir],
                    cellBB.min()[dir]
                );
                limitedMax[dir] = min
                (
                    limitedMax[dir],
                    cellBB.max()[dir]
                );

                if ((limitedMax[dir] - limitedMin[dir]) <= VSMALL)
                {
                    validLimitedBB = false;
                }
            }

            if (validLimitedBB)
            {
                limitedBB = boundBox(limitedMin, limitedMax);
                occupiedVolume = limitedBB.volume();
                occupiedCenter = limitedBB.midpoint();
            }
        }
    }
    else if (cellType == volumeType::unknown)
    {
        // Non-STL geometry fallback.  Arbitrary STL particles use the
        // inside/mixed/outside path above.
        if (cGeomModel_.pointInside(subVolumePropertiesI.center))
        {
            occupiedVolume = cellBB.volume();
        }
    }

    subVolumePropertiesI.cBodyVolume = occupiedVolume;
    subVolumePropertiesI.cBodyCenter = occupiedCenter;
    subVolumePropertiesI.isCBody = occupiedVolume > VSMALL;
    subVolumePropertiesI.toCheck = false;
}
//---------------------------------------------------------------------------//
void virtualMeshWall::resetSubVolume(subVolumeProperties& subVolume)
{
    subVolume.toCheck = true;
    subVolume.isCBody = false;
    subVolume.isTBody = false;
    subVolume.isOnEdge = false;
    subVolume.cBodyVolume = 0;
    subVolume.cBodyCenter = subVolume.center;
}
//---------------------------------------------------------------------------//
label virtualMeshWall::getInternalSV()
{
    autoPtr<DynamicVectorList> nextToCheck(
        new DynamicVectorList);
    autoPtr<DynamicVectorList> auxToCheck(
        new DynamicVectorList);

    nextToCheck->append(bbMatrix_.getSVIndexForPoint(vMeshWallInfo_.getStartingPoint()));
    label innerSVCount(0);
    vectorHashSet octreeSvSet;
    
    while (nextToCheck->size() > 0)
    {
        auxToCheck().clear();
        forAll (nextToCheck(),sV)
        {   
            if (!octreeSvSet.found(nextToCheck()[sV]))
            {
                octreeSvSet.insert(nextToCheck()[sV]);
                subVolumeProperties& cSubVolume = bbMatrix_[nextToCheck()[sV]];
                if (cSubVolume.isCBody)
                {
                    bool isNotOnEdge(true);
                    List<vector> neighbourSubVolumes = bbMatrix_.faceNeighbourSubVolumes(nextToCheck()[sV]);
                    neighbourSubVolumes.append(bbMatrix_.edgeNeighbourSubVolumes(nextToCheck()[sV]));
                    neighbourSubVolumes.append(bbMatrix_.cornerNeighbourSubVolumes(nextToCheck()[sV]));
                    forAll(neighbourSubVolumes,nSV)
                    {
                        isNotOnEdge *= bbMatrix_[neighbourSubVolumes[nSV]].isCBody;
                    }
                    if(!isNotOnEdge)
                    {
                        cSubVolume.isOnEdge = true;
                    }
                    auxToCheck().append(bbMatrix_.faceNeighbourSubVolumes(nextToCheck()[sV]));
                }
                if(!cSubVolume.isOnEdge)
                {
                    innerSVCount++;
                }
            }
        }

        const autoPtr<DynamicVectorList> helpPtr(nextToCheck.ptr());
        nextToCheck.set(auxToCheck.ptr());
        auxToCheck = helpPtr;
    }
    return innerSVCount;
}
// ************************************************************************* //
