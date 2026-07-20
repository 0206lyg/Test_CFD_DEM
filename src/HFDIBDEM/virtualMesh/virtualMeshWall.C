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
#include "wallPlaneInfo.H"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

using namespace Foam;

namespace
{

typedef std::vector<vector> localPolygonVMW;

pointField createJointQuadratureVMW()
{
    const label sampleCount = 64;
    const label bitCount = 32;

    std::uint32_t directions[3][33] = {};

    // Sobol direction numbers for the first three dimensions.
    for (label bit = 1; bit <= bitCount; bit++)
    {
        directions[0][bit] = std::uint32_t(1) << (32 - bit);
    }

    directions[1][1] = std::uint32_t(1) << 31;

    for (label bit = 2; bit <= bitCount; bit++)
    {
        directions[1][bit] =
            directions[1][bit - 1]
          ^ (directions[1][bit - 1] >> 1);
    }

    directions[2][1] = std::uint32_t(1) << 31;
    directions[2][2] = std::uint32_t(3) << 30;

    for (label bit = 3; bit <= bitCount; bit++)
    {
        directions[2][bit] =
            directions[2][bit - 2]
          ^ (directions[2][bit - 2] >> 2)
          ^ directions[2][bit - 1];
    }

    // Fixed digital shifts move the points away from coordinate and diagonal
    // symmetry planes.  The sequence remains deterministic and reproducible.
    const std::uint32_t digitalShift[3] =
    {
        0x12345678u,
        0x87654321u,
        0xdeadbeefu
    };

    pointField quadraturePoints(sampleCount);
    const scalar integerScale = 1.0/4294967296.0;

    for (label sampleI = 0; sampleI < sampleCount; sampleI++)
    {
        const std::uint32_t grayCode =
            std::uint32_t(sampleI ^ (sampleI >> 1));

        point samplePoint(vector::zero);

        for (label dir = 0; dir < 3; dir++)
        {
            std::uint32_t coordinate = digitalShift[dir];

            for (label bit = 1; bit <= bitCount; bit++)
            {
                if (grayCode & (std::uint32_t(1) << (bit - 1)))
                {
                    coordinate ^= directions[dir][bit];
                }
            }

            samplePoint[dir] =
                (scalar(coordinate) + 0.5)*integerScale;
        }

        quadraturePoints[sampleI] = samplePoint;
    }

    return quadraturePoints;
}

const pointField& jointQuadratureVMW()
{
    static const pointField quadraturePoints = createJointQuadratureVMW();
    return quadraturePoints;
}

scalar cross2DVMW
(
    const vector& a,
    const vector& b,
    const vector& c
)
{
    return
        (b[0] - a[0])*(c[1] - a[1])
      - (b[1] - a[1])*(c[0] - a[0]);
}

void appendUniquePointVMW
(
    std::vector<point>& points,
    const point& candidate,
    const scalar tolerance
)
{
    for (std::vector<point>::const_iterator iter = points.begin();
         iter != points.end(); ++iter)
    {
        if (mag(*iter - candidate) <= tolerance)
        {
            return;
        }
    }

    points.push_back(candidate);
}

struct angleLessVMW
{
    vector center_;

    explicit angleLessVMW(const vector& center)
    :
    center_(center)
    {}

    bool operator()(const vector& a, const vector& b) const
    {
        return
            std::atan2(a[1] - center_[1], a[0] - center_[0])
          < std::atan2(b[1] - center_[1], b[0] - center_[0]);
    }
};

localPolygonVMW planeCellPolygonVMW
(
    const boundBox& cellBB,
    const finiteWallGeometry& wallGeometry
)
{
    const pointField boxPoints = cellBB.points();
    const scalar tolerance = wallGeometry.tolerance();
    std::vector<point> planePoints;

    // Derive the twelve AABB edges from coordinate differences instead of
    // relying on a particular boundBox::points() numbering convention.
    for (label pointI = 0; pointI < boxPoints.size(); pointI++)
    {
        for (label pointJ = pointI + 1; pointJ < boxPoints.size(); pointJ++)
        {
            label differingDirections = 0;

            for (label dir = 0; dir < 3; dir++)
            {
                if (mag(boxPoints[pointI][dir] - boxPoints[pointJ][dir]) > tolerance)
                {
                    differingDirections++;
                }
            }

            if (differingDirections != 1)
            {
                continue;
            }

            const scalar distanceI =
                wallGeometry.signedDistance(boxPoints[pointI]);

            const scalar distanceJ =
                wallGeometry.signedDistance(boxPoints[pointJ]);

            if (mag(distanceI) <= tolerance)
            {
                appendUniquePointVMW
                (
                    planePoints,
                    boxPoints[pointI],
                    tolerance
                );
            }

            if (mag(distanceJ) <= tolerance)
            {
                appendUniquePointVMW
                (
                    planePoints,
                    boxPoints[pointJ],
                    tolerance
                );
            }

            if
            (
                (distanceI < -tolerance && distanceJ > tolerance)
             || (distanceJ < -tolerance && distanceI > tolerance)
            )
            {
                const scalar fraction =
                    distanceI/(distanceI - distanceJ);

                appendUniquePointVMW
                (
                    planePoints,
                    boxPoints[pointI]
                  + fraction*(boxPoints[pointJ] - boxPoints[pointI]),
                    tolerance
                );
            }
        }
    }

    localPolygonVMW localPolygon;
    localPolygon.reserve(planePoints.size());

    for (std::vector<point>::const_iterator iter = planePoints.begin();
         iter != planePoints.end(); ++iter)
    {
        localPolygon.push_back(wallGeometry.localCoordinates(*iter));
    }

    if (localPolygon.size() >= 3)
    {
        vector center(vector::zero);

        for (localPolygonVMW::const_iterator iter = localPolygon.begin();
             iter != localPolygon.end(); ++iter)
        {
            center += *iter;
        }

        center /= scalar(localPolygon.size());
        std::sort(localPolygon.begin(), localPolygon.end(), angleLessVMW(center));
    }

    return localPolygon;
}

localPolygonVMW clipToFinitePatchVMW
(
    const localPolygonVMW& inputPolygon,
    const finiteWallGeometry& wallGeometry
)
{
    localPolygonVMW clippedPolygon(inputPolygon);
    const List<vector>& patchVertices = wallGeometry.localVertices();
    const scalar tolerance = wallGeometry.tolerance();

    forAll(patchVertices, edgeI)
    {
        if (clippedPolygon.empty())
        {
            break;
        }

        const vector& edgeStart = patchVertices[edgeI];
        const vector& edgeEnd =
            patchVertices[(edgeI + 1) % patchVertices.size()];

        const scalar edgeLength =
            Foam::sqrt
            (
                sqr(edgeEnd[0] - edgeStart[0])
              + sqr(edgeEnd[1] - edgeStart[1])
            );

        const scalar edgeTolerance = tolerance*edgeLength;

        localPolygonVMW input(clippedPolygon);
        clippedPolygon.clear();

        vector previous = input.back();
        scalar previousValue =
            cross2DVMW(edgeStart, edgeEnd, previous);
        bool previousInside = previousValue >= -edgeTolerance;

        for (localPolygonVMW::const_iterator iter = input.begin();
             iter != input.end(); ++iter)
        {
            const vector current = *iter;
            const scalar currentValue =
                cross2DVMW(edgeStart, edgeEnd, current);
            const bool currentInside = currentValue >= -edgeTolerance;

            if (currentInside != previousInside)
            {
                const scalar denominator =
                    previousValue - currentValue;

                if (mag(denominator) > VSMALL)
                {
                    const scalar fraction = previousValue/denominator;
                    clippedPolygon.push_back
                    (
                        previous + fraction*(current - previous)
                    );
                }
            }

            if (currentInside)
            {
                clippedPolygon.push_back(current);
            }

            previous = current;
            previousValue = currentValue;
            previousInside = currentInside;
        }
    }

    return clippedPolygon;
}

scalar polygonAreaAndCenterVMW
(
    const localPolygonVMW& polygon,
    vector& center
)
{
    center = vector::zero;

    if (polygon.size() < 3)
    {
        return 0;
    }

    scalar twiceArea = 0;

    for (label pointI = 0; pointI < label(polygon.size()); pointI++)
    {
        const vector& a = polygon[pointI];
        const vector& b = polygon[(pointI + 1) % polygon.size()];
        const scalar cross = a[0]*b[1] - b[0]*a[1];

        twiceArea += cross;
        center[0] += (a[0] + b[0])*cross;
        center[1] += (a[1] + b[1])*cross;
    }

    if (mag(twiceArea) <= VSMALL)
    {
        center = vector::zero;
        return 0;
    }

    center /= 3*twiceArea;

    return 0.5*mag(twiceArea);
}

} // End anonymous namespace

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
    vMeshWallInfo.fitToBBox),
finiteWallGeometry_()
{
    if
    (
        vMeshWallInfo_.hasFiniteWallGeometry()
     && !wallPlaneInfo::usesLegacyAxisAlignedFinitePath
        (
            vMeshWallInfo_.wallName
        )
    )
    {
        finiteWallGeometry_ =
            finiteWallGeometry::lookup(vMeshWallInfo_.wallName);
    }
}

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
scalar virtualMeshWall::evaluateContactArea()
{
    contactCenter_ = vector::zero;

    if
    (
        !vMeshWallInfo_.usesFittedBBoxGrid()
     || !finiteWallGeometry_
    )
    {
        return 0;
    }

    scalar contactArea = 0;
    const vector& matrixSize = vMeshWallInfo_.subVolumeNVector;
    const vector cellSize = bbMatrix_.getSubVolumeSize();
    const boundBox matrixBB = bbMatrix_.getBBox();
    const vector& normal = finiteWallGeometry_->normal();

    const label cellCount[3] =
    {
        label(matrixSize[0]),
        label(matrixSize[1]),
        label(matrixSize[2])
    };

    if (cellCount[0] <= 0 || cellCount[1] <= 0 || cellCount[2] <= 0)
    {
        return 0;
    }

    // Sweep columns along the direction with the largest signed-distance
    // change across one cell.  For either of the two remaining cell widths,
    // the plane can then move by no more than one sweep-cell width.  Thus a
    // column contains only a bounded number of plane candidates and the
    // traversal is O(N^2), including for an inclined plane.
    label sweepDirection = 0;
    scalar sweepScale = mag(normal[0])*cellSize[0];

    for (label dir = 1; dir < 3; dir++)
    {
        const scalar directionScale = mag(normal[dir])*cellSize[dir];

        if (directionScale > sweepScale)
        {
            sweepDirection = dir;
            sweepScale = directionScale;
        }
    }

    if (sweepScale <= VSMALL)
    {
        return 0;
    }

    // An axis-aligned finite wall produces a fitted plane grid that is only
    // one cell thick in the sweep direction.  In that case the original
    // dense traversal is already O(N^2) and is cheaper than constructing the
    // sparse column intersections.  Keep a two-cell allowance for planes
    // lying on a grid boundary or for small round-off-induced padding.
    //
    // Both traversal paths call the same evaluateFiniteAreaCell() routine;
    // this branch changes candidate enumeration only, not contact physics.
    const label denseSweepCellLimit = 2;

    if (cellCount[sweepDirection] <= denseSweepCellLimit)
    {
        for (label i = 0; i < cellCount[0]; i++)
        {
            for (label j = 0; j < cellCount[1]; j++)
            {
                for (label k = 0; k < cellCount[2]; k++)
                {
                    const vector subVolumeIndex(i, j, k);
                    const boundBox cellBB =
                        bbMatrix_.getSubVolumeBB(subVolumeIndex);

                    point occupiedCenter(cellBB.midpoint());
                    const scalar occupiedArea =
                        evaluateFiniteAreaCell(cellBB, occupiedCenter);

                    if (occupiedArea > VSMALL)
                    {
                        contactArea += occupiedArea;
                        contactCenter_ += occupiedCenter*occupiedArea;
                    }
                }
            }
        }

        if (contactArea > VSMALL)
        {
            contactCenter_ /= contactArea;
        }

        return contactArea;
    }

    label columnDirections[2] = {-1, -1};
    label columnDirectionI = 0;

    for (label dir = 0; dir < 3; dir++)
    {
        if (dir != sweepDirection)
        {
            columnDirections[columnDirectionI++] = dir;
        }
    }

    const scalar planeConstant =
        normal & finiteWallGeometry_->planePoint();

    const scalar coordinateTolerance =
        finiteWallGeometry_->tolerance()
       /(mag(normal[sweepDirection]) + VSMALL);

    const scalar sweepOrigin = matrixBB.min()[sweepDirection];
    const scalar sweepCellSize = cellSize[sweepDirection];

    for
    (
        label column0 = 0;
        column0 < cellCount[columnDirections[0]];
        column0++
    )
    {
        const scalar column0Min =
            matrixBB.min()[columnDirections[0]]
          + column0*cellSize[columnDirections[0]];

        const scalar column0Max =
            column0Min + cellSize[columnDirections[0]];

        for
        (
            label column1 = 0;
            column1 < cellCount[columnDirections[1]];
            column1++
        )
        {
            const scalar column1Min =
                matrixBB.min()[columnDirections[1]]
              + column1*cellSize[columnDirections[1]];

            const scalar column1Max =
                column1Min + cellSize[columnDirections[1]];

            scalar sweepCoordinateMin = GREAT;
            scalar sweepCoordinateMax = -GREAT;

            // A plane is affine, so the exact sweep-coordinate extrema over
            // this rectangular column occur at its four transverse corners.
            for (label corner0 = 0; corner0 < 2; corner0++)
            {
                const scalar coordinate0 =
                    corner0 == 0 ? column0Min : column0Max;

                for (label corner1 = 0; corner1 < 2; corner1++)
                {
                    const scalar coordinate1 =
                        corner1 == 0 ? column1Min : column1Max;

                    const scalar sweepCoordinate =
                    (
                        planeConstant
                      - normal[columnDirections[0]]*coordinate0
                      - normal[columnDirections[1]]*coordinate1
                    )/normal[sweepDirection];

                    sweepCoordinateMin =
                        min(sweepCoordinateMin, sweepCoordinate);

                    sweepCoordinateMax =
                        max(sweepCoordinateMax, sweepCoordinate);
                }
            }

            // A sweep cell [i, i+1] intersects the closed coordinate interval
            // [a, b] exactly when i >= ceil(a - 1) and i <= floor(b).
            // Expanding a and b by the existing geometry tolerance preserves
            // both cells at a face-coincident plane; the unchanged area
            // evaluator then applies its one-sided ownership rule.
            label firstSweepCell = label
            (
                std::ceil
                (
                    (sweepCoordinateMin - coordinateTolerance - sweepOrigin)
                   /sweepCellSize
                  - 1
                )
            );

            label lastSweepCell = label
            (
                std::floor
                (
                    (sweepCoordinateMax + coordinateTolerance - sweepOrigin)
                   /sweepCellSize
                )
            );

            if
            (
                lastSweepCell < 0
             || firstSweepCell >= cellCount[sweepDirection]
            )
            {
                continue;
            }

            firstSweepCell = max(firstSweepCell, label(0));
            lastSweepCell = min
            (
                lastSweepCell,
                cellCount[sweepDirection] - 1
            );

            for
            (
                label sweepCell = firstSweepCell;
                sweepCell <= lastSweepCell;
                sweepCell++
            )
            {
                vector subVolumeIndex(vector::zero);
                subVolumeIndex[columnDirections[0]] = column0;
                subVolumeIndex[columnDirections[1]] = column1;
                subVolumeIndex[sweepDirection] = sweepCell;

                const boundBox cellBB =
                    bbMatrix_.getSubVolumeBB(subVolumeIndex);

                point occupiedCenter(cellBB.midpoint());
                const scalar occupiedArea =
                    evaluateFiniteAreaCell(cellBB, occupiedCenter);

                if (occupiedArea > VSMALL)
                {
                    contactArea += occupiedArea;
                    contactCenter_ += occupiedCenter*occupiedArea;
                }
            }
        }
    }

    if (contactArea > VSMALL)
    {
        contactCenter_ /= contactArea;
    }

    return contactArea;
}

//---------------------------------------------------------------------------//
scalar virtualMeshWall::evaluateFiniteCell
(
    const boundBox& cellBB,
    point& occupiedCenter
)
{
    if (!finiteWallGeometry_)
    {
        return 0;
    }

    const finiteWallGeometry::cellStatus wallStatus =
        finiteWallGeometry_->classifyContactCell(cellBB);

    if (wallStatus == finiteWallGeometry::cellOutside)
    {
        return 0;
    }

    Foam::subVolume particleCell(cellBB);
    const volumeType particleType =
        cGeomModel_.getVolumeType(particleCell, true);

    if (particleType == volumeType::outside)
    {
        return 0;
    }

    if
    (
        wallStatus == finiteWallGeometry::cellInside
     && particleType == volumeType::inside
    )
    {
        occupiedCenter = cellBB.midpoint();
        return cellBB.volume();
    }

    const pointField& normalizedSamplePoints = jointQuadratureVMW();
    const label sampleCount = normalizedSamplePoints.size();

    pointField samplePoints(sampleCount);

    forAll(samplePoints, sampleI)
    {
        for (label dir = 0; dir < 3; dir++)
        {
            samplePoints[sampleI][dir] =
                cellBB.min()[dir]
              + normalizedSamplePoints[sampleI][dir]*cellBB.span()[dir];
        }
    }

    boolList particleInside(sampleCount, false);

    if (particleType == volumeType::inside)
    {
        forAll(particleInside, pointI)
        {
            particleInside[pointI] = true;
        }
    }
    else
    {
        const boolList sampledInside = cGeomModel_.pointInside(samplePoints);

        if (sampledInside.size() == sampleCount)
        {
            particleInside = sampledInside;
        }
        else
        {
            forAll(samplePoints, pointI)
            {
                particleInside[pointI] =
                    cGeomModel_.pointInside(samplePoints[pointI]);
            }
        }
    }

    label occupiedSamples = 0;
    occupiedCenter = vector::zero;

    forAll(samplePoints, pointI)
    {
        if
        (
            particleInside[pointI]
         && finiteWallGeometry_->pointInContactRegion(samplePoints[pointI])
        )
        {
            occupiedSamples++;
            occupiedCenter += samplePoints[pointI];
        }
    }

    if (occupiedSamples == 0)
    {
        occupiedCenter = cellBB.midpoint();
        return 0;
    }

    occupiedCenter /= scalar(occupiedSamples);

    return
        cellBB.volume()
       *scalar(occupiedSamples)/scalar(sampleCount);
}

//---------------------------------------------------------------------------//
scalar virtualMeshWall::evaluateFiniteAreaCell
(
    const boundBox& cellBB,
    point& occupiedCenter
)
{
    if (!finiteWallGeometry_)
    {
        return 0;
    }

    scalar minDistance = 0;
    scalar maxDistance = 0;
    finiteWallGeometry_->planeDistanceRange
    (
        cellBB,
        minDistance,
        maxDistance
    );

    if
    (
        minDistance > finiteWallGeometry_->tolerance()
     || maxDistance < -finiteWallGeometry_->tolerance()
    )
    {
        return 0;
    }

    // When the wall plane coincides with a Cartesian cell face, only the
    // cell on the penetration side owns that face.  This prevents an exact
    // factor-of-two area error for axis-aligned walls.
    if
    (
        maxDistance <= finiteWallGeometry_->tolerance()
     && minDistance < -finiteWallGeometry_->tolerance()
    )
    {
        return 0;
    }

    localPolygonVMW planePolygon =
        planeCellPolygonVMW(cellBB, *finiteWallGeometry_);

    planePolygon =
        clipToFinitePatchVMW(planePolygon, *finiteWallGeometry_);

    vector localCenter(vector::zero);
    const scalar polygonArea =
        polygonAreaAndCenterVMW(planePolygon, localCenter);

    if (polygonArea <= VSMALL)
    {
        return 0;
    }

    Foam::subVolume particleCell(cellBB);
    const volumeType particleType =
        cGeomModel_.getVolumeType(particleCell, true);

    if (particleType == volumeType::outside)
    {
        return 0;
    }

    if (particleType == volumeType::inside)
    {
        occupiedCenter =
            finiteWallGeometry_->globalCoordinates(localCenter);
        return polygonArea;
    }

    const label triangleCount = label(planePolygon.size());
    const label samplesPerTriangle = 16;
    const pointField& normalizedSamplePoints = jointQuadratureVMW();

    pointField samplePoints(samplesPerTriangle*triangleCount);
    scalarField sampleWeights(samplesPerTriangle*triangleCount, 0);
    label sampleI = 0;

    for (label pointI = 0; pointI < triangleCount; pointI++)
    {
        const vector& a = localCenter;
        const vector& b = planePolygon[pointI];
        const vector& c =
            planePolygon[(pointI + 1) % planePolygon.size()];

        const scalar triangleArea = 0.5*mag(cross2DVMW(a, b, c));

        if (triangleArea <= VSMALL)
        {
            continue;
        }

        // Map the first sixteen low-discrepancy points uniformly from the
        // unit square to this triangle.  This avoids locking the area samples
        // to an inclined particle edge or a Cartesian virtual-mesh diagonal.
        for
        (
            label triangleSampleI = 0;
            triangleSampleI < samplesPerTriangle;
            triangleSampleI++
        )
        {
            const scalar rootU = Foam::sqrt
            (
                normalizedSamplePoints[triangleSampleI][0]
            );

            const scalar v =
                normalizedSamplePoints[triangleSampleI][1];

            const vector localSample =
                (1 - rootU)*a
              + rootU*(1 - v)*b
              + rootU*v*c;

            samplePoints[sampleI] =
                finiteWallGeometry_->globalCoordinates(localSample);

            sampleWeights[sampleI++] =
                triangleArea/samplesPerTriangle;
        }
    }

    samplePoints.setSize(sampleI);
    sampleWeights.setSize(sampleI);

    boolList sampledInside = cGeomModel_.pointInside(samplePoints);

    if (sampledInside.size() != samplePoints.size())
    {
        sampledInside.setSize(samplePoints.size(), false);

        forAll(samplePoints, pointI)
        {
            sampledInside[pointI] =
                cGeomModel_.pointInside(samplePoints[pointI]);
        }
    }

    scalar occupiedArea = 0;
    occupiedCenter = vector::zero;

    forAll(samplePoints, pointI)
    {
        if (sampledInside[pointI])
        {
            occupiedArea += sampleWeights[pointI];
            occupiedCenter += samplePoints[pointI]*sampleWeights[pointI];
        }
    }

    if (occupiedArea <= VSMALL)
    {
        occupiedCenter = cellBB.midpoint();
        return 0;
    }

    occupiedCenter /= occupiedArea;
    return occupiedArea;
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

    if (finiteWallGeometry_)
    {
        point occupiedCenter(cellBB.midpoint());
        const scalar occupiedVolume =
            evaluateFiniteCell(cellBB, occupiedCenter);

        subVolumePropertiesI.cBodyVolume = occupiedVolume;
        subVolumePropertiesI.cBodyCenter = occupiedCenter;
        subVolumePropertiesI.isCBody = occupiedVolume > VSMALL;
        subVolumePropertiesI.toCheck = false;
        return;
    }

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
