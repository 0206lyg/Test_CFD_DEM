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
#include "stlBased.H"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <utility>
#include <vector>

using namespace Foam;

namespace
{

typedef std::vector<vector> localPolygonVMW;
typedef std::vector<point> polygon3DVMW;
typedef std::vector<polygon3DVMW> convexPolyhedronVMW;

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

void cleanPolygon3DVMW
(
    polygon3DVMW& polygon,
    const scalar tolerance
)
{
    polygon3DVMW cleanedPolygon;
    cleanedPolygon.reserve(polygon.size());

    for (polygon3DVMW::const_iterator iter = polygon.begin();
         iter != polygon.end(); ++iter)
    {
        if
        (
            cleanedPolygon.empty()
         || mag(*iter - cleanedPolygon.back()) > tolerance
        )
        {
            cleanedPolygon.push_back(*iter);
        }
    }

    if
    (
        cleanedPolygon.size() > 1
     && mag(cleanedPolygon.front() - cleanedPolygon.back()) <= tolerance
    )
    {
        cleanedPolygon.pop_back();
    }

    if (cleanedPolygon.size() < 3)
    {
        polygon.clear();
        return;
    }

    vector twiceArea(vector::zero);
    const point& origin = cleanedPolygon.front();

    for (label pointI = 1; pointI + 1 < label(cleanedPolygon.size()); pointI++)
    {
        twiceArea +=
            (cleanedPolygon[pointI] - origin)
          ^ (cleanedPolygon[pointI + 1] - origin);
    }

    if (mag(twiceArea) <= sqr(tolerance))
    {
        polygon.clear();
        return;
    }

    polygon.swap(cleanedPolygon);
}

struct planeAngleLessVMW
{
    point center_;
    vector tangent1_;
    vector tangent2_;

    planeAngleLessVMW
    (
        const point& center,
        const vector& tangent1,
        const vector& tangent2
    )
    :
    center_(center),
    tangent1_(tangent1),
    tangent2_(tangent2)
    {}

    bool operator()(const point& a, const point& b) const
    {
        const vector relativeA = a - center_;
        const vector relativeB = b - center_;

        return
            std::atan2(relativeA & tangent2_, relativeA & tangent1_)
          < std::atan2(relativeB & tangent2_, relativeB & tangent1_);
    }
};

convexPolyhedronVMW cellPolyhedronVMW(const boundBox& cellBB)
{
    const scalar x0 = cellBB.min()[0];
    const scalar y0 = cellBB.min()[1];
    const scalar z0 = cellBB.min()[2];
    const scalar x1 = cellBB.max()[0];
    const scalar y1 = cellBB.max()[1];
    const scalar z1 = cellBB.max()[2];

    const point p000(x0, y0, z0);
    const point p001(x0, y0, z1);
    const point p010(x0, y1, z0);
    const point p011(x0, y1, z1);
    const point p100(x1, y0, z0);
    const point p101(x1, y0, z1);
    const point p110(x1, y1, z0);
    const point p111(x1, y1, z1);

    convexPolyhedronVMW polyhedron;
    polyhedron.reserve(6);
    polyhedron.push_back(polygon3DVMW{p000, p001, p011, p010});
    polyhedron.push_back(polygon3DVMW{p100, p110, p111, p101});
    polyhedron.push_back(polygon3DVMW{p000, p100, p101, p001});
    polyhedron.push_back(polygon3DVMW{p010, p011, p111, p110});
    polyhedron.push_back(polygon3DVMW{p000, p010, p110, p100});
    polyhedron.push_back(polygon3DVMW{p001, p101, p111, p011});

    return polyhedron;
}

bool clipPolyhedronVMW
(
    convexPolyhedronVMW& polyhedron,
    const vector& inwardNormal,
    const scalar offset,
    const scalar tolerance
)
{
    convexPolyhedronVMW clippedPolyhedron;
    clippedPolyhedron.reserve(polyhedron.size() + 1);
    std::vector<point> capPoints;

    for (convexPolyhedronVMW::const_iterator faceIter = polyhedron.begin();
         faceIter != polyhedron.end(); ++faceIter)
    {
        if (faceIter->empty())
        {
            continue;
        }

        polygon3DVMW clippedFace;
        clippedFace.reserve(faceIter->size() + 1);

        point previous = faceIter->back();
        scalar previousValue = (inwardNormal & previous) + offset;
        bool previousInside = previousValue >= -tolerance;

        for (polygon3DVMW::const_iterator pointIter = faceIter->begin();
             pointIter != faceIter->end(); ++pointIter)
        {
            const point current = *pointIter;
            const scalar currentValue = (inwardNormal & current) + offset;
            const bool currentInside = currentValue >= -tolerance;

            if (currentInside != previousInside)
            {
                const scalar denominator = previousValue - currentValue;

                if (mag(denominator) > VSMALL)
                {
                    scalar fraction = previousValue/denominator;
                    fraction = max(scalar(0), min(scalar(1), fraction));
                    const point intersection =
                        previous + fraction*(current - previous);

                    clippedFace.push_back(intersection);
                    appendUniquePointVMW
                    (
                        capPoints,
                        intersection,
                        tolerance
                    );
                }
            }

            if (currentInside)
            {
                clippedFace.push_back(current);
            }

            previous = current;
            previousValue = currentValue;
            previousInside = currentInside;
        }

        cleanPolygon3DVMW(clippedFace, tolerance);

        if (!clippedFace.empty())
        {
            clippedPolyhedron.push_back(clippedFace);
        }
    }

    if (clippedPolyhedron.empty())
    {
        polyhedron.clear();
        return false;
    }

    if (capPoints.size() >= 3)
    {
        point capCenter(vector::zero);

        for (std::vector<point>::const_iterator iter = capPoints.begin();
             iter != capPoints.end(); ++iter)
        {
            capCenter += *iter;
        }

        capCenter /= scalar(capPoints.size());

        vector referenceDirection(1, 0, 0);

        if (mag(inwardNormal[0]) > 0.8)
        {
            referenceDirection = vector(0, 1, 0);
        }

        vector tangent1 = inwardNormal ^ referenceDirection;

        if (mag(tangent1) <= VSMALL)
        {
            referenceDirection = vector(0, 0, 1);
            tangent1 = inwardNormal ^ referenceDirection;
        }

        tangent1 /= mag(tangent1);
        vector tangent2 = inwardNormal ^ tangent1;
        tangent2 /= mag(tangent2);

        std::sort
        (
            capPoints.begin(),
            capPoints.end(),
            planeAngleLessVMW(capCenter, tangent1, tangent2)
        );

        polygon3DVMW capPolygon(capPoints.begin(), capPoints.end());
        cleanPolygon3DVMW(capPolygon, tolerance);

        if (!capPolygon.empty())
        {
            clippedPolyhedron.push_back(capPolygon);
        }
    }

    polyhedron.swap(clippedPolyhedron);
    return true;
}

scalar polyhedronVolumeAndCenterVMW
(
    const convexPolyhedronVMW& polyhedron,
    point& center
)
{
    center = vector::zero;
    point interiorPoint(vector::zero);
    label pointCount = 0;

    for (convexPolyhedronVMW::const_iterator faceIter = polyhedron.begin();
         faceIter != polyhedron.end(); ++faceIter)
    {
        for (polygon3DVMW::const_iterator pointIter = faceIter->begin();
             pointIter != faceIter->end(); ++pointIter)
        {
            interiorPoint += *pointIter;
            pointCount++;
        }
    }

    if (pointCount == 0)
    {
        return 0;
    }

    interiorPoint /= scalar(pointCount);
    scalar volume = 0;

    for (convexPolyhedronVMW::const_iterator faceIter = polyhedron.begin();
         faceIter != polyhedron.end(); ++faceIter)
    {
        if (faceIter->size() < 3)
        {
            continue;
        }

        const point& a = faceIter->front();

        for (label pointI = 1; pointI + 1 < label(faceIter->size()); pointI++)
        {
            const point& b = (*faceIter)[pointI];
            const point& c = (*faceIter)[pointI + 1];
            const scalar tetrahedronVolume = mag
            (
                (a - interiorPoint)
              & ((b - interiorPoint) ^ (c - interiorPoint))
            )/6.0;

            if (tetrahedronVolume <= VSMALL)
            {
                continue;
            }

            volume += tetrahedronVolume;
            center +=
                0.25*(interiorPoint + a + b + c)*tetrahedronVolume;
        }
    }

    if (volume > VSMALL)
    {
        center /= volume;
    }

    return volume;
}

scalar clippedConvexVolumeVMW
(
    const boundBox& cellBB,
    const finiteWallGeometry& wallGeometry,
    const DynamicList<vector>& particleNormals,
    const DynamicList<scalar>& particleOffsets,
    const scalar tolerance,
    point& occupiedCenter
)
{
    convexPolyhedronVMW polyhedron = cellPolyhedronVMW(cellBB);

    if
    (
        !clipPolyhedronVMW
        (
            polyhedron,
            wallGeometry.normal(),
            -(wallGeometry.normal() & wallGeometry.planePoint()),
            tolerance
        )
    )
    {
        return 0;
    }

    const List<vector>& edgeNormals = wallGeometry.globalEdgeNormals();
    const List<scalar>& edgeOffsets = wallGeometry.edgeOffsets();

    forAll(edgeNormals, edgeI)
    {
        if
        (
            !clipPolyhedronVMW
            (
                polyhedron,
                edgeNormals[edgeI],
                edgeOffsets[edgeI]
              - (edgeNormals[edgeI] & wallGeometry.planePoint()),
                tolerance
            )
        )
        {
            return 0;
        }
    }

    forAll(particleNormals, planeI)
    {
        if
        (
            !clipPolyhedronVMW
            (
                polyhedron,
                particleNormals[planeI],
                particleOffsets[planeI],
                tolerance
            )
        )
        {
            return 0;
        }
    }

    return polyhedronVolumeAndCenterVMW(polyhedron, occupiedCenter);
}

localPolygonVMW clipToConvexParticleVMW
(
    const localPolygonVMW& inputPolygon,
    const finiteWallGeometry& wallGeometry,
    const DynamicList<vector>& particleNormals,
    const DynamicList<scalar>& particleOffsets,
    const scalar tolerance
)
{
    localPolygonVMW clippedPolygon(inputPolygon);

    forAll(particleNormals, planeI)
    {
        if (clippedPolygon.empty())
        {
            break;
        }

        localPolygonVMW input(clippedPolygon);
        clippedPolygon.clear();

        vector previous = input.back();
        scalar previousValue =
            (particleNormals[planeI]
           & wallGeometry.globalCoordinates(previous))
          + particleOffsets[planeI];
        bool previousInside = previousValue >= -tolerance;

        for (localPolygonVMW::const_iterator iter = input.begin();
             iter != input.end(); ++iter)
        {
            const vector current = *iter;
            const scalar currentValue =
                (particleNormals[planeI]
               & wallGeometry.globalCoordinates(current))
              + particleOffsets[planeI];
            const bool currentInside = currentValue >= -tolerance;

            if (currentInside != previousInside)
            {
                const scalar denominator = previousValue - currentValue;

                if (mag(denominator) > VSMALL)
                {
                    scalar fraction = previousValue/denominator;
                    fraction = max(scalar(0), min(scalar(1), fraction));
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
finiteWallGeometry_(),
exactConvexClipping_(false),
convexParticleNormals_(),
convexParticleOffsets_(),
particleClipTolerance_(1e-12)
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

    if
    (
        finiteWallGeometry_
     && cGeomModel_.getcType() == convex
    )
    {
        initialiseConvexParticleClipping();
    }
}

virtualMeshWall::~virtualMeshWall()
{
}
//---------------------------------------------------------------------------//
void virtualMeshWall::initialiseConvexParticleClipping()
{
    const stlBased* stlGeometry =
        dynamic_cast<const stlBased*>(&cGeomModel_);

    if (!stlGeometry)
    {
        FatalErrorInFunction
            << "Geometry declared as convex for inclined finite-wall "
            << "contact, but it does not provide an STL surface."
            << exit(FatalError);
    }

    const triSurface& surface = stlGeometry->triSurfaceGeometry();
    const pointField& points = surface.points();

    if (surface.empty() || points.size() < 4)
    {
        FatalErrorInFunction
            << "Convex particle used for inclined finite-wall contact must "
            << "contain at least four points and four triangular faces."
            << exit(FatalError);
    }

    typedef std::pair<label, label> edgeKey;
    std::map<edgeKey, label> edgeUseCount;
    std::vector<bool> referencedPoint(points.size(), false);

    forAll(surface, faceI)
    {
        const typename triSurface::FaceType& face = surface[faceI];

        if (face.size() != 3)
        {
            FatalErrorInFunction
                << "Convex particle face " << faceI
                << " is not triangular."
                << exit(FatalError);
        }

        forAll(face, pointI)
        {
            const label vertexI = face[pointI];

            if (vertexI < 0 || vertexI >= points.size())
            {
                FatalErrorInFunction
                    << "Convex particle face " << faceI
                    << " contains invalid point index " << vertexI << "."
                    << exit(FatalError);
            }

            referencedPoint[vertexI] = true;

            const label nextVertexI = face[(pointI + 1) % face.size()];
            const edgeKey edge
            (
                min(vertexI, nextVertexI),
                max(vertexI, nextVertexI)
            );

            edgeUseCount[edge]++;
        }
    }

    for (std::map<edgeKey, label>::const_iterator iter = edgeUseCount.begin();
         iter != edgeUseCount.end(); ++iter)
    {
        if (iter->first.first == iter->first.second || iter->second != 2)
        {
            FatalErrorInFunction
                << "Convex particle STL is not a closed two-manifold: edge ("
                << iter->first.first << ' ' << iter->first.second
                << ") is used " << iter->second << " times."
                << exit(FatalError);
        }
    }

    point interiorPoint(vector::zero);
    label referencedCount = 0;

    forAll(points, pointI)
    {
        if (referencedPoint[pointI])
        {
            interiorPoint += points[pointI];
            referencedCount++;
        }
    }

    if (referencedCount < 4)
    {
        FatalErrorInFunction
            << "Convex particle STL references fewer than four points."
            << exit(FatalError);
    }

    interiorPoint /= scalar(referencedCount);

    const boundBox particleBounds(points, false);
    particleClipTolerance_ = max
    (
        finiteWallGeometry_->tolerance(),
        scalar(1e-10)*max(mag(particleBounds.span()), scalar(1e-3))
    );

    forAll(surface, faceI)
    {
        const typename triSurface::FaceType& face = surface[faceI];
        const point& a = points[face[0]];
        const point& b = points[face[1]];
        const point& c = points[face[2]];
        vector inwardNormal = (b - a) ^ (c - a);
        const scalar twiceFaceArea = mag(inwardNormal);

        if (twiceFaceArea <= sqr(particleClipTolerance_))
        {
            FatalErrorInFunction
                << "Convex particle face " << faceI
                << " is geometrically degenerate."
                << exit(FatalError);
        }

        const scalar centerValue = inwardNormal & (interiorPoint - a);

        if
        (
            mag(centerValue)
         <= particleClipTolerance_*twiceFaceArea
        )
        {
            FatalErrorInFunction
                << "Convex particle has zero three-dimensional extent at "
                << "face " << faceI << "."
                << exit(FatalError);
        }

        if (centerValue < 0)
        {
            inwardNormal = -inwardNormal;
        }

        inwardNormal /= twiceFaceArea;
        const scalar offset = -(inwardNormal & a);

        forAll(points, pointI)
        {
            if
            (
                referencedPoint[pointI]
             && (inwardNormal & points[pointI]) + offset
                < -particleClipTolerance_
            )
            {
                FatalErrorInFunction
                    << "Geometry declared convex is not globally convex: "
                    << "point " << pointI << " lies outside face "
                    << faceI << " by "
                    << -((inwardNormal & points[pointI]) + offset) << "."
                    << exit(FatalError);
            }
        }

        bool duplicatePlane = false;

        forAll(convexParticleNormals_, planeI)
        {
            if
            (
                mag(inwardNormal - convexParticleNormals_[planeI]) <= 1e-9
             && mag(offset - convexParticleOffsets_[planeI])
                <= 10*particleClipTolerance_
            )
            {
                duplicatePlane = true;
                break;
            }
        }

        if (!duplicatePlane)
        {
            convexParticleNormals_.append(inwardNormal);
            convexParticleOffsets_.append(offset);
        }
    }

    if (convexParticleNormals_.size() < 4)
    {
        FatalErrorInFunction
            << "Convex particle provides fewer than four distinct support "
            << "planes."
            << exit(FatalError);
    }

    exactConvexClipping_ = true;
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

    if (exactConvexClipping_)
    {
        const scalar occupiedVolume = clippedConvexVolumeVMW
        (
            cellBB,
            *finiteWallGeometry_,
            convexParticleNormals_,
            convexParticleOffsets_,
            particleClipTolerance_,
            occupiedCenter
        );

        if (occupiedVolume <= VSMALL)
        {
            occupiedCenter = cellBB.midpoint();
            return 0;
        }

        return occupiedVolume;
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

    if (exactConvexClipping_)
    {
        planePolygon = clipToConvexParticleVMW
        (
            planePolygon,
            *finiteWallGeometry_,
            convexParticleNormals_,
            convexParticleOffsets_,
            particleClipTolerance_
        );

        const scalar occupiedArea =
            polygonAreaAndCenterVMW(planePolygon, localCenter);

        if (occupiedArea <= VSMALL)
        {
            occupiedCenter = cellBB.midpoint();
            return 0;
        }

        occupiedCenter =
            finiteWallGeometry_->globalCoordinates(localCenter);
        return occupiedArea;
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
