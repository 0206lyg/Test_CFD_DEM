/*---------------------------------------------------------------------------*\
  Geometry helper for finite planar DEM walls with arbitrary orientation.
\*---------------------------------------------------------------------------*/

#include "finiteWallGeometry.H"

#include "wallPlaneInfo.H"
#include "error.H"

#include <cmath>

using namespace Foam;

namespace
{

scalar maxScalarFWG(const scalar a, const scalar b)
{
    return a > b ? a : b;
}

}

//---------------------------------------------------------------------------//
finiteWallGeometry::finiteWallGeometry(const string& wallName)
:
wallName_(wallName),
planePoint_(vector::zero),
normal_(vector::zero),
tangent1_(vector::zero),
tangent2_(vector::zero),
vertices_(),
localVertices_(),
patchBounds_(boundBox()),
tolerance_(1e-12),
revision_(0)
{
    const HashTable<List<vector>,string,Hash<string>>& planeInfo =
        wallPlaneInfo::getWallPlaneInfo();

    if (!planeInfo.found(wallName_))
    {
        FatalErrorInFunction
            << "No plane geometry is registered for finite wall '"
            << wallName_ << "'."
            << exit(FatalError);
    }

    const HashTable<label,string,Hash<string>>& revisionInfo =
        wallPlaneInfo::getWallGeometryRevisionInfo();

    if (revisionInfo.found(wallName_))
    {
        revision_ = revisionInfo[wallName_];
    }

    normal_ = planeInfo[wallName_][0];
    planePoint_ = planeInfo[wallName_][1];

    if (mag(normal_) <= VSMALL)
    {
        FatalErrorInFunction
            << "Finite wall '" << wallName_ << "' has a zero normal."
            << exit(FatalError);
    }

    normal_ /= mag(normal_);

    const HashTable<pointField,string,Hash<string>>& vertexInfo =
        wallPlaneInfo::getWallPatchVerticesInfo();

    if (vertexInfo.found(wallName_))
    {
        vertices_ = vertexInfo[wallName_];
    }
    else
    {
        createLegacyVertices();
    }

    if (vertices_.size() < 3)
    {
        FatalErrorInFunction
            << "Finite wall '" << wallName_
            << "' needs at least three patch vertices."
            << exit(FatalError);
    }

    const boundBox rawBounds(vertices_, false);
    const scalar patchScale =
        maxScalarFWG(mag(rawBounds.span()), scalar(1e-3));

    // This is a predicate tolerance only.  It never expands the integrated
    // contact volume or area like a physical wall thickness would.
    tolerance_ = maxScalarFWG(scalar(1e-12), scalar(1e-10)*patchScale);

    forAll(vertices_, vertexI)
    {
        const scalar distance = signedDistance(vertices_[vertexI]);

        if (mag(distance) > 100*tolerance_)
        {
            FatalErrorInFunction
                << "Vertex " << vertices_[vertexI] << " of finite wall '"
                << wallName_ << "' is not on its plane. Signed distance: "
                << distance << ", tolerance: " << 100*tolerance_
                << exit(FatalError);
        }

        // Remove harmless dictionary round-off before building the local
        // polygon and all of its half-space tests.
        vertices_[vertexI] -= distance*normal_;
    }

    patchBounds_ = boundBox(vertices_, false);
    initialiseLocalGeometry();
}

//---------------------------------------------------------------------------//
std::shared_ptr<const finiteWallGeometry> finiteWallGeometry::lookup
(
    const string& wallName
)
{
    typedef std::shared_ptr<const finiteWallGeometry> geometryPtr;

    static HashTable<geometryPtr,string,Hash<string>> geometryCache;

    const HashTable<label,string,Hash<string>>& revisionInfo =
        wallPlaneInfo::getWallGeometryRevisionInfo();

    const label currentRevision =
        revisionInfo.found(wallName) ? revisionInfo[wallName] : 0;

    if
    (
        geometryCache.found(wallName)
     && geometryCache[wallName]->revision_ == currentRevision
    )
    {
        return geometryCache[wallName];
    }

    const geometryPtr geometry(new finiteWallGeometry(wallName));

    if (geometryCache.found(wallName))
    {
        geometryCache.erase(wallName);
    }

    geometryCache.insert(wallName, geometry);
    return geometry;
}

//---------------------------------------------------------------------------//
void finiteWallGeometry::createLegacyVertices()
{
    const HashTable<vector,string,Hash<string>>& minInfo =
        wallPlaneInfo::getWallMinBoundInfo();

    const HashTable<vector,string,Hash<string>>& maxInfo =
        wallPlaneInfo::getWallMaxBoundInfo();

    if (!minInfo.found(wallName_) || !maxInfo.found(wallName_))
    {
        FatalErrorInFunction
            << "Finite wall '" << wallName_
            << "' has neither vertices nor legacy minBound/maxBound data."
            << exit(FatalError);
    }

    label normalDirection = -1;
    label nonZeroDirections = 0;

    for (label dir = 0; dir < 3; dir++)
    {
        if (mag(normal_[dir]) > 1e-8)
        {
            normalDirection = dir;
            nonZeroDirections++;
        }
    }

    if (nonZeroDirections != 1)
    {
        FatalErrorInFunction
            << "Inclined finite wall '" << wallName_
            << "' must provide an ordered vertices list.  Legacy "
            << "minBound/maxBound is only unambiguous for axis-aligned walls."
            << exit(FatalError);
    }

    label tangentDirection0 = -1;
    label tangentDirection1 = -1;

    for (label dir = 0; dir < 3; dir++)
    {
        if (dir != normalDirection)
        {
            if (tangentDirection0 < 0)
            {
                tangentDirection0 = dir;
            }
            else
            {
                tangentDirection1 = dir;
            }
        }
    }

    const vector& minBound = minInfo[wallName_];
    const vector& maxBound = maxInfo[wallName_];

    const scalar lo0 = min(minBound[tangentDirection0], maxBound[tangentDirection0]);
    const scalar hi0 = max(minBound[tangentDirection0], maxBound[tangentDirection0]);
    const scalar lo1 = min(minBound[tangentDirection1], maxBound[tangentDirection1]);
    const scalar hi1 = max(minBound[tangentDirection1], maxBound[tangentDirection1]);

    vertices_.setSize(4);

    forAll(vertices_, vertexI)
    {
        vertices_[vertexI] = planePoint_;
    }

    vertices_[0][tangentDirection0] = lo0;
    vertices_[0][tangentDirection1] = lo1;

    vertices_[1][tangentDirection0] = hi0;
    vertices_[1][tangentDirection1] = lo1;

    vertices_[2][tangentDirection0] = hi0;
    vertices_[2][tangentDirection1] = hi1;

    vertices_[3][tangentDirection0] = lo0;
    vertices_[3][tangentDirection1] = hi1;
}

//---------------------------------------------------------------------------//
void finiteWallGeometry::initialiseLocalGeometry()
{
    vector firstEdge(vector::zero);

    for (label vertexI = 1; vertexI < vertices_.size(); vertexI++)
    {
        firstEdge = vertices_[vertexI] - vertices_[0];
        firstEdge -= (firstEdge & normal_)*normal_;

        if (mag(firstEdge) > tolerance_)
        {
            break;
        }
    }

    if (mag(firstEdge) <= tolerance_)
    {
        FatalErrorInFunction
            << "Finite wall '" << wallName_
            << "' has no non-degenerate polygon edge."
            << exit(FatalError);
    }

    tangent1_ = firstEdge/mag(firstEdge);
    tangent2_ = normal_ ^ tangent1_;
    tangent2_ /= mag(tangent2_);

    localVertices_.setSize(vertices_.size());

    forAll(vertices_, vertexI)
    {
        localVertices_[vertexI] = localCoordinates(vertices_[vertexI]);
    }

    scalar twiceArea = 0;

    forAll(localVertices_, vertexI)
    {
        const vector& a = localVertices_[vertexI];
        const vector& b =
            localVertices_[(vertexI + 1) % localVertices_.size()];

        twiceArea += a[0]*b[1] - b[0]*a[1];
    }

    if (mag(twiceArea) <= sqr(tolerance_))
    {
        FatalErrorInFunction
            << "Finite wall '" << wallName_ << "' has zero polygon area."
            << exit(FatalError);
    }

    // All containment tests below use the left half-plane of each edge.
    if (twiceArea < 0)
    {
        pointField reversedVertices(vertices_.size());
        List<vector> reversedLocalVertices(localVertices_.size());

        forAll(vertices_, vertexI)
        {
            reversedVertices[vertexI] =
                vertices_[vertices_.size() - 1 - vertexI];

            reversedLocalVertices[vertexI] =
                localVertices_[localVertices_.size() - 1 - vertexI];
        }

        vertices_ = reversedVertices;
        localVertices_ = reversedLocalVertices;
    }

    forAll(localVertices_, vertexI)
    {
        const vector& a = localVertices_[vertexI];
        const vector& b =
            localVertices_[(vertexI + 1) % localVertices_.size()];
        const vector& c =
            localVertices_[(vertexI + 2) % localVertices_.size()];

        if (mag(b - a) <= tolerance_)
        {
            FatalErrorInFunction
                << "Finite wall '" << wallName_
                << "' has a zero-length polygon edge at vertex "
                << vertexI << "."
                << exit(FatalError);
        }

        if (edgeValue(a, b, c) < -tolerance_)
        {
            FatalErrorInFunction
                << "Finite wall '" << wallName_
                << "' must be an ordered convex planar polygon."
                << exit(FatalError);
        }
    }
}

//---------------------------------------------------------------------------//
scalar finiteWallGeometry::edgeValue
(
    const vector& edgeStart,
    const vector& edgeEnd,
    const vector& checkedPoint
) const
{
    const scalar edgeU = edgeEnd[0] - edgeStart[0];
    const scalar edgeV = edgeEnd[1] - edgeStart[1];
    const scalar pointU = checkedPoint[0] - edgeStart[0];
    const scalar pointV = checkedPoint[1] - edgeStart[1];

    return
        (edgeU*pointV - edgeV*pointU)
       /(Foam::sqrt(sqr(edgeU) + sqr(edgeV)) + VSMALL);
}

//---------------------------------------------------------------------------//
scalar finiteWallGeometry::signedDistance(const point& checkedPoint) const
{
    return (checkedPoint - planePoint_) & normal_;
}

//---------------------------------------------------------------------------//
vector finiteWallGeometry::localCoordinates(const point& checkedPoint) const
{
    const vector relativePoint = checkedPoint - planePoint_;

    return vector
    (
        relativePoint & tangent1_,
        relativePoint & tangent2_,
        0
    );
}

//---------------------------------------------------------------------------//
point finiteWallGeometry::globalCoordinates(const vector& localPoint) const
{
    return
        planePoint_
      + localPoint[0]*tangent1_
      + localPoint[1]*tangent2_;
}

//---------------------------------------------------------------------------//
bool finiteWallGeometry::containsLocalPoint
(
    const vector& localPoint
) const
{
    forAll(localVertices_, edgeI)
    {
        const vector& edgeStart = localVertices_[edgeI];
        const vector& edgeEnd =
            localVertices_[(edgeI + 1) % localVertices_.size()];

        if (edgeValue(edgeStart, edgeEnd, localPoint) < -tolerance_)
        {
            return false;
        }
    }

    return true;
}

//---------------------------------------------------------------------------//
bool finiteWallGeometry::containsProjection
(
    const point& checkedPoint
) const
{
    return containsLocalPoint(localCoordinates(checkedPoint));
}

//---------------------------------------------------------------------------//
bool finiteWallGeometry::pointInContactRegion
(
    const point& checkedPoint
) const
{
    return
        signedDistance(checkedPoint) >= -tolerance_
     && containsProjection(checkedPoint);
}

//---------------------------------------------------------------------------//
bool finiteWallGeometry::planeIntersectsCell(const boundBox& cellBB) const
{
    const pointField cellPoints = cellBB.points();
    scalar minDistance = GREAT;
    scalar maxDistance = -GREAT;

    forAll(cellPoints, pointI)
    {
        const scalar distance = signedDistance(cellPoints[pointI]);
        minDistance = min(minDistance, distance);
        maxDistance = max(maxDistance, distance);
    }

    return
        minDistance <= tolerance_
     && maxDistance >= -tolerance_;
}

//---------------------------------------------------------------------------//
finiteWallGeometry::cellStatus finiteWallGeometry::classifyContactCell
(
    const boundBox& cellBB
) const
{
    const pointField cellPoints = cellBB.points();
    List<vector> localCellPoints(cellPoints.size());

    bool anyPenetratingPoint = false;
    bool allPenetratingPoints = true;

    forAll(cellPoints, pointI)
    {
        const scalar distance = signedDistance(cellPoints[pointI]);

        anyPenetratingPoint =
            anyPenetratingPoint || distance >= -tolerance_;

        allPenetratingPoints =
            allPenetratingPoints && distance >= -tolerance_;

        localCellPoints[pointI] = localCoordinates(cellPoints[pointI]);
    }

    if (!anyPenetratingPoint)
    {
        return cellOutside;
    }

    bool allSupportPoints = true;

    forAll(localVertices_, edgeI)
    {
        const vector& edgeStart = localVertices_[edgeI];
        const vector& edgeEnd =
            localVertices_[(edgeI + 1) % localVertices_.size()];

        bool anyPointInsideEdge = false;
        bool allPointsInsideEdge = true;

        forAll(localCellPoints, pointI)
        {
            const bool insideEdge =
                edgeValue
                (
                    edgeStart,
                    edgeEnd,
                    localCellPoints[pointI]
                ) >= -tolerance_;

            anyPointInsideEdge = anyPointInsideEdge || insideEdge;
            allPointsInsideEdge = allPointsInsideEdge && insideEdge;
        }

        // The polygon prism and the AABB are both convex.  If every AABB
        // corner violates one polygon half-space, they cannot intersect.
        if (!anyPointInsideEdge)
        {
            return cellOutside;
        }

        allSupportPoints = allSupportPoints && allPointsInsideEdge;
    }

    if (allPenetratingPoints && allSupportPoints)
    {
        return cellInside;
    }

    return cellMixed;
}

// ************************************************************************* //
