/*---------------------------------------------------------------------------*\
  Convex polygon-prism geometry and whole-body center sampler for addModels.
\*---------------------------------------------------------------------------*/

#include "polygonPrismGeometry.H"

#include "error.H"

#include <cmath>

using namespace Foam;

namespace
{

scalar maxScalarPPG(const scalar a, const scalar b)
{
    return a > b ? a : b;
}

}

//---------------------------------------------------------------------------//
polygonPrismGeometry::polygonPrismGeometry(const dictionary& coeffs)
:
axis_(coeffs.lookup("axis")),
tangent1_(vector::zero),
tangent2_(vector::zero),
planeOrigin_(vector::zero),
vertices_(coeffs.lookup("vertices")),
localVertices_(),
inwardEdgeNormals_(),
edgeOffsets_(),
minAxial_(readScalar(coeffs.lookup("minAxial"))),
maxAxial_(readScalar(coeffs.lookup("maxAxial"))),
tolerance_(1e-12),
area_(0),
bounds_(boundBox())
{
    if (mag(axis_) <= VSMALL)
    {
        FatalIOErrorInFunction(coeffs)
            << "polygonPrism axis must be non-zero."
            << exit(FatalIOError);
    }

    axis_ /= mag(axis_);

    if (maxAxial_ <= minAxial_)
    {
        FatalIOErrorInFunction(coeffs)
            << "polygonPrism requires maxAxial > minAxial, got "
            << minAxial_ << " and " << maxAxial_ << "."
            << exit(FatalIOError);
    }

    if (vertices_.size() < 3)
    {
        FatalIOErrorInFunction(coeffs)
            << "polygonPrism requires at least three ordered vertices."
            << exit(FatalIOError);
    }

    const boundBox rawBounds(vertices_, false);
    const scalar geometryScale = maxScalarPPG
    (
        mag(rawBounds.span()),
        maxScalarPPG(maxAxial_ - minAxial_, scalar(1e-3))
    );
    tolerance_ = maxScalarPPG
    (
        scalar(1e-12),
        scalar(1e-10)*geometryScale
    );

    planeOrigin_ = vertices_[0];

    vector firstEdge(vector::zero);
    for (label vertexI = 1; vertexI < vertices_.size(); vertexI++)
    {
        firstEdge = vertices_[vertexI] - planeOrigin_;
        firstEdge -= (firstEdge & axis_)*axis_;

        if (mag(firstEdge) > tolerance_)
        {
            break;
        }
    }

    if (mag(firstEdge) <= tolerance_)
    {
        FatalIOErrorInFunction(coeffs)
            << "polygonPrism has no non-degenerate edge."
            << exit(FatalIOError);
    }

    tangent1_ = firstEdge/mag(firstEdge);
    tangent2_ = axis_ ^ tangent1_;
    tangent2_ /= mag(tangent2_);

    localVertices_.setSize(vertices_.size());
    scalar twiceArea = 0;

    forAll(vertices_, vertexI)
    {
        const scalar planeDistance =
            (vertices_[vertexI] - planeOrigin_) & axis_;

        if (mag(planeDistance) > 100*tolerance_)
        {
            FatalIOErrorInFunction(coeffs)
                << "polygonPrism vertex " << vertices_[vertexI]
                << " is not in the cross-section plane. Distance: "
                << planeDistance << "."
                << exit(FatalIOError);
        }

        vertices_[vertexI] -= planeDistance*axis_;
        localVertices_[vertexI] = localCoordinates(vertices_[vertexI]);
    }

    forAll(localVertices_, vertexI)
    {
        const vector& a = localVertices_[vertexI];
        const vector& b =
            localVertices_[(vertexI + 1) % localVertices_.size()];

        twiceArea += a[0]*b[1] - b[0]*a[1];
    }

    if (mag(twiceArea) <= sqr(tolerance_))
    {
        FatalIOErrorInFunction(coeffs)
            << "polygonPrism cross-section has zero area."
            << exit(FatalIOError);
    }

    // All edge predicates use the left half-plane.  Accept either winding
    // from the input but retain the requirement that the order is contiguous.
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
        twiceArea = -twiceArea;
    }

    area_ = 0.5*twiceArea;
    inwardEdgeNormals_.setSize(localVertices_.size());
    edgeOffsets_.setSize(localVertices_.size());

    forAll(localVertices_, edgeI)
    {
        const vector& a = localVertices_[edgeI];
        const vector& b =
            localVertices_[(edgeI + 1) % localVertices_.size()];

        const scalar edgeU = b[0] - a[0];
        const scalar edgeV = b[1] - a[1];
        const scalar edgeLength = Foam::sqrt(sqr(edgeU) + sqr(edgeV));

        if (edgeLength <= tolerance_)
        {
            FatalIOErrorInFunction(coeffs)
                << "polygonPrism has a zero-length edge at vertex "
                << edgeI << "."
                << exit(FatalIOError);
        }

        inwardEdgeNormals_[edgeI] = vector
        (
            -edgeV/edgeLength,
             edgeU/edgeLength,
             0
        );

        edgeOffsets_[edgeI] =
            -(inwardEdgeNormals_[edgeI] & a);

        // Checking every vertex against every edge rejects concave or
        // self-intersecting input before any body addition is attempted.
        forAll(localVertices_, checkedVertexI)
        {
            if
            (
                (inwardEdgeNormals_[edgeI]
               & localVertices_[checkedVertexI])
              + edgeOffsets_[edgeI] < -tolerance_
            )
            {
                FatalIOErrorInFunction(coeffs)
                    << "polygonPrism vertices must form an ordered convex "
                    << "polygon. Edge " << edgeI << " excludes vertex "
                    << checkedVertexI << "."
                    << exit(FatalIOError);
            }
        }
    }

    pointField prismPoints(2*vertices_.size());
    forAll(vertices_, vertexI)
    {
        const scalar vertexAxial = vertices_[vertexI] & axis_;
        prismPoints[2*vertexI] =
            vertices_[vertexI]
          + (minAxial_ - vertexAxial)*axis_;

        prismPoints[2*vertexI + 1] =
            vertices_[vertexI]
          + (maxAxial_ - vertexAxial)*axis_;
    }

    bounds_ = boundBox(prismPoints, false);
}

//---------------------------------------------------------------------------//
vector polygonPrismGeometry::localCoordinates
(
    const point& checkedPoint
) const
{
    const vector relativePoint = checkedPoint - planeOrigin_;

    return vector
    (
        relativePoint & tangent1_,
        relativePoint & tangent2_,
        0
    );
}

//---------------------------------------------------------------------------//
point polygonPrismGeometry::globalCoordinates
(
    const vector& localPoint,
    const scalar axialCoordinate
) const
{
    point inPlane =
        planeOrigin_
      + localPoint[0]*tangent1_
      + localPoint[1]*tangent2_;

    return
        inPlane
      + (axialCoordinate - (inPlane & axis_))*axis_;
}

//---------------------------------------------------------------------------//
scalar polygonPrismGeometry::cross2D
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

//---------------------------------------------------------------------------//
bool polygonPrismGeometry::contains(const point& checkedPoint) const
{
    const scalar axialCoordinate = checkedPoint & axis_;

    if
    (
        axialCoordinate < minAxial_ - tolerance_
     || axialCoordinate > maxAxial_ + tolerance_
    )
    {
        return false;
    }

    const vector localPoint = localCoordinates(checkedPoint);

    forAll(inwardEdgeNormals_, edgeI)
    {
        if
        (
            (inwardEdgeNormals_[edgeI] & localPoint)
          + edgeOffsets_[edgeI] < -tolerance_
        )
        {
            return false;
        }
    }

    return true;
}

//---------------------------------------------------------------------------//
polygonPrismGeometry::localPolygon polygonPrismGeometry::erodedPolygon
(
    const pointField& relativeBodyPoints,
    const scalar isotropicRadius
) const
{
    localPolygon clippedPolygon;
    clippedPolygon.reserve(localVertices_.size());

    List<vector> localRelativeBodyPoints(relativeBodyPoints.size());

    forAll(relativeBodyPoints, pointI)
    {
        localRelativeBodyPoints[pointI] = vector
        (
            relativeBodyPoints[pointI] & tangent1_,
            relativeBodyPoints[pointI] & tangent2_,
            0
        );
    }

    forAll(localVertices_, vertexI)
    {
        clippedPolygon.push_back(localVertices_[vertexI]);
    }

    forAll(inwardEdgeNormals_, edgeI)
    {
        if (clippedPolygon.empty())
        {
            break;
        }

        scalar minBodySupport = 0;

        if (localRelativeBodyPoints.size() > 0)
        {
            minBodySupport = GREAT;

            forAll(localRelativeBodyPoints, pointI)
            {
                minBodySupport = min
                (
                    minBodySupport,
                    inwardEdgeNormals_[edgeI]
                  & localRelativeBodyPoints[pointI]
                );
            }
        }

        const scalar shiftedOffset =
            edgeOffsets_[edgeI]
          + minBodySupport
          - isotropicRadius;

        localPolygon input(clippedPolygon);
        clippedPolygon.clear();

        vector previous = input.back();
        scalar previousValue =
            (inwardEdgeNormals_[edgeI] & previous)
          + shiftedOffset;
        bool previousInside = previousValue >= -tolerance_;

        for
        (
            localPolygon::const_iterator iter = input.begin();
            iter != input.end();
            ++iter
        )
        {
            const vector current = *iter;
            const scalar currentValue =
                (inwardEdgeNormals_[edgeI] & current)
              + shiftedOffset;
            const bool currentInside = currentValue >= -tolerance_;

            if (currentInside != previousInside)
            {
                const scalar denominator = previousValue - currentValue;

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

//---------------------------------------------------------------------------//
point polygonPrismGeometry::sampleCenter
(
    const pointField& relativeBodyPoints,
    const scalar isotropicRadius,
    const scalarField& unitRandomValues,
    bool& valid
) const
{
    valid = false;

    if (unitRandomValues.size() < 4 || isotropicRadius < 0)
    {
        return bounds_.midpoint();
    }

    const localPolygon centerPolygon =
        erodedPolygon(relativeBodyPoints, isotropicRadius);

    if (centerPolygon.size() < 3)
    {
        return bounds_.midpoint();
    }

    scalar minAxialSupport = 0;
    scalar maxAxialSupport = 0;

    if (relativeBodyPoints.size() > 0)
    {
        minAxialSupport = GREAT;
        maxAxialSupport = -GREAT;

        forAll(relativeBodyPoints, pointI)
        {
            const scalar support = relativeBodyPoints[pointI] & axis_;
            minAxialSupport = min(minAxialSupport, support);
            maxAxialSupport = max(maxAxialSupport, support);
        }
    }

    const scalar centerAxialMin =
        minAxial_ - minAxialSupport + isotropicRadius;

    const scalar centerAxialMax =
        maxAxial_ - maxAxialSupport - isotropicRadius;

    if (centerAxialMax < centerAxialMin - tolerance_)
    {
        return bounds_.midpoint();
    }

    scalar totalTriangleArea = 0;
    scalarField triangleAreas(centerPolygon.size() - 2, 0);

    for (label triangleI = 0; triangleI < label(triangleAreas.size()); triangleI++)
    {
        triangleAreas[triangleI] = 0.5*mag
        (
            cross2D
            (
                centerPolygon[0],
                centerPolygon[triangleI + 1],
                centerPolygon[triangleI + 2]
            )
        );

        totalTriangleArea += triangleAreas[triangleI];
    }

    if (totalTriangleArea <= sqr(tolerance_))
    {
        return bounds_.midpoint();
    }

    const scalar selector =
        min(max(unitRandomValues[0], scalar(0)), scalar(1))*totalTriangleArea;

    label selectedTriangle = -1;
    scalar cumulativeArea = 0;

    forAll(triangleAreas, triangleI)
    {
        if (triangleAreas[triangleI] <= sqr(tolerance_))
        {
            continue;
        }

        cumulativeArea += triangleAreas[triangleI];
        selectedTriangle = triangleI;

        if (selector < cumulativeArea)
        {
            break;
        }
    }

    if (selectedTriangle < 0)
    {
        return bounds_.midpoint();
    }

    const scalar rootRandom = Foam::sqrt
    (
        min(max(unitRandomValues[1], scalar(0)), scalar(1))
    );

    const scalar edgeRandom =
        min(max(unitRandomValues[2], scalar(0)), scalar(1));

    const vector& a = centerPolygon[0];
    const vector& b = centerPolygon[selectedTriangle + 1];
    const vector& c = centerPolygon[selectedTriangle + 2];

    const vector localCenter =
        (1 - rootRandom)*a
      + rootRandom*(1 - edgeRandom)*b
      + rootRandom*edgeRandom*c;

    const scalar axialRandom =
        min(max(unitRandomValues[3], scalar(0)), scalar(1));

    const scalar axialCenter =
        centerAxialMin
      + axialRandom*(centerAxialMax - centerAxialMin);

    valid = true;
    return globalCoordinates(localCenter, axialCenter);
}

// ************************************************************************* //
