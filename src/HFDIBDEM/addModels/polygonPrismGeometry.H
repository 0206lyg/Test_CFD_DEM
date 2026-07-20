/*---------------------------------------------------------------------------*\
  Convex polygon-prism geometry and whole-body center sampler for addModels.
\*---------------------------------------------------------------------------*/

#ifndef polygonPrismGeometry_H
#define polygonPrismGeometry_H

#include "boundBox.H"
#include "dictionary.H"
#include "pointField.H"
#include "scalarField.H"

#include <vector>

namespace Foam
{

class polygonPrismGeometry
{
    typedef std::vector<vector> localPolygon;

    vector axis_;
    vector tangent1_;
    vector tangent2_;
    point planeOrigin_;
    pointField vertices_;
    List<vector> localVertices_;
    List<vector> inwardEdgeNormals_;
    List<scalar> edgeOffsets_;
    scalar minAxial_;
    scalar maxAxial_;
    scalar tolerance_;
    scalar area_;
    boundBox bounds_;

    vector localCoordinates(const point& checkedPoint) const;

    point globalCoordinates
    (
        const vector& localPoint,
        const scalar axialCoordinate
    ) const;

    localPolygon erodedPolygon
    (
        const pointField& relativeBodyPoints,
        const scalar isotropicRadius
    ) const;

    static scalar cross2D
    (
        const vector& a,
        const vector& b,
        const vector& c
    );

public:

    explicit polygonPrismGeometry(const dictionary& coeffs);

    const vector& axis() const
    {
        return axis_;
    }

    scalar minAxial() const
    {
        return minAxial_;
    }

    scalar maxAxial() const
    {
        return maxAxial_;
    }

    scalar area() const
    {
        return area_;
    }

    scalar volume() const
    {
        return area_*(maxAxial_ - minAxial_);
    }

    scalar tolerance() const
    {
        return tolerance_;
    }

    const boundBox& bounds() const
    {
        return bounds_;
    }

    const pointField& vertices() const
    {
        return vertices_;
    }

    bool contains(const point& checkedPoint) const;

    point sampleCenter
    (
        const pointField& relativeBodyPoints,
        const scalar isotropicRadius,
        const scalarField& unitRandomValues,
        bool& valid
    ) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
