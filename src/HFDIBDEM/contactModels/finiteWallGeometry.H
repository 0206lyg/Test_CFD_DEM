/*---------------------------------------------------------------------------*\
  Geometry helper for finite planar DEM walls with arbitrary orientation.
\*---------------------------------------------------------------------------*/

#ifndef finiteWallGeometry_H
#define finiteWallGeometry_H

#include "boundBox.H"
#include "pointField.H"
#include "string.H"

#include <memory>

namespace Foam
{

class finiteWallGeometry
{
    public:

    enum cellStatus
    {
        cellOutside,
        cellInside,
        cellMixed
    };

    private:

    string wallName_;
    point planePoint_;
    vector normal_;
    vector tangent1_;
    vector tangent2_;
    pointField vertices_;
    List<vector> localVertices_;
    List<vector> localEdgeNormals_;
    List<vector> globalEdgeNormals_;
    List<scalar> edgeOffsets_;
    boundBox patchBounds_;
    scalar tolerance_;
    label revision_;

    void createLegacyVertices();
    void initialiseLocalGeometry();

    public:

    explicit finiteWallGeometry(const string& wallName);

    static std::shared_ptr<const finiteWallGeometry> lookup
    (
        const string& wallName
    );

    const string& wallName() const
    {
        return wallName_;
    }

    const point& planePoint() const
    {
        return planePoint_;
    }

    const vector& normal() const
    {
        return normal_;
    }

    const vector& tangent1() const
    {
        return tangent1_;
    }

    const vector& tangent2() const
    {
        return tangent2_;
    }

    const pointField& vertices() const
    {
        return vertices_;
    }

    const List<vector>& localVertices() const
    {
        return localVertices_;
    }

    const boundBox& patchBounds() const
    {
        return patchBounds_;
    }

    scalar tolerance() const
    {
        return tolerance_;
    }

    scalar signedDistance(const point& checkedPoint) const;

    vector localCoordinates(const point& checkedPoint) const;

    point globalCoordinates(const vector& localPoint) const;

    bool containsProjection(const point& checkedPoint) const;

    bool containsLocalPoint(const vector& localPoint) const;

    bool pointInContactRegion(const point& checkedPoint) const;

    void planeDistanceRange
    (
        const boundBox& cellBB,
        scalar& minDistance,
        scalar& maxDistance
    ) const;

    bool planeIntersectsCell(const boundBox& cellBB) const;

    cellStatus classifyContactCell(const boundBox& cellBB) const;

    bool cellMayIntersectContactRegion(const boundBox& cellBB) const
    {
        return classifyContactCell(cellBB) != cellOutside;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
