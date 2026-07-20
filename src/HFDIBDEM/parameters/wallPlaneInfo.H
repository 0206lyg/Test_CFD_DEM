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
    contactModel

Description
    Static class for storing wall material information

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef wallPlaneInfo_H
#define wallPlaneInfo_H

#include "pointField.H"
#include "boundBox.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class wallPlaneInfo
{
    static label axisAlignedNormalDirection(const vector& normal)
    {
        label normalDirection = -1;
        label nonZeroDirections = 0;

        for (label dir = 0; dir < 3; dir++)
        {
            if (mag(normal[dir]) > SMALL)
            {
                normalDirection = dir;
                nonZeroDirections++;
            }
        }

        return nonZeroDirections == 1 ? normalDirection : -1;
    }

    static bool isAxisAlignedRectangle
    (
        const vector& normal,
        const point& planePoint,
        const pointField& vertices
    )
    {
        const label normalDirection =
            axisAlignedNormalDirection(normal);

        if (normalDirection < 0 || vertices.size() != 4)
        {
            return false;
        }

        label tangentDirection[2] = {-1, -1};
        label tangentI = 0;

        for (label dir = 0; dir < 3; dir++)
        {
            if (dir != normalDirection)
            {
                tangentDirection[tangentI++] = dir;
            }
        }

        const boundBox bounds(vertices, false);
        const scalar tolerance = max
        (
            scalar(1e-12),
            scalar(1e-10)*max(mag(bounds.span()), scalar(1))
        );

        const scalar lo0 = bounds.min()[tangentDirection[0]];
        const scalar hi0 = bounds.max()[tangentDirection[0]];
        const scalar lo1 = bounds.min()[tangentDirection[1]];
        const scalar hi1 = bounds.max()[tangentDirection[1]];

        if ((hi0 - lo0) <= tolerance || (hi1 - lo1) <= tolerance)
        {
            return false;
        }

        bool cornerFound[4] = {false, false, false, false};

        forAll(vertices, vertexI)
        {
            if
            (
                mag
                (
                    vertices[vertexI][normalDirection]
                  - planePoint[normalDirection]
                ) > tolerance
            )
            {
                return false;
            }

            label side0 = -1;
            label side1 = -1;

            if
            (
                mag(vertices[vertexI][tangentDirection[0]] - lo0)
             <= tolerance
            )
            {
                side0 = 0;
            }
            else if
            (
                mag(vertices[vertexI][tangentDirection[0]] - hi0)
             <= tolerance
            )
            {
                side0 = 1;
            }

            if
            (
                mag(vertices[vertexI][tangentDirection[1]] - lo1)
             <= tolerance
            )
            {
                side1 = 0;
            }
            else if
            (
                mag(vertices[vertexI][tangentDirection[1]] - hi1)
             <= tolerance
            )
            {
                side1 = 1;
            }

            if (side0 < 0 || side1 < 0)
            {
                return false;
            }

            const label cornerI = 2*side0 + side1;

            if (cornerFound[cornerI])
            {
                return false;
            }

            cornerFound[cornerI] = true;
        }

        return
            cornerFound[0]
         && cornerFound[1]
         && cornerFound[2]
         && cornerFound[3];
    }

    public:

    static const HashTable<List<vector>,string,Hash<string>>& getWallPlaneInfo()
    {
        return wallPlaneInfo_;
    }

    static const HashTable<bool,string,Hash<string>>& getWallFiniteInfo()
    {
        return wallFiniteInfo_;
    }

    static const HashTable<vector,string,Hash<string>>& getWallMinBoundInfo()
    {
        return wallMinBoundInfo_;
    }

    static const HashTable<vector,string,Hash<string>>& getWallMaxBoundInfo()
    {
        return wallMaxBoundInfo_;
    }

    static const HashTable<pointField,string,Hash<string>>&
    getWallPatchVerticesInfo()
    {
        return wallPatchVerticesInfo_;
    }

    static const HashTable<label,string,Hash<string>>&
    getWallGeometryRevisionInfo()
    {
        return wallGeometryRevisionInfo_;
    }

    // The Git main fast finite-wall implementation is exact for a
    // coordinate-axis plane whose support is a global-axis rectangle.
    // minBound/maxBound patches qualify directly; vertices patches qualify
    // only after the four corners are proven to match that same rectangle.
    // Axis-aligned trapezoids and other polygons retain the general path.
    static bool usesLegacyAxisAlignedFinitePath(const string& wallPatch)
    {
        return
            wallLegacyFastPathInfo_.found(wallPatch)
         && wallLegacyFastPathInfo_[wallPatch];
    }

    static void touchWallGeometry(const string& wallPatch)
    {
        if (wallGeometryRevisionInfo_.found(wallPatch))
        {
            wallGeometryRevisionInfo_[wallPatch]++;
        }
        else
        {
            wallGeometryRevisionInfo_.insert(wallPatch, 1);
        }
    }

    static void setWallPlanePoint
    (
        const string& wallPatch,
        const vector& planePoint
    )
    {
        if (wallPlaneInfo_.found(wallPatch))
        {
            const vector displacement =
                planePoint - wallPlaneInfo_[wallPatch][1];

            wallPlaneInfo_[wallPatch][1] = planePoint;

            if
            (
                mag(displacement) > VSMALL
             && wallPatchVerticesInfo_.found(wallPatch)
            )
            {
                pointField& vertices = wallPatchVerticesInfo_[wallPatch];

                forAll(vertices, vertexI)
                {
                    vertices[vertexI] += displacement;
                }
            }

            if
            (
                wallFiniteInfo_.found(wallPatch)
             && wallFiniteInfo_[wallPatch]
             && wallMinBoundInfo_.found(wallPatch)
             && wallMaxBoundInfo_.found(wallPatch)
            )
            {
                wallMinBoundInfo_[wallPatch] += displacement;
                wallMaxBoundInfo_[wallPatch] += displacement;
            }

            if (mag(displacement) > VSMALL)
            {
                touchWallGeometry(wallPatch);
            }
        }
    }

    static void setWallFiniteBounds
    (
        const string& wallPatch,
        const vector& minBound,
        const vector& maxBound
    )
    {
        if
        (
            wallMinBoundInfo_.found(wallPatch)
         && wallMaxBoundInfo_.found(wallPatch)
        )
        {
            const vector oldMin = wallMinBoundInfo_[wallPatch];
            const vector oldMax = wallMaxBoundInfo_[wallPatch];
            const bool geometryChanged =
                mag(minBound - oldMin) > VSMALL
             || mag(maxBound - oldMax) > VSMALL;

            if (wallPatchVerticesInfo_.found(wallPatch))
            {
                pointField& vertices = wallPatchVerticesInfo_[wallPatch];
                const scalar tolerance = 1e-10;

                forAll(vertices, vertexI)
                {
                    for (label dir = 0; dir < 3; dir++)
                    {
                        if (mag(vertices[vertexI][dir] - oldMin[dir]) <= tolerance)
                        {
                            vertices[vertexI][dir] = minBound[dir];
                        }
                        else if
                        (
                            mag(vertices[vertexI][dir] - oldMax[dir])
                         <= tolerance
                        )
                        {
                            vertices[vertexI][dir] = maxBound[dir];
                        }
                    }
                }
            }

            wallMinBoundInfo_[wallPatch] = minBound;
            wallMaxBoundInfo_[wallPatch] = maxBound;

            if (geometryChanged)
            {
                touchWallGeometry(wallPatch);
            }
        }
    }

    static void wallPlaneInfo_insert
    (
        const string& wallPatch,
        const vector& patchNVec,
        const vector& planePoint
    )
    {
        List<vector> planeInfo{{patchNVec,planePoint}};

        wallPlaneInfo_.insert(wallPatch,planeInfo);
        wallFiniteInfo_.insert(wallPatch,false);
        wallMinBoundInfo_.insert(wallPatch,vector::zero);
        wallMaxBoundInfo_.insert(wallPatch,vector::zero);
        wallLegacyFastPathInfo_.insert(wallPatch,false);
        touchWallGeometry(wallPatch);
    };

    static void wallPlaneInfo_insert
    (
        const string& wallPatch,
        const vector& patchNVec,
        const vector& planePoint,
        const bool finitePatch,
        const vector& minBound,
        const vector& maxBound
    )
    {
        List<vector> planeInfo{{patchNVec,planePoint}};

        wallPlaneInfo_.insert(wallPatch,planeInfo);
        wallFiniteInfo_.insert(wallPatch,finitePatch);
        wallMinBoundInfo_.insert(wallPatch,minBound);
        wallMaxBoundInfo_.insert(wallPatch,maxBound);
        wallLegacyFastPathInfo_.insert
        (
            wallPatch,
            finitePatch && axisAlignedNormalDirection(patchNVec) >= 0
        );
        touchWallGeometry(wallPatch);
    };

    static void wallPlaneInfo_insert
    (
        const string& wallPatch,
        const vector& patchNVec,
        const vector& planePoint,
        const pointField& patchVertices
    )
    {
        List<vector> planeInfo{{patchNVec,planePoint}};
        const boundBox patchBounds(patchVertices, false);

        wallPlaneInfo_.insert(wallPatch,planeInfo);
        wallFiniteInfo_.insert(wallPatch,true);
        wallMinBoundInfo_.insert(wallPatch,patchBounds.min());
        wallMaxBoundInfo_.insert(wallPatch,patchBounds.max());
        wallPatchVerticesInfo_.insert(wallPatch,patchVertices);
        wallLegacyFastPathInfo_.insert
        (
            wallPatch,
            isAxisAlignedRectangle
            (
                patchNVec,
                planePoint,
                patchVertices
            )
        );
        touchWallGeometry(wallPatch);
    };

    static HashTable<List<vector>,string,Hash<string>> wallPlaneInfo_;

    static HashTable<bool,string,Hash<string>> wallFiniteInfo_;
    static HashTable<vector,string,Hash<string>> wallMinBoundInfo_;
    static HashTable<vector,string,Hash<string>> wallMaxBoundInfo_;
    static HashTable<bool,string,Hash<string>> wallLegacyFastPathInfo_;
    static HashTable<pointField,string,Hash<string>> wallPatchVerticesInfo_;
    static HashTable<label,string,Hash<string>> wallGeometryRevisionInfo_;
};

}

#endif

// ************************************************************************* //
