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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class wallPlaneInfo
{
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

    static const HashTable<scalar,string,Hash<string>>& getWallPadInfo()
    {
        return wallPadInfo_;
    }

    static void setWallPlanePoint
    (
        const string& wallPatch,
        const vector& planePoint
    )
    {
        if (wallPlaneInfo_.found(wallPatch))
        {
            wallPlaneInfo_[wallPatch][1] = planePoint;
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
            wallMinBoundInfo_[wallPatch] = minBound;
            wallMaxBoundInfo_[wallPatch] = maxBound;
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
        wallPadInfo_.insert(wallPatch,0.0);
    };

    static void wallPlaneInfo_insert
    (
        const string& wallPatch,
        const vector& patchNVec,
        const vector& planePoint,
        const bool finitePatch,
        const vector& minBound,
        const vector& maxBound,
        const scalar pad
    )
    {
        List<vector> planeInfo{{patchNVec,planePoint}};

        wallPlaneInfo_.insert(wallPatch,planeInfo);
        wallFiniteInfo_.insert(wallPatch,finitePatch);
        wallMinBoundInfo_.insert(wallPatch,minBound);
        wallMaxBoundInfo_.insert(wallPatch,maxBound);
        wallPadInfo_.insert(wallPatch,pad);
    };

    static HashTable<List<vector>,string,Hash<string>> wallPlaneInfo_;

    static HashTable<bool,string,Hash<string>> wallFiniteInfo_;
    static HashTable<vector,string,Hash<string>> wallMinBoundInfo_;
    static HashTable<vector,string,Hash<string>> wallMaxBoundInfo_;
    static HashTable<scalar,string,Hash<string>> wallPadInfo_;
};

}

#endif

// ************************************************************************* //
