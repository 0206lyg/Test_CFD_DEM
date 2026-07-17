/*---------------------------------------------------------------------------*\
  Moving planar-wall state used only by collision patches that explicitly
  contain a motion sub-dictionary in HFDIBDEMDict.
\*---------------------------------------------------------------------------*/

#ifndef wallMotionInfo_H
#define wallMotionInfo_H

#include "HashTable.H"
#include "label.H"
#include "string.H"
#include "vector.H"
#include "word.H"

namespace Foam
{

class wallMotionState
{
    word sourcePatch_;
    vector previousPlanePoint_;
    vector nextPlanePoint_;
    vector velocity_;
    label timeIndex_;

public:

    wallMotionState()
    :
        sourcePatch_(),
        previousPlanePoint_(vector::zero),
        nextPlanePoint_(vector::zero),
        velocity_(vector::zero),
        timeIndex_(-1)
    {}

    wallMotionState
    (
        const word& sourcePatch,
        const vector& planePoint,
        const label timeIndex
    )
    :
        sourcePatch_(sourcePatch),
        previousPlanePoint_(planePoint),
        nextPlanePoint_(planePoint),
        velocity_(vector::zero),
        timeIndex_(timeIndex)
    {}

    const word& sourcePatch() const
    {
        return sourcePatch_;
    }

    const vector& previousPlanePoint() const
    {
        return previousPlanePoint_;
    }

    const vector& nextPlanePoint() const
    {
        return nextPlanePoint_;
    }

    const vector& velocity() const
    {
        return velocity_;
    }

    label timeIndex() const
    {
        return timeIndex_;
    }

    void setPreviousPlanePoint(const vector& point)
    {
        previousPlanePoint_ = point;
    }

    void setNextPlanePoint(const vector& point)
    {
        nextPlanePoint_ = point;
    }

    void setVelocity(const vector& velocity)
    {
        velocity_ = velocity;
    }

    void setTimeIndex(const label timeIndex)
    {
        timeIndex_ = timeIndex;
    }
};


class wallMotionInfo
{
    static HashTable<wallMotionState,string,Hash<string>> wallMotionInfo_;

public:

    static void clear()
    {
        wallMotionInfo_.clear();
    }

    static void insert
    (
        const string& wallName,
        const word& sourcePatch,
        const vector& planePoint,
        const label timeIndex
    )
    {
        wallMotionInfo_.insert
        (
            wallName,
            wallMotionState(sourcePatch, planePoint, timeIndex)
        );
    }

    static bool isMoving(const string& wallName)
    {
        return wallMotionInfo_.found(wallName);
    }

    static const vector& velocity(const string& wallName)
    {
        if (wallMotionInfo_.found(wallName))
        {
            return wallMotionInfo_[wallName].velocity();
        }

        return vector::zero;
    }

    static wallMotionState& state(const string& wallName)
    {
        return wallMotionInfo_[wallName];
    }

    static const HashTable<wallMotionState,string,Hash<string>>& all()
    {
        return wallMotionInfo_;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
