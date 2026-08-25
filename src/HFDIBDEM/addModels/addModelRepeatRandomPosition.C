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
    Foam

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/
#include "addModelRepeatRandomPosition.H"
#include "meshSearch.H"

#include <cmath>

using namespace Foam;

//---------------------------------------------------------------------------//
addModelRepeatRandomPosition::addModelRepeatRandomPosition
(
    const dictionary& addModelDict,
    const Foam::fvMesh& mesh,
    std::unique_ptr<geomModel> bodyGeomModel,
    List<labelList>& cellPoints
)
:
addModel(mesh, std::move(bodyGeomModel), cellPoints),
addModelDict_(addModelDict),
addMode_(word(addModelDict_.lookup("addModel"))),
bodyAdded_(false),
coeffsDict_(addModelDict_.subDict(addMode_+"Coeffs")),

addDomain_(word(coeffsDict_.lookup("addDomain"))),
scalingMode_(word(coeffsDict_.lookup("scalingMode"))),
rotationMode_(word(coeffsDict_.lookup("rotationMode"))),
addModeI_(word(coeffsDict_.lookup("addMode"))),

addDomainCoeffs_(coeffsDict_.subDict(addDomain_ + "Coeffs")),
scalingModeCoeffs_(coeffsDict_.subDict(scalingMode_ + "Coeffs")),
rotationModeCoeffs_(coeffsDict_.subDict(rotationMode_ + "Coeffs")),
addModeICoeffs_(coeffsDict_.subDict(addModeI_ + "Coeffs")),

useNTimes_(0),
timeBetweenUsage_(0),
partPerAdd_(0),
fieldValue_(0),
addedOnTimeLevel_(0),
partPerAddTemp_(0),

zoneName_(),
minBound_(vector::zero),
maxBound_(vector::zero),

scaleParticles_(false),
minScale_(0),
maxScale_(0),
minScaleFit_(0),
scaleStep_(0),
nTriesBeforeScaling_(0),

rotateParticles_(false),
randomAxis_(false),
axisOfRot_(vector::zero),

bodyAdditionAttemptCounter_(0),
scaleCorrectionCounter_(0),

scaleApplication_(false),
scaleRandomApplication_(false),
rescaleRequirement_(false),
succesfulladition_(false),
scalingFactor_(0),
restartPartCountTemp_(false),
reapeatedAddition_(false),
firstTimeRunning_(true),
cellZoneActive_(false),
boundBoxActive_(false),
octreeField_(mesh_.nCells(), 0),
timeBased_(false),
fieldBased_(false),
fieldCurrentValue_(0),
startAfterTime_(-GREAT),
placementStrategy_(
    coeffsDict_.lookupOrDefault<word>("placementStrategy", "legacy")
),
clearanceScale_(
    coeffsDict_.lookupOrDefault<scalar>("clearanceScale", 1.0)
),
maxCandidateAttemptsPerTimeStep_(
    coeffsDict_.lookupOrDefault<label>
    (
        "maxCandidateAttemptsPerTimeStep",
        50
    )
),
candidateBatchSize_(
    coeffsDict_.lookupOrDefault<label>("candidateBatchSize", 16)
),
attemptsOnTimeLevel_(0),
attemptTimeIndex_(-1),
exactPlacementExhausted_(false),
containingFaces_(6, true),
insertionMinBound_(vector::zero),
insertionMaxBound_(vector::zero),
exactPlacementIsSphere_(false),
exactPlacementRadius_(0),
exactBinSize_(0),
exactBinMin_(vector::zero),
exactBinMax_(vector::zero),
nExactBinsX_(1),
nExactBinsY_(1),
nExactBinsZ_(1),
exactBinsTimeIndex_(-1),
exactBinsBodyCount_(-1),
exactBodyBins_(),
allActiveCellsInMesh_(true),
randGen_(
    placementStrategy_ == "continuousExact"
  ? coeffsDict_.lookupOrDefault<label>("randomSeed", label(12345))
  : label(clock::getTime())
)
{
	init();
}

addModelRepeatRandomPosition::~addModelRepeatRandomPosition()
{
}

//---------------------------------------------------------------------------//
void addModelRepeatRandomPosition::init()
{
    // Optional guard for inlet/replenishment models.  With startAfterTime 0,
    // the model does not alter an independently generated t=0 packing, but it
    // becomes active as soon as simulation time advances beyond zero.
    startAfterTime_ = coeffsDict_.lookupOrDefault<scalar>
    (
        "startAfterTime",
        -GREAT
    );

    if
    (
        placementStrategy_ != "legacy"
     && placementStrategy_ != "continuousExact"
    )
    {
        FatalErrorInFunction
            << "Unknown repeatRandomPosition placementStrategy '"
            << placementStrategy_ << "'. Valid values are legacy and "
            << "continuousExact." << nl << exit(FatalError);
    }

    if (placementStrategy_ == "continuousExact")
    {
        if (addDomain_ != "boundBox")
        {
            FatalErrorInFunction
                << "continuousExact currently requires addDomain boundBox; "
                << "got " << addDomain_ << nl << exit(FatalError);
        }
        if (scalingMode_ != "noScaling")
        {
            FatalErrorInFunction
                << "continuousExact supports only scalingMode noScaling; "
                << "got " << scalingMode_ << nl << exit(FatalError);
        }
        if (mag(clearanceScale_ - 1.0) > SMALL)
        {
            FatalErrorInFunction
                << "continuousExact requires clearanceScale 1.0; got "
                << clearanceScale_ << nl << exit(FatalError);
        }
        if (maxCandidateAttemptsPerTimeStep_ < 1)
        {
            FatalErrorInFunction
                << "maxCandidateAttemptsPerTimeStep must be positive; got "
                << maxCandidateAttemptsPerTimeStep_ << nl
                << exit(FatalError);
        }
        if (candidateBatchSize_ < 1)
        {
            FatalErrorInFunction
                << "candidateBatchSize must be positive; got "
                << candidateBatchSize_ << nl << exit(FatalError);
        }

        if (coeffsDict_.found("nonContainingFaces"))
        {
            const wordList nonContainingFaces
            (
                coeffsDict_.lookup("nonContainingFaces")
            );

            forAll(nonContainingFaces, faceI)
            {
                const word& faceName = nonContainingFaces[faceI];
                label faceIndex(-1);
                if      (faceName == "xMin") faceIndex = 0;
                else if (faceName == "xMax") faceIndex = 1;
                else if (faceName == "yMin") faceIndex = 2;
                else if (faceName == "yMax") faceIndex = 3;
                else if (faceName == "zMin") faceIndex = 4;
                else if (faceName == "zMax") faceIndex = 5;
                else
                {
                    FatalErrorInFunction
                        << "Unknown nonContainingFaces entry '" << faceName
                        << "'. Valid names are xMin xMax yMin yMax zMin "
                        << "zMax." << nl << exit(FatalError);
                }
                containingFaces_[faceIndex] = false;
            }
        }

        InfoH << addModel_Info << "-- addModelMessage-- "
            << "continuousExact placement, clearanceScale 1.0, seed "
            << coeffsDict_.lookupOrDefault<label>("randomSeed", label(12345))
            << ", max candidate attempts per time step "
            << maxCandidateAttemptsPerTimeStep_
            << ", candidate batch size " << candidateBatchSize_ << endl;
    }

    // set sizes to necessary datatypes
    cellsInBoundBox_.setSize(Pstream::nProcs());
    cellZonePoints_.setSize(Pstream::nProcs());


	if (addModeI_ == "timeBased")
	{
        useNTimes_ = (readLabel(addModeICoeffs_.lookup("useNTimes")));
		timeBetweenUsage_ = (readScalar(addModeICoeffs_.lookup("timeBetweenUsage")));
		partPerAdd_ = (readLabel(addModeICoeffs_.lookup("partPerAdd")));
        timeBased_ = true;
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "addModel will control simulation time" << endl;
        InfoH << "-- addModelMessage-- " << "STL will be re-used "
            << useNTimes_ << " times" << endl;
        InfoH << "-- addModelMessage-- " << "STL will be added each "
            << timeBetweenUsage_ << " [T]" << endl;
        InfoH << "-- addModelMessage-- " << "upon each addition, "
            << partPerAdd_ << " bodies will be generated from the given STL"
            << endl;
	}
	else if (addModeI_ == "fieldBased")
	{
		fieldValue_ = (readScalar(addModeICoeffs_.lookup("fieldValue")));
        fieldBased_ = true;
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "addModel will control particles volume fraction" << endl;
		InfoH << "-- addModelMessage-- "
            << "preset volume fraction: " << fieldValue_ << endl;
	}
    else
    {
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "notImplemented, will crash" << endl;
    }

	if (addDomain_ == "cellZone")
	{
		zoneName_ = (word(addDomainCoeffs_.lookup("zoneName")));
		cellZoneActive_ = true;
        initializeCellZone();
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "cellZone based addition zone" << endl;
	}
	else if (addDomain_ == "boundBox")
	{
		minBound_       = (addDomainCoeffs_.lookup("minBound"));
		maxBound_       = (addDomainCoeffs_.lookup("maxBound"));
		boundBoxActive_ = true;
        initializeBoundBox();
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "boundBox based addition zone" << endl;
	}
	else if (addDomain_ == "domain")
	{
		InfoH << addModel_Info << "-- addModelMessage-- "
            << "notImplemented, will crash" << endl;
	}
    else
    {
		InfoH << addModel_Info << "-- addModelMessage-- "
            << "notImplemented, will crash" << endl;
	}

    if (placementStrategy_ == "continuousExact")
    {
        insertionMinBound_ = coeffsDict_.lookupOrDefault<vector>
        (
            "insertionMinBound",
            minBound_
        );
        insertionMaxBound_ = coeffsDict_.lookupOrDefault<vector>
        (
            "insertionMaxBound",
            maxBound_
        );

        for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
        {
            if (insertionMaxBound_[cmpt] <= insertionMinBound_[cmpt])
            {
                FatalErrorInFunction
                    << "Invalid continuousExact insertion bounds: min="
                    << insertionMinBound_ << ", max="
                    << insertionMaxBound_ << nl << exit(FatalError);
            }
        }

        InfoH << addModel_Info << "-- addModelMessage-- "
            << "continuousExact control bounds " << minBound_ << " to "
            << maxBound_ << ", insertion bounds " << insertionMinBound_
            << " to " << insertionMaxBound_ << endl;
    }

    // check, if the whole zone is in the mesh
    scalarList procZoneVols(Pstream::nProcs());
    procZoneVols[Pstream::myProcNo()] = 0;
    forAll (cellsInBoundBox_[Pstream::myProcNo()],cellI)
    {
        procZoneVols[Pstream::myProcNo()]+=mesh_.V()[cellsInBoundBox_[Pstream::myProcNo()][cellI]];
    }
    scalar zoneVol(gSum(procZoneVols));
    scalar zoneBBoxVol(cellZoneBounds_.volume());
    if (zoneVol - zoneBBoxVol > 1e-5*zoneBBoxVol)
    {
        allActiveCellsInMesh_ = false;
        InfoH << addModel_Info << "-- addModelMessage-- "
             << "addition zone NOT completely immersed in mesh "
             << "this computation will be EXPENSIVE" << endl;
        InfoH << zoneVol << " " << zoneBBoxVol << endl;
    }
    else
    {
        InfoH << addModel_Info << "-- addModelMessage-- "
             << "addition zone completely immersed in mesh -> OK" << endl;
    }

	if (scalingMode_ == "noScaling")
	{
		scaleParticles_ = false;
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "all particles will have the same scale" << endl;
	}
	else if (scalingMode_ == "randomScaling")
	{
		minScale_               = (readScalar(scalingModeCoeffs_.lookup("minScale")));
		maxScale_               = (readScalar(scalingModeCoeffs_.lookup("maxScale")));
		scaleParticles_         = false;
		scaleRandomApplication_ = true;
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "particles will be randomly scaled" << endl;
	}
	else if (scalingMode_ == "scaleToFit")
	{
		minScaleFit_        = (readScalar(scalingModeCoeffs_.lookup("minScale")));
		scaleStep_          = (readScalar(scalingModeCoeffs_.lookup("scaleStep")));
		nTriesBeforeScaling_= (readScalar(scalingModeCoeffs_.lookup("nTriesBeforeScaling")));
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "particles will be downscaled to better fill the domain" << endl;
		InfoH << "-- addModelMessage-- "
            << "nTriesBeforeDownScaling: " << nTriesBeforeScaling_ << endl;
		scaleParticles_ = true;
	}
    else
    {
		InfoH << addModel_Info << "-- addModelMessage-- "
            << "notImplemented, will crash" << endl;
	}

	if (rotationMode_ == "noRotation")
	{
		rotateParticles_ = false;
		randomAxis_      = false;
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "source STL will not be rotated upon addition" << endl;
	}
	else if (rotationMode_ == "randomRotation")
	{
		rotateParticles_ = true;
		randomAxis_      = true;
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "source STL will be randomly rotated upon addition" << endl;
	}
	else if (rotationMode_ == "uniformRandomRotation")
	{
		rotateParticles_ = true;
		randomAxis_      = false;
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "source STL orientation will be sampled uniformly in 3D "
            << "upon addition" << endl;
	}
	else if (rotationMode_ == "fixedAxisRandomRotation")
	{
		axisOfRot_       = (rotationModeCoeffs_.lookup("axis"));
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "source STL will be rotated by a random angle around a fixed axis upon addition" << endl;
		InfoH << "-- addModelMessage-- " << "set rotation axis: "
            << axisOfRot_ << endl;
		rotateParticles_ = true;
		randomAxis_      = false;
	}
    else
    {
		InfoH << addModel_Info << "-- addModelMessage-- "
            << "notImplemented, will crash" << endl;
	}

    if (placementStrategy_ == "continuousExact")
    {
        initializeContinuousExactPlacement();
    }

	partPerAddTemp_ = partPerAdd_;
}

//---------------------------------------------------------------------------//
bool addModelRepeatRandomPosition::shouldAddBody(const volScalarField& body)
{
    if (mesh_.time().value() <= startAfterTime_ + SMALL)
    {
        return false;
    }

    if (placementStrategy_ == "continuousExact")
    {
        resetExactPlacementBudget();
        if
        (
            exactPlacementExhausted_
         || attemptsOnTimeLevel_ >= maxCandidateAttemptsPerTimeStep_
        )
        {
            return false;
        }
    }

    if (timeBased_)
    {
        scalar timeVal(mesh_.time().value());
        scalar deltaTime(mesh_.time().deltaT().value());
        scalar tmFrac(timeVal/timeBetweenUsage_);
        tmFrac -=  floor(tmFrac+deltaTime);

        InfoH << addModel_Info << "-- addModelMessage-- "
            << "Time/(Time beween usage) - floor(Time/Time beween usage): "
            << tmFrac << endl;

        InfoH << addModel_Info << "-- addModelMessage-- "
            << "Number of bodies added on this time level: "
            << addedOnTimeLevel_ << endl;

        bool tmLevelOk(tmFrac < deltaTime);

        if (not tmLevelOk)
        {
            addedOnTimeLevel_ = 0;
            return false;
        }

        if (partPerAdd_ <= addedOnTimeLevel_) {return false;}

        return (tmLevelOk and useNTimes_ > 0);
    }

    if (fieldBased_)
    {
        scalar currentLambdaFrac(checkLambdaFraction(body));
        if (currentLambdaFrac < fieldValue_ )
        {
            InfoH << addModel_Info << "-- addModelMessage-- "
                << "Current lambda fraction = " << currentLambdaFrac
                << " < then preset lambda fraction = " << fieldValue_ << endl;
            return true;
        }
    }

    return false;

}

//---------------------------------------------------------------------------//
void addModelRepeatRandomPosition::resetExactPlacementBudget()
{
    const label currentTimeIndex = mesh_.time().timeIndex();
    if (currentTimeIndex != attemptTimeIndex_)
    {
        attemptTimeIndex_ = currentTimeIndex;
        attemptsOnTimeLevel_ = 0;
        exactPlacementExhausted_ = false;
    }
}

//---------------------------------------------------------------------------//
void addModelRepeatRandomPosition::initializeContinuousExactPlacement()
{
    const contactType geometryType = geomModel_->getcType();
    exactPlacementIsSphere_ = geometryType == sphere;

    if (!exactPlacementIsSphere_ && geometryType != convex)
    {
        FatalErrorInFunction
            << "continuousExact supports sphere and convex geometry only; "
            << "the requested geometry has contact type " << geometryType
            << nl << exit(FatalError);
    }

    if (exactPlacementIsSphere_)
    {
        const vector sourceSpan(geomModel_->getBounds().span());
        exactPlacementRadius_ = 0.5*Foam::max
        (
            sourceSpan.x(),
            Foam::max(sourceSpan.y(), sourceSpan.z())
        );
    }
    else
    {
        const pointField sourcePoints(geomModel_->getBodyPoints());
        const vector sourceCentre(geomModel_->getCoM());

        if (sourcePoints.empty())
        {
            FatalErrorInFunction
                << "continuousExact convex geometry returned no support "
                << "points." << nl << exit(FatalError);
        }

        scalar radiusSquared(0);
        forAll(sourcePoints, pointI)
        {
            radiusSquared = Foam::max
            (
                radiusSquared,
                magSqr(sourcePoints[pointI] - sourceCentre)
            );
        }
        exactPlacementRadius_ = Foam::sqrt(radiusSquared);
    }

    if (exactPlacementRadius_ <= VSMALL)
    {
        FatalErrorInFunction
            << "continuousExact obtained a non-positive source radius "
            << exactPlacementRadius_ << nl << exit(FatalError);
    }

    exactBinSize_ = 2.0*exactPlacementRadius_;
    const vector radiusVector
    (
        exactPlacementRadius_,
        exactPlacementRadius_,
        exactPlacementRadius_
    );

    // The grid encloses every candidate support, including support allowed to
    // cross a non-containing insertion face.  Existing bodies are inserted
    // by all AABB-overlapped bins, so larger or differently shaped bodies are
    // still found without assuming a common radius.
    exactBinMin_ = insertionMinBound_ - radiusVector;
    exactBinMax_ = insertionMaxBound_ + radiusVector;
    const vector binSpan(exactBinMax_ - exactBinMin_);

    nExactBinsX_ = Foam::max
    (
        label(1),
        label(std::ceil(binSpan.x()/exactBinSize_))
    );
    nExactBinsY_ = Foam::max
    (
        label(1),
        label(std::ceil(binSpan.y()/exactBinSize_))
    );
    nExactBinsZ_ = Foam::max
    (
        label(1),
        label(std::ceil(binSpan.z()/exactBinSize_))
    );

    InfoH << addModel_Info << "-- addModelMessage-- "
        << "continuousExact dynamic spatial hash: "
        << nExactBinsX_ << 'x' << nExactBinsY_ << 'x' << nExactBinsZ_
        << " bins, bin size " << exactBinSize_ << endl;
}

//---------------------------------------------------------------------------//
label addModelRepeatRandomPosition::exactBinKey
(
    const label ix,
    const label iy,
    const label iz
) const
{
    return ix + nExactBinsX_*(iy + nExactBinsY_*iz);
}

//---------------------------------------------------------------------------//
label addModelRepeatRandomPosition::exactBinCoordinate
(
    const scalar coordinate,
    const direction cmpt
) const
{
    const label nBins = cmpt == 0
      ? nExactBinsX_
      : (cmpt == 1 ? nExactBinsY_ : nExactBinsZ_);

    label binI = label
    (
        std::floor((coordinate - exactBinMin_[cmpt])/exactBinSize_)
    );
    binI = Foam::max(label(0), Foam::min(binI, nBins - 1));
    return binI;
}

//---------------------------------------------------------------------------//
void addModelRepeatRandomPosition::rebuildExactBodyBins
(
    PtrList<immersedBody>& immersedBodies
)
{
    const label currentTimeIndex(mesh_.time().timeIndex());
    if
    (
        exactBinsTimeIndex_ == currentTimeIndex
     && exactBinsBodyCount_ == immersedBodies.size()
    )
    {
        return;
    }

    exactBodyBins_.clear();
    const boundBox hashBounds(exactBinMin_, exactBinMax_);

    forAll(immersedBodies, bodyI)
    {
        const boundBox bodyBounds
        (
            immersedBodies[bodyI].getGeomModel().getBounds()
        );

        if (!bodyBounds.overlaps(hashBounds))
        {
            continue;
        }

        const label minX = exactBinCoordinate(bodyBounds.min().x(), 0);
        const label maxX = exactBinCoordinate(bodyBounds.max().x(), 0);
        const label minY = exactBinCoordinate(bodyBounds.min().y(), 1);
        const label maxY = exactBinCoordinate(bodyBounds.max().y(), 1);
        const label minZ = exactBinCoordinate(bodyBounds.min().z(), 2);
        const label maxZ = exactBinCoordinate(bodyBounds.max().z(), 2);

        for (label iz = minZ; iz <= maxZ; ++iz)
        {
            for (label iy = minY; iy <= maxY; ++iy)
            {
                for (label ix = minX; ix <= maxX; ++ix)
                {
                    const label key(exactBinKey(ix, iy, iz));
                    if (!exactBodyBins_.found(key))
                    {
                        exactBodyBins_.insert(key, DynamicLabelList());
                    }
                    exactBodyBins_[key].append(bodyI);
                }
            }
        }
    }

    exactBinsTimeIndex_ = currentTimeIndex;
    exactBinsBodyCount_ = immersedBodies.size();
}

//---------------------------------------------------------------------------//
bool addModelRepeatRandomPosition::exactCandidateInContact
(
    PtrList<immersedBody>& immersedBodies
) const
{
    const boundBox candidateBounds(geomModel_->getBounds());
    const label minX = exactBinCoordinate(candidateBounds.min().x(), 0);
    const label maxX = exactBinCoordinate(candidateBounds.max().x(), 0);
    const label minY = exactBinCoordinate(candidateBounds.min().y(), 1);
    const label maxY = exactBinCoordinate(candidateBounds.max().y(), 1);
    const label minZ = exactBinCoordinate(candidateBounds.min().z(), 2);
    const label maxZ = exactBinCoordinate(candidateBounds.max().z(), 2);
    pointField candidatePoints;
    if (!exactPlacementIsSphere_)
    {
        candidatePoints = geomModel_->getBodyPoints();
    }

    labelHashSet checkedBodies;
    bool localContact(false);

    for (label iz = minZ; iz <= maxZ && !localContact; ++iz)
    {
        for (label iy = minY; iy <= maxY && !localContact; ++iy)
        {
            for (label ix = minX; ix <= maxX && !localContact; ++ix)
            {
                const label key(exactBinKey(ix, iy, iz));
                if (!exactBodyBins_.found(key))
                {
                    continue;
                }

                const DynamicLabelList& neighbours(exactBodyBins_[key]);
                forAll(neighbours, neighbourI)
                {
                    const label bodyI(neighbours[neighbourI]);
                    if (checkedBodies.found(bodyI))
                    {
                        continue;
                    }
                    checkedBodies.insert(bodyI);

                    geomModel& otherGeometry
                    (
                        immersedBodies[bodyI].getGeomModel()
                    );
                    const boundBox otherBounds(otherGeometry.getBounds());
                    if (!candidateBounds.overlaps(otherBounds))
                    {
                        continue;
                    }

                    if
                    (
                        exactPlacementIsSphere_
                     && otherGeometry.getcType() == sphere
                    )
                    {
                        const vector otherSpan(otherBounds.span());
                        const scalar otherRadius = 0.5*Foam::max
                        (
                            otherSpan.x(),
                            Foam::max(otherSpan.y(), otherSpan.z())
                        );
                        if
                        (
                            magSqr
                            (
                                geomModel_->getCoM()
                              - otherGeometry.getCoM()
                            )
                          <= sqr(exactPlacementRadius_ + otherRadius)
                        )
                        {
                            localContact = true;
                            break;
                        }
                    }
                    else if
                    (
                        !exactPlacementIsSphere_
                     && otherGeometry.getcType() == convex
                    )
                    {
                        const pointField otherPoints
                        (
                            otherGeometry.getBodyPoints()
                        );

                        if (!candidatePoints.empty() && !otherPoints.empty())
                        {
                            scalar otherRadiusSquared(0);
                            const vector otherCentre(otherGeometry.getCoM());
                            forAll(otherPoints, pointI)
                            {
                                otherRadiusSquared = Foam::max
                                (
                                    otherRadiusSquared,
                                    magSqr
                                    (
                                        otherPoints[pointI] - otherCentre
                                    )
                                );
                            }

                            const scalar lengthScale =
                                exactPlacementRadius_
                              + Foam::sqrt(otherRadiusSquared);

                            if
                            (
                                convexBodiesIntersect
                                (
                                    candidatePoints,
                                    otherPoints,
                                    geomModel_->getCoM() - otherCentre,
                                    lengthScale
                                )
                            )
                            {
                                localContact = true;
                                break;
                            }
                        }
                    }
                    // Mixed, cluster and non-convex pairs are deliberately
                    // left to isBodyInContact(), the framework's final guard.
                }
            }
        }
    }

    // All ranks must take the same branch before the next collective random
    // draw or contact check.
    reduce(localContact, orOp<bool>());
    return localContact;
}

//---------------------------------------------------------------------------//
vector addModelRepeatRandomPosition::minkowskiSupport
(
    const pointField& first,
    const pointField& second,
    const vector& direction
)
{
    label firstI(0);
    label secondI(0);
    scalar firstProjection(first[0] & direction);
    scalar secondProjection(second[0] & (-direction));

    for (label pointI = 1; pointI < first.size(); ++pointI)
    {
        const scalar projection(first[pointI] & direction);
        if (projection > firstProjection)
        {
            firstProjection = projection;
            firstI = pointI;
        }
    }
    for (label pointI = 1; pointI < second.size(); ++pointI)
    {
        const scalar projection(second[pointI] & (-direction));
        if (projection > secondProjection)
        {
            secondProjection = projection;
            secondI = pointI;
        }
    }
    return first[firstI] - second[secondI];
}

//---------------------------------------------------------------------------//
bool addModelRepeatRandomPosition::sameDirection
(
    const vector& first,
    const vector& second
)
{
    return (first & second) > 0;
}

//---------------------------------------------------------------------------//
bool addModelRepeatRandomPosition::updateSimplex
(
    DynamicList<vector>& simplex,
    vector& direction
)
{
    if (simplex.size() == 2)
    {
        const vector A(simplex[1]);
        const vector B(simplex[0]);
        const vector AO(-A);
        const vector AB(B - A);

        if (sameDirection(AB, AO))
        {
            direction = (AB ^ AO) ^ AB;
            if (magSqr(direction) <= VSMALL)
            {
                return true;
            }
        }
        else
        {
            simplex.clear();
            simplex.append(A);
            direction = AO;
        }
        return magSqr(direction) <= VSMALL;
    }

    if (simplex.size() == 3)
    {
        const vector A(simplex[2]);
        const vector B(simplex[1]);
        const vector C(simplex[0]);
        const vector AO(-A);
        const vector AB(B - A);
        const vector AC(C - A);
        const vector ABC(AB ^ AC);

        if (sameDirection(ABC ^ AC, AO))
        {
            if (sameDirection(AC, AO))
            {
                simplex.clear();
                simplex.append(C);
                simplex.append(A);
                direction = (AC ^ AO) ^ AC;
            }
            else if (sameDirection(AB, AO))
            {
                simplex.clear();
                simplex.append(B);
                simplex.append(A);
                direction = (AB ^ AO) ^ AB;
            }
            else
            {
                simplex.clear();
                simplex.append(A);
                direction = AO;
            }
        }
        else if (sameDirection(AB ^ ABC, AO))
        {
            if (sameDirection(AB, AO))
            {
                simplex.clear();
                simplex.append(B);
                simplex.append(A);
                direction = (AB ^ AO) ^ AB;
            }
            else
            {
                simplex.clear();
                simplex.append(A);
                direction = AO;
            }
        }
        else if (sameDirection(ABC, AO))
        {
            direction = ABC;
        }
        else
        {
            simplex.clear();
            simplex.append(B);
            simplex.append(C);
            simplex.append(A);
            direction = -ABC;
        }
        return magSqr(direction) <= VSMALL;
    }

    if (simplex.size() == 4)
    {
        const vector A(simplex[3]);
        const vector B(simplex[2]);
        const vector C(simplex[1]);
        const vector D(simplex[0]);
        const vector AO(-A);

        vector faceB(B);
        vector faceC(C);
        vector normal((faceB - A) ^ (faceC - A));
        if ((normal & (D - A)) > 0)
        {
            const vector swapPoint(faceB);
            faceB = faceC;
            faceC = swapPoint;
            normal = -normal;
        }
        if (sameDirection(normal, AO))
        {
            simplex.clear();
            simplex.append(faceC);
            simplex.append(faceB);
            simplex.append(A);
            direction = normal;
            return magSqr(direction) <= VSMALL;
        }

        faceB = C;
        faceC = D;
        normal = (faceB - A) ^ (faceC - A);
        if ((normal & (B - A)) > 0)
        {
            const vector swapPoint(faceB);
            faceB = faceC;
            faceC = swapPoint;
            normal = -normal;
        }
        if (sameDirection(normal, AO))
        {
            simplex.clear();
            simplex.append(faceC);
            simplex.append(faceB);
            simplex.append(A);
            direction = normal;
            return magSqr(direction) <= VSMALL;
        }

        faceB = D;
        faceC = B;
        normal = (faceB - A) ^ (faceC - A);
        if ((normal & (C - A)) > 0)
        {
            const vector swapPoint(faceB);
            faceB = faceC;
            faceC = swapPoint;
            normal = -normal;
        }
        if (sameDirection(normal, AO))
        {
            simplex.clear();
            simplex.append(faceC);
            simplex.append(faceB);
            simplex.append(A);
            direction = normal;
            return magSqr(direction) <= VSMALL;
        }

        return true;
    }

    return false;
}

//---------------------------------------------------------------------------//
bool addModelRepeatRandomPosition::convexBodiesIntersect
(
    const pointField& first,
    const pointField& second,
    const vector& initialDirection,
    const scalar lengthScale
)
{
    vector direction(initialDirection);
    if (magSqr(direction) <= VSMALL)
    {
        direction = vector(1, 0, 0);
    }

    DynamicList<vector> simplex(4);
    simplex.append(minkowskiSupport(first, second, direction));
    direction = -simplex[0];

    const scalar projectionTolerance =
        SMALL*Foam::max(sqr(lengthScale), scalar(1));

    for (label iteration = 0; iteration < 64; ++iteration)
    {
        if (magSqr(direction) <= VSMALL)
        {
            return true;
        }

        const vector support
        (
            minkowskiSupport(first, second, direction)
        );
        if ((support & direction) < -projectionTolerance)
        {
            return false;
        }

        simplex.append(support);
        if (updateSimplex(simplex, direction))
        {
            return true;
        }
    }

    // Degenerate/non-converged poses are rejected conservatively.
    return true;
}

//---------------------------------------------------------------------------//
bool addModelRepeatRandomPosition::moveBodyToContinuousPosition()
{
    const vector centre(geomModel_->getCoM());
    const boundBox bodyBounds(geomModel_->getBounds());
    const vector validDirs((geometricD + vector::one)/2);
    vector translation(vector::zero);

    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        if (validDirs[cmpt] < 0.5)
        {
            translation[cmpt] =
                0.5*(insertionMinBound_[cmpt] + insertionMaxBound_[cmpt])
              - centre[cmpt];
            continue;
        }

        const label minimumFace = 2*cmpt;
        const label maximumFace = minimumFace + 1;

        const scalar minimumTranslation = insertionMinBound_[cmpt] -
        (
            containingFaces_[minimumFace]
          ? bodyBounds.min()[cmpt]
          : centre[cmpt]
        );
        const scalar maximumTranslation = insertionMaxBound_[cmpt] -
        (
            containingFaces_[maximumFace]
          ? bodyBounds.max()[cmpt]
          : centre[cmpt]
        );

        if (maximumTranslation < minimumTranslation)
        {
            return false;
        }

        translation[cmpt] = minimumTranslation
          + randGen_.scalar01()
           *(maximumTranslation - minimumTranslation);
    }

    geomModel_->bodyMovePoints(translation);
    return true;
}

//---------------------------------------------------------------------------//
std::shared_ptr<geomModel>
addModelRepeatRandomPosition::addBodyContinuousExact
(
    const volScalarField& body,
    PtrList<immersedBody>& immersedBodies
)
{
    resetExactPlacementBudget();
    bodyAdded_ = false;
    rebuildExactBodyBins(immersedBodies);

    const label batchEnd = Foam::min
    (
        maxCandidateAttemptsPerTimeStep_,
        attemptsOnTimeLevel_ + candidateBatchSize_
    );

    while (attemptsOnTimeLevel_ < batchEnd)
    {
        attemptsOnTimeLevel_++;
        bodyAdditionAttemptCounter_++;
        geomModel_->resetBody();

        if (rotateParticles_)
        {
            scalar rotAngle(0);

            if (rotationMode_ == "uniformRandomRotation")
            {
                const scalar u1 = randGen_.scalar01();
                const scalar u2 = randGen_.scalar01();
                const scalar u3 = randGen_.scalar01();
                const scalar twoPi =
                    2.0*Foam::constant::mathematical::pi;
                const scalar rootOneMinusU1 = Foam::sqrt(1.0 - u1);
                const scalar rootU1 = Foam::sqrt(u1);

                vector quaternionVector
                (
                    rootOneMinusU1*Foam::sin(twoPi*u2),
                    rootOneMinusU1*Foam::cos(twoPi*u2),
                    rootU1*Foam::sin(twoPi*u3)
                );
                scalar quaternionScalar = rootU1*Foam::cos(twoPi*u3);

                if (quaternionScalar < 0.0)
                {
                    quaternionScalar = -quaternionScalar;
                    quaternionVector = -quaternionVector;
                }

                const scalar sinHalfAngle = mag(quaternionVector);
                if (sinHalfAngle > SMALL)
                {
                    axisOfRot_ = quaternionVector/sinHalfAngle;
                    rotAngle =
                        2.0*Foam::atan2
                        (
                            sinHalfAngle,
                            quaternionScalar
                        );
                }
                else
                {
                    axisOfRot_ = vector(1.0, 0.0, 0.0);
                    rotAngle = 0.0;
                }
            }
            else
            {
                rotAngle =
                    (2.0*randGen_.scalar01() - 1.0)
                   *Foam::constant::mathematical::pi;

                if (randomAxis_)
                {
                    axisOfRot_ = vector
                    (
                        2.0*randGen_.scalar01() - 1.0,
                        2.0*randGen_.scalar01() - 1.0,
                        2.0*randGen_.scalar01() - 1.0
                    );
                    if (mag(axisOfRot_) <= VSMALL)
                    {
                        axisOfRot_ = vector(1.0, 0.0, 0.0);
                    }
                    else
                    {
                        axisOfRot_ /= mag(axisOfRot_);
                    }
                }
            }

            geomModel_->bodyRotatePoints(rotAngle, axisOfRot_);
        }

        if (!moveBodyToContinuousPosition())
        {
            continue;
        }

        if (exactCandidateInContact(immersedBodies))
        {
            continue;
        }

        // Keep the established framework contact test as the final guard for
        // physical walls and all existing body types.  The candidate is at
        // its actual final size: no temporary 1.02 enlargement is applied.
        volScalarField helpBodyField = body;
        geomModel_->createImmersedBody
        (
            helpBodyField,
            octreeField_,
            cellPoints_
        );

        bool canAddBody = !isBodyInContact(immersedBodies);
        reduce(canAddBody, andOp<bool>());

        if (canAddBody)
        {
            bodyAdded_ = true;

            // The exact branch returns before the legacy success block below,
            // so preserve repeatRandomPosition's timeBased accounting here.
            if (timeBased_)
            {
                addedOnTimeLevel_++;
                if (addedOnTimeLevel_ == partPerAdd_)
                {
                    useNTimes_--;
                    reapeatedAddition_ = false;
                }
            }

            InfoH << addModelSummary_Info
                << "[repeatRandomPosition] accepted full-size continuous "
                << "pose after " << attemptsOnTimeLevel_
                << '/' << maxCandidateAttemptsPerTimeStep_
                << " proposals on time index " << attemptTimeIndex_ << endl;
            return geomModel_->getCopy();
        }
    }

    if (attemptsOnTimeLevel_ >= maxCandidateAttemptsPerTimeStep_)
    {
        exactPlacementExhausted_ = true;
        InfoH << addModelSummary_Info
            << "[repeatRandomPosition] no feasible full-size pose in "
            << maxCandidateAttemptsPerTimeStep_
            << " proposals on time index " << attemptTimeIndex_
            << "; insertion is deferred to the next time step" << endl;
    }

    return std::shared_ptr<geomModel>();
}
//---------------------------------------------------------------------------//
std::shared_ptr<geomModel> addModelRepeatRandomPosition::addBody
(
    const volScalarField& body,
    PtrList<immersedBody>& immersedBodies
)
{
    if (placementStrategy_ == "continuousExact")
    {
        return addBodyContinuousExact(body, immersedBodies);
    }

    geomModel_->resetBody();

    bodyAdditionAttemptCounter_++;

    // rotate
    if (rotateParticles_)
    {
        scalar rotAngle(0);
        if (rotationMode_ == "uniformRandomRotation")
        {
            // Shoemake's construction: three uniform random values produce a
            // quaternion distributed uniformly over SO(3).
            const scalar u1 = randGen_.scalar01();
            const scalar u2 = randGen_.scalar01();
            const scalar u3 = randGen_.scalar01();

            const scalar twoPi =
                2.0*Foam::constant::mathematical::pi;
            const scalar rootOneMinusU1 = Foam::sqrt(1.0 - u1);
            const scalar rootU1 = Foam::sqrt(u1);

            vector quaternionVector
            (
                rootOneMinusU1*Foam::sin(twoPi*u2),
                rootOneMinusU1*Foam::cos(twoPi*u2),
                rootU1*Foam::sin(twoPi*u3)
            );
            scalar quaternionScalar = rootU1*Foam::cos(twoPi*u3);

            // q and -q represent the same rotation. Keep the angle in [0, pi].
            if (quaternionScalar < 0.0)
            {
                quaternionScalar = -quaternionScalar;
                quaternionVector = -quaternionVector;
            }

            const scalar sinHalfAngle = mag(quaternionVector);
            if (sinHalfAngle > SMALL)
            {
                axisOfRot_ = quaternionVector/sinHalfAngle;
                rotAngle =
                    2.0*Foam::atan2(sinHalfAngle, quaternionScalar);
            }
            else
            {
                axisOfRot_ = vector(1.0, 0.0, 0.0);
                rotAngle = 0.0;
            }
        }
        else
        {
            rotAngle = returnRandomAngle();
            if (randomAxis_)
            {
                axisOfRot_ = returnRandomRotationAxis();
            }
        }
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "Will rotate by " << rotAngle
            << " rad around axis " << axisOfRot_ << endl;

        geomModel_->bodyRotatePoints(rotAngle,axisOfRot_);
        //~ CoM = gSum(bodySurfMesh.coordinates())/bodySurfMesh.size();
    }

    // scale
    if (scaleApplication_ or scaleRandomApplication_)
    {
        if (scaleRandomApplication_){scaleStep_ = returnRandomScale();}
        geomModel_->bodyScalePoints(scaleStep_);
    }

    geomModel_->bodyScalePoints(1.02);

    vector CoM(geomModel_->getCoM());
    point bBoxCenter = cellZoneBounds_.midpoint();
    geomModel_->bodyMovePoints(bBoxCenter - CoM);

    // translate

    vector randomTrans = geomModel_->addModelReturnRandomPosition(allActiveCellsInMesh_,cellZoneBounds_,randGen_);
    geomModel_->bodyMovePoints(randomTrans);

    // check if the body can be added
    volScalarField helpBodyField_ = body;
    geomModel_->createImmersedBody(
        helpBodyField_,
        octreeField_,
        cellPoints_
    );

    bool canAddBodyI = !isBodyInContact(immersedBodies);

    geomModel_->bodyScalePoints(1.0/1.02);

    reduce(canAddBodyI, andOp<bool>());
    bodyAdded_ = (canAddBodyI);

    if(!bodyAdded_)
	{
		scaleCorrectionCounter_++;
	}

	if(bodyAdded_)
	{
		if(timeBased_)
		{
			InfoH << addModel_Info << "-- addModelMessage-- "
                << "addedOnTimeLevel:  " << addedOnTimeLevel_<< endl;
			addedOnTimeLevel_++;
			InfoH << "-- addModelMessage-- " << "bodyAdded: "
                << bodyAdded_ << " addedOnTimeLevel:  " << addedOnTimeLevel_
                << " useNTimes: " << useNTimes_<<  endl;
			if(addedOnTimeLevel_ == partPerAdd_)
			{
				useNTimes_--;
				InfoH << "-- addModelMessage-- "
                    <<" useNTimes: " << useNTimes_<<  endl;
				reapeatedAddition_ = false;
			}
		}

		scaleCorrectionCounter_ = 0;

	}

	InfoH << addModel_Info << "-- addModelMessage-- "
        << "bodyAdditionAttemptNr  : " << bodyAdditionAttemptCounter_<< endl;
	InfoH << "-- addModelMessage-- " << "sameScaleAttempts      : "
        << scaleCorrectionCounter_<< endl;

	if(scaleCorrectionCounter_ > nTriesBeforeScaling_ && scaleParticles_)
	{
		scaleApplication_ = true;
		rescaleRequirement_ = true;
		scalingFactor_++;
		scaleStep_ = pow(scaleStep_,scalingFactor_);
		if(scaleStep_<minScaleFit_)
		{
			scaleStep_=minScaleFit_;
		}
		scaleCorrectionCounter_ = 0;
	}

    return geomModel_->getCopy();
}
// MODEL SPECIFIC FUNCTIONS==================================================//
//---------------------------------------------------------------------------//
void addModelRepeatRandomPosition::initializeCellZone()
{

	label zoneID = mesh_.cellZones().findZoneID(zoneName_);
	InfoH << addModel_Info << "-- addModelMessage-- "
        << "label of the cellZone " << zoneID << endl;

	const labelList& cellZoneCells = mesh_.cellZones()[zoneID];
    cellsInBoundBox_[Pstream::myProcNo()] = cellZoneCells;

	const pointField& cp = mesh_.C();
	const pointField fCp(cp,cellsInBoundBox_[Pstream::myProcNo()]);
	cellZonePoints_[Pstream::myProcNo()] = fCp;

	updateCellZoneBoundBox();
}
//---------------------------------------------------------------------------//
void addModelRepeatRandomPosition::updateCellZoneBoundBox()
{
		boundBox cellZoneBounds(cellZonePoints_[Pstream::myProcNo()]);

        reduce(cellZoneBounds.min(), minOp<vector>());
        reduce(cellZoneBounds.max(), maxOp<vector>());

        if (Pstream::myProcNo() == 0)
        {
            minBound_ = cellZoneBounds_.min();
            maxBound_ = cellZoneBounds_.max();
            cellZoneBounds_ = boundBox(minBound_,maxBound_);
        }
}
//---------------------------------------------------------------------------//
void addModelRepeatRandomPosition::initializeBoundBox()
{
    octreeField_ *= 0;
    List<DynamicLabelList> bBoxCells(Pstream::nProcs());

    bool isInsideBB(false);
    labelList nextToCheck(1,0);
    label iterCount(0);label iterMax(mesh_.nCells());
    while ((nextToCheck.size() > 0 or not isInsideBB) and iterCount < iterMax)
    {
        iterCount++;
        DynamicLabelList auxToCheck;

        forAll (nextToCheck,cellToCheck)
        {
            auxToCheck.append(
                getBBoxCellsByOctTree(
                    nextToCheck[cellToCheck],
                    isInsideBB,
                    minBound_,maxBound_,bBoxCells
                )
            );
        }
        nextToCheck = auxToCheck;
    }

    cellsInBoundBox_[Pstream::myProcNo()] = bBoxCells[Pstream::myProcNo()];

    cellZoneBounds_ = boundBox(minBound_,maxBound_);
}
//---------------------------------------------------------------------------//
void addModelRepeatRandomPosition::recreateBoundBox()
{
    octreeField_ = Field<label>(mesh_.nCells(), 0);
    cellsInBoundBox_[Pstream::myProcNo()].clear();
    List<DynamicLabelList> bBoxCells(Pstream::nProcs());

    bool isInsideBB(false);
    labelList nextToCheck(1,0);
    label iterCount(0);label iterMax(mesh_.nCells());
    while ((nextToCheck.size() > 0 or not isInsideBB) and iterCount < iterMax)
    {
        iterCount++;
        DynamicLabelList auxToCheck;

        forAll (nextToCheck,cellToCheck)
        {
            auxToCheck.append(
                getBBoxCellsByOctTree(
                    nextToCheck[cellToCheck],
                    isInsideBB,
                    minBound_,maxBound_,bBoxCells
                )
            );
        }
        nextToCheck = auxToCheck;
    }

    cellsInBoundBox_[Pstream::myProcNo()] = bBoxCells[Pstream::myProcNo()];

    InfoH << addModel_Info << "-- addModelMessage-- "
        << "recreated boundBox size "
        << cellsInBoundBox_[Pstream::myProcNo()].size() << endl;

    cellZoneBounds_ = boundBox(minBound_,maxBound_);
}
//---------------------------------------------------------------------------//
labelList addModelRepeatRandomPosition::getBBoxCellsByOctTree
(
    label cellToCheck,
    bool& insideBB,
    vector& bBoxMin,
    vector& bBoxMax,
    List<DynamicLabelList>& bBoxCells
)
{
    labelList retList;

    if (octreeField_[cellToCheck] ==0)
    {
        octreeField_[cellToCheck] = 1;
        vector cCenter = mesh_.C()[cellToCheck];
        label   partCheck(0);
        forAll (bBoxMin,vecI)
        {
            if (cCenter[vecI] >= bBoxMin[vecI] and cCenter[vecI] <= bBoxMax[vecI])
            {
                partCheck++;
            }
        }
        bool cellInside = (partCheck == 3) ? true : false;
        if (cellInside)
        {
            bBoxCells[Pstream::myProcNo()].append(cellToCheck);
            insideBB = true;
        }
        if (not insideBB or cellInside)
        {
            retList = mesh_.cellCells()[cellToCheck];
        }
    }
    return retList;
}
//---------------------------------------------------------------------------//
scalar addModelRepeatRandomPosition::checkLambdaFraction(const volScalarField& body)
{
	scalarList lambdaIntegrate(Pstream::nProcs());
    scalarList volumeIntegrate(Pstream::nProcs());
	scalar lambdaFraction(0);
    forAll (lambdaIntegrate,k)
    {
        lambdaIntegrate[k] = 0;
        volumeIntegrate[k] = 0;
    }
	forAll (cellsInBoundBox_[Pstream::myProcNo()],k)
	{
		label cell = cellsInBoundBox_[Pstream::myProcNo()][k];
		lambdaIntegrate[Pstream::myProcNo()] += mesh_.V()[cell]*body[cell];
		volumeIntegrate[Pstream::myProcNo()] += mesh_.V()[cell];
	}
	lambdaFraction = gSum(lambdaIntegrate)/gSum(volumeIntegrate);
	InfoH << addModel_Info << "-- addModelMessage-- "
        << "lambda fraction in controlled region: " << lambdaFraction<< endl;
	return lambdaFraction;
}
//---------------------------------------------------------------------------//
scalar addModelRepeatRandomPosition::returnRandomAngle()
{
    scalar ranNum = 2.0*randGen_.scalar01() - 1.0;
    scalar angle  = ranNum*Foam::constant::mathematical::pi;
	return angle;
}
//---------------------------------------------------------------------------//
scalar addModelRepeatRandomPosition::returnRandomScale()
{
	scalar ranNum       = randGen_.scalar01();
	scalar scaleDiff    = maxScale_ - minScale_;
    scalar scaleFactor  = minScale_ + ranNum*scaleDiff;
	InfoH << addModel_Info << "-- addModelMessage-- "
        <<"random scaleFactor " << scaleFactor <<endl;
	return scaleFactor;
}
//---------------------------------------------------------------------------//
vector addModelRepeatRandomPosition::returnRandomRotationAxis()
{
	vector  axisOfRotation(vector::zero);
	scalar ranNum = 0;

	for (int i=0;i<3;i++)
	{
		ranNum = randGen_.scalar01();
		axisOfRotation[i] = ranNum;
	}

	axisOfRotation /=mag(axisOfRotation);
	return axisOfRotation;
}
//---------------------------------------------------------------------------//
