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
#include "addModelOnceScatter.H"
#include "meshSearch.H"

#include <cmath>

using namespace Foam;

//---------------------------------------------------------------------------//
addModelOnceScatter::addModelOnceScatter
(
    const dictionary& addModelDict,
    const Foam::fvMesh& mesh,
    const bool startTime0,
    std::unique_ptr<geomModel> bodyGeomModel,
    List<labelList>& cellPoints
)
:
addModel(mesh, std::move(bodyGeomModel), cellPoints),
addModelDict_(addModelDict),
addMode_(word(addModelDict_.lookup("addModel"))),
bodyAdded_(false),
finishedAddition_(false),

coeffsDict_(addModelDict_.subDict(addMode_+"Coeffs")),

addDomain_(word(coeffsDict_.lookup("addDomain"))),
scalingMode_(word(coeffsDict_.lookup("scalingMode"))),
rotationMode_(word(coeffsDict_.lookup("rotationMode"))),
addModeI_(word(coeffsDict_.lookup("addMode"))),

addDomainCoeffs_(coeffsDict_.subDict(addDomain_ + "Coeffs")),
scalingModeCoeffs_(coeffsDict_.subDict(scalingMode_ + "Coeffs")),
rotationModeCoeffs_(coeffsDict_.subDict(rotationMode_ + "Coeffs")),
addModeICoeffs_(coeffsDict_.subDict(addModeI_ + "Coeffs")),

partPerAdd_(0),
fieldValue_(0),
addedOnTimeLevel_(0),
startTimeValue_(mesh_.time().value()),

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
polygonPrismActive_(false),
octreeField_(mesh_.nCells(), 0),
multiBody_(false),
fieldBased_(false),
geometricVolumeBased_(false),
fieldCurrentValue_(0),
targetGeometricFraction_(0),
regionVolume_(0),
particleVolume_(0),
fixedScale_(1),
clearanceScale_(1),
targetBodyCount_(0),
acceptedBodyCount_(0),
totalPackingAttempts_(0),
reportEveryAccepted_(1),
packingPlanReady_(false),
stagnationWindow_(20000),
localRepackAttempts_(200000),
minLocalRepackBodies_(4),
maxLocalRepackBodies_(16),
maxPlannerRestarts_(3),
localRepackCount_(0),
plannerRestartCount_(0),
packingStrategy_("continuousRandomBacktracking"),
uniformRandomRotation_(false),
packingIsSphere_(false),
sourcePoints_(),
sourceCentre_(vector::zero),
sourceSphereRadius_(0),
packingRadius_(0),
packingBinSize_(0),
nPackingBinsX_(1),
nPackingBinsY_(1),
nPackingBinsZ_(1),
acceptedPackingPoints_(),
plannedActualPoints_(),
acceptedPackingBounds_(),
acceptedPackingCentres_(),
acceptedPackingRadii_(),
packingBins_(),
polygonPrismGeometry_(),
allActiveCellsInMesh_(true),
randGen_
(
    coeffsDict_.lookupOrDefault<label>
    (
        "randomSeed",
        addModeI_ == "geometricVolumeBased"
      ? label(12345)
      : label(clock::getTime())
    )
)
{
    if(!startTime0)
    {
        finishedAddition_ = true;
    }
	init();
}

addModelOnceScatter::~addModelOnceScatter()
{
}

//---------------------------------------------------------------------------//
void addModelOnceScatter::init()
{
    // set sizes to necessary datatypes
    cellsInBoundBox_.setSize(Pstream::nProcs());
    cellZonePoints_.setSize(Pstream::nProcs());


	if (addModeI_ == "multiBody")
	{
		partPerAdd_ = (readLabel(addModeICoeffs_.lookup("partPerAdd")));
        multiBody_  = true;
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "addModel will initialize simulation with "
            << partPerAdd_ << " randomly scattered bodies" << endl;
	}
	else if (addModeI_ == "fieldBased")
	{
		fieldValue_ = (readScalar(addModeICoeffs_.lookup("fieldValue")));
        fieldBased_ = true;
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "addModel will fill domain up to preset volume fraction" << endl;
		InfoH << "-- addModelMessage-- "
            << "preset volume fraction: " << fieldValue_ << endl;
	}
    else if (addModeI_ == "geometricVolumeBased")
    {
        targetGeometricFraction_ =
            readScalar(addModeICoeffs_.lookup("fieldValue"));
        geometricVolumeBased_ = true;

        InfoH << addModel_Info << "-- addModelMessage-- "
            << "addModel will place the nearest integer number of bodies "
            << "for geometric volume fraction "
            << targetGeometricFraction_ << endl;
    }
    else
    {
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "notImplemented, will crash" << endl;
    }

	if (addDomain_ == "cellZone")
	{
        if (geometricVolumeBased_)
        {
            FatalErrorInFunction
                << "addMode geometricVolumeBased currently requires "
                << "addDomain boundBox. Exact containment in an arbitrary "
                << "cellZone is not defined." << nl
                << exit(FatalError);
        }
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
	else if (addDomain_ == "polygonPrism")
	{
        polygonPrismGeometry_.reset
        (
            new polygonPrismGeometry(addDomainCoeffs_)
        );
        polygonPrismActive_ = true;
        initializePolygonPrism();
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "convex polygonPrism based addition zone" << endl;
        InfoH << "-- addModelMessage-- polygon area: "
            << polygonPrismGeometry_->area()
            << ", prism volume: " << polygonPrismGeometry_->volume()
            << endl;
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

    if (polygonPrismActive_)
    {
        scalar localZoneVolume = 0;

        forAll(cellsInBoundBox_[Pstream::myProcNo()], cellI)
        {
            localZoneVolume +=
                mesh_.V()
                [cellsInBoundBox_[Pstream::myProcNo()][cellI]];
        }

        scalar zoneVolume = localZoneVolume;
        reduce(zoneVolume, sumOp<scalar>());
        const scalar exactVolume = polygonPrismGeometry_->volume();
        const scalar relativeVolumeError =
            mag(zoneVolume - exactVolume)/(exactVolume + VSMALL);

        allActiveCellsInMesh_ = relativeVolumeError <= 1e-5;

        InfoH << addModel_Info << "-- addModelMessage-- "
            << "polygonPrism selected-cell volume: " << zoneVolume
            << ", exact volume: " << exactVolume
            << ", relative difference: " << relativeVolumeError
            << endl;

        if (!allActiveCellsInMesh_)
        {
            InfoH << addModel_Info << "-- addModelMessage-- "
                << "polygonPrism is not mesh-conformal at cell-center "
                << "resolution; fieldBased volume fraction is approximate."
                << endl;
        }
    }
    else
    {
        // Preserve the existing cellZone/boundBox diagnostic unchanged.
        scalarList procZoneVols(Pstream::nProcs());
        procZoneVols[Pstream::myProcNo()] = 0;
        forAll (cellsInBoundBox_[Pstream::myProcNo()],cellI)
        {
            procZoneVols[Pstream::myProcNo()]
                += mesh_.V()
                [cellsInBoundBox_[Pstream::myProcNo()][cellI]];
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
    }

    if (coeffsDict_.found("randomSeed"))
    {
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "using configured randomSeed "
            << readLabel(coeffsDict_.lookup("randomSeed")) << endl;
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
		InfoH << "-- addModelMessage-- " << "nTriesBeforeDownScaling: "
            << nTriesBeforeScaling_ << endl;
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
        uniformRandomRotation_ = true;
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "source STL orientation will be sampled uniformly in 3D "
            << "upon addition" << endl;
	}
	else if (rotationMode_ == "fixedAxisRandomRotation")
	{
		axisOfRot_       = (rotationModeCoeffs_.lookup("axis"));
        InfoH << addModel_Info << "-- addModelMessage-- "
            << "source STL will be rotated by a random angle around a fixed axis upon addition" << endl;
		InfoH << "-- addModelMessage-- " << "set rotation axis: " << axisOfRot_ << endl;
		rotateParticles_ = true;
		randomAxis_      = false;
	}
    else
    {
		InfoH << addModel_Info << "-- addModelMessage-- "
            << "notImplemented, will crash" << endl;
	}

    if (geometricVolumeBased_)
    {
        initializeGeometricVolumePacking();
    }
}

//---------------------------------------------------------------------------//
bool addModelOnceScatter::shouldAddBody(const volScalarField& body)
{
    if(startTimeValue_ != mesh_.time().value())
    {
        return false;
    }

    if (geometricVolumeBased_)
    {
        return !finishedAddition_ && acceptedBodyCount_ < targetBodyCount_;
    }

    if (!finishedAddition_)
    {
        if (multiBody_)
        {

            InfoH << addModel_Info << "-- addModelMessage-- "
                << "Number of added bodies: " << addedOnTimeLevel_ << endl;

            if (partPerAdd_ < addedOnTimeLevel_)
            {
                finishedAddition_ = false;
            }
            else
            {
                finishedAddition_ = true;
            }
        }

        if (fieldBased_)
        {
            scalar currentLambdaFrac(checkLambdaFraction(body));
            if (currentLambdaFrac < fieldValue_ )
            {
                InfoH << addModel_Info << "-- addModelMessage-- "
                    << "Current lambda fraction = " << currentLambdaFrac
                    << " < then preset lambda fraction = "
                    << fieldValue_ << endl;
                finishedAddition_ = false;
            }
            else
            {
                finishedAddition_ = true;
            }
        }
    }

    return !finishedAddition_;

}
//---------------------------------------------------------------------------//
std::shared_ptr<geomModel> addModelOnceScatter::addBody
(
    const   volScalarField& body,
    PtrList<immersedBody>& immersedBodies
)
{
    if (geometricVolumeBased_)
    {
        return addGeometricVolumeBody(immersedBodies);
    }

    geomModel_->resetBody();

    bodyAdditionAttemptCounter_++;

    // rotate
    if (rotateParticles_)
    {
        scalar rotAngle(0);
        if (rotationMode_ == "uniformRandomRotation")
        {
            // Shoemake's construction: three uniform random values produce
            // a quaternion distributed uniformly over SO(3).
            const scalar u1 = polygonPrismActive_
              ? returnSynchronizedRandom01()
              : randGen_.scalar01();
            const scalar u2 = polygonPrismActive_
              ? returnSynchronizedRandom01()
              : randGen_.scalar01();
            const scalar u3 = polygonPrismActive_
              ? returnSynchronizedRandom01()
              : randGen_.scalar01();

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
    }

    // scale
    if (scaleApplication_ or scaleRandomApplication_)
    {
        if (scaleRandomApplication_){scaleStep_ = returnRandomScale();}
        geomModel_->bodyScalePoints(scaleStep_);
    }

    geomModel_->bodyScalePoints(1.02);

    bool validPosition = true;

    if (polygonPrismActive_)
    {
        scalar isotropicRadius = 0;
        const pointField relativeSupportPoints =
            relativeBodySupportPoints(isotropicRadius);

        scalarField unitRandomValues(4, 0);
        forAll(unitRandomValues, randomI)
        {
            unitRandomValues[randomI] =
                returnSynchronizedRandom01();
        }

        const point sampledCenter =
            polygonPrismGeometry_->sampleCenter
            (
                relativeSupportPoints,
                isotropicRadius,
                unitRandomValues,
                validPosition
            );

        if (validPosition)
        {
            geomModel_->bodyMovePoints
            (
                sampledCenter - geomModel_->getCoM()
            );
        }
        else
        {
            InfoH << addModel_Info << "-- addModelMessage-- "
                << "rotated/scaled body does not fit in polygonPrism"
                << endl;
        }
    }
    else
    {
        // Preserve the existing cellZone/boundBox placement behavior.
        vector CoM(geomModel_->getCoM());
        point bBoxCenter = cellZoneBounds_.midpoint();
        geomModel_->bodyMovePoints(bBoxCenter - CoM);

        vector randomTrans = geomModel_->addModelReturnRandomPosition
        (
            allActiveCellsInMesh_,
            cellZoneBounds_,
            randGen_
        );
        geomModel_->bodyMovePoints(randomTrans);
    }

    // check if the body can be added
    bool canAddBodyI = false;

    if (validPosition)
    {
        volScalarField helpBodyField_ = body;
        geomModel_->createImmersedBody
        (
            helpBodyField_,
            octreeField_,
            cellPoints_
        );

        canAddBodyI = !isBodyInContact(immersedBodies);
    }

    geomModel_->bodyScalePoints(1.0/1.02);

    reduce(canAddBodyI, andOp<bool>());
    bodyAdded_ = (canAddBodyI);

    if(!bodyAdded_)
	{
		scaleCorrectionCounter_++;
	}

	if(bodyAdded_)
	{
		if(multiBody_)
		{
			InfoH << addModel_Info << "-- addModelMessage-- "
                << "addedOnTimeLevel:  " << addedOnTimeLevel_<< endl;
			addedOnTimeLevel_++;
		}

		scaleCorrectionCounter_ = 0;

	}

	InfoH << addModel_Info << "-- addModelMessage-- "
        << "bodyAdditionAttemptNr  : " << bodyAdditionAttemptCounter_<< endl;
	InfoH << "-- addModelMessage-- "
        << "sameScaleAttempts      : " << scaleCorrectionCounter_<< endl;

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

//---------------------------------------------------------------------------//
void addModelOnceScatter::initializeGeometricVolumePacking()
{
    if (!boundBoxActive_)
    {
        FatalErrorInFunction
            << "addMode geometricVolumeBased currently requires "
            << "addDomain boundBox.  Exact containment in an arbitrary "
            << "cellZone is not defined." << nl
            << exit(FatalError);
    }

    if
    (
        targetGeometricFraction_ < 0
     || targetGeometricFraction_ > 1
    )
    {
        FatalErrorInFunction
            << "geometricVolumeBased fieldValue must be in [0,1], got "
            << targetGeometricFraction_ << nl
            << exit(FatalError);
    }

    if (scalingMode_ == "noScaling")
    {
        fixedScale_ = 1.0;
    }
    else if (scalingMode_ == "randomScaling")
    {
        const scalar scaleTolerance =
            SMALL*Foam::max
            (
                scalar(1),
                Foam::max(mag(minScale_), mag(maxScale_))
            );

        if (mag(maxScale_ - minScale_) > scaleTolerance)
        {
            FatalErrorInFunction
                << "geometricVolumeBased needs one known particle volume. "
                << "For randomScaling, minScale and maxScale must be equal; "
                << "got " << minScale_ << " and " << maxScale_ << nl
                << exit(FatalError);
        }
        fixedScale_ = 0.5*(minScale_ + maxScale_);
    }
    else
    {
        FatalErrorInFunction
            << "geometricVolumeBased supports noScaling or fixed "
            << "randomScaling (minScale == maxScale), not "
            << scalingMode_ << nl
            << exit(FatalError);
    }

    if (fixedScale_ <= 0)
    {
        FatalErrorInFunction
            << "Particle scale must be positive, got " << fixedScale_ << nl
            << exit(FatalError);
    }

    clearanceScale_ =
        coeffsDict_.lookupOrDefault<scalar>("clearanceScale", 1.0);
    stagnationWindow_ = coeffsDict_.lookupOrDefault<label>
    (
        "stagnationWindow",
        20000
    );
    localRepackAttempts_ = coeffsDict_.lookupOrDefault<label>
    (
        "localRepackAttempts",
        200000
    );
    minLocalRepackBodies_ = coeffsDict_.lookupOrDefault<label>
    (
        "minLocalRepackBodies",
        4
    );
    maxLocalRepackBodies_ = coeffsDict_.lookupOrDefault<label>
    (
        "maxLocalRepackBodies",
        16
    );
    maxPlannerRestarts_ = coeffsDict_.lookupOrDefault<label>
    (
        "maxPlannerRestarts",
        3
    );
    reportEveryAccepted_ = Foam::max
    (
        label(1),
        coeffsDict_.lookupOrDefault<label>("reportEveryAccepted", 100)
    );

    if (mag(clearanceScale_ - 1.0) > SMALL)
    {
        FatalErrorInFunction
            << "continuousRandomBacktracking requires clearanceScale 1.0; "
            << "got "
            << clearanceScale_ << nl
            << exit(FatalError);
    }
    if
    (
        stagnationWindow_ < 1
     || localRepackAttempts_ < 1
     || minLocalRepackBodies_ < 1
     || maxLocalRepackBodies_ < minLocalRepackBodies_
     || maxPlannerRestarts_ < 0
    )
    {
        FatalErrorInFunction
            << "Invalid continuousRandomBacktracking controls: "
            << "stagnationWindow=" << stagnationWindow_
            << ", localRepackAttempts=" << localRepackAttempts_
            << ", minLocalRepackBodies=" << minLocalRepackBodies_
            << ", maxLocalRepackBodies=" << maxLocalRepackBodies_
            << ", maxPlannerRestarts=" << maxPlannerRestarts_
            << nl << exit(FatalError);
    }

    const contactType geometryType = geomModel_->getcType();
    packingIsSphere_ = geometryType == sphere;
    if (!packingIsSphere_ && geometryType != convex)
    {
        FatalErrorInFunction
            << "Fast geometricVolumeBased packing supports convex and "
            << "sphere geometry only." << nl
            << exit(FatalError);
    }

    sourcePoints_ = geomModel_->getBodyPoints();
    sourceCentre_ = geomModel_->getCoM();
    if (sourcePoints_.empty())
    {
        FatalErrorInFunction
            << "The geometry returned no packing points." << nl
            << exit(FatalError);
    }

    if (packingIsSphere_)
    {
        const vector span(geomModel_->getBounds().span());
        sourceSphereRadius_ = 0.5*Foam::max
        (
            span.x(),
            Foam::max(span.y(), span.z())
        );
        packingRadius_ = sourceSphereRadius_*fixedScale_*clearanceScale_;
    }
    else
    {
        scalar maximumRadiusSquared = 0;
        forAll(sourcePoints_, pointI)
        {
            maximumRadiusSquared = Foam::max
            (
                maximumRadiusSquared,
                magSqr(sourcePoints_[pointI] - sourceCentre_)
            );
        }
        packingRadius_ =
            Foam::sqrt(maximumRadiusSquared)*fixedScale_*clearanceScale_;
    }

    const scalar baseVolume = geomModel_->geometricVolume();
    particleVolume_ = baseVolume*pow(fixedScale_, 3);
    regionVolume_ = cellZoneBounds_.volume();

    if (baseVolume <= VSMALL || particleVolume_ <= VSMALL)
    {
        FatalErrorInFunction
            << "The selected geometry has no positive exact volume "
            << "implementation. Returned base volume: " << baseVolume << nl
            << exit(FatalError);
    }
    if (regionVolume_ <= VSMALL || packingRadius_ <= VSMALL)
    {
        FatalErrorInFunction
            << "Invalid packing region or particle radius. regionVolume="
            << regionVolume_ << ", packingRadius=" << packingRadius_ << nl
            << exit(FatalError);
    }

    const scalar exactBodyCount =
        targetGeometricFraction_*regionVolume_/particleVolume_;

    // Deliberately fixed to nearest-integer rounding; this is not a user
    // option because sub-particle control of volume fraction is impossible.
    targetBodyCount_ = label(std::floor(exactBodyCount + 0.5));

    packingStrategy_ = coeffsDict_.lookupOrDefault<word>
    (
        "packingStrategy",
        "continuousRandomBacktracking"
    );
    if (packingStrategy_ != "continuousRandomBacktracking")
    {
        FatalErrorInFunction
            << "Unknown geometricVolumeBased packingStrategy '"
            << packingStrategy_ << "'. This version supports only "
            << "continuousRandomBacktracking." << nl << exit(FatalError);
    }

    packingBinSize_ = 2.0*packingRadius_;
    const vector regionSpan(cellZoneBounds_.span());
    nPackingBinsX_ = Foam::max
    (
        label(1),
        label(std::ceil(regionSpan.x()/packingBinSize_))
    );
    nPackingBinsY_ = Foam::max
    (
        label(1),
        label(std::ceil(regionSpan.y()/packingBinSize_))
    );
    nPackingBinsZ_ = Foam::max
    (
        label(1),
        label(std::ceil(regionSpan.z()/packingBinSize_))
    );

    InfoH << addModel_Info << "-- addModelMessage-- "
        << "geometricVolumeBased region volume: " << regionVolume_ << endl;
    InfoH << addModel_Info << "-- addModelMessage-- "
        << "scaled particle volume: " << particleVolume_ << endl;
    InfoH << addModel_Info << "-- addModelMessage-- "
        << "nearest target body count: " << targetBodyCount_
        << ", realizable fraction: "
        << targetBodyCount_*particleVolume_/regionVolume_ << endl;
}

//---------------------------------------------------------------------------//
std::shared_ptr<geomModel> addModelOnceScatter::addGeometricVolumeBody
(
    PtrList<immersedBody>& immersedBodies
)
{
    if (immersedBodies.size() != acceptedBodyCount_)
    {
        FatalErrorInFunction
            << "geometricVolumeBased must be the first body model during "
            << "initialization. It found " << immersedBodies.size()
            << " existing bodies but owns " << acceptedBodyCount_
            << " packing entries." << nl
            << exit(FatalError);
    }

    if (!packingPlanReady_)
    {
        prepareContinuousPackingPlan();
    }

    if (acceptedBodyCount_ >= plannedActualPoints_.size())
    {
        FatalErrorInFunction
            << "The continuous packing plan contains only "
            << plannedActualPoints_.size() << " poses, but pose "
            << acceptedBodyCount_ << " was requested." << nl
            << exit(FatalError);
    }

    const pointField& actualPoints =
        plannedActualPoints_[acceptedBodyCount_];

    bodyAdded_ = true;

    if (packingIsSphere_)
    {
        geomModel_->resetBody();
        geomModel_->bodyScalePoints(fixedScale_);
        geomModel_->setBodyPosition(actualPoints);
    }
    else
    {
        geomModel_->setBodyPosition(actualPoints);
    }

    return geomModel_->getCopy();
}

//---------------------------------------------------------------------------//
void addModelOnceScatter::appendPackingPose
(
    const pointField& actualPoints,
    const pointField& collisionPoints,
    const vector& centre,
    const scalar radius,
    const boundBox& bounds
)
{
    const label bodyI = acceptedPackingPoints_.size();

    plannedActualPoints_.append(actualPoints);
    acceptedPackingPoints_.append(collisionPoints);
    acceptedPackingBounds_.append(bounds);
    acceptedPackingCentres_.append(centre);
    acceptedPackingRadii_.append(radius);

    const label key = packingBinKey(centre);
    if (!packingBins_.found(key))
    {
        packingBins_.insert(key, DynamicLabelList());
    }
    packingBins_[key].append(bodyI);
}

//---------------------------------------------------------------------------//
void addModelOnceScatter::rebuildPackingBins()
{
    packingBins_.clear();
    forAll(acceptedPackingCentres_, bodyI)
    {
        const label key = packingBinKey(acceptedPackingCentres_[bodyI]);
        if (!packingBins_.found(key))
        {
            packingBins_.insert(key, DynamicLabelList());
        }
        packingBins_[key].append(bodyI);
    }
}

//---------------------------------------------------------------------------//
bool addModelOnceScatter::tryLocalPackingRebuild
(
    const pointField& triggerActualPoints,
    const pointField& triggerCollisionPoints,
    const vector& triggerCentre,
    const scalar triggerRadius,
    const boundBox& triggerBounds,
    const DynamicLabelList& directBlockers
)
{
    const label oldSize = acceptedPackingPoints_.size();
    if (oldSize == 0 || directBlockers.empty())
    {
        return false;
    }

    // No immersedBody exists yet.  These copies are a complete rollback point
    // for the bounded in-memory displacement transaction.
    const DynamicList<pointField> savedActual(plannedActualPoints_);
    const DynamicList<pointField> savedCollision(acceptedPackingPoints_);
    const DynamicList<boundBox> savedBounds(acceptedPackingBounds_);
    const DynamicList<vector> savedCentres(acceptedPackingCentres_);
    const DynamicList<scalar> savedRadii(acceptedPackingRadii_);

    boolList blockerMask(oldSize, false);
    label directBlockerCount(0);
    forAll(directBlockers, blockerI)
    {
        const label bodyI = directBlockers[blockerI];
        if (bodyI < 0 || bodyI >= oldSize)
        {
            FatalErrorInFunction
                << "Invalid blocker index " << bodyI
                << " for packing plan size " << oldSize << nl
                << exit(FatalError);
        }

        if (!blockerMask[bodyI])
        {
            blockerMask[bodyI] = true;
            directBlockerCount++;
        }
    }

    const label largestRepack = Foam::min(maxLocalRepackBodies_, oldSize);
    if (directBlockerCount > largestRepack)
    {
        return false;
    }

    const label firstRepack = Foam::max
    (
        directBlockerCount,
        Foam::min(minLocalRepackBodies_, largestRepack)
    );

    label repackSize(firstRepack);
    while (repackSize <= largestRepack)
    {
        boolList selected(blockerMask);
        label selectedCount(directBlockerCount);

        // The movable set contains every direct blocker.  Increasing K adds
        // the closest poses, enlarging the local displacement network.
        while (selectedCount < repackSize)
        {
            label nearestBody(-1);
            scalar nearestDistanceSquared(GREAT);

            forAll(savedCentres, bodyI)
            {
                if (selected[bodyI])
                {
                    continue;
                }

                const scalar distanceSquared =
                    magSqr(savedCentres[bodyI] - triggerCentre);
                if (distanceSquared < nearestDistanceSquared)
                {
                    nearestDistanceSquared = distanceSquared;
                    nearestBody = bodyI;
                }
            }

            if (nearestBody < 0)
            {
                break;
            }

            selected[nearestBody] = true;
            selectedCount++;
        }

        vector localMinimum(triggerCentre);
        vector localMaximum(triggerCentre);
        forAll(savedCentres, bodyI)
        {
            if (!selected[bodyI])
            {
                continue;
            }

            for
            (
                direction cmpt = 0;
                cmpt < vector::nComponents;
                ++cmpt
            )
            {
                localMinimum[cmpt] = Foam::min
                (
                    localMinimum[cmpt],
                    savedCentres[bodyI][cmpt]
                );
                localMaximum[cmpt] = Foam::max
                (
                    localMaximum[cmpt],
                    savedCentres[bodyI][cmpt]
                );
            }
        }

        // This is only a centre-proposal window.  Full body support remains
        // constrained by the physical addition boundBox.
        const scalar localPadding = packingRadius_;
        for
        (
            direction cmpt = 0;
            cmpt < vector::nComponents;
            ++cmpt
        )
        {
            localMinimum[cmpt] = Foam::max
            (
                cellZoneBounds_.min()[cmpt],
                localMinimum[cmpt] - localPadding
            );
            localMaximum[cmpt] = Foam::min
            (
                cellZoneBounds_.max()[cmpt],
                localMaximum[cmpt] + localPadding
            );
        }
        const boundBox localCentreBounds(localMinimum, localMaximum);

        // Survivors remain in the ordinary packing arrays and hash.  The local
        // transaction is held separately until all displaced bodies are
        // assigned a conflict-free pose.
        plannedActualPoints_.clear();
        acceptedPackingPoints_.clear();
        acceptedPackingBounds_.clear();
        acceptedPackingCentres_.clear();
        acceptedPackingRadii_.clear();
        packingBins_.clear();

        DynamicList<pointField> localActualPoints;
        DynamicList<pointField> localCollisionPoints;
        DynamicList<boundBox> localBounds;
        DynamicList<vector> localCentres;
        DynamicList<scalar> localRadii;

        for (label bodyI = 0; bodyI < oldSize; ++bodyI)
        {
            if (!selected[bodyI])
            {
                appendPackingPose
                (
                    savedActual[bodyI],
                    savedCollision[bodyI],
                    savedCentres[bodyI],
                    savedRadii[bodyI],
                    savedBounds[bodyI]
                );
            }
            else if (!blockerMask[bodyI])
            {
                localActualPoints.append(savedActual[bodyI]);
                localCollisionPoints.append(savedCollision[bodyI]);
                localBounds.append(savedBounds[bodyI]);
                localCentres.append(savedCentres[bodyI]);
                localRadii.append(savedRadii[bodyI]);
            }
        }

        const label triggerLocalI = localCollisionPoints.size();
        localActualPoints.append(triggerActualPoints);
        localCollisionPoints.append(triggerCollisionPoints);
        localBounds.append(triggerBounds);
        localCentres.append(triggerCentre);
        localRadii.append(triggerRadius);

        bool transactionValid
        (
            packingCandidateBlockers
            (
                triggerCollisionPoints,
                triggerCentre,
                triggerRadius,
                triggerBounds
            ).empty()
        );

        // Removing all direct blockers must leave a conflict-free local set.
        for
        (
            label localI = 0;
            transactionValid && localI < triggerLocalI;
            ++localI
        )
        {
            transactionValid = !packingPosesIntersect
            (
                triggerCollisionPoints,
                triggerCentre,
                triggerRadius,
                triggerBounds,
                localCollisionPoints[localI],
                localCentres[localI],
                localRadii[localI],
                localBounds[localI]
            );
        }

        label unplacedBodies(directBlockerCount);
        label localAttempts(0);

        while
        (
            transactionValid
         && unplacedBodies > 0
         && localAttempts < localRepackAttempts_
        )
        {
            pointField actualPoints;
            pointField collisionPoints;
            vector centre(vector::zero);
            scalar radius(0);
            boundBox bounds;

            const bool candidateFits = createPackingCandidate
            (
                actualPoints,
                collisionPoints,
                centre,
                radius,
                bounds,
                localCentreBounds
            );

            localAttempts++;
            totalPackingAttempts_++;
            bodyAdditionAttemptCounter_++;

            if
            (
                !candidateFits
             || !packingCandidateBlockers
                (
                    collisionPoints,
                    centre,
                    radius,
                    bounds
                ).empty()
            )
            {
                continue;
            }

            label localOverlapCount(0);
            label localOverlapI(-1);
            forAll(localCollisionPoints, localI)
            {
                if
                (
                    packingPosesIntersect
                    (
                        collisionPoints,
                        centre,
                        radius,
                        bounds,
                        localCollisionPoints[localI],
                        localCentres[localI],
                        localRadii[localI],
                        localBounds[localI]
                    )
                )
                {
                    localOverlapCount++;
                    localOverlapI = localI;
                    if (localOverlapCount >= 2)
                    {
                        break;
                    }
                }
            }

            if (localOverlapCount == 0)
            {
                localActualPoints.append(actualPoints);
                localCollisionPoints.append(collisionPoints);
                localBounds.append(bounds);
                localCentres.append(centre);
                localRadii.append(radius);
                unplacedBodies--;
            }
            else if
            (
                localOverlapCount == 1
             && localOverlapI != triggerLocalI
            )
            {
                // Occupy one local slot and displace its previous occupant.
                // The number of unplaced identical bodies is unchanged.
                localActualPoints[localOverlapI] = actualPoints;
                localCollisionPoints[localOverlapI] = collisionPoints;
                localBounds[localOverlapI] = bounds;
                localCentres[localOverlapI] = centre;
                localRadii[localOverlapI] = radius;
            }
        }

        if (transactionValid && unplacedBodies == 0)
        {
            forAll(localCollisionPoints, localI)
            {
                appendPackingPose
                (
                    localActualPoints[localI],
                    localCollisionPoints[localI],
                    localCentres[localI],
                    localRadii[localI],
                    localBounds[localI]
                );
            }

            if (acceptedPackingPoints_.size() == oldSize + 1)
            {
                localRepackCount_++;
                return true;
            }
        }

        // Exact rollback before expanding K or leaving this routine.
        plannedActualPoints_ = savedActual;
        acceptedPackingPoints_ = savedCollision;
        acceptedPackingBounds_ = savedBounds;
        acceptedPackingCentres_ = savedCentres;
        acceptedPackingRadii_ = savedRadii;
        rebuildPackingBins();

        if (repackSize == largestRepack)
        {
            break;
        }
        repackSize = Foam::min(largestRepack, repackSize + 4);
    }

    return false;
}

//---------------------------------------------------------------------------//
void addModelOnceScatter::prepareContinuousPackingPlan()
{
    if (packingPlanReady_)
    {
        return;
    }

    plannedActualPoints_.clear();
    acceptedPackingPoints_.clear();
    acceptedPackingBounds_.clear();
    acceptedPackingCentres_.clear();
    acceptedPackingRadii_.clear();
    packingBins_.clear();

    label failuresSinceProgress(0);
    label bestBlockerCount(targetBodyCount_ + 1);
    DynamicLabelList bestBlockers;
    pointField bestActualPoints;
    pointField bestCollisionPoints;
    vector bestCentre(vector::zero);
    scalar bestRadius(0);
    boundBox bestBounds;

    while (acceptedPackingPoints_.size() < targetBodyCount_)
    {
        pointField actualPoints;
        pointField collisionPoints;
        vector centre(vector::zero);
        scalar radius(0);
        boundBox bounds;

        const bool candidateFits = createPackingCandidate
        (
            actualPoints,
            collisionPoints,
            centre,
            radius,
            bounds,
            cellZoneBounds_
        );

        totalPackingAttempts_++;
        bodyAdditionAttemptCounter_++;

        DynamicLabelList blockers;
        if (candidateFits)
        {
            blockers = packingCandidateBlockers
            (
                collisionPoints,
                centre,
                radius,
                bounds
            );
        }

        if (candidateFits && blockers.empty())
        {
            appendPackingPose
            (
                actualPoints,
                collisionPoints,
                centre,
                radius,
                bounds
            );

            failuresSinceProgress = 0;
            bestBlockerCount = targetBodyCount_ + 1;
            bestBlockers.clear();

            if
            (
                acceptedPackingPoints_.size() % reportEveryAccepted_ == 0
             || acceptedPackingPoints_.size() == targetBodyCount_
            )
            {
                InfoH << addModelSummary_Info
                    << "[onceScatter planner] planned "
                    << acceptedPackingPoints_.size() << '/'
                    << targetBodyCount_ << ", total proposals "
                    << totalPackingAttempts_ << ", local rebuilds "
                    << localRepackCount_ << endl;
            }
            continue;
        }

        failuresSinceProgress++;
        if (candidateFits && blockers.size() < bestBlockerCount)
        {
            bestBlockerCount = blockers.size();
            bestBlockers = blockers;
            bestActualPoints = actualPoints;
            bestCollisionPoints = collisionPoints;
            bestCentre = centre;
            bestRadius = radius;
            bestBounds = bounds;
        }

        if (failuresSinceProgress < stagnationWindow_)
        {
            continue;
        }

        const bool rebuilt =
            !bestBlockers.empty()
         && tryLocalPackingRebuild
            (
                bestActualPoints,
                bestCollisionPoints,
                bestCentre,
                bestRadius,
                bestBounds,
                bestBlockers
            );

        if (rebuilt)
        {
            failuresSinceProgress = 0;
            bestBlockerCount = targetBodyCount_ + 1;
            bestBlockers.clear();

            InfoH << addModelSummary_Info
                << "[onceScatter planner] local rebuild completed, planned "
                << acceptedPackingPoints_.size() << '/'
                << targetBodyCount_ << endl;
            continue;
        }

        plannerRestartCount_++;
        if (plannerRestartCount_ > maxPlannerRestarts_)
        {
            FatalErrorInFunction
                << "continuousRandomBacktracking could not construct all "
                << targetBodyCount_ << " full-size poses after "
                << plannerRestartCount_ << " complete planning attempts. "
                << "No partial packing or structured fallback will be used."
                << nl << "Planned poses before failure: "
                << acceptedPackingPoints_.size()
                << nl << "Total continuous proposals: "
                << totalPackingAttempts_
                << nl << "Completed local rebuilds: "
                << localRepackCount_ << nl
                << exit(FatalError);
        }

        plannedActualPoints_.clear();
        acceptedPackingPoints_.clear();
        acceptedPackingBounds_.clear();
        acceptedPackingCentres_.clear();
        acceptedPackingRadii_.clear();
        packingBins_.clear();
        failuresSinceProgress = 0;
        bestBlockerCount = targetBodyCount_ + 1;
        bestBlockers.clear();

        InfoH << addModelSummary_Info
            << "[onceScatter planner] restarting continuous plan "
            << plannerRestartCount_ << '/' << maxPlannerRestarts_ << endl;
    }

    if
    (
        plannedActualPoints_.size() != targetBodyCount_
     || acceptedPackingPoints_.size() != targetBodyCount_
     || acceptedPackingBounds_.size() != targetBodyCount_
     || acceptedPackingCentres_.size() != targetBodyCount_
     || acceptedPackingRadii_.size() != targetBodyCount_
    )
    {
        FatalErrorInFunction
            << "Internal continuous packing-plan size mismatch." << nl
            << exit(FatalError);
    }

    packingPlanReady_ = true;
    InfoH << addModelSummary_Info
        << "[onceScatter planner] complete: " << targetBodyCount_
        << " continuous full-size poses, " << totalPackingAttempts_
        << " proposals, " << localRepackCount_ << " local rebuilds, "
        << plannerRestartCount_ << " restarts" << endl;
}

//---------------------------------------------------------------------------//
scalar addModelOnceScatter::packingRandom01()
{
    // Every rank constructs this generator from the same dictionary seed.
    // Keeping the exact-packing path branch-identical makes local draws
    // synchronous without an MPI broadcast for every coordinate/quaternion
    // component (which would dominate a high-attempt RSA run).
    return randGen_.scalar01();
}

//---------------------------------------------------------------------------//
tensor addModelOnceScatter::axisAngleRotation
(
    const scalar angle,
    const vector& inputAxis
) const
{
    vector axis(inputAxis);
    const scalar axisMagnitude = mag(axis);
    if (axisMagnitude <= VSMALL)
    {
        return tensor::I;
    }
    axis /= axisMagnitude;

    tensor rotation(Foam::cos(angle)*tensor::I);
    rotation += Foam::sin(angle)*tensor
    (
         0.0,       -axis.z(),  axis.y(),
         axis.z(),   0.0,      -axis.x(),
        -axis.y(),   axis.x(),  0.0
    );
    rotation += (1.0 - Foam::cos(angle))*(axis*axis);
    return rotation;
}

//---------------------------------------------------------------------------//
tensor addModelOnceScatter::randomPackingRotation()
{
    if (!rotateParticles_ || packingIsSphere_)
    {
        return tensor::I;
    }

    if (uniformRandomRotation_)
    {
        // Shoemake's uniform unit-quaternion construction.
        const scalar u1 = geometricVolumeBased_
            ? packingRandom01() : randGen_.scalar01();
        const scalar u2 = geometricVolumeBased_
            ? packingRandom01() : randGen_.scalar01();
        const scalar u3 = geometricVolumeBased_
            ? packingRandom01() : randGen_.scalar01();
        const scalar twoPi = 2.0*Foam::constant::mathematical::pi;

        const scalar qx = Foam::sqrt(1.0 - u1)*Foam::sin(twoPi*u2);
        const scalar qy = Foam::sqrt(1.0 - u1)*Foam::cos(twoPi*u2);
        const scalar qz = Foam::sqrt(u1)*Foam::sin(twoPi*u3);
        const scalar qw = Foam::sqrt(u1)*Foam::cos(twoPi*u3);

        return tensor
        (
            1.0 - 2.0*(qy*qy + qz*qz),
            2.0*(qx*qy - qz*qw),
            2.0*(qx*qz + qy*qw),
            2.0*(qx*qy + qz*qw),
            1.0 - 2.0*(qx*qx + qz*qz),
            2.0*(qy*qz - qx*qw),
            2.0*(qx*qz - qy*qw),
            2.0*(qy*qz + qx*qw),
            1.0 - 2.0*(qx*qx + qy*qy)
        );
    }

    vector axis(axisOfRot_);
    if (randomAxis_)
    {
        axis = vector
        (
            2.0*packingRandom01() - 1.0,
            2.0*packingRandom01() - 1.0,
            2.0*packingRandom01() - 1.0
        );
        if (mag(axis) <= VSMALL)
        {
            axis = vector(1, 0, 0);
        }
    }

    const scalar angle =
        (2.0*packingRandom01() - 1.0)
       *Foam::constant::mathematical::pi;
    return axisAngleRotation(angle, axis);
}

//---------------------------------------------------------------------------//
bool addModelOnceScatter::createPackingCandidate
(
    pointField& actualPoints,
    pointField& collisionPoints,
    vector& centre,
    scalar& radius,
    boundBox& collisionBounds,
    const boundBox& centreSamplingBounds
)
{
    const tensor rotation(randomPackingRotation());

    if (packingIsSphere_)
    {
        radius = sourceSphereRadius_*fixedScale_*clearanceScale_;
        const vector radiusVector(radius, radius, radius);
        vector minimumCentre(cellZoneBounds_.min() + radiusVector);
        vector maximumCentre(cellZoneBounds_.max() - radiusVector);

        for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
        {
            minimumCentre[cmpt] = Foam::max
            (
                minimumCentre[cmpt],
                centreSamplingBounds.min()[cmpt]
            );
            maximumCentre[cmpt] = Foam::min
            (
                maximumCentre[cmpt],
                centreSamplingBounds.max()[cmpt]
            );

            if (maximumCentre[cmpt] < minimumCentre[cmpt])
            {
                return false;
            }
            centre[cmpt] = minimumCentre[cmpt]
                + packingRandom01()
                 *(maximumCentre[cmpt] - minimumCentre[cmpt]);
        }

        actualPoints.setSize(1);
        actualPoints[0] = centre;
        collisionPoints = actualPoints;
        collisionBounds = boundBox
        (
            centre - radiusVector,
            centre + radiusVector
        );
        return true;
    }

    actualPoints = sourcePoints_;
    actualPoints -= sourceCentre_;
    actualPoints *= fixedScale_;
    actualPoints = rotation & actualPoints;
    actualPoints += sourceCentre_;
    centre = sourceCentre_;

    collisionPoints = actualPoints;
    collisionPoints -= centre;
    collisionPoints *= clearanceScale_;
    collisionPoints += centre;

    collisionBounds = boundBox(collisionPoints, false);
    const vector minimumTranslation
    (
        cellZoneBounds_.min() - collisionBounds.min()
    );
    const vector maximumTranslation
    (
        cellZoneBounds_.max() - collisionBounds.max()
    );

    vector minimumCentre(sourceCentre_ + minimumTranslation);
    vector maximumCentre(sourceCentre_ + maximumTranslation);

    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        minimumCentre[cmpt] = Foam::max
        (
            minimumCentre[cmpt],
            centreSamplingBounds.min()[cmpt]
        );
        maximumCentre[cmpt] = Foam::min
        (
            maximumCentre[cmpt],
            centreSamplingBounds.max()[cmpt]
        );

        if (maximumCentre[cmpt] < minimumCentre[cmpt])
        {
            return false;
        }
        centre[cmpt] = minimumCentre[cmpt]
            + packingRandom01()
             *(maximumCentre[cmpt] - minimumCentre[cmpt]);
    }

    const vector translation(centre - sourceCentre_);
    actualPoints += translation;
    collisionPoints += translation;
    collisionBounds = boundBox(collisionPoints, false);

    scalar maximumRadiusSquared = 0;
    forAll(collisionPoints, pointI)
    {
        maximumRadiusSquared = Foam::max
        (
            maximumRadiusSquared,
            magSqr(collisionPoints[pointI] - centre)
        );
    }
    radius = Foam::sqrt(maximumRadiusSquared);
    return true;
}

//---------------------------------------------------------------------------//
label addModelOnceScatter::packingBinKey(const vector& centre) const
{
    const vector relative(centre - cellZoneBounds_.min());
    label ix = label(std::floor(relative.x()/packingBinSize_));
    label iy = label(std::floor(relative.y()/packingBinSize_));
    label iz = label(std::floor(relative.z()/packingBinSize_));

    ix = Foam::max(label(0), Foam::min(ix, nPackingBinsX_ - 1));
    iy = Foam::max(label(0), Foam::min(iy, nPackingBinsY_ - 1));
    iz = Foam::max(label(0), Foam::min(iz, nPackingBinsZ_ - 1));

    return ix + nPackingBinsX_*(iy + nPackingBinsY_*iz);
}

//---------------------------------------------------------------------------//
bool addModelOnceScatter::packingPosesIntersect
(
    const pointField& firstPoints,
    const vector& firstCentre,
    const scalar firstRadius,
    const boundBox& firstBounds,
    const pointField& secondPoints,
    const vector& secondCentre,
    const scalar secondRadius,
    const boundBox& secondBounds
) const
{
    const scalar combinedRadius = firstRadius + secondRadius;
    if (magSqr(firstCentre - secondCentre) > sqr(combinedRadius))
    {
        return false;
    }

    if (!firstBounds.overlaps(secondBounds))
    {
        return false;
    }

    if (packingIsSphere_)
    {
        return true;
    }

    return convexBodiesIntersect
    (
        firstPoints,
        secondPoints,
        firstCentre - secondCentre,
        combinedRadius
    );
}

//---------------------------------------------------------------------------//
addModelOnceScatter::DynamicLabelList
addModelOnceScatter::packingCandidateBlockers
(
    const pointField& points,
    const vector& centre,
    const scalar radius,
    const boundBox& bounds
) const
{
    DynamicLabelList blockers;
    const label centreKey = packingBinKey(centre);
    const label centreX = centreKey % nPackingBinsX_;
    const label centreYZ = centreKey/nPackingBinsX_;
    const label centreY = centreYZ % nPackingBinsY_;
    const label centreZ = centreYZ/nPackingBinsY_;

    for (label dz = -1; dz <= 1; ++dz)
    {
        const label iz = centreZ + dz;
        if (iz < 0 || iz >= nPackingBinsZ_) continue;

        for (label dy = -1; dy <= 1; ++dy)
        {
            const label iy = centreY + dy;
            if (iy < 0 || iy >= nPackingBinsY_) continue;

            for (label dx = -1; dx <= 1; ++dx)
            {
                const label ix = centreX + dx;
                if (ix < 0 || ix >= nPackingBinsX_) continue;

                const label key =
                    ix + nPackingBinsX_*(iy + nPackingBinsY_*iz);
                if (!packingBins_.found(key)) continue;

                const DynamicLabelList& neighbours = packingBins_[key];
                forAll(neighbours, neighbourI)
                {
                    const label bodyI = neighbours[neighbourI];
                    if
                    (
                        packingPosesIntersect
                        (
                            points,
                            centre,
                            radius,
                            bounds,
                            acceptedPackingPoints_[bodyI],
                            acceptedPackingCentres_[bodyI],
                            acceptedPackingRadii_[bodyI],
                            acceptedPackingBounds_[bodyI]
                        )
                    )
                    {
                        blockers.append(bodyI);
                    }
                }
            }
        }
    }
    return blockers;
}

//---------------------------------------------------------------------------//
vector addModelOnceScatter::minkowskiSupport
(
    const pointField& first,
    const pointField& second,
    const vector& direction
)
{
    label firstI = 0;
    label secondI = 0;
    scalar firstProjection = first[0] & direction;
    scalar secondProjection = second[0] & (-direction);

    for (label pointI = 1; pointI < first.size(); ++pointI)
    {
        const scalar projection = first[pointI] & direction;
        if (projection > firstProjection)
        {
            firstProjection = projection;
            firstI = pointI;
        }
    }
    for (label pointI = 1; pointI < second.size(); ++pointI)
    {
        const scalar projection = second[pointI] & (-direction);
        if (projection > secondProjection)
        {
            secondProjection = projection;
            secondI = pointI;
        }
    }
    return first[firstI] - second[secondI];
}

//---------------------------------------------------------------------------//
bool addModelOnceScatter::sameDirection
(
    const vector& first,
    const vector& second
)
{
    return (first & second) > 0;
}

//---------------------------------------------------------------------------//
bool addModelOnceScatter::updateSimplex
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
bool addModelOnceScatter::convexBodiesIntersect
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

    // Degenerate/non-converged cases are treated conservatively as contact.
    return true;
}

//---------------------------------------------------------------------------//
label addModelOnceScatter::maxConsecutiveFailures() const
{
    // Geometric packing is completed atomically before registration.  The
    // driver therefore receives only prevalidated poses, never RSA failures.
    return geometricVolumeBased_ ? -1 : 1000;
}

//---------------------------------------------------------------------------//
label addModelOnceScatter::maxTotalAttempts() const
{
    return geometricVolumeBased_ ? targetBodyCount_ : -1;
}

//---------------------------------------------------------------------------//
bool addModelOnceScatter::strictCompletion() const
{
    return geometricVolumeBased_;
}

//---------------------------------------------------------------------------//
bool addModelOnceScatter::targetReached() const
{
    return !geometricVolumeBased_
        || finishedAddition_
        || acceptedBodyCount_ >= targetBodyCount_;
}

//---------------------------------------------------------------------------//
void addModelOnceScatter::bodyRegistered
(
    const word& bodyName,
    const volScalarField&
)
{
    if (!geometricVolumeBased_)
    {
        return;
    }
    if (!packingPlanReady_ || acceptedBodyCount_ >= targetBodyCount_)
    {
        FatalErrorInFunction
            << "The driver registered unexpected body " << bodyName
            << " at plan index " << acceptedBodyCount_ << '.' << nl
            << exit(FatalError);
    }

    acceptedBodyCount_++;

    if
    (
        acceptedBodyCount_ % reportEveryAccepted_ == 0
     || acceptedBodyCount_ == targetBodyCount_
    )
    {
        InfoH << addModelSummary_Info
            << "[onceScatter " << bodyName << "] accepted "
            << acceptedBodyCount_ << '/' << targetBodyCount_
            << " from prevalidated continuous plan"
            << ", geometric fraction "
            << acceptedBodyCount_*particleVolume_/regionVolume_ << endl;
    }
}

//---------------------------------------------------------------------------//
void addModelOnceScatter::reportAdditionSummary
(
    const word& bodyName,
    const volScalarField& body,
    const word& stopReason
)
{
    if (!geometricVolumeBased_)
    {
        return;
    }
    if (finishedAddition_ && acceptedBodyCount_ == 0)
    {
        // onceScatter is intentionally inactive when starting from a saved
        // non-zero time; the bodies are restored by the normal restart path.
        return;
    }

    const scalar lambdaFraction = checkLambdaFraction(body);
    InfoH << addModelSummary_Info
        << "[onceScatter " << bodyName << "] stopped: " << stopReason
        << ", strategy " << packingStrategy_
        << ", accepted " << acceptedBodyCount_ << '/' << targetBodyCount_
        << ", total pose attempts " << totalPackingAttempts_
        << ", local rebuilds " << localRepackCount_
        << ", planner restarts " << plannerRestartCount_
        << ", geometric fraction "
        << acceptedBodyCount_*particleVolume_/regionVolume_
        << ", lambda fraction " << lambdaFraction << endl;
}
// MODEL SPECIFIC FUNCTIONS==================================================//
//---------------------------------------------------------------------------//
void addModelOnceScatter::initializeCellZone()
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
void addModelOnceScatter::updateCellZoneBoundBox()
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
void addModelOnceScatter::initializeBoundBox()
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
void addModelOnceScatter::initializePolygonPrism()
{
    octreeField_ *= 0;
    cellsInBoundBox_[Pstream::myProcNo()].clear();

    const pointField& cellCenters = mesh_.C();

    forAll(cellCenters, cellI)
    {
        if (polygonPrismGeometry_->contains(cellCenters[cellI]))
        {
            cellsInBoundBox_[Pstream::myProcNo()].append(cellI);
        }
    }

    cellZoneBounds_ = polygonPrismGeometry_->bounds();

    label globalCellCount =
        cellsInBoundBox_[Pstream::myProcNo()].size();
    reduce(globalCellCount, sumOp<label>());

    InfoH << addModel_Info << "-- addModelMessage-- "
        << "polygonPrism selected cells: " << globalCellCount << endl;
}
//---------------------------------------------------------------------------//
labelList addModelOnceScatter::getBBoxCellsByOctTree
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
scalar addModelOnceScatter::checkLambdaFraction(const volScalarField& body)
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
scalar addModelOnceScatter::returnRandomAngle()
{
    const scalar unitRandom = polygonPrismActive_
      ? returnSynchronizedRandom01()
      : randGen_.scalar01();

    scalar ranNum = 2.0*unitRandom - 1.0;
    scalar angle  = ranNum*Foam::constant::mathematical::pi;
	return angle;
}
//---------------------------------------------------------------------------//
scalar addModelOnceScatter::returnRandomScale()
{
	scalar ranNum = polygonPrismActive_
      ? returnSynchronizedRandom01()
      : randGen_.scalar01();
	scalar scaleDiff    = maxScale_ - minScale_;
    scalar scaleFactor  = minScale_ + ranNum*scaleDiff;
	InfoH << addModel_Info << "-- addModelMessage-- "
        <<"random scaleFactor " << scaleFactor <<endl;
	return scaleFactor;
}
//---------------------------------------------------------------------------//
vector addModelOnceScatter::returnRandomRotationAxis()
{
	vector  axisOfRotation(vector::zero);
	scalar ranNum = 0;

	for (int i=0;i<3;i++)
	{
		ranNum = polygonPrismActive_
          ? returnSynchronizedRandom01()
          : randGen_.scalar01();
		axisOfRotation[i] = ranNum;
	}

	axisOfRotation /=mag(axisOfRotation);
	return axisOfRotation;
}
//---------------------------------------------------------------------------//
scalar addModelOnceScatter::returnSynchronizedRandom01()
{
    return randGen_.globalScalar01();
}
//---------------------------------------------------------------------------//
pointField addModelOnceScatter::relativeBodySupportPoints
(
    scalar& isotropicRadius
)
{
    isotropicRadius = 0;
    const point bodyCenter = geomModel_->getCoM();

    if (geomModel_->getcType() == sphere)
    {
        isotropicRadius = 0.5*geomModel_->getDC();
        return pointField(1, vector::zero);
    }

    pointField bodyPoints = geomModel_->getBodyPoints();

    if (bodyPoints.size() == 0)
    {
        // Conservative fallback for geometry types that do not expose their
        // surface points (for example a compound body implementation).
        const boundBox bodyBounds = geomModel_->getBounds();
        const point& minimum = bodyBounds.min();
        const point& maximum = bodyBounds.max();

        bodyPoints.setSize(8);
        bodyPoints[0] = point(minimum.x(), minimum.y(), minimum.z());
        bodyPoints[1] = point(maximum.x(), minimum.y(), minimum.z());
        bodyPoints[2] = point(minimum.x(), maximum.y(), minimum.z());
        bodyPoints[3] = point(maximum.x(), maximum.y(), minimum.z());
        bodyPoints[4] = point(minimum.x(), minimum.y(), maximum.z());
        bodyPoints[5] = point(maximum.x(), minimum.y(), maximum.z());
        bodyPoints[6] = point(minimum.x(), maximum.y(), maximum.z());
        bodyPoints[7] = point(maximum.x(), maximum.y(), maximum.z());
    }

    forAll(bodyPoints, pointI)
    {
        bodyPoints[pointI] -= bodyCenter;
    }

    return bodyPoints;
}
//---------------------------------------------------------------------------//
