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

Description
    class for addition and removal of immersedBodies

SourceFiles
    addModelOnceScatter.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef addModelOnceScatter_H
#define addModelOnceScatter_H

#include "addModel.H"
#include "polygonPrismGeometry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class addModelOnceScatter Declaration
\*---------------------------------------------------------------------------*/
class addModelOnceScatter
:
    public addModel
{

    private:
    typedef DynamicList<label>  DynamicLabelList;
    const dictionary&           addModelDict_;
    word                        addMode_;
    bool                        bodyAdded_;
    bool                        finishedAddition_;

    const dictionary&   coeffsDict_;

	word			    addDomain_;
    word			    scalingMode_;
    word			    rotationMode_;
	word			    addModeI_;

	const dictionary& addDomainCoeffs_;
	const dictionary& scalingModeCoeffs_;
	const dictionary& rotationModeCoeffs_;
	const dictionary& addModeICoeffs_;

    label		partPerAdd_;
	scalar		fieldValue_;
    label		addedOnTimeLevel_;
    const scalar startTimeValue_;

	word		zoneName_;
    vector		minBound_;
    vector		maxBound_;

	bool 		scaleParticles_;
	scalar 		minScale_;
	scalar		maxScale_;
	scalar 		minScaleFit_;
	scalar		scaleStep_;
	scalar		nTriesBeforeScaling_;

	bool 		rotateParticles_;
	bool		randomAxis_;
    vector 		axisOfRot_;

	int 		bodyAdditionAttemptCounter_;
	int 		scaleCorrectionCounter_;

	bool 		scaleApplication_;
	bool 		scaleRandomApplication_;
	bool 		rescaleRequirement_;
	bool 		succesfulladition_;
	scalar		scalingFactor_;
	bool		restartPartCountTemp_;
	bool		reapeatedAddition_;
	bool		firstTimeRunning_;
	bool 		cellZoneActive_;
	bool		boundBoxActive_;
	bool		polygonPrismActive_;
	Field<label>octreeField_;
	bool		multiBody_;
    bool		fieldBased_;
	bool		geometricVolumeBased_;
	scalar		fieldCurrentValue_;

    // Exact-volume packing state.  The complete continuous-pose plan is
    // constructed before the first immersedBody is registered, so initial
    // backtracking never mutates solver/contact state.
    scalar      targetGeometricFraction_;
    scalar      regionVolume_;
    scalar      particleVolume_;
    scalar      fixedScale_;
    scalar      clearanceScale_;
    label       targetBodyCount_;
    label       acceptedBodyCount_;
    label       totalPackingAttempts_;
    label       reportEveryAccepted_;

    bool        packingPlanReady_;
    label       stagnationWindow_;
    label       localRepackAttempts_;
    label       minLocalRepackBodies_;
    label       maxLocalRepackBodies_;
    label       maxPlannerRestarts_;
    label       localRepackCount_;
    label       plannerRestartCount_;

    word        packingStrategy_;

    bool        uniformRandomRotation_;
    bool        packingIsSphere_;

    pointField  sourcePoints_;
    vector      sourceCentre_;
    scalar      sourceSphereRadius_;
    scalar      packingRadius_;
    scalar      packingBinSize_;
    label       nPackingBinsX_;
    label       nPackingBinsY_;
    label       nPackingBinsZ_;

    DynamicList<pointField> acceptedPackingPoints_;
    DynamicList<pointField> plannedActualPoints_;
    DynamicList<boundBox>   acceptedPackingBounds_;
    DynamicList<vector>     acceptedPackingCentres_;
    DynamicList<scalar>     acceptedPackingRadii_;
    HashTable<DynamicLabelList, label, Hash<label>> packingBins_;

    List<DynamicLabelList>  cellsInBoundBox_;
    List<pointField>        cellZonePoints_;
    boundBox                cellZoneBounds_;

    std::unique_ptr<polygonPrismGeometry> polygonPrismGeometry_;

    bool        allActiveCellsInMesh_;

    Random          randGen_;

    scalar returnSynchronizedRandom01();

    pointField relativeBodySupportPoints
    (
        scalar& isotropicRadius
    );

    // Exact-volume packing helpers --------------------------------
    void initializeGeometricVolumePacking();

    void prepareContinuousPackingPlan();

    bool tryLocalPackingRebuild
    (
        const pointField& triggerActualPoints,
        const pointField& triggerCollisionPoints,
        const vector& triggerCentre,
        const scalar triggerRadius,
        const boundBox& triggerBounds,
        const DynamicLabelList& directBlockers
    );

    void appendPackingPose
    (
        const pointField& actualPoints,
        const pointField& collisionPoints,
        const vector& centre,
        const scalar radius,
        const boundBox& bounds
    );

    void rebuildPackingBins();

    std::shared_ptr<geomModel> addGeometricVolumeBody
    (
        PtrList<immersedBody>& immersedBodies
    );

    scalar packingRandom01();

    tensor randomPackingRotation();

    tensor axisAngleRotation(const scalar angle, const vector& axis) const;

    bool createPackingCandidate
    (
        pointField& actualPoints,
        pointField& collisionPoints,
        vector& centre,
        scalar& radius,
        boundBox& collisionBounds,
        const boundBox& centreSamplingBounds
    );

    bool packingPosesIntersect
    (
        const pointField& firstPoints,
        const vector& firstCentre,
        const scalar firstRadius,
        const boundBox& firstBounds,
        const pointField& secondPoints,
        const vector& secondCentre,
        const scalar secondRadius,
        const boundBox& secondBounds
    ) const;

    DynamicLabelList packingCandidateBlockers
    (
        const pointField& points,
        const vector& centre,
        const scalar radius,
        const boundBox& bounds
    ) const;

    label packingBinKey(const vector& centre) const;

    static vector minkowskiSupport
    (
        const pointField& first,
        const pointField& second,
        const vector& direction
    );

    static bool sameDirection(const vector& first, const vector& second);

    static bool updateSimplex(DynamicList<vector>& simplex, vector& direction);

    static bool convexBodiesIntersect
    (
        const pointField& first,
        const pointField& second,
        const vector& initialDirection,
        const scalar lengthScale
    );

    public:

    // Constructors ------------------------------------------------
    addModelOnceScatter
    (
        const dictionary& addModelDict,
        const Foam::fvMesh& mesh,
        const bool startTime0,
        std::unique_ptr<geomModel> bodyGeomModel,
        List<labelList>& cellPoints
    );

    // Destructors -------------------------------------------------
    ~addModelOnceScatter();

    // Public subclasses -------------------------------------------

    // Public variables --------------------------------------------

    // Public methods ----------------------------------------------
    void init();

    bool shouldAddBody(const volScalarField& body);

    std::shared_ptr<geomModel> addBody(
        const volScalarField& body,
        PtrList<immersedBody>& immersedBodies
    );

	scalar returnRandomAngle();

	scalar returnRandomScale();

    vector returnRandomRotationAxis();

	void initializeCellZone();

	void updateCellZoneBoundBox();

    void initializeBoundBox();

    void initializePolygonPrism();

	labelList getBBoxCellsByOctTree
    (
        label cellToCheck,
        bool& insideBB,
        vector& bBoxMin,
        vector& bBoxMax,
        List<DynamicLabelList>& bBoxCells
    );

	scalar checkLambdaFraction(const volScalarField& body);

    bool checkIfAdd()
    {
        return true;
    }

    word getAddMode()
    {
        return addMode_;
    }

    bool getBodyAdded()
    {
        return bodyAdded_;
    }

    label maxConsecutiveFailures() const override;

    label maxTotalAttempts() const override;

    bool strictCompletion() const override;

    bool targetReached() const override;

    void bodyRegistered
    (
        const word& bodyName,
        const volScalarField& body
    ) override;

    void reportAdditionSummary
    (
        const word& bodyName,
        const volScalarField& body,
        const word& stopReason
    ) override;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
