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
#include "prtSubContactInfo.H"
#include "contactModelInfo.H"

using namespace Foam;

//---------------------------------------------------------------------------//
prtSubContactInfo::prtSubContactInfo
(
    const Tuple2<label,label>& contactPair,
    const physicalProperties& physicalProperties
)
:
contactPair_(contactPair),
physicalProperties_(physicalProperties)
{}

prtSubContactInfo::~prtSubContactInfo()
{}
//---------------------------------------------------------------------------//
vector prtSubContactInfo::getVeli(ibContactVars& cVars, vector& lVec)
{
    return (-((lVec-cVars.Axis_*((lVec) & cVars.Axis_))
        ^ cVars.Axis_)*cVars.omega_+ cVars.Vel_);
}
//---------------------------------------------------------------------------//
void prtSubContactInfo::evalVariables(
    ibContactClass& cIb,
    ibContactClass& tIb,
    ibContactVars& cVars,
    ibContactVars& tVars
)
{
    cLVec_ = cIb.getGeomModel().getLVec(prtCntVars_.contactCenter_);
    // cLVec_ = prtCntVars_.contactCenter_ - cIb.getGeomModel().getCoM();
    tLVec_ = tIb.getGeomModel().getLVec(prtCntVars_.contactCenter_);
    // tLVec_ = prtCntVars_.contactCenter_ - tIb.getGeomModel().getCoM();

    cVeli_ = getVeli(cVars, cLVec_);
    tVeli_ = getVeli(tVars, tLVec_);

    Vn_ = -(cVeli_ - tVeli_) & prtCntVars_.contactNormal_;
    Lc_ = (contactModelInfo::getLcCoeff())*mag(cLVec_)*mag(tLVec_)/(mag(cLVec_) + mag(tLVec_));

    physicalProperties_.curAdhN_ = min
    (
        physicalProperties_.maxAdhN_,
        max(physicalProperties_.curAdhN_, physicalProperties_.aY_*prtCntVars_.contactVolume_/(sqr(Lc_)*8
            *Foam::constant::mathematical::pi))
    );
}
//---------------------------------------------------------------------------//
vector prtSubContactInfo::getFNe()
{
    return (physicalProperties_.aY_*prtCntVars_.contactVolume_/(Lc_+SMALL))
        *prtCntVars_.contactNormal_;
}
//---------------------------------------------------------------------------//
vector prtSubContactInfo::getFA()
{
    return ((sqrt(8*Foam::constant::mathematical::pi*physicalProperties_.aY_
        *physicalProperties_.curAdhN_*prtCntVars_.contactVolume_))
        *prtCntVars_.contactNormal_);
}
//---------------------------------------------------------------------------//
vector prtSubContactInfo::getFNd()
{
    return (physicalProperties_.reduceBeta_*sqrt(physicalProperties_.aY_
        *physicalProperties_.reduceM_*prtCntVars_.contactArea_/(Lc_+SMALL))*
        Vn_)*prtCntVars_.contactNormal_;

}
//---------------------------------------------------------------------------//
vector prtSubContactInfo::getFt(scalar deltaT, scalar maxFt)
{
    const scalar normalMag(mag(prtCntVars_.contactNormal_));

    if (normalMag <= SMALL || maxFt <= SMALL || deltaT <= 0)
    {
        FtPrev_ = vector::zero;
        return vector::zero;
    }

    const vector normal(prtCntVars_.contactNormal_/normalMag);

    // Rotate the elastic history into the current contact plane while
    // preserving its magnitude.  Reset it if the projection is singular.
    const scalar FtPrevMag(mag(FtPrev_));
    const vector FtLastP
    (
        FtPrev_
      - (FtPrev_ & normal)*normal
    );
    const scalar FtLastPMag(mag(FtLastP));

    vector FtElastic(vector::zero);
    if
    (
        FtPrevMag > SMALL
     && FtLastPMag > sqrt(SMALL)*FtPrevMag
    )
    {
        FtElastic = FtLastP*(FtPrevMag/FtLastPMag);
    }

    // compute relative tangential velocity
    const vector relVeli(cVeli_ - tVeli_);
    const vector veliNorm
    (
        normal*(relVeli & normal)
    );
    const vector Vt(relVeli - veliNorm);

    scalar kT(0);
    scalar dampingT(0);

    if(contactModelInfo::getUseMindlinRotationalModel())
    {
        if
        (
            physicalProperties_.aG_ > SMALL
         && prtCntVars_.contactArea_ > SMALL
         && Lc_ > SMALL
         && physicalProperties_.reduceM_ > SMALL
        )
        {
            // OpenHFDIB-DEM contact model, Appendix A, Eq. (A.11).
            kT = 8*physicalProperties_.aG_
                *(prtCntVars_.contactArea_/Lc_);

            dampingT = 2*physicalProperties_.reduceBeta_
                *sqrt(kT*physicalProperties_.reduceM_);
        }
    }
    else if(contactModelInfo::getUseChenRotationalModel())
    {
        if
        (
            physicalProperties_.aG_ > SMALL
         && Lc_ > SMALL
         && physicalProperties_.reduceM_ > SMALL
        )
        {
            // Retain the Chen stiffness and damping magnitudes while using
            // the same elastic-history convention as the Mindlin path.
            kT = physicalProperties_.aG_*Lc_;
            dampingT = physicalProperties_.reduceBeta_
                *sqrt(kT*physicalProperties_.reduceM_);
        }
    }

    if (kT <= SMALL)
    {
        FtPrev_ = vector::zero;
        return vector::zero;
    }

    // Only the elastic increment is accumulated.  Damping is instantaneous.
    vector FtElasticTrial(FtElastic - kT*Vt*deltaT);
    const vector FtDamping(-dampingT*Vt);
    vector Ft(FtElasticTrial + FtDamping);

    // Apply the Coulomb slider to the total trial force and back-correct the
    // elastic state so no force above the limit remains hidden in history.
    const scalar frictionLimit(maxFt);
    const scalar FtMag(mag(Ft));
    if (FtMag > frictionLimit && FtMag > SMALL)
    {
        Ft *= frictionLimit/FtMag;
        FtElasticTrial = Ft - FtDamping;
    }

    FtPrev_ = FtElasticTrial;
    return Ft;
}
//---------------------------------------------------------------------------//
void prtSubContactInfo::setVMInfo(boundBox& bBox, scalar subVolumeV)
{
    if (!vmInfo_)
    {
        vmInfo_ = std::make_shared<virtualMeshInfo>(bBox, subVolumeV);
        return;
    }

    vmInfo_->sV = subVolume(bBox);
    vmInfo_->subVolumeV = subVolumeV;
}
//---------------------------------------------------------------------------//
void prtSubContactInfo::setVMInfo(const virtualMeshInfo& vmInfo)
{
    if (!vmInfo_)
    {
        vmInfo_ = std::make_shared<virtualMeshInfo>(vmInfo);
        return;
    }

    vmInfo_->sV = vmInfo.sV;
    vmInfo_->subVolumeV = vmInfo.subVolumeV;
    vmInfo_->startingPoint = vmInfo.startingPoint;
}
//---------------------------------------------------------------------------//
std::shared_ptr<virtualMeshInfo>& prtSubContactInfo::getVMInfo()
{
    return vmInfo_;
}
//---------------------------------------------------------------------------//
void prtSubContactInfo::syncData()
{
    reduce(outForce_.first().F, sumOp<vector>());
    reduce(outForce_.first().T, sumOp<vector>());
    reduce(outForce_.second().F, sumOp<vector>());
    reduce(outForce_.second().T, sumOp<vector>());


    if (vmInfo_)
    {
        point reducePoint = vector::zero;

        if (contactResolved_)
        {
            reducePoint = vmInfo_->getStartingPoint();
        }

        reduce(reducePoint, sumOp<vector>());
        vmInfo_->startingPoint.reset(new point(reducePoint));
    }
}
//---------------------------------------------------------------------------//


// ************************************************************************* //
