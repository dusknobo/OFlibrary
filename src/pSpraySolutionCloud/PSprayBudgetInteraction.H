// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file PSprayBudgetInteraction.H
 * @brief 在原生LocalInteraction上记录离开计算域的液滴累计收支。
 * 记录总液体质量[kg]、非挥发盐[kg]、液体内能[J]和动能[J]。
 * 写出时将MPI全局累计值存入模型属性，重启后接续该账本。
 */
#ifndef PSprayBudgetInteraction_H
#define PSprayBudgetInteraction_H
#include "LocalInteraction.H"
#include "PSpraySolutionThermo.H"

namespace Foam
{
template <class CloudType> class PSprayBudgetInteraction : public LocalInteraction<CloudType>
{
    scalar escapedMass_ = 0, escapedSalt_ = 0, escapedEnergy_ = 0, escapedKinetic_ = 0;

public:
    TypeName("pSprayBudgetInteraction");

    PSprayBudgetInteraction(const dictionary& d, CloudType& c)
        : LocalInteraction<CloudType>(d, c)
    {
    }

    PSprayBudgetInteraction(const PSprayBudgetInteraction& a)
        : LocalInteraction<CloudType>(a)
        , escapedMass_(a.escapedMass_)
        , escapedSalt_(a.escapedSalt_)
        , escapedEnergy_(a.escapedEnergy_)
        , escapedKinetic_(a.escapedKinetic_)
    {
    }

    virtual autoPtr<PatchInteractionModel<CloudType>> clone() const
    {
        return autoPtr<PatchInteractionModel<CloudType>>(new PSprayBudgetInteraction(*this));
    }

    /// 在原生逃逸逻辑重置速度前保存液滴状态，仅对实际删除的逃逸粒子记账。
    virtual bool correct(typename CloudType::parcelType& p, const polyPatch& pp, bool& keep)
    {
        // Film capture is handled first by native MomentumParcel::hitPatch.
        // Snapshot before native escape resets the parcel velocity.
        // mass()是一滴质量；乘nParticle得到本parcel携带的广延质量。
        const scalar m = p.mass() * p.nParticle();
        const scalar x = p.XLiBr(), t = p.T(), kinetic = .5 * m * magSqr(p.U());
        const bool result = LocalInteraction<CloudType>::correct(p, pp, keep);
        if (result && !keep)
        {
            const auto st = pSprayTable(this->owner().constProps().dict()).evaluate(t, x);
            escapedMass_ += m;
            escapedSalt_ += m * x;
            escapedEnergy_ += m * (st.h - st.pEq / st.rho);
            escapedKinetic_ += kinetic;
        }
        return result;
    }

    /// 模型持久化账本加当前分区增量；写出后清零增量，避免重复累计。
    virtual void info(Ostream& os)
    {
        LocalInteraction<CloudType>::info(os);
        const scalar m = this->template getModelProperty<scalar>("escapedLiquidMass") +
                         returnReduce(escapedMass_, sumOp<scalar>());
        const scalar s = this->template getModelProperty<scalar>("escapedSaltMass") +
                         returnReduce(escapedSalt_, sumOp<scalar>());
        const scalar e = this->template getModelProperty<scalar>("escapedLiquidInternalEnergy") +
                         returnReduce(escapedEnergy_, sumOp<scalar>());
        const scalar k = this->template getModelProperty<scalar>("escapedKineticEnergy") +
                         returnReduce(escapedKinetic_, sumOp<scalar>());
        os << "    PSPRAY_ESCAPE cloud=" << this->owner().name() << " mass=" << m << " salt=" << s
           << " E=" << e << " K=" << k << nl;
        if (this->writeTime())
        {
            // Same global ledger on all ranks, including across restarts.
            this->setModelProperty("escapedLiquidMass", m);
            this->setModelProperty("escapedSaltMass", s);
            this->setModelProperty("escapedLiquidInternalEnergy", e);
            this->setModelProperty("escapedKineticEnergy", k);
            escapedMass_ = escapedSalt_ = escapedEnergy_ = escapedKinetic_ = 0;
        }
    }
};
}
#endif
