// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file PSpraySolutionFanInjection.H
 * @brief 在平面扇形注入上设置入口温度、LiBr浓度与一致的液体物性。
 * XLiBr为kg/kg，inletTemperature为K，surfaceTension为N/m。
 * 纯水使用XLiBr=0；p.Y()的单组分占位不代表把LiBr作为可挥发物质。
 */
#ifndef PSpraySolutionFanInjection_H
#define PSpraySolutionFanInjection_H
#include "PSprayFlatFanInjection.H"
#include "PSpraySolutionThermo.H"

namespace Foam
{
template <class C> class PSpraySolutionFanInjection : public PSprayFlatFanInjection<C>
{
    scalar X_, T_, sigma_;

public:
    TypeName("pSpraySolutionFanInjection");

    PSpraySolutionFanInjection(const dictionary& dict, C& owner, const word& name)
        : PSprayFlatFanInjection<C>(dict, owner, name)
        , X_(dict.lookup<scalar>("XLiBr"))
        , T_(dict.lookup<scalar>("inletTemperature"))
        , sigma_(dict.lookup<scalar>("surfaceTension"))
    {
        if (X_ < 0 || X_ > 0.65 || T_ < 278.15 || T_ > 353.15 || sigma_ <= 0)
            FatalIOErrorInFunction(dict) << "Invalid solution inlet state" << exit(FatalIOError);
    }

    PSpraySolutionFanInjection(const PSpraySolutionFanInjection& p)
        : PSprayFlatFanInjection<C>(p)
        , X_(p.X_)
        , T_(p.T_)
        , sigma_(p.sigma_)
    {
    }

    virtual autoPtr<InjectionModel<C>> clone() const
    {
        return autoPtr<InjectionModel<C>>(new PSpraySolutionFanInjection(*this));
    }

    /// 先由父模型确定速度/粒径，再根据入口T、X设置rho、Cp、mu与表面张力。
    virtual void setProperties(const label i,
                               const label n,
                               const scalar t,
                               typename C::parcelType::trackingData& td,
                               typename C::parcelType& p)
    {
        PSprayFlatFanInjection<C>::setProperties(i, n, t, td, p);
        const auto a = pSprayTable(this->owner().constProps().dict()).evaluate(T_, X_);
        p.XLiBr() = X_;
        p.T() = T_;
        p.rho() = a.rho;
        p.Cp() = a.cp;
        p.mu() = a.mu;
        p.sigma() = sigma_;
        p.Y() = scalarField(1, 1);
    }

    /// 所需状态已由注入器完整设置，避免云的默认纯水属性覆盖这些值。
    virtual bool fullyDescribed() const
    {
        return true;
    }
};
}
#endif
