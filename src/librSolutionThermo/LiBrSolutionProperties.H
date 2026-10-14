// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file LiBrSolutionProperties.H
 * @brief OpenFOAM液体热物性适配器：组分相关、不可压缩，以比内能e为能量变量。
 *
 * 热物性表与蒸气能量采用共同参考。压力参数用绝对Pa，温度用K。
 * 浓度范围0≤X≤0.65；温度求逆范围278.15—353.15 K。
 * 实际浓溶液导热还要求T≥293.15 K，数值尾迹的导热参考单独标记。
 */
#ifndef LiBrSolutionProperties_H
#define LiBrSolutionProperties_H
#include "dictionary.H"
#include "fileName.H"
#include "LiBrThermoTable.H"
#include "PSprayConductivity.H"
#include <memory>
#include <algorithm>

namespace Foam
{
class LiBrSolutionProperties
{
    // 所有组分查询共享只读物性表；X_与尾迹状态由mixture查询逐点设置。
    std::shared_ptr<const librSpray::Table> table_;
    scalar X_;
    bool numericalTail_ = false;
    scalar kappaVerification_;
    bool conductivityCorrelation_;

    /// 展开OpenFOAM路径变量；运行算例应引用包内的物性表。
    static std::string tablePath(const dictionary& dict)
    {
        fileName path(dict.lookup<fileName>("tableFile"));
        path.expand();
        return std::string(path.c_str());
    }

    /// 查表越界转换为OpenFOAM FatalError；末尾return仅满足编译器返回路径。
    librSpray::State state(const scalar T) const
    {
        try
        {
            return table_->evaluate(T, X_);
        }
        catch (const std::exception& error)
        {
            FatalErrorInFunction << "T=" << T << ", XLiBr=" << X_ << ": " << error.what()
                                 << exit(FatalError);
        }
        return table_->evaluate(313.15, 0.64);
    }

    /// 固定p、X下对h或e求温度导数[J/(kg·K)]，在表端点缩短差分区间。
    scalar energyDerivative(const scalar p, const scalar T, const bool internal) const
    {
        const scalar lo = std::max(scalar(278.15), T - 1e-4);
        const scalar hi = std::min(scalar(353.15), T + 1e-4);
        if (!(hi > lo))
            FatalErrorInFunction << "Invalid temperature derivative interval" << exit(FatalError);
        return internal ? (he(p, hi) - he(p, lo)) / (hi - lo) : (ha(p, hi) - ha(p, lo)) / (hi - lo);
    }

public:
    static constexpr bool incompressible = true;
    static constexpr bool isochoric = false;

    static word typeName()
    {
        return "librSolutionProperties";
    }

    static word heName()
    {
        return "e";
    }

    static bool enthalpy()
    {
        return false;
    }

    explicit LiBrSolutionProperties(const dictionary& dict)
        : table_(std::make_shared<const librSpray::Table>(tablePath(dict)))
        , X_(0)
        , kappaVerification_(dict.lookupOrDefault<scalar>("verificationThermalConductivity", 0))
        , conductivityCorrelation_(dict.lookup<word>("thermalConductivityMode") == "ASHRAE526RP")
    {
        if (!conductivityCorrelation_ &&
            (dict.lookup<word>("thermalConductivityMode") != "constantForVerification" ||
             !(kappaVerification_ > 0)))
            FatalIOErrorInFunction(dict)
                << "This prototype requires a positive explicitly labelled verification thermal conductivity"
                << exit(FatalIOError);
    }

    /// 设置当前质量分数；越界直接报错，不以截断浓度掩盖盐质量超限。
    void composition(const scalar X, const bool numericalTail = false)
    {
        if (!(X >= 0 && X <= 0.65))
            FatalErrorInFunction << "XLiBr outside [0,0.65]: " << X << exit(FatalError);
        X_ = X;
        numericalTail_ = numericalTail;
    }

    /// 混合物摩尔质量[kg/kmol]，由水和LiBr的质量分数换算。
    scalar W() const
    {
        return 1 / ((1 - X_) / 18.015268 + X_ / 86.85);
    }

    scalar rho(const scalar, const scalar T) const
    {
        return state(T).rho;
    }

    /// 不可压缩模型的∂rho/∂p为零；rho仍随T和X变化。
    scalar psi(const scalar, const scalar) const
    {
        return 0;
    }

    scalar mu(const scalar, const scalar T) const
    {
        return state(T).mu;
    }

    /**
     * @brief 导热系数[W/(m·K)]。
     * constantForVerification仅用于显式标记的组件核验。
     * 真实浓溶液低于关联式有效温度时报错；数值尾迹仅使用有界导热参考，
     * 不修改该位置真实的T、e、rho或X。
     */
    scalar kappa(const scalar, const scalar T) const
    {
        if (!conductivityCorrelation_)
            return kappaVerification_;
        if (X_ >= .3 && T < 293.15 && !numericalTail_)
            FatalErrorInFunction
                << "ASHRAE conductivity for resolved X>=0.3 liquid requires T>=293.15 K; T=" << T
                << ", X=" << X_ << exit(FatalError);
        // There is no continuum liquid layer below the native film's minimum-thickness numerical
        // resolution. Use a bounded transport reference there, never change its
        // thermodynamic T/e/rho or claim a low-T concentrated-liquid fit.
        return librSpray::conductivity(numericalTail_ && X_ >= .3 ? max(T, scalar(293.15)) : T, X_);
    }

    // Common HEOS-water reference. There is no separate formation-enthalpy shift.
    /// 表内焓位于平衡压力pEq，按h(p)=hEq+(p-pEq)/rho修正至调用压力。
    scalar ha(const scalar p, const scalar T) const
    {
        const auto s = state(T);
        return s.h + (p - s.pEq) / s.rho;
    }

    scalar hs(const scalar p, const scalar T) const
    {
        return ha(p, T);
    }

    scalar hf() const
    {
        return 0;
    }

    /// 比内能[J/kg]：e=h(p)-p/rho=hEq-pEq/rho，显式压力项相消。
    scalar he(const scalar, const scalar T) const
    {
        const auto s = state(T);
        // The table enthalpy is evaluated on the solution saturation curve.
        // For an incompressible liquid, h(p)=h_eq+(p-p_eq)/rho;
        // subtracting p/rho therefore leaves u=h_eq-p_eq/rho.
        // Omitting the enthalpy pressure correction creates a spurious
        // negative du/dp and destabilises pressure/temperature iterations.
        return s.h - s.pEq / s.rho;
    }

    /// Cp与Cv分别从本适配器实际采用的h(p,T,X)、e(T,X)求导。
    scalar Cp(const scalar p, const scalar T) const
    {
        return energyDerivative(p, T, false);
    }

    scalar Cv(const scalar p, const scalar T) const
    {
        return energyDerivative(p, T, true);
    }

    scalar Cpv(const scalar p, const scalar T) const
    {
        return Cv(p, T);
    }

    /**
     * @brief 从比内能反求温度[K]。
     * 采用有界Newton迭代；试探点落在括区外时退回二分。
     * 超出能量可表示范围或60次内未收敛时报错，保证热物性定义一致。
     */
    scalar The(const scalar energy, const scalar p, const scalar guess) const
    {
        scalar lo = 278.15, hi = 353.15;
        const scalar lowEnergy = he(p, lo), highEnergy = he(p, hi);
        if (energy < lowEnergy - 1e-6 || energy > highEnergy + 1e-6)
            FatalErrorInFunction << "Energy outside tabulated temperature range: e=" << energy
                                 << ", X=" << X_ << exit(FatalError);
        scalar T = std::max(lo, std::min(hi, guess));
        for (label iter = 0; iter < 60; ++iter)
        {
            const scalar residual = he(p, T) - energy;
            if (mag(residual) < 1e-6)
                return T;
            if (residual > 0)
                hi = T;
            else
                lo = T;
            const scalar trial = T - residual / Cv(p, T);
            T = trial > lo && trial < hi ? trial : (lo + hi) / 2;
        }
        FatalErrorInFunction << "Energy-temperature inversion failed" << exit(FatalError);
        return T;
    }
};
}
#endif
