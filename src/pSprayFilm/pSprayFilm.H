// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file pSprayFilm.H
 * @brief 在Foundation 14原生film上增加非挥发LiBr输运与守恒诊断。
 * 主变量为alpha、U、比内能e和XLiBr；膜厚delta=alpha*VbyA。
 * T和rho由e、X反求。该模型使用厚度平均状态，不求深液池自由液面。
 */
#ifndef pSprayFilm_H
#define pSprayFilm_H
#include "film.H"
#include "LiBrThermoTable.H"

namespace Foam
{
namespace solvers
{
class pSprayFilm : public film
{
    // X引用热物性已注册的浓度场；互扩散系数单位m²/s。
    volScalarField& X_;
    dimensionedScalar diffusivity_;
    librSpray::Table table_;
    // 尾迹约束和连续性匹配项的人工贡献单独输出，不从总误差中扣除。
    scalar tailEnergyPower_ = 0;
    scalar tailLiquidMass_ = 0;
    scalar tailSaltRate_ = 0;
    scalar tailSaltMass_ = 0;
    scalar continuitySaltRate_ = 0;
    scalar continuitySaltAbsoluteRate_ = 0;
    /// 配合最新质量通量求盐方程，再构造同一守恒能量方程的初始迭代。
    void correctComposition();
    /// 湿区连通面上的D*alpha*rho，采用有界调和插值。
    tmp<surfaceScalarField> saltDiffusionCoefficient() const;

public:
    TypeName("pSprayFilm");
    pSprayFilm(fvMesh& mesh);
    /// 原生质量预测后同步组分，供随后的热物性查询使用。
    virtual void prePredictor();
    /// 质量修正后同步盐，并重解相同时间步的能量方程。
    virtual void pressureCorrector();
    /// 包含盐扩散携带焓；求解后按保存质量重算alpha和delta。
    virtual void thermophysicalPredictor();
    /// 输出质量、盐、能量和壁面热流，供累计收支核验。
    virtual void postSolve();
};
}
}
#endif
