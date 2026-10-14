// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file librSolutionThermo.H
 * @brief 注册浓度相关液体热物性，并为每个单元/边界提供真实组分。
 * XLiBr为必读、自动写出的无量纲质量分数字段；尾迹标记允许从重启场恢复。
 */
#ifndef librSolutionThermo_H
#define librSolutionThermo_H
#include "rhoFluidThermo.H"
#include "LiBrSolutionMixture.H"

namespace Foam
{
// Retain the fluid-thermo interface, but supply actual solution composition
// to the thermodynamic mixture instead of the pureThermo nil composition.
class LiBrSolutionComposite : public rhoFluidThermo::composite
{
    // 物性与pSprayFilm盐输运访问同一个注册字段，避免重复的组分状态。
    volScalarField XLiBr_;
    volScalarField conductivityTail_;

public:
    template <class MixtureType>
    LiBrSolutionComposite(const dictionary& dict,
                          const MixtureType& mixture,
                          const fvMesh& mesh,
                          const word& phaseName)
        : rhoFluidThermo::composite(dict, mixture, mesh, phaseName)
        , XLiBr_(IOobject(dict.lookupOrDefault<word>("concentrationField",
                                                     IOobject::groupName("XLiBr", phaseName)),
                          mesh.time().name(),
                          mesh,
                          IOobject::MUST_READ,
                          IOobject::AUTO_WRITE),
                 mesh)
        , conductivityTail_(IOobject("pSprayConductivityTail",
                                     mesh.time().name(),
                                     mesh,
                                     IOobject::READ_IF_PRESENT,
                                     IOobject::AUTO_WRITE),
                            mesh,
                            dimensionedScalar(dimless, 0),
                            "zeroGradient")
    {
    }

    /// 热物性框架的nil占位输入由该单元的XLiBr与尾迹状态替换。
    LiBrSolutionComposition cellComposition(const nil, const label celli) const
    {
        return {XLiBr_[celli], conductivityTail_[celli] > .5};
    }

    LiBrSolutionComposition
    patchFaceComposition(const nil, const label patchi, const label facei) const
    {
        return {XLiBr_.boundaryField()[patchi][facei],
                conductivityTail_.boundaryField()[patchi][facei] > .5};
    }

    /// 未提供独立溶液组分源的体积注入不能使用此接口，避免默认按纯水注入。
    LiBrSolutionComposition sourceCellComposition(const nil, const label) const
    {
        FatalErrorInFunction
            << "Volume injection requires an explicit solution-composition source; use configured inlet patches"
            << exit(FatalError);
        return {0, false};
    }
};

typedef RhoFluidThermo<BasicThermo<LiBrSolutionMixture, LiBrSolutionComposite>>
    LiBrSolutionThermoBase;

class librSolutionThermo : public LiBrSolutionThermoBase
{
public:
    TypeName("librSolutionThermo");

    librSolutionThermo(const fvMesh& mesh, const word& phaseName)
        : LiBrSolutionThermoBase(mesh, phaseName)
    {
    }
};
}
#endif
