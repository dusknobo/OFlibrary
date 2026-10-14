// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file pSpraySmoothFilmPhaseChange.H
 * @brief 每次外迭代更新的实验性气膜相变源项。
 * 在Implicit本构源项上对极薄膜截止区增加smoothstep激活；过渡宽度单位m。
 * massRate>0为蒸发；膜侧减、气侧加同一质量/动量/焓源，盐不转移。
 * pressureDerivative保存-∂massRate/∂p，供气相压力/密度方程隐式修正。
 * alphaSink、energySink分别用于液膜质量和能量方程的稳定隐式修正。
 * 本类不调用共享filmReaction预测器；须核验外迭代和时间步精度。
 */
#ifndef pSpraySmoothFilmPhaseChange_H
#define pSpraySmoothFilmPhaseChange_H
#include "fvModel.H"
#include "isothermalFilm.H"
#include "PSprayFilmReaction.H"

namespace Foam
{
namespace fv
{
class pSpraySmoothFilmPhaseChange : public fvModel
{
    const solvers::isothermalFilm& film_;
    librSpray::Table table_;
    scalar accommodation_, wetThickness_, tailTransitionWidth_, maxFraction_, cpv_, href_, Tref_;
    bool phaseEnabled_;
    scalar phaseStartTime_, phaseRampDuration_;
    label lastTime_ = -1;
    volScalarField massRate_, energyRate_, pressureDerivative_, alphaSink_, energySink_;
    volVectorField momentumRate_;
    void requireUpdated() const;

public:
    TypeName("pSpraySmoothFilmPhaseChange");
    pSpraySmoothFilmPhaseChange(const word&, const word&, const fvMesh&, const dictionary&);

    const solvers::isothermalFilm& film() const
    {
        return film_;
    }

    const volScalarField& massRate() const
    {
        return massRate_;
    }

    const volScalarField& energyRate() const
    {
        return energyRate_;
    }

    const volScalarField& pressureDerivative() const
    {
        return pressureDerivative_;
    }

    const volVectorField& momentumRate() const
    {
        return momentumRate_;
    }

    void update(const scalarField& allocatedGasVolume);
    virtual wordList addSupFields() const;
    virtual void addSup(const volScalarField&, const volScalarField&, fvMatrix<scalar>&) const;
    virtual void addSup(const volScalarField&,
                        const volScalarField&,
                        const volScalarField&,
                        fvMatrix<scalar>&) const;
    virtual void addSup(const volScalarField&,
                        const volScalarField&,
                        const volVectorField&,
                        fvMatrix<vector>&) const;

    virtual bool movePoints()
    {
        return false;
    }

    virtual void topoChange(const polyTopoChangeMap&)
    {
        NotImplemented;
    }

    virtual void mapMesh(const polyMeshMap&)
    {
        NotImplemented;
    }

    virtual void distribute(const polyDistributionMap&)
    {
        NotImplemented;
    }
};

class pSpraySmoothVapourFilmSource : public fvModel
{
    wordList filmRegions_;
    word filmModel_;
    label lastTime_ = -1;
    volScalarField massRate_, energyRate_, pressureDerivative_;
    volVectorField momentumRate_;

public:
    TypeName("pSpraySmoothVapourFilmSource");
    pSpraySmoothVapourFilmSource(const word&, const word&, const fvMesh&, const dictionary&);
    virtual wordList addSupFields() const;
    virtual void correct();
    virtual void addSup(const volScalarField&, fvMatrix<scalar>&) const;
    virtual void addSup(const volScalarField&, const volScalarField&, fvMatrix<scalar>&) const;
    virtual void addSup(const volScalarField&, const volVectorField&, fvMatrix<vector>&) const;

    virtual bool movePoints()
    {
        return false;
    }

    virtual void topoChange(const polyTopoChangeMap&)
    {
        NotImplemented;
    }

    virtual void mapMesh(const polyMeshMap&)
    {
        NotImplemented;
    }

    virtual void distribute(const polyDistributionMap&)
    {
        NotImplemented;
    }
};
}
}
#endif
