// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file pSprayFilmPhaseChange.H
 * @brief 基线气膜相变：局部反应预测后，在一个时间步内冻结配对源项。
 * massRate正为蒸发[kg/(m³·s)]；energyRate为向气相转移的焓功率[W/m³]。
 * 液膜使用负源，气相使用正源；LiBr不挥发，因此不向盐方程加入相变源。
 * pSprayVapourFilmSource先更新，再由膜侧addSup读取同一组转移量。
 */
#ifndef pSprayFilmPhaseChange_H
#define pSprayFilmPhaseChange_H
#include "fvModel.H"
#include "isothermalFilm.H"
#include "PSprayFilmReaction.H"

namespace Foam
{
namespace fv
{
class pSprayFilmPhaseChange : public fvModel
{
    const solvers::isothermalFilm& film_;
    librSpray::Table table_;
    scalar accommodation_, wetThickness_, maxFraction_, cpv_, href_, Tref_;
    bool phaseEnabled_;
    scalar phaseStartTime_, phaseRampDuration_;
    label lastTime_ = -1;
    volScalarField massRate_, energyRate_;
    volVectorField momentumRate_;
    void requireUpdated() const;

public:
    TypeName("pSprayFilmPhaseChange");
    pSprayFilmPhaseChange(const word&, const word&, const fvMesh&, const dictionary&);

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

class pSprayVapourFilmSource : public fvModel
{
    wordList filmRegions_;
    word filmModel_;
    label lastTime_ = -1;
    volScalarField massRate_, energyRate_;
    volVectorField momentumRate_;

public:
    TypeName("pSprayVapourFilmSource");
    pSprayVapourFilmSource(const word&, const word&, const fvMesh&, const dictionary&);
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
