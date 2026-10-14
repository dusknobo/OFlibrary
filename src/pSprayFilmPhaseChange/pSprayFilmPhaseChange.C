// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file pSprayFilmPhaseChange.C
 * @brief 将局部双温度反应预测量映射为成对的气相/膜有限体积源项。
 * update在同一timeIndex只计算一次，因此存在与空间输运的时间分裂误差。
 * 与Implicit/Smooth版本每次外迭代重新计算的本构源项实现应明确区分。
 */
#include "pSprayFilmPhaseChange.H"
#include "mappedFvPatchBaseBase.H"
#include "fvModels.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace fv
{
defineTypeNameAndDebug(pSprayFilmPhaseChange, 0);
addToRunTimeSelectionTable(fvModel, pSprayFilmPhaseChange, dictionary);
defineTypeNameAndDebug(pSprayVapourFilmSource, 0);
addToRunTimeSelectionTable(fvModel, pSprayVapourFilmSource, dictionary);

static std::string tablePath(const dictionary& d)
{
    fileName p(d.lookup<fileName>("tableFile"));
    p.expand();
    return std::string(p.c_str());
}

pSprayFilmPhaseChange::pSprayFilmPhaseChange(const word& name,
                                             const word& type,
                                             const fvMesh& mesh,
                                             const dictionary& d)
    : fvModel(name, type, mesh, d)
    , film_(mesh.lookupObject<solvers::isothermalFilm>(solver::typeName))
    , table_(tablePath(d))
    , accommodation_(d.lookup<scalar>("accommodation"))
    , wetThickness_(d.lookup<scalar>("wetThickness"))
    , maxFraction_(d.lookupOrDefault<scalar>("maxFractionChange", .1))
    , cpv_(d.lookupOrDefault<scalar>("vapourCp", 1860))
    , href_(d.lookupOrDefault<scalar>("vapourHRef", 2547276.823116537))
    , Tref_(d.lookupOrDefault<scalar>("vapourTRef", 298.15))
    , phaseEnabled_(d.lookupOrDefault<Switch>("phaseEnabled", true))
    , phaseStartTime_(d.lookupOrDefault<scalar>("phaseStartTime", 0))
    , phaseRampDuration_(d.lookupOrDefault<scalar>("phaseRampDuration", 0))
    , massRate_(
          IOobject(
              "phaseMassRate", mesh.time().name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
          mesh,
          dimensionedScalar(dimDensity / dimTime, 0))
    , energyRate_(
          IOobject(
              "phaseEnergyRate", mesh.time().name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
          mesh,
          dimensionedScalar(dimPower / dimVolume, 0))
    , momentumRate_(IOobject("phaseMomentumRate",
                             mesh.time().name(),
                             mesh,
                             IOobject::NO_READ,
                             IOobject::AUTO_WRITE),
                    mesh,
                    dimensionedVector(dimDensity * dimVelocity / dimTime, Zero))
{
    if (!(accommodation_ > 0 && accommodation_ <= 1 && wetThickness_ > 0 && maxFraction_ > 0 &&
          maxFraction_ <= .5 && phaseRampDuration_ >= 0))
        FatalIOErrorInFunction(d) << "Invalid kinetic coefficient, wet scale or mass-fraction limit"
                                  << exit(FatalIOError);
}

wordList pSprayFilmPhaseChange::addSupFields() const
{
    return {film_.alpha.name(), film_.thermo.he().name(), film_.U.name()};
}

void pSprayFilmPhaseChange::requireUpdated() const
{
    if (lastTime_ != mesh().time().timeIndex())
        FatalErrorInFunction << "Paired gas pSprayVapourFilmSource must be corrected first"
                             << exit(FatalError);
}

void pSprayFilmPhaseChange::update(const scalarField& allocatedGasVolume)
{
    // 冻结整步源项，避免同一时间步多次消耗局部可用液体/蒸气库存。
    if (lastTime_ == mesh().time().timeIndex())
        return;
    lastTime_ = mesh().time().timeIndex();
    massRate_ = dimensionedScalar(massRate_.dimensions(), 0);
    energyRate_ = dimensionedScalar(energyRate_.dimensions(), 0);
    momentumRate_ = dimensionedVector(momentumRate_.dimensions(), Zero);
    const scalar age = mesh().time().value() - phaseStartTime_;
    const scalar r = !phaseEnabled_ || age < 0 ? scalar(0)
                     : phaseRampDuration_ > 0  ? min(age / phaseRampDuration_, scalar(1))
                                               : scalar(1);
    // smoothstep只安排启动渐启；物理适应系数由accommodation单独给出。
    const scalar phaseScale = r * r * (3 - 2 * r);
    if (phaseScale <= 0)
    {
        Info << "Film phase " << mesh().name()
             << " evaporation kg/s=0 absorption kg/s=0 active=0 limited=0"
             << " limitedTransfer kg=0 energyToGasW=0 predictorResidual kg=0"
             << " startupScale=0" << nl;
        return;
    }
    const auto& map = film_.surfacePatchMap();
    const auto& gm = map.nbrMesh();
    const label gp = map.nbrFvPatch().index();
    const auto& fp = film_.surfacePatch();
    const scalarField tg(map.fromNeighbour(
        gm.lookupObject<volScalarField>("T").boundaryField()[gp].patchInternalField()()));
    const scalarField pg(map.fromNeighbour(
        gm.lookupObject<volScalarField>("p").boundaryField()[gp].patchInternalField()()));
    const scalarField hg(map.fromNeighbour(
        gm.lookupObject<volScalarField>("h").boundaryField()[gp].patchInternalField()()));
    const vectorField ug(map.fromNeighbour(
        gm.lookupObject<volVectorField>("U").boundaryField()[gp].patchInternalField()()));
    const scalarField ag(map.fromNeighbour(gm.boundary()[gp].magSf()));
    const scalarField vg(map.fromNeighbour(allocatedGasVolume));
    const auto& X = mesh().lookupObject<volScalarField>("XLiBr");
    const auto& tl = film_.thermo.T();
    const scalar dt = mesh().time().deltaTValue();
    constexpr scalar Rv = 8.31446261815324 / .018015268;
    const scalar tailThickness =
        film_.thermo.properties().lookupOrDefault<scalar>("numericalTailThickness", 1e-8);
    scalar evap = 0, absorb = 0, energyToGas = 0, limitedMass = 0, maxResidual = 0;
    label limited = 0, active = 0;
    forAll (fp, i)
    {
        const label c = fp.faceCells()[i];
        const scalar ml = film_.alpha[c] * film_.rho[c] * mesh().V()[c];
        const scalar area =
            phaseScale * ag[i] * min(max(film_.delta[c] / wetThickness_, scalar(0)), scalar(1));
        if (ml <= 1e-20 || area <= 0 || film_.delta[c] < tailThickness)
            continue;
        if (mag(hg[i] - (href_ + cpv_ * (tg[i] - Tref_))) > 1e-5)
            FatalErrorInFunction << "Gas must use the declared COMMON reference hConst enthalpy; h="
                                 << hg[i] << ", T=" << tg[i] << exit(FatalError);
        // 纯水蒸气按理想气体求质量；vg是按接触面积分配的气相体积。
        const scalar mg = pg[i] / (Rv * tg[i]) * vg[i];
        try
        {
            librSpray::FilmReactionInput a{ml,
                                           mg,
                                           vg[i],
                                           area,
                                           tl[c],
                                           tg[i],
                                           X[c],
                                           dt,
                                           accommodation_,
                                           maxFraction_,
                                           cpv_,
                                           href_,
                                           Tref_};
            // 预测器返回的是dm[kg]，先除dt再除膜体积，得到体积质量源。
            const auto z = librSpray::filmReaction(table_, a);
            const scalar rate = z.dm / dt;
            massRate_[c] = rate / mesh().V()[c];
            energyRate_[c] = massRate_[c] * z.hv;
            energyToGas += rate * z.hv;
            // Momentum carried by transferred water; the film solver projects out
            // normal velocity and applies the corresponding wall reaction.
            const vector donorU = rate >= 0 ? film_.U[c] : ug[i];
            momentumRate_[c] = massRate_[c] * donorU;
            evap += max(rate, scalar(0));
            absorb += max(-rate, scalar(0));
            if (z.limited)
            {
                ++limited;
                limitedMass += mag(z.dm);
            }
            maxResidual = max(maxResidual, mag(z.residual));
            ++active;
        }
        catch (const std::exception& e)
        {
            FatalErrorInFunction << "Film reaction cell=" << c << ", T=" << tl[c] << ", X=" << X[c]
                                 << ", ml=" << ml << ", mg=" << mg << ", area=" << area << ": "
                                 << e.what() << exit(FatalError);
        }
    }
    Info << "Film phase " << mesh().name()
         << " evaporation kg/s=" << returnReduce(evap, sumOp<scalar>())
         << " absorption kg/s=" << returnReduce(absorb, sumOp<scalar>())
         << " active=" << returnReduce(active, sumOp<label>())
         << " limited=" << returnReduce(limited, sumOp<label>())
         << " limitedTransfer kg=" << returnReduce(limitedMass, sumOp<scalar>())
         << " energyToGasW=" << returnReduce(energyToGas, sumOp<scalar>())
         << " predictorResidual kg=" << returnReduce(maxResidual, maxOp<scalar>())
         << " startupScale=" << phaseScale << nl;
}

void pSprayFilmPhaseChange::addSup(const volScalarField&,
                                   const volScalarField& alpha,
                                   fvMatrix<scalar>& eqn) const
{
    requireUpdated();
    if (&alpha != &film_.alpha)
        NotImplemented;
    // mDot>0代表蒸发，因此膜连续性中的源为负。
    eqn -= massRate_();
}

void pSprayFilmPhaseChange::addSup(const volScalarField&,
                                   const volScalarField&,
                                   const volScalarField& he,
                                   fvMatrix<scalar>& eqn) const
{
    requireUpdated();
    if (&he != &film_.thermo.he())
        NotImplemented;
    eqn -= energyRate_();
}

void pSprayFilmPhaseChange::addSup(const volScalarField&,
                                   const volScalarField&,
                                   const volVectorField&,
                                   fvMatrix<vector>& eqn) const
{
    requireUpdated();
    eqn -= momentumRate_();
}

pSprayVapourFilmSource::pSprayVapourFilmSource(const word& name,
                                               const word& type,
                                               const fvMesh& mesh,
                                               const dictionary& d)
    : fvModel(name, type, mesh, d)
    , filmRegions_(d.lookup<wordList>("filmRegions"))
    , filmModel_(d.lookupOrDefault<word>("filmModel", "waterPhaseChange"))
    , massRate_(
          IOobject(
              "filmMassRate", mesh.time().name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
          mesh,
          dimensionedScalar(dimDensity / dimTime, 0))
    , energyRate_(
          IOobject(
              "filmEnergyRate", mesh.time().name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
          mesh,
          dimensionedScalar(dimPower / dimVolume, 0))
    , momentumRate_(IOobject("filmMomentumRate",
                             mesh.time().name(),
                             mesh,
                             IOobject::NO_READ,
                             IOobject::AUTO_WRITE),
                    mesh,
                    dimensionedVector(dimDensity * dimVelocity / dimTime, Zero))
{
}

wordList pSprayVapourFilmSource::addSupFields() const
{
    return {"rho", "h", "U"};
}

void pSprayVapourFilmSource::correct()
{
    if (lastTime_ == mesh().time().timeIndex())
        return;
    lastTime_ = mesh().time().timeIndex();
    massRate_ = dimensionedScalar(massRate_.dimensions(), 0);
    energyRate_ = dimensionedScalar(energyRate_.dimensions(), 0);
    momentumRate_ = dimensionedVector(momentumRate_.dimensions(), Zero);
    DynamicList<pSprayFilmPhaseChange*> models;
    scalarField sumArea(mesh().nCells(), 0);
    forAll (filmRegions_, ri)
    {
        const auto& fm = mesh().time().lookupObject<fvMesh>(filmRegions_[ri]);
        auto& fvms = fvModels::New(fm);
        pSprayFilmPhaseChange* found = nullptr;
        forAll (fvms, j)
            if (fvms[j].name() == filmModel_)
                found = &refCast<pSprayFilmPhaseChange>(fvms[j]);
        if (!found)
            FatalErrorInFunction << "Missing " << filmModel_ << " in " << fm.name()
                                 << exit(FatalError);
        models.append(found);
        const auto& patch = found->film().surfacePatchMap().nbrFvPatch();
        forAll (patch, i)
            sumArea[patch.faceCells()[i]] += patch.magSf()[i];
    }
    scalar totalFilmM = 0, totalFilmE = 0, absoluteFilmM = 0, absoluteFilmE = 0;
    forAll (models, j)
    {
        auto& model = *models[j];
        const auto& film = model.film();
        const auto& map = film.surfacePatchMap();
        const auto& gp = map.nbrFvPatch();
        const auto& fc = film.surfacePatch().faceCells();
        scalarField volume(gp.size(), 0);
        // 同一气相单元可能接触多个膜面；按面积分配其蒸气体积，避免重复库存。
        forAll (gp, i)
        {
            const label c = gp.faceCells()[i];
            volume[i] = mesh().V()[c] * gp.magSf()[i] / sumArea[c];
        }
        model.update(volume);
        // Send extensive per-face transfers. Normalise by the mapped GAS face
        // area, not the extruded film wall area, to preserve the paired source.
        const scalarField ag(map.fromNeighbour(gp.magSf()));
        scalarField m(fc.size()), e(fc.size());
        vectorField u(fc.size());
        forAll (fc, i)
        {
            const label c = fc[i];
            const scalar v = model.mesh().V()[c];
            m[i] = model.massRate()[c] * v / ag[i];
            e[i] = model.energyRate()[c] * v / ag[i];
            u[i] = model.momentumRate()[c] * v / ag[i];
            totalFilmM += model.massRate()[c] * v;
            totalFilmE += model.energyRate()[c] * v;
            absoluteFilmM += mag(model.massRate()[c] * v);
            absoluteFilmE += mag(model.energyRate()[c] * v);
        }
        const scalarField gm(map.toNeighbour(m)), ge(map.toNeighbour(e));
        const vectorField gu(map.toNeighbour(u));
        forAll (gp, i)
        {
            const label c = gp.faceCells()[i];
            const scalar av = gp.magSf()[i] / mesh().V()[c];
            massRate_[c] += gm[i] * av;
            energyRate_[c] += ge[i] * av;
            momentumRate_[c] += gu[i] * av;
        }
    }
    totalFilmM = returnReduce(totalFilmM, sumOp<scalar>());
    totalFilmE = returnReduce(totalFilmE, sumOp<scalar>());
    absoluteFilmM = returnReduce(absoluteFilmM, sumOp<scalar>());
    absoluteFilmE = returnReduce(absoluteFilmE, sumOp<scalar>());
    // 检查映射前后广延质量率[kg/s]和焓功率[W]，而非只比较面密度。
    const scalar errM = gSum(massRate_.primitiveField() * mesh().V()) - totalFilmM;
    const scalar errE = gSum(energyRate_.primitiveField() * mesh().V()) - totalFilmE;
    if (mag(errM) > 1e-10 * max(absoluteFilmM, scalar(1e-15)) ||
        mag(errE) > 1e-10 * max(absoluteFilmE, scalar(1e-10)))
        FatalErrorInFunction << "Nonconservative film mapping: dm/dt=" << errM << " dE/dt=" << errE
                             << exit(FatalError);
    Info << "Paired film-vapour source mapping residual kg/s=" << errM << " W=" << errE << nl;
}

void pSprayVapourFilmSource::addSup(const volScalarField&, fvMatrix<scalar>& eqn) const
{
    eqn += massRate_();
}

void pSprayVapourFilmSource::addSup(const volScalarField&,
                                    const volScalarField& he,
                                    fvMatrix<scalar>& eqn) const
{
    if (he.name() != "h")
        NotImplemented;
    eqn += energyRate_();
}

void pSprayVapourFilmSource::addSup(const volScalarField&,
                                    const volVectorField&,
                                    fvMatrix<vector>& eqn) const
{
    eqn += momentumRate_();
}
}
}
