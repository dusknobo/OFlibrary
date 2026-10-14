// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file pSpraySmoothFilmPhaseChange.C
 * @brief 实验性相变本构通量、Jacobian近似和守恒气膜源项映射。
 * 每次外迭代按当前场更新源项；系数在该次矩阵求解中固定。
 * correction(Sp)围绕当前迭代线性化，在收敛点消失，不额外指定相变量。
 */
#include "pSpraySmoothFilmPhaseChange.H"
#include "mappedFvPatchBaseBase.H"
#include "fvModels.H"
#include "fvMatrix.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace fv
{
defineTypeNameAndDebug(pSpraySmoothFilmPhaseChange, 0);
addToRunTimeSelectionTable(fvModel, pSpraySmoothFilmPhaseChange, dictionary);
defineTypeNameAndDebug(pSpraySmoothVapourFilmSource, 0);
addToRunTimeSelectionTable(fvModel, pSpraySmoothVapourFilmSource, dictionary);

static std::string tablePath(const dictionary& d)
{
    fileName p(d.lookup<fileName>("tableFile"));
    p.expand();
    return std::string(p.c_str());
}

pSpraySmoothFilmPhaseChange::pSpraySmoothFilmPhaseChange(const word& name,
                                                         const word& type,
                                                         const fvMesh& mesh,
                                                         const dictionary& d)
    : fvModel(name, type, mesh, d)
    , film_(mesh.lookupObject<solvers::isothermalFilm>(solver::typeName))
    , table_(tablePath(d))
    , accommodation_(d.lookup<scalar>("accommodation"))
    , wetThickness_(d.lookup<scalar>("wetThickness"))
    , tailTransitionWidth_(d.lookupOrDefault<scalar>("tailTransitionWidth", 0))
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
    , pressureDerivative_(IOobject("phasePressureDerivative",
                                   mesh.time().name(),
                                   mesh,
                                   IOobject::NO_READ,
                                   IOobject::NO_WRITE),
                          mesh,
                          dimensionedScalar(dimDensity / dimTime / dimPressure, 0))
    , alphaSink_(
          IOobject(
              "phaseAlphaSink", mesh.time().name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
          mesh,
          dimensionedScalar(dimDensity / dimTime, 0))
    , energySink_(
          IOobject(
              "phaseEnergySink", mesh.time().name(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
          mesh,
          dimensionedScalar(dimDensity / dimTime, 0))
    , momentumRate_(IOobject("phaseMomentumRate",
                             mesh.time().name(),
                             mesh,
                             IOobject::NO_READ,
                             IOobject::AUTO_WRITE),
                    mesh,
                    dimensionedVector(dimDensity * dimVelocity / dimTime, Zero))
{
    if (!(accommodation_ > 0 && accommodation_ <= 1 && wetThickness_ > 0 && maxFraction_ > 0 &&
          maxFraction_ <= .5 && phaseRampDuration_ >= 0 && tailTransitionWidth_ >= 0))
        FatalIOErrorInFunction(d) << "Invalid kinetic coefficient, wet scale or mass-fraction limit"
                                  << exit(FatalIOError);
    const scalar tail =
        film_.thermo.properties().lookupOrDefault<scalar>("numericalTailThickness", 1e-8);
    if (!(tail >= 0 && tail + tailTransitionWidth_ <= wetThickness_))
        FatalIOErrorInFunction(d) << "Tail transition must lie below the resolved wet-film scale"
                                  << exit(FatalIOError);
    Info << "Phase tail regularization " << mesh.name() << " tail=" << tail
         << " transitionWidth=" << tailTransitionWidth_ << nl;
}

wordList pSpraySmoothFilmPhaseChange::addSupFields() const
{
    return {film_.alpha.name(), film_.thermo.he().name(), film_.U.name()};
}

void pSpraySmoothFilmPhaseChange::requireUpdated() const
{
    if (lastTime_ != mesh().time().timeIndex())
        FatalErrorInFunction << "Paired gas pSpraySmoothVapourFilmSource must be corrected first"
                             << exit(FatalError);
}

void pSpraySmoothFilmPhaseChange::update(const scalarField& allocatedGasVolume)
{
    // Experimental fully coupled constitutive source. Update on every outer
    // iteration; implicit correction terms vanish only when that iteration
    // converges. Independent M/S/E and source-change checks are essential.
    const scalarField previousRate(massRate_.primitiveField());
    lastTime_ = mesh().time().timeIndex();
    massRate_ = dimensionedScalar(massRate_.dimensions(), 0);
    energyRate_ = dimensionedScalar(energyRate_.dimensions(), 0);
    momentumRate_ = dimensionedVector(momentumRate_.dimensions(), Zero);
    pressureDerivative_ = dimensionedScalar(pressureDerivative_.dimensions(), 0);
    alphaSink_ = dimensionedScalar(alphaSink_.dimensions(), 0);
    energySink_ = dimensionedScalar(energySink_.dimensions(), 0);
    const scalar age = mesh().time().value() - phaseStartTime_;
    const scalar r = !phaseEnabled_ || age < 0 ? scalar(0)
                     : phaseRampDuration_ > 0  ? min(age / phaseRampDuration_, scalar(1))
                                               : scalar(1);
    // 相变启动时间和渐启时长只控制初始化安排，仍需独立核验成熟相变阶段。
    const scalar phaseScale = r * r * (3 - 2 * r);
    if (phaseScale <= 0)
    {
        if (film_.pimple.finalIter())
            Info << "Film phase " << mesh().name()
                 << " evaporation kg/s=0 absorption kg/s=0 active=0 limited=0"
                 << " limitedTransfer kg=0 energyToGasW=0 startupScale=0" << nl;
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
    const scalar tail =
        film_.thermo.properties().lookupOrDefault<scalar>("numericalTailThickness", 1e-8);
    scalar evap = 0, absorb = 0, energyToGas = 0, limitedMass = 0, change = 0, absolute = 0;
    label limited = 0, active = 0;
    forAll (fp, i)
    {
        const label cell = fp.faceCells()[i];
        const scalar V = mesh().V()[cell];
        const scalar ml = film_.alpha[cell] * film_.rho[cell] * V;
        const scalar delta = film_.delta[cell];
        if (ml <= 1e-20 || delta < tail)
            continue;
        // Smooth only the numerical dry-tail switch. The original wet-area factor,
        // Schrage law, transport properties and startup schedule remain unchanged
        // for delta >= tail + tailTransitionWidth. A zero width reproduces the
        // original sharp cutoff exactly.
        scalar activation = 1;
        if (tailTransitionWidth_ > 0)
        {
            const scalar s = min(max((delta - tail) / tailTransitionWidth_, scalar(0)), scalar(1));
            activation = s * s * (3 - 2 * s);
        }
        const scalar area =
            phaseScale * ag[i] * min(max(delta / wetThickness_, scalar(0)), scalar(1)) * activation;
        if (area <= 0)
            continue;
        if (mag(hg[i] - (href_ + cpv_ * (tg[i] - Tref_))) > 1e-5)
            FatalErrorInFunction << "Gas enthalpy reference mismatch" << exit(FatalError);
        try
        {
            const auto state = table_.evaluate(tl[cell], X[cell]);
            const scalar K = area * 2 * accommodation_ / (2 - accommodation_) /
                             sqrt(2 * 3.14159265358979323846 * Rv);
            // raw[kg/s]正为蒸发、负为吸收；分别使用液体/气体供体温度。
            const scalar raw = K * (state.pEq / sqrt(tl[cell]) - pg[i] / sqrt(tg[i]));
            const scalar mg = pg[i] / (Rv * tg[i]) * vg[i];
            // 上限同时限制可用液态水和最高X=0.65；下限限制可用蒸气质量。
            scalar upper = maxFraction_ * ml * (1 - X[cell]) / dt;
            if (X[cell] > 0)
                upper = min(upper, max(scalar(0), ml - ml * X[cell] / .65) * (1 - 1e-10) / dt);
            const scalar rate = max(-maxFraction_ * mg / dt, min(upper, raw));
            const bool clipped = mag(rate - raw) > 1e-12 * max(mag(raw), scalar(1e-30));
            const scalar hv = href_ + cpv_ * ((rate >= 0 ? tl[cell] : tg[i]) - Tref_);
            massRate_[cell] = rate / V;
            energyRate_[cell] = rate * hv / V;
            momentumRate_[cell] = rate / V * (rate >= 0 ? film_.U[cell] : ug[i]);
            alphaSink_[cell] = max(rate, scalar(0)) / (V * max(film_.alpha[cell], scalar(1e-30)));
            // 限幅切换点不沿原始通量求Jacobian，以免引入不一致的线性导数。
            if (!clipped)
            {
                pressureDerivative_[cell] = K / (V * sqrt(tg[i]));
                const scalar tmin = X[cell] >= .3 ? 293.15 : 278.15;
                const scalar lo = max(tmin, tl[cell] - .001),
                             hi = min(scalar(353.15), tl[cell] + .001);
                if (!(hi > lo))
                    throw std::runtime_error("No valid thermal Jacobian bracket");
                const auto a = table_.evaluate(lo, X[cell]);
                const auto b = table_.evaluate(hi, X[cell]);
                // 以实际e(T,X)导数将温度敏感性换算为能量方程的隐式系数。
                const scalar dedt = ((b.h - b.pEq / b.rho) - (a.h - a.pEq / a.rho)) / (hi - lo);
                if (!(dedt > 0))
                    throw std::runtime_error("Non-positive liquid energy slope");
                const scalar dmdt = K * (b.pEq / sqrt(hi) - a.pEq / sqrt(lo)) / (hi - lo);
                energySink_[cell] =
                    max(scalar(0), (dmdt * hv + max(rate, scalar(0)) * cpv_) / (V * dedt));
            }
            evap += max(rate, scalar(0));
            absorb += max(-rate, scalar(0));
            energyToGas += rate * hv;
            if (clipped)
            {
                ++limited;
                limitedMass += mag(rate) * dt;
            }
            ++active;
        }
        catch (const std::exception& e)
        {
            FatalErrorInFunction << "Implicit phase cell=" << cell << ", T=" << tl[cell]
                                 << ", X=" << X[cell] << ": " << e.what() << exit(FatalError);
        }
    }
    forAll (previousRate, cell)
    {
        change += mag(massRate_[cell] - previousRate[cell]) * mesh().V()[cell];
        absolute += mag(massRate_[cell]) * mesh().V()[cell];
    }
    change = returnReduce(change, sumOp<scalar>());
    absolute = returnReduce(absolute, sumOp<scalar>());
    Info << "Implicit phase outer " << mesh().name()
         << " sourceRelativeChange=" << change / max(absolute, scalar(1e-30)) << nl;
    if (film_.pimple.finalIter())
        Info << "Film phase " << mesh().name()
             << " evaporation kg/s=" << returnReduce(evap, sumOp<scalar>())
             << " absorption kg/s=" << returnReduce(absorb, sumOp<scalar>())
             << " active=" << returnReduce(active, sumOp<label>())
             << " limited=" << returnReduce(limited, sumOp<label>())
             << " limitedTransfer kg=" << returnReduce(limitedMass, sumOp<scalar>())
             << " energyToGasW=" << returnReduce(energyToGas, sumOp<scalar>())
             << " startupScale=" << phaseScale << nl;
}

void pSpraySmoothFilmPhaseChange::addSup(const volScalarField&,
                                         const volScalarField& alpha,
                                         fvMatrix<scalar>& eqn) const
{
    requireUpdated();
    if (&alpha != &film_.alpha)
        NotImplemented;
    // 显式配对质量源与围绕当前alpha的隐式修正共同组装。
    eqn -= massRate_() + correction(fvm::Sp(alphaSink_(), eqn.psi()));
}

void pSpraySmoothFilmPhaseChange::addSup(const volScalarField&,
                                         const volScalarField&,
                                         const volScalarField& he,
                                         fvMatrix<scalar>& eqn) const
{
    requireUpdated();
    if (&he != &film_.thermo.he())
        NotImplemented;
    eqn -= energyRate_() + correction(fvm::Sp(energySink_(), eqn.psi()));
}

void pSpraySmoothFilmPhaseChange::addSup(const volScalarField&,
                                         const volScalarField&,
                                         const volVectorField&,
                                         fvMatrix<vector>& eqn) const
{
    requireUpdated();
    eqn -= momentumRate_() + correction(fvm::Sp(posPart(massRate_()), eqn.psi()));
}

pSpraySmoothVapourFilmSource::pSpraySmoothVapourFilmSource(const word& name,
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
    , pressureDerivative_(IOobject("filmPressureDerivative",
                                   mesh.time().name(),
                                   mesh,
                                   IOobject::NO_READ,
                                   IOobject::NO_WRITE),
                          mesh,
                          dimensionedScalar(dimDensity / dimTime / dimPressure, 0))
    , momentumRate_(IOobject("filmMomentumRate",
                             mesh.time().name(),
                             mesh,
                             IOobject::NO_READ,
                             IOobject::AUTO_WRITE),
                    mesh,
                    dimensionedVector(dimDensity * dimVelocity / dimTime, Zero))
{
}

wordList pSpraySmoothVapourFilmSource::addSupFields() const
{
    return {"rho", "h", "U"};
}

void pSpraySmoothVapourFilmSource::correct()
{
    lastTime_ = mesh().time().timeIndex();
    massRate_ = dimensionedScalar(massRate_.dimensions(), 0);
    pressureDerivative_ = dimensionedScalar(pressureDerivative_.dimensions(), 0);
    energyRate_ = dimensionedScalar(energyRate_.dimensions(), 0);
    momentumRate_ = dimensionedVector(momentumRate_.dimensions(), Zero);
    DynamicList<pSpraySmoothFilmPhaseChange*> models;
    scalarField sumArea(mesh().nCells(), 0);
    forAll (filmRegions_, ri)
    {
        const auto& fm = mesh().time().lookupObject<fvMesh>(filmRegions_[ri]);
        auto& fvms = fvModels::New(fm);
        pSpraySmoothFilmPhaseChange* found = nullptr;
        forAll (fvms, j)
            if (fvms[j].name() == filmModel_)
                found = &refCast<pSpraySmoothFilmPhaseChange>(fvms[j]);
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
        // 按所有接触膜面的面积分配气相单元体积，避免多膜面重复使用蒸气库存。
        forAll (gp, i)
        {
            const label c = gp.faceCells()[i];
            volume[i] = mesh().V()[c] * gp.magSf()[i] / sumArea[c];
        }
        model.update(volume);
        // Send extensive per-face transfers. Normalise by the mapped GAS face
        // area, not the extruded film wall area, to preserve the paired source.
        const scalarField ag(map.fromNeighbour(gp.magSf()));
        scalarField m(fc.size()), e(fc.size()), kp(fc.size());
        vectorField u(fc.size());
        forAll (fc, i)
        {
            const label c = fc[i];
            const scalar v = model.mesh().V()[c];
            m[i] = model.massRate()[c] * v / ag[i];
            e[i] = model.energyRate()[c] * v / ag[i];
            u[i] = model.momentumRate()[c] * v / ag[i];
            kp[i] = model.pressureDerivative()[c] * v / ag[i];
            totalFilmM += model.massRate()[c] * v;
            totalFilmE += model.energyRate()[c] * v;
            absoluteFilmM += mag(model.massRate()[c] * v);
            absoluteFilmE += mag(model.energyRate()[c] * v);
        }
        const scalarField gm(map.toNeighbour(m)), ge(map.toNeighbour(e)), gkp(map.toNeighbour(kp));
        const vectorField gu(map.toNeighbour(u));
        forAll (gp, i)
        {
            const label c = gp.faceCells()[i];
            const scalar av = gp.magSf()[i] / mesh().V()[c];
            massRate_[c] += gm[i] * av;
            energyRate_[c] += ge[i] * av;
            momentumRate_[c] += gu[i] * av;
            pressureDerivative_[c] += gkp[i] * av;
        }
    }
    totalFilmM = returnReduce(totalFilmM, sumOp<scalar>());
    totalFilmE = returnReduce(totalFilmE, sumOp<scalar>());
    absoluteFilmM = returnReduce(absoluteFilmM, sumOp<scalar>());
    absoluteFilmE = returnReduce(absoluteFilmE, sumOp<scalar>());
    // 累积正/负相变可能抵消；误差阈值以绝对转移总量缩放。
    const scalar errM = gSum(massRate_.primitiveField() * mesh().V()) - totalFilmM;
    const scalar errE = gSum(energyRate_.primitiveField() * mesh().V()) - totalFilmE;
    if (mag(errM) > 1e-10 * max(absoluteFilmM, scalar(1e-15)) ||
        mag(errE) > 1e-10 * max(absoluteFilmE, scalar(1e-10)))
        FatalErrorInFunction << "Nonconservative film mapping: dm/dt=" << errM << " dE/dt=" << errE
                             << exit(FatalError);
    Info << "Paired film-vapour source mapping residual kg/s=" << errM << " W=" << errE << nl;
}

void pSpraySmoothVapourFilmSource::addSup(const volScalarField&, fvMatrix<scalar>& eqn) const
{
    eqn += massRate_();
    // 气相连续性可装配到压力或密度变量；rho分支用p=rho*Rv*T换算导数。
    if (eqn.psi().dimensions() == dimPressure)
        eqn -= correction(fvm::Sp(pressureDerivative_(), eqn.psi()));
    else if (eqn.psi().dimensions() == dimDensity)
    {
        const dimensionedScalar R(
            "Rv", dimPressure / dimDensity / dimTemperature, 8.31446261815324 / .018015268);
        const volScalarField::Internal K(pressureDerivative_() * R *
                                         mesh().lookupObject<volScalarField>("T")());
        eqn -= correction(fvm::Sp(K, eqn.psi()));
    }
    else
        NotImplemented;
}

void pSpraySmoothVapourFilmSource::addSup(const volScalarField&,
                                          const volScalarField& he,
                                          fvMatrix<scalar>& eqn) const
{
    if (he.name() != "h")
        NotImplemented;
    eqn += energyRate_() + correction(fvm::Sp(negPart(massRate_()), eqn.psi()));
}

void pSpraySmoothVapourFilmSource::addSup(const volScalarField&,
                                          const volVectorField&,
                                          fvMatrix<vector>& eqn) const
{
    eqn += momentumRate_() + correction(fvm::Sp(negPart(massRate_()), eqn.psi()));
}
}
}
