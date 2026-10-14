// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file pSprayFilm.C
 * @brief 壁膜盐输运、内能求解和干湿前沿数值处理。
 * 每次物性更新保留alpha*rho表示的质量。尾迹约束、连续性匹配项及
 * 浮点舍入产生的盐/能量贡献明确记录，供外部独立审计。
 */
#include "pSprayFilm.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvModels.H"
#include "fvc.H"
#include "mappedFvPatchBaseBase.H"
#include "fvConstraints.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace solvers
{
defineTypeNameAndDebug(pSprayFilm, 0);
addToRunTimeSelectionTable(solver, pSprayFilm, fvMesh);

// Fourier只处理温度导热；盐扩散携带的偏焓通量由本求解器加入。
pSprayFilm::pSprayFilm(fvMesh& mesh)
    : film(mesh)
    , X_(mesh.lookupObjectRef<volScalarField>("XLiBr"))
    , diffusivity_("mutualDiffusivity", sqr(dimLength) / dimTime, thermo_.properties())
    , table_(std::string(
          thermo_.properties().subDict("mixture").lookup<fileName>("tableFile").c_str()))
{
    if (thermophysicalTransport->type() != "Fourier")
        FatalErrorInFunction
            << "Variable-X LiBr film requires constant/<filmRegion>/thermophysicalTransport: "
            << "laminar { model Fourier; }. Salt-diffusion enthalpy is handled explicitly."
            << exit(FatalError);
    thermo_.he().writeOpt() = IOobject::AUTO_WRITE;
}

void pSprayFilm::prePredictor()
{
    film::prePredictor();
    correctComposition();
}

void pSprayFilm::pressureCorrector()
{
    film::pressureCorrector();
    correctComposition();
    // Final mass flux and salt are now consistent. Re-solve the SAME
    // conservative energy equation (same old-time state and frozen sources),
    // rather than interpreting old-composition energy at the new X.
    thermophysicalPredictor();
}

// 湿-干面扩散通量为零；MPI耦合面也采用同一调和插值规则。
tmp<surfaceScalarField> pSprayFilm::saltDiffusionCoefficient() const
{
    const volScalarField gamma(diffusivity_ * alpha_ * rho);
    auto result = fvc::interpolate(gamma);
    auto& gf = result.ref();
    const auto& w = mesh.weights();
    const auto harmonic = [](const scalar a, const scalar b, const scalar weight)
    {
        // Zero liquid on either side means no connected diffusion path.
        // Reverse-linear reciprocal weighting is the native harmonic rule;
        // use ratios to avoid dividing by zero or underflowing tiny films.
        if (a <= 0 || b <= 0)
            return scalar(0);
        return a < b ? a / ((1 - weight) + weight * (a / b))
                     : b / ((1 - weight) * (b / a) + weight);
    };
    forAll (gf, i)
        gf[i] = harmonic(gamma[mesh.owner()[i]], gamma[mesh.neighbour()[i]], w[i]);
    forAll (gf.boundaryField(), pi)
    {
        if (!gamma.boundaryField()[pi].coupled())
            continue;
        const scalarField own(gamma.boundaryField()[pi].patchInternalField());
        const scalarField nbr(gamma.boundaryField()[pi].patchNeighbourField());
        auto& pf = gf.boundaryFieldRef()[pi];
        forAll (pf, j)
            pf[j] = harmonic(own[j], nbr[j], w.boundaryField()[pi][j]);
    }
    return result;
}

void pSprayFilm::correctComposition()
{
    // M为单位挤出膜体积的液体质量[kg/m³]，不是整单元质量。
    // 盐对流使用原生膜alphaRhoPhi，与液体质量使用同一通量。
    // Both the continuity predictor and final conservative film mass flux
    // require an accompanying salt update before thermophysical inversion.
    const scalarField M(alpha_.primitiveField() * rho.primitiveField());
    const scalarField previousX(X_.primitiveField());
    // Relative to this region's inventory, even one affected tail per cell
    // remains far below the independent salt-conservation audit tolerance.
    const scalar saltRoundoffMass = max(scalar(1e-20), 1e-12 * gSum(mag(M) * mesh.V()));
    fvScalarMatrix saltEqn(
        fvm::ddt(alpha_, rho, X_) +
            fvm::div(alphaRhoPhi_, X_)
            // Match the native film momentum/energy treatment of continuity
            // truncation and alpha limiting. This preserves a uniform passive
            // concentration at wet/dry fronts instead of dividing two differently
            // converged tiny inventories. Its extra salt is independently audited.
            - fvm::Sp(contErr(), X_) - fvm::laplacian(saltDiffusionCoefficient(), X_) ==
        fvModels().source(alpha_, rho, X_));
    // Exactly dry disconnected rows carry no mass and no salt: retain their
    // reference composition without adding a fictitious mass source.
    forAll (M, i)
    {
        if (saltEqn.diag()[i] == 0)
        {
            if (mag(saltEqn.source()[i]) > 1e-25)
                FatalErrorInFunction << "Salt source in a zero-mass isolated film cell " << i
                                     << exit(FatalError);
            saltEqn.diag()[i] = 1;
            saltEqn.source()[i] = previousX[i];
        }
    }
    volScalarField::Internal rowMagnitude(
        IOobject("saltEquationRowMagnitude", runTime.name(), mesh),
        mesh,
        dimensionedScalar(dimless, 1));
    const auto rowNorm = [&]()
    {
        scalarField norm(mag(saltEqn.D()));
        if (saltEqn.hasUpper())
        {
            const auto& upper = saltEqn.upper();
            const auto& lower = saltEqn.lower();
            forAll (upper, fi)
            {
                norm[mesh.owner()[fi]] += mag(upper[fi]);
                norm[mesh.neighbour()[fi]] += mag(lower[fi]);
            }
        }
        forAll (mesh.boundary(), pi)
        {
            if (!mesh.boundary()[pi].coupled())
                continue;
            const auto& bc = saltEqn.boundaryCoeffs()[pi];
            const auto& fc = mesh.boundary()[pi].faceCells();
            forAll (bc, j)
                norm[fc[j]] += mag(bc[j]);
        }
        return max(norm, mag(saltEqn.source()))();
    };
    const scalarField physicalRowNorm(rowNorm());
    DynamicList<label> saltTailCells;
    DynamicList<scalar> saltTailValues;
    scalar tailMass = 0;
    forAll (M, i)
    {
        if (mag(M[i] * mesh.V()[i]) < saltRoundoffMass &&
            physicalRowNorm[i] * runTime.deltaTValue() < saltRoundoffMass)
        {
            saltTailCells.append(i);
            saltTailValues.append(previousX[i]);
            tailMass += mag(M[i] * mesh.V()[i]);
        }
    }
    // 保存约束前方程，用其残差计算尾迹约束额外产生的盐速率。
    const fvScalarMatrix physicalSaltEquation(saltEqn);
    saltEqn.setValues(labelList(saltTailCells), scalarField(saltTailValues));
    rowMagnitude.primitiveFieldRef() = max(rowNorm(), scalar(1e-300));
    tailSaltMass_ = returnReduce(tailMass, sumOp<scalar>());
    Info << "Salt matrix diagonal range before scaling = " << gMin(saltEqn.diag()) << " "
         << gMax(saltEqn.diag()) << nl;
    // Algebraically equivalent row division avoids reciprocal overflow in
    // vanishing-film cells without adding any fictitious liquid or salt.
    saltEqn /= rowMagnitude;
    saltEqn.solve();
    const scalarField saltResidual(physicalSaltEquation.residual());
    scalar saltPower = 0;
    forAll (saltTailCells, j)
        saltPower -= saltResidual[saltTailCells[j]];
    tailSaltRate_ = returnReduce(saltPower, sumOp<scalar>());
    // 有实际库存的越界浓度报错；舍入误差和极小尾迹的有界处理单独记账。
    scalar roundoffSalt = 0, tailSaltCorrection = 0, tailLiquidMass = 0;
    forAll (X_, i)
    {
        const scalar x = X_[i];
        const scalar mass = M[i] * mesh.V()[i];
        const bool outside = x < -1e-12 || x > .65 + 1e-12;
        if (!std::isfinite(x) || (outside && mag(mass) >= saltRoundoffMass))
            FatalErrorInFunction << "Nonvolatile salt fraction outside table: X=" << x
                                 << ", cell=" << i << ", liquid mass kg=" << mass
                                 << exit(FatalError);
        X_[i] = max(scalar(0), min(scalar(.65), x));
        if (outside)
        {
            tailSaltCorrection += (X_[i] - x) * mass;
            tailLiquidMass += mag(mass);
        }
        else
            roundoffSalt += (X_[i] - x) * mass;
    }
    X_.correctBoundaryConditions();
    const scalarField saltContinuityRate(contErr().primitiveField() * X_.primitiveField() *
                                         mesh.V());
    continuitySaltRate_ = gSum(saltContinuityRate);
    continuitySaltAbsoluteRate_ = gSum(mag(saltContinuityRate));
    auto& tail = mesh.lookupObjectRef<volScalarField>("pSprayConductivityTail");
    const scalar tailThickness =
        thermo_.properties().lookupOrDefault<scalar>("numericalTailThickness", 1e-8);
    forAll (tail, i)
        tail[i] = alpha_[i] * VbyA_[i] < tailThickness ? 1 : 0;
    tail.correctBoundaryConditions();
    auto& he = thermo_.he();
    // Preserve physical time-n energy. The following is ONLY an initial
    // nonlinear iterate for the next conservative energy solve, not a heat
    // source or a temperature clamp on the converged solution.
    he.oldTime();
    auto& T = thermo_.T();
    T.oldTime();
    forAll (he, i)
    {
        const scalar lower = X_[i] >= .3 && tail[i] < .5 ? scalar(293.15) : scalar(278.15);
        // A newly wetted cell can inherit an invalid dry-state temperature.
        // Start that nonlinear iteration inside the valid interval, so a
        // coupled boundary does not fall outside it before the first solve.
        const scalar guess = T[i] < lower ? scalar(313.15) : min(T[i], scalar(353.15));
        const auto state = table_.evaluate(guess, X_[i]);
        // Supply a consistent T as well as e. At the conductivity-domain
        // endpoint, inverting e from a stale colder T can undershoot by
        // roundoff before the energy equation has even been assembled.
        T[i] = guess;
        he[i] = state.h - state.pEq / state.rho;
    }
    forAll (he.boundaryField(), patchi)
    {
        auto& pT = T.boundaryFieldRef()[patchi];
        scalarField guess(pT);
        const auto& xb = X_.boundaryField()[patchi];
        // RhoFluidThermo::calculate treats mixed/coupledTemperature as
        // fixesValue(): it recomputes e from pT and ignores an e-only guess.
        // Refresh their current iterate too. A prescribed fixed temperature
        // remains the physical boundary value, including its domain check.
        if (pT.type() != "fixedValue")
        {
            forAll (guess, j)
            {
                const scalar lower = xb[j] >= .3 && tail.boundaryField()[patchi][j] < .5
                                         ? scalar(293.15)
                                         : scalar(278.15);
                guess[j] = guess[j] < lower ? scalar(313.15) : min(guess[j], scalar(353.15));
            }
            pT == guess;
        }
        he.boundaryFieldRef()[patchi] == thermo_.he(guess, patchi);
    }
    thermo_.correct();
    // 新rho对应更新后的组分/迭代温度，按原M重算alpha保持质量。
    forAll (M, i)
        alpha_[i] = M[i] / rho[i];
    alpha_.correctBoundaryConditions();
    delta_ = alpha_ * VbyA_;
    delta_.correctBoundaryConditions();
    Info << "LiBr film X min/max = " << gMin(X_) << " " << gMax(X_)
         << "; roundoff salt correction kg = " << returnReduce(roundoffSalt, sumOp<scalar>())
         << "; numerical-tail salt correction kg="
         << returnReduce(tailSaltCorrection, sumOp<scalar>())
         << "; affected liquid mass kg=" << returnReduce(tailLiquidMass, sumOp<scalar>())
         << "; salt roundoff mass threshold kg=" << saltRoundoffMass << nl;
}

void pSprayFilm::thermophysicalPredictor()
{
    // he在本模型中为比内能e[J/kg]；物性修正前保存液体质量密度。
    volScalarField& he = thermo_.he();
    const scalarField oldReferenceEnergy(he.primitiveField());
    const scalarField conservedMassDensity(alpha.primitiveField() * rho.primitiveField());
    volScalarField hX(IOobject("diffusivePartialEnthalpy", runTime.name(), mesh),
                      mesh,
                      dimensionedScalar(dimEnergy / dimMass, 0),
                      "zeroGradient");
    // hX=∂h(p,T,X)/∂X，含表焓、平衡压力和密度的组分导数。
    // 乘盐扩散质量通量后得到其携带的能量通量[W]。
    forAll (hX, i)
    {
        const scalar x = X_[i], t = thermo_.T()[i];
        const auto st = table_.evaluate(t, x);
        const auto pe = table_.pressureData(t, x);
        const scalar xp = min(x + 1e-6, scalar(.65)), xm = max(x - 1e-6, scalar(0));
        const scalar rhoX = (table_.density(t, xp) - table_.density(t, xm)) / (xp - xm);
        hX[i] = st.dh_dX - st.pEq * pe[2] / st.rho - (p[i] - st.pEq) * rhoX / sqr(st.rho);
    }
    hX.correctBoundaryConditions();
    const surfaceScalarField diffusiveEnergyFlux("diffusiveEnergyFlux",
                                                 -saltDiffusionCoefficient() * fvc::snGrad(X_) *
                                                     mesh.magSf() * fvc::interpolate(hX));
    // Fourier handles -k grad(T). Nonvolatile mutual diffusion separately
    // carries the partial enthalpy difference h_X times the LiBr mass flux.
    fvScalarMatrix eqn(fvm::ddt(alpha, rho, he) + fvm::div(alphaRhoPhi, he) +
                           fvc::div(diffusiveEnergyFlux) - fvm::Sp(contErr(), he) +
                           thermophysicalTransport->divq(he) ==
                       fvModels().source(alpha, rho, he));
    eqn.relax();
    fvConstraints().constrain(eqn);
    scalarField effectiveDiagonal(eqn.D());
    DynamicList<label> unoccupied;
    DynamicList<scalar> referenceValues;
    const scalar dryMassFloor = 1e-20;
    const scalar rateFloor = dryMassFloor / runTime.deltaTValue();
    scalar constrainedLiquidMass = 0, constrainedMaxDelta = 0;
    forAll (effectiveDiagonal, i)
    {
        const scalar mass = alpha[i] * rho[i] * mesh.V()[i];
        if (mag(mass) < dryMassFloor && mag(effectiveDiagonal[i]) < rateFloor &&
            mag(eqn.source()[i]) < rateFloor * 1e6)
        {
            unoccupied.append(i);
            referenceValues.append(oldReferenceEnergy[i]);
            constrainedLiquidMass += mag(mass);
            constrainedMaxDelta = max(constrainedMaxDelta, delta_[i]);
        }
    }
    // A specific temperature is undefined where the liquid inventory is at
    // floating-point numerical tails. Retain its reference value only if
    // mass, thermal coupling and heat source are all below explicit bounds.
    // 使用约束前方程残差核算尾迹参考值引入的人工功率[W]。
    const fvScalarMatrix physicalEquation(eqn);
    eqn.setValues(labelList(unoccupied), scalarField(referenceValues));
    effectiveDiagonal = eqn.D();
    tailLiquidMass_ = returnReduce(constrainedLiquidMass, sumOp<scalar>());
    Info << "Energy reference constrained numerical-tail liquid mass kg = " << tailLiquidMass_
         << "; maxDelta=" << returnReduce(constrainedMaxDelta, maxOp<scalar>())
         << "; massThreshold=" << dryMassFloor << nl;
    forAll (effectiveDiagonal, i)
    {
        if (effectiveDiagonal[i] == 0)
        {
            if (mag(eqn.source()[i]) > 1e-20)
                FatalErrorInFunction << "Heat source in an isolated zero-mass adiabatic film cell "
                                     << i << exit(FatalError);
            eqn.diag()[i] = 1;
            eqn.source()[i] = oldReferenceEnergy[i];
            effectiveDiagonal[i] = 1;
        }
    }
    volScalarField::Internal rowMagnitude(
        IOobject("energyEquationRowMagnitude", runTime.name(), mesh),
        mesh,
        dimensionedScalar(dimless, 1));
    rowMagnitude.primitiveFieldRef() = max(mag(effectiveDiagonal), scalar(1e-300));
    eqn /= rowMagnitude;
    eqn.solve();
    fvConstraints().constrain(he);
    const scalarField physicalResidual(physicalEquation.residual());
    scalar artificialPower = 0;
    forAll (unoccupied, j)
        artificialPower -= physicalResidual[unoccupied[j]];
    tailEnergyPower_ = returnReduce(artificialPower, sumOp<scalar>());
    // 检查求解后真实浓溶液温度；此处输出诊断，不截断温度。
    const auto predictedT = thermo_.The(he, p, thermo_.T());
    const auto& tail = mesh.lookupObject<volScalarField>("pSprayConductivityTail");
    forAll (he, i)
        if (X_[i] >= .3 && predictedT()[i] < 293.15 && tail[i] < .5)
            Pout << "PSPRAY_THERMO_DOMAIN cell=" << i << " C=" << mesh.C()[i] << " X=" << X_[i]
                 << " Tsolve=" << predictedT()[i] << " Titerate=" << thermo_.T()[i]
                 << " e=" << he[i] << " mass=" << alpha[i] * rho[i] * mesh.V()[i]
                 << " diag=" << effectiveDiagonal[i] << " source=" << eqn.source()[i] << endl;
    forAll (he.boundaryField(), pi)
    {
        const auto& tb = thermo_.T().boundaryField()[pi];
        const auto& xb = X_.boundaryField()[pi];
        forAll (xb, j)
        {
            const scalar actualT = tb.fixesValue() ? tb[j] : predictedT().boundaryField()[pi][j];
            if (xb[j] >= .3 && actualT < 293.15 && tail.boundaryField()[pi][j] < .5)
                Pout << "PSPRAY_THERMO_DOMAIN patch=" << mesh.boundary()[pi].name() << " face=" << j
                     << " X=" << xb[j] << " Tsolve=" << actualT << " fixesValue=" << tb.fixesValue()
                     << " massCell="
                     << alpha[tb.patch().faceCells()[j]] * rho[tb.patch().faceCells()[j]] *
                            mesh.V()[tb.patch().faceCells()[j]]
                     << endl;
        }
    }
    thermo_.correct();
    forAll (alpha_, i)
        alpha_[i] = conservedMassDensity[i] / rho[i];
    alpha_.correctBoundaryConditions();
    delta_ = alpha_ * VbyA_;
    delta_.correctBoundaryConditions();
}

void pSprayFilm::postSolve()
{
    // 广延库存单位为kg、kg、J；Qout正为向域外放热。
    // gSum/returnReduce跨MPI分区求总量，后处理不能重复累加同一总量。
    Info << "PSPRAY_CONTINUITY_SALT time=" << runTime.value() << " region=" << mesh.name()
         << " artificialSaltRate=" << continuitySaltRate_
         << " absoluteSaltRate=" << continuitySaltAbsoluteRate_ << nl;

    film::postSolve();
    Info << "PSPRAY_TAIL_SALT time=" << runTime.value() << " region=" << mesh.name()
         << " artificialSaltRate=" << tailSaltRate_ << " constrainedLiquidMass=" << tailSaltMass_
         << nl;
    Info << "PSPRAY_TAIL_ENERGY time=" << runTime.value() << " region=" << mesh.name()
         << " artificialPower=" << tailEnergyPower_ << " constrainedLiquidMass=" << tailLiquidMass_
         << nl;
    const scalarField mass(alpha.primitiveField() * rho.primitiveField() * mesh.V());
    const auto& tail = mesh.lookupObject<volScalarField>("pSprayConductivityTail");
    scalar referenceMass = 0, referenceMaxDelta = 0;
    label referenceCells = 0;
    forAll (mass, i)
    {
        if (tail[i] > .5 && X_[i] >= .3 && thermo_.T()[i] < 293.15)
        {
            referenceMass += mass[i];
            referenceMaxDelta = max(referenceMaxDelta, delta_[i]);
            ++referenceCells;
        }
    }
    Info << "PSPRAY_TRANSPORT_REFERENCE time=" << runTime.value() << " region=" << mesh.name()
         << " mass=" << returnReduce(referenceMass, sumOp<scalar>())
         << " maxDelta=" << returnReduce(referenceMaxDelta, maxOp<scalar>())
         << " cells=" << returnReduce(referenceCells, sumOp<label>()) << nl;
    Info << "PSPRAY_BUDGET time=" << runTime.value() << " region=" << mesh.name()
         << " mass=" << gSum(mass) << " salt=" << gSum(mass * X_.primitiveField())
         << " E=" << gSum(mass * thermo_.he().primitiveField()) << " maxDelta=" << gMax(delta_)
         << nl;
    const auto tq = thermophysicalTransport->q();
    forAll (mesh.boundary(), pi)
    {
        if (mesh.boundary()[pi].coupled() || mesh.boundary()[pi].type() == "empty")
            continue;
        Info << "PSPRAY_HEAT time=" << runTime.value() << " region=" << mesh.name()
             << " patch=" << mesh.boundary()[pi].name()
             << " Qout=" << gSum(tq().boundaryField()[pi] * mesh.magSf().boundaryField()[pi])
             << " area=" << gSum(mesh.magSf().boundaryField()[pi]) << nl;
    }
    if (mesh.name() == "tubeFilm" || mesh.name() == "casingFilm")
    {
        // Spatial diagnostics only. Bottom bands overlap their bank totals.
        const auto& centres = mesh.C().primitiveField();
        const scalar bottomY = gMin(centres.component(vector::Y)) + .01;
        scalarField wallHeat(mesh.nCells(), 0);
        forAll (mesh.boundary(), pi)
        {
            if (mesh.boundary()[pi].name() != "wall")
                continue;
            const auto& fc = mesh.boundary()[pi].faceCells();
            const auto& q = tq().boundaryField()[pi];
            const auto& a = mesh.magSf().boundaryField()[pi];
            forAll (fc, j)
                wallHeat[fc[j]] += q[j] * a[j];
        }
        const volScalarField* phaseM = mesh.foundObject<volScalarField>("phaseMassRate")
                                           ? &mesh.lookupObject<volScalarField>("phaseMassRate")
                                           : nullptr;
        const volScalarField* phaseE = mesh.foundObject<volScalarField>("phaseEnergyRate")
                                           ? &mesh.lookupObject<volScalarField>("phaseEnergyRate")
                                           : nullptr;
        const bool physicalPools =
            mesh.name() == "casingFilm" && thermo_.properties().found("poolGeometry");
        const dictionary pool(physicalPools ? thermo_.properties().subDict("poolGeometry")
                                            : dictionary::null);
        const scalar slope = pool.lookupOrDefault<scalar>("slope", 0);
        const scalar intercept = pool.lookupOrDefault<scalar>("intercept", 0);
        const scalar verticalX = pool.lookupOrDefault<scalar>("verticalX", 0);
        const scalar dividerTopY = pool.lookupOrDefault<scalar>("dividerTopY", 0);
        const scalar leftFloorY = pool.lookupOrDefault<scalar>("leftFloorY", 0);
        const scalar rightFloorY = pool.lookupOrDefault<scalar>("rightFloorY", 0);
        const scalar bandHeight = pool.lookupOrDefault<scalar>("bandHeight", .01);
        const label groups = physicalPools ? 5 : (mesh.name() == "casingFilm" ? 4 : 2);
        for (label zi = 0; zi < groups; ++zi)
        {
            const bool left = physicalPools ? (zi == 0 || zi == 3) : (zi % 2 == 0);
            const bool bottom = physicalPools ? zi >= 3 : zi >= 2;
            const word zone = physicalPools
                                  ? word(zi == 0   ? "leftLowerChamber"
                                         : zi == 1 ? "rightLowerChamber"
                                         : zi == 2 ? "upperCasing"
                                         : zi == 3 ? "leftBottom10mm"
                                                   : "rightBottom10mm")
                                  : (bottom ? (left ? "leftBottom10mm" : "rightBottom10mm")
                                            : (left ? "left" : "right"));
            const scalar zoneBottomY =
                physicalPools ? (left ? leftFloorY : rightFloorY) + bandHeight : bottomY;
            scalar m = 0, s = 0, e = 0, v = 0, a = 0, aw = 0, d = 0, mt = 0, qe = 0, ab = 0, ev = 0,
                   qw = 0;
            forAll (mass, i)
            {
                if (physicalPools)
                {
                    const bool lower = centres[i].y() < dividerTopY;
                    const scalar dividerX = min(verticalX, slope * centres[i].y() + intercept);
                    const bool onLeft = centres[i].x() < dividerX;
                    if (zi == 2 ? lower : (!lower || onLeft != left))
                        continue;
                    if (bottom && centres[i].y() > zoneBottomY)
                        continue;
                }
                else if ((centres[i].x() < 0) != left || (bottom && centres[i].y() > bottomY))
                    continue;
                const scalar area = mesh.V()[i] / VbyA_[i];
                m += mass[i];
                s += mass[i] * X_[i];
                e += mass[i] * thermo_.he()[i];
                v += alpha[i] * mesh.V()[i];
                a += area;
                if (delta_[i] >= 1e-6)
                    aw += area;
                d = max(d, delta_[i]);
                mt += mass[i] * thermo_.T()[i];
                qw += wallHeat[i];
                if (phaseM)
                {
                    const scalar rate = (*phaseM)[i] * mesh.V()[i];
                    ev += max(rate, scalar(0));
                    ab += max(-rate, scalar(0));
                    qe += (*phaseE)[i] * mesh.V()[i];
                }
            }
            m = returnReduce(m, sumOp<scalar>());
            s = returnReduce(s, sumOp<scalar>());
            e = returnReduce(e, sumOp<scalar>());
            v = returnReduce(v, sumOp<scalar>());
            a = returnReduce(a, sumOp<scalar>());
            aw = returnReduce(aw, sumOp<scalar>());
            d = returnReduce(d, maxOp<scalar>());
            mt = returnReduce(mt, sumOp<scalar>());
            ev = returnReduce(ev, sumOp<scalar>());
            ab = returnReduce(ab, sumOp<scalar>());
            qe = returnReduce(qe, sumOp<scalar>());
            qw = returnReduce(qw, sumOp<scalar>());
            Info << "PSPRAY_ZONE time=" << runTime.value() << " region=" << mesh.name()
                 << " zone=" << zone << " mass=" << m << " salt=" << s << " E=" << e
                 << " volume=" << v << " area=" << a << " wetAreaAbove1um=" << aw
                 << " maxDelta=" << d << " massWeightedT=" << (m > 0 ? mt / m : 0)
                 << " evaporation=" << ev << " absorption=" << ab << " phaseEnergyToGas=" << qe
                 << " wallQout=" << qw << " bottomBandYmax=" << zoneBottomY
                 << " physicalPoolGeometry=" << physicalPools << nl;
        }
    }

    // The first film region also reports the common carrier budget.
    if (mesh.name() != "film" && mesh.name() != "tubeFilm")
        return;
    const auto& gm = surfacePatchMap().nbrMesh();
    const auto& gt = gm.lookupObject<fluidThermo>(physicalProperties::typeName);
    const auto tr = gt.rho();
    const auto& gr = tr();
    const auto& gu = gm.lookupObject<volVectorField>("U");
    const auto& gp = gm.lookupObject<volScalarField>("p");
    const auto& gh = gt.he();
    const auto& phi = gm.lookupObject<surfaceScalarField>("phi");
    const auto& gg = gm.lookupObject<uniformDimensionedVectorField>("g");
    const scalarField gasMass(gr.primitiveField() * gm.V());
    const scalar gasE = gSum(gasMass * (gh.primitiveField() + .5 * magSqr(gu.primitiveField())) -
                             gp.primitiveField() * gm.V());
    Info << "PSPRAY_GAS time=" << runTime.value() << " mass=" << gSum(gasMass) << " E=" << gasE
         << " gravityPower=" << gSum(gasMass * (gu.primitiveField() & gg.value())) << nl;
    forAll (gm.boundary(), pi)
    {
        const auto& patch = gm.boundary()[pi];
        if (patch.coupled() || patch.type() == "empty")
            continue;
        // Carrier is explicitly pure H2O, laminar; its hConst enthalpy
        // gradient and temperature gradient produce identical Fourier heat.
        const scalarField qout(-gt.kappa().boundaryField()[pi] *
                               gt.T().boundaryField()[pi].snGrad());
        Info << "PSPRAY_GAS_PATCH time=" << runTime.value() << " patch=" << patch.name()
             << " massOut=" << gSum(phi.boundaryField()[pi]) << " energyOut="
             << gSum(phi.boundaryField()[pi] *
                     (gh.boundaryField()[pi] + .5 * magSqr(gu.boundaryField()[pi])))
             << " Qout=" << gSum(qout * gm.magSf().boundaryField()[pi])
             << " area=" << gSum(gm.magSf().boundaryField()[pi]) << nl;
    }
}
}
}
