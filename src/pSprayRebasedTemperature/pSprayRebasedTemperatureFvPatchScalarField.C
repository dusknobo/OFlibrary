// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file pSprayRebasedTemperatureFvPatchScalarField.C
 * @brief mixed温度边界的等价线性表示。
 * Tb=f*Tr+(1-f)*(Ti+g/d)，d为deltaCoeffs[1/m]，g为法向温度梯度[K/m]。
 * 配对调整Tr与g，使任意Ti下的真实Tb及法向梯度保持不变。
 */
#include "pSprayRebasedTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "Time.H"
#include <cmath>

namespace Foam
{
pSprayRebasedTemperatureFvPatchScalarField::pSprayRebasedTemperatureFvPatchScalarField(
    const fvPatch& p, const DimensionedField<scalar, fvMesh>& field, const dictionary& dict)
    : pSprayConservativeTemperatureFvPatchScalarField(p, field, dict)
    , referenceRebase_(dict.lookupOrDefault<Switch>("referenceRebase", false))
    , referenceMin_(dict.lookupOrDefault<scalar>("referenceMin", 278.15))
    , referenceMax_(dict.lookupOrDefault<scalar>("referenceMax", 353.15))
{
    if (!(referenceMin_ < referenceMax_))
        FatalIOErrorInFunction(dict) << "Invalid reference-temperature range" << exit(FatalIOError);
}

pSprayRebasedTemperatureFvPatchScalarField::pSprayRebasedTemperatureFvPatchScalarField(
    const pSprayRebasedTemperatureFvPatchScalarField& other,
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& field,
    const fieldMapper& mapper)
    : pSprayConservativeTemperatureFvPatchScalarField(other, p, field, mapper)
    , referenceRebase_(other.referenceRebase_)
    , referenceMin_(other.referenceMin_)
    , referenceMax_(other.referenceMax_)
{
}

pSprayRebasedTemperatureFvPatchScalarField::pSprayRebasedTemperatureFvPatchScalarField(
    const pSprayRebasedTemperatureFvPatchScalarField& other,
    const DimensionedField<scalar, fvMesh>& field)
    : pSprayConservativeTemperatureFvPatchScalarField(other, field)
    , referenceRebase_(other.referenceRebase_)
    , referenceMin_(other.referenceMin_)
    , referenceMax_(other.referenceMax_)
{
}

void pSprayRebasedTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
        return;
    pSprayConservativeTemperatureFvPatchScalarField::updateCoeffs();
    const scalarField inside(patchInternalField());
    const scalarField& d = patch().deltaCoeffs();
    bool reported = false;
    forAll (refValue(), i)
    {
        const scalar tr = refValue()[i];
        if (tr >= referenceMin_ && tr <= referenceMax_)
            continue;
        // 原系数求得的真实温度和梯度作为等价变换检查基准。
        const scalar f = valueFraction()[i], g = refGrad()[i], ti = inside[i];
        const scalar tb = f * tr + (1 - f) * (ti + g / d[i]);
        const scalar normalGradient = f * d[i] * (tr - ti) + (1 - f) * g;
        if (!reported)
        {
            Pout << "PSPRAY_T_REFERENCE time=" << db().time().value()
                 << " region=" << patch().mesh().name() << " patch=" << patch().name()
                 << " face=" << i << " cell=" << patch().faceCells()[i] << " Cf=" << patch().Cf()[i]
                 << " refValue=" << tr << " actualMixedT=" << tb << " cellT=" << ti
                 << " valueFraction=" << f << " refGradient=" << g << " deltaCoeffs=" << d[i]
                 << " rebase=" << referenceRebase_ << endl;
            reported = true;
        }
        // Diagnostic mode preserves every coefficient and lets native thermo
        // produce the original error and stack.
        if (!referenceRebase_)
            continue;
        // 真实温度越界或f≈1使梯度表达奇异时，报错而不修补物理状态。
        if (!std::isfinite(tb) || tb < referenceMin_ || tb > referenceMax_ || !(f >= 0 && f <= 1) ||
            !(d[i] > 0) || !(1 - f > 1e-12))
            FatalErrorInFunction << "Cannot rebase a mixed reference while retaining a valid "
                                 << "physical interface temperature: patch=" << patch().name()
                                 << " face=" << i << " reference=" << tr << " actualMixedT=" << tb
                                 << " fraction=" << f << exit(FatalError);
        // The mixed temperature and its normal gradient are invariant for
        // EVERY owner-cell temperature, not only the current iterate:
        // f*Tr + (1-f)*(Ti+g/d)
        // = f*TrNew + (1-f)*(Ti+gNew/d).
        // Do not clamp a physical T, energy, heat flux, mass or concentration.
        refValue()[i] = tb;
        refGrad()[i] = g + f * d[i] * (tr - tb) / (1 - f);
        // 同时检查温度与梯度；只匹配当前Tb不足以证明传热算子不变。
        const scalar checkedT = f * refValue()[i] + (1 - f) * (ti + refGrad()[i] / d[i]);
        const scalar checkedGradient = f * d[i] * (refValue()[i] - ti) + (1 - f) * refGrad()[i];
        if (mag(checkedT - tb) > 1e-10 * max(scalar(1), mag(tb)) ||
            mag(checkedGradient - normalGradient) > 1e-10 * max(scalar(1), mag(normalGradient)))
            FatalErrorInFunction << "Mixed-temperature reparameterisation changed temperature "
                                 << "or normal heat-gradient operator on face " << i
                                 << exit(FatalError);
    }
}

void pSprayRebasedTemperatureFvPatchScalarField::write(Ostream& os) const
{
    pSprayConservativeTemperatureFvPatchScalarField::write(os);
    writeEntry(os, "referenceRebase", referenceRebase_);
    writeEntry(os, "referenceMin", referenceMin_);
    writeEntry(os, "referenceMax", referenceMax_);
}

makePatchTypeField(fvPatchScalarField, pSprayRebasedTemperatureFvPatchScalarField);
}
