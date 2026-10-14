// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file pSprayConservativeTemperatureFvPatchScalarField.C
 * @brief 修正映射温度边界的面积基准，保持气膜界面总热量守恒。
 * 这里处理气膜接触界面；传热管固体侧wall的codedMixed在算例中配置。
 */
#include "pSprayConservativeTemperatureFvPatchScalarField.H"
#include "mappedFvPatchBaseBase.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
pSprayConservativeTemperatureFvPatchScalarField::pSprayConservativeTemperatureFvPatchScalarField(
    const fvPatch& patch, const DimensionedField<scalar, fvMesh>& field, const dictionary& dict)
    : coupledTemperatureFvPatchScalarField(patch, field, dict)
    , neighbourTemperatureName_(dict.lookupOrDefault<word>("Tnbr", "T"))
{
}

pSprayConservativeTemperatureFvPatchScalarField::pSprayConservativeTemperatureFvPatchScalarField(
    const pSprayConservativeTemperatureFvPatchScalarField& other,
    const fvPatch& patch,
    const DimensionedField<scalar, fvMesh>& field,
    const fieldMapper& mapper)
    : coupledTemperatureFvPatchScalarField(other, patch, field, mapper)
    , neighbourTemperatureName_(other.neighbourTemperatureName_)
{
}

pSprayConservativeTemperatureFvPatchScalarField::pSprayConservativeTemperatureFvPatchScalarField(
    const pSprayConservativeTemperatureFvPatchScalarField& other,
    const DimensionedField<scalar, fvMesh>& field)
    : coupledTemperatureFvPatchScalarField(other, field)
    , neighbourTemperatureName_(other.neighbourTemperatureName_)
{
}

void pSprayConservativeTemperatureFvPatchScalarField::getNbr(tmp<scalarField>& sumKappaTByDeltaNbr,
                                                             tmp<scalarField>& sumKappaByDeltaNbr,
                                                             tmp<scalarField>& qNbr) const
{
    coupledTemperatureFvPatchScalarField::getNbr(sumKappaTByDeltaNbr, sumKappaByDeltaNbr, qNbr);

    const mappedFvPatchBaseBase& mapper = mappedFvPatchBaseBase::getMap(patch());
    const fvPatch& neighbour = mapper.nbrFvPatch();
    const auto& neighbourTemperature =
        neighbour.lookupPatchField<volScalarField, scalar>(neighbourTemperatureName_);
    if (!isA<pSprayConservativeTemperatureFvPatchScalarField>(neighbourTemperature))
    {
        FatalErrorInFunction << "Both sides of " << patch().name()
                             << " must use pSprayConservativeTemperature" << exit(FatalError);
    }

    // This case has one-to-one mapped gas/film faces. The extruded film's
    // outer face is larger than its gas counterpart on curved surfaces.
    // Export conductance and correction in W/K and W, divided by the
    // receiving face area. Both sides then share one interface temperature
    // and exchange equal, opposite integrated heat at convergence.
    const scalarField receivingArea(mapper.fromNeighbour(neighbour.magSf()));
    forAll (receivingArea, face)
    {
        if (!(receivingArea[face] > VSMALL))
        {
            FatalErrorInFunction << "Invalid mapped receiving area on " << patch().name()
                                 << " face " << face << exit(FatalError);
        }
    }
    // 每面积热导/修正通量乘A发送/A接收，三项必须使用同一面积比例。
    const scalarField areaRatio(patch().magSf() / receivingArea);
    sumKappaTByDeltaNbr = areaRatio * sumKappaTByDeltaNbr();
    sumKappaByDeltaNbr = areaRatio * sumKappaByDeltaNbr();
    if (qNbr.valid())
    {
        qNbr = areaRatio * qNbr();
    }
}

makePatchTypeField(fvPatchScalarField, pSprayConservativeTemperatureFvPatchScalarField);
}
