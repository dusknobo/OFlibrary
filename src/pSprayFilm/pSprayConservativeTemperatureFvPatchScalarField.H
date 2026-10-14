// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file pSprayConservativeTemperatureFvPatchScalarField.H
 * @brief 当前一一对应气膜面的面积守恒温度边界。
 * 曲面挤出膜与气相面面积不同；按接收面积换算总热导[W/K]与热量[W]。
 * Tnbr指定对侧温度字段，两侧采用本类型或其派生类型。
 */
#ifndef pSprayConservativeTemperatureFvPatchScalarField_H
#define pSprayConservativeTemperatureFvPatchScalarField_H
#include "coupledTemperatureFvPatchScalarField.H"

namespace Foam
{
class pSprayConservativeTemperatureFvPatchScalarField : public coupledTemperatureFvPatchScalarField
{
    const word neighbourTemperatureName_;

protected:
    /// 先获取原生热耦合数据，再按发送/接收面面积换算。
    virtual void getNbr(tmp<scalarField>& sumKappaTByDeltaNbr,
                        tmp<scalarField>& sumKappaByDeltaNbr,
                        tmp<scalarField>& qNbr) const override;

public:
    TypeName("pSprayConservativeTemperature");

    pSprayConservativeTemperatureFvPatchScalarField(const fvPatch&,
                                                    const DimensionedField<scalar, fvMesh>&,
                                                    const dictionary&);
    pSprayConservativeTemperatureFvPatchScalarField(
        const pSprayConservativeTemperatureFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, fvMesh>&,
        const fieldMapper&);
    pSprayConservativeTemperatureFvPatchScalarField(
        const pSprayConservativeTemperatureFvPatchScalarField&,
        const DimensionedField<scalar, fvMesh>&);

    virtual tmp<fvPatchScalarField> clone(const DimensionedField<scalar, fvMesh>& iF) const override
    {
        return tmp<fvPatchScalarField>(
            new pSprayConservativeTemperatureFvPatchScalarField(*this, iF));
    }
};
}
#endif
