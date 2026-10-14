// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file pSprayRebasedTemperatureFvPatchScalarField.H
 * @brief 保持真实混合温度和梯度算子不变的参考值重参数化。
 * refValue是mixed线性表达的参考温度，不必等于真实界面温度。
 * referenceRebase=false仅诊断；Min/Max[K]不是物理场的截断器。
 */
#ifndef pSprayRebasedTemperatureFvPatchScalarField_H
#define pSprayRebasedTemperatureFvPatchScalarField_H
#include "pSprayConservativeTemperatureFvPatchScalarField.H"
#include "Switch.H"

namespace Foam
{
class pSprayRebasedTemperatureFvPatchScalarField
    : public pSprayConservativeTemperatureFvPatchScalarField
{
    Switch referenceRebase_;
    scalar referenceMin_, referenceMax_;

public:
    TypeName("pSprayRebasedTemperature");
    pSprayRebasedTemperatureFvPatchScalarField(const fvPatch&,
                                               const DimensionedField<scalar, fvMesh>&,
                                               const dictionary&);
    pSprayRebasedTemperatureFvPatchScalarField(const pSprayRebasedTemperatureFvPatchScalarField&,
                                               const fvPatch&,
                                               const DimensionedField<scalar, fvMesh>&,
                                               const fieldMapper&);
    pSprayRebasedTemperatureFvPatchScalarField(const pSprayRebasedTemperatureFvPatchScalarField&,
                                               const DimensionedField<scalar, fvMesh>&);

    virtual tmp<fvPatchScalarField>
    clone(const DimensionedField<scalar, fvMesh>& field) const override
    {
        return tmp<fvPatchScalarField>(
            new pSprayRebasedTemperatureFvPatchScalarField(*this, field));
    }

    /// 先更新面积守恒热耦合，再重表示refValue/refGrad并核对算子不变性。
    virtual void updateCoeffs() override;
    /// 写出开关和参考范围，供重启恢复。
    virtual void write(Ostream&) const override;
};
}
#endif
