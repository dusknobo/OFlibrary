// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file LiBrSolutionMixture.H
 * @brief 将单元/边界的LiBr质量分数传给OpenFOAM热物性和输运接口。
 * numericalTail仅标记低于壁膜厚度分辨率的数值尾迹，供导热参考处理使用。
 */
#ifndef LiBrSolutionMixture_H
#define LiBrSolutionMixture_H
#include "LiBrSolutionProperties.H"

namespace Foam
{
struct LiBrSolutionComposition
{
    // X[kg/kg]；尾迹标记不代表另一种液体组分。
    scalar X;
    bool numericalTail;
};

class LiBrSolutionMixture
{
    // 按查询位置更新并复用物性对象；返回引用依赖当前组分查询状态。
    mutable LiBrSolutionProperties mixture_;

public:
    typedef LiBrSolutionProperties thermoType;
    typedef LiBrSolutionProperties thermoMixtureType;
    typedef LiBrSolutionProperties transportMixtureType;

    static word typeName()
    {
        return "librSolutionMixture";
    }

    explicit LiBrSolutionMixture(const dictionary& dict)
        : mixture_(dict.subDict("mixture"))
    {
    }

    /// 每次求物性前设置该位置的X和尾迹标记，避免误用纯水nil组分。
    const thermoMixtureType& thermoMixture(const LiBrSolutionComposition& c) const
    {
        mixture_.composition(c.X, c.numericalTail);
        return mixture_;
    }

    const transportMixtureType& transportMixture(const LiBrSolutionComposition& c) const
    {
        return thermoMixture(c);
    }

    const transportMixtureType& transportMixture(const LiBrSolutionComposition& c,
                                                 const thermoMixtureType&) const
    {
        return thermoMixture(c);
    }

    /// 重新读取mixture字典，使查表与导热模式随物性配置更新。
    void read(const dictionary& dict)
    {
        mixture_ = LiBrSolutionProperties(dict.subDict("mixture"));
    }
};
}
#endif
