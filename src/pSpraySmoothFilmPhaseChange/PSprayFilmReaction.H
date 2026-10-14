// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file PSprayFilmReaction.H
 * @brief 一个时间步内的局部气膜相变反应预测器。
 *
 * 只预测水的相间转移与对应能量/浓度变化，不求空间输运或管壁热惯性。
 * 返回的广延转移量在同一时间步冻结一次，供气相和液膜有限体积方程配对使用。
 * 各相变库保留相同副本；修改物理闭合时必须同步并重新核验时间分裂误差。
 */
#ifndef PSprayFilmReaction_H
#define PSprayFilmReaction_H
#include "LiBrThermoTable.H"

namespace librSpray
{
// SI. This local implicit reaction predictor supplies ONE frozen extensive
// transfer to the native gas and film FV equations over the same time step.
// Advection, wall heat and conduction are solved by OpenFOAM, outside this
// predictor. Temporal splitting error must be checked by time refinement.
struct FilmReactionInput
{
    // ml、mg[kg]，vg[m³]，area[m²]，Tl/Tg[K]，X[kg/kg]，dt[s]。
    // beta为Schrage适应系数，maxFraction为单步允许转移的可用质量比例。
    double ml, mg, vg, area, Tl, Tg, X, dt, beta, maxFraction;
    // 蒸气焓参考：hv=href+cpv*(T-Tref)，cpv[J/(kg·K)]、href[J/kg]。
    double cpv, href, Tref;
};

struct FilmReactionResult
{
    // dm[kg]正为蒸发、负为吸收；hv[J/kg]取转移水的供体温度。
    // Tl/Tg/X是局部预测终态，不直接替代原生方程求得的输运场。
    // residual[kg]为dm与双温度Schrage隐式转移量的差。
    double dm = 0, hv = 0, Tl = 0, Tg = 0, X = 0, residual = 0;
    int iterations = 0;
    bool limited = false;
};

/**
 * @brief 同时约束盐质量、液体能量与气相能量后，求解隐式转移质量dm。
 * 先为候选dm求Tl/Tg/X，再对Schrage质量残差使用带括区的割线/二分。
 * limited表示最优根超出可用水、蒸气或有效物性状态边界，应单独计入限幅审计。
 */
inline FilmReactionResult filmReaction(const Table& table, const FilmReactionInput& a)
{
    constexpr double Rv = 8.31446261815324 / 0.018015268;
    const double cv = a.cpv - Rv, c = a.href - a.cpv * a.Tref;
    const auto e = [&](double T, double X)
    {
        const auto s = table.evaluate(T, X);
        return s.h - s.pEq / s.rho;
    };
    FilmReactionResult out;
    out.Tl = a.Tl;
    out.Tg = a.Tg;
    out.X = a.X;
    if (a.ml <= 1e-20 || a.mg <= 1e-25 || a.area <= 0)
        return out;
    // LiBr不挥发：反应前后盐质量salt保持不变；只改变水质量。
    const double el = e(a.Tl, a.X), salt = a.ml * a.X;
    const double A =
        a.area * 2 * a.beta / (2 - a.beta) / std::sqrt(2 * 3.14159265358979323846 * Rv);
    // 候选状态必须保留两相正质量，并满足浓度与导热关联式的温度范围。
    const auto trial = [&](double dm, FilmReactionResult& z) -> bool
    {
        z.dm = dm;
        z.X = salt / (a.ml - dm);
        if (!(a.ml - dm > 0 && a.mg + dm > 0 && z.X >= 0 && z.X <= .65))
            return false;
        const double tmin = z.X >= .3 ? 293.15 : 278.15, tmax = 353.15;
        double tg = a.Tg;
        if (dm < 0)
            tg = a.mg * cv * a.Tg / (a.mg * cv - dm * Rv);
        // 液体能量残差[J]：新液体内能 + 转移水焓 - 原液体内能。
        // 蒸发供体为液体Tl，吸收供体为气体Tg；同一hv供两相配对使用。
        const auto f = [&](double tl)
        {
            const double hv = c + a.cpv * (dm >= 0 ? tl : tg);
            return (a.ml - dm) * e(tl, z.X) + dm * hv - a.ml * el;
        };
        double lo = tmin, hi = tmax, flo = f(lo), fhi = f(hi);
        const double etol =
            2e-13 * a.ml * std::max(std::abs(el), table.evaluate(a.Tl, a.X).cp * a.Tl);
        if (flo > etol || fhi < -etol)
            return false;
        double tl = std::max(lo, std::min(hi, a.Tl));
        for (int k = 0; k < 24; ++k)
        {
            const double ft = f(tl);
            if (std::abs(ft) < etol)
                break;
            if (ft > 0)
                hi = tl;
            else
                lo = tl;
            const double cp = table.evaluate(tl, z.X).cp;
            double tn = tl - ft / ((a.ml - dm) * cp + (dm >= 0 ? dm * a.cpv : 0));
            if (!(tn > lo && tn < hi))
                tn = .5 * (lo + hi);
            tl = tn;
        }
        if (dm >= 0)
            tg = (a.mg * cv * a.Tg + dm * a.cpv * tl) / ((a.mg + dm) * cv);
        if (!(tg > 150 && tg < 1000))
            return false;
        z.Tl = tl;
        z.Tg = tg;
        z.hv = c + a.cpv * (dm >= 0 ? tl : tg);
        // 局部纯水蒸气压力由预测质量和温度求得，为绝对压力。
        // 双温度Schrage项分别使用pEq(Tl,X)/sqrt(Tl)与pg/sqrt(Tg)。
        const double pg = (a.mg + dm) * Rv * tg / a.vg;
        z.residual =
            dm - a.dt * A * (table.evaluate(tl, z.X).pEq / std::sqrt(tl) - pg / std::sqrt(tg));
        return true;
    };
    // 吸收不能消耗过多可用蒸气；蒸发不能超过可用水或使X高于0.65。
    double lo = -a.maxFraction * a.mg;
    double hi = a.maxFraction * a.ml * (1 - a.X);
    if (salt > 0)
        hi = std::min(hi, std::max(0., a.ml - salt / .65) * (1 - 1e-10));
    FilmReactionResult zl, zh;
    // 如端点不满足能量/物性范围，向零转移量二分收缩到有效端点。
    const auto bound = [&](double v, FilmReactionResult& z)
    {
        if (trial(v, z))
            return v;
        double good = 0, bad = v;
        for (int k = 0; k < 48; ++k)
        {
            const double mid = .5 * (good + bad);
            FilmReactionResult q;
            if (trial(mid, q))
                good = mid;
            else
                bad = mid;
        }
        trial(good, z);
        return good;
    };
    lo = bound(lo, zl);
    hi = bound(hi, zh);
    // 根在有效括区之外时返回端点，并明确标记限幅，供累计收支核验。
    if (zl.residual >= 0)
    {
        zl.limited = true;
        return zl;
    }
    if (zh.residual <= 0)
    {
        zh.limited = true;
        return zh;
    }
    // Bracketed secant, periodically bisected to avoid false-position stalling.
    for (int k = 0; k < 96; ++k)
    {
        double dm = (lo * zh.residual - hi * zl.residual) / (zh.residual - zl.residual);
        if (k % 4 == 3 || !(dm > lo && dm < hi))
            dm = .5 * (lo + hi);
        if (!trial(dm, out))
            throw std::runtime_error("Invalid interior film reaction state");
        out.iterations = k + 1;
        if (std::abs(out.residual) < 1e-12 * std::min(a.ml, a.mg) ||
            hi - lo < 1e-12 * std::min(a.ml, a.mg))
            return out;
        if (out.residual > 0)
        {
            hi = dm;
            zh = out;
        }
        else
        {
            lo = dm;
            zl = out;
        }
    }
    throw std::runtime_error("Film reaction predictor did not converge");
}
}
#endif
