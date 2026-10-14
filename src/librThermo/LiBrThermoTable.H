// Local model kernel. SI units; no OpenFOAM or Python dependency at runtime.
/**
 * @file LiBrThermoTable.H
 * @brief 水/LiBr溶液的共享查表内核，供液滴、壁膜和相变源项使用。
 *
 * 温度用K，X为LiBr质量分数[kg/kg]，压力为绝对压力[Pa]。
 * 输入表规定有效范围；越界抛出异常，不对高浓度或低温状态外推。
 * 与OpenFOAM相关的异常处理由调用层负责，本文件也可独立用于物性核验。
 */
#ifndef LiBrThermoTable_H
#define LiBrThermoTable_H
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace librSpray
{
struct State
{
    // 字段顺序与表一致：rho[kg/m³]、cp[J/(kg·K)]、mu[Pa·s]、pEq[Pa]、
    // 水相变焓差latentWater[J/kg]、溶液比焓h[J/kg]、dh_dX[J/kg]。
    double rho, cp, mu, pEq, latentWater, h, dh_dX;
};

class Table
{
    // T-X规则网格；每行的黏度和平衡蒸气压以自然对数存储。
    int nT_, nX_;
    double Tmin_, dT_, Xmin_, dX_;
    std::vector<std::array<double, 7>> rows_;

public:
    /// 读取LIBR_WATER_THERMO_V2格式，检查网格、有限数值和尾部多余数据。
    explicit Table(const std::string& path)
    {
        std::ifstream file(path);
        std::string version;
        if (!(file >> version) || version != "LIBR_WATER_THERMO_V2")
            throw std::runtime_error("Invalid or missing LiBr water table: " + path);
        if (!(file >> nT_ >> nX_ >> Tmin_ >> dT_ >> Xmin_ >> dX_) || nT_ < 2 || nX_ < 2 ||
            nT_ > 10000 || nX_ > 10000 || dT_ <= 0 || dX_ <= 0)
            throw std::runtime_error("Invalid LiBr table grid");
        rows_.resize(std::size_t(nT_) * nX_);
        for (auto& row : rows_)
            for (auto& v : row)
                if (!(file >> v) || !std::isfinite(v))
                    throw std::runtime_error("Invalid LiBr table value");
        double extra;
        if (file >> extra)
            throw std::runtime_error("Unexpected extra LiBr table data");
    }

    /**
     * @brief 计算指定T、X下的物性，并保持h与其温度/浓度导数一致。
     * 密度等字段先双线性插值；焓和导数再由同一个Hermite曲面计算。
     * 黏度在稀溶液区采用解析浓度项，以解析X→0附近的陡峭变化。
     */
    State evaluate(double T, double X) const
    {
        if (!std::isfinite(T) || !std::isfinite(X) || T < Tmin_ || T > Tmin_ + (nT_ - 1) * dT_ ||
            X < Xmin_ || X > Xmin_ + (nX_ - 1) * dX_)
            throw std::out_of_range("LiBr state outside table; extrapolation is disabled");
        const double u = (T - Tmin_) / dT_, v = (X - Xmin_) / dX_;
        const int i = std::min(int(u), nT_ - 2), j = std::min(int(v), nX_ - 2);
        const double a = u - i, b = v - j;
        std::array<double, 7> out{};
        for (int k = 0; k < 7; ++k)
            out[k] =
                (1 - a) * ((1 - b) * rows_[i * nX_ + j][k] + b * rows_[i * nX_ + j + 1][k]) +
                a * ((1 - b) * rows_[(i + 1) * nX_ + j][k] + b * rows_[(i + 1) * nX_ + j + 1][k]);
        // Tensor-product Hermite interpolation of the SAME solution enthalpy.
        // Nodes provide h, dh/dT and dh/dX; the cross derivative is computed
        // from the heat-capacity table with centred/one-sided finite differences.
        const auto hermite = [](double q, double y0, double y1, double d0, double d1, double width)
        {
            const double q2 = q * q, q3 = q2 * q;
            return (2 * q3 - 3 * q2 + 1) * y0 + (q3 - 2 * q2 + q) * width * d0 +
                   (-2 * q3 + 3 * q2) * y1 + (q3 - q2) * width * d1;
        };
        const auto hermiteD = [](double q, double y0, double y1, double d0, double d1, double width)
        {
            const double q2 = q * q;
            return ((6 * q2 - 6 * q) * y0 + (-6 * q2 + 6 * q) * y1) / width +
                   (3 * q2 - 4 * q + 1) * d0 + (3 * q2 - 2 * q) * d1;
        };
        // 混合导数∂²h/(∂T∂X)由cp表沿X差分获得；边界采用单侧公式。
        const auto crossDerivative = [&](int ti, int xi)
        {
            const auto cp = [&](int k) { return rows_[ti * nX_ + k][1]; };
            if (nX_ == 2)
                return (cp(1) - cp(0)) / dX_;
            if (xi == 0)
                return (-3 * cp(0) + 4 * cp(1) - cp(2)) / (2 * dX_);
            if (xi == nX_ - 1)
                return (3 * cp(xi) - 4 * cp(xi - 1) + cp(xi - 2)) / (2 * dX_);
            if (xi >= 2 && xi < nX_ - 2)
                return (cp(xi - 2) - 8 * cp(xi - 1) + 8 * cp(xi + 1) - cp(xi + 2)) / (12 * dX_);
            return (cp(xi + 1) - cp(xi - 1)) / (2 * dX_);
        };
        const auto interpolateRow = [&](int ti)
        {
            const auto& s0 = rows_[ti * nX_ + j];
            const auto& s1 = rows_[ti * nX_ + j + 1];
            const double mixed0 = crossDerivative(ti, j), mixed1 = crossDerivative(ti, j + 1);
            return std::array<double, 4>{{hermite(b, s0[5], s1[5], s0[6], s1[6], dX_),
                                          hermite(b, s0[1], s1[1], mixed0, mixed1, dX_),
                                          hermiteD(b, s0[5], s1[5], s0[6], s1[6], dX_),
                                          hermiteD(b, s0[1], s1[1], mixed0, mixed1, dX_)}};
        };
        const auto low = interpolateRow(i), high = interpolateRow(i + 1);
        out[5] = hermite(a, low[0], high[0], low[1], high[1], dT_);
        out[1] = hermiteD(a, low[0], high[0], low[1], high[1], dT_);
        out[6] = hermite(a, low[2], high[2], low[3], high[3], dT_);
        // Evaluate concentration dependence analytically to resolve the steep
        // dilute-limit viscosity curve; interpolate only log(mu_water) in T.
        const double logMuWater = (1 - a) * rows_[i * nX_][2] + a * rows_[(i + 1) * nX_][2];
        const double mu = std::exp(logMuWater + 2.477 * std::pow(T / 647.096, 3.609 - 5.871 * X) *
                                                    std::pow(X, 0.2365));
        return {out[0], out[1], mu, std::exp(out[3]), out[4], out[5], out[6]};
    }

    /// 返回{ln(pEq), ∂ln(pEq)/∂T, ∂ln(pEq)/∂X}；不是压力本身。
    std::array<double, 3> pressureData(double T, double X) const
    {
        if (T < Tmin_ || T > Tmin_ + (nT_ - 1) * dT_ || X < Xmin_ || X > Xmin_ + (nX_ - 1) * dX_)
            throw std::out_of_range("Pressure state outside table");
        const double u = (T - Tmin_) / dT_, v = (X - Xmin_) / dX_;
        const int ti = std::min(int(u), nT_ - 2), xi = std::min(int(v), nX_ - 2);
        const double a = u - ti, b = v - xi;
        const double q00 = rows_[ti * nX_ + xi][3], q01 = rows_[ti * nX_ + xi + 1][3];
        const double q10 = rows_[(ti + 1) * nX_ + xi][3], q11 = rows_[(ti + 1) * nX_ + xi + 1][3];
        return {{(1 - a) * ((1 - b) * q00 + b * q01) + a * ((1 - b) * q10 + b * q11),
                 ((1 - b) * (q10 - q00) + b * (q11 - q01)) / dT_,
                 ((1 - a) * (q01 - q00) + a * (q11 - q10)) / dX_}};
    }

    /// 只求密度，用于盐密度反求浓度，避免重复计算整组物性。
    double density(double T, double X) const
    {
        if (!std::isfinite(T) || !std::isfinite(X) || T < Tmin_ || T > Tmin_ + (nT_ - 1) * dT_ ||
            X < Xmin_ || X > Xmin_ + (nX_ - 1) * dX_)
            throw std::out_of_range("Density state outside table");
        const double u = (T - Tmin_) / dT_, v = (X - Xmin_) / dX_;
        const int ti = std::min(int(u), nT_ - 2), xi = std::min(int(v), nX_ - 2);
        const double a = u - ti, b = v - xi;
        return (1 - a) * ((1 - b) * rows_[ti * nX_ + xi][0] + b * rows_[ti * nX_ + xi + 1][0]) +
               a * ((1 - b) * rows_[(ti + 1) * nX_ + xi][0] +
                    b * rows_[(ti + 1) * nX_ + xi + 1][0]);
    }

    /**
     * @brief 由单位液体体积中的盐质量[kg/m³]反求X。
     * 解X*rho(T,X)=saltPerLiquidVolume，使更新密度后仍保持盐质量。
     * 在浓度网格内二分定位，再利用分段线性rho解析求根；超过表上限报错。
     */
    double massFractionFromSaltDensity(double T, double saltPerLiquidVolume) const
    {
        const double xmax = Xmin_ + (nX_ - 1) * dX_;
        if (!std::isfinite(saltPerLiquidVolume) || saltPerLiquidVolume < 0)
            throw std::out_of_range("Negative or nonfinite salt concentration");
        if (saltPerLiquidVolume == 0)
            return 0;
        const double maxSalt = xmax * density(T, xmax);
        if (saltPerLiquidVolume > maxSalt)
            throw std::out_of_range("Conserved salt would require XLiBr above the table range");
        int lo = 0, hi = nX_ - 1;
        while (hi - lo > 1)
        {
            const int mid = (lo + hi) / 2;
            const double x = Xmin_ + mid * dX_;
            if (x * density(T, x) < saltPerLiquidVolume)
                lo = mid;
            else
                hi = mid;
        }
        const double x0 = Xmin_ + lo * dX_, x1 = Xmin_ + hi * dX_;
        const double r0 = density(T, x0), r1 = density(T, x1);
        const double slope = (r1 - r0) / (x1 - x0), intercept = r0 - slope * x0;
        // 二次方程采用避免两个接近数相减的根表达式，减少消减误差。
        return std::abs(slope) < 1e-12
                   ? saltPerLiquidVolume / r0
                   : 2 * saltPerLiquidVolume /
                         (intercept +
                          std::sqrt(intercept * intercept + 4 * slope * saltPerLiquidVolume));
    }
};

/// 等温Schrage辅助通量[kg/(m²·s)]：正值为液体向蒸气蒸发。
/// 双温度局部反应闭合见PSprayFilmReaction.H，此函数不代替该预测器。
inline double schrageFlux(const State& state, double T, double pV, double accommodation)
{
    if (!std::isfinite(T) || T <= 0 || !std::isfinite(pV) || pV <= 0 ||
        !std::isfinite(accommodation) || accommodation <= 0 || accommodation > 1)
        throw std::invalid_argument("Schrage closure requires T>0, p_abs>0, 0<accommodation<=1");
    constexpr double Rv = 8.31446261815324 / 0.018015268;
    constexpr double pi = 3.14159265358979323846;
    return 2 * accommodation / (2 - accommodation) * (state.pEq - pV) / std::sqrt(2 * pi * Rv * T);
}

struct Sources
{
    // 两相质量源[kg/(m³·s)]、盐源[kg/(m³·s)]、显热源[W/m³]。
    double liquidMass, vapourMass, saltMass, sensibleHeat;
};

/// 相变只转移水：液体和蒸气质量源相反，LiBr盐源严格为零。
inline Sources phaseSources(double mDotPerVolume, double differentialLatent)
{
    // mDot > 0 transfers WATER from solution/water to vapour; LiBr is nonvolatile.
    // Sensible heat is negative on evaporation, positive on absorption.
    return {-mDotPerVolume, mDotPerVolume, 0.0, -mDotPerVolume * differentialLatent};
}
}
#endif
