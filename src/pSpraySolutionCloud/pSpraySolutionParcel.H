// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file pSpraySolutionParcel.H
 * @brief 具有独立LiBr质量分数的原生sprayParcel扩展。
 * 单个parcel代表nParticle个物理液滴；每滴质量与parcel总质量应分别处理。
 * 非挥发LiBr随液滴追踪、克隆、MPI传送和重启保存，不注册为气相组分。
 */
// Uses native sprayParcel tracking and momentum/thermal integration.
// Nonvolatile LiBr is an explicit conserved parcel property, never a gas species.
// In-flight phase change, breakup and atomisation are disabled in this first extension.
#ifndef pSpraySolutionParcel_H
#define pSpraySolutionParcel_H
#include "sprayParcel.H"
#include "PSpraySolutionThermo.H"

namespace Foam
{
class pSpraySolutionParcel : public sprayParcel
{
    scalar XLiBr_ = 0;

public:
    TypeName("pSpraySolutionParcel");
    using sprayParcel::sprayParcel;

    pSpraySolutionParcel(const pSpraySolutionParcel& p)
        : sprayParcel(p)
        , XLiBr_(p.XLiBr_)
    {
    }

    explicit pSpraySolutionParcel(const sprayParcel& p)
        : sprayParcel(p)
        , XLiBr_(0)
    {
    }

    pSpraySolutionParcel(Istream& is, bool read = true)
        : sprayParcel(is, read)
    {
        if (read)
        {
            if (is.format() == IOstream::ASCII)
                is >> XLiBr_;
            else
                is.read(reinterpret_cast<char*>(&XLiBr_), sizeof(scalar));
        }
        is.check("pSpraySolutionParcel(Istream)");
    }

    virtual autoPtr<particle> clone() const
    {
        return autoPtr<particle>(new pSpraySolutionParcel(*this));
    }

    static autoPtr<pSpraySolutionParcel> New(Istream& is)
    {
        return autoPtr<pSpraySolutionParcel>(new pSpraySolutionParcel(is));
    }

    scalar XLiBr() const
    {
        return XLiBr_;
    }

    scalar& XLiBr()
    {
        return XLiBr_;
    }

    /**
     * @brief 先刷新溶液物性，再使用原生动量/热积分并保持液滴质量。
     * rho变化时按相同每滴质量重算直径；原生定比热温升随后换算到共同焓参考。
     * 破碎、雾化、飞行相变和辐射在当前扩展中必须关闭。
     */
    template <class TrackCloudType>
    void calc(TrackCloudType& cloud, trackingData& td, const scalar dt)
    {
        if (cloud.breakup().type() != "none" || cloud.atomisation().type() != "none" ||
            cloud.phaseChange().type() != "none" || cloud.radiation())
            FatalErrorInFunction
                << "This prototype requires no breakup, atomisation, in-flight phase change or radiation"
                << exit(FatalError);
        const auto& table = pSprayTable(cloud.constProps().dict());
        const scalar Tstart = T(), pres = td.pc(), m = mass();
        auto a = table.evaluate(Tstart, XLiBr_);
        rho() = a.rho;
        Cp() = a.cp;
        mu() = a.mu;
        d() = cbrt(6 * m / (constant::mathematical::pi * rho()));
        const scalar h0 = pSprayEnthalpy(table, pres, Tstart, XLiBr_), cp0 = Cp();
        cloud.constProps().setTMax(353.15);
        // Explicitly select native no-phase-change thermal/momentum integration.
        ThermoParcel<MomentumParcel<particle>>::calc(cloud, td, dt);
        // 将原生积分得到的显热增量转为目标焓，再用真实组分相关物性反求T。
        const scalar target = h0 + cp0 * (T() - Tstart);
        scalar tt = T();
        for (label i = 0; i < 15; ++i)
        {
            const scalar res = pSprayEnthalpy(table, pres, tt, XLiBr_) - target;
            if (mag(res) < 1e-7)
                break;
            const scalar lo = max(scalar(278.15), tt - 1e-4), hi = min(scalar(353.15), tt + 1e-4);
            const scalar cp = (pSprayEnthalpy(table, pres, hi, XLiBr_) -
                               pSprayEnthalpy(table, pres, lo, XLiBr_)) /
                              (hi - lo);
            tt -= res / cp;
            if (tt < 278.15 || tt > 353.15)
                FatalErrorInFunction << "Parcel energy outside table" << exit(FatalError);
        }
        T() = tt;
        a = table.evaluate(T(), XLiBr_);
        rho() = a.rho;
        Cp() = a.cp;
        mu() = a.mu;
        d() = cbrt(6 * m / (constant::mathematical::pi * rho()));
    }

    /// 单独读取XLiBr粒子场；空云由size控制读取，避免要求不存在的粒子文件。
    template <class C> static void readSalt(C& c)
    {
        IOField<scalar> x(c.fieldIOobject("XLiBr", IOobject::MUST_READ), bool(c.size()));
        c.checkFieldIOobject(c, x);
        label i = 0;
        forAllIter (typename C, c, iter)
            iter().XLiBr_ = x[i++];
    }

    /// 原生spray字段之外写出XLiBr，维持粒子遍历顺序与原生字段一致。
    template <class C> static void writeSalt(const C& c)
    {
        IOField<scalar> x(c.fieldIOobject("XLiBr", IOobject::NO_READ), c.size());
        label i = 0;
        forAllConstIter (typename C, c, iter)
            x[i++] = iter().XLiBr_;
        x.write(c.size() > 0);
    }

    template <class C> static void readFields(C& c)
    {
        sprayParcel::readFields(c);
        readSalt(c);
    }

    template <class C, class Comp> static void readFields(C& c, const Comp& comp)
    {
        sprayParcel::readFields(c, comp);
        readSalt(c);
    }

    template <class C> static void writeFields(const C& c)
    {
        sprayParcel::writeFields(c);
        writeSalt(c);
    }

    template <class C, class Comp> static void writeFields(const C& c, const Comp& comp)
    {
        sprayParcel::writeFields(c, comp);
        writeSalt(c);
    }

    /// MPI/流序列化在原生parcel数据后追加XLiBr，与读取构造函数配对。
    friend Ostream& operator<<(Ostream& os, const pSpraySolutionParcel& p)
    {
        os << static_cast<const sprayParcel&>(p);
        if (os.format() == IOstream::ASCII)
            os << token::SPACE << p.XLiBr_;
        else
            os.write(reinterpret_cast<const char*>(&p.XLiBr_), sizeof(scalar));
        return os;
    }
};

template <> inline bool contiguous<pSpraySolutionParcel>()
{
    return false;
}
}
#endif
