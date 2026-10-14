// GPL-3.0-or-later. Extends OpenFOAM Foundation 14 ConeInjection.
// Point-source flat fan: uniform mass per angle is an explicit modelling assumption.
/**
 * @file PSprayFlatFanInjection.H
 * @brief 点源平面扇形注入；按扇面内角度均匀抽样属于当前模型假设。
 * fanDirection给出中心轴，fanTangent给出扇面方向，fanSpeed单位m/s。
 * fanHalfAngle输入为角度，内部转换为弧度；45°对应总夹角90°。
 * 流量、粒径分布等仍由原生ConeInjection及算例输入负责。
 */
#ifndef PSprayFlatFanInjection_H
#define PSprayFlatFanInjection_H
#include "ConeInjection.H"
#include "mathematicalConstants.H"

namespace Foam
{
template <class CloudType> class PSprayFlatFanInjection : public ConeInjection<CloudType>
{
    vector axis_, tangent_;
    scalar halfAngle_, speed_;

public:
    TypeName("pSprayFlatFanInjection");

    PSprayFlatFanInjection(const dictionary& dict, CloudType& owner, const word& name)
        : ConeInjection<CloudType>(dict, owner, name)
        , axis_(normalised(dict.lookup<vector>("fanDirection")))
        , tangent_(dict.lookup<vector>("fanTangent"))
        , halfAngle_(dict.lookup<scalar>("fanHalfAngle") * constant::mathematical::pi / 180)
        , speed_(dict.lookup<scalar>("fanSpeed"))
    {
        // Gram-Schmidt正交化使扇面方向垂直于中心轴，然后归一化。
        tangent_ -= axis_ * (axis_ & tangent_);
        if (mag(tangent_) < small || speed_ <= 0 || halfAngle_ <= 0 ||
            halfAngle_ > constant::mathematical::pi / 2)
            FatalIOErrorInFunction(dict)
                << "Invalid flat-fan axes, speed or angle" << exit(FatalIOError);
        tangent_ = normalised(tangent_);
    }

    PSprayFlatFanInjection(const PSprayFlatFanInjection& other)
        : ConeInjection<CloudType>(other)
        , axis_(other.axis_)
        , tangent_(other.tangent_)
        , halfAngle_(other.halfAngle_)
        , speed_(other.speed_)
    {
    }

    virtual autoPtr<InjectionModel<CloudType>> clone() const
    {
        return autoPtr<InjectionModel<CloudType>>(new PSprayFlatFanInjection(*this));
    }

    virtual void setProperties(const label parcelI,
                               const label nParcels,
                               const scalar time,
                               typename CloudType::parcelType::trackingData& td,
                               typename CloudType::parcelType& parcel)
    {
        ConeInjection<CloudType>::setProperties(parcelI, nParcels, time, td, parcel);
        // 在[-halfAngle,+halfAngle]均匀抽样；此处不采用空间圆锥立体角分布。
        const scalar angle = halfAngle_ * (2 * this->owner().rndGen().scalar01() - 1);
        parcel.U() = speed_ * (cos(angle) * axis_ + sin(angle) * tangent_);
    }
};
}
#endif
