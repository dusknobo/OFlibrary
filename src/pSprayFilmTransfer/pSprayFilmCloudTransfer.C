/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2026 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

/**
 * @file pSprayFilmCloudTransfer.C
 * @brief 膜侧捕集、回滴与守恒源项组装。
 * 显式累计广延量，使用同一实际deltaT换算回滴和撞膜方程源。
 * correct()保持无操作；云cacheFilmFields调用resetFromCloudFields完成准备。
 */
#include "pSprayFilmCloudTransfer.H"
#include "mappedFvPatchBaseBase.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
defineTypeNameAndDebug(pSprayFilmCloudTransfer, 0);

addToRunTimeSelectionTable(fvModel, pSprayFilmCloudTransfer, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::pSprayFilmCloudTransfer::pSprayFilmCloudTransfer(const word& sourceName,
                                                           const word& modelType,
                                                           const fvMesh& mesh,
                                                           const dictionary& dict)
    : fvModel(sourceName, modelType, mesh, dict)
    , film_(mesh.lookupObject<solvers::isothermalFilm>(solver::typeName))
    , cloudFieldsTransferred_(false)
    , ejectionCloud_(dict.lookupOrDefault<word>("ejectionCloud", "solutionSpray"))
    , maxEjectionFraction_(dict.lookupOrDefault<scalar>("maxEjectionFraction", .25))
    , ejectionMass_(
          IOobject(
              "ejectionMass", mesh.time().name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
          mesh,
          dimensionedScalar(dimMass, 0))
    , ejectionSalt_(
          IOobject(
              "ejectionSalt", mesh.time().name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
          mesh,
          dimensionedScalar(dimMass, 0))
    , ejectionEnergy_(
          IOobject(
              "ejectionEnergy", mesh.time().name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
          mesh,
          dimensionedScalar(dimEnergy, 0))
    , ejectionMomentum_(IOobject("ejectionMomentum",
                                 mesh.time().name(),
                                 mesh,
                                 IOobject::NO_READ,
                                 IOobject::AUTO_WRITE),
                        mesh,
                        dimensionedVector(dimMass * dimVelocity, Zero))
    , correctEjection_(false)
    , ejection_(dict.found("ejection") ? ejectionModel::New(dict.subDict("ejection"), film_)
                                       : autoPtr<ejectionModel>(nullptr))
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::pSprayFilmCloudTransfer::addSupFields() const
{
    return wordList{film_.alpha.name(), film_.thermo.he().name(), film_.U.name(), "XLiBr"};
}

void Foam::fv::pSprayFilmCloudTransfer::correct()
{
    // Ejection is prepared when the return cloud caches film properties.
    // Both the created parcels and the film sinks use the same frozen
    // extensive quantities in THIS time step, including variable deltaT.
}

void Foam::fv::pSprayFilmCloudTransfer::prepareEjection()
{
    // 本步全部回滴广延量先置零，随后与新parcel生成量一一对应。
    ejectionMass_ = dimensionedScalar(dimMass, 0);
    ejectionSalt_ = dimensionedScalar(dimMass, 0);
    ejectionEnergy_ = dimensionedScalar(dimEnergy, 0);
    ejectionMomentum_ = dimensionedVector(dimMass * dimVelocity, Zero);
    if (!ejection_.valid())
        return;
    if (!(maxEjectionFraction_ > 0 && maxEjectionFraction_ <= .5))
        FatalErrorInFunction << "Require 0<maxEjectionFraction<=0.5" << exit(FatalError);
    ejection_->correct();
    const auto& X = mesh().lookupObject<volScalarField>("XLiBr");
    const auto& he = film_.thermo.he();
    const scalar dt = mesh().time().deltaTValue();
    const volScalarField* phase = mesh().foundObject<volScalarField>("phaseMassRate")
                                      ? &mesh().lookupObject<volScalarField>("phaseMassRate")
                                      : nullptr;
    scalar total = 0, limited = 0;
    label parcels = 0;
    forAll (ejectionMass_, i)
    {
        const scalar m = film_.alpha[i] * film_.rho[i] * mesh().V()[i];
        // 预留已请求的蒸发水量，防止相变与回滴同时消耗同一液体库存。
        const scalar reserved = phase ? max((*phase)[i], scalar(0)) * dt * mesh().V()[i] : 0;
        const scalar request = max(m * ejection_->rate()[i] * dt, scalar(0));
        scalar dm = min(request, maxEjectionFraction_ * max(m - reserved, scalar(0)));
        const scalar d = ejection_->diameter()[i];
        const scalar dropMass = constant::mathematical::pi / 6 * film_.rho[i] * pow3(d);
        // Native SurfaceFilmModel drops parcels of nParticle <= 0.001.
        // Retain that inventory in the film instead of losing it.
        if (d <= 0 || dm <= .001 * dropMass)
            dm = 0;
        // 回滴继承膜的X、e和U；新云质量与膜负源使用同一个dm。
        ejectionMass_[i] = dm;
        ejectionSalt_[i] = dm * X[i];
        ejectionEnergy_[i] = dm * he[i];
        ejectionMomentum_[i] = dm * film_.U[i];
        total += dm;
        if (dm > 0)
            ++parcels;
        if (dm < request)
            limited += request - dm;
    }
    Info << "Same-step ejection " << mesh().name()
         << " mass kg=" << returnReduce(total, sumOp<scalar>())
         << " cells=" << returnReduce(parcels, sumOp<label>())
         << " retained kg=" << returnReduce(limited, sumOp<scalar>()) << nl;
}

template <class Type>
Foam::tmp<Foam::VolInternalField<Type>> inline Foam::fv::pSprayFilmCloudTransfer::
    CloudToFilmTransferRate(const Field<Type>& prop, const dimensionSet& dimProp) const
{
    tmp<VolInternalField<Type>> tSu(VolInternalField<Type>::New(
        "Su", mesh(), dimensioned<Type>(dimProp / dimVolume / dimTime, Zero)));

    if (cloudFieldsTransferred_)
    {
        const fvMesh& cloudMesh = film_.surfacePatchMap().nbrMesh();
        const label cloudPatchi = film_.surfacePatchMap().nbrFvPatch().index();

        // 先按气相面面积换算再映射，随后按原生膜几何VbyA和dt求体积源。
        UIndirectList<Type>(tSu.ref(), film_.surfacePatch().faceCells()) =
            film_.surfacePatchMap().fromNeighbour(prop / cloudMesh.boundary()[cloudPatchi].magSf());

        tSu.ref().primitiveFieldRef() /= film_.VbyA;
        tSu.ref().primitiveFieldRef() /= mesh().time().deltaTValue();
    }

    return tSu;
}

void Foam::fv::pSprayFilmCloudTransfer::addSup(const volScalarField& rho,
                                               const volScalarField& alpha,
                                               fvMatrix<scalar>& eqn) const
{
    if (debug)
    {
        Info << type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (&alpha == &film_.alpha && &eqn.psi() == &film_.alpha)
    {
        eqn += CloudToFilmTransferRate<scalar>(massFromCloud_, dimMass);

        if (ejection_.valid())
        {
            eqn -= ejectionMass_() / (mesh().V() * mesh().time().deltaT());
        }
    }
    else
    {
        FatalErrorInFunction << "Support for field " << alpha.name() << " is not implemented"
                             << exit(FatalError);
    }
}

void Foam::fv::pSprayFilmCloudTransfer::addSup(const volScalarField& alpha,
                                               const volScalarField& rho,
                                               const volScalarField& he,
                                               fvMatrix<scalar>& eqn) const
{
    if (debug)
    {
        Info << type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (&he == &film_.thermo.he() && &eqn.psi() == &film_.thermo.he())
    {
        eqn += CloudToFilmTransferRate<scalar>(energyFromCloud_, dimEnergy);

        if (ejection_.valid())
        {
            eqn -= ejectionEnergy_() / (mesh().V() * mesh().time().deltaT());
        }
    }
    // 盐输运使用撞膜/回滴的真实盐质量，而不是仅指定某个目标浓度。
    else if (he.name() == "XLiBr")
    {
        eqn += CloudToFilmTransferRate<scalar>(saltFromCloud_, dimMass);
        if (ejection_.valid())
            eqn -= ejectionSalt_() / (mesh().V() * mesh().time().deltaT());
    }
    else
    {
        FatalErrorInFunction << "Support for field " << he.name() << " is not implemented"
                             << exit(FatalError);
    }
}

void Foam::fv::pSprayFilmCloudTransfer::addSup(const volScalarField& alpha,
                                               const volScalarField& rho,
                                               const volVectorField& U,
                                               fvMatrix<vector>& eqn) const
{
    if (debug)
    {
        Info << type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (&U == &film_.U && &U == &film_.U)
    {
        eqn += CloudToFilmTransferRate<vector>(momentumFromCloud_, dimensions::momentum);

        if (ejection_.valid())
        {
            eqn -= ejectionMomentum_() / (mesh().V() * mesh().time().deltaT());
        }
    }
    else
    {
        FatalErrorInFunction << "Support for field " << U.name() << " is not implemented"
                             << exit(FatalError);
    }
}

void Foam::fv::pSprayFilmCloudTransfer::resetFromCloudFields()
{
    // 两种云共享一次清零，保持同一步全部撞膜贡献。
    if (lastCloudReset_ == mesh().time().timeIndex())
        return;
    lastCloudReset_ = mesh().time().timeIndex();
    const fvMesh& cloudMesh = film_.surfacePatchMap().nbrMesh();
    const label cloudPatchi = film_.surfacePatchMap().nbrFvPatch().index();
    const label nCloudPatchFaces = cloudMesh.boundary()[cloudPatchi].size();

    if (massFromCloud_.size() != nCloudPatchFaces)
    {
        massFromCloud_.setSize(nCloudPatchFaces);
        saltFromCloud_.setSize(nCloudPatchFaces);
        momentumFromCloud_.setSize(nCloudPatchFaces);
        energyFromCloud_.setSize(nCloudPatchFaces);
    }

    massFromCloud_ = 0;
    saltFromCloud_ = 0;
    momentumFromCloud_ = Zero;
    energyFromCloud_ = 0;

    cloudFieldsTransferred_ = true;

    // 保留原生状态标记；本实现立即准备同一步回滴，不等待空的correct()。
    correctEjection_ = true;
    prepareEjection();
}

void Foam::fv::pSprayFilmCloudTransfer::parcelFromCloud(const label facei,
                                                        const scalar mass,
                                                        const vector& momentum,
                                                        const scalar energy,
                                                        const scalar saltMass)
{
    // 此处只累计一个parcel事件的广延量，不再乘nParticle或除dt。
    massFromCloud_[facei] += mass;
    saltFromCloud_[facei] += saltMass;
    momentumFromCloud_[facei] += momentum;
    energyFromCloud_[facei] += energy;
}

template <class Type>
Foam::tmp<Foam::Field<Type>> inline Foam::fv::pSprayFilmCloudTransfer::filmToCloudTransfer(
    const VolInternalField<Type>& prop) const
{
    return film_.surfacePatchMap().toNeighbour(Field<Type>(prop, film_.surfacePatch().faceCells()));
}

bool Foam::fv::pSprayFilmCloudTransfer::ejecting(const word& cloudName) const
{
    return ejection_.valid() && cloudName == ejectionCloud_;
}

Foam::tmp<Foam::Field<Foam::scalar>> Foam::fv::pSprayFilmCloudTransfer::ejectedMassToCloud() const
{
    return filmToCloudTransfer<scalar>(ejectionMass_());
}

Foam::tmp<Foam::Field<Foam::scalar>>
Foam::fv::pSprayFilmCloudTransfer::ejectedDiameterToCloud() const
{
    return filmToCloudTransfer<scalar>(ejection_->diameter());
}

Foam::tmp<Foam::Field<Foam::scalar>> Foam::fv::pSprayFilmCloudTransfer::deltaToCloud() const
{
    return filmToCloudTransfer<scalar>(film_.delta);
}

Foam::tmp<Foam::Field<Foam::vector>> Foam::fv::pSprayFilmCloudTransfer::UToCloud() const
{
    return filmToCloudTransfer<vector>(film_.U);
}

Foam::tmp<Foam::Field<Foam::scalar>> Foam::fv::pSprayFilmCloudTransfer::rhoToCloud() const
{
    return filmToCloudTransfer<scalar>(film_.rho);
}

Foam::tmp<Foam::Field<Foam::scalar>> Foam::fv::pSprayFilmCloudTransfer::TToCloud() const
{
    return filmToCloudTransfer<scalar>(film_.thermo.T());
}

Foam::tmp<Foam::Field<Foam::scalar>> Foam::fv::pSprayFilmCloudTransfer::CpToCloud() const
{
    return filmToCloudTransfer<scalar>(film_.thermo.Cp());
}

void Foam::fv::pSprayFilmCloudTransfer::topoChange(const polyTopoChangeMap& map)
{
    // Set the cloud field state to false, will be updated by the cloud tracking
    // If the film is evaluated before the cloud it would be better
    // if the cloud fields were mapped
    cloudFieldsTransferred_ = false;

    if (ejection_.valid())
    {
        ejection_->topoChange(map);
    }
}

void Foam::fv::pSprayFilmCloudTransfer::mapMesh(const polyMeshMap& map)
{
    // Set the cloud field state to false, will be updated by the cloud tracking
    // If the film is evaluated before the cloud it would be better
    // if the cloud fields were mapped
    cloudFieldsTransferred_ = false;

    if (ejection_.valid())
    {
        ejection_->mapMesh(map);
    }
}

void Foam::fv::pSprayFilmCloudTransfer::distribute(const polyDistributionMap& map)
{
    // Set the cloud field state to false, will be updated by the cloud tracking
    // If the film is evaluated before the cloud it would be better
    // if the cloud fields were mapped
    cloudFieldsTransferred_ = false;

    if (ejection_.valid())
    {
        ejection_->distribute(map);
    }
}

bool Foam::fv::pSprayFilmCloudTransfer::movePoints()
{
    if (ejection_.valid())
    {
        ejection_->movePoints();
    }

    return true;
}

// ************************************************************************* //

Foam::tmp<Foam::Field<Foam::scalar>> Foam::fv::pSprayFilmCloudTransfer::XToCloud() const
{
    return filmToCloudTransfer<scalar>(mesh().lookupObject<volScalarField>("XLiBr")());
}
