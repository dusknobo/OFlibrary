/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

Class
    Foam::fv::pSprayFilmCloudTransfer

Description
    Film<->cloud transfer model

Usage
    Example usage:
    \verbatim
    pSprayFilmCloudTransfer
    {
        type    pSprayFilmCloudTransfer;
    }
    \endverbatim

SourceFiles
    pSprayFilmCloudTransfer.C

\*---------------------------------------------------------------------------*/

/**
 * @file pSprayFilmCloudTransfer.H
 * @brief 液膜侧的液滴捕集/回滴源项与同一步广延转移账本。
 * cloud入账为kg、kg、kg·m/s、J；施加方程源时才除以单元体积和dt。
 * 回滴质量、盐、动量与内能在同一时间步冻结，供新parcel和膜侧负源共用。
 */
#ifndef pSprayFilmCloudTransfer_H
#define pSprayFilmCloudTransfer_H

#include "fvModel.H"
#include "isothermalFilm.H"
#include "ejectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class pSprayFilmCloudTransfer Declaration
\*---------------------------------------------------------------------------*/

class pSprayFilmCloudTransfer : public fvModel
{
    // Private Data

    //- The film model
    const solvers::isothermalFilm& film_;

    // Transfers from cloud

    //- Switch to indicate the cloud fields have been transferred
    bool cloudFieldsTransferred_;

    // 两种云在气相面上的撞膜增量共同累加；每个timeIndex只清零一次。
    scalarField massFromCloud_;
    scalarField saltFromCloud_;
    label lastCloudReset_ = -1;
    word ejectionCloud_;
    scalar maxEjectionFraction_;
    // 单元本步回滴广延量，不是速率；重启审计需区分此字段和phase*Rate。
    volScalarField ejectionMass_, ejectionSalt_, ejectionEnergy_;
    volVectorField ejectionMomentum_;
    /// 预留相变水量、限制单步回滴比例，并排除原生会删除的过小parcel。
    void prepareEjection();
    vectorField momentumFromCloud_;
    scalarField energyFromCloud_;

    //- Switch to ensure the ejection rate is not updated until
    //  the last set of ejected parcels have been included in the cloud
    bool correctEjection_;

    //- Pointer to the ejection model
    autoPtr<ejectionModel> ejection_;

    // Private Member Functions

    //- Cloud to film rate transfer function
    template <class Type>
    tmp<VolInternalField<Type>> inline CloudToFilmTransferRate(const Field<Type>& prop,
                                                               const dimensionSet& dimProp) const;

    //- Film to cloud transfer function
    template <class Type>
    tmp<Field<Type>> inline filmToCloudTransfer(const VolInternalField<Type>& prop) const;

public:
    //- Runtime type information
    TypeName("pSprayFilmCloudTransfer");

    // Constructors

    //- Construct from explicit source name and mesh
    pSprayFilmCloudTransfer(const word& sourceName,
                            const word& modelType,
                            const fvMesh& mesh,
                            const dictionary& dict);

    //- Disallow default bitwise copy construction
    pSprayFilmCloudTransfer(const pSprayFilmCloudTransfer&) = delete;

    // Member Functions

    // Checks

    //- Return the list of fields for which the option adds source term
    //  to the transport equation
    virtual wordList addSupFields() const;

    // Correct

    //- 本实现不在correct中求解液膜；回滴在云缓存膜状态时准备。
    virtual void correct();

    // Sources

    //- Add source to phase continuity equation
    virtual void
    addSup(const volScalarField& rho, const volScalarField& alpha, fvMatrix<scalar>& eqn) const;

    //- Add source to phase energy equation
    virtual void addSup(const volScalarField& alpha,
                        const volScalarField& rho,
                        const volScalarField& he,
                        fvMatrix<scalar>& eqn) const;

    //- Add source to mixture momentum equation
    virtual void addSup(const volScalarField& alpha,
                        const volScalarField& rho,
                        const volVectorField& U,
                        fvMatrix<vector>& eqn) const;

    // Transfer from cloud

    //- Reset the fields accumulated cloud transfer fields
    /// 同一步第二种云再次访问时保留第一种云已累计的转移。
    void resetFromCloudFields();

    //- Transfer parcel properties from cloud to the film
    //  accumulated fields
    void parcelFromCloud(const label facei,
                         const scalar mass,
                         const vector& momentum,
                         const scalar energy,
                         const scalar saltMass);

    // Transfer to cloud

    //- Return true if the film is ejecting to the cloud
    /// 只向ejectionCloud指定的一种云返回液滴，避免两种云重复取走膜质量。
    bool ejecting(const word& cloudName) const;
    tmp<Field<scalar>> XToCloud() const;

    //- Transfer the ejected mass to the cloud
    tmp<Field<scalar>> ejectedMassToCloud() const;

    //- Transfer the ejected droplet diameter to the cloud
    tmp<Field<scalar>> ejectedDiameterToCloud() const;

    //- Transfer the film delta field to the cloud
    tmp<Field<scalar>> deltaToCloud() const;

    //- Transfer the film velocity field to the cloud
    tmp<Field<vector>> UToCloud() const;

    //- Transfer the film density field to the cloud
    tmp<Field<scalar>> rhoToCloud() const;

    //- Transfer the film temperature field to the cloud
    tmp<Field<scalar>> TToCloud() const;

    //- Transfer the film heat capacity field to the cloud
    tmp<Field<scalar>> CpToCloud() const;

    // Mesh changes

    //- Update topology using the given map
    virtual void topoChange(const polyTopoChangeMap&);

    //- Update from another mesh using the given map
    virtual void mapMesh(const polyMeshMap&);

    //- Redistribute or update using the given distribution map
    virtual void distribute(const polyDistributionMap&);

    //- Update for mesh motion
    virtual bool movePoints();

    // Member Operators

    //- Disallow default bitwise assignment
    void operator=(const pSprayFilmCloudTransfer&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
