/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::pSpraySolutionCloud

Description
    Cloud class to introduce reacting spray parcels

\*---------------------------------------------------------------------------*/

#ifndef pSpraySolutionCloud_H
#define pSpraySolutionCloud_H

/**
 * @file pSpraySolutionCloud.H
 * @brief 以pSpraySolutionParcel组成原生Momentum/Thermo/Reacting/Spray云。
 * 两侧液体使用同一云类型、不同入口X和T；子模型是否启用由算例字典决定。
 * LiBr保存在自定义parcel属性中，气相仍只有水蒸气。
 */
#include "ParcelCloudBase.H"
#include "MomentumCloud.H"
#include "ThermoCloud.H"
#include "ReactingCloud.H"
#include "SprayCloud.H"
#include "ParcelCloud.H"
#include "pSpraySolutionParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// 逐层保留原生cloud继承链，便于核对各层parcel能力。
// clang-format off
typedef ParcelCloud
    <
        SprayCloud
        <
            ReactingCloud
            <
                ThermoCloud
                <
                    MomentumCloud
                    <
                        ParcelCloudBase
                        <
                            pSpraySolutionParcel
                        >
                    >
                >
            >
        >
    > pSpraySolutionCloud;
// clang-format on
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
