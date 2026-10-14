/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
 * @file makePSpraySolutionSubmodels.C
 * @brief 为自定义cloud实例化并注册原生及溶液扩展子模型。
 * 注册某项能力不等于在当前算例中启用它；parcel会检查并禁止当前不支持的
 * 破碎、雾化、飞行相变及辐射。膜转移模板在独立FilmTransfer库中注册。
 */
#include "pSpraySolutionCloud.H"

#include "makeParcelCloudFunctionObjects.H"

// Momentum
#include "makeThermoParcelForces.H"
#include "makeParcelDispersionModels.H"
#include "makeSprayParcelInjectionModels.H"
#include "makeParcelPatchInteractionModels.H"
#include "makeSprayParcelStochasticCollisionModels.H"

// Thermodynamic
#include "makeParcelHeatTransferModels.H"
#include "makeParcelCompositionModels.H"

// Reacting
#include "makeReactingParcelPhaseChangeModels.H"
#include "makeParcelSurfaceFilmModels.H"

// Spray
#include "DistortedSphereDragForce.H"
#include "makeSprayParcelAtomisationModels.H"
#include "makeSprayParcelBreakupModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeParcelCloudFunctionObjects(pSpraySolutionCloud);

// Momentum sub-models
makeThermoParcelForces(pSpraySolutionCloud);
makeParcelDispersionModels(pSpraySolutionCloud);
makeSprayParcelInjectionModels(pSpraySolutionCloud);
makeParcelPatchInteractionModels(pSpraySolutionCloud);
makeSprayParcelStochasticCollisionModels(pSpraySolutionCloud);

// Thermo sub-models
makeParcelHeatTransferModels(pSpraySolutionCloud);
makeParcelCompositionModels(pSpraySolutionCloud);

// Reacting sub-models
makeReactingParcelPhaseChangeModels(pSpraySolutionCloud);
makeParcelSurfaceFilmModels(pSpraySolutionCloud);

// Spray sub-models
makeParticleForceModelType(DistortedSphereDragForce, pSpraySolutionCloud);
makeSprayParcelAtomisationModels(pSpraySolutionCloud);
makeSprayParcelBreakupModels(pSpraySolutionCloud);

// ************************************************************************* //

#include "PSpraySolutionFanInjection.H"
makeInjectionModelType(PSprayFlatFanInjection, pSpraySolutionCloud);
makeInjectionModelType(PSpraySolutionFanInjection, pSpraySolutionCloud);

#include "PSprayBudgetInteraction.H"
makePatchInteractionModelType(PSprayBudgetInteraction, pSpraySolutionCloud);
