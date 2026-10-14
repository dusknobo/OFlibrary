// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file librSolutionThermo.C
 * @brief 将librSolutionThermo注册到basic/fluid/rhoFluid三个热物性选择表。
 * 类型名与算例physicalProperties中选择的热物性名称保持一致。
 */
#include "librSolutionThermo.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
defineTypeNameAndDebug(librSolutionThermo, 0);
addToRunTimeSelectionTable(basicThermo, librSolutionThermo, fvMesh);
addToRunTimeSelectionTable(fluidThermo, librSolutionThermo, fvMesh);
addToRunTimeSelectionTable(rhoFluidThermo, librSolutionThermo, fvMesh);
}
