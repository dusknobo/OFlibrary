// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file PSpraySolutionThermo.H
 * @brief 液滴与撞膜模型共用的物性表缓存和压力修正焓。
 * 与液膜采用相同表和能量参考，避免转移时人为产生混合/相变热量。
 */
#ifndef PSpraySolutionThermo_H
#define PSpraySolutionThermo_H
#include "LiBrThermoTable.H"
#include "dictionary.H"
#include "fileName.H"
#include <map>
#include <memory>

namespace Foam
{
/// 按展开后的solutionTable路径缓存只读表；每个MPI进程具有自己的缓存。
inline const librSpray::Table& pSprayTable(const dictionary& dict)
{
    fileName path(dict.lookup<fileName>("solutionTable"));
    path.expand();
    static std::map<std::string, std::shared_ptr<librSpray::Table>> tables;
    const std::string key(path.c_str());
    auto it = tables.find(key);
    if (it == tables.end())
        it = tables.emplace(key, std::make_shared<librSpray::Table>(key)).first;
    return *it->second;
}

/// 比焓[J/kg]：h(p)=hEq+(p-pEq)/rho，p为绝对压力[Pa]。
inline scalar pSprayEnthalpy(const librSpray::Table& table, scalar p, scalar T, scalar X)
{
    const auto a = table.evaluate(T, X);
    return a.h + (p - a.pEq) / a.rho;
}
}
#endif
