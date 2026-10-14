/**
 * @file table_audit.cpp
 * @brief Python库存核验调用的独立表查询工具。
 * argv[1]为物性表路径；标准输入每行T[K]、X[kg/kg]，
 * 输出rho[kg/m³]、e[J/kg]、pEq[Pa]，以17位精度供累计收支核验。
 */
#include "LiBrThermoTable.H"
#include <iostream>
#include <iomanip>

int main(int argc, char** argv)
{
    librSpray::Table tab(argv[1]);
    double T, X;
    std::cout << std::setprecision(17);
    while (std::cin >> T >> X)
    {
        const auto s = tab.evaluate(T, X);
        std::cout << s.rho << " " << s.h - s.pEq / s.rho << " " << s.pEq << "\n";
    }
}
