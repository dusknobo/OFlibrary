/**
 * @file thermoProbe.C
 * @brief 独立物性/等温Schrage辅助核验程序，不推进CFD时间。
 * 输入每行T[K]、X[kg/kg]、p[Pa]、适应系数；输出物性、通量及配对源。
 * 面积/体积比固定为1 m⁻¹用于内核归一化核验，不代表实际设备几何。
 */
#include "LiBrThermoTable.H"
#include <iostream>
#include <iomanip>

int main(int argc, char** argv)
{
    try
    {
        if (argc != 2)
            throw std::runtime_error("Usage: thermoProbe table.dat < T_X_p_accommodation_rows");
        librSpray::Table table(argv[1]);
        double T, X, p, coefficient;
        std::cout << std::setprecision(17);
        while (std::cin >> T >> X >> p >> coefficient)
        {
            const auto s = table.evaluate(T, X);
            const double flux = librSpray::schrageFlux(s, T, p, coefficient);
            const double areaPerVolume = 1.0; // m^-1, source-kernel verification normalization
            const auto source = librSpray::phaseSources(flux * areaPerVolume, s.latentWater);
            std::cout << s.rho << ' ' << s.cp << ' ' << s.mu << ' ' << s.pEq << ' ' << s.latentWater
                      << ' ' << flux << ' ' << source.liquidMass << ' ' << source.vapourMass << ' '
                      << source.saltMass << ' ' << source.sensibleHeat << ' ' << s.h << ' '
                      << s.dh_dX << '\n';
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
