// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file reaction_audit.cpp
 * @brief 局部隐式反应预测器的独立核验入口。
 * 输入依次为ml、mg、vg、area、Tl、Tg、X、dt、beta、maxFraction，均为SI。
 * 输出dm、hv、Tl、Tg、X、质量残差、限幅标记和迭代次数，不推进整机CFD。
 */
#include "PSprayFilmReaction.H"
#include <iostream>
#include <iomanip>

int main(int argc, char** argv)
{
    if (argc != 2)
        return 2;
    librSpray::Table table(argv[1]);
    librSpray::FilmReactionInput a;
    a.cpv = 1860.;
    a.href = 2547276.823116537;
    a.Tref = 298.15;
    std::cout << std::setprecision(17);
    while (std::cin >> a.ml >> a.mg >> a.vg >> a.area >> a.Tl >> a.Tg >> a.X >> a.dt >> a.beta >>
           a.maxFraction)
    {
        try
        {
            auto z = librSpray::filmReaction(table, a);
            std::cout << z.dm << " " << z.hv << " " << z.Tl << " " << z.Tg << " " << z.X << " "
                      << z.residual << " " << z.limited << " " << z.iterations << "\n";
        }
        catch (const std::exception& e)
        {
            std::cerr << e.what() << "\n";
            return 1;
        }
    }
}
