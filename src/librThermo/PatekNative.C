// Patek-Klomfar original correlation code preserved by CoolProp.
// Source: CoolProp src/Backends/Incompressible/IncompressibleLibrary.cpp.
// DOI: 10.1016/j.ijrefrig.2005.10.007. See LICENSE.CoolProp and provenance.
#include <cmath>
namespace patekNative {
using std::pow; using std::exp; using std::log;
double ps_H2O(double); double rho_H2O(double); double cp_H2O(double);
double h_H2O(double); double s_H2O(double);
double const M_H2O  = 0.018015268; /* kg/mol, molar mass of H2O */
double const M_LiBr = 0.08685; /* kg/mol, molar mass of LiBr */
double const T0     = 221; /* K, constant */

/* Critical point of H2O */
double const Tc_H2O   = 647.096; /* K, temperature  */
double const pc_H2O   = 22.064; /* MPa, pressure */
double const rhoc_H2O = 17873; /* mol/m^3 (322 kg/m^3), molar density */
double const hc_H2O   = 37548.5; /* J/mol, molar enthalpy */
double const sc_H2O   = 79.3933; /* J/(mol.K) molar entropy*/

/*Triple point of H2O */
double const Tt_H2O   = 273.16; /* K, temperature */
double const cpt_H2O  = 76.0226; /* J/(mol.K), molar isobaric heat capacity*/

double ps_mix(double T, double x)
/* Equation (1) */
{
    static double m[8] = { 3.0, 4.0, 4.0, 8.0, 1.0, 1.0, 4.0, 6.0 };
    static double n[8] = { 0.0, 5.0, 6.0, 3.0, 0.0, 2.0, 6.0, 0.0 };
    static double t[8] = { 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 };
    static double a[8] = { -2.41303e2, 1.91750e7, -1.75521e8, 3.25430e7,
            3.92571e2, -2.12626e3, 1.85127e8, 1.91216e3 };
    double tau, suma = 0.0;
    int i;

    tau = T / Tc_H2O;
    for (i = 0; i <= 7; i++)
        suma += a[i] * pow(x, m[i]) * pow(0.4 - x, n[i]) * pow(tau, t[i]);
    return (ps_H2O(T - suma));

} /* end function ps_mix */

double rho_mix(double T, double x)
/* Equation (2) */
{
    static double m[2] = { 1.0, 1.0 };
    static double n[2] = { 0.0, 6.0 };
    static double a[2] = { 1.746, 4.709 };

    double tau, suma = 0.0;
    int i;

    tau = T / Tc_H2O;
    for (i = 0; i <= 1; i++)
        suma += a[i] * pow(x, m[i]) * pow(tau, n[i]);

    return ((1.0 - x) * rho_H2O(T) + rhoc_H2O * suma);

} /* end function rho_mix */

double cp_mix(double T, double x)
/* Equation (3) */
{
    static double m[8] = { 2.0, 3.0, 3.0, 3.0, 3.0, 2.0, 1.0, 1.0 };
    static double n[8] = { 0.0, 0.0, 1.0, 2.0, 3.0, 0.0, 3.0, 2.0 };
    static double t[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 3.0, 4.0 };
    static double a[8] = { -1.42094e1, 4.04943e1, 1.11135e2, 2.29980e2,
            1.34526e3, -1.41010e-2, 1.24977e-2, -6.83209e-4 };

    double tau, suma = 0.0;
    int i;

    tau = Tc_H2O / (T - T0);
    for (i = 0; i <= 7; i++)
        suma += a[i] * pow(x, m[i]) * pow(0.4 - x, n[i]) * pow(tau, t[i]);

    return ((1.0 - x) * cp_H2O(T) + cpt_H2O * suma);

} /* end function cp_mix */

double h_mix(double T, double x)
/* Equation (4) */
{
    static double m[30] = { 1.0, 1.0, 2.0, 3.0, 6.0, 1.0, 3.0, 5.0, 4.0,
            5.0, 5.0, 6.0, 6.0, 1.0, 2.0, 2.0, 2.0, 5.0, 6.0, 7.0, 1.0, 1.0,
            2.0, 2.0, 2.0, 3.0, 1.0, 1.0, 1.0, 1.0 };

    static double n[30] = { 0.0, 1.0, 6.0, 6.0, 2.0, 0.0, 0.0, 4.0, 0.0,
            4.0, 5.0, 5.0, 6.0, 0.0, 3.0, 5.0, 7.0, 0.0, 3.0, 1.0, 0.0, 4.0,
            2.0, 6.0, 7.0, 0.0, 0.0, 1.0, 2.0, 3.0 };

    static double t[30] = { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0,
            2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0,
            4.0, 4.0, 4.0, 4.0, 5.0, 5.0, 5.0, 5.0 };

    static double a[30] = { 2.27431, -7.99511, 3.85239e2, -1.63940e4,
            -4.22562e2, 1.13314e-1, -8.33474, -1.73833e4, 6.49763,
            3.24552e3, -1.34643e4, 3.99322e4, -2.58877e5, -1.93046e-3,
            2.80616, -4.04479e1, 1.45342e2, -2.74873, -4.49743e2,
            -1.21794e1, -5.83739e-3, 2.33910e-1, 3.41888e-1, 8.85259,
            -1.78731e1, 7.35179e-2, -1.79430e-4, 1.84261e-3, -6.24282e-3,
            6.84765e-3 };

    double tau, suma = 0.0;
    int i;

    tau = Tc_H2O / (T - T0);
    for (i = 0; i <= 29; i++)
        suma += a[i] * pow(x, m[i]) * pow(0.4 - x, n[i]) * pow(tau, t[i]);

    return ((1.0 - x) * h_H2O(T) + hc_H2O * suma);

} /* end function h_mix */

double s_mix(double T, double x)
/* Equation (5) */
{
    static double m[29] = { 1.0, 1.0, 2.0, 3.0, 6.0, 1.0, 3.0, 5.0, 1.0,
            2.0, 2.0, 4.0, 5.0, 5.0, 6.0, 6.0, 1.0, 3.0, 5.0, 7.0, 1.0, 1.0,
            1.0, 2.0, 3.0, 1.0, 1.0, 1.0, 1.0 };

    static double n[29] = { 0.0, 1.0, 6.0, 6.0, 2.0, 0.0, 0.0, 4.0, 0.0,
            0.0, 4.0, 0.0, 4.0, 5.0, 2.0, 5.0, 0.0, 4.0, 0.0, 1.0, 0.0, 2.0,
            4.0, 7.0, 1.0, 0.0, 1.0, 2.0, 3.0 };

    static double t[29] = { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0,
            2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0,
            4.0, 4.0, 4.0, 5.0, 5.0, 5.0, 5.0 };

    static double a[29] = { 1.53091, -4.52564, 6.98302e+2, -2.16664e+4,
            -1.47533e+3, 8.47012e-2, -6.59523, -2.95331e+4, 9.56314e-3,
            -1.88679e-1, 9.31752, 5.78104, 1.38931e+4, -1.71762e+4,
            4.15108e+2, -5.55647e+4, -4.23409e-3, 3.05242e+1, -1.67620,
            1.48283e+1, 3.03055e-3, -4.01810e-2, 1.49252e-1, 2.59240,
            -1.77421e-1, -6.99650e-5, 6.05007e-4, -1.65228e-3, 1.22966e-3 };

    double tau, suma = 0.0;
    int i;

    tau = Tc_H2O / (T - T0);
    for (i = 0; i <= 28; i++)
        suma += a[i] * pow(x, m[i]) * pow(0.4 - x, n[i]) * pow(tau, t[i]);

    return ((1.0 - x) * s_H2O(T) + sc_H2O * suma);

} /* end function s_mix */

double ps_H2O(double T)
/* Equation (28) */
{
    static double a[7] = { 0.0, -7.85951783, 1.84408259, -11.7866497,
            22.6807411, -15.9618719, 1.80122502 };

    double tau, ps;

    tau = 1 - T / Tc_H2O;

    ps = pc_H2O
            * exp(
                    Tc_H2O / T
                            * (a[1] * tau + a[2] * pow(tau, 1.5)
                                    + a[3] * pow(tau, 3.0)
                                    + a[4] * pow(tau, 3.5)
                                    + a[5] * pow(tau, 4.0)
                                    + a[6] * pow(tau, 7.5)));

    return (ps * 1.0e6);

} /* end function ps_H2O */

double rho_H2O(double T)
/* Equation (29) */
{
    static double b[7] = { 0.0, 1.99274064, 1.09965342, -0.510839303,
            -1.75493479, -45.5170352, -6.7469445e5 };
    double theta, rho;

    theta = 1.0 - T / Tc_H2O;

    rho = rhoc_H2O
            * (1.0 + b[1] * pow(theta, 1.0 / 3.0)
                    + b[2] * pow(theta, 2.0 / 3.0)
                    + b[3] * pow(theta, 5.0 / 3.0)
                    + b[4] * pow(theta, 16.0 / 3.0)
                    + b[5] * pow(theta, 43.0 / 3.0)
                    + b[6] * pow(theta, 110.0 / 3.0));

    return (rho);
} /* end function rho_H2O */

double cp_H2O(double T)
/* Equation (30) */
{
    static double a[5] =
            { 1.38801, -2.95318, 3.18721, -0.645473, 9.18946e5 };
    static double b[5] = { 0.0, 2.0, 3.0, 6.0, 34.0 };
    static double c[5] = { 0.0, 2.0, 3.0, 5.0, 0.0 };

    double suma = 0.0;
    int i;

    for (i = 0; i <= 4; i++)
        suma += a[i] * exp(b[i] * log(1.0 - T / Tc_H2O))
                * exp(c[i] * log(T / Tt_H2O));

    return (cpt_H2O * suma);

} /* end function cp_H2O */

double h_H2O(double T)
/* Equation (31) */
{
    static double a[4] = { -4.37196e-1, 3.03440e-1, -1.29582, -1.76410e-1 };
    static double alpha[4] = { 1.0 / 3.0, 2.0 / 3.0, 5.0 / 6.0, 21.0 / 6.0 };

    double suma = 0.0;
    int i;

    for (i = 0; i <= 3; i++)
        suma += a[i] * exp(alpha[i] * log(1.0 - T / Tc_H2O));

    return (hc_H2O * (1.0 + suma));

} /* end function h_H2O */

double s_H2O(double T)
/* Equation (32)  */
{
    static double a[4] = { -3.34112e-1, -8.47987e-1, -9.11980e-1, -1.64046 };
    static double alpha[4] = { 1.0 / 3.0, 3.0 / 3.0, 8.0 / 3.0, 24.0 / 3.0 };

    double suma = 0.0;
    int i;

    for (i = 0; i <= 3; i++)
        suma += a[i] * exp(alpha[i] * log(1.0 - T / Tc_H2O));

    return (sc_H2O * (1.0 + suma));

} /* end function s_H2O */


/** Finished with the code from the paper. Now we need to
 *  convert the molar values to mass-based units. */
double massToMole(double w)
/* Equation (7)  */
{
    return (w/M_LiBr)/(w/M_LiBr+(1.-w)/M_H2O);
    //return (w*M_LiBr)/(w*M_LiBr+(1.-w)*M_H2O);
}

double molarToSpecific(double w, double value)
/* Equation (7,8)  */
{
    double x = massToMole(w);
    //return w/(x*M_LiBr) * value;
    return 1. / ( x*M_LiBr + (1.-x)*M_H2O ) * value;
}

bool const debug = false;



}
extern "C" double pk_pressure(double T,double w) {return patekNative::ps_mix(T,patekNative::massToMole(w));}
extern "C" double pk_density(double T,double w) {return patekNative::rho_mix(T,patekNative::massToMole(w))/patekNative::molarToSpecific(w,1.0);}
extern "C" double pk_cp(double T,double w) {return patekNative::molarToSpecific(w,patekNative::cp_mix(T,patekNative::massToMole(w)));}
extern "C" double pk_enthalpy(double T,double w) {return patekNative::molarToSpecific(w,patekNative::h_mix(T,patekNative::massToMole(w)));}
