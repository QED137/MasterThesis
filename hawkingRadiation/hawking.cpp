#include "hawking.hpp"
#include <cmath>

double hawkingTemperature(double M, double a){
    double rp = M + sqrt(M*M - a*a);
    double rm = M - sqrt(M*M - a*a);
    return (rp - rm)/(4*M_PI*(rp*rp + a*a));
}

// double hawkingSpectrum(const Params& P, double Gamma){
//     double TH = hawkingTemperature(P.M, P.a);
//     double x = (P.omega - P.m*P.a/(2*P.M*(P.M+sqrt(P.M*P.M-P.a*P.a))))/TH;

//     double denominator = exp(x) - pow(-1.0,2*P.s);
//     return (P.omega * Gamma) / (2*M_PI*denominator);
// }

double hawkingSpectrum(const Params& P, double Gamma){
    double TH = hawkingTemperature(P.M, P.a);

    double OmegaH = P.a / (2*P.M*(P.M + sqrt(P.M*P.M - P.a*P.a)));
    double x = (P.omega - P.m*OmegaH) / TH;

    // Fermi vs Bose
    double denom;
    if (P.s == -0.5)  // fermion
        denom = exp(x) + 1.0;
    else
        denom = exp(x) - 1.0;

    return P.polarizations * (P.omega * Gamma) / (2*M_PI * denom);
}
