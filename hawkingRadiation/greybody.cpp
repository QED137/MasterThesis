#include "greybody.hpp"
#include <complex>

using cdouble = std::complex<double>;

double computeGreybody(const TeukolskySolution& sol, const Params& P) {
    int N = sol.r.size();
    int i0 = N - 50;               // choose a region near infinity
    int i1 = N - 1;

    cdouble A(0,0), B(0,0);

    for (int i=i0; i<i1; i++) {
        double r = sol.r[i];
        cdouble R = sol.R[i];

        cdouble up = std::exp(cdouble(0, +P.omega*r));
        cdouble dn = std::exp(cdouble(0, -P.omega*r));

        // fit R = A * up + B * dn
        A += R * std::conj(up);
        B += R * std::conj(dn);
    }

    double Zin = std::norm(B);
    double Zout = std::norm(A);

    double Gamma = 1.0 - Zout/Zin;
    if (Gamma < 0) Gamma = 0;
    return Gamma;
}

