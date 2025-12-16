#pragma once
#include <complex>
#include <vector>
#include <cmath>

// struct Params {
//     double M;
//     double a;
//     int s;
//     int l, m;
//     double omega;
// };
struct Params {
    double M;
    double a;
    double omega;

    int s;      // Teukolsky spin weight (0, -1/2, -1, -2)
    double ell; // orbital mode (can be half-integer)
    int m;      // azimuthal mode

    int polarizations; // g: number of degrees of freedom
};


inline double Delta(double r, const Params& p) {
    return r*r - 2*p.M*r + p.a*p.a;
}

inline double K(double r, const Params& p) {
    return (r*r + p.a*p.a)*p.omega - p.a*p.m;
}

// inline double lambda(const Params& p) {
//     return p.l*(p.l+1) - p.s*(p.s+1);
// }
inline double lambda(const Params& p) {
    return p.ell*(p.ell+1.0) - p.s*(p.s+1.0);
}


struct TeukolskySolution {
    std::vector<std::complex<double>> R;
    std::vector<double> r;
    std::complex<double> Z_in, Z_out;
};

class RadialTeukolsky {
public:
    double r_min, r_max;
    int N;

    RadialTeukolsky(double rmin, double rmax, int Npts)
        : r_min(rmin), r_max(rmax), N(Npts) {}

    TeukolskySolution integrate(const Params& P);
};

