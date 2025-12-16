#include "teukolsky_radial.hpp"
#include <complex>
#include <vector>

using cdouble = std::complex<double>;

static cdouble RHS(double r, cdouble R, cdouble dR, const Params& P) {
    double D = Delta(r, P);
    double k = K(r, P);
    double lam = lambda(P);
    int s = P.s;

    cdouble V =
        (k*k - 2.0*cdouble(0,1)*double(s)*(r-P.M)*k)/D +
        4.0*cdouble(0,1)*double(s)*P.omega*r -
        lam;

    return -( (P.s+1.0)*D* dR/(r*r + P.a*P.a) + V*R );
}

TeukolskySolution RadialTeukolsky::integrate(const Params& P) {
    TeukolskySolution sol;
    sol.r.resize(N);
    sol.R.resize(N);

    double dr = (r_max - r_min)/double(N-1);

    for (int i=0; i<N; i++)
        sol.r[i] = r_min + i*dr;

    // initial boundary conditions near horizon
    double r0 = sol.r[0];
    double rp = P.M + sqrt(P.M*P.M - P.a*P.a);
    double wm = P.omega - P.m*P.a/(2*P.M*rp);

    cdouble R = pow(Delta(r0,P), -double(P.s)/2.0) *
                std::exp(cdouble(0,-wm*r0));
    cdouble dR = cdouble(0,0);

    sol.R[0] = R;

    // RK4 integration
    for (int i=0; i<N-1; i++) {
        double r = sol.r[i];
        cdouble k1 = dr * dR;
        cdouble l1 = dr * RHS(r, R, dR, P);

        cdouble k2 = dr * (dR + 0.5*l1);
        cdouble l2 = dr * RHS(r+0.5*dr, R+0.5*k1, dR+0.5*l1, P);

        cdouble k3 = dr * (dR + 0.5*l2);
        cdouble l3 = dr * RHS(r+0.5*dr, R+0.5*k2, dR+0.5*l2, P);

        cdouble k4 = dr * (dR + l3);
        cdouble l4 = dr * RHS(r+dr, R+k3, dR+l3, P);

        R  += (k1 + 2.0*(k2+k3) + k4)/6.0;
        dR += (l1 + 2.0*(l2+l3) + l4)/6.0;

        sol.R[i+1] = R;
    }

    return sol;
}

