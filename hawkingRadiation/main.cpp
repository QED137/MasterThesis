#include <iostream>
#include <fstream>

#include "hawking.hpp"
#include "greybody.hpp"
#include "teukolsky_radial.hpp"

int main(){
    RadialTeukolsky solver(1.0001, 80.0, 6000);

    std::ofstream f("page_spectrum_all_species.csv");
    f << "species,omega,ell,m,Gamma,dEdtdw\n";

    for (auto s : {0, -1, -2}) {
        int pol = (s == 0 ? 1 : 2);

        for (double ell = std::abs(s); ell < 7; ell += 1.0) {
            for (double m = -ell; m <= ell; m += 1.0) {

                for (double w = 0.01; w < 2.0; w += 0.01) {
                    Params P;
                    P.M = 1.0;
                    P.a = 0.7;
                    P.s = s;
                    P.ell = ell;
                    P.m = m;
                    P.omega = w;
                    P.polarizations = pol;

                    auto sol = solver.integrate(P);
                    double Gamma = computeGreybody(sol, P);
                    double dEdtdw = hawkingSpectrum(P, Gamma);

                    f << s << "," << w << "," << ell << "," << m << "," << Gamma << "," << dEdtdw << "\n";
                    std::cout << "s=" << s << " ell=" << ell << " w=" << w << "\n";
                }
            }
        }
    }

    // Handle fermions (s = −1/2) separately
    for (double ell = 0.5; ell < 5.0; ell += 1.0) {
        for (double m = -ell; m <= ell; m += 1.0) {
            for (double w = 0.01; w < 2.0; w += 0.01) {

                Params P;
                P.M = 1.0;
                P.a = 0.7;
                P.s = -0.5;
                P.ell = ell;
                P.m = m;
                P.omega = w;
                P.polarizations = 2;

                auto sol = solver.integrate(P);
                double Gamma = computeGreybody(sol, P);
                double dEdtdw = hawkingSpectrum(P, Gamma);

                f << "fermion" << "," << w << "," << ell << "," << m << "," << Gamma << "," << dEdtdw << "\n";
                std::cout << "fermion ell="<<ell<<" w="<<w<<"\n";
            }
        }
    }

    std::cout << "Done.\n";
}

