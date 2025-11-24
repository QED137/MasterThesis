#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

// -----------------------------
// Data structures
// -----------------------------

struct ComplexPsi {
    double real;
    double imag;
};

// -----------------------------
// Physical Parameters
// -----------------------------

struct Params {
    double M = 1.0;       // Mass
    double spin = 0.5;    // a/M
    int s = -2;           // spin weight
    int l = 2, m = 2;     // mode numbers
    double omega = 0.5;   // frequency
};

// -----------------------------
// Simulation grid
// -----------------------------

struct Grid {
    int Nr = 200;
    double rMin = 1.5;
    double rMax = 50.0;
    double dr;

    vector<double> r;
    vector<double> potential;
    vector<ComplexPsi> psi;
};

// ------------------------------------------------------
// Compute effective radial Teukolsky potential
// ------------------------------------------------------

double teukolsky_potential(double r, const Params& P) {
    double M = P.M;
    double a = P.spin * M;
    double s = P.s;

    double Delta = r*r - 2*M*r + a*a;
    double K = (r*r + a*a) * P.omega - a * P.m;
    double lambda = P.l * (P.l + 1) - s * (s + 1);

    double term1 = lambda;
    double term2 = (2 * s * (r - M) / Delta) * (K - a * P.omega);
    double term3 = (K*K) / Delta;

    return (Delta / (r*r + a*a)) * (term1 + term2 + term3);
}

// ------------------------------------------------------
// Initialize grid, potential, and wavefunction
// ------------------------------------------------------

void initialize(Grid& G, const Params& P) {
    G.dr = (G.rMax - G.rMin) / (G.Nr - 1);
    G.r.resize(G.Nr);
    G.potential.resize(G.Nr);
    G.psi.resize(G.Nr);

    // Build r-grid
    for (int i = 0; i < G.Nr; i++) {
        G.r[i] = G.rMin + i * G.dr;
    }

    // Potential
    for (int i = 0; i < G.Nr; i++) {
        G.potential[i] = teukolsky_potential(G.r[i], P);
    }

    // Initial Gaussian wavepacket
    double r0 = (G.rMin + G.rMax) * 0.5;
    double sigma = 5.0;
    double k0 = P.omega;

    for (int i = 0; i < G.Nr; i++) {
        double gaussian = exp(-pow((G.r[i] - r0) / sigma, 2));
        double phase = k0 * (G.r[i] - r0);

        G.psi[i].real = gaussian * cos(phase);
        G.psi[i].imag = gaussian * sin(phase);
    }

    cout << "Initialization complete." << endl;
}

// ------------------------------------------------------
// Time evolution step
// ------------------------------------------------------

void evolve(Grid& G, const Params& P, double dt) {
    vector<ComplexPsi> newPsi = G.psi;

    for (int i = 1; i < G.Nr - 1; i++) {
        double d2_real =
            (G.psi[i+1].real - 2*G.psi[i].real + G.psi[i-1].real) / (G.dr * G.dr);
        double d2_imag =
            (G.psi[i+1].imag - 2*G.psi[i].imag + G.psi[i-1].imag) / (G.dr * G.dr);

        double V = G.potential[i];

        // i ∂ψ/∂t = (-∂²ψ/∂r*² + V) ψ
        newPsi[i].real =
            G.psi[i].real - dt * (-d2_imag + V * G.psi[i].imag);

        newPsi[i].imag =
            G.psi[i].imag + dt * (-d2_real + V * G.psi[i].real);
    }

    // Absorbing boundaries
    int dampingWidth = 10;
    for (int i = 0; i < dampingWidth; i++) {
        double damp = exp(-pow((dampingWidth - i) / 5.0, 2));
        newPsi[i].real *= damp;
        newPsi[i].imag *= damp;
        newPsi[G.Nr-1 - i].real *= damp;
        newPsi[G.Nr-1 - i].imag *= damp;
    }

    G.psi = newPsi;
}

// ------------------------------------------------------
// Save wavefunction to CSV
// ------------------------------------------------------

void save_csv(const Grid& G, int step) {
    string filename = "psi_step_" + to_string(step) + ".csv";
    ofstream f(filename);

    f << "r,real,imag\n";
    for (int i = 0; i < G.Nr; i++) {
        f << G.r[i] << ","
          << G.psi[i].real << ","
          << G.psi[i].imag << "\n";
    }

    f.close();
    cout << "Saved: " << filename << endl;
}

// ------------------------------------------------------
// Main Simulation Loop
// ------------------------------------------------------

int main() {
    Params P;
    Grid G;

    initialize(G, P);

    double dt = 0.05;
    int totalSteps = 2000;

    for (int step = 0; step < totalSteps; step++) {
        evolve(G, P, dt);

        if (step % 100 == 0) {
            cout << "Step " << step << endl;
            save_csv(G, step);
        }
    }

    cout << "Simulation complete." << endl;
    return 0;
}

