#ifndef BLACK_HOLE_H
#define BLACK_HOLE_H

#include <complex>
#include <vector>
#include <functional>

using Complex = std::complex<double>;

// Particle type enum
enum class ParticleType {
    BOSON,      // Integer spin: 0, 1, 2
    FERMION     // Half-integer spin: 1/2, 3/2, ...
};

// Configuration structure
struct Config {
    double M;           // Black hole mass
    double w;           // Frequency
    int l;              // Angular momentum number
    int s;              // Spin (0, 1, 2 for bosons; use s_half for fermions)
    double s_half;      // Spin for fermions (0.5, 1.5, etc.)
    ParticleType type;  // BOSON or FERMION
    double r_ex;        // Extraction radius
    double r_horizon;   // Horizon radius
    int max_steps;      // Maximum ODE solver steps
    
    // Constructor with defaults for bosons
    Config() : M(1.0), w(0.1), l(0), s(0), s_half(0.0), 
               type(ParticleType::BOSON), r_ex(100.0), 
               r_horizon(2.0), max_steps(1000000) {}
};

// Asymptotic expansion coefficients
struct AsymptoticCoeffs {
    Complex a_in, b_in, c_in, d_in;   // Ingoing wave coefficients
    Complex a_out, b_out, c_out, d_out; // Outgoing wave coefficients
};

// ODE System for perturbation equations
class PerturbationODE {
public:
    PerturbationODE(const Config& cfg);
    
    // Schwarzschild case - bosons
    void schwarzschild_equation(double r, const std::vector<Complex>& y, 
                               std::vector<Complex>& dydr);
    
    // Schwarzschild case - fermions (Dirac equation)
    void schwarzschild_fermion_equation(double r, const std::vector<Complex>& y,
                                       std::vector<Complex>& dydr);
    
    // Extreme Kerr case
    void extreme_kerr_equation(double r, const std::vector<Complex>& y,
                              std::vector<Complex>& dydr);
    
    // Compute horizon boundary conditions
    Complex horizon_bc_value(double r);
    Complex horizon_bc_derivative(double r);
    
    // Fermion boundary conditions
    Complex horizon_bc_value_fermion(double r);
    Complex horizon_bc_derivative_fermion(double r);
    
    // Compute asymptotic coefficients
    AsymptoticCoeffs compute_asymptotic_coeffs();
    AsymptoticCoeffs compute_asymptotic_coeffs_fermion();
    
private:
    Config config_;
    Complex compute_horizon_expansion_a();
    Complex compute_horizon_expansion_b();
    Complex compute_horizon_expansion_a_fermion();
    Complex compute_horizon_expansion_b_fermion();
};

// ODE Solver using Runge-Kutta 4th order
class ODESolver {
public:
    using DerivFunc = std::function<void(double, const std::vector<Complex>&, 
                                        std::vector<Complex>&)>;
    
    ODESolver(DerivFunc deriv, double r_start, double r_end, 
             const std::vector<Complex>& y0, int max_steps = 1000000);
    
    bool solve();
    Complex value_at_end() const { return y_end_[0]; }
    Complex derivative_at_end() const { return y_end_[1]; }
    
private:
    void rk4_step(double r, double h, std::vector<Complex>& y);
    
    DerivFunc deriv_;
    double r_start_, r_end_;
    std::vector<Complex> y0_;
    std::vector<Complex> y_end_;
    int max_steps_;
};

// Main calculator class
class AbsorptionCalculator {
public:
    AbsorptionCalculator(const Config& cfg);
    
    // Compute absorption cross-section (returns 1 - |reflection_coeff|^2)
    double compute_absorption();
    
    // Compute thermal flux at given frequency
    double compute_thermal_flux(double temperature, int l_max);
    
    // Get reflection and transmission coefficients
    std::pair<Complex, Complex> get_coefficients() const {
        return {C_in_, C_out_};
    }
    
private:
    Config config_;
    Complex C_in_;  // Incoming wave amplitude
    Complex C_out_; // Outgoing wave amplitude
    
    bool solve_radial_equation();
    bool solve_radial_equation_fermion();
    void match_asymptotic_solution();
};

// Neutrino calculator using Page's analytical formulas
class NeutrinoCalculator {
public:
    NeutrinoCalculator(double M = 1.0);
    
    // Compute neutrino absorption using Page's analytical formula
    double compute_absorption_analytical(double omega, int l);
    
    // Compute full neutrino spectrum
    double compute_spectrum(double omega, int l_max);
    
    // Get absorption for all l modes
    std::vector<double> compute_all_modes(double omega, int l_max);
    
private:
    double M_;
    
    // Page's Gamma function for neutrinos (Eq. 16)
    Complex compute_gamma_neutrino(double omega, int l, int m, int p);
    
private:
    Config config_;
    Complex C_in_;  // Incoming wave amplitude
    Complex C_out_; // Outgoing wave amplitude
    
    bool solve_radial_equation();
    void match_asymptotic_solution();
};

#endif // BLACK_HOLE_H
