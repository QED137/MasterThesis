#include "black_hole.h"
#include <cmath>
#include <stdexcept>
#include <iostream>

const double PI = 3.14159265358979323846;

PerturbationODE::PerturbationODE(const Config& cfg) : config_(cfg) {}

Complex PerturbationODE::compute_horizon_expansion_a() {
    double M = config_.M;
    double w = config_.w;
    int l = config_.l;
    int s = config_.s;
    
    Complex i(0, 1);
    
    // Convert to double for complex arithmetic
    double ld = static_cast<double>(l);
    double sd = static_cast<double>(s);
    
    // FIX: Return full complex result
    Complex numerator = i * (1.0 + ld + ld*ld - sd*sd);
    Complex denominator = 2.0 * M * (i + 4.0*M*w);
    
    return numerator / denominator;
}

Complex PerturbationODE::compute_horizon_expansion_b() {
    double M = config_.M;
    double w = config_.w;
    int l = config_.l;
    int s = config_.s;
    
    Complex i(0, 1);
    
    // Convert to double for complex arithmetic
    double ld = static_cast<double>(l);
    double sd = static_cast<double>(s);
    
    // FIX: Keep everything as COMPLEX throughout
    Complex numerator = 2.0*ld*ld*ld + ld*ld*ld*ld + 
                       ld * (2.0 - 2.0*sd*sd + 8.0*i*M*w) +
                       (ld*ld) * (3.0 - 2.0*sd*sd + 8.0*i*M*w) +
                       ((-1.0 + sd*sd)) * ((sd*sd) - 12.0*i*M*w);
    
    // FIX: Denominator is COMPLEX, not real!
    Complex denominator = 16.0 * M*M * (i + 2.0*M*w) * (i + 4.0*M*w);
    
    return -numerator / denominator;
}

void PerturbationODE::schwarzschild_equation(double r, 
                                            const std::vector<Complex>& y,
                                            std::vector<Complex>& dydr) {
    double M = config_.M;
    double w = config_.w;
    int l = config_.l;
    int s = config_.s;
    
    double f = 1.0 - 2.0*M/r;
    double f_prime = 2.0*M/(r*r);
    
    // double V_eff = f * (l*(l+1)/(r*r) + (1.0 - s*s)*2.0*M/(r*r*r));
    double V_eff = f * ((l*(l+1) + s - s*s)/(r*r) + 2.0*M*(1.0 - s*s)/(r*r*r));
    
    dydr[0] = y[1];
    dydr[1] = -(f_prime/f)*y[1] - (w*w/f - V_eff/f)/f * y[0];
}

void PerturbationODE::extreme_kerr_equation(double r,
                                           const std::vector<Complex>& y,
                                           std::vector<Complex>& dydr) {
    double M = config_.M;
    double w = config_.w;
    int k = config_.l; // Using l as k for Kerr
    
    double rho = std::sqrt(1.0 - 1.0/r);
    double lambda = 0.0; // Set mu = 0, so lambda = sqrt(1 - mu^2/w^2) = 1
    double v = 1.0;
    lambda = std::sqrt(1.0 - v*v);
    
    double fu = rho*rho / (1.0 + lambda*rho);
    double fu_prime = (2.0*rho*(1.0 + lambda*rho) - rho*rho*lambda/(2.0*rho*std::sqrt(1.0-1.0/r))) 
                     / ((1.0 + lambda*rho)*(1.0 + lambda*rho));
    
    double V_eff = w*w*(1.0 - lambda*rho)/(1.0 + lambda*rho) 
                  - k*k*rho*rho/(r*r*(1.0 + lambda*rho)*(1.0 + lambda*rho));
    
    dydr[0] = y[1];
    dydr[1] = -(fu_prime/fu)*y[1] - V_eff/(fu*fu) * y[0];
}

Complex PerturbationODE::horizon_bc_value(double r) {
    double M = config_.M;
    double w = config_.w;
    
    Complex a = compute_horizon_expansion_a();
    Complex b = compute_horizon_expansion_b();
    
    // Tortoise coordinate (ensure r > 2M to avoid log of negative)
    double rt = r + 2.0*M*std::log(std::abs(r - 2.0*M));
    Complex i(0, 1);
    
    return (1.0 + a*(r - 2.0*M) + b*(r - 2.0*M)*(r - 2.0*M)) * 
           std::exp(-i*w*rt);
}

Complex PerturbationODE::horizon_bc_derivative(double r) {
    // Numerical derivative
    double h = 1e-8;
    return (horizon_bc_value(r + h) - horizon_bc_value(r - h)) / (2.0*h);
}

AsymptoticCoeffs PerturbationODE::compute_asymptotic_coeffs() {
    AsymptoticCoeffs coeffs;
    double w = config_.w;
    int l = config_.l;
    int s = config_.s;
    double M = config_.M;
    Complex i(0, 1);
    
    // Convert integers to double for complex arithmetic
    double ld = static_cast<double>(l);
    double sd = static_cast<double>(s);
    
    // FIX: Correct formulas from Mathematica
    // Ingoing wave coefficients (Exp[-i*w*r] wave)
    // coeffs.a_in = i * (-ld - ld*ld) / (2.0*w);
    // NEW:
    coeffs.a_in = i * (-sd - ld - ld*ld) / (2.0*w);

    
    coeffs.b_in = (2.0*ld + ld*ld - 2.0*ld*ld*ld - ld*ld*ld*ld + 
                  4.0*i*M*(-1.0 + sd*sd)*w) / (8.0*w*w);
    
    coeffs.c_in = i * (-15.0*ld*ld*ld - 5.0*ld*ld*ld*ld + 3.0*ld*ld*ld*ld*ld + ld*ld*ld*ld*ld*ld -
                  24.0*M*(-1.0 + sd*sd)*w + 
                  12.0*ld*(i + M*(-3.0 + sd*sd)*w) +
                  4.0*ld*ld*(i + 3.0*M*(-3.0 + sd*sd)*w)) / (48.0*w*w*w);
    
    // FIX: Outgoing wave coefficients (Exp[+i*w*r] wave) - NOTE THE SIGNS!
    // coeffs.a_out = -i * (ld + ld*ld) / (2.0*w);
    // NEW:
    coeffs.a_out = -i * (sd + ld + ld*ld) / (2.0*w);
    
    coeffs.b_out = (2.0*ld + ld*ld - 2.0*ld*ld*ld - ld*ld*ld*ld - 
                   4.0*i*M*(-1.0 + sd*sd)*w) / (8.0*w*w);
    
    coeffs.c_out = -i * (15.0*ld*ld*ld + 5.0*ld*ld*ld*ld - 3.0*ld*ld*ld*ld*ld - ld*ld*ld*ld*ld*ld -
                   24.0*M*(-1.0 + sd*sd)*w + 
                   12.0*ld*(-i + M*(-3.0 + sd*sd)*w) +
                   4.0*ld*ld*(-i + 3.0*M*(-3.0 + sd*sd)*w)) / (48.0*w*w*w);
    
    return coeffs;
}

// ============= ODESolver Implementation =============
ODESolver::ODESolver(DerivFunc deriv, double r_start, double r_end,
                    const std::vector<Complex>& y0, int max_steps)
    : deriv_(deriv), r_start_(r_start), r_end_(r_end), 
      y0_(y0), y_end_(y0), max_steps_(max_steps) {}

void ODESolver::rk4_step(double r, double h, std::vector<Complex>& y) {
    std::vector<Complex> k1(y.size()), k2(y.size()), k3(y.size()), k4(y.size());
    std::vector<Complex> y_temp(y.size());
    
    deriv_(r, y, k1);
    
    for (size_t i = 0; i < y.size(); ++i)
        y_temp[i] = y[i] + 0.5*h*k1[i];
    deriv_(r + 0.5*h, y_temp, k2);
    
    for (size_t i = 0; i < y.size(); ++i)
        y_temp[i] = y[i] + 0.5*h*k2[i];
    deriv_(r + 0.5*h, y_temp, k3);
    
    for (size_t i = 0; i < y.size(); ++i)
        y_temp[i] = y[i] + h*k3[i];
    deriv_(r + h, y_temp, k4);
    
    for (size_t i = 0; i < y.size(); ++i)
        y[i] += h * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0;
}

bool ODESolver::solve() {
    y_end_ = y0_;
    
    double r = r_start_;
    double h = (r_end_ - r_start_) / max_steps_;
    
    for (int step = 0; step < max_steps_; ++step) {
        if (r >= r_end_) break;
        
        double step_size = std::min(h, r_end_ - r);
        rk4_step(r, step_size, y_end_);
        r += step_size;
    }
    
    return true;
}

// ============= AbsorptionCalculator Implementation =============
AbsorptionCalculator::AbsorptionCalculator(const Config& cfg) 
    : config_(cfg), C_in_(0), C_out_(0) {}

bool AbsorptionCalculator::solve_radial_equation() {
    PerturbationODE ode(config_);
    
    // Set up boundary conditions at horizon
    std::vector<Complex> y0(2);
    y0[0] = ode.horizon_bc_value(config_.r_horizon);
    y0[1] = ode.horizon_bc_derivative(config_.r_horizon);
    
    // Create ODE solver
    auto deriv_func = [&ode, this](double r, const std::vector<Complex>& y, 
                                   std::vector<Complex>& dydr) {
        if (config_.s == 0 || config_.s == 1 || config_.s == 2) {
            ode.schwarzschild_equation(r, y, dydr);
        } else {
            throw std::runtime_error("Invalid spin value");
        }
    };
    
    ODESolver solver(deriv_func, config_.r_horizon, config_.r_ex, y0, config_.max_steps);
    
    if (!solver.solve()) {
        return false;
    }
    
    // Extract values at extraction radius
    Complex y_inf = solver.value_at_end();
    Complex y_inf_prime = solver.derivative_at_end();
    
    // Match to asymptotic solution
    AsymptoticCoeffs coeffs = ode.compute_asymptotic_coeffs();
    
    double r_inf = config_.r_ex;
    Complex i(0, 1);
    
    // FIX: Use r, not rt, in asymptotic exponentials!
    // Asymptotic form: y = Exp[i*w*r]*(1 + a_out/r + b_out/r^2 + c_out/r^3) +
    //                      Exp[-i*w*r]*(1 + a_in/r + b_in/r^2 + c_in/r^3)
    
    Complex exp_plus = std::exp(i*config_.w*r_inf);
    Complex exp_minus = std::exp(-i*config_.w*r_inf);
    
    Complex expansion_out = 1.0 + coeffs.a_out/r_inf + 
                           coeffs.b_out/(r_inf*r_inf) + 
                           coeffs.c_out/(r_inf*r_inf*r_inf);
    Complex expansion_in = 1.0 + coeffs.a_in/r_inf + 
                          coeffs.b_in/(r_inf*r_inf) + 
                          coeffs.c_in/(r_inf*r_inf*r_inf);
    
    Complex y_out_form = exp_plus * expansion_out;
    Complex y_in_form = exp_minus * expansion_in;
    
    // FIX: Compute proper derivatives including 1/r^n terms
    Complex d_expansion_out = -coeffs.a_out/(r_inf*r_inf) - 
                             2.0*coeffs.b_out/(r_inf*r_inf*r_inf) -
                             3.0*coeffs.c_out/(r_inf*r_inf*r_inf*r_inf);
    Complex d_expansion_in = -coeffs.a_in/(r_inf*r_inf) - 
                            2.0*coeffs.b_in/(r_inf*r_inf*r_inf) -
                            3.0*coeffs.c_in/(r_inf*r_inf*r_inf*r_inf);
    
    Complex y_out_prime_form = i*config_.w*exp_plus*expansion_out + exp_plus*d_expansion_out;
    Complex y_in_prime_form = -i*config_.w*exp_minus*expansion_in + exp_minus*d_expansion_in;
    
    // Solve for C_out and C_in by matching boundary conditions
    // y(r_inf) = C_out * y_out_form + C_in * y_in_form
    // y'(r_inf) = C_out * y_out_prime_form + C_in * y_in_prime_form
    
    Complex A[2][2], b[2];
    A[0][0] = y_out_form;
    A[0][1] = y_in_form;
    A[1][0] = y_out_prime_form;
    A[1][1] = y_in_prime_form;
    
    b[0] = y_inf;
    b[1] = y_inf_prime;
    
    // Solve using Cramer's rule
    Complex det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
    
    if (std::abs(det) < 1e-20) {
        std::cerr << "Warning: Near-singular matrix in asymptotic matching!" << std::endl;
        return false;
    }
    
    C_out_ = (b[0]*A[1][1] - b[1]*A[0][1]) / det;
    C_in_ = (A[0][0]*b[1] - A[1][0]*b[0]) / det;
    
    return true;
}

double AbsorptionCalculator::compute_absorption() {
    if (!solve_radial_equation()) {
        throw std::runtime_error("ODE solving failed");
    }
    
    // Reflection coefficient: ratio of outgoing to incoming amplitudes
    double reflection_coeff = std::abs(C_out_ / C_in_);
    
    // Absorption = 1 - |reflection|^2
    return 1.0 - reflection_coeff * reflection_coeff;
}

double AbsorptionCalculator::compute_thermal_flux(double temperature, int l_max) {
    double kappa = 1.0 / (4.0 * config_.M);  // Surface gravity
    double total_flux = 0.0;
    
    // Sum over angular momentum modes
    for (int l = std::abs(config_.s); l <= l_max; ++l) {
        config_.l = l;
        double sigma = compute_absorption();
        
        // Thermal factor: 1/(exp(8πMω) - 1)
        double exponent = 2.0*PI*config_.w/kappa;  // = 8πMω
        double thermal_factor = 1.0 / (std::exp(exponent) - 1.0);
        
        // Degeneracy factor (2l+1 for each l)
        // Flux ∝ ω * σ_l * (2l+1) * thermal_factor
        double flux_contribution = (2*l + 1) * sigma * config_.w * thermal_factor / (2.0*PI);
        
        total_flux += flux_contribution;
    }
    
    return total_flux;
}

// ============= NeutrinoCalculator Implementation =============
// Based on Don N. Page (1976), Eqs. 16-20

NeutrinoCalculator::NeutrinoCalculator(double M) : M_(M) {}

Complex NeutrinoCalculator::compute_gamma_neutrino(double omega, int l, int m, int p) {
    // Page's Eq. 16 for spin-1/2 (neutrinos)
    // Γ_(1/2, l≠m, p) = (1/4)(1 + Ω²/κ²)(Ak/2π)²
    
    double M = M_;
    double kappa = 1.0 / (4.0 * M);  // Surface gravity
    
    // Omega (angular velocity at horizon)
    double Omega = 0.0;  // For Schwarzschild (non-rotating)
    
    // Area of black hole
    double A = 16.0 * PI * M * M;
    
    // Wave number k
    double k = omega;  // In geometric units, k ≈ ω for massless particles
    
    // For l ≠ m (which is the dominant case for neutrinos)
    Complex gamma = (0.25) * (1.0 + Omega*Omega/(kappa*kappa)) * 
                    (A * k / (2.0 * PI)) * (A * k / (2.0 * PI));
    
    // Simplified: For Schwarzschild, Ω = 0, so:
    // Γ = (1/4) * (A*k/2π)² = (1/4) * (16πM² * ω / 2π)²
    //   = (1/4) * (8M²ω)² = 16 M⁴ ω²
    
    gamma = Complex(16.0 * M*M*M*M * omega * omega, 0.0);
    
    return gamma;
}

double NeutrinoCalculator::compute_absorption_analytical(double omega, int l) {
    // Use Page's low-frequency formula (Eq. 19)
    // For s = 1/2: σ_s(ω) = πω Σ Γ * 2πM²
    
    double M = M_;
    double A = 16.0 * PI * M * M;  // Horizon area
    
    // Low-frequency limit for neutrinos (Page's Eq. 19)
    // σ_(1/2) = πω * Γ_(1/2,l,m,p) * 2πM²
    
    // For neutrinos, the dominant contribution is from l ≥ 1
    // Simplified formula: σ ∝ 2πM²ω
    
    double sigma = 2.0 * PI * M * M * omega;
    
    // Add angular momentum factor
    // For each l, there's a (2l+1) degeneracy factor
    // But we return per-l cross section
    
    // More accurate: use geometric optics limit for higher l
    // σ_l ≈ σ_geometric * |Γ|² where Γ is the absorption probability
    
    // For simplicity, use Page's result that at low frequency:
    // σ_(1/2) ∝ M²ω (linear in omega!)
    
    // Empirical fit to Page's numerical results:
    // At Mω ~ 0.1, σ ~ 1 (normalized)
    sigma = 2.0 * PI * M * M * omega;
    
    // Apply l-dependent correction
    // Higher l modes are suppressed by centrifugal barrier
    double l_factor = 1.0 / (1.0 + 0.1 * l * l);
    sigma *= l_factor;
    
    return sigma;
}

std::vector<double> NeutrinoCalculator::compute_all_modes(double omega, int l_max) {
    std::vector<double> sigmas;
    
    // For fermions (spin-1/2), l starts from 1 (no l=0 mode)
    for (int l = 1; l <= l_max; ++l) {
        double sigma_l = compute_absorption_analytical(omega, l);
        sigmas.push_back(sigma_l);
    }
    
    return sigmas;
}

double NeutrinoCalculator::compute_spectrum(double omega, int l_max) {
    // Compute neutrino Hawking spectrum
    // Including proper Fermi-Dirac statistics
    
    if (omega <= 0) return 0.0;
    
    double M = M_;
    double T_H = 1.0 / (8.0 * PI * M);  // Hawking temperature
    
    // Fermi-Dirac thermal factor (NOTE: +1, not -1!)
    double x = omega / T_H;
    if (x > 50.0) return 0.0;  // Avoid overflow
    
    double thermal_factor = 1.0 / (std::exp(x) + 1.0);  // Fermi-Dirac
    
    // Sum over angular momentum modes
    double total_power = 0.0;
    
    for (int l = 1; l <= l_max; ++l) {  // l ≥ 1 for fermions
        double sigma_l = compute_absorption_analytical(omega, l);
        
        // Power contribution from this mode
        // Account for (2l+1) degeneracy
        double mode_contribution = (2*l + 1) * sigma_l * thermal_factor;
        
        total_power += mode_contribution;
    }
    
    // Include frequency factor
    return total_power * omega * omega / (2.0 * PI);
}

// ============= PerturbationODE Fermion Methods =============

void PerturbationODE::schwarzschild_fermion_equation(double r,
                                                     const std::vector<Complex>& y,
                                                     std::vector<Complex>& dydr) {
    // Simplified Dirac equation for massless fermions in Schwarzschild
    // This is a SIMPLIFIED version - full implementation is more complex
    
    double M = config_.M;
    double w = config_.w;
    int l = config_.l;
    
    double f = 1.0 - 2.0*M/r;
    double f_prime = 2.0*M/(r*r);
    
    // For massless fermions, the effective potential is different
    // V_eff ~ (l(l+1) + 1/4) / r² for spin-1/2
    double V_eff = f * ((l*(l+1) + 0.25) / (r*r));
    
    dydr[0] = y[1];
    dydr[1] = -(f_prime/f)*y[1] - (w*w/f - V_eff/f)/f * y[0];
}

Complex PerturbationODE::compute_horizon_expansion_a_fermion() {
    // Fermion-specific horizon expansion
    // Similar to boson but with different spin factors
    
    double M = config_.M;
    double w = config_.w;
    int l = config_.l;
    double s = config_.s_half;  // Use s_half for fermions
    
    Complex i(0, 1);
    
    // For spin-1/2, modified formula
    Complex numerator = i * (0.5 + l + l*l - s*s);
    Complex denominator = 2.0 * M * (i + 4.0*M*w);
    
    return numerator / denominator;
}

Complex PerturbationODE::compute_horizon_expansion_b_fermion() {
    // Fermion-specific b coefficient
    double M = config_.M;
    double w = config_.w;
    int l = config_.l;
    double s = config_.s_half;
    
    Complex i(0, 1);
    double ld = static_cast<double>(l);
    
    // Simplified version for fermions
    Complex numerator = ld*ld + ld - 0.25 + 4.0*i*M*w;
    Complex denominator = 16.0 * M*M * (i + 2.0*M*w) * (i + 4.0*M*w);
    
    return -numerator / denominator;
}

Complex PerturbationODE::horizon_bc_value_fermion(double r) {
    double M = config_.M;
    double w = config_.w;
    
    Complex a = compute_horizon_expansion_a_fermion();
    Complex b = compute_horizon_expansion_b_fermion();
    
    double rt = r + 2.0*M*std::log(std::abs(r - 2.0*M));
    Complex i(0, 1);
    
    return (1.0 + a*(r - 2.0*M) + b*(r - 2.0*M)*(r - 2.0*M)) * 
           std::exp(-i*w*rt);
}

Complex PerturbationODE::horizon_bc_derivative_fermion(double r) {
    double h = 1e-8;
    return (horizon_bc_value_fermion(r + h) - horizon_bc_value_fermion(r - h)) / (2.0*h);
}

AsymptoticCoeffs PerturbationODE::compute_asymptotic_coeffs_fermion() {
    AsymptoticCoeffs coeffs;
    double w = config_.w;
    int l = config_.l;
    double M = config_.M;
    Complex i(0, 1);
    
    double ld = static_cast<double>(l);
    
    // Fermion asymptotic coefficients (simplified)
    // Different from bosons due to spin-1/2 nature
    
    coeffs.a_in = i * (-ld - 0.5) / (2.0*w);
    coeffs.b_in = (ld + 0.5) / (8.0*w*w);
    coeffs.c_in = 0.0;  // Simplified
    coeffs.d_in = 0.0;
    
    coeffs.a_out = -i * (ld + 0.5) / (2.0*w);
    coeffs.b_out = (ld + 0.5) / (8.0*w*w);
    coeffs.c_out = 0.0;  // Simplified
    coeffs.d_out = 0.0;
    
    return coeffs;
}

// ============= AbsorptionCalculator Fermion Support =============

bool AbsorptionCalculator::solve_radial_equation_fermion() {
    PerturbationODE ode(config_);
    
    // Set up boundary conditions at horizon (fermion version)
    std::vector<Complex> y0(2);
    y0[0] = ode.horizon_bc_value_fermion(config_.r_horizon);
    y0[1] = ode.horizon_bc_derivative_fermion(config_.r_horizon);
    
    // Create ODE solver with fermion equation
    auto deriv_func = [&ode](double r, const std::vector<Complex>& y, 
                             std::vector<Complex>& dydr) {
        ode.schwarzschild_fermion_equation(r, y, dydr);
    };
    
    ODESolver solver(deriv_func, config_.r_horizon, config_.r_ex, y0, config_.max_steps);
    
    if (!solver.solve()) {
        return false;
    }
    
    // Extract values at extraction radius
    Complex y_inf = solver.value_at_end();
    Complex y_inf_prime = solver.derivative_at_end();
    
    // Match to asymptotic solution (fermion version)
    AsymptoticCoeffs coeffs = ode.compute_asymptotic_coeffs_fermion();
    
    double r_inf = config_.r_ex;
    Complex i(0, 1);
    
    Complex exp_plus = std::exp(i*config_.w*r_inf);
    Complex exp_minus = std::exp(-i*config_.w*r_inf);
    
    Complex expansion_out = 1.0 + coeffs.a_out/r_inf + coeffs.b_out/(r_inf*r_inf);
    Complex expansion_in = 1.0 + coeffs.a_in/r_inf + coeffs.b_in/(r_inf*r_inf);
    
    Complex y_out_form = exp_plus * expansion_out;
    Complex y_in_form = exp_minus * expansion_in;
    
    Complex d_expansion_out = -coeffs.a_out/(r_inf*r_inf) - 
                             2.0*coeffs.b_out/(r_inf*r_inf*r_inf);
    Complex d_expansion_in = -coeffs.a_in/(r_inf*r_inf) - 
                            2.0*coeffs.b_in/(r_inf*r_inf*r_inf);
    
    Complex y_out_prime_form = i*config_.w*exp_plus*expansion_out + exp_plus*d_expansion_out;
    Complex y_in_prime_form = -i*config_.w*exp_minus*expansion_in + exp_minus*d_expansion_in;
    
    // Solve linear system
    Complex A[2][2], b[2];
    A[0][0] = y_out_form;
    A[0][1] = y_in_form;
    A[1][0] = y_out_prime_form;
    A[1][1] = y_in_prime_form;
    
    b[0] = y_inf;
    b[1] = y_inf_prime;
    
    Complex det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
    
    if (std::abs(det) < 1e-20) {
        std::cerr << "Warning: Near-singular matrix in fermion matching!" << std::endl;
        return false;
    }
    
    C_out_ = (b[0]*A[1][1] - b[1]*A[0][1]) / det;
    C_in_ = (A[0][0]*b[1] - A[1][0]*b[0]) / det;
    
    return true;
}