#include "black_hole.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

const double PI = 3.14159265358979323846;

// Helper function to compute proper extraction radius
double get_extraction_radius(double w, double M) {
    // Need r >> M for asymptotic expansion to be valid
    // But also need reasonable computational cost
    // Rule: r_ex = max(100M, min(1000M, c/w)) where c ~ 50-100
    double r_min = 100.0 * M;
    double r_max = 1000.0 * M;
    double r_wave = 50.0 / std::abs(w);  // ~wavelength dependent
    
    return std::max(r_min, std::min(r_max, r_wave));
}

// Compute absorption cross-section for a single (l, s, w) mode
double compute_single_mode_absorption(double w, int l, int s, double M) {
    Config cfg;
    cfg.M = M;
    cfg.w = w;
    cfg.l = l;
    cfg.s = s;
    cfg.r_horizon = 2.0 * M * (1.0 + 1e-4);
    cfg.r_ex = get_extraction_radius(w, M);
    cfg.max_steps = 1000000;
    
    AbsorptionCalculator calc(cfg);
    return calc.compute_absorption();
}

// Compute Hawking radiation spectrum at given frequency
double compute_hawking_spectrum(double w, int s, double M, int l_max) {
    if (w <= 0) return 0.0;
    
    // Hawking temperature for Schwarzschild: T_H = 1/(8πM)
    double T_H = 1.0 / (8.0 * PI * M);
    
    // Thermal factor: 1/(exp(ω/T_H) - 1)
    double x = w / T_H;  // = 8πMω
    
    // Avoid numerical overflow for large x
    if (x > 50.0) return 0.0;
    
    double thermal_factor = 1.0 / (std::exp(x) - 1.0);
    
    // Sum over angular momentum modes
    double total_power = 0.0;
    int l_min = std::abs(s);  // Minimum l for given spin
    
    for (int l = l_min; l <= l_max; ++l) {
        try {
            double sigma_l = compute_single_mode_absorption(w, l, s, M);
            
            // Power contribution from this mode:
            // dP/dω ∝ (2l+1) * σ_l * ω^2 * thermal_factor
            // (The ω^2 comes from phase space; some formulations use ω^3)
            double mode_contribution = (2*l + 1) * sigma_l * thermal_factor;
            
            total_power += mode_contribution;
            
            std::cout << "    l=" << l << ": sigma=" << std::setw(10) << sigma_l 
                     << ", contribution=" << std::setw(12) << mode_contribution << "\n";
        } catch (const std::exception& e) {
            std::cerr << "    Error at l=" << l << ": " << e.what() << "\n";
        }
    }
    
    // Include frequency factor
    return total_power * w * w / (2.0 * PI);
}

void compute_graviton_spectrum() {
    std::cout << "==============================================\n";
    std::cout << "Computing Graviton (spin-2) Hawking Spectrum\n";
    std::cout << "==============================================\n\n";
    
    double M = 1.0;
    int s = 2;
    int l_max = 5;  // Sum l=2,3,4,5 (minimum is l=2 for s=2)
    
    std::ofstream out("graviton_hawking_spectrum.dat");
    out << "# omega    hawking_power    temperature=" << 1.0/(8.0*PI*M) << "\n";
    out << std::setprecision(12);
    
    // Scan frequency range
    // Hawking temperature: T_H ~ 0.04 for M=1
    // Peak around ω ~ 2.5*T_H ~ 0.1
    
    std::vector<double> omega_values;
    
    // Low frequency region (detailed)
    for (int i = 1; i <= 20; ++i) {
        omega_values.push_back(i * 0.01);  // 0.01 to 0.20
    }
    // Medium frequency
    for (int i = 1; i <= 20; ++i) {
        omega_values.push_back(0.20 + i * 0.02);  // 0.22 to 0.60
    }
    // High frequency
    for (int i = 1; i <= 10; ++i) {
        omega_values.push_back(0.60 + i * 0.05);  // 0.65 to 1.10
    }
    
    for (double w : omega_values) {
        std::cout << "Computing omega = " << std::setw(8) << w << "\n";
        
        try {
            double power = compute_hawking_spectrum(w, s, M, l_max);
            out << w << "    " << power << "\n";
            
            std::cout << "  Total power: " << power << "\n\n";
        } catch (const std::exception& e) {
            std::cerr << "Error at omega = " << w << ": " << e.what() << "\n\n";
        }
    }
    
    out.close();
    std::cout << "Results saved to graviton_hawking_spectrum.dat\n\n";
}

void compute_photon_spectrum() {
    std::cout << "==============================================\n";
    std::cout << "Computing Photon (spin-1) Hawking Spectrum\n";
    std::cout << "==============================================\n\n";
    
    double M = 1.0;
    int s = 1;
    int l_max = 20;  // Sum l=1,2,3,4,5
    
    std::ofstream out("photon_hawking_spectrum.dat");
    out << "# omega    hawking_power    temperature=" << 1.0/(8.0*PI*M) << "\n";
    out << std::setprecision(12);
    
    std::vector<double> omega_values;
    for (int i = 1; i <= 20; ++i) {
        omega_values.push_back(i * 0.01);
    }
    for (int i = 1; i <= 20; ++i) {
        omega_values.push_back(0.20 + i * 0.02);
    }
    for (int i = 1; i <= 10; ++i) {
        omega_values.push_back(0.60 + i * 0.05);
    }
    
    for (double w : omega_values) {
        std::cout << "Computing omega = " << std::setw(8) << w << "\n";
        
        try {
            double power = compute_hawking_spectrum(w, s, M, l_max);
            out << w << "    " << power << "\n";
            
            std::cout << "  Total power: " << power << "\n\n";
        } catch (const std::exception& e) {
            std::cerr << "Error at omega = " << w << ": " << e.what() << "\n\n";
        }
    }
    
    out.close();
    std::cout << "Results saved to photon_hawking_spectrum.dat\n\n";
}

void compute_scalar_spectrum() {
    std::cout << "==============================================\n";
    std::cout << "Computing Scalar (spin-0) Hawking Spectrum\n";
    std::cout << "==============================================\n\n";
    
    double M = 1.0;
    int s = 0;
    int l_max = 5;  // Sum l=0,1,2,3,4,5 (monopole is important!)
    
    std::ofstream out("scalar_hawking_spectrum.dat");
    out << "# omega    hawking_power    temperature=" << 1.0/(8.0*PI*M) << "\n";
    out << std::setprecision(12);
    
    std::vector<double> omega_values;
    for (int i = 1; i <= 20; ++i) {
        omega_values.push_back(i * 0.01);
    }
    for (int i = 1; i <= 20; ++i) {
        omega_values.push_back(0.20 + i * 0.02);
    }
    for (int i = 1; i <= 10; ++i) {
        omega_values.push_back(0.60 + i * 0.05);
    }
    
    for (double w : omega_values) {
        std::cout << "Computing omega = " << std::setw(8) << w << "\n";
        
        try {
            double power = compute_hawking_spectrum(w, s, M, l_max);
            out << w << "    " << power << "\n";
            
            std::cout << "  Total power: " << power << "\n\n";
        } catch (const std::exception& e) {
            std::cerr << "Error at omega = " << w << ": " << e.what() << "\n\n";
        }
    }
    
    out.close();
    std::cout << "Results saved to scalar_hawking_spectrum.dat\n\n";
}

// Quick test function to verify code is working
void quick_test() {
    std::cout << "==============================================\n";
    std::cout << "Quick Test: Single Point Calculation\n";
    std::cout << "==============================================\n\n";
    
    double M = 1.0;
    double w = 0.1;  // Near peak of spectrum
    
    std::cout << "Testing scalar field at omega = " << w << "\n";
    std::cout << "Hawking temperature T_H = " << 1.0/(8.0*PI*M) << "\n\n";
    
    for (int l = 0; l <= 3; ++l) {
        try {
            double sigma = compute_single_mode_absorption(w, l, 0, M);
            std::cout << "l=" << l << ": sigma = " << sigma << "\n";
        } catch (const std::exception& e) {
            std::cerr << "l=" << l << ": ERROR - " << e.what() << "\n";
        }
    }
    std::cout << "\n";
}

int main(int argc, char* argv[]) {
    std::cout << std::setprecision(10);
    
    std::cout << "==============================================\n";
    std::cout << "Black Hole Hawking Radiation Calculator\n";
    std::cout << "==============================================\n\n";
    
    // Run quick test first
    quick_test();
    
    // If test passes, run full calculations
    std::cout << "Starting full spectrum calculations...\n";
    std::cout << "(This may take a while - each frequency requires solving ODEs for multiple l values)\n\n";
    
    // User can choose which to compute
    if (argc > 1) {
        std::string arg = argv[1];
        if (arg == "scalar") {
            compute_scalar_spectrum();
        } else if (arg == "photon") {
            compute_photon_spectrum();
        } else if (arg == "graviton") {
            compute_graviton_spectrum();
        } else if (arg == "all") {
            compute_scalar_spectrum();
            compute_photon_spectrum();
            compute_graviton_spectrum();
        } else if (arg == "quick") {
            // Only quick test, no full computation
            std::cout << "Quick test only - skipping full calculations\n";
            return 0;
        } else {
            std::cout << "Usage: " << argv[0] << " [scalar|photon|graviton|all|quick]\n";
            std::cout << "  scalar   - Compute scalar field spectrum only\n";
            std::cout << "  photon   - Compute photon spectrum only\n";
            std::cout << "  graviton - Compute graviton spectrum only\n";
            std::cout << "  all      - Compute all three spectra (default)\n";
            std::cout << "  quick    - Quick test only, no full computation\n";
            return 1;
        }
    } else {
        // Default: compute ALL THREE spectra
        std::cout << "Computing all three spectra (scalar, photon, graviton)...\n";
        std::cout << "This will take approximately 2-3 hours.\n";
        std::cout << "To compute only one type, use: " << argv[0] << " scalar|photon|graviton\n\n";
        
        compute_scalar_spectrum();
        compute_photon_spectrum();
        compute_graviton_spectrum();
    }
    
    std::cout << "All computations complete!\n";
    std::cout << "\nTo plot results, use gnuplot:\n";
    std::cout << "  gnuplot> plot 'scalar_hawking_spectrum.dat' with lines\n";
    
    return 0;
}