#include "black_hole.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

const double PI = 3.14159265358979323846;

// Helper function to compute proper extraction radius
double get_extraction_radius(double w, double M) {
    double r_min = 100.0 * M;
    double r_max = 1000.0 * M;
    double r_wave = 50.0 / std::abs(w);
    return std::max(r_min, std::min(r_max, r_wave));
}

// Compute neutrino spectrum using Page's analytical formula
void compute_neutrino_spectrum_analytical() {
    std::cout << "==============================================\n";
    std::cout << "Computing Neutrino (spin-1/2) Spectrum\n";
    std::cout << "Using Page (1976) Analytical Formulas\n";
    std::cout << "==============================================\n\n";
    
    double M = 1.0;
    int l_max = 20;  // Sum l=1,2,...,20 (fermions: l_min=1, no l=0 mode!)
    
    NeutrinoCalculator neutrino_calc(M);
    
    std::ofstream out("neutrino_hawking_spectrum.dat");
    out << "# omega    hawking_power    temperature=" << 1.0/(8.0*PI*M) << "\n";
    out << std::setprecision(12);
    
    // Frequency sampling (same as other spectra)
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
            double power = neutrino_calc.compute_spectrum(w, l_max);
            out << w << "    " << power << "\n";
            
            std::cout << "  Total power: " << power << "\n\n";
        } catch (const std::exception& e) {
            std::cerr << "Error at omega = " << w << ": " << e.what() << "\n\n";
        }
    }
    
    out.close();
    std::cout << "Results saved to neutrino_hawking_spectrum.dat\n\n";
    std::cout << "Note: Using Page (1976) analytical formulas\n";
    std::cout << "      Fermi-Dirac statistics applied (1/(exp(x)+1))\n\n";
}

// Original boson spectrum functions
double compute_single_mode_absorption(double w, int l, int s, double M) {
    Config cfg;
    cfg.M = M;
    cfg.w = w;
    cfg.l = l;
    cfg.s = s;
    cfg.type = ParticleType::BOSON;
    cfg.r_horizon = 2.0 * M * (1.0 + 1e-4);
    cfg.r_ex = get_extraction_radius(w, M);
    cfg.max_steps = 1000000;
    
    AbsorptionCalculator calc(cfg);
    return calc.compute_absorption();
}

double compute_hawking_spectrum(double w, int s, double M, int l_max) {
    if (w <= 0) return 0.0;
    
    double T_H = 1.0 / (8.0 * PI * M);
    double x = w / T_H;
    
    if (x > 50.0) return 0.0;
    
    // BOSE-EINSTEIN for bosons
    double thermal_factor = 1.0 / (std::exp(x) - 1.0);
    
    double total_power = 0.0;
    // CRITICAL: l_min = |s| ensures correct angular momentum selection rules
    // Scalar (s=0): l = 0,1,2,... | Photon (s=1): l = 1,2,3,... | Graviton (s=2): l = 2,3,4,...
    int l_min = std::abs(s);
    
    for (int l = l_min; l <= l_max; ++l) {
        try {
            double sigma_l = compute_single_mode_absorption(w, l, s, M);
            double mode_contribution = (2*l + 1) * sigma_l * thermal_factor;
            total_power += mode_contribution;
            
            std::cout << "    l=" << l << ": sigma=" << std::setw(10) << sigma_l 
                     << ", contribution=" << std::setw(12) << mode_contribution << "\n";
        } catch (const std::exception& e) {
            std::cerr << "    Error at l=" << l << ": " << e.what() << "\n";
        }
    }
    
    return total_power * w * w / (2.0 * PI);
}

void compute_graviton_spectrum() {
    std::cout << "==============================================\n";
    std::cout << "Computing Graviton (spin-2) Hawking Spectrum\n";
    std::cout << "==============================================\n\n";
    
    double M = 1.0;
    int s = 2;
    int l_max = 20;  // Need l_max=20 to capture peaks correctly (Page 1976)
    
    std::ofstream out("graviton_hawking_spectrum.dat");
    out << "# omega    hawking_power    temperature=" << 1.0/(8.0*PI*M) << "\n";
    out << std::setprecision(12);
    
    std::vector<double> omega_values;
    for (int i = 1; i <= 20; ++i) omega_values.push_back(i * 0.01);
    for (int i = 1; i <= 20; ++i) omega_values.push_back(0.20 + i * 0.02);
    for (int i = 1; i <= 10; ++i) omega_values.push_back(0.60 + i * 0.05);
    
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
    int l_max = 20;  // Need l_max=20 to capture peaks correctly (Page 1976)
    
    std::ofstream out("photon_hawking_spectrum.dat");
    out << "# omega    hawking_power    temperature=" << 1.0/(8.0*PI*M) << "\n";
    out << std::setprecision(12);
    
    std::vector<double> omega_values;
    for (int i = 1; i <= 20; ++i) omega_values.push_back(i * 0.01);
    for (int i = 1; i <= 20; ++i) omega_values.push_back(0.20 + i * 0.02);
    for (int i = 1; i <= 10; ++i) omega_values.push_back(0.60 + i * 0.05);
    
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
    int l_max = 20;  // Need l_max=20 to capture peaks correctly (Page 1976)
    
    std::ofstream out("scalar_hawking_spectrum.dat");
    out << "# omega    hawking_power    temperature=" << 1.0/(8.0*PI*M) << "\n";
    out << std::setprecision(12);
    
    std::vector<double> omega_values;
    for (int i = 1; i <= 20; ++i) omega_values.push_back(i * 0.01);
    for (int i = 1; i <= 20; ++i) omega_values.push_back(0.20 + i * 0.02);
    for (int i = 1; i <= 10; ++i) omega_values.push_back(0.60 + i * 0.05);
    
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

void quick_test() {
    std::cout << "==============================================\n";
    std::cout << "Quick Test: Single Point Calculation\n";
    std::cout << "==============================================\n\n";
    
    double M = 1.0;
    double w = 0.1;
    
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
    std::cout << "Based on Don N. Page (1976)\n";
    std::cout << "==============================================\n\n";
    
    // Run quick test first
    quick_test();
    
    std::cout << "Starting full spectrum calculations...\n";
    std::cout << "(This may take a while)\n\n";
    
    if (argc > 1) {
        std::string arg = argv[1];
        if (arg == "scalar") {
            compute_scalar_spectrum();
        } else if (arg == "photon") {
            compute_photon_spectrum();
        } else if (arg == "graviton") {
            compute_graviton_spectrum();
        } else if (arg == "neutrino") {
            compute_neutrino_spectrum_analytical();
        } else if (arg == "all") {
            compute_scalar_spectrum();
            compute_photon_spectrum();
            compute_graviton_spectrum();
            compute_neutrino_spectrum_analytical();
        } else if (arg == "page") {
            // Page's three main particles: photons, gravitons, neutrinos
            compute_photon_spectrum();
            compute_graviton_spectrum();
            compute_neutrino_spectrum_analytical();
        } else if (arg == "quick") {
            std::cout << "Quick test only - skipping full calculations\n";
            return 0;
        } else {
            std::cout << "Usage: " << argv[0] << " [command]\n\n";
            std::cout << "Commands:\n";
            std::cout << "  scalar    - Scalar field (spin-0)\n";
            std::cout << "  photon    - Photons (spin-1)\n";
            std::cout << "  graviton  - Gravitons (spin-2)\n";
            std::cout << "  neutrino  - Neutrinos (spin-1/2, analytical)\n";
            std::cout << "  page      - Page's three: photon + graviton + neutrino\n";
            std::cout << "  all       - Everything (scalar + photon + graviton + neutrino)\n";
            std::cout << "  quick     - Quick test only\n";
            return 1;
        }
    } else {
        // Default: compute Page's main results (photons + gravitons + neutrinos)
        std::cout << "Computing Page (1976) main results...\n";
        std::cout << "  - Photons (spin-1, bosons)\n";
        std::cout << "  - Gravitons (spin-2, bosons)\n";
        std::cout << "  - Neutrinos (spin-1/2, fermions, analytical)\n\n";
        
        compute_photon_spectrum();
        compute_graviton_spectrum();
        compute_neutrino_spectrum_analytical();
    }
    
    std::cout << "\nAll computations complete!\n";
    std::cout << "\nTo plot results in Page (1976) style:\n";
    std::cout << "  python3 plot_page_1976_exact.py\n";
    
    return 0;
}