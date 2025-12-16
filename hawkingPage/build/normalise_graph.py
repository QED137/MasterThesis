#!/usr/bin/env python3
"""
Normalized Hawking Radiation Spectrum Plotter
Following Don Page's normalization: M × dE/(dω dt)
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Try to use nice plotting style
try:
    plt.style.use('seaborn-v0_8-darkgrid')
except:
    try:
        plt.style.use('seaborn-darkgrid')
    except:
        pass

def load_spectrum(filename):
    """Load spectrum data from file"""
    if not os.path.exists(filename):
        print(f"Error: File {filename} not found!")
        return None
    
    data = np.loadtxt(filename, comments='#')
    
    if data.size == 0:
        print(f"Error: File {filename} is empty!")
        return None
    
    return data

def normalize_page_style(omega, power, M=1.0):
    """
    Normalize spectrum in Don Page style: M × dE/(dω dt)
    
    The raw power output from our code is essentially dE/(dω dt).
    To match Page's plots, we multiply by M.
    
    For M=1: this is just the power itself
    For other M: scale by M
    """
    # The normalization is: M × (power spectrum)
    # This makes the y-axis dimensionless and independent of M
    normalized_power = M * power
    
    return normalized_power

def plot_page_style_normalized():
    """Plot all spectra with Don Page normalization"""
    
    files = {
        'scalar_hawking_spectrum.dat': ('Scalar (s=0)', 'blue', '-'),
        'photon_hawking_spectrum.dat': ('Photons (s=1)', 'green', '-'),
        'graviton_hawking_spectrum.dat': ('Gravitons (s=2)', 'red', '-')
    }
    
    # Check which files exist
    available_files = {f: info for f, info in files.items() if os.path.exists(f)}
    
    if not available_files:
        print("Error: No spectrum files found!")
        print("Expected: scalar_hawking_spectrum.dat, photon_hawking_spectrum.dat, graviton_hawking_spectrum.dat")
        return
    
    # Physical parameters
    M = 1.0  # Black hole mass
    T_H = 1.0 / (8.0 * np.pi * M)  # Hawking temperature
    
    # Create figure matching Page's style
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Track total (sum of all species)
    omega_total = None
    power_total = None
    
    # Plot each spectrum
    for filename, (label, color, linestyle) in available_files.items():
        data = load_spectrum(filename)
        if data is not None:
            omega = data[:, 0]
            power_raw = data[:, 1]
            
            # Apply Page normalization: M × dE/(dω dt)
            power_normalized = normalize_page_style(omega, power_raw, M)
            
            # Plot
            ax.plot(omega, power_normalized, color=color, linewidth=2.5,
                   linestyle=linestyle, label=label)
            
            # Accumulate for total
            if omega_total is None:
                omega_total = omega.copy()
                power_total = power_normalized.copy()
            else:
                # Interpolate to common omega grid if needed
                if len(omega) == len(omega_total) and np.allclose(omega, omega_total):
                    power_total += power_normalized
    
    # Plot total if we have multiple species
    if len(available_files) > 1 and power_total is not None:
        ax.plot(omega_total, power_total, 'black', linewidth=3,
               linestyle='-', label='Total', zorder=10)
    
    # Also add "Thermal Radiation" (sum of all) with dashed line for comparison
    # This would be the theoretical perfect blackbody
    if omega_total is not None:
        # Theoretical thermal: proportional to omega^2 / (exp(omega/T_H) - 1)
        thermal_factor = omega_total**2 / (np.exp(omega_total / T_H) - 1.0)
        # Normalize to match peak height of actual total
        if np.max(thermal_factor) > 0:
            thermal_normalized = thermal_factor * (np.max(power_total) / np.max(thermal_factor))
            ax.plot(omega_total, thermal_normalized, 'gray', linewidth=2,
                   linestyle='--', label='Thermal (ideal)', alpha=0.6)
    
    # Labels matching Page's style
    ax.set_xlabel(r'$M\omega$', fontsize=16)
    ax.set_ylabel(r'$M \times \frac{dE}{d\omega \, dt}$', fontsize=16)
    ax.set_title('Hawking Radiation Spectrum (Page Normalization)', 
                 fontsize=18, fontweight='bold')
    
    # Legend
    ax.legend(fontsize=13, loc='upper right', framealpha=0.95)
    
    # Grid
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Set reasonable axis limits
    ax.set_xlim(0, omega_total[-1])
    ax.set_ylim(0, None)
    
    plt.tight_layout()
    
    # Save
    plt.savefig('hawking_spectrum_page_normalized.png', dpi=300, bbox_inches='tight')
    print("Saved normalized plot to: hawking_spectrum_page_normalized.png")
    
    plt.show()

def plot_page_style_with_neutrinos():
    """
    Plot similar to Page's figure, including a neutrino curve
    Note: Our code doesn't compute neutrinos (fermions), but we can show
    where they would go for comparison
    """
    
    files = {
        'scalar_hawking_spectrum.dat': ('Scalar (s=0)', 'blue', '-'),
        'photon_hawking_spectrum.dat': ('Photons (s=1)', 'green', '-'),
        'graviton_hawking_spectrum.dat': ('Gravitons (s=2)', 'red', '-')
    }
    
    available_files = {f: info for f, info in files.items() if os.path.exists(f)}
    
    if not available_files:
        print("Error: No spectrum files found!")
        return
    
    M = 1.0
    T_H = 1.0 / (8.0 * np.pi * M)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    omega_total = None
    power_total = None
    
    # Plot computed spectra
    for filename, (label, color, linestyle) in available_files.items():
        data = load_spectrum(filename)
        if data is not None:
            omega = data[:, 0]
            power_raw = data[:, 1]
            power_normalized = normalize_page_style(omega, power_raw, M)
            
            ax.plot(omega, power_normalized, color=color, linewidth=2.5,
                   linestyle=linestyle, label=label)
            
            if omega_total is None:
                omega_total = omega.copy()
                power_total = power_normalized.copy()
            else:
                if len(omega) == len(omega_total) and np.allclose(omega, omega_total):
                    power_total += power_normalized
    
    # Add approximate neutrino curve (fermions would be about 7/8 of photons)
    if 'photon_hawking_spectrum.dat' in available_files:
        data = load_spectrum('photon_hawking_spectrum.dat')
        if data is not None:
            omega = data[:, 0]
            power_raw = data[:, 1]
            power_normalized = normalize_page_style(omega, power_raw, M)
            
            # Approximate neutrino contribution
            # For massless fermions: ~(7/8) × (number of species) × photon spectrum
            # Assuming 3 neutrino types (but simplified)
            neutrino_approx = 0.7 * power_normalized
            
            ax.plot(omega, neutrino_approx, 'purple', linewidth=2.5,
                   linestyle='--', label='Neutrinos (approx.)', alpha=0.7)
            
            if omega_total is not None and len(omega) == len(omega_total):
                power_total += neutrino_approx
    
    # Plot total
    if power_total is not None:
        ax.plot(omega_total, power_total, 'black', linewidth=3,
               linestyle='-', label='Total', zorder=10)
    
    # Labels
    ax.set_xlabel(r'$M\omega$', fontsize=16)
    ax.set_ylabel(r'$M \times \frac{dE}{d\omega \, dt}$', fontsize=16)
    ax.set_title('Hawking Radiation Spectrum (Don Page Style)', 
                 fontsize=18, fontweight='bold')
    
    ax.legend(fontsize=13, loc='upper right', framealpha=0.95)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_xlim(0, omega_total[-1])
    ax.set_ylim(0, None)
    
    plt.tight_layout()
    plt.savefig('hawking_spectrum_page_with_neutrinos.png', dpi=300, bbox_inches='tight')
    print("Saved Page-style plot to: hawking_spectrum_page_with_neutrinos.png")
    
    plt.show()

def plot_comparison_raw_vs_normalized():
    """Compare raw output vs Page-normalized"""
    
    filename = 'scalar_hawking_spectrum.dat'
    data = load_spectrum(filename)
    
    if data is None:
        print(f"Need {filename} for comparison plot")
        return
    
    omega = data[:, 0]
    power_raw = data[:, 1]
    
    M = 1.0
    power_normalized = normalize_page_style(omega, power_raw, M)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Raw plot
    ax1.plot(omega, power_raw, 'b-', linewidth=2.5)
    ax1.set_xlabel(r'$\omega$', fontsize=14)
    ax1.set_ylabel(r'$\frac{dE}{d\omega \, dt}$ (raw)', fontsize=14)
    ax1.set_title('Raw Spectrum', fontsize=16)
    ax1.grid(True, alpha=0.3)
    
    # Normalized plot
    ax2.plot(omega, power_normalized, 'r-', linewidth=2.5)
    ax2.set_xlabel(r'$M\omega$', fontsize=14)
    ax2.set_ylabel(r'$M \times \frac{dE}{d\omega \, dt}$', fontsize=14)
    ax2.set_title('Page-Normalized Spectrum', fontsize=16)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('comparison_raw_vs_normalized.png', dpi=300, bbox_inches='tight')
    print("Saved comparison to: comparison_raw_vs_normalized.png")
    
    plt.show()

def print_page_statistics():
    """Print statistics in Page normalization"""
    
    files = ['scalar_hawking_spectrum.dat', 
             'photon_hawking_spectrum.dat',
             'graviton_hawking_spectrum.dat']
    
    M = 1.0
    T_H = 1.0 / (8.0 * np.pi * M)
    
    print("\n" + "="*70)
    print("Statistics in Don Page Normalization")
    print("="*70)
    print(f"Black hole mass M = {M}")
    print(f"Hawking temperature T_H = {T_H:.6f}")
    print(f"Dimensionless temperature: M × T_H = {M * T_H:.6f}")
    print("="*70)
    
    for filename in files:
        if not os.path.exists(filename):
            continue
            
        data = load_spectrum(filename)
        if data is None:
            continue
        
        omega = data[:, 0]
        power_raw = data[:, 1]
        power_normalized = normalize_page_style(omega, power_raw, M)
        
        # Find peak
        peak_idx = np.argmax(power_normalized)
        peak_omega = omega[peak_idx]
        peak_power = power_normalized[peak_idx]
        
        # Dimensionless peak frequency
        Mw_peak = M * peak_omega
        
        # Integrated power
        total_power = np.trapz(power_normalized, omega)
        
        print(f"\n{filename}:")
        print(f"  Peak frequency (Mω): {Mw_peak:.4f}")
        print(f"  Peak power (M×dE/dωdt): {peak_power:.6f}")
        print(f"  Integrated total (M×dE/dt): {total_power:.6f}")
        print(f"  Peak/T_H ratio: {peak_omega/T_H:.3f}")
    
    print("="*70 + "\n")

def main():
    """Main function"""
    print("\n" + "="*70)
    print("Hawking Radiation Plotter - Don Page Normalization")
    print("="*70 + "\n")
    
    if len(sys.argv) > 1:
        command = sys.argv[1].lower()
        
        if command == 'page':
            plot_page_style_normalized()
            print_page_statistics()
        
        elif command == 'neutrinos':
            plot_page_style_with_neutrinos()
        
        elif command == 'compare':
            plot_comparison_raw_vs_normalized()
        
        elif command == 'stats':
            print_page_statistics()
        
        elif command == 'help':
            print("Usage: python3 plot_page_normalized.py [command]")
            print("\nCommands:")
            print("  page       - Plot with Page normalization (default)")
            print("  neutrinos  - Include approximate neutrino curve")
            print("  compare    - Compare raw vs normalized")
            print("  stats      - Print statistics only")
            print("  help       - Show this help")
        
        else:
            print(f"Unknown command: {command}")
            print("Use 'python3 plot_page_normalized.py help' for usage")
    
    else:
        # Default: Page-style normalized plot
        plot_page_style_normalized()
        print_page_statistics()

if __name__ == '__main__':
    main()