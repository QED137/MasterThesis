#!/usr/bin/env python3
"""
Plotting script for Hawking Radiation Spectrum
Plots the output from the C++ Hawking radiation calculator
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Check if matplotlib has nice style
try:
    plt.style.use('seaborn-v0_8-darkgrid')
except:
    try:
        plt.style.use('seaborn-darkgrid')
    except:
        pass  # Use default style

def load_spectrum(filename):
    """Load spectrum data from file"""
    if not os.path.exists(filename):
        print(f"Error: File {filename} not found!")
        return None
    
    # Load data, skipping comment lines
    data = np.loadtxt(filename, comments='#')
    
    if data.size == 0:
        print(f"Error: File {filename} is empty!")
        return None
    
    return data

def plot_single_spectrum(filename, title="Hawking Radiation Spectrum"):
    """Plot a single spectrum"""
    data = load_spectrum(filename)
    if data is None:
        return
    
    omega = data[:, 0]
    power = data[:, 1]
    
    # Calculate Hawking temperature for M=1
    M = 1.0
    T_H = 1.0 / (8.0 * np.pi * M)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot spectrum
    ax.plot(omega, power, 'b-', linewidth=2.5, label='Hawking Power')
    
    # Mark the Hawking temperature
    ax.axvline(T_H, color='red', linestyle='--', linewidth=1.5, 
               alpha=0.7, label=f'$T_H = {T_H:.4f}$')
    
    # Mark peak location
    peak_idx = np.argmax(power)
    peak_omega = omega[peak_idx]
    peak_power = power[peak_idx]
    ax.plot(peak_omega, peak_power, 'ro', markersize=10, 
            label=f'Peak: ω = {peak_omega:.3f}')
    
    # Labels and title
    ax.set_xlabel('Frequency ω', fontsize=14)
    ax.set_ylabel('Power', fontsize=14)
    ax.set_title(title + f' (M = {M})', fontsize=16, fontweight='bold')
    ax.legend(fontsize=12, loc='best')
    ax.grid(True, alpha=0.3)
    
    # Make it look nice
    plt.tight_layout()
    
    # Save figure
    output_name = filename.replace('.dat', '.png')
    plt.savefig(output_name, dpi=300, bbox_inches='tight')
    print(f"Saved plot to {output_name}")
    
    plt.show()

def plot_all_spectra():
    """Plot all three spectra on the same graph"""
    # File names and labels
    files = {
        'scalar_hawking_spectrum.dat': ('Scalar (s=0)', 'blue', 'o'),
        'photon_hawking_spectrum.dat': ('Photon (s=1)', 'green', 's'),
        'graviton_hawking_spectrum.dat': ('Graviton (s=2)', 'red', '^')
    }
    
    # Check which files exist
    available_files = {f: info for f, info in files.items() if os.path.exists(f)}
    
    if not available_files:
        print("Error: No spectrum files found!")
        print("Expected files: scalar_hawking_spectrum.dat, photon_hawking_spectrum.dat, graviton_hawking_spectrum.dat")
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 7))
    
    M = 1.0
    T_H = 1.0 / (8.0 * np.pi * M)
    
    # Plot each spectrum
    for filename, (label, color, marker) in available_files.items():
        data = load_spectrum(filename)
        if data is not None:
            omega = data[:, 0]
            power = data[:, 1]
            
            # Plot with both line and markers
            ax.plot(omega, power, color=color, linewidth=2.5, 
                   label=label, marker=marker, markevery=5, markersize=6)
    
    # Mark Hawking temperature
    ax.axvline(T_H, color='black', linestyle='--', linewidth=1.5, 
               alpha=0.5, label=f'$T_H = {T_H:.4f}$')
    
    # Labels and title
    ax.set_xlabel('Frequency ω', fontsize=14)
    ax.set_ylabel('Power', fontsize=14)
    ax.set_title('Hawking Radiation Spectrum - All Spins (M=1)', 
                 fontsize=16, fontweight='bold')
    ax.legend(fontsize=12, loc='best')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save
    plt.savefig('all_hawking_spectra.png', dpi=300, bbox_inches='tight')
    print("Saved combined plot to all_hawking_spectra.png")
    
    plt.show()

def plot_log_scale():
    """Plot spectrum with logarithmic y-axis to see exponential decay"""
    files = {
        'scalar_hawking_spectrum.dat': ('Scalar (s=0)', 'blue'),
        'photon_hawking_spectrum.dat': ('Photon (s=1)', 'green'),
        'graviton_hawking_spectrum.dat': ('Graviton (s=2)', 'red')
    }
    
    available_files = {f: info for f, info in files.items() if os.path.exists(f)}
    
    if not available_files:
        print("Error: No spectrum files found!")
        return
    
    fig, ax = plt.subplots(figsize=(12, 7))
    
    M = 1.0
    T_H = 1.0 / (8.0 * np.pi * M)
    
    for filename, (label, color) in available_files.items():
        data = load_spectrum(filename)
        if data is not None:
            omega = data[:, 0]
            power = data[:, 1]
            
            # Filter out zeros for log plot
            mask = power > 0
            ax.semilogy(omega[mask], power[mask], color=color, 
                       linewidth=2.5, label=label)
    
    ax.axvline(T_H, color='black', linestyle='--', linewidth=1.5, 
               alpha=0.5, label=f'$T_H = {T_H:.4f}$')
    
    ax.set_xlabel('Frequency ω', fontsize=14)
    ax.set_ylabel('Power (log scale)', fontsize=14)
    ax.set_title('Hawking Radiation Spectrum - Logarithmic Scale', 
                 fontsize=16, fontweight='bold')
    ax.legend(fontsize=12, loc='best')
    ax.grid(True, alpha=0.3, which='both')
    
    plt.tight_layout()
    plt.savefig('hawking_spectrum_log.png', dpi=300, bbox_inches='tight')
    print("Saved log-scale plot to hawking_spectrum_log.png")
    
    plt.show()

def plot_comparison_with_planck():
    """Compare Hawking spectrum with Planck distribution"""
    # Load scalar spectrum as example
    filename = 'scalar_hawking_spectrum.dat'
    data = load_spectrum(filename)
    
    if data is None:
        print(f"Error: Need {filename} to make comparison plot")
        return
    
    omega = data[:, 0]
    power = data[:, 1]
    
    M = 1.0
    T_H = 1.0 / (8.0 * np.pi * M)
    
    # Calculate theoretical Planck distribution (normalized)
    # Planck: ~ omega^3 / (exp(omega/T) - 1)
    planck = omega**3 / (np.exp(omega / T_H) - 1.0)
    
    # Normalize both to have same peak for comparison
    power_normalized = power / np.max(power)
    planck_normalized = planck / np.max(planck)
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
    
    # Plot 1: Actual spectrum
    ax1.plot(omega, power, 'b-', linewidth=2.5, label='Computed Hawking Spectrum')
    ax1.axvline(T_H, color='red', linestyle='--', alpha=0.7, label=f'$T_H = {T_H:.4f}$')
    ax1.set_xlabel('Frequency ω', fontsize=12)
    ax1.set_ylabel('Power', fontsize=12)
    ax1.set_title('Computed Hawking Radiation', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Comparison with Planck (normalized)
    ax2.plot(omega, power_normalized, 'b-', linewidth=2.5, 
            label='Hawking Spectrum (normalized)')
    ax2.plot(omega, planck_normalized, 'r--', linewidth=2, 
            label='Planck $\\omega^3/(e^{\\omega/T_H}-1)$')
    ax2.axvline(T_H, color='gray', linestyle=':', alpha=0.7)
    ax2.set_xlabel('Frequency ω', fontsize=12)
    ax2.set_ylabel('Normalized Power', fontsize=12)
    ax2.set_title('Comparison with Planck Distribution', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('hawking_planck_comparison.png', dpi=300, bbox_inches='tight')
    print("Saved comparison plot to hawking_planck_comparison.png")
    
    plt.show()

def print_statistics(filename):
    """Print statistics about the spectrum"""
    data = load_spectrum(filename)
    if data is None:
        return
    
    omega = data[:, 0]
    power = data[:, 1]
    
    M = 1.0
    T_H = 1.0 / (8.0 * np.pi * M)
    
    # Find peak
    peak_idx = np.argmax(power)
    peak_omega = omega[peak_idx]
    peak_power = power[peak_idx]
    
    # Calculate integrated power (trapezoidal rule)
    total_power = np.trapz(power, omega)
    
    print(f"\n{'='*60}")
    print(f"Statistics for {filename}")
    print(f"{'='*60}")
    print(f"Hawking Temperature (T_H):     {T_H:.6f}")
    print(f"Peak frequency (ω_peak):        {peak_omega:.6f}")
    print(f"Peak power:                     {peak_power:.6f}")
    print(f"ω_peak / T_H ratio:             {peak_omega/T_H:.3f}")
    print(f"Integrated total power:         {total_power:.6f}")
    print(f"Frequency range:                {omega[0]:.3f} - {omega[-1]:.3f}")
    print(f"Number of frequency points:     {len(omega)}")
    print(f"{'='*60}\n")

def main():
    """Main function"""
    print("\n" + "="*60)
    print("Hawking Radiation Spectrum Plotting Tool")
    print("="*60 + "\n")
    
    # Check command line arguments
    if len(sys.argv) > 1:
        command = sys.argv[1].lower()
        
        if command == 'scalar':
            plot_single_spectrum('scalar_hawking_spectrum.dat', 
                               'Scalar Field Hawking Spectrum')
            print_statistics('scalar_hawking_spectrum.dat')
        
        elif command == 'photon':
            plot_single_spectrum('photon_hawking_spectrum.dat', 
                               'Photon Hawking Spectrum')
            print_statistics('photon_hawking_spectrum.dat')
        
        elif command == 'graviton':
            plot_single_spectrum('graviton_hawking_spectrum.dat', 
                               'Graviton Hawking Spectrum')
            print_statistics('graviton_hawking_spectrum.dat')
        
        elif command == 'all':
            plot_all_spectra()
            for f in ['scalar_hawking_spectrum.dat', 
                     'photon_hawking_spectrum.dat',
                     'graviton_hawking_spectrum.dat']:
                if os.path.exists(f):
                    print_statistics(f)
        
        elif command == 'log':
            plot_log_scale()
        
        elif command == 'compare':
            plot_comparison_with_planck()
        
        elif command == 'help':
            print("Usage: python3 plot_hawking.py [command]")
            print("\nCommands:")
            print("  scalar    - Plot scalar spectrum")
            print("  photon    - Plot photon spectrum")
            print("  graviton  - Plot graviton spectrum")
            print("  all       - Plot all spectra together")
            print("  log       - Plot with logarithmic y-axis")
            print("  compare   - Compare with Planck distribution")
            print("  help      - Show this help message")
            print("\nDefault (no command): Plot all available spectra")
        
        else:
            print(f"Unknown command: {command}")
            print("Use 'python3 plot_hawking.py help' for usage")
    
    else:
        # Default: plot all available spectra
        print("Plotting all available spectra...")
        plot_all_spectra()
        
        # Print statistics for all files
        for f in ['scalar_hawking_spectrum.dat', 
                 'photon_hawking_spectrum.dat',
                 'graviton_hawking_spectrum.dat']:
            if os.path.exists(f):
                print_statistics(f)

if __name__ == '__main__':
    main()