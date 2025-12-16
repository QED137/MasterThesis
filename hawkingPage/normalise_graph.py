#!/usr/bin/env python3
"""
Hawking Radiation Spectrum Plotter - Don Page (1976) Style

Following Page's 1976 paper "Particle Emission Rates from a Black Hole"

Page's plots show:
- X-axis: Dimensionless frequency Mω (where M is black hole mass)
- Y-axis: Particle emission rate per species (dimensionless)
- Each curve represents emission for different spin species (s=0,1,2,...)

Key physics:
1. Hawking temperature: T_H = 1/(8πM) for Schwarzschild BH
2. Thermal spectrum: ∝ 1/(e^(ω/T) - 1) for bosons
3. Greybody factors suppress emission (especially at low ω)
4. Peak emission around Mω ~ 0.1-0.4 depending on spin

Note: Page's figures show that:
- Low-spin particles (s=0) have higher emission
- High-spin particles (s=2) have lower emission (stronger suppression)
- All spectra peak around similar frequencies but different amplitudes
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

def normalize_page_style(omega, power, M=1.0, T_H=None):
    """
    Normalize spectrum following Don Page (1976)
    
    Page (1976) equation (2.16): The dimensionless emission rate is
        Γ_s(ω) = (27/4π) × M² × ω² × σ_abs(ω) / (e^(ω/T_H) - 1)
    
    Our C++ code computes:
        power = [Σ_l (2l+1) σ_l / (e^(ω/T) - 1)] × ω² / (2π)
    
    To match Page's plots, we want dimensionless quantity.
    Page's y-axis is essentially: M² × (emission rate)
    
    For comparison across different M, use: power (already dimensionless for M=1)
    For strict Page normalization: M² × power
    
    Actually, looking at Page's Fig 1-4, the y-axis scale is O(0.001-0.01)
    which matches our power values O(10^-4), so we're close.
    
    The main issue is NOT normalization but the spectrum shape (peak location).
    """
    if T_H is None:
        T_H = 1.0 / (8.0 * np.pi * M)
    
    # For M=1 (our case), just return power as-is
    # This gives dimensionless emission rate matching Page's scale
    # (Our power already has the right dimensions)
    
    # Optional: Scale by M² to make truly dimensionless across different masses
    normalized_power = power  # For M=1, this is already correct
    
    return normalized_power

def plot_page_style_normalized():
    """Plot all spectra with Don Page normalization"""
    
    # Degeneracy factors for each spin type
    # spin 0: 1 DOF
    # spin 1: 2 polarizations  
    # spin 2: 2 polarizations
    # spin 1/2: 2 DOF × (7/8) = 1.75 per neutrino species
    
    files = {
        'scalar_hawking_spectrum.dat': ('Scalar (s=0)', 'blue', '-', 1),
        'photon_hawking_spectrum.dat': ('Photons (s=1)', 'green', '-', 2),
        'graviton_hawking_spectrum.dat': ('Gravitons (s=2)', 'red', '-', 2)
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
    total_degeneracy = 0  # Sum of all degeneracy factors
    
    # Plot each spectrum
    for filename, (label, color, linestyle, degeneracy) in available_files.items():
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
                total_degeneracy = degeneracy
            else:
                # Interpolate to common omega grid if needed
                if len(omega) == len(omega_total) and np.allclose(omega, omega_total):
                    power_total += power_normalized
                    total_degeneracy += degeneracy
    
    # Plot total if we have multiple species
    if len(available_files) > 1 and power_total is not None:
        ax.plot(omega_total, power_total, 'black', linewidth=3,
               linestyle='-', label='Total', zorder=10)
    
    # Add "Ideal Thermal" spectrum for comparison
    # This is the theoretical blackbody spectrum accounting for all species
    if omega_total is not None and total_degeneracy > 0:
        # Ideal blackbody: (degeneracy) × omega^2 / (exp(omega/T_H) - 1)
        # For bosons: g_i × ω² / (e^(ω/T) - 1)
        ideal_thermal = total_degeneracy * omega_total**2 / (np.exp(omega_total / T_H) - 1.0)
        
        # Apply Page normalization to ideal spectrum
        ideal_normalized = M * ideal_thermal
        
        # Scale to match actual total radiated power (accounts for greybody suppression)
        actual_total_power = np.trapz(power_total, omega_total)
        ideal_total_power = np.trapz(ideal_normalized, omega_total)
        
        if ideal_total_power > 0:
            # Renormalize ideal to match actual (this shows the greybody suppression)
            ideal_scaled = ideal_normalized * (actual_total_power / ideal_total_power)
            ax.plot(omega_total, ideal_scaled, 'gray', linewidth=2.5,
                   linestyle='--', label=f'Ideal Blackbody (g={total_degeneracy})', alpha=0.7)
    
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
    
    # Add text box with info
    if omega_total is not None:
        actual_power = np.trapz(power_total, omega_total)
        textstr = f'Total DOF: {total_degeneracy}\n'
        textstr += f'$T_H = {T_H:.5f}$\n'
        textstr += f'Total Power: {actual_power:.4f}'
        ax.text(0.02, 0.98, textstr, transform=ax.transAxes,
                fontsize=11, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
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

def plot_with_greybody_factors():
    """Plot showing the greybody suppression for each species"""
    
    files = {
        'scalar_hawking_spectrum.dat': ('Scalar (s=0)', 'blue', '-', 1),
        'photon_hawking_spectrum.dat': ('Photons (s=1)', 'green', '-', 2),
        'graviton_hawking_spectrum.dat': ('Gravitons (s=2)', 'red', '-', 2)
    }
    
    M = 1.0
    T_H = 1.0 / (8.0 * np.pi * M)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Left plot: Actual vs Ideal spectra
    # Right plot: Greybody factors
    
    for filename, (label, color, linestyle, degeneracy) in files.items():
        if not os.path.exists(filename):
            continue
            
        data = load_spectrum(filename)
        if data is None:
            continue
        
        omega = data[:, 0]
        power_raw = data[:, 1]
        power_normalized = normalize_page_style(omega, power_raw, M)
        
        # Ideal blackbody for this species
        ideal_bb = degeneracy * omega**2 / (np.exp(omega / T_H) - 1.0)
        ideal_normalized = M * ideal_bb
        
        # Greybody factor = actual / ideal
        # Avoid division by zero
        greybody = np.zeros_like(omega)
        mask = ideal_normalized > 1e-10
        greybody[mask] = power_normalized[mask] / ideal_normalized[mask]
        
        # Plot spectra
        ax1.plot(omega, power_normalized, color=color, linewidth=2.5,
                linestyle='-', label=f'{label} (actual)')
        ax1.plot(omega, ideal_normalized, color=color, linewidth=2,
                linestyle='--', alpha=0.5, label=f'{label} (ideal)')
        
        # Plot greybody factors
        ax2.plot(omega, greybody, color=color, linewidth=2.5,
                label=label)
    
    # Configure left plot
    ax1.set_xlabel(r'$M\omega$', fontsize=14)
    ax1.set_ylabel(r'$M \times \frac{dE}{d\omega \, dt}$', fontsize=14)
    ax1.set_title('Actual vs Ideal Spectra', fontsize=16, fontweight='bold')
    ax1.legend(fontsize=10, loc='upper right', ncol=2)
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.set_xlim(0, None)
    ax1.set_ylim(0, None)
    
    # Configure right plot
    ax2.set_xlabel(r'$M\omega$', fontsize=14)
    ax2.set_ylabel(r'Greybody Factor $\Gamma(\omega)$', fontsize=14)
    ax2.set_title('Greybody Suppression Factors', fontsize=16, fontweight='bold')
    ax2.legend(fontsize=12, loc='best')
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.set_xlim(0, None)
    ax2.set_ylim(0, 1.1)
    ax2.axhline(y=1.0, color='gray', linestyle=':', linewidth=1, alpha=0.5)
    
    plt.tight_layout()
    plt.savefig('greybody_factors_analysis.png', dpi=300, bbox_inches='tight')
    print("Saved greybody analysis to: greybody_factors_analysis.png")
    
    plt.show()

def diagnose_spectrum():
    """Diagnose potential issues with the spectrum calculation"""
    
    print("\n" + "="*70)
    print("SPECTRUM DIAGNOSTIC")
    print("="*70)
    
    M = 1.0
    T_H = 1.0 / (8.0 * np.pi * M)
    
    files = ['scalar_hawking_spectrum.dat', 'photon_hawking_spectrum.dat', 'graviton_hawking_spectrum.dat']
    
    for filename in files:
        if not os.path.exists(filename):
            continue
        
        data = load_spectrum(filename)
        if data is None:
            continue
        
        omega = data[:, 0]
        power = data[:, 1]
        
        print(f"\n{filename}:")
        print(f"  Frequency range: Mω = {omega[0]:.4f} to {omega[-1]:.4f}")
        print(f"  Number of points: {len(omega)}")
        
        # Find peak
        peak_idx = np.argmax(power)
        print(f"  Peak at: Mω = {omega[peak_idx]:.4f}, Power = {power[peak_idx]:.6e}")
        
        # Check if monotonic decrease (WRONG - should have a peak!)
        if peak_idx <= 2:
            print(f"  ⚠️  WARNING: Peak at lowest frequencies!")
            print(f"     This suggests:")
            print(f"     - Insufficient l_max in calculation (need higher angular momenta)")
            print(f"     - Or wrong low-frequency behavior in absorption cross-section")
            print(f"     - Page's spectra peak around Mω ~ 0.1-0.4, not at Mω→0")
        
        # Check thermal factor behavior
        thermal = 1.0 / (np.exp(omega / T_H) - 1.0)
        print(f"  Thermal factor at peak: {thermal[peak_idx]:.6e}")
        print(f"  Thermal factor at Mω=0.1: {thermal[np.argmin(np.abs(omega-0.1))]:.6e}")
        
        # Estimate if ω² × thermal has right behavior
        expected_shape = omega**2 * thermal  # Without greybody
        expected_peak_idx = np.argmax(expected_shape)
        print(f"  Expected peak (no greybody): Mω = {omega[expected_peak_idx]:.4f}")
        
    print("="*70 + "\n")

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
        
        elif command == 'greybody':
            plot_with_greybody_factors()
        
        elif command == 'diagnose':
            diagnose_spectrum()
        
        elif command == 'help':
            print("Usage: python3 normalise_graph.py [command]")
            print("\nCommands:")
            print("  page       - Plot with Page normalization (default)")
            print("  neutrinos  - Include approximate neutrino curve")
            print("  compare    - Compare raw vs normalized")
            print("  greybody   - Show greybody suppression factors")
            print("  diagnose   - Diagnose spectrum calculation issues")
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