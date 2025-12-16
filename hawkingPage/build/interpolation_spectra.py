#!/usr/bin/env python3
"""
Interpolate Hawking spectra using cubic splines
Mimics Mathematica's Interpolation[..., InterpolationOrder->3]
"""

import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import os
import sys

def interpolate_spectrum(input_file, output_file, n_points=1000, plot=False):
    """
    Interpolate spectrum data using cubic splines
    
    Args:
        input_file: Original discrete data file
        output_file: Output file for interpolated data
        n_points: Number of points in interpolated curve
        plot: Whether to create comparison plot
    
    Returns:
        (omega_dense, power_dense) arrays
    """
    
    if not os.path.exists(input_file):
        print(f"WARNING: {input_file} not found, skipping")
        return None, None
    
    # Load data
    data = np.loadtxt(input_file, comments='#')
    if len(data) == 0:
        print(f"WARNING: {input_file} is empty, skipping")
        return None, None
    
    omega = data[:, 0]
    power = data[:, 1]
    
    n_orig = len(omega)
    print(f"  Original points: {n_orig}")
    
    # Create cubic spline
    # bc_type='natural' gives natural boundary conditions (2nd derivative = 0 at ends)
    # This matches Mathematica's default behavior
    cs = CubicSpline(omega, power, bc_type='natural')
    
    # Create dense grid
    omega_dense = np.linspace(omega[0], omega[-1], n_points)
    power_dense = cs(omega_dense)
    
    # Ensure no negative values (can happen with interpolation at edges)
    power_dense = np.maximum(power_dense, 0.0)
    
    # Save interpolated data
    with open(output_file, 'w') as f:
        f.write(f"# Cubic spline interpolation of {input_file}\n")
        f.write(f"# Original: {n_orig} points → Interpolated: {n_points} points\n")
        f.write(f"# Method: SciPy CubicSpline with natural boundary conditions\n")
        f.write(f"# Mimics Mathematica InterpolationOrder->3\n")
        f.write("# omega    power\n")
        for w, p in zip(omega_dense, power_dense):
            f.write(f"{w:.8f}    {p:.12e}\n")
    
    print(f"  Interpolated points: {n_points}")
    print(f"  Saved to {output_file}")
    
    # Optional: create comparison plot
    if plot:
        plot_filename = output_file.replace('.dat', '_comparison.png')
        
        plt.figure(figsize=(10, 6))
        plt.plot(omega, power, 'o', markersize=6, label=f'Original ({n_orig} points)', alpha=0.7)
        plt.plot(omega_dense, power_dense, '-', linewidth=2, label=f'Cubic Spline ({n_points} points)')
        
        plt.xlabel('ω', fontsize=12)
        plt.ylabel('Power', fontsize=12)
        plt.title(f'Interpolation: {os.path.basename(input_file)}', fontsize=14)
        plt.legend(fontsize=11)
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(plot_filename, dpi=200)
        print(f"  Comparison plot: {plot_filename}")
        plt.close()
    
    return omega_dense, power_dense

def interpolate_all_spectra(n_points=1000, make_plots=False):
    """Interpolate all available spectrum files"""
    
    files = {
        'scalar_hawking_spectrum.dat': 'Scalar (s=0)',
        'photon_hawking_spectrum.dat': 'Photon (s=1)',
        'graviton_hawking_spectrum.dat': 'Graviton (s=2)',
        'neutrino_hawking_spectrum.dat': 'Neutrino (s=1/2)'
    }
    
    print("\n" + "="*70)
    print("Cubic Spline Interpolation of Hawking Spectra")
    print("="*70 + "\n")
    
    results = {}
    
    for filename, description in files.items():
        if os.path.exists(filename):
            output = filename.replace('.dat', '_interpolated.dat')
            print(f"{description}:")
            omega, power = interpolate_spectrum(filename, output, n_points, plot=make_plots)
            if omega is not None:
                results[filename] = (omega, power)
            print()
        else:
            print(f"{description}: Not found, skipping\n")
    
    return results

def create_comparison_plot(results):
    """Create a plot comparing original and interpolated data"""
    
    if not results:
        print("No data to plot")
        return
    
    fig, axes = plt.subplots(len(results), 1, figsize=(10, 4*len(results)))
    if len(results) == 1:
        axes = [axes]
    
    colors = {
        'scalar_hawking_spectrum.dat': 'blue',
        'photon_hawking_spectrum.dat': 'green',
        'graviton_hawking_spectrum.dat': 'red',
        'neutrino_hawking_spectrum.dat': 'purple'
    }
    
    labels = {
        'scalar_hawking_spectrum.dat': 'Scalar (s=0)',
        'photon_hawking_spectrum.dat': 'Photon (s=1)',
        'graviton_hawking_spectrum.dat': 'Graviton (s=2)',
        'neutrino_hawking_spectrum.dat': 'Neutrino (s=1/2)'
    }
    
    for ax, (filename, (omega_interp, power_interp)) in zip(axes, results.items()):
        # Load original data
        data_orig = np.loadtxt(filename, comments='#')
        omega_orig = data_orig[:, 0]
        power_orig = data_orig[:, 1]
        
        color = colors.get(filename, 'black')
        label = labels.get(filename, filename)
        
        # Plot
        ax.plot(omega_orig, power_orig, 'o', color=color, markersize=5, 
               alpha=0.6, label=f'Original ({len(omega_orig)} pts)')
        ax.plot(omega_interp, power_interp, '-', color=color, linewidth=2,
               label=f'Interpolated ({len(omega_interp)} pts)')
        
        ax.set_xlabel('ω', fontsize=11)
        ax.set_ylabel('Power', fontsize=11)
        ax.set_title(label, fontsize=13, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('all_spectra_interpolation_comparison.png', dpi=250, bbox_inches='tight')
    print("="*70)
    print("Saved combined comparison: all_spectra_interpolation_comparison.png")
    print("="*70 + "\n")

def verify_interpolation_quality(filename, n_test=10):
    """
    Verify interpolation quality by checking intermediate points
    """
    
    if not os.path.exists(filename):
        return
    
    # Load original data
    data = np.loadtxt(filename, comments='#')
    omega = data[:, 0]
    power = data[:, 1]
    
    # Create spline
    cs = CubicSpline(omega, power, bc_type='natural')
    
    # Check at midpoints between original points
    print(f"Quality check for {filename}:")
    print(f"  Checking {n_test} random midpoints...")
    
    max_error = 0.0
    for i in range(min(n_test, len(omega)-1)):
        idx = np.random.randint(0, len(omega)-1)
        
        # Midpoint
        w_mid = 0.5 * (omega[idx] + omega[idx+1])
        
        # True value (linear interpolation as reference)
        p_true = 0.5 * (power[idx] + power[idx+1])
        
        # Spline value
        p_spline = cs(w_mid)
        
        # Relative error
        if p_true > 0:
            error = abs(p_spline - p_true) / p_true
            max_error = max(max_error, error)
    
    print(f"  Max relative error: {100*max_error:.2f}%")
    if max_error < 0.05:
        print(f"  ✓ Good quality (error < 5%)")
    elif max_error < 0.10:
        print(f"  ⚠ Acceptable (error < 10%)")
    else:
        print(f"  ⚠ Warning: High error, may need more points")
    print()

def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Interpolate Hawking spectra using cubic splines'
    )
    parser.add_argument('--points', type=int, default=1000,
                       help='Number of points in interpolated curve (default: 1000)')
    parser.add_argument('--plots', action='store_true',
                       help='Create comparison plots')
    parser.add_argument('--verify', action='store_true',
                       help='Verify interpolation quality')
    
    args = parser.parse_args()
    
    # Interpolate all spectra
    results = interpolate_all_spectra(n_points=args.points, make_plots=args.plots)
    
    # Create combined comparison plot
    if results and args.plots:
        create_comparison_plot(results)
    
    # Verify quality
    if args.verify:
        print("="*70)
        print("Interpolation Quality Verification")
        print("="*70 + "\n")
        for filename in results.keys():
            verify_interpolation_quality(filename)
    
    # Summary
    print("="*70)
    print("Summary")
    print("="*70)
    print(f"\nInterpolated {len(results)} spectra with {args.points} points each")
    print("\nOutput files:")
    for filename in results.keys():
        output = filename.replace('.dat', '_interpolated.dat')
        print(f"  ✓ {output}")
    
    print("\nUsage:")
    print("  - Use *_interpolated.dat files for smooth plotting")
    print("  - Curves will match Mathematica's InterpolationOrder->3")
    print("  - Compatible with all plotting scripts")
    print("\nExample:")
    print("  python3 plot_page_1976_corrected.py  # (modify to use _interpolated.dat)")
    print("="*70 + "\n")

if __name__ == '__main__':
    main()