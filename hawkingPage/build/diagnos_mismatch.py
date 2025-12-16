#!/usr/bin/env python3
"""
Compare your computed values with Page (1976) Table I exact values
Table I shows Γ values (absorption rates) for different particles and l modes
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# Page's Table I - Selected values (from paper page 6)
# Format: {particle: {l: Gamma_value}}
# These are Γ values (dimensionless absorption rates) at specific Mω values

# At Mω = 0.1:
PAGE_TABLE_I_0p1 = {
    'photon': {
        1: 0.0065,   # l=1
        2: 0.00037,  # l=2
        3: 0.000039  # l=3
    },
    'graviton': {
        2: 0.00015,  # l=2
        3: 0.000017, # l=3
        4: 0.0000025 # l=4
    }
}

# At Mω = 0.2:
PAGE_TABLE_I_0p2 = {
    'photon': {
        1: 0.039,    # l=1
        2: 0.0060,   # l=2
        3: 0.0012    # l=3
    },
    'graviton': {
        2: 0.0033,   # l=2
        3: 0.00063,  # l=3
        4: 0.00015   # l=4
    }
}

def load_your_data(filename):
    """Load your computed spectrum"""
    if not os.path.exists(filename):
        print(f"File not found: {filename}")
        return None, None
    
    data = np.loadtxt(filename, comments='#')
    return data[:, 0], data[:, 1]

def find_closest_omega(omega_array, target_omega):
    """Find index of closest omega value"""
    idx = np.argmin(np.abs(omega_array - target_omega))
    return idx, omega_array[idx]

def compare_with_page_table(omega, power, particle_name, target_mw=0.1):
    """Compare values at specific Mω with Page's Table I"""
    
    print(f"\n{'='*70}")
    print(f"Comparing {particle_name.upper()} at Mω ≈ {target_mw}")
    print(f"{'='*70}")
    
    # Find closest point in your data
    idx, actual_mw = find_closest_omega(omega, target_mw)
    your_power = power[idx]
    
    print(f"\nYour data:")
    print(f"  Closest Mω = {actual_mw:.4f}")
    print(f"  Power = {your_power:.6e}")
    
    # Get Page's values for this Mω
    if target_mw == 0.1:
        page_table = PAGE_TABLE_I_0p1
    elif target_mw == 0.2:
        page_table = PAGE_TABLE_I_0p2
    else:
        print(f"No Page data for Mω = {target_mw}")
        return
    
    if particle_name not in page_table:
        print(f"No Page data for {particle_name}")
        return
    
    print(f"\nPage's Table I values:")
    page_values = page_table[particle_name]
    for l, gamma in page_values.items():
        print(f"  l={l}: Γ = {gamma:.6f}")
    
    # Note: Direct comparison difficult because Page shows Γ (absorption)
    # while we compute Power = Σ(2l+1) × Γ × thermal × ω²/(2π)
    print(f"\nNote: Direct comparison requires converting Power → Γ_l")
    print(f"      This needs decomposing your total into individual l contributions")

def diagnose_mismatch():
    """Diagnose why spectra don't match"""
    
    print("\n" + "="*70)
    print("SPECTRUM MISMATCH DIAGNOSIS")
    print("="*70)
    
    files = {
        'photon': 'photon_hawking_spectrum.dat',
        'graviton': 'graviton_hawking_spectrum.dat',
        'neutrino': 'neutrino_hawking_spectrum.dat'
    }
    
    data = {}
    for name, filename in files.items():
        omega, power = load_your_data(filename)
        if omega is not None:
            data[name] = (omega, power)
            print(f"✓ Loaded {filename}")
        else:
            print(f"✗ Missing {filename}")
    
    if len(data) < 2:
        print("\nERROR: Need at least 2 spectra files")
        return
    
    # Check 1: Are curves smooth?
    print("\n" + "="*70)
    print("CHECK 1: Smoothness")
    print("="*70)
    
    for name, (omega, power) in data.items():
        # Compute derivative to check smoothness
        if len(omega) > 2:
            d2power = np.diff(power, 2)
            roughness = np.std(d2power)
            print(f"\n{name.capitalize()}:")
            print(f"  Points: {len(omega)}")
            print(f"  Roughness: {roughness:.6e}")
            if roughness > 1e-6:
                print(f"  ⚠ Curve may be rough - consider interpolation")
            else:
                print(f"  ✓ Curve is smooth")
    
    # Check 2: Are peaks in right locations?
    print("\n" + "="*70)
    print("CHECK 2: Peak Locations")
    print("="*70)
    
    expected_peaks = {
        'neutrino': 0.18,
        'photon': 0.15,
        'graviton': 0.12
    }
    
    for name, (omega, power) in data.items():
        peak_idx = np.argmax(power)
        peak_omega = omega[peak_idx]
        peak_power = power[peak_idx]
        expected = expected_peaks.get(name, 0.15)
        
        print(f"\n{name.capitalize()}:")
        print(f"  Your peak: Mω = {peak_omega:.4f}, Power = {peak_power:.6e}")
        print(f"  Expected:  Mω ≈ {expected:.2f}")
        
        error = abs(peak_omega - expected)
        if error < 0.03:
            print(f"  ✓ Peak location correct (error = {error:.4f})")
        elif error < 0.05:
            print(f"  ⚠ Peak location acceptable (error = {error:.4f})")
        else:
            print(f"  ❌ Peak location wrong (error = {error:.4f})")
    
    # Check 3: Are amplitudes correct?
    print("\n" + "="*70)
    print("CHECK 3: Relative Amplitudes")
    print("="*70)
    
    if 'photon' in data and 'graviton' in data:
        omega_p, power_p = data['photon']
        omega_g, power_g = data['graviton']
        
        # At peak
        peak_p = np.max(power_p)
        peak_g = np.max(power_g)
        ratio = peak_p / peak_g
        expected_ratio = 16.7 / 1.9  # From percentages
        
        print(f"\nPhoton/Graviton peak ratio:")
        print(f"  Your value: {ratio:.2f}:1")
        print(f"  Expected:   {expected_ratio:.2f}:1")
        
        error = abs(ratio - expected_ratio)
        if error < 2:
            print(f"  ✓ Ratio correct")
        else:
            print(f"  ❌ Ratio wrong (error = {error:.2f})")
            print(f"\n  Possible causes:")
            print(f"    - Wrong l_min or l_max for one of the particles")
            print(f"    - Absorption calculation error")
            print(f"    - Missing (2l+1) factors")
    
    # Check 4: Low-frequency behavior
    print("\n" + "="*70)
    print("CHECK 4: Low-Frequency Scaling")
    print("="*70)
    print("\nFor BOSONS, at low ω:")
    print("  Photons (s=1):   σ ∝ ω²")
    print("  Gravitons (s=2): σ ∝ ω⁴")
    
    for name, (omega, power) in data.items():
        if name in ['photon', 'graviton']:
            # Check first 5 points
            omega_low = omega[:5]
            power_low = power[:5]
            
            print(f"\n{name.capitalize()} (first 5 points):")
            
            if name == 'photon':
                # Should scale as ω²
                ratios = power_low / (omega_low**2 + 1e-20)
                std_ratio = np.std(ratios) / np.mean(ratios)
                print(f"  Power/ω² = {ratios}")
                print(f"  Variation: {100*std_ratio:.1f}%")
                if std_ratio < 0.3:
                    print(f"  ✓ Correct ω² scaling")
                else:
                    print(f"  ⚠ Scaling may be wrong")
                    
            elif name == 'graviton':
                # Should scale as ω⁴
                ratios = power_low / (omega_low**4 + 1e-20)
                std_ratio = np.std(ratios) / np.mean(ratios)
                print(f"  Power/ω⁴ = {ratios}")
                print(f"  Variation: {100*std_ratio:.1f}%")
                if std_ratio < 0.3:
                    print(f"  ✓ Correct ω⁴ scaling")
                else:
                    print(f"  ⚠ Scaling may be wrong")
    
    # Check 5: Integrated totals
    print("\n" + "="*70)
    print("CHECK 5: Integrated Powers and Percentages")
    print("="*70)
    
    totals = {}
    for name, (omega, power) in data.items():
        total = np.trapz(power, omega)
        totals[name] = total
        print(f"\n{name.capitalize()}: ∫ Power dω = {total:.6e}")
    
    grand_total = sum(totals.values())
    print(f"\nGrand total: {grand_total:.6e}")
    
    expected_pct = {
        'neutrino': 81.4,
        'photon': 16.7,
        'graviton': 1.9
    }
    
    print(f"\nPercentage contributions:")
    for name, total in totals.items():
        pct = 100 * total / grand_total
        exp_pct = expected_pct.get(name, 0)
        
        print(f"  {name.capitalize()}:")
        print(f"    Your value:  {pct:>6.2f}%")
        print(f"    Page's value: {exp_pct:>6.2f}%")
        
        error = abs(pct - exp_pct)
        if error < 2:
            print(f"    ✓ Match!")
        elif error < 5:
            print(f"    ⚠ Close (error = {error:.1f}%)")
        else:
            print(f"    ❌ Mismatch (error = {error:.1f}%)")
    
    # Final recommendations
    print("\n" + "="*70)
    print("RECOMMENDATIONS")
    print("="*70)
    
    # Analyze what's wrong
    if 'photon' in data and 'graviton' in data:
        omega_p, power_p = data['photon']
        omega_g, power_g = data['graviton']
        ratio = np.max(power_p) / np.max(power_g)
        
        if ratio < 5:
            print("\n❌ CRITICAL: Photon/Graviton ratio too small!")
            print("   This means gravitons are too strong relative to photons.")
            print("\n   Most likely cause:")
            print("   → Gravitons are including forbidden l=0 or l=1 modes")
            print("\n   Fix:")
            print("   → In compute_hawking_spectrum(), verify that:")
            print("      int l_min = std::abs(s);")
            print("   → For gravitons (s=2), l must start at 2, not 0 or 1!")
            
        elif ratio > 15:
            print("\n❌ CRITICAL: Photon/Graviton ratio too large!")
            print("   This means photons are too strong relative to gravitons.")
            print("\n   Possible causes:")
            print("   → Missing l=2 mode for gravitons")
            print("   → Wrong absorption calculation")
            
        else:
            print("\n✓ Photon/Graviton ratio looks reasonable")
    
    # Check if interpolation needed
    has_interpolated = any(os.path.exists(f.replace('.dat', '_interpolated.dat')) 
                          for f in files.values())
    
    if not has_interpolated:
        print("\n⚠ No interpolated data found")
        print("   Run: python3 interpolate_spectra.py")
        print("   This will smooth your curves to match Page's figure")

def main():
    diagnose_mismatch()
    
    print("\n" + "="*70)
    print("Next steps:")
    print("="*70)
    print("\n1. Review the checks above")
    print("2. If Photon/Graviton ratio is wrong:")
    print("   → Check l_min in C++ code")
    print("   → Verify for gravitons: l starts at 2, not 0")
    print("3. If peaks in wrong location:")
    print("   → Add more frequency points near peaks")
    print("4. If curves rough:")
    print("   → Run: python3 interpolate_spectra.py")
    print("5. After fixes:")
    print("   → Recompile: cd build && make")
    print("   → Rerun: ./bin/hawking_radiation page")
    print("   → Check: python3 diagnose_detailed.py")
    print("="*70 + "\n")

if __name__ == '__main__':
    main()