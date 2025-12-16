# How to Run

## 1. Build
```bash
mkdir build && cd build
cmake ..
make -j4
```

## 2. Run Simulations
The executable computes Hawking radiation spectra for different particle types.

### Default (Page's main results: photon + graviton + neutrino)
```bash
./bin/hawking_radiation
```

### Individual particle types
```bash
./bin/hawking_radiation scalar    # Scalar field (spin-0)
./bin/hawking_radiation photon    # Photons (spin-1)
./bin/hawking_radiation graviton  # Gravitons (spin-2)
./bin/hawking_radiation neutrino  # Neutrinos (spin-1/2, analytical)
```

### All particle types
```bash
./bin/hawking_radiation all       # Compute all four spectra
./bin/hawking_radiation page      # Page's three: photon + graviton + neutrino
./bin/hawking_radiation quick     # Quick test only (no full computation)
```

## 3. Generate Smooth Interpolated Data
```bash
python3 interpolation_spectra.py
```
This creates interpolated versions of the spectrum files:
- `scalar_hawking_spectrum_interpolated.dat`
- `photon_hawking_spectrum_interpolated.dat`
- `graviton_hawking_spectrum_interpolated.dat`
- `neutrino_hawking_spectrum_interpolated.dat`

## 4. Plot Results (Page 1976 Style)
```bash
python3 plot_page_1976_exact.py page
```

This generates a plot matching Don Page (1976) Figure 1 with:
- Dual y-axes (dimensionless and physical units)
- Top x-axis in MeV
- Color-coded particle spectra:
  - Green: Scalar (s=0)
  - Blue: Photons (s=1)
  - Red: Gravitons (s=2)
  - Magenta (dashed): Neutrinos (s=1/2)
  - Black (thick): Total spectrum