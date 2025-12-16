# Master Thesis: Black Hole Physics Simulations

This repository contains the computational implementation of my master thesis work on black hole physics, focusing on Hawking radiation, the Page curve, and solutions to the Teukolsky equation.

## Overview

This project provides numerical simulations and computational tools for studying:
- **Hawking Radiation**: Calculation of thermal emission from black holes

- **Teukolsky Equation**: Solutions for perturbations in Kerr black hole spacetimes

## Project Structure

### `/hawkingRadiation`
Simulation of Hawking radiation emission and greybody factors:
- `hawking.cpp/hpp`: Core Hawking radiation calculations
- `greybody.cpp/hpp`: Greybody factor computations
- `teukolsky_radial.cpp/hpp`: Radial solutions to the Teukolsky equation
- `page_spectrum_all_species.csv`: Calculated particle spectra data

### `/hawkingPage`
Page curve calculations for black hole evaporation:
- `black_hole.cpp/h`: Black hole thermodynamics implementation
- `bh_absorption`: Absorption cross-section calculations

### `/teukolskyEquation`
Numerical solutions to the Teukolsky equation:
- `teukolsky1D.cpp`: 1D radial equation solver
- `teukolsky2D.cpp`: 2D time-evolution solver
- `visualize.py`: Visualization tools for solution data
- `/data`: Time-series solution snapshots

### `/hawkingRadiation/plottingData`
Visualization scripts for analysis:
- `plot_page_spectrum.py`: Page curve spectrum plots
- `evaporation_curve.py`: Black hole evaporation timeline
- `qnm_peaks.py`: Quasi-normal mode analysis

## Dependencies

- **C++ Compiler**: C++11 or later
- **CMake**: Version 3.10+
- **Eigen3**: Linear algebra library
- **Python 3**: For visualization scripts
- **NumPy/Matplotlib**: Python plotting dependencies

## Building the Project

```bash
# Clone the repository
git clone https://github.com/QED137/MasterThesis.git
cd MasterThesis

# Build Hawking radiation simulations
cd hawkingRadiation
mkdir -p build && cd build
cmake ..
make -j

# Build Teukolsky equation solvers
cd ../../teukolskyEquation
g++ -std=c++11 -O3 teukolsky1D.cpp -o teukolsky1D
g++ -std=c++11 -O3 teukolsky2D.cpp -o teukolsky2D
```

## Usage

```bash
# Run Hawking radiation calculation
cd hawkingRadiation/build
./hawking_radiation

# Generate Page curve data
cd ../../hawkingPage/build
./page_curve

# Solve Teukolsky equation
cd ../../teukolskyEquation
./teukolsky1D
./teukolsky2D

# Visualize results
python visualize.py
cd ../hawkingRadiation/plottingData
python plot_page_spectrum.py
```

## Future Development

The goal is to develop this into a full interactive browser-based visualization tool using modern web frameworks, allowing real-time exploration of black hole physics simulations.

## License

See [LICENSE](LICENSE) for details.

## Contact

For questions or collaboration, please open an issue on this repository.
