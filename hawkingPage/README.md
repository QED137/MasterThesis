# How to run 
# 1. Build
mkdir build && cd build
cmake .. && make -j4

# 2. Run (computes ALL THREE spectra automatically)
./bin/hawking_radiation

# 3. Plot
cd ..
# Standard Page normalization (shows total matching ideal shape)
python3 ../normalise_graph.py page

# Show greybody suppression factors
python3 ../normalise_graph.py greybody

# Other options
python3 ../normalise_graph.py neutrinos  # Add neutrino approximation
python3 ../normalise_graph.py stats      # Print statistics only