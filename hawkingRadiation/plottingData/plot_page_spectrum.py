#page_spectrum.py
#energy spectrum of Hawking radiation from C++ output
#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load C++ spectrum
df = pd.read_csv("page_spectrum_all_species.csv")

# Unique species
species = sorted(df["species"].unique(), key=str)

plt.figure(figsize=(10,6))

for sp in species:
    sub = df[df["species"] == sp]
    
    # Page spectrum summed over (ℓ, m)
    grouped = sub.groupby("omega")["dEdtdw"].sum()

    plt.plot(grouped.index, grouped.values, label=f"s={sp}")

plt.xlabel(r"$\omega M$")
plt.ylabel(r"$dE/dt\, d\omega$ (dimensionless)")
plt.title("Hawking Emission Spectrum (Reproducing Page 1976)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

