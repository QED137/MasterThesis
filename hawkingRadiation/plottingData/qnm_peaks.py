#Hawking emission often shows resonance bumps near quasinormal mode frequencies.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

df = pd.read_csv("page_spectrum_all_species.csv")

# Choose gravitational waves only, since QNMs strongest there
grav = df[df["species"] == -2]
spec = grav.groupby("omega")["dEdtdw"].sum()

omega = spec.index.values
power = spec.values

# Find peaks
peaks, _ = find_peaks(power, height=np.max(power)*0.2, distance=10)

print("Detected QNM-like peaks at:")
for p in peaks:
    print("ω =", omega[p])

plt.figure(figsize=(10,6))
plt.plot(omega, power, label="Gravitational Spectrum")
plt.plot(omega[peaks], power[peaks], "ro", label="Detected peaks")
plt.xlabel(r"$\omega M$")
plt.ylabel("Emission power")
plt.title("QNM Peak Detection in Hawking Emission Spectrum")
plt.grid(True)
plt.legend()
plt.show()
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

df = pd.read_csv("page_spectrum_all_species.csv")

# Choose gravitational waves only, since QNMs strongest there
grav = df[df["species"] == -2]
spec = grav.groupby("omega")["dEdtdw"].sum()

omega = spec.index.values
power = spec.values

# Find peaks
peaks, _ = find_peaks(power, height=np.max(power)*0.2, distance=10)

print("Detected QNM-like peaks at:")
for p in peaks:
    print("ω =", omega[p])

plt.figure(figsize=(10,6))
plt.plot(omega, power, label="Gravitational Spectrum")
plt.plot(omega[peaks], power[peaks], "ro", label="Detected peaks")
plt.xlabel(r"$\omega M$")
plt.ylabel("Emission power")
plt.title("QNM Peak Detection in Hawking Emission Spectrum")
plt.grid(True)
plt.legend()
plt.show()

