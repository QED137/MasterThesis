# #evaporation_curve.py
# #black hole evaporation curve using Page's model
# #!/usr/bin/env python3
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt

# df = pd.read_csv("page_spectrum_all_species.csv")

# # integrate over omega to get total luminosity
# L = df.groupby("omega")["dEdtdw"].sum()

# total_power = np.trapz(L.values, L.index)

# print("Total luminosity dM/dt =", total_power)

# # Now integrate evaporation
# M0 = 1.0   # initial mass
# M = M0
# t = 0
# dt = 0.001  # small time step

# Ts = [t]
# Ms = [M]

# while M > 0.05:
#     L_M = total_power / (M**2)   # Hawking power scales ~ 1/M^2
#     M -= L_M * dt
#     t += dt
#     Ts.append(t)
#     Ms.append(M)

# plt.figure(figsize=(10,6))
# plt.plot(Ts, Ms)
# plt.xlabel("t (in units of M)")
# plt.ylabel("Mass M(t)")
# plt.title("Black Hole Evaporation Curve (Don Page Model)")
# plt.grid(True)
# plt.show()


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv("page_spectrum_all_species.csv")

# integrate over omega to get total luminosity (dimensionless units)
L = df.groupby("omega")["dEdtdw"].sum()
total_power = np.trapz(L.values, L.index)  # Page's spectrum
print("Total luminosity dM/dt =", total_power)

# initial mass
M = 1.0
t = 0.0

# store only 200 points
times = []
masses = []

N_points = 200
for _ in range(N_points):
    times.append(t)
    masses.append(M)

    # Hawking mass loss ~ L/M^2
    dMdt = -total_power / (M**2)

    # adaptive timestep: larger when BH is large, tiny near the end
    dt = min(0.01 * M**2, 0.005)

    M += dMdt * dt
    t += dt

    if M <= 0.05:
        break

plt.figure(figsize=(10,6))
plt.plot(times, masses)
plt.xlabel("t (in units of M)")
plt.ylabel("Mass M(t)")
plt.title("Black Hole Evaporation Curve (Fast Version)")
plt.grid(True)
plt.show()
