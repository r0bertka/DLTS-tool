import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.signal import savgol_filter
from scipy.constants import elementary_charge as e
from scipy.constants import epsilon_0 as eps0
from scipy.constants import k as kB

from src.calculate_doping_profile import calculate_doping_profile

path = os.path.join('..', 'data', 'F_DLTSprofile_600K_Obelix_-19V.txt')
Vr, C, Vp, S = np.loadtxt(path, unpack=True, skiprows=2)

# parameters
eps = 9.7
A = 7.85e-7# in m2
T = 600  # in K
NC = 3.25e21 * T**(3/2) # in m-3
Et = 1.55  # trap energy in eV

number_of_data_points = len(Vr)

# calculate doping profile and bulk values of C, w and n
w, n_smooth = calculate_doping_profile(Vr, C, eps, A)
w = w[:-1]
n_smooth = n_smooth[:-1]
n_bulk = n_smooth[-1]
w0 = w[-1]
C0 = C[-1]

EF = - kB * T / e * np.log(n_bulk/NC)  # calculate Fermi level
Lambda = np.sqrt(2 * eps0 * eps / (e * n_bulk) * (Et-EF))  # calculate Lambda
x = w - Lambda  # calculate spatial coordinate of defect response

# flip DLTS measurement order s.th. response from deep inside semiconductor comes last, same as in CV measurements. Necessary for Obelix
Vp = Vp[:-1]
S = S[:-1]
reversed_Vp = Vp[::-1]
reversed_S = S[::-1]

# now the trap profile
smoothed_reversed_S = savgol_filter(reversed_S, 5, 3)
delta_C_C = smoothed_reversed_S / (C0 * 1e12)
diff_delta_C_C = np.diff(delta_C_C) / np.diff(Vp)
diff_delta_C_C = np.append(diff_delta_C_C, diff_delta_C_C[-1])
N_trap = -(e * w0**2 * n_bulk / (eps0 * eps)) * n_smooth * diff_delta_C_C / 2.449 * w / x  # 2.449 is the correction factor for td = 5ms, ti = 640ms, lock-in correlation
smoothed_N_trap = savgol_filter(N_trap, 7, 3)

plt.plot(w / 1e-6, n_smooth / 1e6, x / 1e-6, N_trap / 1e6)
plt.xlabel('w (µm)')
plt.ylabel('density (cm$^{-3}$)')
plt.yscale('log')
plt.show()

