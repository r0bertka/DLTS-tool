import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
from scipy.constants import elementary_charge as e
from scipy.constants import epsilon_0 as eps0
from scipy.constants import k as kB


path = os.path.join('..', 'data', 'He4_#1_after1000C_6h_Z12profile_296K_0.1Vstep_to+1Vp_newObelix.txt')
# names = ['Vr', 'C', 'Vp', 'S']
# col_index = range(4)
# data = pd.read_csv(path, delimiter='\t', skiprows=2, header=None, names=names)
Vr, C, Vp, S = np.loadtxt(path, unpack=True, skiprows=2)
# Vr, C, Vp, S = data['Vr'], data['C'], data['Vp'], data['S']

# some constants
eps = 9.7
A = 7.85e-7# in m2
T = 295  # in K
NC = 3.25e21 * T**(3/2) # in m-3
Et = 0.67  # trap energy in eV

number_of_data_points = len(Vr)

# calculate doping profile
w = A * eps0 * eps / C
# w = w[:-1]
inverse_square_of_cap = 1 / C**2
diff_inverse_square_of_cap = np.diff(inverse_square_of_cap) / np.diff(Vr)
diff_inverse_square_of_cap = np.append(diff_inverse_square_of_cap, diff_inverse_square_of_cap[-1])
smooth_diff_inverse_square_of_cap = savgol_filter(diff_inverse_square_of_cap, 9, 3)
# Vr = Vr[:-1]
# n = -2/(e * A**2 * eps0 * eps * diff_invsqrC)
n_smooth = -2/(e * A**2 * eps0 * eps * smooth_diff_inverse_square_of_cap)

# calculate Fermi level
n_bulk = n_smooth[-1]
w0 = w[-1]
C0 = C[-1]
EF = - kB * T / e * np.log(n_bulk/NC)

# calculate Lambda
Lambda = np.sqrt(2 * eps0 * eps / (e * n_bulk) * (Et-EF))

# calculate spatial coordinate of defect response
x = w - Lambda

# flip DLTS measurement order s.th. response from deep inside semiconductor comes last, same as in CV measurements. Necessary for Obelix
reversed_Vp = Vp[::-1]
reversed_S = S[::-1]

# now the trap profile
delta_C_C = reversed_S / (C0 * 1e12)
diff_delta_C_C = np.diff(delta_C_C) / np.diff(Vp)
diff_delta_C_C = np.append(diff_delta_C_C, diff_delta_C_C[-1])
N_trap = -(e * w0**2 * n_bulk / (eps0 * eps)) * n_smooth * diff_delta_C_C / 2.449 * w / x  # 2.449 is the correction factor for td = 5ms, ti = 640ms, lock-in correlation
smoothed_N_trap = savgol_filter(N_trap, 7, 3)

plt.plot(w / 1e-6, n_smooth / 1e6, x / 1e-6, smoothed_N_trap / 1e6)
plt.xlabel('w (µm)')
plt.ylabel('density (cm$^{-3}$)')
plt.yscale('log')
plt.show()

