import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.signal import savgol_filter

path = os.path.join('..', 'data', 'He4_#1_after1000C_6h_Z12profile_296K_0.1Vstep_to+1Vp_newObelix.txt')
names = ['Vr', 'C', 'Vp', 'S']
# col_index = range(4)
data = pd.read_csv(path, delimiter='\t', skiprows=2, header=None, names=names)
Vr, C, Vp, S = data['Vr'], data['C'], data['Vp'], data['S']

# some constants
eps0 = 8.85e-12 # in As/Vm
eps = 9.7
e = 1.602e-19 # in C
kB = 1.38e-23 # in J/K
A = 7.85e-7# in m2
T = 295 # in K
NC = 3.25e21 * T**(3/2) # in m-3
Et = 0.67 # trap energy in eV

# calculate doping profile
w = A*eps0*eps/C
w = w[:-1]
invsqrC = 1/C**2
diff_invsqrC = np.diff(invsqrC)/np.diff(Vr)
smooth_diff_invsqrC = savgol_filter(diff_invsqrC, 9, 3)
Vr = Vr[:-1]
# n = -2/(e * A**2 * eps0 * eps * diff_invsqrC)
n_smooth = -2/(e * A**2 * eps0 * eps * smooth_diff_invsqrC)

# calculate Fermi level
n_bulk = n_smooth[-1]
w0 = w[-1]
EF = - kB * T / e * np.log(n_bulk/NC)

# calculate Lambda
Lambda = np.sqrt(2 * eps0 * eps / (e * n_bulk) * (Et-EF))

# calculate spatial coordinate of defect response
x = w - Lambda
x = x[:-1]

# now the trap profile
diff_S = np.diff(S)/np.diff(Vp)
n_smooth = n_smooth[:-1]
N_trap = -(e * w0**2 * n_bulk / (eps0 * eps)) * n_smooth * diff_S * (w0**2 / ((w0 - Lambda)**2 - (w - Lambda)**2))

plt.plot(x, N_trap)
#plt.xlabel('w (µm)')
#plt.ylabel('net doping density $N_D$ (cm$^{-3}$)')
plt.show()

# data.plot(x='Vp', y='S')
# plt.show()


