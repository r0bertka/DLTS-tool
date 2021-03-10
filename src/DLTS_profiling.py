import numpy as np
import matplotlib.pyplot as plt
import os

from src.calculate_doping_profile import calculate_doping_profile
from src.calculate_trap_profile import calculate_trap_profile

path = os.path.join('..', 'data', '1e10_all_connected_after300C_wirebonded_295K_DLTS_profile_to2.2V.txt')
reverse_voltage, capacitance, pulse_voltage, dlts_signal = np.loadtxt(path, unpack=True, skiprows=2)

# parameters
eps = 9.7
area = 8.73e-7  # in m2
temperature = 295  # in K
effective_density_of_states_cb = 3.25e21 * temperature**(3/2)  # in m-3
trap_energy = 0.74  # trap energy in eV
setup = 'Obelix'

# calculate doping profile
scr_depth, doping_density = calculate_doping_profile(reverse_voltage, capacitance, eps, area)


pulse_voltage = pulse_voltage[:-1]  # cut off last "data point" that is per definition NaN because the Vp-S series has one point less than the Vr-C series
dlts_signal = dlts_signal[:-1]
if setup == 'Obelix':
    pulse_voltage = pulse_voltage[::-1]  # flip DLTS measurement order s.th. response from deep inside semiconductor comes last, same as in CV measurements. Necessary for Obelix
    dlts_signal = dlts_signal[::-1]

# now the trap profile
depth_trap_response, trap_density = calculate_trap_profile(capacitance, pulse_voltage, dlts_signal, scr_depth, doping_density, eps, temperature, effective_density_of_states_cb, trap_energy)

fig, doping_axis = plt.subplots()

color = 'blue'
doping_axis.set_xlabel('spatial position (µm)')
doping_axis.set_ylabel('doping density (cm$^{-3}$)', color = color)
doping_axis.set_yscale('log')
doping_axis.plot(scr_depth / 1e-6, doping_density / 1e6, color=color)
doping_axis.tick_params(axis='y', labelcolor = color)

trap_axis = doping_axis.twinx()

color = 'red'
trap_axis.set_ylabel('trap density (cm$^{-3}$)', color = color)
trap_axis.set_yscale('log')
trap_axis.plot(depth_trap_response / 1e-6, trap_density / 1e6, color=color)
trap_axis.tick_params(axis='y', labelcolor = color)

fig.tight_layout()

plt.show()

