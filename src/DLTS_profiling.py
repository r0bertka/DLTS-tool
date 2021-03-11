import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

from src.calculate_carrier_profile import calculate_carrier_profile
from src.calculate_trap_profile import calculate_trap_profile

raw_data_filename = 'B_#1_DLTS_profilling_650K_EH67_secondAttempt'
raw_data_filename_extension = raw_data_filename + '.txt'
raw_data_path = os.path.join('..', 'data', raw_data_filename_extension)
reverse_voltage, capacitance, pulse_voltage, dlts_signal = np.loadtxt(raw_data_path, unpack=True, skiprows=2)

# parameters
eps = 9.7
area = 7.85e-7  # in m2
temperature = 600  # in K
effective_density_of_states_cb = 3.25e21 * temperature**(3/2)  # in m-3
trap_energy = 1.55  # trap energy in eV
setup = 'Obelix'

# calculate doping profile
scr_depth, carrier_density = calculate_carrier_profile(reverse_voltage, capacitance, eps, area)
scr_depth = scr_depth[:-1]
carrier_density = carrier_density[:-1]

pulse_voltage = pulse_voltage[:-1]  # cut off last "data point" that is per definition NaN because the Vp-S series has one point less than the Vr-C series
dlts_signal = dlts_signal[:-1]

if setup == 'Obelix':
    pulse_voltage = pulse_voltage[::-1]  # flip DLTS measurement order s.th. response from deep inside semiconductor comes last, same as in CV measurements. Necessary for Obelix
    dlts_signal = dlts_signal[::-1]

# now the trap profile
depth_trap_response, trap_density = calculate_trap_profile(capacitance, pulse_voltage, dlts_signal, scr_depth, carrier_density, eps, temperature, effective_density_of_states_cb, trap_energy)

# pd.DataFrame((scr_depth / 1e-6, carrier_density / 1e6, depth_trap_response / 1e-6, trap_density / 1e6), columns=['SCR depth (um)', 'carrier density (cm-3)', 'depth trap response (um)', 'trap density (cm-3)'])
eval_data = pd.DataFrame({'SCR depth (um)': scr_depth / 1e-6,
                          'carrier density (cm-3)': carrier_density / 1e6,
                          'depth trap response (um)': depth_trap_response / 1e-6,
                          'trap density (cm-3)': trap_density / 1e6}, columns = ['SCR depth (um)', 'carrier density (cm-3)', 'depth trap response (um)', 'trap density (cm-3)'])
eval_data_filename = raw_data_filename + '_eval.txt'
eval_data_path = os.path.join('..', 'data', 'eval', eval_data_filename)

eval_data.to_csv(eval_data_path, sep='\t', float_format='%1.6E', index=False)


# fig, doping_axis = plt.subplots()
#
# color = 'blue'
# doping_axis.set_xlabel('spatial position (µm)')
# doping_axis.set_ylabel('doping density (cm$^{-3}$)', color = color)
# doping_axis.set_yscale('log')
# doping_axis.plot(scr_depth / 1e-6, carrier_density / 1e6, color=color)
# doping_axis.tick_params(axis='y', labelcolor = color)
#
# trap_axis = doping_axis.twinx()
#
# color = 'red'
# trap_axis.set_ylabel('trap density (cm$^{-3}$)', color = color)
# trap_axis.set_yscale('log')
# trap_axis.plot(depth_trap_response / 1e-6, trap_density / 1e6, color=color)
# trap_axis.tick_params(axis='y', labelcolor = color)
#
# fig.tight_layout()
#
# plt.show()

