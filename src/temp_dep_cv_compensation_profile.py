import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.constants import elementary_charge as e
from scipy.constants import k as kB

from src.calculate_carrier_profile import calculate_carrier_profile

eps = 9.7
area = 7.85e-7  # in m2
temperature = 295 # in K
trap_degeneracy = 2
trap_binding_energy =  0.072 # in eV
effective_density_of_states_cb = 3.25e21 * temperature**(3/2)  # in m-3
donor_density = 3e20  # in per cubic meter

raw_data_filename = 'B_as-rec_CV_RT_dark'
raw_data_filename_extension = raw_data_filename + '.txt'

path = os.path.join('..', 'data', raw_data_filename_extension)
reverse_voltage, capacitance, scr_width, carrier_density, conductivity = np.loadtxt(path, unpack=True, skiprows=2)

# calculated_scr_width, calculated_carrier_density = calculate_carrier_profile(reverse_voltage, capacitance, eps, area)

acceptor_density = (effective_density_of_states_cb * np.exp(-e * trap_binding_energy / (kB * temperature)) * (donor_density - carrier_density * 1e6) - trap_degeneracy * (carrier_density * 1e6) **2) / ((carrier_density * 1e6 * trap_degeneracy) + effective_density_of_states_cb * np.exp(-(trap_binding_energy * e) / (kB * temperature)))

eval_data = pd.DataFrame({'SCR depth (um)': scr_width / 1e-4,
                          'carrier density (cm-3)': carrier_density,
                          'acceptor density (cm-3)': acceptor_density / 1e6}, columns = ['SCR depth (um)', 'carrier density (cm-3)', 'acceptor density (cm-3)'])
eval_data_filename = raw_data_filename + '_eval.txt'
eval_data_path = os.path.join('..', 'data', 'eval', eval_data_filename)

eval_data.to_csv(eval_data_path, sep='\t', float_format='%1.6E', index=False)

plt.plot(scr_width / 1e-4, carrier_density, scr_width / 1e-4, acceptor_density / 1e6)
plt.yscale('log')
plt.show()