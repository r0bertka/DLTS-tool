import numpy as np
from scipy.signal import savgol_filter
from scipy.constants import elementary_charge as e
from scipy.constants import epsilon_0 as eps0
from scipy.constants import k as kB


def calculate_trap_profile(capacitance, pulse_voltage, DLTS_signal,
                           scr_depth, doping_density, epsilon, temperature, effective_density_of_states_cb,
                           trap_energy):
    scr_depth = scr_depth[:-1]
    doping_density = doping_density[:-1]
    bulk_doping_density = doping_density[-1]
    bulk_depth_of_scr = scr_depth[-1]
    bulk_capacitance = capacitance[-1]

    EF = - kB * temperature / e * np.log(bulk_doping_density / effective_density_of_states_cb)  # calculate Fermi level
    Lambda = np.sqrt(2 * eps0 * epsilon / (e * bulk_doping_density) * (trap_energy - EF))  # calculate Lambda
    depth_trap_response = scr_depth - Lambda  # calculate spatial coordinate of defect response

    smoothed_DLTS_signal = savgol_filter(DLTS_signal, 7, 3)
    delta_C_C = smoothed_DLTS_signal / (bulk_capacitance * 1e12)
    diff_delta_C_C = np.diff(delta_C_C) / np.diff(pulse_voltage)
    diff_delta_C_C = np.append(diff_delta_C_C, diff_delta_C_C[-1])
    N_trap = (e * bulk_depth_of_scr ** 2 * bulk_doping_density / (
                eps0 * epsilon)) * doping_density * diff_delta_C_C / 2.449 * scr_depth / depth_trap_response  # 2.449 is the correction factor for td = 5ms, ti = 640ms, lock-in correlation
    smoothed_trap_density = savgol_filter(N_trap, 7, 3)

    return depth_trap_response, smoothed_trap_density
