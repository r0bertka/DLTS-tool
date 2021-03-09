import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.signal import savgol_filter
from scipy.constants import elementary_charge as e
from scipy.constants import epsilon_0 as eps0


def calculate_doping_profile(voltage, capacitance, eps, area):
    w = area * eps0 * eps / capacitance
    inverse_square_of_cap = 1 / capacitance ** 2
    diff_inverse_square_of_cap = np.diff(inverse_square_of_cap) / np.diff(voltage)
    diff_inverse_square_of_cap = np.append(diff_inverse_square_of_cap, diff_inverse_square_of_cap[-1])
    smooth_diff_inverse_square_of_cap = savgol_filter(diff_inverse_square_of_cap, 9, 3)
    n_smooth = -2 / (e * area ** 2 * eps0 * eps * smooth_diff_inverse_square_of_cap)

    return w, n_smooth
