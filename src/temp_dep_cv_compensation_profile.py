import numpy as np
import matplotlib.pyplot as plt
import os

path = os.path.join('..', 'data', 'B_as-rec_CV_RT_dark.txt')
reverse_voltage, capacitance, scr_width, carrier_density, conductivity = np.loadtxt(path, unpack=True, skiprows=2)

plt.plot(scr_width, carrier_density)
plt.show()