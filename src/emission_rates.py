import numpy as np
import math as mt
import matplotlib.pyplot as plt

number_points = 64
emission_rate_spectrum = np.zeros(number_points)
DeltaC = np.zeros(number_points)
time_points = np.zeros(number_points)
transient = np.zeros(number_points)
transient_grid = np.zeros((number_points, number_points))
min_rate = 1 / 640e-3
max_rate = 1 / 20e-3

scale_factor = (max_rate / min_rate) ** (1 / (number_points - 1))

for rate_point in range(number_points):
    emission_rate_spectrum[rate_point] = min_rate * scale_factor ** rate_point
    inv_point = number_points - rate_point - 1
    time_points[inv_point] = 1 / emission_rate_spectrum[rate_point]
    if emission_rate_spectrum[rate_point] >= max_rate:
        break

DeltaC[10] = 2e-12
DeltaC[32] = 4e-12

for time_point in range(number_points):
    transient_sum_over_rates = 0
    for rate_point in range(number_points):
        transient_grid[time_point, rate_point] = DeltaC[rate_point] * mt.exp(
            -emission_rate_spectrum[rate_point] * time_points[time_point])
        transient_sum_over_rates += transient_grid[time_point, rate_point]
    transient[time_point] = transient_sum_over_rates

plt.plot(time_points, transient)
plt.show()
