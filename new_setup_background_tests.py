import numpy as np
import matplotlib.pyplot as plt
from file_reader import plot_runs, Run

path = "new setup background tests/"

no_walls = Run(path + "2018_05_25__13_44_05.txt")
same_side = Run(path + "2018_05_25__13_57_52.txt")
opposite_side = Run(path + "2018_05_25__14_06_33.txt")

min_oil = Run(path + "2018_05_29__16_22_24.txt")

runs = [no_walls, same_side, opposite_side, min_oil]

plot_runs(runs, labels=["no walls in air", "walls opposite laser and detector in air", "walls opposite laser and on same side as detector in air", "new setup in mineral oil"], show=False, rot=True)

# copied from 10_16_background_tests.py, using data from background_measurements/background in mineral oil no sample old walls new walls lens tubes.txt
distance_from_sample_to_photodiode = 5.435
photodiode_radius = (9 / 2.0) / 25.4
photodiode_solid_angle = np.pi * np.power(photodiode_radius, 2) / np.power(distance_from_sample_to_photodiode, 2)
flux_i = 0.005940 * 100e-6
sensitivity = 100 * 1e-9
intensity_factor = sensitivity / (photodiode_solid_angle * flux_i)

total_data = np.loadtxt("old_setup_data/background in mineral oil no sample old walls new walls lens tubes.txt", skiprows=1)
total_data_x = [180 - 90 - np.round(d[0], 2) for d in total_data]
total_data_y = [intensity_factor * d[1] for d in total_data]

old_data_x = total_data_x[0:302]
old_data_y = total_data_y[0:302]
plt.scatter(old_data_x, old_data_y, marker="x", c="c", s=5, label="old setup in mineral oil")

plt.legend()
plt.xlim(0, 170)

plt.show()
