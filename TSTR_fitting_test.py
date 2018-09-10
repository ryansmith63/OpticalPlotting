import matplotlib.pyplot as plt
import numpy as np
from file_reader import Run, get_independent_variables_and_relative_intensities
from plotting import plot_runs, plot_TSTR_fit
from TSTR_fit_new import fit_parameters, BRIDF_plotter, reflectance_diffuse, reflectance_specular

m18path = "vuv_height_comparison_and_first_data/M18 turn/center of sample/"
m1830 = Run(m18path + "2018_08_30__15_16_18.txt")
m1845 = Run(m18path + "2018_08_30__15_22_02.txt")
m1860 = Run(m18path + "2018_08_30__15_29_09.txt")
m1875 = Run(m18path + "2018_08_30__15_35_27.txt")

runs = [m1830, m1845, m1860, m1875]


plot_runs(runs, title="M18 turn, 175nm", log=True, labels=["30 degrees", "45 degrees", "60 degrees", "75 degrees"])

fit_params = fit_parameters(get_independent_variables_and_relative_intensities(runs))
plot_TSTR_fit(30., 1., fit_params, label="fitted", color="k")
plot_TSTR_fit(45., 1., fit_params, color="k")
plot_TSTR_fit(60., 1., fit_params, color="k")
plot_TSTR_fit(75., 1., fit_params, color="k")



plt.figure()

x = [30, 45, 60, 75]

y_diffuse = [reflectance_diffuse(theta, 1., 0.5, fit_params) for theta in x]
y_specular = [reflectance_specular(theta, 1., 0.5, fit_params) for theta in x]
y_total = [y_diffuse[i] + y_specular[i] for i in range(len(y_specular))]

plt.plot(x, y_diffuse, label="diffuse")
plt.plot(x, y_specular, label="specular")
plt.plot(x, y_total, label="total")

plt.xlabel("incident angle (degrees)")
plt.ylabel("reflectance (fraction)")
plt.legend()

plt.title("Fitted M18 Reflectance, 175 nm")

plt.show()
