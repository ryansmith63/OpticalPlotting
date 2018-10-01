import matplotlib.pyplot as plt
import numpy as np
from file_reader import Run, get_independent_variables_and_relative_intensities
from plotting import plot_runs, plot_TSTR_fit
from TSTR_fit_new import fit_parameters, fit_parameters_and_angle, BRIDF_plotter, reflectance_diffuse, reflectance_specular, BRIDF
import time

# Plots BRIDF at a few angles w/ and w/o specular spike, as well as hemispherical reflectances
# Used for testing effects of the specular spike model params

fit_params=[.74,1.45,0.049] # From Coimbra fit
fit_params_spec=[.74,1.45,0.049,1.7] # From Coimbra fit
plot_TSTR_fit(45., 1., fit_params, color="k", label="no spike", average_angle=4., precision=0.25)
plot_TSTR_fit(65., 1., fit_params, color="k", average_angle=4., precision=0.25)
plot_TSTR_fit(80., 1., fit_params, color="k", average_angle=4., precision=0.25)
plot_TSTR_fit(45., 1., fit_params_spec, color="r", label="with spike", average_angle=4., precision=0.25)
plot_TSTR_fit(65., 1., fit_params_spec, color="r", average_angle=4., precision=0.25)
plot_TSTR_fit(80., 1., fit_params_spec, color="r", average_angle=4., precision=0.25)
plt.text(0.1,0.2,r"Fit: $\rho_L$={0:.3f}, n={1:.2f}, $\gamma$={2:.3f}, K={3:.3f}".format(*fit_params_spec),transform=plt.gca().transAxes,fontsize=13)
plt.yscale("log")

plt.figure()

x = [0,10,20,30, 45, 55, 60, 65, 70, 75, 80, 85]#

#fit_params=[0.58,1.49,0.064] # Coimbra fitted parameters for skived sample from thesis
#fit_params=[0.52,1.51,0.057] # Coimbra fitted parameters for unpolished molded sample from paper
y_diffuse = [reflectance_diffuse(theta, 1., 0.5, fit_params) for theta in x]
y_specular = [reflectance_specular(theta, 1., 0.5, fit_params) for theta in x]
y_total = [y_diffuse[i] + y_specular[i] for i in range(len(y_specular))]

plt.plot(x, y_diffuse, label="diffuse")
plt.plot(x, y_specular, label="specular")
plt.plot(x, y_total, label="total")

y_specular_spike = [reflectance_specular(theta, 1., 0.5, fit_params_spec) for theta in x]
y_total_spike = [y_diffuse[i] + y_specular_spike[i] for i in range(len(y_specular))]
plt.plot(x, y_specular_spike, label="specular w/ spike")
plt.plot(x, y_total_spike, label="total w/ spike")

plt.xlabel("incident angle (degrees)")
plt.ylabel("reflectance (fraction)")
plt.legend()

plt.title("Reflectance from Coimbra specular spike fits, 175 nm")
plt.show()
