import matplotlib.pyplot as plt
import numpy as np
from file_reader import Run, get_independent_variables_and_relative_intensities
from plotting import plot_runs, plot_TSTR_fit, plot_double_TSTR_fit
from TSTR_fit_new import fit_parameters, fit_parameters_and_angle, BRIDF_plotter, reflectance_diffuse, reflectance_specular, BRIDF_specular_plotter, BRIDF_diffuse_plotter, BRIDF_plotter
import time

path = "First Xe Run Measurements/first measurements with no bubbles in cell 11-01-2/Sample 1/"

s1_30 = Run(path + "2018_11_02__18_28_43.txt")
s1_45 = Run(path + "2018_11_02__18_35_18.txt")
s1_52 = Run(path + "2018_11_02__18_41_20.txt")
s1_60 = Run(path + "2018_11_02__18_47_03.txt")
s1_67 = Run(path + "2018_11_02__18_51_55.txt")
s1_75 = Run(path + "2018_11_02__18_57_35.txt")
#runs = [s1_30, s1_45, s1_52, s1_60, s1_67, s1_75]
runs = [s1_30, s1_45, s1_60, s1_75]

#labels=["30 degrees", "45 degrees", "52 degrees", "60 degrees", "67 degrees", "75 degrees"]
labels=["30 degrees", "45 degrees", "60 degrees", "75 degrees"]

# Plot BRIDF data
sample_name="Sample 1 (LUX)"
#plot_runs(runs, title=sample_name+", 178nm", log=True, labels=labels, errorbars=True)
#plt.show()

# In a loop: ask for input params, plot model overlaid on data
get_input=True
params=[0.24,1.27,0.1,5.0,
        0.1,0.5]

plot_TSTR_fit(30., 1.69, params, label="total", color="k", average_angle=4., precision=0.25, sigma_theta_i=2.)

# these parts of the plot do not currently average over the aperture
theta_r_in_degrees_array = np.linspace(0., 85., 100)
#plt.plot(theta_r_in_degrees_array, BRIDF_plotter(theta_r_in_degrees_array, 0., 30., 1.69, 0.5, params[:5], average_angle=4, precision=0.25, sigma_theta_i=2.), label="total, averaged over aperture")
plt.plot(theta_r_in_degrees_array, BRIDF_plotter(theta_r_in_degrees_array, 0., 30., 1.69, 0.5, params[:5]), label="total, not averaged over aperture")
plt.plot(theta_r_in_degrees_array, BRIDF_specular_plotter(theta_r_in_degrees_array, 0., 30., 1.69, 0.5, params[:5]), label="specular")
plt.plot(theta_r_in_degrees_array, BRIDF_diffuse_plotter(theta_r_in_degrees_array, 0., 30., 1.69, 0.5, params[:5]), label="diffuse")

plt.legend()
plt.show()

