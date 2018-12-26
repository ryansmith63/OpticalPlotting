import matplotlib.pyplot as plt
import numpy as np
from file_reader import Run, get_independent_variables_and_relative_intensities
from plotting import plot_runs, plot_TSTR_fit
from TSTR_fit_new import fit_parameters, fit_parameters_and_angle, BRIDF_plotter, reflectance_diffuse, reflectance_specular, BRIDF
import time

path = "vuv_height_comparison_and_first_data\\M18 turn\\center of sample\\"#\\M18 turn\\center of sample\\"
#path = "vuv_height_comparison_and_first_data/M18 turn polished/center of sample/"
specpath = "vuv_height_comparison_and_first_data\\Spectralon\\Centered, 500 nm\\"
# m18pol30 = Run(path + "2018_08_30__11_23_42.txt") 
# m18pol45 = Run(path + "2018_08_30__11_29_08.txt") 
# m18pol60 = Run(path + "2018_08_30__11_34_35.txt") 
# m18pol75 = Run(path + "2018_08_30__11_39_45.txt") 
# nxt8530 = Run(path + "2018_08_30__14_13_18.txt") 
# nxt8545 = Run(path + "2018_08_30__14_19_14.txt") 
# nxt8560 = Run(path + "2018_08_30__14_24_51.txt") 
# nxt8575 = Run(path + "2018_08_30__14_30_46.txt")
# skived30 = Run(path + "2018_08_29__17_29_38.txt") 
# skived45 = Run(path + "2018_08_29__17_35_12.txt") 
# skived60 = Run(path + "2018_08_29__17_40_57.txt") 
# skived75 = Run(path + "2018_08_29__17_46_39.txt") 
m1830 = Run(path + "2018_08_30__15_16_18.txt") # Fit to 36 deg
m1830.change_theta_i(33)
m1845 = Run(path + "2018_08_30__15_22_02.txt") # Fit to 46 deg
m1860 = Run(path + "2018_08_30__15_29_09.txt") # Fit to 60 deg
m1860.change_theta_i(61)
m1875 = Run(path + "2018_08_30__15_35_27.txt") # Fit to 85 deg
m1875.change_theta_i(76)
# m1730 = Run(path + "2018_08_30__16_18_57.txt")
# m1745 = Run(path + "2018_08_30__16_11_27.txt")
# m1760 = Run(path + "2018_08_30__16_06_11.txt")
# m1775 = Run(path + "2018_08_30__15_59_33.txt")
# nx30 = Run(path + "2018_08_29__15_57_32.txt")
# nx45 = Run(path + "2018_08_29__16_02_55.txt")
# nx60 = Run(path + "2018_08_29__16_08_16.txt")
# nx75 = Run(path + "2018_08_29__16_13_39.txt")
# spec8 = Run(specpath + "2018_09_04__16_38_49.txt")
# spec30 = Run(specpath + "2018_09_04__16_45_48.txt")
# spec45 = Run(specpath + "2018_09_04__16_51_29.txt")
# spec60 = Run(specpath + "2018_09_04__16_56_58.txt")
#spec8 = Run(specpath + "2018_09_04__16_38_49.txt")
#spec30 = Run(specpath + "2018_09_04__11_31_30.txt")
#spec35 = Run(specpath + "2018_10_02__16_55_53.txt")
#spec45 = Run(specpath + "2018_09_04__11_39_44.txt")
#spec50 = Run(specpath + "2018_10_02__17_01_27.txt")
#spec60 = Run(specpath + "2018_09_04__11_45_54.txt")

runs = [m1830,m1845,m1860,m1875]#[m18pol30, m18pol45, m18pol60, m18pol75]#[m1875]#
labels=["30 degrees","45 degrees","60 degrees","75 degrees"]#["30 degrees", "45 degrees", "60 degrees", "75 degrees"]#["75 degrees"]#

# Plot BRIDF data
sample_name="M18 turn"
plot_runs(runs, title=sample_name+", 175 nm", log=True, labels=labels, errorbars=True,legend_loc=6)
t0=time.time()

# Fit data
fit_params = fit_parameters(get_independent_variables_and_relative_intensities(runs),p0=[0.75,1.57,0.2],average_angle=4., precision=-1,use_errs=True,use_spike=False)
#fit_params_ang = fit_parameters_and_angle(get_independent_variables_and_relative_intensities(runs),average_angle=4.)
#fit_ang = fit_params_ang[0]
#fit_params = fit_params_ang[1:]
print("Fit parameters (rho_L, n, gamma): "+str(fit_params))
#print("Fit angle: "+str(fit_ang))
t1=time.time()
print("Fitting time: {0}".format(t1-t0))

# Plot BRIDF model from fits
#plot_TSTR_fit(8., 1., fit_params, label="fitted", color="k", average_angle=4., precision=-1)
plot_TSTR_fit(33., 1., fit_params, color="k", average_angle=4., precision=-1)
plot_TSTR_fit(45., 1., fit_params, color="k", average_angle=4., precision=-1)
plot_TSTR_fit(61., 1., fit_params, color="k", average_angle=4., precision=-1)
plot_TSTR_fit(76, 1., fit_params, color="k", average_angle=4., precision=-1)
#plt.text(0.1,0.1,r"Fit: $\rho_L$={0:.3f}, n={1:.2f}, $\gamma$={2:.3f}".format(*fit_params),transform=plt.gca().transAxes,fontsize=13)
#plt.text(0.1,0.2,r"Fit: $\rho_L$={0:.3f}, n={1:.2f}, $\gamma$={2:.3f}, K={3:.3f}".format(*fit_params),transform=plt.gca().transAxes,fontsize=13)
t2=time.time()
print("Plotting time: {0}".format(t2-t1))

# Now calculate hemispherical reflectance
# plt.figure()

# # Different incident angles to calculate for
# x = [8,30,45,60]#[30, 45, 60, 75]#[0,10,20,30, 45, 55, 60, 65, 70, 75, 80, 85]#

# y_diffuse = [reflectance_diffuse(theta, 1., 0.5, fit_params) for theta in x]
# y_specular = [reflectance_specular(theta, 1., 0.5, fit_params) for theta in x]
# y_total = [y_diffuse[i] + y_specular[i] for i in range(len(y_specular))]

# plt.plot(x, y_diffuse, label="diffuse")
# plt.plot(x, y_specular, label="specular")
# plt.plot(x, y_total, label="total")

# plt.xlabel("incident angle (degrees)")
# plt.ylabel("reflectance (fraction)")
# plt.legend()

# plt.title("Fitted "+sample_name+" Reflectance, 500 nm")
# t3=time.time()
# print("Reflectance calc time: {0}".format(t3-t2))
plt.show()
