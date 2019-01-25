import matplotlib.pyplot as plt
import numpy as np
from file_reader import Run, get_independent_variables_and_relative_intensities
from plotting import plot_runs, plot_TSTR_fit
from TSTR_fit_new import fit_parameters, fit_parameters_new, fit_parameters_and_angle, BRIDF_plotter, reflectance_diffuse, reflectance_specular, BRIDF
import time

#path = "vuv_height_comparison_and_first_data\\M18 turn\\center of sample\\"#\\M18 turn\\center of sample\\"
path = "Vacuum measurements after 3rd xenon run/Jan 14/LZ Skived/400nm/"
#specpath = "vuv_height_comparison_and_first_data\\Spectralon\\Centered, 500 nm\\"

###August m18 polished 175nm
# m18pol30 = Run(path + "2018_08_30__11_23_42.txt") 
# m18pol45 = Run(path + "2018_08_30__11_29_08.txt") 
# m18pol60 = Run(path + "2018_08_30__11_34_35.txt") 
# m18pol75 = Run(path + "2018_08_30__11_39_45.txt") 

###August NXT85 175nm
# nxt8530 = Run(path + "2018_08_30__14_13_18.txt") 
# nxt8545 = Run(path + "2018_08_30__14_19_14.txt") 
# nxt8560 = Run(path + "2018_08_30__14_24_51.txt") 
# nxt8575 = Run(path + "2018_08_30__14_30_46.txt")
###Before 2nd Xe run 175nm NXT85 
#nxt8530 = Run(path + "2018_11_26__16_50_58.txt")
#nxt8545 = Run(path + "2018_11_26__16_55_51.txt")
#nxt8552 = Run(path + "2018_11_26__17_00_30.txt")
#nxt8560 = Run(path + "2018_11_26__17_05_19.txt")
#nxt8567 = Run(path + "2018_11_26__17_09_57.txt")
#nxt8575 = Run(path + "2018_11_26__17_14_45.txt")

### August LUX 178nm
### 1/14/19 LUX 400nm  Vacuum measurements after 3rd xenon run/Jan 14/LUX/400nm/
#lux30 = Run(path + "2019_01_14__15_11_37.txt")
#lux45 = Run(path + "2019_01_14__15_06_53.txt")
#lux52 = Run(path + "2019_01_14__15_02_04.txt")
#lux60 = Run(path + "2019_01_14__14_56_59.txt")
#lux67 = Run(path + "2019_01_14__14_52_14.txt")
#lux75 = Run(path + "2019_01_14__14_44_03.txt")
### 1/14/19 LUX 500nm  Vacuum measurements after 3rd xenon run/Jan 14/LUX/500nm/
#lux30 = Run(path + "2019_01_14__17_33_40.txt")
#lux45 = Run(path + "2019_01_14__17_29_06.txt")
#lux52 = Run(path + "2019_01_14__17_24_26.txt")
#lux60 = Run(path + "2019_01_14__17_19_41.txt")
#lux67 = Run(path + "2019_01_14__17_14_52.txt")
#lux75 = Run(path + "2019_01_14__17_09_39.txt")

### August LZ Skived 175nm  2nd Xenon Run Measurements/
#skived30 = Run(path + "2018_08_29__17_29_38.txt") 
#skived30.change_theta_i(29)
#skived45 = Run(path + "2018_08_29__17_35_12.txt") 
#skived45.change_theta_i(43)
#skived60 = Run(path + "2018_08_29__17_40_57.txt") 
#skived60.change_theta_i(58)
#skived75 = Run(path + "2018_08_29__17_46_39.txt") 
#skived75.change_theta_i(74)
###Following first Xe LZ Skived 175nm  First Xe Run Measurements/follow-up vacuum measurements/Sample9
#skived30 = Run(path + "2019_11_09__15_00_15.txt")
#skived30.change_theta_i(20)
#skived45 = Run(path + "2019_11_09__14_53_56.txt")
#skived45.change_theta_i(36)
#skived52 = Run(path + "2019_11_09__14_47_23.txt")
#skived52.change_theta_i(43)
#skived60 = Run(path + "2019_11_09__14_40_52.txt")
#skived60.change_theta_i(52)
#skived67 = Run(path + "2019_11_09__14_33_41.txt")
#skived67.change_theta_i(63)
#skived75 = Run(path + "2019_11_09__14_22_32.txt")
#skived75.change_theta_i(71)
### 1/14/19 LZ Skived 400nm  Vacuum measurements after 3rd xenon run/Jan 14/LZ Skived/400nm/
skived30 = Run(path + "2019_01_14__15_43_06.txt")
skived30.change_theta_i(29)
skived45 = Run(path + "2019_01_14__15_38_26.txt")
skived45.change_theta_i(44)
skived52 = Run(path + "2019_01_14__15_33_45.txt")
skived52.change_theta_i(51)
skived60 = Run(path + "2019_01_14__15_29_06.txt")
skived60.change_theta_i(59)
skived67 = Run(path + "2019_01_14__15_24_17.txt")
skived75 = Run(path + "2019_01_14__15_19_30.txt")
skived75.change_theta_i(74)
### 1/14/19 LZ Skived 500nm  Vacuum measurements after 3rd xenon run/Jan 14/LZ Skived/500nm/
#skived30 = Run(path + "2019_01_14__16_27_42.txt")
#skived30.change_theta_i(27)
#skived45 = Run(path + "2019_01_14__16_23_09.txt")
#skived45.change_theta_i(43)
#skived52 = Run(path + "2019_01_14__16_18_20.txt")
#skived52.change_theta_i(50)
#skived60 = Run(path + "2019_01_14__16_13_34.txt")
#skived60.change_theta_i(58)
#skived67 = Run(path + "2019_01_14__16_08_49.txt")
#skived67.change_theta_i(65.5)
#skived75 = Run(path + "2019_01_14__16_03_48.txt")
#skived75.change_theta_i(73)

###August M18 175nm
#m1830 = Run(path + "2018_08_30__15_16_18.txt") # Fit to 36 deg
#m1830.change_theta_i(33)
#m1845 = Run(path + "2018_08_30__15_22_02.txt") # Fit to 46 deg
#m1860 = Run(path + "2018_08_30__15_29_09.txt") # Fit to 60 deg
#m1860.change_theta_i(61)
#m1875 = Run(path + "2018_08_30__15_35_27.txt") # Fit to 85 deg
#m1875.change_theta_i(76)

###August M17 175nm
#m1730 = Run(path + "2018_08_30__16_18_57.txt")
#m1745 = Run(path + "2018_08_30__16_11_27.txt")
#m1760 = Run(path + "2018_08_30__16_06_11.txt")
#m1760.change_theta_i(62)
#m1775 = Run(path + "2018_08_30__15_59_33.txt")
#m1775.change_theta_i(76)

###August 807NX 175nm
# nx30 = Run(path + "2018_08_29__15_57_32.txt")
# nx45 = Run(path + "2018_08_29__16_02_55.txt")
# nx60 = Run(path + "2018_08_29__16_08_16.txt")
# nx75 = Run(path + "2018_08_29__16_13_39.txt")

###September Spectralon
# spec8 = Run(specpath + "2018_09_04__16_38_49.txt")
# spec30 = Run(specpath + "2018_09_04__16_45_48.txt")
# spec45 = Run(specpath + "2018_09_04__16_51_29.txt")
# spec60 = Run(specpath + "2018_09_04__16_56_58.txt")
###
#spec8 = Run(specpath + "2018_09_04__16_38_49.txt")
#spec30 = Run(specpath + "2018_09_04__11_31_30.txt")
#spec35 = Run(specpath + "2018_10_02__16_55_53.txt")
#spec45 = Run(specpath + "2018_09_04__11_39_44.txt")
#spec50 = Run(specpath + "2018_10_02__17_01_27.txt")
#spec60 = Run(specpath + "2018_09_04__11_45_54.txt")
### Spectralon 1/14/14 400nm
### Spectralon 1/14/14 500nm
#spec8 = Run(path + "2019_01_14__17_03_36.txt")
#spec20 = Run(path + "2019_01_14__16_59_12.txt")
#spec30 = Run(path + "2019_01_14__16_54_11.txt")
#spec45 = Run(path + "2019_01_14__16_49_27.txt")
#spec60 = Run(path + "2019_01_14__16_44_45.txt")
#spec60.change_theta_i(59)
#spec67 = Run(path + "2019_01_14__16_40_00.txt")
#spec67.change_theta_i(66)
#spec75 = Run(path + "2019_01_14__16_35_18.txt")
#spec75.change_theta_i(74)

runs = [skived30,skived45,skived52,skived60,skived67,skived75]#[m18pol30, m18pol45, m18pol60, m18pol75]#[m1875]#
labels= ["30 degrees","45 degrees","52 degrees","60 degrees","67 degrees","75 degrees"]

# Plot BRIDF data
sample_name="LZ Skived"
plot_runs(runs, title=sample_name+", 400 nm", log=False, labels=labels, errorbars=False)
t0=time.time()

# Fit data
fit_params = fit_parameters(get_independent_variables_and_relative_intensities(runs),p0=[0.75,1.57,0.05],average_angle=4, precision=0.25,use_errs=True,use_spike=False)
#fit_params = fit_parameters_new(get_independent_variables_and_relative_intensities(runs)[0], get_independent_variables_and_relative_intensities(runs)[1], get_independent_variables_and_relative_intensities(runs)[2], 0, 1, 30, 1, 2, 30, .01, .3, 30, plot=False, show=False)
#fit_params_ang = fit_parameters_and_angle(get_independent_variables_and_relative_intensities(runs),average_angle=4.)
#fit_ang = fit_params_ang[0]
#fit_params = fit_params_ang[1:]
print("Fit parameters (rho_L, n, gamma): "+str(fit_params))
#print("Fit angle: "+str(fit_ang))
t1=time.time()
print("Fitting time: {0}".format(t1-t0))

# Plot BRIDF model from fits
#plot_TSTR_fit(8., 1., fit_params, color="r", average_angle=4., precision=-1)
#plot_TSTR_fit(20., 1., fit_params, color="g", average_angle=4., precision=-1)
plot_TSTR_fit(29., 1., fit_params, color="r", average_angle=4, sigma_theta_i =2, precision=0.25)
plot_TSTR_fit(44., 1., fit_params, color="g", average_angle=4, sigma_theta_i =2, precision=0.25)
plot_TSTR_fit(51., 1., fit_params, color="b", average_angle=4, sigma_theta_i =2, precision=0.25)
plot_TSTR_fit(59, 1., fit_params, color="m", average_angle=4, sigma_theta_i =2, precision=0.25)
plot_TSTR_fit(67, 1., fit_params, color="c", average_angle=4, sigma_theta_i =2, precision=0.25)
plot_TSTR_fit(74, 1., fit_params, color="y", average_angle=4, sigma_theta_i =2, precision=0.25)
#plt.text(0.1,0.1,r"Fit: $\rho_L$={0:.3f}, n={1:.2f}, $\gamma$={2:.3f}".format(*fit_params),transform=plt.gca().transAxes,fontsize=13)
#plt.text(0.1,0.2,r"Fit: $\rho_L$={0:.3f}, n={1:.2f}, $\gamma$={2:.3f}, K={3:.3f}".format(*fit_params),transform=plt.gca().transAxes,fontsize=13)
t2=time.time()
print("Plotting time: {0}".format(t2-t1))

# Now calculate hemispherical reflectance
plt.figure()

# # Different incident angles to calculate for
x = [8,30,45,60]#[0,10,20,30, 45, 55, 60, 65, 70, 75, 80, 85] #[30, 45, 60, 75]

y_diffuse = [reflectance_diffuse(theta, 1., 0.5, fit_params) for theta in x]
y_specular = [reflectance_specular(theta, 1., 0.5, fit_params) for theta in x]
y_total = [y_diffuse[i] + y_specular[i] for i in range(len(y_specular))]

plt.plot(x, y_diffuse, label="diffuse")
plt.plot(x, y_specular, label="specular")
plt.plot(x, y_total, label="total")

plt.xlabel("incident angle (degrees)")
plt.ylabel("reflectance (fraction)")
plt.legend()

plt.title("Fitted "+sample_name+" Reflectance")
t3=time.time()
print("Reflectance calc time: {0}".format(t3-t2))
plt.show()
