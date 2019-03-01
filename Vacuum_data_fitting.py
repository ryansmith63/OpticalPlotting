import matplotlib.pyplot as plt
import numpy as np
from file_reader import Run, get_independent_variables_and_relative_intensities
from plotting import plot_runs, plot_TSTR_fit
from TSTR_fit_new import fit_parameters, fit_parameters_grid, fit_parameters_and_angle, BRIDF_plotter, reflectance_diffuse, reflectance_specular, BRIDF
import time

#path = "vuv_height_comparison_and_first_data\\M18 turn\\center of sample\\"#\\M18 turn\\center of sample\\"
path = "Vacuum measurements after 3rd xenon run/Jan 9-12/M17 turn/Blue height/500nm/"
#specpath = "vuv_height_comparison_and_first_data\\Spectralon\\Centered, 500 nm\\"

###August m18 polished 175nm  vuv_height_comparison_and_first_data/M18 turn polished/center of sample/
#m18pol30 = Run(path + "2018_08_30__11_23_42.txt") 
#m18pol30.change_theta_i(29)
#m18pol45 = Run(path + "2018_08_30__11_29_08.txt") 
#m18pol45.change_theta_i(43)
#m18pol60 = Run(path + "2018_08_30__11_34_35.txt") 
#m18pol75 = Run(path + "2018_08_30__11_39_45.txt") 
#m18pol75.change_theta_i(76)
### M18 polished October 1  vuv_height_comparison_and_first_data/M18 turn polished/center of sample repeat oct 1/
#m18pol30 = Run(path + "2018_10_01__14_06_14.txt")
#m18pol45 = Run(path + "2018_10_01__13_59_27.txt")
#m18pol45.change_theta_i(46)
#m18pol60 = Run(path + "2018_10_01__13_52_32.txt")
#m18pol75 = Run(path + "2018_10_01__12_09_58.txt")
### M18 polished 2/20/19 178nm  Vacuum measurements 2-20-19/M18 polished/178nm/
#m18pol30 = Run(path + "2019_02_20__16_13_06.txt")
#m18pol30.change_theta_i(28.5)
#m18pol45 = Run(path + "2019_02_20__16_07_47.txt")
#m18pol45.change_theta_i(43)
#m18pol52 = Run(path + "2019_02_20__16_00_46.txt")
#m18pol52.change_theta_i(50)
#m18pol60 = Run(path + "2019_02_20__15_55_37.txt")
#m18pol60.change_theta_i(58.5)
#m18pol67 = Run(path + "2019_02_20__15_50_27.txt")
#m18pol67.change_theta_i(65.5)
#m18pol75 = Run(path + "2019_02_20__15_44_29.txt")
#m18pol75.change_theta_i(73.5)
### M18 polished 2/20/19 178nm  Vacuum measurements 2-20-19/M18 polished/165nm/
#m18pol30 = Run(path + "2019_02_20__16_52_02.txt")
#m18pol30.change_theta_i(29)
#m18pol45 = Run(path + "2019_02_20__16_46_30.txt")
#m18pol45.change_theta_i(44)
#m18pol52 = Run(path + "2019_02_20__16_41_01.txt")
#m18pol52.change_theta_i(50.5)
#m18pol60 = Run(path + "2019_02_20__16_35_46.txt")
#m18pol60.change_theta_i(58)
#m18pol67 = Run(path + "2019_02_20__16_30_19.txt")
#m18pol67.change_theta_i(66.5)
#m18pol75 = Run(path + "2019_02_20__16_25_21.txt")
#m18pol75.change_theta_i(74.5)

###August NXT85 175nm  vuv_height_comparison_and_first_data/NXT85 turn/center of sample/
#nxt8530 = Run(path + "2018_08_30__14_13_18.txt") 
#nxt8545 = Run(path + "2018_08_30__14_19_14.txt") 
#nxt8560 = Run(path + "2018_08_30__14_24_51.txt") 
#nxt8575 = Run(path + "2018_08_30__14_30_46.txt")
###Before 2nd Xe run 175nm NXT85  2nd Xenon Run Measurements/Vacuum Test before Condensing, plus mirror tests in air-vacuum/Sample 5 measurements/
#nxt8530 = Run(path + "2018_11_26__17_14_45.txt")
#nxt8545 = Run(path + "2018_11_26__17_09_57.txt")
#nxt8552 = Run(path + "2018_11_26__17_05_19.txt")
#nxt8560 = Run(path + "2018_11_26__17_00_30.txt")
#nxt8567 = Run(path + "2018_11_26__16_55_51.txt")
#nxt8575 = Run(path + "2018_11_26__16_50_58.txt")

### August LUX 175nm  vuv_height_comparison_and_first_data/old sample/Center of sample
#lux30 = Run(path + "2018_08_29__15_34_41.txt")
#lux45 = Run(path + "2018_08_29__15_25_14.txt")
#lux60 = Run(path + "2018_08_29__15_15_22.txt")
#lux75 = Run(path + "2018_08_29__15_04_19.txt")
### LUX after 1st run  First Xe Run Measurements/follow-up vacuum measurements/Sample 1/
#lux30 = Run(path + "2018_11_09__14_14_55.txt")
#lux30.change_theta_i(27)
#lux45 = Run(path + "2018_11_09__14_08_54.txt")
#lux45.change_theta_i(42)
#lux52 = Run(path + "2018_11_09__14_02_45.txt")
#lux52.change_theta_i(49)
#lux60 = Run(path + "2018_11_09__13_55_54.txt")
#lux60.change_theta_i(58)
#lux67 = Run(path + "2018_11_09__13_47_53.txt")
#lux67.change_theta_i(65)
#lux75 = Run(path + "2018_11_09__13_41_18.txt")
#lux75.change_theta_i(73)
### 1/14/19 LUX 400nm  Vacuum measurements after 3rd xenon run/Jan 14/LUX/400nm/
#lux30 = Run(path + "2019_01_14__15_11_37.txt")
#lux45 = Run(path + "2019_01_14__15_06_53.txt")
#lux52 = Run(path + "2019_01_14__15_02_04.txt")
#lux60 = Run(path + "2019_01_14__14_56_59.txt")
#lux60.change_theta_i(59)
#lux67 = Run(path + "2019_01_14__14_52_14.txt")
#lux67.change_theta_i(66)
#lux75 = Run(path + "2019_01_14__14_44_03.txt")
#lux75.change_theta_i(72.5)
### 1/14/19 LUX 500nm  Vacuum measurements after 3rd xenon run/Jan 14/LUX/500nm/
#lux30 = Run(path + "2019_01_14__17_33_40.txt")
#lux30.change_theta_i(29)
#lux45 = Run(path + "2019_01_14__17_29_06.txt")
#lux45.change_theta_i(44)
#lux52 = Run(path + "2019_01_14__17_24_26.txt")
#lux52.change_theta_i(51)
#lux60 = Run(path + "2019_01_14__17_19_41.txt")
#lux60.change_theta_i(59)
#lux67 = Run(path + "2019_01_14__17_14_52.txt")
#lux67.change_theta_i(66)
#lux75 = Run(path + "2019_01_14__17_09_39.txt")
#lux75.change_theta_i(73)

### August LZ Skived 175nm  vuv_height_comparison_and_first_data/Skived LZ/center of sample/
#skived30 = Run(path + "2018_08_29__17_29_38.txt") 
#skived30.change_theta_i(29)
#skived45 = Run(path + "2018_08_29__17_35_12.txt") 
#skived45.change_theta_i(43)
#skived60 = Run(path + "2018_08_29__17_40_57.txt") 
#skived60.change_theta_i(58)
#skived75 = Run(path + "2018_08_29__17_46_39.txt") 
#skived75.change_theta_i(74)
###October 3 LZ Skived 255nm  vuv_height_comparison_and_first_data/Skived LZ/center of sample/255nm/
#skived30 = Run(path + "2018_10_03__13_18_54.txt")
#skived30.change_theta_i(28)
#skived45 = Run(path + "2018_10_03__13_24_45.txt")
#skived45.change_theta_i(42)
#skived60 = Run(path + "2018_10_03__13_30_19.txt")
#skived60.change_theta_i(56.5)
#skived75 = Run(path + "2018_10_03__13_36_22.txt")
#skived75.change_theta_i(71)
###October 3 LZ Skived 310 nm  vuv_height_comparison_and_first_data/Skived LZ/center of sample/310nm/
#skived30 = Run(path + "2018_10_03__11_18_57.txt")
#skived30.change_theta_i(28)
#skived45 = Run(path + "2018_10_03__11_23_54.txt")
#skived45.change_theta_i(43)
#skived60 = Run(path + "2018_10_03__11_28_48.txt")
#skived60.change_theta_i(58)
#skived75 = Run(path + "2018_10_03__12_57_41.txt")
#skived75.change_theta_i(73)
###October 4 LZ Skived 500nm  vuv_height_comparison_and_first_data/Skived LZ/center of sample/500nm/
#skived30 = Run(path + "2018_10_04__16_35_27.txt")
#skived30.change_theta_i(32)
#skived45 = Run(path + "2018_10_04__16_41_08.txt")
#skived45.change_theta_i(47)
#skived60 = Run(path + "2018_10_04__16_50_12.txt")
#skived60.change_theta_i(63)
#skived75 = Run(path + "2018_10_04__16_55_04.txt")
#skived75.change_theta_i(77)
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
### After 3rd run LZ Skived 178nm Vacuum measurements after 3rd xenon run/Jan 9-12/LZ Skived/Blue height/178nm/
#skived30 = Run(path + "2019_01_10__10_52_19.txt")
#skived30.change_theta_i(28)
#skived45 = Run(path + "2019_01_10__10_45_18.txt")
#skived52 = Run(path + "2019_01_10__10_39_26.txt")
#skived60 = Run(path + "2019_01_10__10_34_19.txt")
#skived67 = Run(path + "2019_01_10__10_29_11.txt")
#skived75 = Run(path + "2019_01_10__10_23_22.txt")
### 1/14/19 LZ Skived 400nm  Vacuum measurements after 3rd xenon run/Jan 14/LZ Skived/400nm/
#skived30 = Run(path + "2019_01_14__15_43_06.txt")
#skived30.change_theta_i(28)
#skived45 = Run(path + "2019_01_14__15_38_26.txt")
#skived45.change_theta_i(43)
#skived52 = Run(path + "2019_01_14__15_33_45.txt")
#skived52.change_theta_i(51)
#skived60 = Run(path + "2019_01_14__15_29_06.txt")
#skived60.change_theta_i(59)
#skived67 = Run(path + "2019_01_14__15_24_17.txt")
#skived67.change_theta_i(66)
#skived75 = Run(path + "2019_01_14__15_19_30.txt")
#skived75.change_theta_i(73)
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

###August M18 turn  175nm  vuv_height_comparison_and_first_data/M18 turn/center of sample/
#m1830 = Run(path + "2018_08_30__15_16_18.txt")
#m1830.change_theta_i(32)
#m1845 = Run(path + "2018_08_30__15_22_02.txt") 
#m1860 = Run(path + "2018_08_30__15_29_09.txt") 
#m1860.change_theta_i(61)
#m1875 = Run(path + "2018_08_30__15_35_27.txt") 
#m1875.change_theta_i(77)
### M18 turn after first run  First Xe Run Measurements/follow-up vacuum measurements/Sample 3/
#m1830 = Run(path + "2018_11_09__12_13_51.txt")
#m1830.change_theta_i(28)
#m1845 = Run(path + "2018_11_09__12_04_53.txt")
#m1845.change_theta_i(44)
#m1852 = Run(path + "2018_11_09__12_00_07.txt")
#m1852.change_theta_i(51)
#m1860 = Run(path + "2018_11_09__11_54_58.txt")
#m1867 = Run(path + "2018_11_09__11_49_54.txt")
#m1867.change_theta_i(66)
#m1875 = Run(path + "2018_11_09__11_44_37.txt")
#m1875.change_theta_i(73)
### M18 turn 2/20/19  Vacuum measurements 2-20-19/M18 turn/
#m1830 = Run(path + "2019_02_20__15_37_40.txt")
#m1845 = Run(path + "2019_02_20__15_32_29.txt")
#m1852 = Run(path + "2019_02_20__15_26_58.txt")
#m1860 = Run(path + "2019_02_20__15_21_51.txt")
#m1867 = Run(path + "2019_02_20__15_16_34.txt")
#m1875 = Run(path + "2019_02_20__15_11_22.txt")

###August M17 175nm  vuv_height_comparison_and_first_data/M17 turn/Centered on sample/
#m1730 = Run(path + "2018_08_30__16_18_57.txt")
#m1745 = Run(path + "2018_08_30__16_11_27.txt")
#m1760 = Run(path + "2018_08_30__16_06_11.txt")
#m1760.change_theta_i(62)
#m1775 = Run(path + "2018_08_30__15_59_33.txt")
#m1775.change_theta_i(76)
###September 21 no cell M17 175nm  vuv_height_comparison_and_first_data/M17 turn/Centered on sample, no cell/
#m1730 = Run(path + "2018_09_21__17_50_04.txt")
#m1745 = Run(path + "2018_09_21__17_44_58.txt")
#m1760 = Run(path + "2018_09_21__17_32_43.txt")
#m1775 = Run(path + "2018_09_21__17_40_00.txt")
###After 3rd run M17 178nm  Vacuum measurements after 3rd xenon run/Jan 9-12/M17 turn/Blue height/178nm
#m1730 = Run(path + "2019_01_10__12_00_27.txt")
#m1745 = Run(path + "2019_01_10__11_55_45.txt")
#m1752 = Run(path + "2019_01_10__11_50_15.txt")
#m1760 = Run(path + "2019_01_10__11_45_29.txt")
#m1767 = Run(path + "2019_01_10__11_40_29.txt")
#m1775 = Run(path + "2019_01_10__11_35_25.txt")
#m1775.change_theta_i(74)
###After 3rd run M17 400nm  Vacuum measurements after 3rd xenon run/Jan 9-12/M17 turn/Blue height/400nm/
#m1730 = Run(path + "2019_01_10__16_39_24.txt")
#m1745 = Run(path + "2019_01_10__16_34_37.txt")
#m1752 = Run(path + "2019_01_10__16_29_15.txt")
#m1760 = Run(path + "2019_01_10__16_23_48.txt")
#m1767 = Run(path + "2019_01_10__16_13_08.txt")
#m1775 = Run(path + "2019_01_10__16_08_12.txt")
###After 3rd run M17 500nm  Vacuum measurements after 3rd xenon run/Jan 9-12/M17 turn/Blue height/500nm/
m1730 = Run(path + "2019_01_10__17_10_44.txt")
m1745 = Run(path + "2019_01_10__17_05_52.txt")
m1752 = Run(path + "2019_01_10__17_01_07.txt")
m1760 = Run(path + "2019_01_10__16_56_14.txt")
m1767 = Run(path + "2019_01_10__16_51_34.txt")
m1775 = Run(path + "2019_01_10__16_46_53.txt")

###August 807NX 175nm  vuv_height_comparison_and_first_data/807NX turn/Center of sample/
#nx30 = Run(path + "2018_08_29__15_57_32.txt")
#nx45 = Run(path + "2018_08_29__16_02_55.txt")
#nx60 = Run(path + "2018_08_29__16_08_16.txt")
#nx75 = Run(path + "2018_08_29__16_13_39.txt")
###807NX 2/20/19  Vacuum measurements 2-20-19/807NX/
#nx30 = Run(path + "2019_02_20__15_02_12.txt")
#nx30.change_theta_i(29)
#nx45 = Run(path + "2019_02_20__14_56_50.txt")
#nx45.change_theta_i(44)
#nx52 = Run(path + "2019_02_20__14_51_22.txt")
#nx52.change_theta_i(51)
#nx60 = Run(path + "2019_02_20__14_44_13.txt")
#nx67 = Run(path + "2019_02_20__14_29_45.txt")
#nx75 = Run(path + "2019_02_20__14_24_17.txt")
#nx75.change_theta_i(74)

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
#spec8 = Run(path + "2019_01_14__14_38_24.txt")
#spec20 = Run(path + "2019_01_14__14_33_56.txt")
#spec30 = Run(path + "2019_01_14__14_29_01.txt")
#spec45 = Run(path + "2019_01_14__14_24_06.txt")
#spec52 = Run(path + "2019_01_14__14_18_57.txt")
#spec60 = Run(path + "2019_01_14__14_14_15.txt")
#spec67 = Run(path + "2019_01_14__14_09_32.txt")
#spec75 = Run(path + "2019_01_14__14_00_52.txt")
#spec75.change_theta_i(73)
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

runs = [m1730, m1745, m1752, m1760, m1767]#[m1875]#
labels= ["30 degrees","45 degrees","52 degrees","60 degrees","67 degrees"]

# Plot BRIDF data
sample_name="M17 turn"
plot_runs(runs, title=sample_name+", 500nm", log=True, labels=labels, errorbars=True)
t0=time.time()

# Fit data
fit_params = fit_parameters(get_independent_variables_and_relative_intensities(runs),p0=[0.75,1.57,0.05],average_angle=4, precision=0.25,sigma_theta_i =2, use_errs=True,use_spike=False)
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
plot_TSTR_fit(30., 1., fit_params, color="r", average_angle=4, sigma_theta_i =2, precision=0.25)
plot_TSTR_fit(45., 1., fit_params, color="g", average_angle=4, sigma_theta_i =2, precision=0.25)
plot_TSTR_fit(52., 1., fit_params, color="b", average_angle=4, sigma_theta_i =2, precision=0.25)
plot_TSTR_fit(60., 1., fit_params, color="m", average_angle=4, sigma_theta_i =2, precision=0.25)
plot_TSTR_fit(67., 1., fit_params, color="c", average_angle=4, sigma_theta_i =2, precision=0.25)
#plot_TSTR_fit(75., 1., fit_params, color="y", average_angle=4, sigma_theta_i =2, precision=0.25)
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
