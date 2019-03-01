import matplotlib.pyplot as plt
import numpy as np
from file_reader import Run, get_independent_variables_and_relative_intensities
from plotting import plot_runs, plot_TSTR_fit
from TSTR_fit_new import fit_parameters, fit_parameters_and_angle, BRIDF_plotter, reflectance_diffuse, reflectance_specular, BRIDF
import time

path_run1 = "First Xe Run Measurements//first measurements with no bubbles in cell 11-01-2//"
path = "2nd Xenon Run Measurements//"
#path = "vuv_height_comparison_and_first_data/M18 turn polished/center of sample/"
# m18pol30 = Run(path + "2018_08_30__11_23_42.txt") 
# m18pol45 = Run(path + "2018_08_30__11_29_08.txt") 
# m18pol60 = Run(path + "2018_08_30__11_34_35.txt") 
# m18pol75 = Run(path + "2018_08_30__11_39_45.txt") 
# nxt8530 = Run(path + "2018_08_30__14_13_18.txt") 
# nxt8545 = Run(path + "2018_08_30__14_19_14.txt") 
# nxt8560 = Run(path + "2018_08_30__14_24_51.txt") 
# nxt8575 = Run(path + "2018_08_30__14_30_46.txt")
s5_30 = Run(path + "Sample 5 with bubbles//2018_11_30__14_24_00.txt")
s5_45 = Run(path + "Sample 5 with bubbles//2018_11_30__14_29_18.txt")
s5_52 = Run(path + "Sample 5 with bubbles//2018_11_30__14_38_44.txt") # Two files exist for 52 deg; first probably got redone bc the angle range was too small
s5_60 = Run(path + "Sample 5 with bubbles//2018_11_30__14_43_44.txt") 
s5_67 = Run(path + "Sample 5 with bubbles//2018_11_30__14_48_38.txt")
s5_75 = Run(path + "Sample 5 with bubbles//2018_11_30__14_54_19.txt")
s8_30 = Run(path + "Sample 8 no bubbles//2018_12_03__14_52_24.txt")
s8_45 = Run(path + "Sample 8 no bubbles//2018_12_03__14_46_48.txt")
s8_52 = Run(path + "Sample 8 no bubbles//2018_12_03__14_42_00.txt")
s8_60 = Run(path + "Sample 8 no bubbles//2018_12_03__14_57_12.txt") 
s8_67 = Run(path + "Sample 8 no bubbles//2018_12_03__15_01_51.txt")
s8_75 = Run(path + "Sample 8 no bubbles//2018_12_03__15_06_25.txt") 
s9_first_30 = Run(path_run1 + "Sample 9//2018_11_02__20_38_01.txt")
s9_first_45 = Run(path_run1 + "Sample 9//2018_11_02__20_32_34.txt")
s9_first_52 = Run(path_run1 + "Sample 9//2018_11_02__20_27_27.txt")
s9_first_60 = Run(path_run1 + "Sample 9//2018_11_02__20_22_13.txt")
s9_first_67 = Run(path_run1 + "Sample 9//2018_11_02__20_17_11.txt")
s9_first_75 = Run(path_run1 + "Sample 9//2018_11_02__20_11_44.txt")
s9_bubbles_30 = Run(path + "Sample 9 with bubbles//2018_11_30__15_26_05.txt")
s9_bubbles_45 = Run(path + "Sample 9 with bubbles//2018_11_30__15_21_25.txt")
s9_bubbles_52 = Run(path + "Sample 9 with bubbles//2018_11_30__15_16_32.txt")
s9_bubbles_60 = Run(path + "Sample 9 with bubbles//2018_11_30__15_11_46.txt") 
s9_bubbles_67 = Run(path + "Sample 9 with bubbles//2018_11_30__15_07_06.txt")
s9_bubbles_75 = Run(path + "Sample 9 with bubbles//2018_11_30__15_02_17.txt") 
s9_nobubbles_30 = Run(path + "Sample 9 no bubbles//Before getter//2018_11_30__17_05_41.txt")
s9_nobubbles_45 = Run(path + "Sample 9 no bubbles//Before getter//2018_11_30__17_01_01.txt")
s9_nobubbles_52 = Run(path + "Sample 9 no bubbles//Before getter//2018_11_30__16_56_06.txt")
s9_nobubbles_60 = Run(path + "Sample 9 no bubbles//Before getter//2018_11_30__16_51_00.txt") 
s9_nobubbles_67 = Run(path + "Sample 9 no bubbles//Before getter//2018_11_30__16_46_13.txt")
s9_nobubbles_75 = Run(path + "Sample 9 no bubbles//Before getter//2018_11_30__16_40_45.txt") 
s9_getter_30 = Run(path + "Sample 9 no bubbles//After getter//2018_12_03__11_13_32.txt")
s9_getter_45 = Run(path + "Sample 9 no bubbles//After getter//2018_12_03__11_08_40.txt")
s9_getter_52 = Run(path + "Sample 9 no bubbles//After getter//2018_12_03__11_03_50.txt")
s9_getter_60 = Run(path + "Sample 9 no bubbles//After getter//2018_12_03__10_58_47.txt") 
s9_getter_67 = Run(path + "Sample 9 no bubbles//After getter//2018_12_03__10_53_50.txt")
s9_getter_75 = Run(path + "Sample 9 no bubbles//After getter//2018_12_03__10_48_50.txt") 
s9_lowp_30 = Run(path + "Sample 9 lower pressure//2018_12_05__13_52_04.txt")
s9_lowp_45 = Run(path + "Sample 9 lower pressure//2018_12_05__13_47_09.txt")
s9_lowp_52 = Run(path + "Sample 9 lower pressure//2018_12_05__13_42_21.txt")
s9_lowp_60 = Run(path + "Sample 9 lower pressure//2018_12_05__13_37_31.txt") 
s9_lowp_67 = Run(path + "Sample 9 lower pressure//2018_12_05__13_32_47.txt")
s9_lowp_75 = Run(path + "Sample 9 lower pressure//2018_12_05__13_27_55.txt")
s9_lowp2_30 = Run(path + "Sample 9 lower pressure 2//2018_12_07__13_28_41.txt")
s9_lowp2_45 = Run(path + "Sample 9 lower pressure 2//2018_12_07__13_23_51.txt")
s9_lowp2_52 = Run(path + "Sample 9 lower pressure 2//2018_12_07__13_18_50.txt")
s9_lowp2_60 = Run(path + "Sample 9 lower pressure 2//2018_12_07__13_14_04.txt") 
s9_lowp2_67 = Run(path + "Sample 9 lower pressure 2//2018_12_07__13_09_15.txt")
s9_lowp2_75 = Run(path + "Sample 9 lower pressure 2//2018_12_07__13_04_19.txt")
s9_medp_30 = Run(path + "Sample 9 medium pressure//2018_12_05__15_56_28.txt")
s9_medp_45 = Run(path + "Sample 9 medium pressure//2018_12_05__15_51_07.txt")
s9_medp_52 = Run(path + "Sample 9 medium pressure//2018_12_05__15_44_43.txt")
s9_medp_60 = Run(path + "Sample 9 medium pressure//2018_12_05__15_39_51.txt") 
s9_medp_67 = Run(path + "Sample 9 medium pressure//2018_12_05__15_35_13.txt")
s9_medp_75 = Run(path + "Sample 9 medium pressure//2018_12_05__15_29_31.txt")  
s9_medp2_30 = Run(path + "Sample 9 medium pressure 2//2018_12_07__14_42_49.txt")
s9_medp2_45 = Run(path + "Sample 9 medium pressure 2//2018_12_07__14_36_41.txt")
s9_medp2_52 = Run(path + "Sample 9 medium pressure 2//2018_12_07__14_31_00.txt")
s9_medp2_60 = Run(path + "Sample 9 medium pressure 2//2018_12_07__14_25_38.txt") 
s9_medp2_67 = Run(path + "Sample 9 medium pressure 2//2018_12_07__14_20_48.txt")
s9_medp2_75 = Run(path + "Sample 9 medium pressure 2//2018_12_07__14_14_47.txt") 
s9_hip_30 = Run(path + "Sample 9 higher pressure//2018_12_05__17_33_42.txt")
s9_hip_45 = Run(path + "Sample 9 higher pressure//2018_12_05__17_28_40.txt")
s9_hip_52 = Run(path + "Sample 9 higher pressure//2018_12_05__17_19_59.txt")
s9_hip_60 = Run(path + "Sample 9 higher pressure//2018_12_05__17_15_21.txt") 
s9_hip_67 = Run(path + "Sample 9 higher pressure//2018_12_05__17_10_17.txt")
s9_hip_75 = Run(path + "Sample 9 higher pressure//2018_12_05__17_03_53.txt") 
s9_hip2_30 = Run(path + "Sample 9 higher pressure 2//2018_12_07__16_59_52.txt")
s9_hip2_45 = Run(path + "Sample 9 higher pressure 2//2018_12_07__16_55_11.txt")
s9_hip2_52 = Run(path + "Sample 9 higher pressure 2//2018_12_07__16_50_31.txt")
s9_hip2_60 = Run(path + "Sample 9 higher pressure 2//2018_12_07__16_45_32.txt") 
s9_hip2_67 = Run(path + "Sample 9 higher pressure 2//2018_12_07__16_40_40.txt")
s9_hip2_75 = Run(path + "Sample 9 higher pressure 2//2018_12_07__16_35_58.txt") 
s9_lowp_above_30 = Run(path + "Sample 9 1-8_ above center//2018_12_07__11_37_47.txt")
s9_lowp_above_45 = Run(path + "Sample 9 1-8_ above center//2018_12_07__11_32_48.txt")
s9_lowp_above_52 = Run(path + "Sample 9 1-8_ above center//2018_12_07__11_27_11.txt")
s9_lowp_above_60 = Run(path + "Sample 9 1-8_ above center//2018_12_07__11_22_25.txt") 
s9_lowp_above_67 = Run(path + "Sample 9 1-8_ above center//2018_12_07__11_17_42.txt")
s9_lowp_above_75 = Run(path + "Sample 9 1-8_ above center//2018_12_07__11_12_10.txt")


runs = [s9_lowp_30,s9_lowp2_30,s9_medp_30,s9_medp2_30,s9_hip_30,s9_hip2_30]#[s9_medp2_30,s9_medp2_45,s9_medp2_52,s9_medp2_60,s9_medp2_67,s9_medp2_75,s9_hip2_30,s9_hip2_45,s9_hip2_52,s9_hip2_60,s9_hip2_67,s9_hip2_75]#[s9_lowp2_30,s9_lowp2_45,s9_lowp2_52,s9_lowp2_60,s9_lowp2_67,s9_lowp2_75,s9_medp2_30,s9_medp2_45,s9_medp2_52,s9_medp2_60,s9_medp2_67,s9_medp2_75,s9_hip2_30,s9_hip2_45,s9_hip2_52,s9_hip2_60,s9_hip2_67,s9_hip2_75]#[s9_lowp_30,s9_lowp_45,s9_lowp_52,s9_lowp_60,s9_lowp_67,s9_lowp_75,s9_lowp_above_30,s9_lowp_above_45,s9_lowp_above_52,s9_lowp_above_60,s9_lowp_above_67,s9_lowp_above_75]#[s9_lowp_75,s9_medp_75,s9_hip_75]#[s9_medp_30,s9_medp_45,s9_medp_52,s9_medp_60,s9_medp_67,s9_medp_75,s9_hip_30,s9_hip_45,s9_hip_52,s9_hip_60,s9_hip_67,s9_hip_75]#[s9_lowp_30,s9_lowp_45,s9_lowp_52,s9_lowp_60,s9_lowp_67,s9_lowp_75,s9_medp_30,s9_medp_45,s9_medp_52,s9_medp_60,s9_medp_67,s9_medp_75,s9_hip_30,s9_hip_45,s9_hip_52,s9_hip_60,s9_hip_67,s9_hip_75]#[s5_30,s5_45,s5_52,s5_60,s5_67,s5_75]#[s9_nobubbles_30,s9_nobubbles_45,s9_nobubbles_52,s9_nobubbles_60,s9_nobubbles_67,s9_nobubbles_75,s9_getter_30,s9_getter_45,s9_getter_52,s9_getter_60,s9_getter_67,s9_getter_75]#[s9_first_30,s9_first_45,s9_first_52,s9_first_60,s9_first_67,s9_first_75,s9_nobubbles_30,s9_nobubbles_45,s9_nobubbles_52,s9_nobubbles_60,s9_nobubbles_67,s9_nobubbles_75]#
labels=["0.2 barg","0.2 barg fixed","1.0 barg","1.0 barg fixed","1.5 barg","1.5 barg fixed"]#["1.0 barg, 30 deg","1.0 barg, 45 deg","1.0 barg, 52 deg","1.0 barg, 60 deg","1.0 barg, 67 deg","1.0 barg, 75 deg","1.5 barg, 30 deg","1.5 barg, 45 deg","1.5 barg, 52 deg","1.5 barg, 60 deg","1.5 barg, 67 deg","1.5 barg, 75 deg"]#["0.2 barg, 30 deg","0.2 barg, 45 deg","0.2 barg, 52 deg","0.2 barg, 60 deg","0.2 barg, 67 deg","0.2 barg, 75 deg","1.0 barg, 30 deg","1.0 barg, 45 deg","1.0 barg, 52 deg","1.0 barg, 60 deg","1.0 barg, 67 deg","1.0 barg, 75 deg","1.5 barg, 30 deg","1.5 barg, 45 deg","1.5 barg, 52 deg","1.5 barg, 60 deg","1.5 barg, 67 deg","1.5 barg, 75 deg"]#["center, 30 deg","center, 45 deg","center, 52 deg","center, 60 deg","center, 67 deg","center, 75 deg","1/8/" above, 30 deg","1/8/" above, 45 deg","1/8/" above, 52 deg","1/8/" above, 60 deg","1/8/" above, 67 deg","1/8/" above, 75 deg"]#["0.2 barg","1.0 barg","1.5 barg"]#["1.0 barg, 30 deg","1.0 barg, 45 deg","1.0 barg, 52 deg","1.0 barg, 60 deg","1.0 barg, 67 deg","1.0 barg, 75 deg","1.5 barg, 30 deg","1.5 barg, 45 deg","1.5 barg, 52 deg","1.5 barg, 60 deg","1.5 barg, 67 deg","1.5 barg, 75 deg"]#["0.2 barg, 30 deg","0.2 barg, 45 deg","0.2 barg, 52 deg","0.2 barg, 60 deg","0.2 barg, 67 deg","0.2 barg, 75 deg","1.0 barg, 30 deg","1.0 barg, 45 deg","1.0 barg, 52 deg","1.0 barg, 60 deg","1.0 barg, 67 deg","1.0 barg, 75 deg","1.5 barg, 30 deg","1.5 barg, 45 deg","1.5 barg, 52 deg","1.5 barg, 60 deg","1.5 barg, 67 deg","1.5 barg, 75 deg"]#["30 deg","45 deg","52 deg","60 deg","67 deg","75 deg"]#["before getter, 30 deg", "before getter, 45 deg", "before getter, 52 deg", "before getter, 60 deg", "before getter, 67 deg", "before getter, 75 deg", "after getter, 30 deg", "after getter, 45 deg", "after getter, 52 deg", "after getter, 60 deg", "after getter, 67 deg", "after getter, 75 deg"]#["75 degrees"]#["run 1, 30 deg", "run 1, 45 deg", "run 1, 52 deg", "run 1, 60 deg", "run 1, 67 deg", "run 1, 75 deg", "run 2, 30 deg", "run 2, 45 deg", "run 2, 52 deg", "run 2, 60 deg", "run 2, 67 deg", "run 2, 75 deg"]#

# Plot BRIDF data
sample_name="LZ Skived, 30 deg"
plot_runs(runs, title=sample_name+", 178nm", log=True, labels=labels, errorbars=True, legend_loc=0)
t0=time.time()

# Fit data
#fit_params = fit_parameters(get_independent_variables_and_relative_intensities(runs),p0=[0.24,1.2,.1,5.0],average_angle=4., precision=0.25,use_errs=True,use_spike=True)
#fit_params_ang = fit_parameters_and_angle(get_independent_variables_and_relative_intensities(runs),average_angle=4.)
#fit_ang = fit_params_ang[0]
#fit_params = fit_params_ang[1:]
#print("Fit parameters (rho_L, n, gamma): "+str(fit_params))
#print("Fit angle: "+str(fit_ang))
t1=time.time()
print("Fitting time: {0}".format(t1-t0))

# Plot BRIDF model from fits
#plot_TSTR_fit(30., 1., fit_params, label="fitted", color="k", average_angle=4., precision=0.25)
#plot_TSTR_fit(45., 1., fit_params, color="k", average_angle=4., precision=0.25)
#plot_TSTR_fit(60., 1., fit_params, color="k", average_angle=4., precision=0.25)
#plot_TSTR_fit(75., 1., fit_params, color="k", average_angle=4., precision=0.25)
#plt.text(0.1,0.1,r"Fit: $/rho_L$={0:.3f}, n={1:.2f}, $/gamma$={2:.3f}".format(*fit_params),transform=plt.gca().transAxes,fontsize=13)
#plt.text(0.1,0.2,r"Fit: $/rho_L$={0:.3f}, n={1:.2f}, $/gamma$={2:.3f}, K={3:.3f}".format(*fit_params),transform=plt.gca().transAxes,fontsize=13)
t2=time.time()
print("Plotting time: {0}".format(t2-t1))

# Now calculate hemispherical reflectance
# plt.figure()

# # Different incident angles to calculate for
# x = [30, 45, 60, 75]#[0,10,20,30, 45, 55, 60, 65, 70, 75, 80, 85]#

# y_diffuse = [reflectance_diffuse(theta, 1., 0.5, fit_params) for theta in x]
# y_specular = [reflectance_specular(theta, 1., 0.5, fit_params) for theta in x]
# y_total = [y_diffuse[i] + y_specular[i] for i in range(len(y_specular))]

# plt.plot(x, y_diffuse, label="diffuse")
# plt.plot(x, y_specular, label="specular")
# plt.plot(x, y_total, label="total")

# plt.xlabel("incident angle (degrees)")
# plt.ylabel("reflectance (fraction)")
# plt.legend()

# plt.title("Fitted "+sample_name+" Reflectance, 175 nm")
# t3=time.time()
# print("Reflectance calc time: {0}".format(t3-t2))
plt.show()
