import matplotlib.pyplot as plt
import numpy as np
from file_reader import Run, get_independent_variables_and_relative_intensities
from plotting import plot_runs, plot_TSTR_fit
from TSTR_fit_new import fit_parameters, fit_parameters_and_angle, BRIDF_plotter, reflectance_diffuse, reflectance_specular, BRIDF, chi_squared
import time

#path = "vuv_height_comparison_and_first_data/M18 turn polished/center of sample/"#\\M18 turn\\center of sample\\"
#path = "vuv_height_comparison_and_first_data\\Skived LZ\\center of sample\\"#\\M18 turn\\center of sample\\"
path = "2nd Xenon Run Measurements\\"
# m18pol30 = Run(path + "2018_08_30__11_23_42.txt") 
# m18pol45 = Run(path + "2018_08_30__11_29_08.txt") 
# m18pol60 = Run(path + "2018_08_30__11_34_35.txt") 
# m18pol75 = Run(path + "2018_08_30__11_39_45.txt") 
# nxt8530 = Run(path + "2018_08_30__14_13_18.txt") 
# nxt8545 = Run(path + "2018_08_30__14_19_14.txt") 
# nxt8560 = Run(path + "2018_08_30__14_24_51.txt") 
# nxt8575 = Run(path + "2018_08_30__14_30_46.txt")
# skived30 = Run(path + "2018_08_29__17_29_38.txt") 
# skived30.change_theta_i(29)
# skived45 = Run(path + "2018_08_29__17_35_12.txt") 
# skived45.change_theta_i(43)
# skived60 = Run(path + "2018_08_29__17_40_57.txt") 
# skived60.change_theta_i(59)
# skived75 = Run(path + "2018_08_29__17_46_39.txt") 
# skived75.change_theta_i(74)
#m1830 = Run(path + "2018_08_30__15_16_18.txt") # Fit to 36 deg
#m1830.change_theta_i(33)
#m1845 = Run(path + "2018_08_30__15_22_02.txt") # Fit to 46 deg
#m1860 = Run(path + "2018_08_30__15_29_09.txt") # Fit to 60 deg
#m1860.change_theta_i(61)
#m1875 = Run(path + "2018_08_30__15_35_27.txt") # Fit to 85 deg
#m1875.change_theta_i(76)
# m1730 = Run(path + "2018_08_30__16_18_57.txt")
# m1745 = Run(path + "2018_08_30__16_11_27.txt")
# m1760 = Run(path + "2018_08_30__16_06_11.txt")
# m1775 = Run(path + "2018_08_30__15_59_33.txt")
# nx30 = Run(path + "2018_08_29__15_57_32.txt")
# nx45 = Run(path + "2018_08_29__16_02_55.txt")
# nx60 = Run(path + "2018_08_29__16_08_16.txt")
# nx75 = Run(path + "2018_08_29__16_13_39.txt")
s9_lowp_30 = Run(path + "Sample 9 lower pressure\\2018_12_05__13_52_04.txt")
#s9_lowp_30.change_theta_i(30)
s9_lowp_45 = Run(path + "Sample 9 lower pressure\\2018_12_05__13_47_09.txt")
#s9_lowp_45.change_theta_i(45)
s9_lowp_52 = Run(path + "Sample 9 lower pressure\\2018_12_05__13_42_21.txt")
#s9_lowp_52.change_theta_i(52)
s9_lowp_60 = Run(path + "Sample 9 lower pressure\\2018_12_05__13_37_31.txt") 
#s9_lowp_60.change_theta_i(60)
s9_lowp_67 = Run(path + "Sample 9 lower pressure\\2018_12_05__13_32_47.txt")
#s9_lowp_67.change_theta_i(67)
s9_lowp_75 = Run(path + "Sample 9 lower pressure\\2018_12_05__13_27_55.txt")
#s9_lowp_75.change_theta_i(75)

vacuum_path = "First Xe Run Measurements\\follow-up vacuum measurements\\"

vacuum_9_path = vacuum_path + "Sample 9\\"
v9_75 = Run(vacuum_9_path + "2018_11_09__14_22_32.txt")
v9_75.change_theta_i(74)
v9_67 = Run(vacuum_9_path + "2018_11_09__14_33_41.txt")
v9_67.change_theta_i(65)
v9_60 = Run(vacuum_9_path + "2018_11_09__14_40_54.txt")
v9_60.change_theta_i(56)
v9_52 = Run(vacuum_9_path + "2018_11_09__14_47_23.txt")
v9_52.change_theta_i(48)
v9_45 = Run(vacuum_9_path + "2018_11_09__14_53_56.txt")
v9_45.change_theta_i(40)
v9_30 = Run(vacuum_9_path + "2018_11_09__15_00_15.txt")
v9_30.change_theta_i(25)
vacuum_9_runs = [v9_30, v9_45, v9_52, v9_60, v9_67, v9_75]

runs = vacuum_9_runs#[s9_lowp_30,s9_lowp_45,s9_lowp_52,s9_lowp_60,s9_lowp_67,s9_lowp_75]#[m1875]#
labels=["30 degrees", "45 degrees", "52 degrees", "60 degrees", "67 degrees", "75 degrees"]#["75 degrees"]#

run_data = get_independent_variables_and_relative_intensities(runs)
# Plot BRIDF data
sample_name="LZ skived"
plot_runs(runs, title=sample_name+", 178 nm", log=True, labels=labels, errorbars=True)
#plt.show()

# In a loop: ask for input params, plot model overlaid on data
get_input=True
params=[.8,1.55,0.1,-1,2.0]
param_names=["rho_L","n","gamma","K","sigma_theta_i"]
fig_num=0
while get_input:
    fig_num += 1
    for ii,name in enumerate(param_names):
        while True:
            param = input('Input '+name+' (or enter to use '+str(params[ii])+'): ')
            try:
                params[ii] = float(param)
                break
            except ValueError: 
                if param == '': break
                print("Error: input a number")
    sigma_theta_i = params[-1]
    fit_params = params[:-1]
    #plt.figure(fig_num)
    plot_runs(runs, title=sample_name+", 175 nm", log=True, labels=labels, errorbars=True)
	# Trying setting incident angles by hand, for better fits to specular spike
    n_LXe_178 = 1.69
    n=1.0#n_LXe_178
    #sigma_theta_i=2.0
    precision=0.25
    average_angle=4.
    plot_TSTR_fit(25., n, fit_params, label="fitted", color="r", average_angle=average_angle, precision=precision, sigma_theta_i=sigma_theta_i)
    plot_TSTR_fit(40., n, fit_params, color="g", average_angle=average_angle, precision=precision, sigma_theta_i=sigma_theta_i)
    plot_TSTR_fit(48., n, fit_params, color="b", average_angle=average_angle, precision=precision, sigma_theta_i=sigma_theta_i)
    plot_TSTR_fit(56., n, fit_params, color="m", average_angle=average_angle, precision=precision, sigma_theta_i=sigma_theta_i)
    plot_TSTR_fit(65., n, fit_params, color="c", average_angle=average_angle, precision=precision, sigma_theta_i=sigma_theta_i)
    plot_TSTR_fit(74., n, fit_params, color="y", average_angle=average_angle, precision=precision, sigma_theta_i=sigma_theta_i)

    #plt.text(0.05,0.1,r"Fit: $\rho_L$={0:.3f}, n={1:.2f}, $\gamma$={2:.3f}, K={3:.3f}, $\sigma (\theta_i)$={4:.3f}".format(*params),transform=plt.gca().transAxes,fontsize=11)
    chi_sq = chi_squared(run_data[0], run_data[1], run_data[2], fit_params, average_angle=average_angle, precision=precision)
    print("Reduced chi^2: ",chi_sq)
    plt.show()
    get_input=input("New parameters? (y/n) ").lower()=='y'
