import matplotlib.pyplot as plt
import numpy as np
from file_reader import Run, get_independent_variables_and_relative_intensities
from plotting import plot_runs, plot_TSTR_fit
from TSTR_fit_new import fit_parameters, fit_parameters_and_angle, BRIDF_plotter, reflectance_diffuse, reflectance_specular, BRIDF
import time

path = "vuv_height_comparison_and_first_data/M18 turn polished/center of sample/"#\\M18 turn\\center of sample\\"
#path = "vuv_height_comparison_and_first_data\\M18 turn polished\\center of sample\\"#\\M18 turn\\center of sample\\"
m18pol30 = Run(path + "2018_08_30__11_23_42.txt") 
m18pol45 = Run(path + "2018_08_30__11_29_08.txt") 
m18pol60 = Run(path + "2018_08_30__11_34_35.txt") 
m18pol75 = Run(path + "2018_08_30__11_39_45.txt") 
# nxt8530 = Run(path + "2018_08_30__14_13_18.txt") 
# nxt8545 = Run(path + "2018_08_30__14_19_14.txt") 
# nxt8560 = Run(path + "2018_08_30__14_24_51.txt") 
# nxt8575 = Run(path + "2018_08_30__14_30_46.txt")
# skived30 = Run(path + "2018_08_29__17_29_38.txt") 
# skived45 = Run(path + "2018_08_29__17_35_12.txt") 
# skived60 = Run(path + "2018_08_29__17_40_57.txt") 
# skived75 = Run(path + "2018_08_29__17_46_39.txt") 
# m1830 = Run(path + "2018_08_30__15_16_18.txt") # Fit to 36 deg
# m1845 = Run(path + "2018_08_30__15_22_02.txt") # Fit to 46 deg
# m1860 = Run(path + "2018_08_30__15_29_09.txt") # Fit to 60 deg
# m1875 = Run(path + "2018_08_30__15_35_27.txt") # Fit to 85 deg
# m1730 = Run(path + "2018_08_30__16_18_57.txt")
# m1745 = Run(path + "2018_08_30__16_11_27.txt")
# m1760 = Run(path + "2018_08_30__16_06_11.txt")
# m1775 = Run(path + "2018_08_30__15_59_33.txt")
# nx30 = Run(path + "2018_08_29__15_57_32.txt")
# nx45 = Run(path + "2018_08_29__16_02_55.txt")
# nx60 = Run(path + "2018_08_29__16_08_16.txt")
# nx75 = Run(path + "2018_08_29__16_13_39.txt")

runs = [m18pol30, m18pol45, m18pol60, m18pol75]#[m1875]#
labels=["30 degrees", "45 degrees", "60 degrees", "75 degrees"]#["75 degrees"]#

# Plot BRIDF data
sample_name="M18 turn polished"
plot_runs(runs, title=sample_name+", 175nm", log=True, labels=labels, errorbars=True)
#plt.show()

# In a loop: ask for input params, plot model overlaid on data
get_input=True
params=[.24,1.27,0.1,5.0,2.0]
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
    plot_runs(runs, title=sample_name+", 175nm", log=True, labels=labels, errorbars=True)
	# Trying setting incident angles by hand, for better fits to specular spike
    plot_TSTR_fit(29., 1., fit_params, label="fitted", color="k", average_angle=4., precision=0.25)
    plot_TSTR_fit(42., 1., fit_params, color="k", average_angle=4., precision=0.25)
    plot_TSTR_fit(60., 1., fit_params, color="k", average_angle=4., precision=0.25)
    plot_TSTR_fit(77., 1., fit_params, color="k", average_angle=4., precision=0.25)
    plt.text(0.05,0.1,r"Fit: $\rho_L$={0:.3f}, n={1:.2f}, $\gamma$={2:.3f}, K={3:.3f}, $\sigma (\theta_i)$={4:.3f}".format(*params),transform=plt.gca().transAxes,fontsize=11)
    plt.show()
    get_input=input("New parameters? (y/n) ").lower()=='y'
