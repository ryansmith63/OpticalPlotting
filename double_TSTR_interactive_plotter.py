import matplotlib.pyplot as plt
import numpy as np
from file_reader import Run, get_independent_variables_and_relative_intensities
from plotting import plot_runs, plot_TSTR_fit, plot_double_TSTR_fit
from TSTR_fit_new import fit_parameters, fit_parameters_and_angle, BRIDF_plotter, reflectance_diffuse, reflectance_specular, BRIDF
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
param_names=["rho_L_1","n_1","gamma_1","K_1",
             "gamma_2","frac_1"]
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
    
    #plt.figure(fig_num)
    plot_runs(runs, title=sample_name+", 175nm", log=True, labels=labels, errorbars=True)
	# Trying setting incident angles by hand, for better fits to specular spike
    plot_double_TSTR_fit(30., 1.69, params, labels=["fitted total", "fitted 1 (LXe-PTFE)", "fitted 2 (LXe-gas)"], colors=["k", "r", "b"], average_angle=4., precision=0.25)
    print("done with 30...")
    plot_double_TSTR_fit(45., 1.69, params, colors=["k", "r", "b"], average_angle=4., precision=0.25)
    print("done with 45...")
    #plot_double_TSTR_fit(52., 1.69, params, colors=["k", "r", "b"], average_angle=4., precision=0.25)
    #print("done with 52...")
    plot_double_TSTR_fit(60., 1.69, params, colors=["k", "r", "b"], average_angle=4., precision=0.25)
    print("done with 60...")
    #plot_double_TSTR_fit(67., 1.69, params, colors=["k", "r", "b"], average_angle=4., precision=0.25)
    #print("done with 67...")
    plot_double_TSTR_fit(75., 1.69, params, colors=["k", "r", "b"], average_angle=4., precision=0.25)
    print("done with 75...")

    #plt.text(0.05,0.1,r"Fit: $\rho_L$={0:.3f}, n={1:.2f}, $\gamma$={2:.3f}, K={3:.3f}, $\sigma (\theta_i)$={4:.3f}".format(*params),transform=plt.gca().transAxes,fontsize=11)
    plt.show()
    get_input=input("New parameters? (y/n) ").lower()=='y'
