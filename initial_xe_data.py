import matplotlib.pyplot as plt
import numpy as np
from file_reader import Run, get_independent_variables_and_relative_intensities
from plotting import plot_runs, plot_TSTR_fit
from TSTR_fit_new import fit_parameters, BRIDF_plotter, reflectance_diffuse, reflectance_specular

initial_xe_path = "First Xe Run Measurements/first measurements with no bubbles in cell 11-01-2/"

sample_1_path = initial_xe_path + "Sample 1/"
s1_30 = Run(sample_1_path + "2018_11_02__18_28_43.txt")
s1_45 = Run(sample_1_path + "2018_11_02__18_35_18.txt")
s1_52 = Run(sample_1_path + "2018_11_02__18_41_20.txt")
s1_60 = Run(sample_1_path + "2018_11_02__18_47_03.txt")
s1_67 = Run(sample_1_path + "2018_11_02__18_51_55.txt")
s1_75 = Run(sample_1_path + "2018_11_02__18_57_35.txt")
sample_1_runs = [s1_30, s1_45, s1_52, s1_60, s1_67, s1_75]

sample_3_path = initial_xe_path + "Sample 3/"
s3_30 = Run(sample_3_path + "2018_11_02__19_37_08.txt")
s3_45 = Run(sample_3_path + "2018_11_02__19_42_42.txt")
s3_52 = Run(sample_3_path + "2018_11_02__19_47_37.txt")
s3_60 = Run(sample_3_path + "2018_11_02__19_52_30.txt")
s3_67 = Run(sample_3_path + "2018_11_02__19_59_24.txt")
s3_75 = Run(sample_3_path + "2018_11_02__20_04_42.txt")
sample_3_runs = [s3_30, s3_45, s3_52, s3_60, s3_67, s3_75]

sample_9_path = initial_xe_path + "Sample 9/"
s9_75 = Run(sample_9_path + "2018_11_02__20_11_44.txt")
s9_67 = Run(sample_9_path + "2018_11_02__20_17_11.txt")
s9_60 = Run(sample_9_path + "2018_11_02__20_22_13.txt")
s9_52 = Run(sample_9_path + "2018_11_02__20_27_27.txt")
s9_45 = Run(sample_9_path + "2018_11_02__20_32_34.txt")
s9_30 = Run(sample_9_path + "2018_11_02__20_38_01.txt")
sample_9_runs = [s9_30, s9_45, s9_52, s9_60, s9_67, s9_75]

vacuum_path = "First Xe Run Measurements/follow-up vacuum measurements/"

vacuum_1_path = vacuum_path + "Sample 1/"
v1_75 = Run(vacuum_1_path + "2018_11_09__13_41_18.txt")
v1_67 = Run(vacuum_1_path + "2018_11_09__13_47_53.txt")
v1_60 = Run(vacuum_1_path + "2018_11_09__13_55_54.txt")
v1_52 = Run(vacuum_1_path + "2018_11_09__14_02_45.txt")
v1_45 = Run(vacuum_1_path + "2018_11_09__14_08_54.txt")
v1_30 = Run(vacuum_1_path + "2018_11_09__14_14_55.txt")
vacuum_1_runs = [v1_30, v1_45, v1_52, v1_60, v1_67, v1_75]

vacuum_3_path = vacuum_path + "Sample 3/"
v3_75 = Run(vacuum_3_path + "2018_11_09__11_44_37.txt")
v3_67 = Run(vacuum_3_path + "2018_11_09__11_49_54.txt")
v3_60 = Run(vacuum_3_path + "2018_11_09__11_54_58.txt")
v3_52 = Run(vacuum_3_path + "2018_11_09__12_00_07.txt")
v3_45 = Run(vacuum_3_path + "2018_11_09__12_04_53.txt")
v3_30 = Run(vacuum_3_path + "2018_11_09__12_13_51.txt")
vacuum_3_runs = [v3_30, v3_45, v3_52, v3_60, v3_67, v3_75]

vacuum_9_path = vacuum_path + "Sample 9/"
v9_75 = Run(vacuum_9_path + "2018_11_09__14_22_32.txt")
v9_67 = Run(vacuum_9_path + "2018_11_09__14_33_41.txt")
v9_60 = Run(vacuum_9_path + "2018_11_09__14_40_54.txt")
v9_52 = Run(vacuum_9_path + "2018_11_09__14_47_23.txt")
v9_45 = Run(vacuum_9_path + "2018_11_09__14_53_56.txt")
v9_30 = Run(vacuum_9_path + "2018_11_09__15_00_15.txt")
vacuum_9_runs = [v9_30, v9_45, v9_52, v9_60, v9_67, v9_75]

all_1_runs = sample_1_runs + vacuum_1_runs
all_3_runs = sample_3_runs + vacuum_3_runs
all_9_runs = sample_9_runs + vacuum_9_runs

# vacuum individual fits

fit_params = fit_parameters(get_independent_variables_and_relative_intensities(v1_30), p0=[0.24,1.2,.1,5.0],average_angle=4., precision=0.25, use_errs=True, use_spike=True)
plot_runs(v1_30, title="Sample 1 (LUX), 178nm, Vacuum, 30 degrees", label="data", log=True)
plot_TSTR_fit(30., 1., fit_params, label="fitted", color="k", average_angle=4., precision=0.25)
print("done with 30")

fit_params = fit_parameters(get_independent_variables_and_relative_intensities(v1_45), p0=[0.24,1.2,.1,5.0],average_angle=4., precision=0.25, use_errs=True, use_spike=True)
plot_runs(v1_45, title="Sample 1 (LUX), 178nm, Vacuum, 45 degrees", label="data", log=True)
plot_TSTR_fit(45., 1., fit_params, label="fitted", color="k", average_angle=4., precision=0.25)
print("done with 45")

fit_params = fit_parameters(get_independent_variables_and_relative_intensities(v1_52), p0=[0.24,1.2,.1,5.0],average_angle=4., precision=0.25, use_errs=True, use_spike=True)
plot_runs(v1_52, title="Sample 1 (LUX), 178nm, Vacuum, 52 degrees", label="data", log=True)
plot_TSTR_fit(52., 1., fit_params, label="fitted", color="k", average_angle=4., precision=0.25)
print("done with 52")

fit_params = fit_parameters(get_independent_variables_and_relative_intensities(v1_60), p0=[0.24,1.2,.1,5.0],average_angle=4., precision=0.25, use_errs=True, use_spike=True)
plot_runs(v1_60, title="Sample 1 (LUX), 178nm, Vacuum, 60 degrees", label="data", log=True)
plot_TSTR_fit(60., 1., fit_params, label="fitted", color="k", average_angle=4., precision=0.25)
print("done with 60")

fit_params = fit_parameters(get_independent_variables_and_relative_intensities(v1_67), p0=[0.24,1.2,.1,5.0],average_angle=4., precision=0.25, use_errs=True, use_spike=True)
plot_runs(v1_67, title="Sample 1 (LUX), 178nm, Vacuum, 67 degrees", label="data", log=True)
plot_TSTR_fit(67., 1., fit_params, label="fitted", color="k", average_angle=4., precision=0.25)
print("done with 67")

fit_params = fit_parameters(get_independent_variables_and_relative_intensities(v1_75), p0=[0.24,1.2,.1,5.0],average_angle=4., precision=0.25, use_errs=True, use_spike=True)
plot_runs(v1_75, title="Sample 1 (LUX), 178nm, Vacuum, 75 degrees", label="data", log=True)
plot_TSTR_fit(75., 1., fit_params, label="fitted", color="k", average_angle=4., precision=0.25)
print("done with 75")


"""
# LXe individual fits
fit_params = fit_parameters(get_independent_variables_and_relative_intensities(s1_30), p0=[0.24,1.2,.1,5.0],average_angle=4., precision=0.25, use_errs=True, use_spike=True)
plot_runs(s1_30, title="Sample 1 (LUX), 178nm, LXe, 30 degrees", label="data", log=True)
plot_TSTR_fit(30., 1.69, fit_params, label="fitted", color="k", average_angle=4., precision=0.25)
print("done with 30")

fit_params = fit_parameters(get_independent_variables_and_relative_intensities(s1_45), p0=[0.24,1.2,.1,5.0],average_angle=4., precision=0.25, use_errs=True, use_spike=True)
plot_runs(s1_45, title="Sample 1 (LUX), 178nm, LXe, 45 degrees", label="data", log=True)
plot_TSTR_fit(45., 1.69, fit_params, label="fitted", color="k", average_angle=4., precision=0.25)
print("done with 45")

fit_params = fit_parameters(get_independent_variables_and_relative_intensities(s1_52), p0=[0.24,1.2,.1,5.0],average_angle=4., precision=0.25, use_errs=True, use_spike=True)
plot_runs(s1_52, title="Sample 1 (LUX), 178nm, LXe, 52 degrees", label="data", log=True)
plot_TSTR_fit(52., 1.69, fit_params, label="fitted", color="k", average_angle=4., precision=0.25)
print("done with 52")

fit_params = fit_parameters(get_independent_variables_and_relative_intensities(s1_60), p0=[0.24,1.2,.1,5.0],average_angle=4., precision=0.25, use_errs=True, use_spike=True)
plot_runs(s1_60, title="Sample 1 (LUX), 178nm, LXe, 60 degrees", label="data", log=True)
plot_TSTR_fit(60., 1.69, fit_params, label="fitted", color="k", average_angle=4., precision=0.25)
print("done with 60")

fit_params = fit_parameters(get_independent_variables_and_relative_intensities(s1_67), p0=[0.24,1.2,.1,5.0],average_angle=4., precision=0.25, use_errs=True, use_spike=True)
plot_runs(s1_67, title="Sample 1 (LUX), 178nm, LXe, 67 degrees", label="data", log=True)
plot_TSTR_fit(67., 1.69, fit_params, label="fitted", color="k", average_angle=4., precision=0.25)
print("done with 67")

fit_params = fit_parameters(get_independent_variables_and_relative_intensities(s1_75), p0=[0.24,1.2,.1,5.0],average_angle=4., precision=0.25, use_errs=True, use_spike=True)
plot_runs(s1_75, title="Sample 1 (LUX), 178nm, LXe, 75 degrees", label="data", log=True)
plot_TSTR_fit(75., 1.69, fit_params, label="fitted", color="k", average_angle=4., precision=0.25)
print("done with 75")
"""


#fit_params_1 = [0.172, 1.33, 0.355, 14.938] # from line below
#fit_params_1 = fit_parameters(get_independent_variables_and_relative_intensities(vacuum_1_runs), p0=[0.24,1.2,.1,5.0],average_angle=4., precision=0.25, use_errs=True, use_spike=True)

"""
fit_params_1 = fit_parameters(get_independent_variables_and_relative_intensities(sample_1_runs), p0=[0.24,1.2,.1,5.0],average_angle=4., precision=0.25, use_errs=True, use_spike=True)


plot_runs(vacuum_1_runs, title="Sample 1 (LUX), 175nm, Vacuum\n Fit from sample 1 LXe data", log=True, 
	labels=["30 degrees", "45 degrees", "52 degrees", "60 degrees", "67 degrees", "75 degrees"])
plot_TSTR_fit(30., 1., fit_params_1, label="fitted", color="k", average_angle=4., precision=0.25)
plot_TSTR_fit(45., 1., fit_params_1, color="k", average_angle=4., precision=0.25)
plot_TSTR_fit(52., 1., fit_params_1, color="k", average_angle=4., precision=0.25)
plot_TSTR_fit(60., 1., fit_params_1, color="k", average_angle=4., precision=0.25)
plot_TSTR_fit(67., 1., fit_params_1, color="k", average_angle=4., precision=0.25)
plot_TSTR_fit(75., 1., fit_params_1, color="k", average_angle=4., precision=0.25)

plot_runs(sample_1_runs, title="Sample 1 (LUX), 175nm, LXe\n Fit from this data", log=True, 
	labels=["30 degrees", "45 degrees", "52 degrees", "60 degrees", "67 degrees", "75 degrees"])
plot_TSTR_fit(30., 1.69, fit_params_1, label="fitted", color="k", average_angle=4., precision=0.25)
plot_TSTR_fit(45., 1.69, fit_params_1, color="k", average_angle=4., precision=0.25)
plot_TSTR_fit(52., 1.69, fit_params_1, color="k", average_angle=4., precision=0.25)
plot_TSTR_fit(60., 1.69, fit_params_1, color="k", average_angle=4., precision=0.25)
plot_TSTR_fit(67., 1.69, fit_params_1, color="k", average_angle=4., precision=0.25)
plot_TSTR_fit(75., 1.69, fit_params_1, color="k", average_angle=4., precision=0.25)
"""

# this plots the individual plots
"""
plot_runs(sample_1_runs, title="Sample 1 (LUX), 178nm, LXe", log=True, 
	labels=["30 degrees", "45 degrees", "52 degrees", "60 degrees", "67 degrees", "75 degrees"])
plot_runs(sample_3_runs, title="Sample 3 (M18 Turn), 178nm, LXe", log=True, 
	labels=["30 degrees", "45 degrees", "52 degrees", "60 degrees", "67 degrees", "75 degrees"])
plot_runs(sample_9_runs, title="Sample 9 (Skived LZ), 178nm, LXe", log=True, 
	labels=["30 degrees", "45 degrees", "52 degrees", "60 degrees", "67 degrees", "75 degrees"])


plot_runs(vacuum_1_runs, title="Sample 1 (LUX), 178nm, Vacuum", log=True, 
	labels=["75 degrees", "67 degrees", "60 degrees", "52 degrees", "45 degrees", "30 degrees"])
plot_runs(vacuum_3_runs, title="Sample 3 (M18 Turn), 178nm, Vacuum", log=True, 
	labels=["75 degrees", "67 degrees", "60 degrees", "52 degrees", "45 degrees", "30 degrees"])
plot_runs(vacuum_9_runs, title="Sample 9 (Skived LZ), 178nm, Vacuum", log=True, 
	labels=["75 degrees", "67 degrees", "60 degrees", "52 degrees", "45 degrees", "30 degrees"])
"""

# this plots the LXe vacuum comparisons
"""
plot_runs(vacuum_1_runs, title="Sample 1 (LUX), 178nm", log=True, label="Vacuum", color="b")
plot_runs(sample_1_runs, figure=False, label="LXe", color="r")

plot_runs(vacuum_3_runs, title="Sample 3 (M18 Turn), 178nm", log=True, label="Vacuum", color="b")
plot_runs(sample_3_runs, figure=False, label="LXe", color="r")

plot_runs(vacuum_9_runs, title="Sample 9 (Skived LZ), 178nm", log=True, label="Vacuum", color="b")
plot_runs(sample_9_runs, figure=False, label="LXe", color="r")
"""

plt.show()
