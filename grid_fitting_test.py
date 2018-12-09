import matplotlib.pyplot as plt
import numpy as np
from file_reader import Run, get_independent_variables_and_relative_intensities
from plotting import plot_runs, plot_TSTR_fit
from TSTR_fit_new import fit_parameters, BRIDF_plotter, reflectance_diffuse, reflectance_specular, fit_parameters_new, chi_squared

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



fit_params_old = fit_parameters(get_independent_variables_and_relative_intensities(v1_45), p0=[0.24,1.2,.1,5.0], 
	average_angle=0, precision=-1, use_errs=False, use_spike=False)
chi_squared_old = chi_squared(v1_45.independent_variables_array, v1_45.relative_intensities, v1_45.relative_std, fit_params_old)

fit_params_new = fit_parameters_new(v1_45.independent_variables_array, v1_45.relative_intensities, v1_45.relative_std,
	#rho_start=0.15, rho_end=0.25, rho_num=10, n_start=1., n_end=1.2, n_num=10, gamma_start=0.05, gamma_end=0.2, gamma_num=10, plot=True, show=False) # fine
	rho_start=0.05, rho_end=0.5, rho_num=10, n_start=1., n_end=1.5, n_num=10, gamma_start=0.01, gamma_end=0.5, gamma_num=10, plot=True, show=False) # coarse

chi_squared_new = chi_squared(v1_45.independent_variables_array, v1_45.relative_intensities, v1_45.relative_std, fit_params_new)

plot_runs(v1_45, title="Sample 1 (LUX), 178nm, Vacuum, 45 degrees", label="data", log=True)
plot_TSTR_fit(45., 1., fit_params_old, label="fitted with old fitter", color="b", average_angle=0, precision=-1, 
	fit_text="Old Fit, chi^2=" + str(round(chi_squared_old,4)) + ": ", fit_text_offset=0.05)
plot_TSTR_fit(45., 1., fit_params_new, label="fitted with new fitter", color="g", average_angle=0, precision=-1, 
	fit_text="New Fit, chi^2=" + str(round(chi_squared_new,4)) + ": ")

plt.show()



