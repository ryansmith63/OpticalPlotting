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
"""
fit_params_1 = fit_parameters(get_independent_variables_and_relative_intensities([s1_30]),
	 p0=[0.24,1.2,.1,5.0],average_angle=4., precision=0.25,use_errs=True,use_spike=True)
plot_runs(sample_1_runs, title="Sample 1 (LUX), 175nm, LXe", log=True, 
	labels=["30 degrees", "45 degrees", "52 degrees", "60 degrees", "67 degrees", "75 degrees"])
plot_TSTR_fit(30., 1., fit_params, label="fitted", color="k", average_angle=4., precision=0.25)
plot_TSTR_fit(45., 1., fit_params, color="k", average_angle=4., precision=0.25)
plot_TSTR_fit(60., 1., fit_params, color="k", average_angle=4., precision=0.25)
plot_TSTR_fit(75., 1., fit_params, color="k", average_angle=4., precision=0.25)
"""
plot_runs(sample_1_runs, title="Sample 1 (LUX), 175nm, LXe", log=True, 
	labels=["30 degrees", "45 degrees", "52 degrees", "60 degrees", "67 degrees", "75 degrees"])
plot_runs(sample_3_runs, title="Sample 3 (M18 Turn), 175nm, LXe", log=True, 
	labels=["30 degrees", "45 degrees", "52 degrees", "60 degrees", "67 degrees", "75 degrees"])
plot_runs(sample_9_runs, title="Sample 9 (Skived LZ), 175nm, LXe", log=True, 
	labels=["30 degrees", "45 degrees", "52 degrees", "60 degrees", "67 degrees", "75 degrees"])


plt.show()
