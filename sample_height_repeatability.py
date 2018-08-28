import matplotlib.pyplot as plt
from file_reader import plot_runs, Run

path = "sample_height_repeatability/"

runs_no_pole = [
	Run(path + "2018_06_01__15_24_11.txt"),
	Run(path + "2018_06_01__15_28_01.txt"),
	Run(path + "2018_06_01__15_19_56.txt"),
	Run(path + "2018_06_01__15_31_34.txt")
]

runs_pole = [
	Run(path + "2018_06_04__14_48_21.txt"),
	Run(path + "2018_06_04__14_53_05.txt"),
	Run(path + "2018_06_04__14_57_33.txt"),
	Run(path + "2018_06_04__15_03_01.txt"),
	Run(path + "2018_06_04__15_07_27.txt")
]

plot_runs(runs_no_pole, title="NXT85 in air, 45 degrees, trying to repeat same sample height with no pole", 
	labels=["test 1", "test 2", "test 3", "test 4"], smooth=True, show=False, ylabel="relative intensity")
plot_runs(runs_pole, title="NXT85 in air, 45 degrees, trying to repeat same sample height with pole", 
	labels=["test 1", "test 2", "test 3", "test 4", "test 5"], smooth=True, show=False, ylabel="relative intensity")

plt.show()
