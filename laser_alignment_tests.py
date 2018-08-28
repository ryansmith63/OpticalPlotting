import matplotlib.pyplot as plt
from file_reader import plot_runs, Run

path = "laser_alignment_tests/"

runs = [
	Run(path + "2018_06_07__15_02_01.txt"),
	Run(path + "2018_06_07__15_15_04.txt"),
	Run(path + "2018_06_07__15_20_18.txt"),
]


plot_runs(runs, title="Laser alignment tests", 
	labels=["alignment 1", "alignment 2", "alignment 3"], smooth=True, show=False, ylabel="relative intensity")

plt.show()
