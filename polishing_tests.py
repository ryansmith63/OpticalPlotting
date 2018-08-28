import matplotlib.pyplot as plt
from file_reader import plot_runs, Run

path = "polishing_tests/"

runs_unpolished = [
	Run(path + "2018_06_07__15_53_10.txt"),
	Run(path + "2018_06_07__15_57_42.txt"),
	Run(path + "2018_06_07__16_01_32.txt"),
	Run(path + "2018_06_07__16_05_22.txt"),
	Run(path + "2018_06_07__16_09_12.txt"),
	Run(path + "2018_06_07__16_13_01.txt"),
	Run(path + "2018_06_07__16_17_38.txt"),
	Run(path + "2018_06_07__16_17_38.txt")
]

runs_polished = [
	Run(path + "2018_06_08__14_30_48.txt"),
	Run(path + "2018_06_08__14_34_56.txt"),
	Run(path + "2018_06_08__14_38_38.txt"),
	Run(path + "2018_06_08__14_42_37.txt"),
	Run(path + "2018_06_08__14_47_53.txt"),
	Run(path + "2018_06_08__14_52_12.txt")
]


plot_runs(runs_unpolished, title="M18 turned, adjusting sample height", 
	labels=["height 1 \n(beam hitting center of sample)", "height 2", "height 3", "height 4", 
	"height 5", "height 6", "height 7", "height 8\n(beam hitting near top of sample)"], smooth=True, show=False)
plot_runs(runs_polished, title="M18 polished, adjusting sample height", 
	labels=["height 1 \n(beam hitting center of sample)", "height 2", "height 3", "height 4", 
	"height 5", "height 6\n(beam hitting near top of sample)"], smooth=True, show=False)


plt.show()
