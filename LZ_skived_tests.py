import matplotlib.pyplot as plt
from file_reader import plot_runs, Run

path = "LZ_skived_tests/"

runs_unmarked = [
	Run(path + "2018_06_13__16_11_58.txt"),
	Run(path + "2018_06_13__16_19_19.txt"),
	Run(path + "2018_06_13__16_25_51.txt"),
	Run(path + "2018_06_13__16_33_40.txt")
]

# reversed in time because I went outer to inner
runs_unmarked_adjusted = [
	Run(path + "2018_06_13__17_02_41.txt"),
	Run(path + "2018_06_13__16_59_03.txt"),
	Run(path + "2018_06_13__16_55_14.txt"),
	Run(path + "2018_06_13__16_50_48.txt"),
	Run(path + "2018_06_13__16_47_03.txt")
]

runs_marked = [
	Run(path + "2018_06_14__12_32_29.txt"),
	Run(path + "2018_06_14__12_36_26.txt"),
	Run(path + "2018_06_14__12_40_10.txt"),
	Run(path + "2018_06_14__12_44_20.txt"),
	Run(path + "2018_06_14__12_48_09.txt")
]

plot_runs(runs_unmarked, title="LZ skived, unmarked side, not pushed all the way in, adjusting sample height", 
	labels=["height 1 \n(beam hitting near center of sample)", "height 2", "height 3", "height 4\n(beam hitting near top of sample)"], smooth=True, show=False)
plot_runs(runs_unmarked_adjusted, title="LZ skived, unmarked side, adjusting sample height", 
	labels=["height 1 \n(beam hitting near center of sample)", "height 2", "height 3", "height 4", 
	"height 5\n(beam hitting near top of sample)"], smooth=True, show=False)
plot_runs(runs_marked, title="LZ skived, marked side, adjusting sample height", 
	labels=["height 1 \n(beam hitting near center of sample)", "height 2", "height 3", "height 4", 
	"height 5\n(beam hitting near top of sample)"], smooth=True, show=False)

plt.show()
