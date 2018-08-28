import matplotlib.pyplot as plt
from file_reader import plot_runs, Run

path = "sample_height_tests/"

runs_nxt_85 = [
	Run(path + "2018_05_30__13_11_32.txt"),
	Run(path + "2018_05_30__13_15_19.txt"),
	Run(path + "2018_05_30__13_18_58.txt"),
	Run(path + "2018_05_30__13_22_31.txt"),
	Run(path + "2018_05_30__13_25_58.txt"),
	Run(path + "2018_05_30__13_29_46.txt"),
	Run(path + "2018_05_30__13_33_24.txt"),
	Run(path + "2018_05_30__13_36_54.txt"),
	Run(path + "2018_05_30__13_40_30.txt")
]

runs_m18_turn = [
	Run(path + "2018_05_31__12_44_30.txt"),
	Run(path + "2018_05_31__12_49_09.txt"),
	Run(path + "2018_05_31__13_00_57.txt"),
	Run(path + "2018_05_31__13_04_27.txt"),
	Run(path + "2018_05_31__13_07_57.txt"),
	Run(path + "2018_05_31__13_11_28.txt"),
	Run(path + "2018_05_31__13_14_54.txt")
]

runs_old_sample = [
	Run(path + "2018_05_31__13_25_38.txt"),
	Run(path + "2018_05_31__14_14_42.txt"),
	Run(path + "2018_05_31__14_18_09.txt"),
	Run(path + "2018_05_31__14_21_59.txt"),
	Run(path + "2018_05_31__14_25_21.txt"),
	Run(path + "2018_05_31__14_28_42.txt"),
	Run(path + "2018_05_31__14_32_14.txt")
]

plot_runs(runs_nxt_85, title="NXT85 in air, adjusting sample height", 
	labels=["height 1 \n(beam hitting center of sample)", "height 2", "height 3", "height 4", 
	"height 5", "height 6", "height 7", "height 8", "height 9 \n(beam hitting near top of sample)"], smooth=True, show=False)
plot_runs(runs_m18_turn, title="M18 Turn in air, adjusting sample height", 
	labels=["height 1 \n(beam hitting center of sample)", "height 2", "height 3", "height 4", 
	"height 5", "height 6", "height 7 \n(beam hitting near top of sample)"], smooth=True, show=False)
plot_runs(runs_old_sample, title="Old Sample in air, adjusting sample height", 
	labels=["height 1 \n(beam hitting center of sample)", "height 2", "height 3", "height 4", 
	"height 5", "height 6", "height 7 \n(beam hitting near top of sample)"], smooth=True, show=False)

plt.show()
