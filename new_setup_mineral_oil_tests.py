import matplotlib.pyplot as plt
from file_reader import Run
from plotting import plot_runs

path = "new_setup_mineral_oil_tests/"

_M18_turn_45 = [
	Run(path + "2018_06_25__12_40_19.txt"),
	Run(path + "2018_06_25__12_44_54.txt"),
	Run(path + "2018_06_25__12_49_03.txt")
]

_M18_turn_60 = [
	Run(path + "2018_06_25__12_56_14.txt"),
	Run(path + "2018_06_25__13_00_24.txt"),
	Run(path + "2018_06_25__13_04_31.txt")
]

_M18_skived_60 = [
	Run(path + "2018_06_25__13_21_55.txt"),
	Run(path + "2018_06_25__13_25_51.txt"),
	Run(path + "2018_06_25__13_30_07.txt")
]

_M18_skived_45 = [
	Run(path + "2018_06_25__13_35_43.txt"),
	Run(path + "2018_06_25__13_39_21.txt"),
	Run(path + "2018_06_25__13_43_02.txt")
]

_M18_polished_45 = [
	Run(path + "2018_06_25__13_51_03.txt"),
	Run(path + "2018_06_25__13_54_51.txt"),
	Run(path + "2018_06_25__13_58_33.txt"),
	Run(path + "2018_06_25__14_02_37.txt")
]

_M18_polished_60 = [
	Run(path + "2018_06_25__14_08_19.txt"),
	Run(path + "2018_06_25__14_13_18.txt"),
	Run(path + "2018_06_25__14_17_21.txt")
]

_M17_turn_60 = [
	Run(path + "2018_06_25__16_28_04.txt"),
	Run(path + "2018_06_25__16_32_18.txt"),
	Run(path + "2018_06_25__16_36_18.txt")
]

_M17_turn_45 = [
	Run(path + "2018_06_25__16_41_10.txt"),
	Run(path + "2018_06_25__16_45_05.txt"),
	Run(path + "2018_06_25__16_48_57.txt")
]

_NXT85_turn_45 = [
	Run(path + "2018_06_26__10_49_37.txt"),
	Run(path + "2018_06_26__10_53_55.txt"),
	Run(path + "2018_06_26__11_00_51.txt")
]

_NXT85_turn_60 = [
	Run(path + "2018_06_26__11_05_45.txt"),
	Run(path + "2018_06_26__11_09_27.txt"),
	Run(path + "2018_06_26__11_15_26.txt")
]

_807NX_turn_60 = [
	Run(path + "2018_06_26__11_26_29.txt"),
	Run(path + "2018_06_26__11_35_07.txt"),
	Run(path + "2018_06_26__11_39_32.txt")
]

_807NX_turn_45 = [
	Run(path + "2018_06_26__11_47_29.txt"),
	Run(path + "2018_06_26__11_51_14.txt"),
	Run(path + "2018_06_26__11_57_17.txt")
]

_old_sample_45 = [
	Run(path + "2018_06_26__12_11_45.txt"),
	Run(path + "2018_06_26__12_18_36.txt"),
	Run(path + "2018_06_26__12_22_10.txt")
]

_old_sample_60 = [
	Run(path + "2018_06_26__12_26_52.txt"),
	Run(path + "2018_06_26__12_30_44.txt"),
	Run(path + "2018_06_26__12_34_17.txt")
]

all_runs_45 = [_M18_turn_45, _M18_skived_45, _M18_polished_45, _M17_turn_45, _NXT85_turn_45, _807NX_turn_45, _old_sample_45]
all_runs_60 = [_M18_turn_60, _M18_skived_60, _M18_polished_60, _M17_turn_60, _NXT85_turn_60, _807NX_turn_60, _old_sample_60]

all_runs = all_runs_45 + all_runs_60

plot_runs(all_runs, title="mineral oil at different sample heights", labels=["M18 turn, 45 degrees", "M18 skived, 45 degrees", "M18 polished, 45 degrees", "M17 turn, 45 degrees", "NXT85 turn, 45 degrees", "807NX turn, 45 degrees", "old sample, 45 degrees", "M18 turn, 60 degrees", "M18 skived, 60 degrees", "M18 polished, 60 degrees", "M17 turn, 60 degrees", "NXT85 turn, 60 degrees", "807NX turn, 60 degrees", "old sample, 60 degrees"], smooth=True, show=False, legend_loc=2)

plot_runs(all_runs_45, title="45 degrees in mineral oil at different sample heights", labels=["M18 turn", "M18 skived", "M18 polished", "M17 turn", "NXT85 turn", "807NX turn", "old sample"], smooth=True, show=False, legend_loc=2)
plot_runs(all_runs_60, title="60 degrees in mineral oil at different sample heights", labels=["M18 turn", "M18 skived", "M18 polished", "M17 turn", "NXT85 turn", "807NX turn", "old sample"], smooth=True, show=False, legend_loc=2)

plot_runs([_M18_turn_45], title="M18 turn at 45 degrees in mineral oil at different sample heights", smooth=True, show=False, include_legend=False)
plot_runs([_M18_turn_60], title="M18 turn at 60 degrees in mineral oil at different sample heights", smooth=True, show=False, include_legend=False)
plot_runs([_M18_skived_45], title="M18 skived at 45 degrees in mineral oil at different sample heights", smooth=True, show=False, include_legend=False)
plot_runs([_M18_skived_60], title="M18 skived at 60 degrees in mineral oil at different sample heights", smooth=True, show=False, include_legend=False)
plot_runs([_M18_polished_45], title="M18 polished at 45 degrees in mineral oil at different sample heights", smooth=True, show=False, include_legend=False)
plot_runs([_M18_polished_60], title="M18 polished at 60 degrees in mineral oil at different sample heights", smooth=True, show=False, include_legend=False)
plot_runs([_M17_turn_45], title="M17 turn at 45 degrees in mineral oil at different sample heights", smooth=True, show=False, include_legend=False)
plot_runs([_M17_turn_60], title="M17 turn at 60 degrees in mineral oil at different sample heights", smooth=True, show=False, include_legend=False)
plot_runs([_NXT85_turn_45], title="NXT85 turn at 45 degrees in mineral oil at different sample heights", smooth=True, show=False, include_legend=False)
plot_runs([_NXT85_turn_60], title="NXT85 turn at 60 degrees in mineral oil at different sample heights", smooth=True, show=False, include_legend=False)
plot_runs([_807NX_turn_45], title="807NX turn at 45 degrees in mineral oil at different sample heights", smooth=True, show=False, include_legend=False)
plot_runs([_807NX_turn_60], title="807NX turn at 60 degrees in mineral oil at different sample heights", smooth=True, show=False, include_legend=False)
plot_runs([_old_sample_45], title="old sample at 45 degrees in mineral oil at different sample heights", smooth=True, show=False, include_legend=False)
plot_runs([_old_sample_60], title="old sample at 60 degrees in mineral oil at different sample heights", smooth=True, show=False, include_legend=False)

plt.show()
