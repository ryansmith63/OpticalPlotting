import matplotlib.pyplot as plt
from file_reader import Run
from plotting import plot_runs

# the path to the txt files
path = "initial_PMT_tests_M18_turn/"

# giving each run of data a variable name

power_405 = Run(path + "2018_08_21__10_16_41.txt")

_45_405 = Run(path + "2018_08_21__10_35_50.txt")

_60_405 = Run(path + "2018_08_21__10_42_54.txt")

_60_175 = Run(path + "2018_08_21__10_51_23.txt")

_45_175 = Run(path + "2018_08_21__10_57_15.txt")

power_175 = Run(path + "2018_08_21__11_02_58.txt")

# this step is done to normalize the data (numbers gotten from manually looking at text file)

max_405 = 135464

max_175 = 640232

# changing the relative intensities to the normalized ones 
# (normally done automatically, but not input into LabVIEW for these tests)
# by default, intensities are plotted, not relative intensities

_45_405.relative_intensities = [x / max_405 for x in _45_405.intensities]
_60_405.relative_intensities = [x / max_405 for x in _60_405.intensities]
_45_175.relative_intensities = [x / max_175 for x in _45_175.intensities]
_60_175.relative_intensities = [x / max_175 for x in _60_175.intensities]

# doing the plotting

plot_runs(power_405, title="Power Measurement, 405 nm", 
	smooth=True, show=False, voltage=True, ylabel="rate (Hz)", rot=True)
plot_runs(power_175, title="Power Measurement, 175 nm", 
	smooth=True, show=False, voltage=True, ylabel="rate (Hz)", rot=True)

plot_runs([_45_405, _60_405], title="M18 turn, 405 nm", 
	labels=["45 degrees", "60 degrees"], smooth=True, show=False, voltage=True, ylabel="rate (Hz)")
plot_runs([_45_175, _60_175], title="M18 turn, 175 nm", 
	labels=["45 degrees", "60 degrees"], smooth=True, show=False, voltage=True, ylabel="rate (Hz)")

plot_runs([_45_405, _60_405, _45_175, _60_175], title="M18 turn", 
	labels=["45 degrees, 405 nm", "60 degrees, 405 nm", "45 degrees, 175 nm", "60 degrees, 175 nm"], smooth=True, show=False, ylabel="relative intensity")

plt.show()
