import matplotlib.pyplot as plt
import numpy as np
from file_reader import Run
from plotting import plot_runs
#TODO compare with 30 rather than with power measurement, don't change within 5 percent
mirror_path = "Vacuum measurements after 3rd xenon run/Jan 9-12/Mirror alignment/Blue height 2/"

# must be in ascending angle order
mirror_filenames_and_angles = [
[30., mirror_path + "2019_01_10__16_03_43.txt"],
[45., mirror_path + "2019_01_10__16_01_50.txt"],
[52., mirror_path + "2019_01_10__16_00_07.txt"],
[60., mirror_path + "2019_01_10__15_58_07.txt"],
[67., mirror_path + "2019_01_10__15_56_24.txt"],
[75., mirror_path + "2019_01_10__15_53_05.txt"]
]

path = "Vacuum measurements after 3rd xenon run/Jan 9-12/LUX/178nm/"

v_75 = Run(path + "2019_01_11__11_31_13.txt")
v_67 = Run(path + "2019_01_11__11_36_08.txt")
v_60 = Run(path + "2019_01_11__11_40_55.txt")
v_52 = Run(path + "2019_01_11__11_45_38.txt")
v_45 = Run(path + "2019_01_11__11_50_42.txt")
v_30 = Run(path + "2019_01_11__11_55_48.txt")

v_75c = Run(path + "2019_01_11__11_31_13.txt", mirror_filenames_and_angles)
v_67c = Run(path + "2019_01_11__11_36_08.txt", mirror_filenames_and_angles)
v_60c = Run(path + "2019_01_11__11_40_55.txt", mirror_filenames_and_angles)
v_52c = Run(path + "2019_01_11__11_45_38.txt", mirror_filenames_and_angles)
v_45c = Run(path + "2019_01_11__11_50_42.txt", mirror_filenames_and_angles)
v_30c = Run(path + "2019_01_11__11_55_48.txt", mirror_filenames_and_angles)

all_runs = [v_30, v_45, v_52, v_60, v_67, v_75]
all_runsc = [v_30c, v_45c, v_52c, v_60c, v_67c, v_75c]

labels=["30 degrees", "45 degrees", "52 degrees", "60 degrees", "67 degrees", "75 degrees"]

plot_runs(all_runs, title="not corrected", include_legend=False)
plot_runs(all_runsc, title="corrected using mirror data", include_legend=False)

plt.show()

