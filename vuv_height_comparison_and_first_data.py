import matplotlib.pyplot as plt
from file_reader import Run
from plotting import plot_runs

log = True
voltage = False

nxpath_off = "vuv_height_comparison_and_first_data/807NX turn/One eighth inch above center of sample/"
nx30_off = Run(nxpath_off + "2018_08_29__16_49_57.txt")
nx45_off = Run(nxpath_off + "2018_08_29__16_55_41.txt")
nx60_off = Run(nxpath_off + "2018_08_29__17_01_55.txt")
# no 75 degree measurement taken

nxpath = "vuv_height_comparison_and_first_data/807NX turn/Center of sample/"
nx30 = Run(nxpath + "2018_08_29__15_57_32.txt")
nx45 = Run(nxpath + "2018_08_29__16_02_55.txt")
nx60 = Run(nxpath + "2018_08_29__16_08_16.txt")
nx75 = Run(nxpath + "2018_08_29__16_13_39.txt")

m18path = "vuv_height_comparison_and_first_data/M18 turn/center of sample/"
m1830 = Run(m18path + "2018_08_30__15_16_18.txt")
m1845 = Run(m18path + "2018_08_30__15_22_02.txt")
m1860 = Run(m18path + "2018_08_30__15_29_09.txt")
m1875 = Run(m18path + "2018_08_30__15_35_27.txt")

m18path_off = "vuv_height_comparison_and_first_data/M18 turn/old/Good data/"
m1830_off = Run(m18path_off + "2018_08_28__17_26_14.txt")
m1845_off = Run(m18path_off + "2018_08_28__17_34_54.txt")
m1860_off = Run(m18path_off + "2018_08_28__17_43_46.txt")
m1875_off = Run(m18path_off + "2018_08_28__17_52_54.txt")

nxtpath = "vuv_height_comparison_and_first_data/NXT85 turn/center of sample/"
nxt30 = Run(nxtpath + "2018_08_30__14_13_18.txt")
nxt45 = Run(nxtpath + "2018_08_30__14_19_14.txt")
nxt60 = Run(nxtpath + "2018_08_30__14_24_51.txt")
nxt75 = Run(nxtpath + "2018_08_30__14_30_46.txt")

"""
plt.figure()
plot_runs([nx30, nx45, nx60], title="height comparison, 807NX turn, 175nm", smooth=True, figure=False, show=False, voltage=voltage, log=log,
	labels=["30 degrees centered", "45 degrees centered", "60 degrees centered"])
plot_runs([nx30_off, nx45_off, nx60_off], title="height comparison, 807NX turn, 175nm", smooth=True, figure=False, show=False, voltage=voltage, linestyle="--", log=log,
	labels=["30 degrees 1/8\" above center", "45 degrees 1/8\" above center", "60 degrees 1/8\" above center"])
"""

plt.figure()
plot_runs([m1830, m1845, m1860, m1875], title="height comparison, M18 turn, 175nm", smooth=True, figure=False, show=False, voltage=voltage, log=log,
	labels=["30 degrees centered", "45 degrees centered", "60 degrees centered", "75 degrees centered"])
plot_runs([m1830_off, m1845_off, m1860_off, m1875_off], title="height comparison, M18 turn, 175nm", smooth=True, figure=False, show=False, voltage=voltage, linestyle="--", log=log,
	labels=["30 degrees 1/8\" above center", "45 degrees 1/8\" above center", "60 degrees 1/8\" above center", "75 degrees 1/8\" above center"])

plot_runs([m1830, m1845, m1860, m1875], title="M18 turn, 175 nm", smooth=True, show=False, voltage=voltage, log=log,
	labels=["30 degrees", "45 degrees", "60 degrees", "75 degrees"])

plot_runs([nxt30, nxt45, nxt60, nxt75], title="NXT85 turn, 175 nm", smooth=True, show=False, voltage=voltage, log=log,
	labels=["30 degrees", "45 degrees", "60 degrees", "75 degrees"])

plt.show()
