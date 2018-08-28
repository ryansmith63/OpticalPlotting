import numpy as np
import matplotlib.pyplot as plt
from file_reader import Run
from plotting import plot_runs, plot_fit, plot_fit_by_params, get_total_reflection

path = "sample_height_tests/"

run = Run(path + "2018_05_31__13_04_27.txt") # random M18_turn test

plot_runs(run, title="NXT85 in air, adjusting sample height", labels=["data"])

#plot_fit_by_params(0., 85., [100., 100., 10., 50.])

plot_fit(run, mu=49.)

print(get_total_reflection(run))

plt.legend(loc=3)
plt.show()
