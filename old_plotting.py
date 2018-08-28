import numpy as np

def make_data_by_run(filename, lower_cutoff=-90, upper_cutoff=90, intensity_factor=1):

    data = np.loadtxt(filename, skiprows=1)
    data_by_run = []

    for i in range(len(data)):
        # if the x_data decreases (new run)
        if i == 0 or data[i][0] < data[i - 1][0]:
            data_by_run.append([[], []])
        if lower_cutoff <= data[i][0] <= upper_cutoff:
            data_by_run[-1][0].append(np.round(data[i][0], 3))
            data_by_run[-1][1].append(intensity_factor * data[i][1])
    return data_by_run