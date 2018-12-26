import numpy as np
import matplotlib.pyplot as plt

#path = "/Users/leehagaman/Desktop/OpticalPlotting/"
path = "C:\\Users\\swkra\\OneDrive\\Documents\\GitHub\\OpticalPlotting\\"

class Run:
	def __init__(self, filename):
		file = open(path + filename)
		lines = file.readlines()
		data = np.loadtxt(path + filename, skiprows=12)
		self.angles = [datum[0] for datum in data]
		bkg = 100 # constant background from e.g. dark rate; 500 for initial Spectralon 400 nm data, 100 for 500 nm data
		# chosen to be slightly less than lowest rate during background measurement in LXe
		self.intensities = [datum[1]-bkg for datum in data]
		self.intensity_std = [datum[2] for datum in data]

		self.date_time = lines[1][15:]
		self.name = lines[2][6:]
		self.description = lines[3][13:]
		self.substance = lines[4][11:]
		self.sampleid = lines[5][11:-1]
		self.preampsensitivity = float(lines[6][20:-1]) * 1e-9
		self.incidentangle = float(lines[7][16:-1])
		self.incidentpower = float(lines[8][16:-1])
		self.wavelength = float(lines[9][12:-1])
		self.temperature = float(lines[10][13:-1])
		self.pressure = float(lines[11][10:-1])
		self.relative_intensities = [intensity*intensity_factor(self.incidentpower) for intensity in self.intensities]
		const_err = 100 # error to add to std from e.g. error on background
		self.relative_std = [(std+const_err)*intensity_factor(self.incidentpower) for std in self.intensity_std]

		self.rot_angles = [180. - self.incidentangle - a for a in self.angles]

		if self.substance[:3].lower() == "vac" or self.substance[:3].lower() == "air":
			self.n = 1.
		if self.substance[:3].lower() == "lxe":
			self.n = 1.69
			# from https://arxiv.org/pdf/physics/0307044

		# independent_variables_array is a list where each element is of the form [theta_r_in_degrees, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization]
		self.independent_variables_array = [[angle, 0, self.incidentangle, self.n, 0.5] for angle in self.angles]

	def change_theta_i(self, new_theta_i):
		# Changes incidentangle, angles, and angles in independent_variables_array; does not change rot_angles (which should be fixed)
		delta_theta = new_theta_i-self.incidentangle
		self.angles = [angle-delta_theta for angle in self.angles] # If incident angle increases, angles relative to normal decrease
		self.independent_variables_array = [[list[0]-delta_theta]+list[1:] for list in self.independent_variables_array]
		self.incidentangle = new_theta_i

def intensity_factor(incidentpower):
	
	# Angular diameter is arcsin(0.385"/5.4") = 4.1 deg; to get to angular radius of 2.44 for same aperture, distance must be ~4.5" (~4.35" for 2.54)
	# Remeasured, got 0.385-0.390" aperture, ~5.2" distance = 4.3 degrees (2.15 degrees radius)
	# Spectralon calibrated reflectance: 0.9745 at 400 nm, 0.9784 at 500 nm; index is listed as 1.35 (no wavelength specified...)
	# fit to data matches this for angular_radius of 2.67 (400 nm), 2.68 (500 nm) when using bkg of 50, err of 50
	# for 400 nm, bkg of 500, err of 100 (closer match to background measurement at 405 nm, similar power), get 2.44
	# for 500 nm, bkg of 100, err of 100 (closer match to bkg at 500 nm when scaled by power, but taken a month later), get 2.61
	# for 500 nm, bkg of 150, err of 100 (background at 405 nm from same day, but scaled by power), get 2.54
	photodiode_angular_radius = 2.54 #degrees; old plots used 4.0 incorrectly
	photodiode_angular_size = photodiode_angular_radius * np.pi / 180.

	photodiode_solid_angle = np.pi * np.power(photodiode_angular_size, 2)

	# product of this and measured rate is (rate/str)/rate_i
	# intensity_factor * rate = (rate / photodiode_solid_angle) / flux_i
	intensity_factor = 1. / (photodiode_solid_angle * incidentpower)
	
	return intensity_factor

def get_independent_variables_and_relative_intensities(runs):
	if type(runs) != type([]): # if runs is a single run
		runs = [runs] # make it a list
		
	independent_variables_array = []
	relative_intensities = []
	relative_stds = []

	for run in runs:
		for i in range(len(run.relative_intensities)):
			independent_variables_array.append(run.independent_variables_array[i])
			relative_intensities.append(run.relative_intensities[i])
			relative_stds.append(run.relative_std[i])
	return [independent_variables_array, relative_intensities, relative_stds]





