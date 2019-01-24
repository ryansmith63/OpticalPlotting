import numpy as np
import matplotlib.pyplot as plt

path = "/Users/leehagaman/Desktop/OpticalPlotting/"
#path = "C:\\Users\\swkra\\OneDrive\\Documents\\GitHub\\OpticalPlotting\\"

class Run:
	def __init__(self, filename, mirror_filenames_and_angles=[]):
		file = open(path + filename)
		lines = file.readlines()
		data = np.loadtxt(path + filename, skiprows=12)
		self.angles = [datum[0] for datum in data]
		bkg = 70 # constant background from e.g. dark rate; 500 for initial Spectralon 400 nm data, 100 for 500 nm data
		# chosen to be slightly less than lowest rate during background measurement in LXe
		# Have been using 70 for 1st/2nd run LXe data at 178 nm, 150 for vacuum data at 178/175 nm
		self.intensities = [datum[1]-bkg for datum in data]
		self.intensity_std = [datum[2] for datum in data]

		self.date_time = lines[1][15:]
		self.name = lines[2][6:]
		self.description = lines[3][13:]
		self.substance = lines[4][11:]
		self.sampleid = lines[5][11:-1]
		self.preampsensitivity = float(lines[6][20:-1]) * 1e-9
		self.incidentangle = float(lines[7][16:-1])
		self.uncorrectedincidentpower = float(lines[8][16:-1])
		self.incidentpower = self.uncorrectedincidentpower * incident_power_factor_function(mirror_filenames_and_angles)(self.incidentangle)
		self.wavelength = float(lines[9][12:-1])
		self.temperature = float(lines[10][13:-1])
		self.pressure = float(lines[11][10:-1])
		self.relative_intensities = [intensity*intensity_factor(self.incidentpower) for intensity in self.intensities]
		const_err = 300 # error to add to std from e.g. error on background; using 300 for 1st/2nd run LXe data at 178 nm, 100 for vacuum at 178/175 nm
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

# this uses a mirror test to see how much of the beam power is lost at a certain incident angle, presumably
# due to the beam missing part of the sample and hitting the sample rack or back of the chamber. So this should be
# a different file for each incident angle and change in alignment of the sample rack. This file must have an accurate
# incident power stored in the file.
def incident_power_factor(mirror_incident_power, mirror_filename):
	mirror_data = np.loadtxt(path + mirror_filename, skiprows=12)

	mirror_intensities = [datum[1] for datum in mirror_data]
	max_mirror_intensity = np.max(mirror_intensities)

	return max_mirror_intensity / mirror_incident_power


# takes in an incident power file and a list of pairs of incident angles and mirror test filenames,
# returns a list of pairs of incident angles and incident power factors
def incident_power_factor_list(mirror_incident_power_filename, mirror_incident_angle_and_filename_list):
	mirror_incident_power_data = np.loadtxt(path + mirror_incident_power_filename, skiprows=12)
	mirror_incident_power_intensities = [datum[1] for datum in mirror_incident_power_data]
	mirror_incident_power = np.max(mirror_incident_power_intensities)

	lst = []

	for i in range(len(mirror_incident_angle_and_filename_list)):
		lst.append([mirror_incident_angle_and_filename_list[i][0], 
			incident_power_factor(mirror_incident_power, mirror_incident_angle_and_filename_list[i][1])])
	return lst

# retuns a function which gives incident power factor as a function of incident angle by interpolating between the points in incident_power_factor_list
def incident_power_factor_function(mirror_filenames_and_angles):
	if mirror_filenames_and_angles == []:
		return lambda x : 1
	mirror_incident_power_filename = mirror_filenames_and_angles[0]
	mirror_incident_angle_and_filename_list = mirror_filenames_and_angles[1]

	lst = incident_power_factor_list(mirror_incident_power_filename, mirror_incident_angle_and_filename_list)
	angles = [pair[0] for pair in lst]
	factors = [pair[1] for pair in lst]
	def ret(theta):
		if theta <= np.min(angles):
			return factors[0]
		if theta >= np.max(angles):
			return factors[-1]
		for i in range(len(lst)):
			if angles[i] <= theta < angles[i + 1]:
				return factors[i] + (angles[i + 1] - angles[i]) * (factors[i + 1] - factors[i]) / (angles[i + 1] - angles[i]) # interpolating
		print("error, failed to interpolate")
		return 0 #failed to interpolate, should never happen
	return ret


def intensity_factor(incidentpower):
	
	# Angular diameter is arcsin(0.385"/5.4") = 4.1 deg; to get to angular radius of 2.44 for same aperture, distance must be ~4.5" (~4.35" for 2.54)
	# Remeasured, got 0.385-0.390" aperture, ~5.2" distance = 4.3 degrees (2.15 degrees radius)
	# Spectralon calibrated reflectance: 0.9745 at 400 nm, 0.9784 at 500 nm; index is listed as 1.35 (no wavelength specified...)
	# fit to data matches this for angular_radius of 2.67 (400 nm), 2.68 (500 nm) when using bkg of 50, err of 50
	# for 400 nm, bkg of 500, err of 100 (closer match to background measurement at 405 nm, similar power), get 2.44
	# for 500 nm, bkg of 100, err of 100 (closer match to bkg at 500 nm when scaled by power, but taken a month later), get 2.61
	# for 500 nm, bkg of 150, err of 100 (background at 405 nm from same day, but scaled by power), get 2.54
	photodiode_angular_radius = 2.44 #degrees; old plots used 4.0 incorrectly
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





