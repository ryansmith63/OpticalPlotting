import numpy as np
import matplotlib.pyplot as plt

path = "/Users/leehagaman/Desktop/OpticalPlotting/"

class Run:
	def __init__(self, filename):
		file = open(path + filename)
		lines = file.readlines()
		data = np.loadtxt(path + filename, skiprows=12)
		self.angles = [datum[0] for datum in data]
		self.intensities = [datum[1] for datum in data]

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
		self.relative_intensities = [relative_intensity(intensity, self.incidentpower) for intensity in self.intensities]

		self.rot_angles = [180. - self.incidentangle - a for a in self.angles]

		

def relative_intensity(intensity, incidentpower):
	
	photodiode_angular_size = 4. * np.pi / 180.

	photodiode_solid_angle = np.pi * np.power(photodiode_angular_size, 2)

	# product of this and measured rate is (rate/str)/rate_i
	# intensity_factor * rate = (rate / photodiode_solid_angle) / flux_i
	intensity_factor = 1. / (photodiode_solid_angle * incidentpower)
	
	return intensity_factor * intensity

