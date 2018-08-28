import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import dblquad

def fit_gaussian(x, y):
	def func(x, A, B, sigma, mu):
		return np.abs(A) * np.exp(-(x - mu) * (x - mu) / (2 * sigma * sigma)) + np.abs(B) * np.cos(x * np.pi / 180.)
	return curve_fit(func, np.array(x), np.array(y))[0]

def fit_gaussian_with_mu(x, y, mu):
	def func(x, A, B, sigma):
		return np.abs(A) * np.exp(-(x - mu) * (x - mu) / (2 * sigma * sigma)) + np.abs(B) * np.cos(x * np.pi / 180.)
	ret = curve_fit(func, np.array(x), np.array(y))[0]
	ret = np.append(ret, mu)
	return ret

def plot_gaussian_fit_by_params(xmin, xmax, params, label):
	def func(x, A, B, sigma, mu):
		return np.abs(A) * np.exp(-(x - mu) * (x - mu) / (2 * sigma * sigma)) + np.abs(B) * np.cos(x * np.pi / 180.)
	x = np.linspace(xmin, xmax, 1000)
	y = [func(x_, params[0], params[1], params[2], params[3]) for x_ in x]
	plt.plot(x, y, label=label)

def gaussian_get_params(run, mu=1000):
	x = run.angles
	y = run.relative_intensities

	if mu == 1000: #placeholder value, means that mu is unspecified
		params = fit_gaussian(x, y)
	else:
		params = fit_gaussian_with_mu(x, y, mu)
	return params


def gaussian_total_reflectance(params):
	A, B, sigma, mu = params
	# r is angular distance from specular reflection to viewing angle
	# theta is angular distance from sample normal to viewing angle
	def get_r(theta1, phi1, theta2, phi2): # angular distance
		v1 = [np.sin(theta1) * np.cos(phi1), np.sin(theta1) * np.sin(phi1), np.cos(theta1)]
		v2 = [np.sin(theta2) * np.cos(phi2), np.sin(theta2) * np.sin(phi2), np.cos(theta2)]
		return 180. / np.pi * np.arccos(np.dot(v1, v2))

	def f(theta, phi):
		r = get_r(theta, phi, mu, 0)
		return np.abs(A) * np.exp(-(r - mu) * (r - mu) / (2 * sigma * sigma)) + np.abs(B) * np.cos(theta * np.pi / 180.)

	def f_sin_theta(theta, phi):
		return f(theta, phi) * np.sin(theta)

	return dblquad(f_sin_theta, 0, 2 * np.pi, lambda x: 0, lambda x: np.pi / 2.)[0]

def gaussian_get_total_reflection(run):
	return gaussian_total_reflectance(get_params(run))
