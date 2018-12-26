# TSTR_fit

# specular follows a torrance-sparrow model
# diffuse follows a trowbridge-reitz distribution of surface normals combined with lambertian diffusion

# described in Reflectance of Polytetrafluoroethylene (PTFE) for Xenon_Scintillation Light, LIP-Coimbra

# called geometrical optical approximation (GOA) model in other Coimbra paper,
# but this doesn't include the correction factor N

import numpy as np
import scipy.optimize
import copy
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.colors

# following equations could be refined for mu !~ mu_0, but teflon has mu~mu_0


# this function is used in calculating the total reflectance
def G_calc(theta_r, phi_r, theta_i, n_0, polarization, parameters):
    rho_L = parameters[0]
    n = parameters[1]
    gamma = parameters[2]

    theta_i_prime = 0.5 * np.arccos(np.cos(theta_i) * np.cos(theta_r) -
                                    np.sin(theta_i) * np.sin(theta_r) * np.cos(phi_r))

    theta_r_prime = 0.5 * np.arccos(np.cos(theta_i) * np.cos(theta_r) -
                                    np.sin(theta_i) * np.sin(theta_r))

    def G_prime(theta):
        return 2 / (1 + np.sqrt(1 + np.power(gamma * np.tan(theta), 2)))

        # this has negative of the inside of the H functions from the paper, I think it was a typo

    G = H(np.pi / 2 - theta_i_prime) * H(np.pi / 2 - theta_r_prime) * \
        G_prime(theta_i) * G_prime(theta_r)

    return G

def F_s(theta_i, n_0, n):
    if np.size(theta_i)>1: return F_s_array(theta_i, n_0, n)
    else: return F_s_single(theta_i, n_0, n)
    
# polarized electric field perpendicular to plane of incidence (vertical)
def F_s_single(theta_i, n_0, n):
    if np.abs(n_0 / n * np.sin(theta_i)) > 1:
        # total internal reflection
        return 1
    theta_t = np.arcsin(n_0 / n * np.sin(theta_i))
    return np.power((n_0 * np.cos(theta_i) - n * np.cos(theta_t)) /
                    (n_0 * np.cos(theta_i) + n * np.cos(theta_t)), 2)
    
# polarized electric field perpendicular to plane of incidence (vertical)
# designed to work when theta_i is an array; more efficient for arrays, but much 
# slower for single values due to masked arrays being slow
def F_s_array(theta_i, n_0, n):
    mask = np.abs(n_0 / n * np.sin(theta_i)) > 1
    #print(type(mask))
    theta_i_mask = np.ma.masked_array(theta_i, mask)
    theta_t = np.arcsin(n_0 / n * np.sin(theta_i_mask))
    F_calc = np.power((n_0 * np.cos(theta_i) - n * np.cos(theta_t)) /
                    (n_0 * np.cos(theta_i) + n * np.cos(theta_t)), 2)
    #print(F_calc)
    return np.ma.filled(F_calc,fill_value=1)

def F_p(theta_i, n_0, n):
    if np.size(theta_i)>1: return F_p_array(theta_i, n_0, n)
    else: return F_p_single(theta_i, n_0, n)
    
# polarized electric field parallel to plane of incidence (horizontal)
def F_p_single(theta_i, n_0, n):
    if np.abs(n_0 / n * np.sin(theta_i)) > 1:
        # total internal reflection
        return 1
    theta_t = np.arcsin(n_0 / n * np.sin(theta_i))
    return np.power((n_0 * np.cos(theta_t) - n * np.cos(theta_i)) /
                    (n_0 * np.cos(theta_t) + n * np.cos(theta_i)), 2)

# polarized electric field parallel to plane of incidence (horizontal)
# designed to work when theta_i is an array; more efficient for arrays, but much 
# slower for single values due to masked arrays being slow
def F_p_array(theta_i, n_0, n):
    # Have to check for arcsin errors; using masked array to handle multiple theta_i values
    mask = np.abs(n_0 / n * np.sin(theta_i)) > 1
    #print(type(mask))
    theta_i_mask = np.ma.masked_array(theta_i, mask)
    theta_t = np.arcsin(n_0 / n * np.sin(theta_i_mask))
    #print(type(theta_t))
    F_calc = np.power((n_0 * np.cos(theta_t) - n * np.cos(theta_i)) /
                    (n_0 * np.cos(theta_t) + n * np.cos(theta_i)), 2)
    #print(type(F_calc))
    return np.ma.filled(F_calc,fill_value=1)


def F_unpolarized(theta_i, n_0, n):
    return 0.5 * (F_s(theta_i, n_0, n) + F_p(theta_i, n_0, n))


def F(theta_i, n_0, n, polarization):
    return polarization * F_p(theta_i, n_0, n) + (1 - polarization) * F_s(theta_i, n_0, n)


# heaviside step function
def H(x):
    return 0.5 * (np.sign(x) + 1)


# BRIDF as described in paper
# takes an array, easy to plot
# parameters has form [rho_L, n, gamma]; optional fourth parameter K
# All angles are in degrees
# average_angle sets aperture size (diameter) to average over (no avg if 0)
# precision sets grid spacing for 2D averaging (1D average used if < 0)
def BRIDF_plotter(theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, parameters, average_angle=0, precision=-1, sigma_theta_i=2):
    phi_r = phi_r_in_degrees * np.pi / 180
    theta_i = theta_i_in_degrees * np.pi / 180
    precision_rad = precision * np.pi / 180
    if precision < 0: # Skip 2D averaging
        return_array = []
        for theta_r_in_degrees in theta_r_in_degrees_array:
            theta_r = np.pi * theta_r_in_degrees / 180
            return_array.append(BRIDF(theta_r, phi_r, theta_i, n_0, polarization, parameters, precision=precision_rad))
        if average_angle == 0:
            return return_array
        else: # Just average adjacent points in plane
            return average_by_angle(theta_r_in_degrees_array, return_array, average_angle)
    else: # Do a full 2D average across both theta_r and phi_r, with grid spacing given by precision
        avgs = []
        log_rho_L = np.log(parameters[0])
        log_n_minus_one = np.log(parameters[1]-1)
        log_gamma = np.log(parameters[2])
        if len(parameters)>3: 
            if parameters[3]>0: log_K = np.log(parameters[3])
            else: log_K = None
        else: log_K = None
        for theta_r_in_degrees in theta_r_in_degrees_array: # For each point, calculate the BRIDF along a grid nearby and average it
            avgs.append(average_BRIDF([theta_r_in_degrees,phi_r_in_degrees,theta_i_in_degrees,n_0,polarization], log_rho_L, log_n_minus_one, log_gamma, log_K, average_angle=average_angle, precision=precision, sigma_theta_i=sigma_theta_i))
        return avgs

# Takes a set of viewing angles and intensities and returns the average, for each point, of adjacent points w/in +/- angle/2
# Assumes angles are sorted; if arrays include multiple runs, will not average across runs
def average_by_angle(theta_r_in_degrees_array, intensity_array, angle):
    #print("avg angle: "+str(angle))
    #print("theta array: "+str(theta_r_in_degrees_array))
    grouped = [[] for i in intensity_array]
    angles = [[] for i in intensity_array]
    for i in range(len(intensity_array)):
        for j in range(i,len(intensity_array)):
            if np.abs(theta_r_in_degrees_array[j] - theta_r_in_degrees_array[i]) < angle / 2.:
                grouped[i].append(intensity_array[j])
                angles[i].append(theta_r_in_degrees_array[j])
            else: break # if we've left the averaging angle range, stop looking forward to avoid combining runs
        for j in range(i,0,-1):
            if np.abs(theta_r_in_degrees_array[j] - theta_r_in_degrees_array[i]) < angle / 2.:
                grouped[i].append(intensity_array[j])
                angles[i].append(theta_r_in_degrees_array[j])
            else: break # if we've left the averaging angle range, stop looking backward to avoid combining runs
    averaged = [np.average(grouped[i]) for i in range(len(intensity_array))]
    return averaged

# Not updated yet to allow averaging
def BRIDF_radiance_plotter(theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization, parameters):
    phi_r = phi_r_in_degrees * np.pi / 180
    theta_i = theta_i_in_degrees * np.pi / 180
    return_array = []
    for theta_r_in_degrees in theta_r_in_degrees_array:
        theta_r = np.pi * theta_r_in_degrees / 180
        return_array.append(BRIDF(theta_r, phi_r, theta_i, n_0, polarization, parameters) / np.cos(theta_r))
    return return_array

# Not updated yet to allow averaging
def BRIDF_specular_plotter(theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees,
                           n_0, polarization, parameters):
    phi_r = phi_r_in_degrees * np.pi / 180
    theta_i = theta_i_in_degrees * np.pi / 180
    return_array = []
    for theta_r_in_degrees in theta_r_in_degrees_array:
        theta_r = np.pi * theta_r_in_degrees / 180
        return_array.append(BRIDF_specular(theta_r, phi_r, theta_i, n_0, polarization, parameters))
    return return_array

# Not updated yet to allow averaging
def BRIDF_diffuse_plotter(theta_r_in_degrees_array, phi_r_in_degrees, theta_i_in_degrees,
                          n_0, polarization, parameters):
    phi_r = phi_r_in_degrees * np.pi / 180
    theta_i = theta_i_in_degrees * np.pi / 180
    return_array = []
    for theta_r_in_degrees in theta_r_in_degrees_array:
        theta_r = np.pi * theta_r_in_degrees / 180
        return_array.append(BRIDF_diffuse(theta_r, phi_r, theta_i, n_0, polarization, parameters))
    return return_array

# Analytical calculation of BRIDF function, using Coimbra model
# If parameters includes K, specular spike is included; precision sets length scale of delta function in specular spike
# All angles should be in radians at this stage
def BRIDF_pair(theta_r, phi_r, theta_i, n_0, polarization, parameters, precision=-1):

    rho_L = parameters[0]
    n = parameters[1]
    gamma = parameters[2]

    # Local angle, relative to microfacet; from page 85 of Claudio's thesis
    theta_prime = 0.5 * np.arccos(np.cos(theta_i) * np.cos(theta_r) -
        np.sin(theta_i) * np.sin(theta_r) * np.cos(phi_r))

    theta_i_prime = theta_prime
    theta_r_prime = theta_prime

    # Trowbridge-Reitz micro-facet angle distribution
    def P(alpha_):
        return np.power(gamma, 2) / \
               (np.pi * np.power(np.cos(alpha_), 4) *
                np.power(np.power(gamma, 2) + np.power(np.tan(alpha_), 2), 2))

    # Check probably unnecessary; if cos errors pop up, will have to fix to work for arrays 
    # if theta_i == theta_r:
        # alpha_specular = 0
    # else:
        # alpha_specular = np.arccos((np.cos(theta_i) + np.cos(theta_r)) / (2 * np.cos(theta_i_prime)))
    # Micro-facet angle for reflection into angle theta_r, also from p. 85
    alpha_specular = np.arccos((np.cos(theta_i) + np.cos(theta_r)) / (2 * np.cos(theta_i_prime)))

    P_ = P(alpha_specular)

    # Wolff correction factor, p. 114 of thesis
    # Operating on masked arrays is slow; only use them if needed (inputs are arrays)
    if np.size(theta_r) > 1 or np.size(theta_i) > 1 or np.size(phi_r) > 1: # Have to use masks
        mask = np.abs(n_0 / n * np.sin(theta_r)) > 1
        theta_r_mask = np.ma.masked_array(theta_r, mask)
        W = (1 - F(theta_i, n_0, n, polarization)) * (1 - F_unpolarized(np.arcsin(n_0 / n * np.sin(theta_r_mask)), n, n_0))
        W = np.ma.filled(W, fill_value=0)
    else:
        if np.abs(n_0 / n * np.sin(theta_r)) > 1:
            # the W expression would have an arcsin error
            W = 0
        else:
            # There was a typo here (?), I changed by moving parentheses and taking a reciprocal
            W = (1 - F(theta_i, n_0, n, polarization)) * (1 - F_unpolarized(np.arcsin(n_0 / n * np.sin(theta_r)), n, n_0))

    # Oren-Nayar correction factor(s) (1-A+B)*cos(theta_i);
    # From p. 8 of arXiv:0910.1056v1 (Reflectance of PTFE for Xenon Scintillation Light), using Trowbridge-Reitz
    # (For comparizon, Claudio's thesis p. 117 has Oren-Nayar correction factor, N, for Torrance-Sparrow)
    A = 0.5 * np.power(gamma, 2) / (np.power(gamma, 2) + 0.92)

    theta_m = np.minimum(theta_i, theta_r)

    theta_M = np.maximum(theta_i, theta_r)

    B = 0.45 * np.power(gamma, 2) / (np.power(gamma, 2) + 0.25) * \
        H(np.cos(phi_r)) * np.cos(phi_r) * np.sin(theta_M) * np.tan(theta_m)

    # Shadowing and masking, for Trowbridge-Reitz distribution; p. 108 of thesis
    def G_prime(theta):
        return 2 / (1 + np.sqrt(1 + np.power(gamma * np.tan(theta), 2)))

    # this has negative of the inside of the H functions from the paper, I think it was a typo
    G = H(np.pi / 2 - theta_i_prime) * H(np.pi / 2 - theta_r_prime) * \
        G_prime(theta_i) * G_prime(theta_r)

    # Use Fresnel factor relative to local angle theta_i_prime
    F_ = F(theta_i_prime, n_0, n, polarization)

    # add specular spike, if parameters length >3
    specular_delta=0
    if len(parameters)>3:    
        K = parameters[3]
        # this is where it differs from semi-empirical fit, Lambda was slightly different 
        # see Silva Oct. 2009 "A model of the reflection distribution in the vacuum ultra violet region"
        C = specular_spike_norm(theta_r, theta_i, K)
        # integrated over the photodiode solid angle, the delta functions go away; make sure spike isn't double counted
        is_theta_spec = np.logical_or(np.abs(theta_r - theta_i) < precision/2., theta_r==theta_i+precision/2.) 
        is_phi_spec = np.logical_or(np.abs(phi_r) < precision/2.,phi_r==precision/2.) 
        is_specular = np.logical_and(is_theta_spec, is_phi_spec)
        #print(is_specular)
        specular_delta = is_specular/precision**2 # only one point in grid is non-zero; magnitude scales w/ grid size to preserve average
    else:
        C = 0
    

    # Semi-empirical formula from p. 121
    # specular component is split into specular lobe w/ normalization=1-C and specular spike w/ normalization=C
    return [(1-C) * F_ * G * P_ / (4 * np.cos(theta_i)) + C * F_ * G * specular_delta,
        rho_L / np.pi * W * (1 - A + B) * np.cos(theta_r)]


def BRIDF_specular(theta_r, phi_r, theta_i, n_0, polarization, parameters, precision=-1):
    return BRIDF_pair(theta_r, phi_r, theta_i, n_0, polarization, parameters, precision=precision)[0]


def BRIDF_diffuse(theta_r, phi_r, theta_i, n_0, polarization, parameters, precision=-1):
    return BRIDF_pair(theta_r, phi_r, theta_i, n_0, polarization, parameters, precision=precision)[1]


def BRIDF(theta_r, phi_r, theta_i, n_0, polarization, parameters, precision=-1):
    pair = BRIDF_pair(theta_r, phi_r, theta_i, n_0, polarization, parameters, precision=precision)
    return pair[0] + pair[1]


def specular_spike_norm(theta_r, theta_i, K):
    return np.exp(-K / 2 * (np.cos(theta_i) + np.cos(theta_r))) # Version from Coimbra paper
    #return np.exp(-K *np.cos(theta_i)) # Version from Claudio's thesis
    
# Calculates BRIDF as a function of parameters; inputs are in degrees, converted to radians; can be arrays
# independent variables has the form [theta_r_in_degrees, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization]
def unvectorized_fitter(independent_variables, log_rho_L, log_n_minus_one, log_gamma, log_K=None, precision=-1):
    theta_r = independent_variables[0] * np.pi / 180
    phi_r = independent_variables[1] * np.pi / 180
    theta_i = independent_variables[2] * np.pi / 180
    n_0 = independent_variables[3]
    polarization = independent_variables[4]
    precision = precision * np.pi / 180
    rho_L = np.exp(log_rho_L)
    n = np.exp(log_n_minus_one) + 1
    gamma = np.exp(log_gamma)
    parameters = [rho_L, n, gamma]
    if log_K is not None:
        parameters.append(np.exp(log_K))
    return BRIDF(theta_r, phi_r, theta_i, n_0, polarization, parameters, precision=precision)


# independent variables without angle has the form [theta_r_in_degrees, phi_r_in_degrees, n_0, polarization]
def unvectorized_fitter_with_angle(independent_variables_without_angle, theta_i, log_rho_L, log_n_minus_one, log_gamma, log_K=None, precision=-1):
    theta_r = independent_variables_without_angle[0] * np.pi / 180
    phi_r = independent_variables_without_angle[1] * np.pi / 180
    n_0 = independent_variables_without_angle[2]
    polarization = independent_variables_without_angle[2]
    precision = precision * np.pi / 180
    rho_L = np.exp(log_rho_L)
    n = np.exp(log_n_minus_one) + 1
    gamma = np.exp(log_gamma)
    parameters = [rho_L, n, gamma]
    if log_K is not None:
        parameters.append(np.exp(log_K))
    return BRIDF(theta_r, phi_r, theta_i, n_0, polarization, parameters, precision=precision)


# gets weights to use for the gaussian distribution of incident angles
# x is the array of incident angles to get weights for
# outputs weights that should correspond to a list of N points separated by 
def get_relative_gaussian_weights(x, mu, sigma):
    return [np.exp(-(x_ - mu) * (x_ - mu) / (2 * sigma * sigma)) for x_ in x]


def average_BRIDF(independent_variables, log_rho_L, log_n_minus_one, log_gamma, log_K, average_angle, precision=0.25, sigma_theta_i=2.0):
    # For a single point, calculate BRIDF values on a grid within a circle of diameter average_angle from center
    # and average them all; grid spacing is given by precision; all angles in degrees
    # Timing tests suggest that there's almost no speed up between a precision of 0.5 and 1 (for average_angle=4), and
    # only a slight speed up for 0.5 vs 0.25, so a precision of 0.5 or 0.25 seems reasonable
    # If sigma_theta_i>0, also add smearing in theta_i with width sigma_theta_i, to simulate the finite width of the incoming light
	# 0.14" collimator ID, 7" length -> arcsin(0.14/7) ~= 1.1 deg full opening angle for cone (if all light starts at a point before collimator)
    theta_r0 = independent_variables[0]
    phi_r0 = independent_variables[1]
    theta_i = independent_variables[2]
    theta_i_scalar = theta_i
    
    #print(theta_i,theta_r0,phi_r0)
    if sigma_theta_i > 0: # Expand theta_i to multiple values, drawn from a distribution centered on theta_i with width +/- sigma_theta_i
        n_rand = 20 # number of random values to draw
        #theta_i_prime = np.repeat(theta_i, n_rand)
        #delta_theta_i = np.random.normal(scale=sigma_theta_i, size=n_rand)
        eps=0.0001
        delta_theta_i = np.linspace(-3 * sigma_theta_i-eps,sigma_theta_i * 3+eps,n_rand) # offset by eps to avoid errors related to floating point precision
        #print(delta_theta_i)
        theta_i = theta_i+delta_theta_i
    


    # Make a grid of points in theta_r and phi_r
    #print(average_angle)
    theta_r_list = np.linspace(theta_r0-average_angle/2.,theta_r0+average_angle/2.,average_angle/precision+1)
    #print(theta_r_list)
    phi_r_list = np.linspace(phi_r0-average_angle/2.,phi_r0+average_angle/2.,average_angle/precision+1)
    if np.size(theta_i)>1: # If smearing theta_i, make the grid 3D
        tt, pp, ti = np.meshgrid(theta_r_list, phi_r_list, theta_i)
    else:
        tt, pp = np.meshgrid(theta_r_list, phi_r_list)
    tt=tt.flatten()
    pp=pp.flatten()
    # Include only points w/in circle
    in_circle=(tt-theta_r0)**2+(pp-phi_r0)**2 <= (average_angle/2)**2
    tt = tt[in_circle]
    pp = pp[in_circle]
    if np.size(theta_i)>1:
        ti=ti.flatten()
        theta_i=ti[in_circle]
        #print(theta_i[:25])
        weights = get_relative_gaussian_weights(theta_i, theta_i_scalar, sigma_theta_i)
    else:
        weights = [1.0]*len(tt)
    #print(tt[:25],pp[:25])
    return np.sum(weights * unvectorized_fitter([tt,pp,theta_i]+list(independent_variables[3:]),log_rho_L, log_n_minus_one, log_gamma, log_K, precision=precision)) / np.sum(weights)

    
# takes arrays of independent variable lists; last two independent variables are precision and average_angle (assumed to be the same for all points)
def fitter(independent_variables_array, log_rho_L, log_n_minus_one, log_gamma, log_K=None):
    average_angle = independent_variables_array[0][-1]
    precision = independent_variables_array[0][-2]
    # Could be fitting many runs, w/ different params like incident angle or medium index

    if precision < 0: 
        arr = []
        for independent_variables in independent_variables_array:
            arr.append(unvectorized_fitter(independent_variables, log_rho_L, log_n_minus_one, log_gamma, log_K, precision=precision))
        if average_angle > 0: # Just average adjacent points in plane
            theta_r_in_degrees_array = [ind_vars[0] for ind_vars in independent_variables_array]
            # need to have theta_r in degrees! ind. vars is a list, each entry of which contains a list including theta_r (first)
            return average_by_angle(theta_r_in_degrees_array, arr, average_angle)
        else: # No averaging
            return arr
    else: # Do a full 2D average across both theta_r and phi_r, with grid spacing given by precision
        # t1=time.time()
        # print("Average by angle time: {0}".format(t1-t0))
        # t2=time.time()
        avgs = []
        for independent_variables in independent_variables_array: # For each point, calculate the BRIDF along a grid nearby and average it
            avgs.append(average_BRIDF(independent_variables, log_rho_L, log_n_minus_one, log_gamma, log_K, average_angle, precision=precision))
        # t3=time.time()
        # print("Average BRIDF time: {0}".format(t3-t2))
        return avgs

def fitter_with_angle(independent_variables_without_angle_array, theta_i, log_rho_L, log_n_minus_one, log_gamma, log_K=None, precision=-1):
    arr = []
    for independent_variables_without_angle in independent_variables_without_angle_array:
        fit_val = unvectorized_fitter_with_angle(independent_variables_without_angle, theta_i, log_rho_L, log_n_minus_one, log_gamma, log_K, precision=precision)
        arr.append(fit_val)
    return arr

# independent variables array has as elements lists of independent variables for each point
# at each point, it has the form [theta_r_in_degrees, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization]
# average_angle sets the angle in degrees over which to average the BRIDF for fitting to the data
# if precision is positive, the averaging is a full 2D average over theta_r, phi_r with grid size given by precision; else just a 1D average in plane
# use_errs sets whether to use the run's relative_std values to weight the fitting by or to use uniform errors
# use_spike sets whether to include the specular spike in the BRIDF model
def fit_parameters(independent_variables_array_intensity_array, p0=[0.5, 1.5, 0.05], average_angle=0, precision=-1, use_errs=True, use_spike=False):
    independent_variables_array = independent_variables_array_intensity_array[0]
    independent_variables_array = [list+[precision,average_angle] for list in independent_variables_array]
    #print(independent_variables_array)
    intensity_array = independent_variables_array_intensity_array[1]
    if use_errs: std_array = independent_variables_array_intensity_array[2]
    else: std_array = None
    p0_log = [np.log(p0[0]), np.log(p0[1]-1), np.log(p0[2])]
    if use_spike:
        if len(p0)>3: p0_log.append(np.log(p0[3]))
        else: p0_log.append(np.log(2.0)) # default K value
    # initial parameters are the ones found in the paper
    # fitter_p0 = fitter(independent_variables_array, p0_log[0], p0_log[1], p0_log[2])
    # res_fit_p0 = np.sum((np.array(intensity_array) - np.array(fitter_p0))**2)
    fit_params_log = scipy.optimize.curve_fit(fitter, np.array(independent_variables_array), np.array(intensity_array),
                                          p0=p0_log, sigma=std_array)[0]
    # fitter_popt = fitter(independent_variables_array, fit_params[0], fit_params[1], fit_params[2])
    # res_fit_popt = np.sum((np.array(intensity_array) - np.array(fitter_popt))**2)
    # print(res_fit_p0, res_fit_popt)
    params = [np.exp(fit_params_log[0]), np.exp(fit_params_log)[1] + 1, np.exp(fit_params_log[2])]
    if len(fit_params_log)>3: params.append(np.exp(fit_params_log[3]))
    return params


def change_theta_i(vars_intensity_array, new_theta_i):
    new_points = copy.deepcopy(vars_intensity_array)
    for point in new_points[0]:
        point[2] = new_theta_i
    return new_points


# returns [fitted_theta_i, rho_L, n, gamma]
# not intended for running data with more than one incident angle at a time
def fit_parameters_and_angle(independent_variables_array_intensity_array, p0=[0.5, 1.5, 0.05], average_angle=0, precision=-1, use_errs=True, use_spike=False):
    independent_variables_array = independent_variables_array_intensity_array[0]
    independent_variables_array = [list+[precision, average_angle] for list in independent_variables_array]
    intensities = independent_variables_array_intensity_array[1]
    
    # determine if there is a sharp peak
    max_ind = intensities.index(max(intensities))
    theta_r_peak = independent_variables_array[max_ind][0]
    theta_i_peak = independent_variables_array[max_ind][2]
    # peak is within 15 degrees of where it's expected
    if np.abs(theta_r_peak - theta_i_peak) < 15:
        # use peak as theta_i
        points_1 = change_theta_i(independent_variables_array_intensity_array, theta_r_peak)
        parameters_with_theta_i_peak = fit_parameters(points_1,p0,average_angle, use_errs=use_errs, use_spike=use_spike)
        theta_r_in_degrees_array = [point[0] for point in independent_variables_array]
        point0 = points_1[0][0]
        fit_array = BRIDF_plotter(theta_r_in_degrees_array, *point0[1:], parameters_with_theta_i_peak)
        fit_peak_index = fit_array.index(max(fit_array))
        fit_peak = independent_variables_array[fit_peak_index][0]
        peak_offset = fit_peak - theta_r_peak
        points_2 = change_theta_i(independent_variables_array_intensity_array, theta_r_peak - peak_offset)
        return [theta_r_peak - peak_offset] + fit_parameters(points_2,p0,average_angle, use_errs=use_errs, use_spike=use_spike)

    else:
        # use experimental as theta_i
        return [theta_i_peak] + fit_parameters(independent_variables_array_intensity_array,p0,average_angle, use_spike=use_spike)


# returns [fitted_theta_i, rho_L, n, gamma]
# not intended for running data with more than one incident angle at a time
# Note: not yet updated for newer data format (Run)
def fit_parameters_and_angle_with_starting_theta_i(points, starting_theta_i, average_angle=0):
    independent_variables_without_angle_array = []
    intensity_array = []
    for point in points:
        independent_variables_without_angle_array.append([point.theta_r_in_degrees, point.phi_r_in_degrees,
                                            point.n_0, point.polarization, average_angle])
        intensity_array.append(point.intensity)
    # initial parameters are the ones found in the paper
    fit_params = scipy.optimize.curve_fit(fitter_with_angle, independent_variables_without_angle_array, intensity_array,
                                          p0=[starting_theta_i, np.log(0.5), np.log(1.5 - 1), np.log(0.05)])[0]
    return [fit_params[0], np.exp(fit_params[1]), np.exp(fit_params)[2] + 1, np.exp(fit_params[3])]

# Calculates hemispherical reflectance for a given incident angle
# Note that specular lobe depends on whether specular spike is included or not
def reflectance(theta_i_in_degrees, n_0, polarization, parameters):

    theta_i = theta_i_in_degrees * np.pi / 180

    a = lambda x: 0
    b = lambda x: np.pi / 2

    def BRIDF_int(theta_r, phi_r):
        # solid angle not used here
        G = G_calc(theta_r, phi_r, theta_i, n_0, polarization, parameters)
        # precision=-1 ensures that specular spike isn't included, though it may modify specular lobe normalization if parameters include K
        # will add specular spike by hand, since integration won't have a fixed grid, so ensuring a real delta function is difficult
        # sin(theta_r) is from solid angle, G is to correct for the fact that shadowed/masked light ends up reflected, see eq. 16 of Silva 2009 Xe scintillation
        return BRIDF(theta_r, phi_r, theta_i, n_0, polarization, parameters, precision=-1) * np.sin(theta_r) / G 

    spec_spike = 0
    if len(parameters)>3:
        n=parameters[1]
        K=parameters[3]
        C = specular_spike_norm(theta_i, theta_i, K)
        F_ = F(theta_i, n_0, n, polarization)
        spec_spike = C * F_
        
    return scipy.integrate.dblquad(BRIDF_int, -np.pi, np.pi, a, b)[0] + spec_spike

# TODO: check if 1/G is supposed to go here
def reflectance_diffuse(theta_i_in_degrees, n_0, polarization, parameters):

    theta_i = theta_i_in_degrees * np.pi / 180

    a = lambda x: 0
    b = lambda x: np.pi / 2

    def BRIDF_int(theta_r, phi_r):
        # solid angle not used here
        G = G_calc(theta_r, phi_r, theta_i, n_0, polarization, parameters)
        # sin(theta_r) is from solid angle, G is to correct for the fact that shadowed/masked light ends up reflected, see eq. 16 of Silva 2009 Xe scintillation
        return BRIDF_diffuse(theta_r, phi_r, theta_i, n_0, polarization, parameters) * np.sin(theta_r) / G

    return scipy.integrate.dblquad(BRIDF_int, -np.pi, np.pi, a, b)[0]

# Note that specular lobe depends on whether specular spike is included or not
def reflectance_specular(theta_i_in_degrees, n_0, polarization, parameters):

    theta_i = theta_i_in_degrees * np.pi / 180

    a = lambda x: 0
    b = lambda x: np.pi / 2

    def BRIDF_int(theta_r, phi_r):
        G = G_calc(theta_r, phi_r, theta_i, n_0, polarization, parameters)
        # precision=-1 ensures that specular spike isn't included, though it may modify specular lobe normalization if parameters include K
        # will add specular spike by hand, since integration won't have a fixed grid, so ensuring a real delta function is difficult
        # sin(theta_r) is from solid angle, G is to correct for the fact that shadowed/masked light ends up reflected, see eq. 16 of Silva 2009 Xe scintillation
        return BRIDF_specular(theta_r, phi_r, theta_i, n_0, polarization, parameters, precision=-1) * np.sin(theta_r) / G 

    spec_spike = 0
    if len(parameters)>3:
        n=parameters[1]
        K=parameters[3]
        C = specular_spike_norm(theta_i, theta_i, K)
        F_ = F(theta_i, n_0, n, polarization)
        spec_spike = C * F_
        
    return scipy.integrate.dblquad(BRIDF_int, -np.pi, np.pi, a, b)[0] + spec_spike

# uses 3d grid of parameters
# independent variables array has as elements lists of independent variables for each point
# at each point, it has the form [theta_r_in_degrees, phi_r_in_degrees, theta_i_in_degrees, n_0, polarization]
def fit_parameters_new(independent_variables_array, intensity_array, std_array, rho_start, rho_end, rho_num, n_start, n_end, n_num, gamma_start, gamma_end, gamma_num, plot=True, show=True):
    rho_array = np.linspace(rho_start, rho_end, rho_num)
    n_array = np.linspace(n_start, n_end, n_num)
    gamma_array = np.linspace(gamma_start, gamma_end, gamma_num)

    grid_points = []
    for rho in rho_array:
        for n in n_array:
            for gamma in gamma_array:
                grid_points.append([rho, n, gamma])

    length = len(grid_points)

    chi2 = []
    for i in range(len(grid_points)):
        grid_point = grid_points[i]
        chi2.append(chi_squared(independent_variables_array, intensity_array, std_array, grid_point))
        if i%10 == 0:
            print(str(i) + "/" + str(length))

    min_index = np.argmin(chi2)
    min_params = grid_points[min_index]

    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        p0 = []
        p1 = []
        p2 = []
        color = []
        for i in range(len(grid_points)):
            p0.append(grid_points[i][0])
            p1.append(grid_points[i][1])
            p2.append(grid_points[i][2])
            color.append(chi2[i])

        a = ax.scatter(p0, p1, p2, c=color, norm=matplotlib.colors.LogNorm(), s=15)

        fig.colorbar(a)

        ax.text2D(0., 0.95, "optimal parameters: rho_L = " + str(round(min_params[0], 5)) + ", n = " + str(round(min_params[1], 5)) + ", gamma = " + str(round(min_params[2], 5)), transform=ax.transAxes)

        ax.set_xlabel('rho_L')
        ax.set_ylabel('n')
        ax.set_zlabel('gamma')

        if show:
            plt.show()

    return grid_points[min_index]

def chi_squared(independent_variables_array, intensity_array, std_array, grid_point):
    # assumes entire array is at one incident angle and external index of refraction
    theta_i_in_degrees = independent_variables_array[0][2] 
    n_0 = independent_variables_array[0][3]

    # assumes data is 85 degrees to 0 degrees in 1 degree increments
    theta_r_in_degrees_array = np.linspace(85, 0, 86)

    fit = BRIDF_plotter(theta_r_in_degrees_array, 0, theta_i_in_degrees, n_0, 0.5, grid_point, average_angle=0, precision=-1, sigma_theta_i=0)

    chi2 = 0
    for i in range(len(fit)):
        chi2 += (fit[i] - intensity_array[i]) * (fit[i] - intensity_array[i]) / (std_array[i] * std_array[i])

    return chi2






