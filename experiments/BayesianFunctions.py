'''

2024-08-29
Sanskriti Chitransh

File containing all functions needed for Bayesian adaptive T1 

'''

import time
import numpy as np
from scipy.interpolate import interp1d
import sympy
from sympy import Symbol
from sympy import lambdify

class Symbols:

    tau_plus = Symbol("tau_-")
    tau_minus = Symbol("tau_+")
    gamma_plus = Symbol("Gamma_+")
    gamma_minus = Symbol("Gamma_-")

    G = sympy.sqrt(gamma_plus**2 + gamma_minus**2 - gamma_plus * gamma_minus)
    beta_plus = gamma_plus + gamma_minus + G
    beta_minus = gamma_plus + gamma_minus - G

    M_tilde_plus = (G + gamma_plus) * sympy.exp(-tau_plus * beta_plus) + (
        G - gamma_plus
    ) * sympy.exp(-tau_plus * beta_minus)
    M_tilde_minus = (G + gamma_minus) * sympy.exp(-tau_minus * beta_plus) + (
        G - gamma_minus
    ) * sympy.exp(-tau_minus * beta_minus)
    M_tilde_plus = M_tilde_plus / (2 * G)
    M_tilde_minus = M_tilde_minus / (2 * G)

    pp = lambdify(
        [gamma_plus, gamma_minus, tau_plus],
        sympy.diff(M_tilde_plus, gamma_plus).simplify(),
    )
    pm = lambdify(
        [gamma_plus, gamma_minus, tau_plus],
        sympy.diff(M_tilde_plus, gamma_minus).simplify(),
    )
    mm = lambdify(
        [gamma_plus, gamma_minus, tau_minus],
        sympy.diff(M_tilde_minus, gamma_minus).simplify(),
    )
    mp = lambdify(
        [gamma_plus, gamma_minus, tau_minus],
        sympy.diff(M_tilde_minus, gamma_plus).simplify(),
    )
    Mtplus = lambdify([gamma_plus, gamma_minus, tau_plus], M_tilde_plus.simplify())
    Mtminus = lambdify([gamma_plus, gamma_minus, tau_minus], M_tilde_minus.simplify())


class BF:

	def trap_normalize_2D_pdf(pdf, gamma_grid):

		gamma_plus, gamma_minus = gamma_grid[0][0], np.transpose(gamma_grid[1])[0]
		norm = np.trapz(np.trapz(pdf, gamma_minus), gamma_plus)

		return pdf / norm


	def calculate_conf_intervals(prior_gamma, gamma_grid, delta_gamma):

	    gamma_plus = gamma_grid[0][0]
	    gamma_minus = np.transpose(gamma_grid[1])[0]

	    # print('Normalizing pdfs')
	    gamma_plus_distr = np.sum(prior_gamma, 0) * delta_gamma
	    norm_gamma_plus_distr = gamma_plus_distr / (np.sum(gamma_plus_distr) * delta_gamma)
	    # print('normalizing_factor_gamma_plus_distr ', (np.sum(gamma_plus_distr) * delta_gamma))
	    gamma_minus_distr = np.sum(prior_gamma, 1) * delta_gamma
	    norm_gamma_minus_distr = gamma_minus_distr / (np.sum(gamma_minus_distr) * delta_gamma)
	    # print('normalizing_factor_gamma_minus_distr ', (np.sum(gamma_minus_distr) * delta_gamma))

	    # print('Calculating cdfs')
	    gamma_plus_cdf = np.cumsum(norm_gamma_plus_distr) * delta_gamma
	    gamma_minus_cdf = np.cumsum(norm_gamma_minus_distr) * delta_gamma

	    # print('Interpolating')
	    gamma_plus_cdf_interp = interp1d(gamma_plus_cdf, gamma_plus)
	    gamma_minus_cdf_interp = interp1d(gamma_minus_cdf, gamma_minus)

	    gamma_plus_vals = gamma_plus_cdf_interp(np.array([0.05, 0.5, 0.95]))
	    gamma_minus_vals = gamma_minus_cdf_interp(np.array([0.05, 0.5, 0.95]))

	    return gamma_plus_vals, gamma_minus_vals


	def nob_calculate_tau_opt(tau_grid, repetitions, gamma_plus, gamma_minus, T_overhead):
	    """
	    Cost function is only a function of tau, and we want to optimize it
	    to find the 'best' value of tau. This is done by finding the minimum
	    of the cost function evaluated over a grid of tau values.

	    On the other hand, the cost function is only a function of one value
	    of gamma_plus and gamma_minus. These values are the average value of
	    gamma under the gamma prior. gamma_plus and gamma_minus are single
	    float values, not arrays.
	    """
	    start_tau_optimization = time.time()
	    tau_plus, tau_minus = tau_grid

	    num = (
	        (gamma_minus * Symbols.mm(gamma_plus, gamma_minus, tau_minus)) ** 2
	        + (gamma_minus * Symbols.pm(gamma_plus, gamma_minus, tau_plus)) ** 2
	        + (gamma_plus * Symbols.mp(gamma_plus, gamma_minus, tau_minus)) ** 2
	        + (gamma_plus * Symbols.pp(gamma_plus, gamma_minus, tau_plus)) ** 2
	    )

	    den = (
	         Symbols.pm(gamma_plus, gamma_minus, tau_plus) *  Symbols.mp(gamma_plus, gamma_minus, tau_minus)
	    ) -  Symbols.pp(gamma_plus, gamma_minus, tau_plus) *  Symbols.mm(gamma_plus, gamma_minus, tau_minus)
	    den = den**2

	    T = 2 * repetitions * (tau_plus + tau_minus) + T_overhead

	    # actual cost function to compare to paper
	    cost_function = ((T * num / den) ** 0.5) / (gamma_minus * gamma_plus)

	    # clean up cost function by changing all nans to an arbitrarily large value
	    flag = 0
	    for nni in range(len(tau_plus)):
	        for nnj in range(len(tau_plus)):
	            if np.isnan(cost_function[nni][nnj]) or np.isinf(cost_function[nni][nnj]):
	                # print('AHA!')
	                cost_function[nni][nnj] = 1e100
	                flag += 1

	    tp, tm = tau_plus.flatten(), tau_minus.flatten()
	    min_cost_idx = np.argmin(cost_function.flatten())

	    print('Time taken to optimize tau ', time.time() - start_tau_optimization)

	    return np.array([tp[min_cost_idx], tm[min_cost_idx]])


	def calculate_likelihood(M_measured_mean, M_measured_sigma, tau_opt, gamma_grid):
		
	    start_likelihood_calculation = time.time()
	    M_plus_measured, M_minus_measured = M_measured_mean
	    sigma_M_plus_measured, sigma_M_minus_measured = M_measured_sigma
	    

	    M_plus_tilde, M_minus_tilde = BF.calculate_M_tildes(gamma_grid, tau_opt)

	    chi_plus = (M_plus_measured - M_plus_tilde) / ((2**0.5) * (sigma_M_plus_measured))
	    chi_minus = (M_minus_measured - M_minus_tilde) / (
	        (2**0.5) * (sigma_M_minus_measured)
	    )

	    chi_sq = (chi_plus**2) + (chi_minus**2)
	    chi_sq_final = chi_sq - np.min(chi_sq)
	    likelihood = np.exp(-(chi_sq_final))

	    print('Time taken to calculate likelihood ', time.time() - start_likelihood_calculation)

	    return likelihood


	def calculate_M_tildes(gamma_grid, tau_opt):

	    gamma_plus, gamma_minus = gamma_grid
	    tau_plus, tau_minus = tau_opt
	    # pp, pm, mm, mp, Mtplus, Mtminus = BF.give_sympy_functions()
	    
	    M_plus_tilde = Symbols.Mtplus(gamma_plus, gamma_minus, tau_plus)
	    M_minus_tilde = Symbols.Mtminus(gamma_plus, gamma_minus, tau_minus)

	    return [M_plus_tilde, M_minus_tilde]




