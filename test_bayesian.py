"""
TESTING BAYESIAN OPERATION WITH FAKE MEASUREMENT DATA
"""

"""
Questions:
1. Should gammas be log spaced? Should taus? 
2. Figure out the cost function 
3. How to get fake tau-dependent counts data?
4. Is the covariance/std calculation valid? Paper requires sigma_gamma_plus but how to get that from the joint pdf prior? 
"""


import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import derivative

def BayesianT1(N, T1):

    # decay rates grid
    gamma_lower = 0.1   # in ms^-1
    gamma_upper = 10    # in ms^-1
    n_gamma = 500
    gamma_plus_arr = np.linspace(gamma_lower, gamma_upper, n_gamma) / 1000     # for us^-1
    # print(gamma_plus_arr)
    gamma_minus_arr = np.linspace(gamma_lower, gamma_upper, n_gamma) / 1000    # for us^-1
    gamma_grid = np.meshgrid(gamma_plus_arr, gamma_minus_arr, indexing = 'ij')
    # print('gamma_grid', gamma_grid)
    delta_gamma = gamma_plus_arr[1] - gamma_plus_arr[0]

    # decay rates distribution
    gamma_distr = np.ones((n_gamma, n_gamma)) / (n_gamma**2)

    # relaxometry delay tau grid
    tau_lower = 3       # in us
    tau_upper = 5500    # in us
    n_tau = 1000
    tau_plus_arr = np.linspace(tau_lower, tau_upper, n_tau)
    tau_minus_arr = np.linspace(tau_lower, tau_upper, n_tau)
    tau_grid = np.meshgrid(tau_plus_arr, tau_minus_arr, indexing = 'ij')

    # begin with a flat prior in gammas
    prior_gamma = gamma_distr
    repetitions = 10

    for rep in range(N):

        print("Doing rep number ", rep)

        # calculate gamma_plus_mean and gamma_minus_mean from prior_gamma_distr
        gamma_mean, gamma_sigma = calculate_gamma_params(
            N, rep, prior_gamma, gamma_grid, gamma_plus_arr, gamma_minus_arr, n_gamma, delta_gamma
        )

        # use mean gamma values to find optimized taus from cost function
        tau_opt = calculate_tau_opt(tau_grid, gamma_mean, gamma_sigma, repetitions)

        # use taus in measurement
        M_measured = np.zeros((2, repetitions))
        for t in range(2):
            # FAKE COUNTING
            M_measured[t] = fake_counts(tau_opt[t], repetitions, T1[t])

        # calculate likelihood from measurement result
        likelihood = calculate_likelihood(M_measured, tau_opt, gamma_grid, n_gamma)

        # print('prior_gamma ', prior_gamma)
        # print('likelihood ', likelihood)
        # calculate posterior
        posterior_gamma_unnorm = np.multiply(likelihood, prior_gamma)

        # normalize posterior and update prior
        prior_gamma = posterior_gamma_unnorm / np.sum(posterior_gamma_unnorm)
        # print('new prior_gamma ', prior_gamma)


def calculate_gamma_params(
    N, rep, prior_gamma, gamma_grid, gamma_plus_arr, gamma_minus_arr, n_gamma, delta_gamma
):

    gamma_plus, gamma_minus = gamma_grid

    gamma_plus_distr = np.sum(prior_gamma, 0) 
    # print('gamma_plus_distr ', gamma_plus_distr)
    gamma_minus_distr = np.sum(prior_gamma, 1) 
    # print('gamma_minus_distr ', gamma_minus_distr)

    gamma_mean_plus = np.sum(gamma_plus * gamma_plus_distr) / n_gamma
    gamma_mean_minus = np.sum(gamma_minus * gamma_minus_distr) / n_gamma

    gamma_sigma_plus = (np.sum((gamma_plus - gamma_mean_plus) ** 2) / n_gamma) ** 0.5
    gamma_sigma_minus = (np.sum((gamma_minus - gamma_mean_minus) ** 2) / n_gamma) ** 0.5

    gamma_mean = [gamma_mean_plus, gamma_mean_minus]
    gamma_sigma = [gamma_sigma_plus,gamma_sigma_minus]
    print('gamma_mean ', gamma_mean)
    # print('gamma_sigma ', gamma_sigma)

    # PRINT GAMMA_PLUS AND GAMMA_MINUS FROM CURRENT PDF
    # print("T1_plus estimate in us: ", 1 / (gamma_plus_arr[np.argmax(gamma_plus_distr)]))
    # print("T1_minus estimate in us: ", 1 / (gamma_minus_arr[np.argmax(gamma_minus_distr)]))

    # PLOT GAMMA_PLU AND GAMMA_MINUS PDFs
    # if rep == 0 or rep == (N - 1):
    #     fig, axs = plt.subplots(2)
    #     fig.suptitle("Gamma (ms^-1) distributions for " + str(rep) + "th iteration")
    #     axs[0].plot(1000 * gamma_plus_arr, gamma_plus_distr)
    #     axs[1].plot(1000 * gamma_minus_arr, gamma_minus_distr)
    #     plt.show()

    return gamma_mean, gamma_sigma


def calculate_tau_opt(tau_grid, gamma_mean, gamma_sigma, repetitions):

    tau_plus, tau_minus = tau_grid
    tp, tm = tau_plus.flatten(), tau_minus.flatten()

    cost_function = (
        (gamma_sigma[0] / gamma_mean[0]) ** 2
        + (gamma_sigma[1] / gamma_mean[1] ** 2) ** 0.5
    ) * (2 * repetitions * (tau_plus + tau_minus)) ** 0.5
    # print('cost_function ', cost_function)

    tau_optimized = [tp[np.argmin(cost_function)], tm[np.argmin(cost_function)]]
    # print('tau_optimized ', tau_optimized)

    return tau_optimized


def nob_g(gp, gm):
    
    return (gp ** 2 + gm ** 2 - gp * gm) ** 0.5

    
def nob_bp(gp, gm):

    return (gp + gm + nob_g(gp, gm))


def nob_bm(gp, gm):

    return (gp + gm - nob_g(gp, gm))


def nob_mtp(gp, gm):

    g = nob_g(gp, gm)
    bp = nob_bp(gp, gm)
    bm = nob_bm(gp, gm)
    
    return ((g + gp) * np.exp(-1 * bp * tp) + (g - gp) * np.exp(-1 * bm * tp)) / (2 * g)


def nob_mtp(gp, gm):

    g = nob_g(gp, gm)
    bp = nob_bp(gp, gm)
    bm = nob_bm(gp, gm)
    
    return ((g + gp) * np.exp(-1 * bp * tp) + (g - gp) * np.exp(-1 * bm * tp)) / (2 * g)   


# def nob_cost_function(gp, gm):





def fake_counts(tau, num_samples, T1):

    return np.random.normal((np.exp(-(tau / T1))), 1e-4, num_samples)


def calculate_likelihood(M_measured, tau_opt, gamma_grid, n_gamma):

    M_plus_measured = np.mean(M_measured[0]) * np.ones((n_gamma, n_gamma))
    sigma_M_plus_measured = np.std(M_measured[0])
    # print('sigma_M_plus_measured ', sigma_M_plus_measured)

    M_minus_measured = np.mean(M_measured[1]) * np.ones((n_gamma, n_gamma))
    sigma_M_minus_measured = np.std(M_measured[1])

    gamma_plus, gamma_minus = gamma_grid
    tau_plus, tau_minus = tau_opt[0], tau_opt[1]

    g = ((gamma_plus**2) + (gamma_minus**2) - (gamma_plus * gamma_minus)) ** 0.5
    beta_plus = gamma_plus + gamma_minus + g
    beta_minus = gamma_plus + gamma_minus - g

    M_plus_tilde = (
        ((g + gamma_plus) * np.exp(-beta_plus * tau_plus))
        + ((g - gamma_plus) * np.exp(-beta_minus * tau_plus))
    ) / (2 * g)
    M_minus_tilde = (
        ((g + gamma_minus) * np.exp(-beta_plus * tau_minus))
        + ((g - gamma_minus) * np.exp(-beta_minus * tau_minus))
    ) / (2 * g)
    # print('M_plus_measured ', M_plus_measured)
    # print('M_plus_tilde ', M_plus_tilde)

    chi_plus = (M_plus_measured - M_plus_tilde) / ((2**0.5) * (sigma_M_plus_measured))
    chi_minus = (M_minus_measured - M_minus_tilde) / (
        (2**0.5) * (sigma_M_minus_measured)
    )

    chi_sq = (chi_plus**2) + (chi_minus**2)
    chi_sq_final = chi_sq - np.min(chi_sq) 
    likelihood = np.exp(-(chi_sq_final))

    return likelihood
