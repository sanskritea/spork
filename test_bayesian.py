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


def BayesianT1(N, T1):

    # decay rates grid
    gamma_lower = 0.1  # in us^-1
    gamma_upper = 2  # in us^-1
    n_gamma = 2000
    gamma_plus_arr = np.linspace(gamma_lower, gamma_upper, n_gamma)
    gamma_minus_arr = np.linspace(gamma_lower, gamma_upper, n_gamma)
    delta_gamma = gamma_plus_arr[1] - gamma_plus_arr[0]

    # decay rates distribution
    gamma_distr = np.ones((n_gamma, n_gamma)) / (n_gamma ** 2)

    # relevant experimental times
    n_tau = 1000
    tau_plus_arr = np.linspace(1, 5500, n_tau)  # in us
    tau_minus_arr = np.linspace(1, 5500, n_tau)  # in us

    # begin with a flat prior
    prior_gamma_distr = gamma_distr  # starting with flat prior in gamma_plus
    repetitions = 10

    for rep in range(N):

        print("Doing rep number ", rep)

        # calculate gamma_plus_mean and gamma_minus_mean from prior_gamma_distr
        gamma_mean, gamma_sigma = calculate_gamma_params(
            N,
            rep,
            prior_gamma_distr,
            gamma_plus_arr,
            gamma_minus_arr,
            n_gamma,
            delta_gamma,
        )

        # use mean gamma values to find optimized taus from cost function
        tau_plus_opt, tau_minus_opt = calculate_tau_opt(
            tau_plus_arr, tau_minus_arr, gamma_mean, gamma_sigma, n_tau, repetitions)
        tau_vals = [tau_plus_opt, tau_minus_opt]
        # print("tau_vals ", tau_vals)

        # use taus in measurement
        M_measured = np.zeros((2, repetitions))

        for t in range(2):

            # FAKE COUNTING
            M_measured[t] = fake_counts(tau_vals[t], repetitions, T1[t])

        # calculate likelihood from measurement result
        likelihood = calculate_likelihood(
            M_measured, tau_vals, gamma_plus_arr, gamma_minus_arr, n_gamma
        )

        # calculate posterior
        posterior_unnorm = np.zeros((n_gamma, n_gamma))
        for n in range(n_gamma):
            for nn in range(n_gamma):
                posterior_unnorm[n][nn] = likelihood[n][nn] * prior_gamma_distr[n][nn]

        # print("np.sum(prior_gamma_distr)", np.sum(prior_gamma_distr))
        # print("np.sum(posterior_unnorm)", np.sum(posterior_unnorm))

        # normalize posterior and update prior
        prior_gamma_distr = posterior_unnorm / (
            (delta_gamma**2) * np.sum(posterior_unnorm)
        )


def calculate_gamma_params(
    N, rep, prior_gamma_distr, gamma_plus_arr, gamma_minus_arr, n_gamma, delta_gamma
):

    gamma_plus_distr = []
    gamma_minus_distr = []
    gamma_mean = []
    gamma_sigma = []
    sigma_sums = []

    for n in range(n_gamma):
        sum_all_rows = 0
        sum_all_cols = 0
        for nn in range(n_gamma):
            sum_all_rows += prior_gamma_distr[n][nn]
            sum_all_cols += prior_gamma_distr[nn][n]
        gamma_plus_distr.append(sum_all_rows)
        gamma_minus_distr.append(sum_all_cols)

    if rep == 0 or rep == (N-1):
        fig, axs = plt.subplots(2)
        fig.suptitle('Gamma distributions for '+ str(rep) + 'th iteration')
        axs[0].plot(gamma_plus_arr, gamma_plus_distr)
        axs[1].plot(gamma_minus_arr, gamma_minus_distr)
        plt.show()

    gamma_mean_plus = 0
    gamma_mean_minus = 0
    for n in range(n_gamma):
        gamma_mean_plus += gamma_plus_arr[n] * gamma_plus_distr[n] 
        gamma_mean_minus += gamma_minus_arr[n] * gamma_minus_distr[n] 

    sigma_sums_plus = 0
    sigma_sums_minus = 0
    for n in range(n_gamma):
        sigma_sums_plus += (gamma_plus_arr[n] - gamma_mean_plus) ** 2
        sigma_sums_minus += (gamma_plus_arr[n] - gamma_mean_minus) ** 2

    gamma_mean = [gamma_mean_plus, gamma_mean_minus]
    gamma_sigma = [
        (sigma_sums_plus / n_gamma) ** 0.5,
        (sigma_sums_minus / n_gamma) ** 0.5,
    ]

    # print('gamma_mean ', gamma_mean)
    # print('gamma_sigma ', gamma_sigma)

    return gamma_mean, gamma_sigma


def calculate_tau_opt(tau_plus_arr, tau_minus_arr, gamma_mean, gamma_sigma, n_tau, repetitions):

    cost_function = np.zeros((n_tau, n_tau))

    for i in range(n_tau):
        for j in range(n_tau):
            cost_function[i][j] = (
                (
                    (gamma_sigma[0] / gamma_mean[0]) ** 2
                    + (gamma_sigma[1] / gamma_mean[1]) ** 2
                )
                ** 0.5
            ) * 2 * repetitions * (tau_plus_arr[i] + tau_minus_arr[j])

    # print('cost_function ', cost_function)

    tau_plus_optimized = tau_plus_arr[np.argmin(cost_function, 0)[0]]
    tau_minus_optimized = tau_minus_arr[np.argmin(cost_function, 1)[0]]

    return tau_plus_optimized, tau_minus_optimized


def fake_counts(tau, num_samples, T1):

    return np.random.normal((np.exp(-(tau / T1))), 1e-3, num_samples)


def calculate_likelihood(
    M_measured, tau_vals, gamma_plus_arr, gamma_minus_arr, n_gamma
):

    M_plus_measured = np.mean(M_measured[0]) * np.ones((n_gamma, n_gamma))
    sigma_M_plus_measured = np.std(M_measured[0])
    # print('sigma_M_plus_measured ', sigma_M_plus_measured)

    M_minus_measured = np.mean(M_measured[1]) * np.ones((n_gamma, n_gamma))
    sigma_M_minus_measured = np.std(M_measured[1])

    tau_plus = tau_vals[0]
    tau_minus = tau_vals[1]

    M_plus_tilde = np.zeros((n_gamma, n_gamma))
    M_minus_tilde = np.zeros((n_gamma, n_gamma))

    for n in range(n_gamma):

        gamma_plus = gamma_plus_arr[n]

        for nn in range(n_gamma):

            gamma_minus = gamma_minus_arr[nn]

            g = ((gamma_plus**2) + (gamma_minus**2) - (gamma_plus * gamma_minus)) ** 0.5

            beta_plus = gamma_plus + gamma_minus + g
            beta_minus = gamma_plus + gamma_minus - g

            M_plus_tilde[n][nn] = (
                ((g + gamma_plus) * np.exp(-beta_plus * tau_plus))
                + ((g - gamma_plus) * np.exp(-beta_minus * tau_plus))
            ) / (2 * g)

            M_minus_tilde[n][nn] = (
                ((g + gamma_minus) * np.exp(-beta_plus * tau_minus))
                + ((g - gamma_minus) * np.exp(-beta_minus * tau_minus))
            ) / (2 * g)

    # print('M_plus_measured ', M_plus_measured)
    # print('M_plus_tilde ', M_plus_tilde)

    chi_plus = (M_plus_measured - M_plus_tilde) / ((2**0.5) * (sigma_M_plus_measured))
    # print('chi_plus ', chi_plus)
    chi_minus = (M_minus_measured - M_minus_tilde) / ((2**0.5) * (sigma_M_minus_measured))
    # print('chi_minus ', chi_minus)

    chi_sq = (chi_plus**2) + (chi_minus**2)
    chi_sq_final = chi_sq - np.min(chi_sq) * np.ones((n_gamma, n_gamma))

    likelihood = np.exp(-(chi_sq_final))
    # print('likelihood ', likelihood)

    return likelihood
