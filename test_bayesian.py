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
from findiff import FinDiff


def BayesianT1(N, T1):

    # decay rates grid
    gamma_lower = 0.0002    # in us^-1
    gamma_upper = 0.5       # in us^-1
    n_gamma = 500
    gamma_plus_arr = np.linspace(gamma_lower, gamma_upper, n_gamma) * 1000
    # print(gamma_plus_arr)
    gamma_minus_arr = np.linspace(gamma_lower, gamma_upper, n_gamma)    
    gamma_grid = np.meshgrid(gamma_plus_arr, gamma_minus_arr, indexing="ij")
    # print('gamma_grid', gamma_grid)
    delta_gamma = gamma_plus_arr[1] - gamma_plus_arr[0]

    # decay rates distribution
    gamma_distr = np.ones((n_gamma, n_gamma)) / (n_gamma**2)

    # relaxometry delay tau grid
    tau_lower = 2           # in us
    tau_upper = 5000        # in us
    n_tau = 500
    tau_plus_arr = np.linspace(tau_lower, tau_upper, n_tau)
    tau_minus_arr = np.linspace(tau_lower, tau_upper, n_tau)
    tau_grid = np.meshgrid(tau_plus_arr, tau_minus_arr, indexing="ij")

    # begin with a flat prior in gammas
    prior_gamma = gamma_distr
    repetitions = N
    printing_and_plotting(gamma_grid, prior_gamma, gamma_plus_arr, gamma_minus_arr)

    for rep in range(repetitions):

        print("Doing measurement number ", (rep + 1))

        # find optimized taus from NOB cost function
        tau_opt = nob_calculate_tau_opt(tau_grid, repetitions, gamma_grid, delta_gamma)

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

    printing_and_plotting(gamma_grid, prior_gamma, gamma_plus_arr, gamma_minus_arr)


def nob_calculate_tau_opt(tau_grid, repetitions, gamma_grid, delta_gamma):

    tau_plus, tau_minus = tau_grid
    tp, tm = tau_plus.flatten(), tau_minus.flatten()
    gamma_plus, gamma_minus = gamma_grid
    M_tilde_plus, M_tilde_minus = calculate_M_tildes(gamma_grid, tau_grid)

    pp = (
        np.exp(-2 * (gamma_minus + gamma_plus) * tau_plus)
        * (
            (
                -1
                * np.exp(
                    (
                        gamma_minus
                        + gamma_plus
                        - (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2)
                        ** 0.5
                    )
                    * tau_plus
                )
                - 1
                * np.exp(
                    (
                        gamma_minus
                        + gamma_plus
                        + (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2)
                        ** 0.5
                    )
                    * tau_plus
                )
            )
            * (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 2.5
            * tau_plus
            + np.exp(
                (
                    gamma_minus
                    + gamma_plus
                    + (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 0.5
                )
                * tau_plus
            )
            * (
                (-0.25 * gamma_minus + 0.5 * gamma_plus)
                * gamma_plus
                * (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 1
                + (0.25 * gamma_minus - 0.5 * gamma_plus)
                * gamma_plus
                * (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 1.5
                * tau_plus
                + (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 2
                * (-0.5 - 0.5 * gamma_minus * tau_plus + 1.5 * gamma_plus * tau_plus)
            )
            + np.exp(
                (
                    gamma_minus
                    + gamma_plus
                    - (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 0.5
                )
                * tau_plus
            )
            * (
                (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 2
                * (0.5 + 0.5 * gamma_minus * tau_plus - 1.5 * gamma_plus * tau_plus)
                + (0.25 * gamma_minus - 0.5 * gamma_plus)
                * gamma_plus
                * (
                    (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 1
                    + (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 1.5
                    * tau_plus
                )
            )
        )
    ) / (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 2.5

    pm = (
        (1 / ((gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 2.5))
        * np.exp(-2 * (gamma_minus + gamma_plus) * tau_plus)
        * (
            (
                -1
                * np.exp(
                    (
                        gamma_minus
                        + gamma_plus
                        - (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2)
                        ** 0.5
                    )
                    * tau_plus
                )
                + 1
                * np.exp(
                    (
                        gamma_minus
                        + gamma_plus
                        + (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2)
                        ** 0.5
                    )
                    * tau_plus
                )
            )
            * gamma_minus
            * (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 2
            * tau_plus
            + (
                -1
                * np.exp(
                    (
                        gamma_minus
                        + gamma_plus
                        - (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2)
                        ** 0.5
                    )
                    * tau_plus
                )
                - 1
                * np.exp(
                    (
                        gamma_minus
                        + gamma_plus
                        + (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2)
                        ** 0.5
                    )
                    * tau_plus
                )
            )
            * (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 2.5
            * tau_plus
            + (1 * gamma_minus - 0.5 * gamma_plus)
            * gamma_plus
            * (
                np.exp(
                    (
                        gamma_minus
                        + gamma_plus
                        - (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2)
                        ** 0.5
                    )
                    * tau_plus
                )
                * (
                    -0.5
                    * (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 1
                    - 0.5
                    * (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 1.5
                    * tau_plus
                )
                + np.exp(
                    (
                        gamma_minus
                        + gamma_plus
                        + (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2)
                        ** 0.5
                    )
                    * tau_plus
                )
                * (
                    0.5
                    * (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 1
                    - 0.5
                    * (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 1.5
                    * tau_plus
                )
            )
        )
    )

    mp = (
        (1 / ((gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 2.5))
        * np.exp(-2 * (gamma_minus + gamma_plus) * tau_minus)
        * (
            (
                -1
                * np.exp(
                    (
                        gamma_minus
                        + gamma_plus
                        - (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2)
                        ** 0.5
                    )
                    * tau_minus
                )
                + 1
                * np.exp(
                    (
                        gamma_minus
                        + gamma_plus
                        + (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2)
                        ** 0.5
                    )
                    * tau_minus
                )
            )
            * gamma_plus
            * (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 2
            * tau_minus
            + (
                -1
                * np.exp(
                    (
                        gamma_minus
                        + gamma_plus
                        - (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2)
                        ** 0.5
                    )
                    * tau_minus
                )
                - 1
                * np.exp(
                    (
                        gamma_minus
                        + gamma_plus
                        + (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2)
                        ** 0.5
                    )
                    * tau_minus
                )
            )
            * (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 2.5
            * tau_minus
            + gamma_minus
            * (1 * gamma_minus - 2 * gamma_plus)
            * (
                np.exp(
                    (
                        gamma_minus
                        + gamma_plus
                        + (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2)
                        ** 0.5
                    )
                    * tau_minus
                )
                * (
                    -0.25
                    * (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 1
                    + 0.25
                    * (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 1.5
                    * tau_minus
                )
                + np.exp(
                    (
                        gamma_minus
                        + gamma_plus
                        - (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2)
                        ** 0.5
                    )
                    * tau_minus
                )
                * (
                    0.25
                    * (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 1
                    + 0.25
                    * (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 1.5
                    * tau_minus
                )
            )
        )
    )

    mm = (
        (1 / (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 2.5)
        * np.exp(-2 * (gamma_minus + gamma_plus) * tau_minus)
        * (
            (
                -1
                * np.exp(
                    (
                        gamma_minus
                        + gamma_plus
                        - (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2)
                        ** 0.5
                    )
                    * tau_minus
                )
                - 1
                * np.exp(
                    (
                        gamma_minus
                        + gamma_plus
                        + (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2)
                        ** 0.5
                    )
                    * tau_minus
                )
            )
            * (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 2.5
            * tau_minus
            + np.exp(
                (
                    gamma_minus
                    + gamma_plus
                    + (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 0.5
                )
                * tau_minus
            )
            * (
                gamma_minus
                * (0.5 * gamma_minus - 0.25 * gamma_plus)
                * (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2)
                + gamma_minus
                * (-0.5 * gamma_minus + 0.25 * gamma_plus)
                * (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 1.5
                * tau_minus
                + (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 2
                * (-0.5 + 1.5 * gamma_minus * tau_minus - 0.5 * gamma_plus * tau_minus)
            )
            + np.exp(
                (
                    gamma_minus
                    + gamma_plus
                    - (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 0.5
                )
                * tau_minus
            )
            * (
                (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 2
                * (0.5 - 1.5 * gamma_minus * tau_minus + 0.5 * gamma_plus * tau_minus)
                + gamma_minus
                * (-0.5 * gamma_minus + 0.25 * gamma_plus)
                * (
                    (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2)
                    + (gamma_minus**2 - gamma_minus * gamma_plus + gamma_plus**2) ** 1.5
                    * tau_minus
                )
            )
        )
    )

    num = (
        (gamma_minus * mm) ** 2
        + (gamma_minus * pm) ** 2
        + (gamma_plus * mp) ** 2
        + (gamma_plus * pp) ** 2
    )
    den = (pm * mp - pp * mm) ** 2
    T = 2 * repetitions * (tau_plus + tau_minus)

    # take log of cost function to maintain computational sanity
    log_cost_function = (
        0.5 * (np.log(T) + np.log(num) - np.log(den))
        - np.log(gamma_minus)
        - np.log(gamma_plus)
    )

    tau_optimized = [tp[np.argmin(log_cost_function)], tm[np.argmin(log_cost_function)]]
    print("tau_optimized ", tau_optimized)

    return tau_optimized


def calculate_M_tildes(gamma_grid, tau_grid):

    gamma_plus, gamma_minus = gamma_grid
    tau_plus, tau_minus = tau_grid

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

    return [M_plus_tilde, M_plus_tilde]


def fake_counts(tau, num_samples, T1):

    M = np.random.normal((np.exp(-(tau / T1))), 1e-4, num_samples)

    return M


def calculate_likelihood(M_measured, tau_opt, gamma_grid, n_gamma):

    M_plus_measured = np.mean(M_measured[0]) * np.ones((n_gamma, n_gamma))
    sigma_M_plus_measured = np.std(M_measured[0])
    # print('sigma_M_plus_measured ', sigma_M_plus_measured)

    M_minus_measured = np.mean(M_measured[1]) * np.ones((n_gamma, n_gamma))
    sigma_M_minus_measured = np.std(M_measured[1])

    tau_plus, tau_minus = tau_opt[0], tau_opt[1]
    M_plus_tilde, M_minus_tilde = calculate_M_tildes(gamma_grid, tau_opt)
    # print('M_plus_measured ', M_plus_measured)
    # print('M_plus_tilde ', M_plus_tilde)

    chi_plus = (M_plus_measured - M_plus_tilde) / ((2**0.5) * (sigma_M_plus_measured))
    chi_minus = (M_minus_measured - M_minus_tilde) / (
        (2**0.5) * (sigma_M_minus_measured)
    )

    chi_sq = (chi_plus**2) + (chi_minus**2)
    chi_sq_final = chi_sq - np.min(chi_sq)
    likelihood = np.exp(-(chi_sq_final))
    # print('likelihood ', likelihood)

    return likelihood


def printing_and_plotting(gamma_grid, prior_gamma, gamma_plus_arr, gamma_minus_arr):

    gamma_plus_distr = np.sum(prior_gamma, 0)
    # print('gamma_plus_distr ', gamma_plus_distr)
    gamma_minus_distr = np.sum(prior_gamma, 1)
    # print('gamma_minus_distr ', gamma_minus_distr)

    # PRINT GAMMA_PLUS AND GAMMA_MINUS FROM CURRENT PDF
    print("T1_plus estimate in us: ", 1000 / (gamma_plus_arr[np.argmax(gamma_plus_distr)]))
    print(
        "T1_minus estimate in us: ", 1000 / (gamma_minus_arr[np.argmax(gamma_minus_distr)])
    )

    # PLOT GAMMA_PLUS AND GAMMA_MINUS PDFs
    fig, axs = plt.subplots(2)
    fig.suptitle("Gamma_plus and Gamma_minus (in us^-1) distributions")
    axs[0].plot(gamma_plus_arr, gamma_plus_distr)
    axs[1].plot(gamma_minus_arr, gamma_minus_distr)
    plt.show()


"""
OBE type estimation functions

def calculate_gamma_params(
    prior_gamma,
    gamma_grid,
    gamma_plus_arr,
    gamma_minus_arr,
    n_gamma,
    delta_gamma,
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
    gamma_sigma = [gamma_sigma_plus, gamma_sigma_minus]
    print('gamma_mean ', gamma_mean)
    # print('gamma_sigma ', gamma_sigma)

    return gamma_mean, gamma_sigma


def zeroth_calculate_tau_opt(tau_grid, gamma_params, repetitions):

    gamma_mean, gamma_sigma = gamma_params
    tau_plus, tau_minus = tau_grid
    tp, tm = tau_plus.flatten(), tau_minus.flatten()

    cost_function = (
        (gamma_sigma[0] / gamma_mean[0]) ** 2
        + (gamma_sigma[1] / gamma_mean[1] ** 2) ** 0.5
    ) * (2 * repetitions * (tau_plus + tau_minus)) ** 0.5
    # print('cost_function ', cost_function)

    tau_optimized = [tp[np.argmin(cost_function)], tm[np.argmin(cost_function)]]
    print('tau_optimized ', tau_optimized)

    return tau_optimized

"""
