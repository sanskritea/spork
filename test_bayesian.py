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


import numpy as np
import matplotlib.pyplot as plt
from sympy import Symbol
import sympy
from sympy import lambdify

GAMMA_SIM = [6, 8]  # [gamma_plus_sim, gamma_minus_sim]


def give_sympy_functions():
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

    return pp, pm, mm, mp, Mtplus, Mtminus


pp, pm, mm, mp, Mtplus, Mtminus = give_sympy_functions()


def BayesianT1(
    N_bayesian,
    gamma_lower=1,  # in ms^-1
    gamma_upper=20,  # in ms^-1
    n_gamma=1000,
    tau_lower=0.001,  # in ms
    tau_upper=0.7,  # in ms
    n_tau=1000,
    repetitions=1000,
):

    gamma_plus_arr = np.linspace(gamma_lower, gamma_upper, n_gamma)
    # print(gamma_plus_arr)
    gamma_minus_arr = np.linspace(gamma_lower, gamma_upper, n_gamma)
    gamma_grid = np.meshgrid(gamma_plus_arr, gamma_minus_arr)
    # print('gamma_grid', gamma_grid)
    delta_gamma = gamma_plus_arr[1] - gamma_plus_arr[0]

    # decay rates distribution
    gamma_distr = np.ones((n_gamma, n_gamma))
    gamma_distr = normalize_2D_pdf(gamma_distr, delta_gamma, delta_gamma)

    # relaxometry delay tau grid
    tau_plus_arr = np.geomspace(tau_lower, tau_upper, n_tau)
    tau_minus_arr = np.geomspace(tau_lower, tau_upper, n_tau)
    tau_grid = np.meshgrid(tau_plus_arr, tau_minus_arr)

    # begin with a flat prior in gammas
    prior_gamma = gamma_distr.copy()
    # printing_and_plotting(gamma_grid, prior_gamma, gamma_plus_arr, gamma_minus_arr)

    for num in range(N_bayesian):

        print("Doing adaptive cycle ", num)

        # find optimized taus from NOB cost function
        gamma_plus, gamma_minus = calc_mean_gammas(prior_gamma, gamma_grid, delta_gamma)
        tau_opt = nob_calculate_tau_opt(tau_grid, repetitions, gamma_plus, gamma_minus)

        # use taus in measurement
        M_measured = fake_counts(tau_opt, repetitions, GAMMA_SIM)

        # for t in range(2):

        #     # FAKE COUNTING, REPLACE WITH REAL DATA
        #     M_measured[t] = fake_counts(tau_opt[t], repetitions, GAMMA_SIM[t])

        # calculate likelihood from measurement result
        M_measured_mean = [np.mean(M_measured[0]), np.mean(M_measured[1])]
        M_measured_sigma = [np.std(M_measured[0]), np.std(M_measured[1])]

        likelihood = calculate_likelihood(
            M_measured_mean, M_measured_sigma, tau_opt, gamma_grid
        )

        # print('prior_gamma ', prior_gamma)
        # print('likelihood ', likelihood)
        # calculate posterior
        posterior_gamma_unnorm = likelihood * prior_gamma

        # normalize posterior and update prior
        prior_gamma = normalize_2D_pdf(posterior_gamma_unnorm, delta_gamma, delta_gamma)
        # print('new prior_gamma ', prior_gamma)

    printing_and_plotting(gamma_grid, prior_gamma, gamma_plus_arr, gamma_minus_arr)


def normalize_2D_pdf(pdf, delta_x, delta_y):
    return pdf / (np.sum(pdf) * delta_x * delta_y)


def normalize_1D_pdf(pdf, delta_x):
    return pdf / (np.sum(pdf) * delta_x)


def nob_calculate_tau_opt(tau_grid, repetitions, gamma_plus, gamma_minus):
    """
    Cost function is only a function of tau, and we want to optimize it
    to find the `best' value of tau. This is done by finding the minimum
    of the cost function evaluated over a grid of tau values.

    On the other hand, the cost function is only a function of one value
    of gamma_plus and gamma_minus. These values are the average value of
    gamma under the gamma prior. gamma_plus and gamma_minus are single
    float values, not arrays.
    """

    tau_plus, tau_minus = tau_grid

    num = (gamma_minus * mm(gamma_plus, gamma_minus, tau_minus)) ** 2
    num += (gamma_minus * pm(gamma_plus, gamma_minus, tau_plus)) ** 2
    num += (gamma_plus * mp(gamma_plus, gamma_minus, tau_minus)) ** 2
    num += (gamma_plus * pp(gamma_plus, gamma_minus, tau_plus)) ** 2

    den = pm(gamma_plus, gamma_minus, tau_plus) * mp(gamma_plus, gamma_minus, tau_minus)
    den -= pp(gamma_plus, gamma_minus, tau_plus) * mm(
        gamma_plus, gamma_minus, tau_minus
    )
    den = den**2

    T = 2 * repetitions * (tau_plus + tau_minus)

    # take log of cost function to maintain computational sanity
    log_cost_function = (
        0.5 * (np.log(T) + np.log(num) - np.log(den))
        - np.log(gamma_minus)
        - np.log(gamma_plus)
    )

    # plt.pcolor(tau_plus, tau_minus, log_cost_function)
    # plt.colorbar()
    # plt.show()

    print(len(log_cost_function[np.isnan(num)]))

    tp, tm = tau_plus.flatten(), tau_minus.flatten()
    min_cost_idx = np.argmin(log_cost_function.flatten())
    print("min cost idx", min_cost_idx)

    return tp[min_cost_idx], tm[min_cost_idx]


def calculate_M_tildes(gamma_grid, tau_opt):

    gamma_plus, gamma_minus = gamma_grid
    tau_plus, tau_minus = tau_opt

    M_plus_tilde = Mtplus(gamma_plus, gamma_minus, tau_plus)
    M_minus_tilde = Mtminus(gamma_plus, gamma_minus, tau_minus)

    return [M_plus_tilde, M_minus_tilde]


def fake_counts(tau, num_samples, gamma):

    means = calculate_M_tildes(gamma, tau)
    print(gamma, tau, means)
    stds = [1e-1, 1e-1]

    M = np.random.normal(means, stds, (num_samples, 2))

    # M = np.random.normal((np.exp(-(tau * gamma * 0.001))), 1e-4, num_samples)
    # # print('M ', M)
    # # print('np.std(M) ', np.std(M))

    return M


def calculate_likelihood(M_measured_mean, M_measured_sigma, tau_opt, gamma_grid):

    M_plus_measured, M_minus_measured = M_measured_mean
    sigma_M_plus_measured, sigma_M_minus_measured = M_measured_sigma

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


def calc_mean_gammas(prior, gamma_grid, delta_gamma):
    return (
        np.sum(prior * gamma_grid[0]) * delta_gamma**2,
        np.sum(prior * gamma_grid[1]) * delta_gamma**2,
    )


def printing_and_plotting(gamma_grid, prior_gamma, gamma_plus_arr, gamma_minus_arr):

    gamma_plus_distr = np.sum(prior_gamma, 0)
    # print('gamma_plus_distr ', gamma_plus_distr)
    gamma_minus_distr = np.sum(prior_gamma, 1)
    # print('gamma_minus_distr ', gamma_minus_distr)

    # PRINT GAMMA_PLUS AND GAMMA_MINUS FROM CURRENT PDF
    # print("T1_plus estimate in us: ", 1 / (gamma_plus_arr[np.argmax(gamma_plus_distr)]))
    # print(
    #     "T1_minus estimate in us: ", 1 / (gamma_minus_arr[np.argmax(gamma_minus_distr)])
    # )

    # PLOT GAMMA_PLUS AND GAMMA_MINUS PDFs
    fig, axs = plt.subplots(2)
    fig.suptitle("Gamma_plus and Gamma_minus (in ms^-1) distributions")
    axs[0].plot(gamma_plus_arr, gamma_plus_distr)
    axs[0].axvline(x=GAMMA_SIM[0], color="r", linestyle="--")
    axs[0].set_xlabel("Gamma plus")
    axs[1].plot(gamma_minus_arr, gamma_minus_distr)
    axs[1].axvline(x=GAMMA_SIM[1], color="r", linestyle="--")
    axs[1].set_xlabel("Gamma minus")
    fig.tight_layout()
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
