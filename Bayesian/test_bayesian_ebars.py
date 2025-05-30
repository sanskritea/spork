"""

Use this code with a two-measurement pulse sequence to perform bayesian adaptive relaxometry

Reference: Caouette, Childress et al, DOI: 10.1103/PhysRevApplied.17.064031

Sanskriti Chitransh, Aditya Vijaykumar (CITA)

Command line statement:
python -c 'import test_bayesian; test_bayesian.BayesianT1(N_bayesian)'

"""

import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d
import sympy
from sympy import Symbol
from sympy import lambdify


GAMMA_SIM = np.array([1, 2])  # [gamma_plus_sim, gamma_minus_sim]


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
    gamma_lower=0.055,  # in ms^-1
    gamma_upper=10,     # in ms^-1, decided after calculating const function many times, original value 32
    n_gamma=1001,       # original value 3001
    tau_lower=0.003,    # in ms
    tau_upper=10,       # in ms, original value 5.5
    n_tau=2000,         # original value 1000
    repetitions=1000,
):

    # decay rates grid
    gamma_plus_arr = np.linspace(gamma_lower, gamma_upper, n_gamma)
    gamma_minus_arr = np.linspace(gamma_lower, gamma_upper, n_gamma)
    gamma_grid = np.meshgrid(gamma_plus_arr, gamma_minus_arr)
    delta_gamma = gamma_plus_arr[1] - gamma_plus_arr[0]

    # decay rates distribution
    gamma_distr = np.ones((n_gamma, n_gamma))
    # gamma_distr = normalize_2D_pdf(gamma_distr, delta_gamma, delta_gamma)
    gamma_distr = trap_normalize_2D_pdf(gamma_distr, gamma_grid)

    # relaxometry delay tau grid
    tau_plus_arr = np.geomspace(tau_lower, tau_upper, n_tau)
    tau_minus_arr = np.geomspace(tau_lower, tau_upper, n_tau)
    tau_grid = np.meshgrid(tau_plus_arr, tau_minus_arr)

    # begin with a flat prior in gammas
    prior_gamma = gamma_distr.copy()
    tau_plus, tau_minus = 0, 0

    # plotting arrays
    cycles_range = np.array(np.linspace(1, N_bayesian, N_bayesian).astype(int))
    gamma_plus_est = np.zeros(N_bayesian)
    gamma_minus_est = np.zeros(N_bayesian)
    gamma_plus_pointzerofive = np.zeros(N_bayesian)
    gamma_plus_pointninefive = np.zeros(N_bayesian)
    gamma_minus_pointzerofive = np.zeros(N_bayesian)
    gamma_minus_pointninefive = np.zeros(N_bayesian)
    tau_plus_final = np.zeros(N_bayesian)
    tau_minus_final = np.zeros(N_bayesian)
    time_taken = np.zeros(N_bayesian)

    for num in range(N_bayesian):

        start_time = time.time()

        print("Doing adaptive cycle", num)

        ## find optimized taus from NOB cost function
        # calculate mean gammas
        # gamma_plus, gamma_minus = calc_mean_gammas(prior_gamma, gamma_grid, delta_gamma)

        # use cdf to calculate mean gammas from asymmetric priors
        gamma_vals = calculate_conf_intervals(prior_gamma, gamma_grid, delta_gamma)
        gamma_plus = gamma_vals[0][1]
        gamma_minus = gamma_vals[1][1]

        tau_opt = nob_calculate_tau_opt(tau_grid, repetitions, gamma_plus, gamma_minus)
        tau_plus, tau_minus = tau_opt

        ## use taus in measurement
        M_measured = fake_counts(tau_opt, repetitions, GAMMA_SIM)
        M_measured_mean = [np.mean(np.transpose(M_measured)[0]), np.mean(np.transpose(M_measured)[1])]
        M_measured_sigma = [np.std(np.transpose(M_measured)[0]), np.std(np.transpose(M_measured)[1])]

        # calculate likelihood from measurement result
        likelihood = calculate_likelihood(
            M_measured_mean, M_measured_sigma, tau_opt, gamma_grid
        )

        # calculate posterior
        posterior_gamma_unnorm = likelihood * prior_gamma
        # posterior = normalize_2D_pdf(posterior_gamma_unnorm, delta_gamma, delta_gamma)

        # normalize posterior
        posterior = trap_normalize_2D_pdf(posterior_gamma_unnorm, gamma_grid)

        # update prior
        prior_gamma = posterior

        # prepare plotting arrays
        gamma_plus_est[num] = gamma_plus
        gamma_plus_pointzerofive[num] = gamma_vals[0][0]
        gamma_plus_pointninefive[num] = gamma_vals[0][2]
        gamma_minus_est[num] = gamma_minus
        gamma_minus_pointzerofive[num] = gamma_vals[1][0]
        gamma_minus_pointninefive[num] = gamma_vals[1][2]
        tau_plus_final[num] = tau_plus
        tau_minus_final[num] = tau_minus
        time_taken[num] = time.time() - start_time

        # if num == 0 or num == N_bayesian - 1:
        #     printing_and_plotting(gamma_grid, prior_gamma)

    # plot trends
    fig, axes = plt.subplots(1, 5, figsize=(10, 2.5))
    fig.suptitle("Trends with Adaptive Cycles")

    axes[0].plot(cycles_range, gamma_plus_est)
    axes[0].fill_between(
        cycles_range,
        gamma_plus_pointzerofive,
        gamma_plus_pointninefive,
        alpha=0.2,
    )
    axes[0].axhline(y=GAMMA_SIM[1], color="black", linestyle="--")
    axes[0].axhline(y=GAMMA_SIM[0], color="red", linestyle="--")
    axes[0].set_title("gamma_plus", fontsize=10)
    axes[0].set_xlabel("cycles", fontsize=9)
    axes[0].set_ylabel("ms^-1", fontsize=9)

    axes[1].plot(cycles_range, gamma_minus_est)
    axes[1].fill_between(
        cycles_range,
        gamma_minus_pointzerofive,
        gamma_minus_pointninefive,
        alpha=0.2,
    )
    axes[1].axhline(y=GAMMA_SIM[0], color="black", linestyle="--")
    axes[1].axhline(y=GAMMA_SIM[1], color="red", linestyle="--")
    axes[1].set_title("gamma_minus", fontsize=10)
    axes[1].set_xlabel("cycles", fontsize=9)
    axes[1].set_ylabel("ms^-1", fontsize=9)

    axes[2].plot(cycles_range, tau_plus_final * 1000)
    axes[2].set_title("tau_plus", fontsize=10)
    axes[2].set_xlabel("cycles", fontsize=9)
    axes[2].set_ylabel("us", fontsize=9)

    axes[3].plot(cycles_range, tau_minus_final * 1000)
    axes[3].set_title("tau_minus", fontsize=10)
    axes[3].set_xlabel("cycles", fontsize=9)
    axes[3].set_ylabel("us", fontsize=9)

    time_taken = np.cumsum(time_taken)
    axes[4].plot(cycles_range, time_taken)
    axes[4].set_title("Time taken", fontsize=10)
    axes[4].set_xlabel("cycles", fontsize=9)
    axes[4].set_ylabel("s", fontsize=9)

    plt.tight_layout()
    plt.show()


def normalize_2D_pdf(pdf, delta_x, delta_y):

    return pdf / (np.sum(pdf) * delta_x * delta_y)


def trap_normalize_2D_pdf(pdf, gamma_grid):

    gamma_plus, gamma_minus = gamma_grid[0][0], np.transpose(gamma_grid[1])[0]
    norm = np.trapz(np.trapz(pdf, gamma_minus), gamma_plus)
    print('norm ', norm)
    return pdf / norm


def calculate_conf_intervals(prior_gamma, gamma_grid, delta_gamma):

    gamma_plus = gamma_grid[0][0]
    gamma_minus = np.transpose(gamma_grid[1])[0]

    gamma_plus_distr = np.sum(prior_gamma, 0) * delta_gamma
    norm_gamma_plus_distr = gamma_plus_distr / (np.sum(gamma_plus_distr) * delta_gamma)
    gamma_minus_distr = np.sum(prior_gamma, 1) * delta_gamma
    norm_gamma_minus_distr = gamma_minus_distr / (
        np.sum(gamma_minus_distr) * delta_gamma
    )

    gamma_plus_cdf = np.cumsum(norm_gamma_plus_distr) * delta_gamma
    gamma_minus_cdf = np.cumsum(norm_gamma_minus_distr) * delta_gamma

    gamma_plus_cdf_interp = interp1d(gamma_plus_cdf, gamma_plus)
    gamma_minus_cdf_interp = interp1d(gamma_minus_cdf, gamma_minus)

    gamma_plus_vals = gamma_plus_cdf_interp(np.array([0.05, 0.5, 0.95]))

    gamma_minus_vals = gamma_minus_cdf_interp(np.array([0.05, 0.5, 0.95]))

    return gamma_plus_vals, gamma_minus_vals


def calc_mean_gammas(prior_gamma, gamma_grid, delta_gamma):

    return (
        np.sum(prior_gamma * gamma_grid[0]) * delta_gamma**2,
        np.sum(prior_gamma * gamma_grid[1]) * delta_gamma**2,
    )


def calc_std_gammas(prior_gamma, gamma_grid, delta_gamma):

    gamma_plus_mean, gamma_minus_mean = calc_mean_gammas(
        prior_gamma, gamma_grid, delta_gamma
    )

    gamma_plus_distr = np.sum(prior_gamma, 0)
    gamma_minus_distr = np.sum(prior_gamma, 1)

    gamma_plus = gamma_grid[0][0]
    gamma_minus = np.transpose(gamma_grid[1])[0]

    gamma_plus_std = np.sqrt(
        np.sum(np.square(gamma_plus - gamma_plus_mean) * gamma_plus_distr * delta_gamma)
    )
    gamma_minus_std = np.sqrt(
        np.sum(
            np.square(gamma_minus - gamma_minus_mean) * gamma_minus_distr * delta_gamma
        )
    )

    return (gamma_plus_std, gamma_minus_std)


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
    T_overhead = 0

    num = (
        (gamma_minus * mm(gamma_plus, gamma_minus, tau_minus)) ** 2
        + (gamma_minus * pm(gamma_plus, gamma_minus, tau_plus)) ** 2
        + (gamma_plus * mp(gamma_plus, gamma_minus, tau_minus)) ** 2
        + (gamma_plus * pp(gamma_plus, gamma_minus, tau_plus)) ** 2
    )

    den = (
        pm(gamma_plus, gamma_minus, tau_plus) * mp(gamma_plus, gamma_minus, tau_minus)
    ) - pp(gamma_plus, gamma_minus, tau_plus) * mm(gamma_plus, gamma_minus, tau_minus)
    den = den**2
    print('den ', den)
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

    # plot cost function
    # plot_cost_function(cost_function, tau_minus, tau_plus)

    return np.array([tp[min_cost_idx], tm[min_cost_idx]])


def fake_counts(tau, num_samples, gamma):

    means = calculate_M_tildes(gamma, tau)
    stds = [1e-2, 1e-2]

    M = np.random.normal(means, stds, (num_samples, 2))

    return M


def calculate_likelihood(M_measured_mean, M_measured_sigma, tau_opt, gamma_grid):

    M_plus_measured, M_minus_measured = M_measured_mean
    sigma_M_plus_measured, sigma_M_minus_measured = M_measured_sigma

    M_plus_tilde, M_minus_tilde = calculate_M_tildes(gamma_grid, tau_opt)

    chi_plus = (M_plus_measured - M_plus_tilde) / ((2**0.5) * (sigma_M_plus_measured))
    chi_minus = (M_minus_measured - M_minus_tilde) / (
        (2**0.5) * (sigma_M_minus_measured)
    )

    chi_sq = (chi_plus**2) + (chi_minus**2)
    chi_sq_final = chi_sq - np.min(chi_sq)
    likelihood = np.exp(-(chi_sq_final))

    return likelihood


def calculate_M_tildes(gamma_grid, tau_opt):

    gamma_plus, gamma_minus = gamma_grid
    tau_plus, tau_minus = tau_opt

    M_plus_tilde = Mtplus(gamma_plus, gamma_minus, tau_plus)
    M_minus_tilde = Mtminus(gamma_plus, gamma_minus, tau_minus)

    return [M_plus_tilde, M_minus_tilde]


def printing_and_plotting(gamma_grid, prior_gamma):

    gamma_plus_arr = gamma_grid[0][0]
    gamma_minus_arr = np.transpose(gamma_grid[1])[0]

    gamma_plus_distr = np.sum(prior_gamma, 0)
    # print('final gamma_plus distr ', gamma_plus_distr)
    gamma_minus_distr = np.sum(prior_gamma, 1)

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
    axs[0].axvline(x=GAMMA_SIM[1], color="b", linestyle="--")
    axs[0].set_xlabel("Gamma plus")
    axs[1].plot(gamma_minus_arr, gamma_minus_distr)
    axs[1].axvline(x=GAMMA_SIM[0], color="r", linestyle="--")
    axs[1].axvline(x=GAMMA_SIM[1], color="b", linestyle="--")
    axs[1].set_xlabel("Gamma minus")
    fig.tight_layout()
    plt.show()


##############################################################################################
##############################################################################################


"""
TEST FUNCTIONS


def good_fake_counts(
    tau,
):
    
    # 1. Evaluate using R later to find the optimum value for the pillar

    # 2. Check performance for small Rs with Cauchy distribution (ratio of two uncorrelated normal distributions)
    

    tau_plus, tau_minus = tau
    gamma_plus, gamma_minus = GAMMA_SIM
    contrast = 0.5
    R = int(10)

    mean_S1_zero = 1e6
    mean_S2_zero = mean_S1_zero * (1 - contrast)
    var_S1_zero = mean_S1_zero
    var_S2_zero = mean_S2_zero
    mean_delta = mean_S1_zero - mean_S2_zero
    var_delta = var_S1_zero + var_S2_zero
    std_delta = np.sqrt(var_delta)

    Z_max = (1 / (4 * var_delta)) * (
        np.sqrt((mean_delta**2) + (8 * var_delta)) - mean_delta
    )
    std_Z = (Z_max**2) * std_delta / np.sqrt(2 - (Z_max * mean_delta))

    mean_S1_tau_plus = mean_S1_zero * np.exp(-tau_plus * gamma_plus)
    mean_S2_tau_plus = mean_S1_tau_plus * (1 - contrast)
    var_S1_tau_plus = mean_S1_tau_plus
    var_S2_tau_plus = mean_S2_tau_plus
    mean_A_n_plus = mean_S1_tau_plus - mean_S2_tau_plus
    var_A_n_plus = var_S1_tau_plus + var_S2_tau_plus
    std_A_n_plus = np.sqrt(var_A_n_plus)

    mean_S1_tau_minus = mean_S1_zero * np.exp(-tau_minus * gamma_minus)
    mean_S2_tau_minus = mean_S1_tau_minus * (1 - contrast)
    var_S1_tau_minus = mean_S1_tau_minus
    var_S2_tau_minus = mean_S2_tau_minus
    mean_A_n_minus = mean_S1_tau_minus - mean_S2_tau_minus
    var_A_n_minus = var_S1_tau_minus + var_S2_tau_minus
    std_A_n_minus = np.sqrt(var_A_n_minus)

    # # let's try drawing the numerator and denominator and taking mean and std from there
    # num_plus = np.random.normal(mean_A_n_plus, std_A_n_plus, R)
    # num_minus = np.random.normal(mean_A_n_minus, std_A_n_minus, R)
    # one_by_den = np.random.normal(Z_max, std_Z, R)

    # drawn_M_plus = num_plus * one_by_den
    # drawn_M_minus = num_minus * one_by_den

    # drawn_mean_M_plus = np.mean(drawn_M_plus)
    # drawn_mean_M_minus = np.mean(drawn_M_minus)

    # drawn_std_M_plus = np.std(drawn_M_plus)
    # drawn_std_M_minus = np.std(drawn_M_minus)

    # print('drawn_M_plus ', drawn_M_plus)
    # print('drawn_M_minus ', drawn_M_minus)
    # print('drawn_mean_M_plus ', drawn_mean_M_plus)
    # print('drawn_mean_M_minus ', drawn_mean_M_minus)
    # print('drawn_std_M_plus ', drawn_std_M_plus)
    # print('drawn_std_M_minus ', drawn_std_M_minus)

    calc_mean_M_plus = mean_A_n_plus * Z_max
    calc_mean_M_minus = mean_A_n_minus * Z_max

    calc_std_M_plus = calc_mean_M_plus * np.sqrt(
        ((std_A_n_plus / mean_A_n_plus) ** 2) + ((std_Z / Z_max) ** 2)
    )
    calc_std_M_minus = calc_mean_M_minus * np.sqrt(
        ((std_A_n_minus / mean_A_n_minus) ** 2) + ((std_Z / Z_max) ** 2)
    )

    # print('drawn_mean_M_plus ', drawn_mean_M_plus)
    # print('calc_mean_M_plus ', calc_mean_M_plus)
    # print('drawn_std_M_plus ', drawn_std_M_plus)
    # print('calc_std_M_plus ', calc_std_M_plus)

    # mean_M = [drawn_mean_M_plus, drawn_mean_M_minus]
    # std_M = [drawn_std_M_plus, drawn_std_M_minus]  

    mean_M = [calc_mean_M_plus, calc_mean_M_minus]
    std_M = [calc_std_M_plus, calc_std_M_minus] 

    # fig, axes = plt.subplots(1, 2)

    # axes[0].hist(drawn_M_plus)
    # axes[0].axvline(x = drawn_mean_M_plus, color = "black")
    # axes[0]. axvline(x = drawn_mean_M_plus - drawn_std_M_plus, color = "black", linestyle = "--")
    # axes[0]. axvline(x = drawn_mean_M_plus + drawn_std_M_plus, color = "black", linestyle = "--")
    # axes[0].axvline(x = calc_mean_M_plus, color = "red")
    # axes[0]. axvline(x = calc_mean_M_plus - calc_std_M_plus, color = "red", linestyle = "--")
    # axes[0]. axvline(x = calc_mean_M_plus + calc_std_M_plus, color = "red", linestyle = "--")

    # axes[1].hist(drawn_M_minus)
    # axes[1].axvline(x = drawn_mean_M_minus, color = "black")
    # axes[1]. axvline(x = drawn_mean_M_minus - drawn_std_M_minus, color = "black", linestyle = "--")
    # axes[1]. axvline(x = drawn_mean_M_minus + drawn_std_M_minus, color = "black", linestyle = "--")
    # axes[1].axvline(x = calc_mean_M_minus, color = "red")
    # axes[1]. axvline(x = calc_mean_M_minus - calc_std_M_minus, color = "red", linestyle = "--")
    # axes[1]. axvline(x = calc_mean_M_minus + calc_std_M_minus, color = "red", linestyle = "--")

    # plt.tight_layout()
    # plt.show()

    return mean_M, std_M


def gamma_spreads(cycle_num, rep_number):

    gamma_plus_means = np.zeros(rep_number)
    gamma_plus_stds = np.zeros(rep_number)
    gamma_minus_means = np.zeros(rep_number)
    gamma_minus_stds = np.zeros(rep_number)

    # find mean and std of gammas after 100 bayesian cycles
    for r in range(rep_number):

        gamma_pdf, gamma_grid, delta_gamma = BayesianT1(cycle_num)
        (
            gamma_plus_means[r],
            gamma_minus_means[r],
            gamma_plus_stds[r],
            gamma_minus_stds[r],
        ) = calc_std_gammas(gamma_pdf, gamma_grid, delta_gamma)

    # print('gamma_plus_means ', gamma_plus_means)
    # print('gamma_minus_means ', gamma_minus_means)
    # print('gamma_plus_stds ', gamma_plus_stds)
    # print('gamma_minus_stds ', gamma_minus_stds)

    gamma_plus_function = (gamma_plus_means - GAMMA_SIM[0]) / gamma_plus_stds
    gamma_minus_function = (gamma_minus_means - GAMMA_SIM[1]) / gamma_minus_stds
    # print('gamma_plus_function ', gamma_plus_function)
    # print('gamma_minus_function ', gamma_minus_function)

    # plot gamma distributions
    fig, axes = plt.subplots(2)
    fig.suptitle("Gamma distributions after " + str(rep_number) + " bayesian cycles")

    axes[0].hist(gamma_plus_function)
    axes[0].set_title("Gamma_plus")

    axes[1].hist(gamma_minus_function)
    axes[1].set_title("Gamma_minus")

    plt.show()


##############################################################################################


def tau_trends(max_cycle):

    cycles = len(cycles_range)
    tau_plus_final = np.zeros(cycles)
    tau_minus_final = np.zeros(cycles)

    #  find gamma_plus and gamma_minus for different cycle numbers
    for n in range(cycles):

        tau_plus_final[n], tau_minus_final[n] = BayesianT1(cycles_range[n])

    # plot tau trends
    plt.plot(cycles_range, tau_plus_final, label="tau_plus")
    plt.plot(cycles_range, tau_minus_final, label="tau_minus")
    plt.xlabel("Number of Bayesian cycles")
    plt.ylabel("tau in ms")
    plt.legend()
    plt.show()


##############################################################################################


def trends(max_cycle):

    cycles_range = np.array(np.linspace(1, max_cycle, max_cycle).astype(int))
    cycles = len(cycles_range)
    gamma_plus_est = np.zeros(cycles)
    gamma_minus_est = np.zeros(cycles)
    gamma_plus_delta = np.zeros(cycles)
    gamma_minus_delta = np.zeros(cycles)
    gamma_plus_pointzerofive = np.zeros(cycles)
    gamma_plus_pointninefive = np.zeros(cycles)
    gamma_minus_pointzerofive = np.zeros(cycles)
    gamma_minus_pointninefive = np.zeros(cycles)
    tau_plus_final = np.zeros(cycles)
    tau_minus_final = np.zeros(cycles)
    time_taken = np.zeros(cycles)

    #  find gamma_plus and gamma_minus for different cycle numbers
    for n in range(cycles):

        bayesian_cycle_num = cycles_range[n]
        print("Running ", bayesian_cycle_num, "bayesian cycles")

        start_time = time.time()
        (
            gamma_pdf,
            gamma_grid,
            delta_gamma,
            tau_plus_final[n],
            tau_minus_final[n],
            n_gamma,
        ) = BayesianT1(bayesian_cycle_num)
        time_taken[n] = time.time() - start_time

        # gamma_plus_est[n] = gamma_grid[0][0][np.argmax(np.sum(gamma_pdf, 0))]
        # gamma_minus_est[n] = np.transpose(gamma_grid[1])[0][
        #     np.argmax(np.sum(gamma_pdf, 1))
        # ]
        # gamma_plus_delta[n], gamma_minus_delta[n] = calc_std_gammas(
        #     gamma_pdf, gamma_grid, delta_gamma
        # )

        [
            [
                gamma_plus_pointzerofive[n],
                gamma_plus_est[n],
                gamma_plus_pointninefive[n],
            ],
            [
                gamma_minus_pointzerofive[n],
                gamma_minus_est[n],
                gamma_minus_pointninefive[n],
            ],
        ] = calculate_conf_intervals(gamma_pdf, gamma_grid, delta_gamma)

        print('gamma_plus_pointzerofive[n] ', gamma_plus_pointzerofive[n])
        print('gamma_plus_est[n] ', gamma_plus_est[n])
        print('gamma_plus_pointzerofive[n] / gamma_plus_est[n] ', gamma_plus_pointzerofive[n] / gamma_plus_est[n])

        if gamma_plus_pointzerofive[n] / gamma_plus_est[n] > 0.91and gamma_minus_pointzerofive[n] / gamma_minus_est[n] > 0.91:

            print('BREAKING NOW')
            # print('n ', n)
            # print('len(cycles_range) ', len(cycles_range))
            # print('n max ', len(cycles_range) - 1)
            cycles_range = np.delete(cycles_range, (n+1, cycles - 1))
            gamma_plus_est = np.delete(gamma_plus_est, (n+1, cycles - 1)) 
            gamma_minus_est = np.delete(gamma_minus_est, (n+1, cycles - 1)) 
            gamma_plus_delta = np.delete(gamma_plus_delta, (n+1, cycles - 1)) 
            gamma_minus_delta= np.delete(gamma_minus_delta, (n+1, cycles - 1)) 
            gamma_plus_pointzerofive = np.delete(gamma_plus_pointzerofive, (n+1, cycles - 1)) 
            gamma_plus_pointninefive = np.delete(gamma_plus_pointninefive, (n+1, cycles - 1))
            gamma_minus_pointzerofive = np.delete(gamma_minus_pointzerofive, (n+1, cycles - 1))
            gamma_minus_pointninefive = np.delete(gamma_minus_pointninefive, (n+1, cycles - 1)) 
            tau_plus_final = np.delete(tau_plus_final, (n+1, cycles - 1))
            tau_minus_final = np.delete(tau_minus_final, (n+1, cycles - 1))
            time_taken = np.delete(time_taken, (n+1, cycles - 1))

            break


    # store results in txt file
    fname = (
        str(np.datetime64("today")) + "_BayesianTrends_Ngamma(" + str(n_gamma) + ")_"
    )
    np.savetxt(
        fname + "GammaPlus.txt",
        gamma_plus_est,
        delimiter=" ",
        newline="\n",
        header="",
        footer="",
    )
    np.savetxt(
        fname + "GammaMinus.txt",
        gamma_minus_est,
        delimiter=" ",
        newline="\n",
        header="",
        footer="",
    )
    np.savetxt(
        fname + "GammaPlusDelta.txt",
        gamma_plus_delta,
        delimiter=" ",
        newline="\n",
        header="",
        footer="",
    )
    np.savetxt(
        fname + "GammaMinusDelta.txt",
        gamma_minus_delta,
        delimiter=" ",
        newline="\n",
        header="",
        footer="",
    )
    np.savetxt(
        fname + "TauPlus.txt",
        tau_plus_final,
        delimiter=" ",
        newline="\n",
        header="",
        footer="",
    )
    np.savetxt(
        fname + "TauMinus.txt",
        tau_minus_final,
        delimiter=" ",
        newline="\n",
        header="",
        footer="",
    )
    np.savetxt(
        fname + "TimeTaken.txt",
        time_taken,
        delimiter=" ",
        newline="\n",
        header="",
        footer="",
    )

    # plot gamma trends
    fig, axes = plt.subplots(1, 5, figsize=(10, 2.5))
    fig.suptitle("Trends with Adaptive Cycles")

    axes[0].plot(cycles_range, gamma_plus_est)
    # axes[0].fill_between(
    #     cycles_range,
    #     gamma_plus_est - gamma_plus_delta,
    #     gamma_plus_est + gamma_plus_delta,
    #     alpha=0.2,
    # )
    axes[0].fill_between(
        cycles_range,
        gamma_plus_pointzerofive,
        gamma_plus_pointninefive,
        alpha=0.2,
    )
    axes[0].axhline(y=GAMMA_SIM[0], color="black", linestyle="--")
    axes[0].set_title("gamma_plus", fontsize=10)
    axes[0].set_xlabel("cycles", fontsize=9)
    axes[0].set_ylabel("ms^-1", fontsize=9)
    # axes[0].axis('square')

    axes[1].plot(cycles_range, gamma_minus_est)
    # axes[1].fill_between(
    #     cycles_range,
    #     gamma_minus_est - gamma_minus_delta,
    #     gamma_minus_est + gamma_minus_delta,
    #     alpha=0.2,
    # )
    axes[1].fill_between(
        cycles_range,
        gamma_minus_pointzerofive,
        gamma_minus_pointninefive,
        alpha=0.2,
    )
    axes[1].axhline(y=GAMMA_SIM[1], color="black", linestyle="--")
    axes[1].set_title("gamma_minus", fontsize=10)
    axes[1].set_xlabel("cycles", fontsize=9)
    axes[1].set_ylabel("ms^-1", fontsize=9)
    # axes[1].axis('square')

    axes[2].plot(cycles_range, tau_plus_final * 1000)
    axes[2].set_title("tau_plus", fontsize=10)
    axes[2].set_xlabel("cycles", fontsize=9)
    axes[2].set_ylabel("us", fontsize=9)
    # axes[2].axis('square')

    axes[3].plot(cycles_range, tau_minus_final * 1000)
    axes[3].set_title("tau_minus", fontsize=10)
    axes[3].set_xlabel("cycles", fontsize=9)
    axes[3].set_ylabel("us", fontsize=9)
    # axes[3].axis('square')

    axes[4].plot(cycles_range, time_taken)
    axes[4].set_title("Time taken", fontsize=10)
    axes[4].set_xlabel("cycles", fontsize=9)
    axes[4].set_ylabel("s", fontsize=9)
    # axes[4].axis('square')

    plt.tight_layout()
    plt.show()


##############################################################################################


def plot_pdfs(prior_gamma, likelihood, posterior):

    # fig , axes = plt.subplots(4)
    # fig.suptitle("PDFs")

    # prior_plot = axes[0].pcolor(prior)
    # plt.colorbar(prior_plot, ax = axes[0])
    # axes[0].set_title("Prior")

    # likelihood_plot = axes[1].pcolor(likelihood)
    # plt.colorbar(likelihood_plot, ax = axes[1])
    # axes[1].set_title("Likelihood")

    # posterior_plot = axes[2].pcolor(posterior)
    # plt.colorbar(posterior_plot, ax = axes[2])
    # axes[2].set_title("Posterior")

    # gamma_plus_prior = axes[3].plot(gamma_plus_arr, np.sum(prior_gamma, 0))
    # gamma_minus_prior = axes[3].plot(gamma_plus_arr, np.sum(prior_gamma, 1))
    # axes[3].set_title("Gamma priors")

    # gamma_plus_posterior = axes[4].plot(gamma_plus_arr, np.sum(posterior, 0))
    # gamma_minus_posterior = axes[4].plot(gamma_plus_arr, np.sum(posterior, 1))
    # axes[4].set_title("Gamma posteriors")

    # plt.show()


##############################################################################################


def calculate_gamma_opt():

    # Finding decent range of gammas to prevent overflow errors

    tau_plus = 5.5  # in ms, lower limit is 0.003 ms
    tau_minus = 5.5  # in ms
    gamma_lower = 1 / 5.5  # in ms^-1
    gamma_upper = 32  # in ms^-1
    n_gamma = 1000
    gamma_plus_arr = np.geomspace(gamma_lower, gamma_upper, n_gamma)
    gamma_minus_arr = np.geomspace(gamma_lower, gamma_upper, n_gamma)
    gamma_plus, gamma_minus = np.meshgrid(gamma_plus_arr, gamma_minus_arr)

    num = (
        (gamma_minus * mm(gamma_plus, gamma_minus, tau_minus)) ** 2
        + (gamma_minus * pm(gamma_plus, gamma_minus, tau_plus)) ** 2
        + (gamma_plus * mp(gamma_plus, gamma_minus, tau_minus)) ** 2
        + (gamma_plus * pp(gamma_plus, gamma_minus, tau_plus)) ** 2
    )

    den = (
        pm(gamma_plus, gamma_minus, tau_plus) * mp(gamma_plus, gamma_minus, tau_minus)
    ) - pp(gamma_plus, gamma_minus, tau_plus) * mm(gamma_plus, gamma_minus, tau_minus)
    den = den**2

    T = 2 * 1000 * (tau_plus + tau_minus)

    # actual cost function to compare to paper
    cost_function = ((T * num / den) ** 0.5) / (gamma_minus * gamma_plus)

    print("let us turn all nans to 1e-7")
    flag = 0
    for nni in range(n_gamma):
        for nnj in range(n_gamma):
            if np.isnan(cost_function[nni][nnj]) or np.isinf(cost_function[nni][nnj]):
                # print('AHA!')
                cost_function[nni][nnj] = 0.0000001
                flag += 1

    print("number of changes made ", flag)

    print("np.max(cost_function) ", np.max(cost_function))
    print("np.min(cost_function) ", np.min(cost_function))
    plt.title("Evaluated Cost Function")
    plt.pcolor(
        gamma_minus,
        gamma_plus,
        cost_function,
        norm=LogNorm(vmin=np.min(cost_function), vmax=np.max(cost_function)),
    )
    # ,
    #     cmap="viridis_r",
    # )
    plt.xscale("log")
    plt.yscale("log")
    plt.colorbar()
    plt.show()


##############################################################################################


def plot_M_tildes():

    tau_plus_arr = np.geomspace(0.001, 10, 1000)
    tau_minus_arr = np.geomspace(0.001, 10, 1000)
    tau_grid = np.meshgrid(tau_plus_arr, tau_minus_arr)
    tau_plus, tau_minus = tau_grid

    gamma_plus = 20
    gamma_minus = 10

    M_plus_tilde = Mtplus(gamma_plus, gamma_minus, tau_plus)
    M_minus_tilde = Mtminus(gamma_plus, gamma_minus, tau_minus)

    print("np.shape(M_plus_tilde) ", np.shape(M_plus_tilde))

    plt.title("M_tildes calculated from flattened tau meshgrid")
    _, unique_id = np.unique(tau_plus.flatten(), return_index=True)
    plt.plot(tau_plus.flatten()[unique_id], M_plus_tilde.flatten()[unique_id], "o")
    _, unique_id = np.unique(tau_minus.flatten(), return_index=True)
    plt.plot(tau_minus.flatten()[unique_id], M_minus_tilde.flatten()[unique_id], "o")
    plt.xscale("log")
    plt.xlabel("taus")
    plt.show()


##############################################################################################


def plot_cost_function(cost_function, tau_minus, tua_plus):
   
    print("np.max(cost_function) ", np.max(cost_function))
    print("np.min(cost_function) ", np.min(cost_function))
    plt.title("Evaluated Cost Function")
    plt.pcolor(
        tau_minus,
        tau_plus,
        cost_function,
        norm=LogNorm(vmin=np.min(cost_function), vmax=np.max(cost_function)),
        cmap="viridis_r",
    )
    plt.xscale("log")
    plt.yscale("log")
    plt.colorbar()
    plt.show()


##############################################################################################


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


##############################################################################################


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



##############################################################################################
##############################################################################################

"""
