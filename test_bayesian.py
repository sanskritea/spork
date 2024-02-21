"""

Use this code with a two-measurement pulse sequence to perform bayesian adaptive relaxometry

Reference: Caouette, Childress et al, DOI: 10.1103/PhysRevApplied.17.064031

Sanskriti Chitransh, Aditya Vijaykumar (CITA)

"""

import time

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import sympy
from sympy import Symbol
from sympy import lambdify



GAMMA_SIM = [4, 11]  # [gamma_plus_sim, gamma_minus_sim]


def find_bias(max_cycle, num_cycle):

    cycles_range = np.unique(np.geomspace(1, max_cycle, num_cycle).astype(int))
    cycles = len(cycles_range)
    gamma_plus_est = np.zeros(cycles)
    gamma_minus_est = np.zeros(cycles)
    gamma_plus_delta = np.zeros(cycles)
    gamma_minus_delta = np.zeros(cycles)
    time_taken = np.zeros(cycles)

    #  find gamma_plus and gamma_minus for different cycle numbers
    for n in range(cycles):

        bayesian_cycle_num = cycles_range[n]
        print('Running ', bayesian_cycle_num, 'bayesian cycles')
        start_time = time.time()

        gamma_arr, gamma_pdf = BayesianT1(bayesian_cycle_num)
        time_taken[n]   = time.time() - start_time

        gamma_plus = gamma_arr[np.argmax(np.sum(gamma_pdf, 0))]
        gamma_plus_est[n] = gamma_plus
        gamma_plus_delta[n] = np.abs(GAMMA_SIM[0] - gamma_plus)

        gamma_minus = gamma_arr[np.argmax(np.sum(gamma_pdf, 1))]
        gamma_minus_est[n] = gamma_minus
        gamma_minus_delta[n] = np.abs(GAMMA_SIM[1] - gamma_minus)

    # plot gamma trends
    fig, axes = plt.subplots(3)
    fig.suptitle("Gamma Trends with Adaptive Cycles")

    axes[0].plot(cycles_range, gamma_plus_est, label = 'estimated gamma_plus')
    axes[0].plot(cycles_range, gamma_minus_est, label = 'estimated gamma_minus')
    axes[0].axhline(x=GAMMA_SIM[0], color="r", linestyle="--", label = 'real gamma_plus')
    axes[0].axhline(x=GAMMA_SIM[1], color="b", linestyle="--", label = 'real gamma_minus')
    axes[0].set_title('Gamma VS Bayesian Cycles')
    axes[0].legend()

    axes[1].plot(cycles_range, gamma_plus_delta, label = 'gamma_plus deviation')
    axes[1].plot(cycles_range, gamma_minus_delta, label = 'gamma_minus deviation')
    axes[1].axhline(x=0, color="r", linestyle="--", label = 'zero')
    axes[1].set_title('Deviation from Real Gamma VS Bayesian Cycles')
    axes[1].legend()

    axes[2].plot(cycles_range, time_taken)
    axes[2].set_title('Time taken for Bayesian Cycles')

    plt.show()


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
    gamma_upper=32,  # in ms^-1, decided after calculating const function many times
    n_gamma=1000,
    tau_lower=0.003,  # in ms
    tau_upper=5.5,  # in ms
    n_tau=1000,
    repetitions=1000,
):

    # decay rates grid
    gamma_plus_arr = np.linspace(gamma_lower, gamma_upper, n_gamma)
    gamma_minus_arr = np.linspace(gamma_lower, gamma_upper, n_gamma)
    gamma_grid = np.meshgrid(gamma_plus_arr, gamma_minus_arr)
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

        # print("Doing adaptive cycle", num)

        # find optimized taus from NOB cost function
        gamma_plus, gamma_minus = calc_mean_gammas(prior_gamma, gamma_grid, delta_gamma)

        tau_opt = nob_calculate_tau_opt(tau_grid, repetitions, gamma_plus, gamma_minus)

        # use taus in measurement
        M_measured = fake_counts(tau_opt, repetitions, GAMMA_SIM)

        # calculate likelihood from measurement result
        M_measured_mean = [np.mean(M_measured[0]), np.mean(M_measured[1])]
        M_measured_sigma = [np.std(M_measured[0]), np.std(M_measured[1])]

        likelihood = calculate_likelihood(
            M_measured_mean, M_measured_sigma, tau_opt, gamma_grid
        )

        # calculate posterior
        posterior_gamma_unnorm = likelihood * prior_gamma
        posterior = normalize_2D_pdf(posterior_gamma_unnorm, delta_gamma, delta_gamma)

        # PLOT PRIOR, LIKELIHOOD, POSTERIOR
        # plot_pdfs(prior_gamma, likelihood, posterior)

        # normalize posterior and update prior
        prior_gamma = posterior

    # printing_and_plotting(gamma_grid, prior_gamma, gamma_plus_arr, gamma_minus_arr)
    return (gamma_plus_arr, prior_gamma)


def normalize_2D_pdf(pdf, delta_x, delta_y):

    return pdf / (np.sum(pdf) * delta_x * delta_y)


def normalize_1D_pdf(pdf, delta_x):

    return pdf / (np.sum(pdf) * delta_x)


def calc_mean_gammas(prior, gamma_grid, delta_gamma):

    return (
        np.sum(prior * gamma_grid[0]) * delta_gamma**2,
        np.sum(prior * gamma_grid[1]) * delta_gamma**2,
    )


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

    T = 2 * repetitions * (tau_plus + tau_minus)

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

    return tp[min_cost_idx], tm[min_cost_idx]


def fake_counts(tau, num_samples, gamma):

    means = calculate_M_tildes(gamma, tau)
    stds = [1e-3, 1e-3]

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


def printing_and_plotting(gamma_grid, prior_gamma, gamma_plus_arr, gamma_minus_arr):

    gamma_plus_distr = np.sum(prior_gamma, 0)
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


"""
TEST FUNCTIONS


def plot_pdfs(prior, likelihood, posterior):

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

    # gamma_plus_prior = axes[3].plot(gamma_plus_arr, np.sum(prior, 0))
    # gamma_minus_prior = axes[3].plot(gamma_plus_arr, np.sum(prior, 1))
    # axes[3].set_title("Gamma priors")

    # gamma_plus_posterior = axes[4].plot(gamma_plus_arr, np.sum(posterior, 0))
    # gamma_minus_posterior = axes[4].plot(gamma_plus_arr, np.sum(posterior, 1))
    # axes[4].set_title("Gamma posteriors")

    # plt.show()


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
