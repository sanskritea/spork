import numpy as np
import emcee
import matplotlib.pyplot as plt
import sympy
from sympy import Symbol
from sympy import lambdify
import pandas as pd


# Get M and M_err from Two Meaured T1 datasets
df = pd.read_json('C://Users/awschlab/Desktop/data/250122/Bayesian T1/Bayesian T1213203final.json')

tau_plus_list = np.array(df['datasets']['TauPlus']) # convert taus to ms
tau_minus_list = np.array(df['datasets']['TauPlus']) # convert taus to ms
# print('tau_plus ', tau_plus)
M_plus = np.array(df['datasets']['M_plus'])
M_minus = np.array(df['datasets']['M_minus'])
M_plus_err = np.array(df['datasets']['M_plus_err'])
M_minus_err = np.array(df['datasets']['M_minus_err'])

S000plus = np.array(df['datasets']['S000Plus'])
S000minus = np.array(df['datasets']['S000Minus'])
Slevel00plus = np.array(df['datasets']['Slevel00Plus'])
Slevel00minus = np.array(df['datasets']['Slevel00Minus'])
Slevel0tauplus = np.array(df['datasets']['Slevel0tauPlus'])
Slevel0tauminus = np.array(df['datasets']['Slevel0tauMinus'])
S00tauplus = np.array(df['datasets']['S00tauPlus'])
S00tauminus = np.array(df['datasets']['S00tauMinus'])

print('trying skellam distribution to calculate M and errors')
print(np.shape(S000plus))

num_meas = len(tau_plus_list)

# Symbols to calculate likelihoods
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
Gval = lambdify([gamma_plus, gamma_minus], G.simplify())


def calculate_tau_opt(gamma_plus, gamma_minus, tau, repetitions = 100000, T_overhead = 0):

        # run a 2d optimizing sampler for calculating tau opt also where the cost function becomes the likelihood

        tau_plus, tau_minus = np.meshgrid(tau, tau)

        num = (
            (gamma_minus * mm(gamma_plus, gamma_minus, tau_minus)) ** 2
            + (gamma_minus * pm(gamma_plus, gamma_minus, tau_plus)) ** 2
            + (gamma_plus * mp(gamma_plus, gamma_minus, tau_minus)) ** 2
            + (gamma_plus * pp(gamma_plus, gamma_minus, tau_plus)) ** 2
        )

        den = (
             pm(gamma_plus, gamma_minus, tau_plus) *  mp(gamma_plus, gamma_minus, tau_minus)
        ) -  pp(gamma_plus, gamma_minus, tau_plus) *  mm(gamma_plus, gamma_minus, tau_minus)
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

        # print('Time taken to optimize tau ', time.time() - start_tau_optimization)

        return np.array([tp[min_cost_idx], tm[min_cost_idx]])


def M_tildes(gamma, tau_plus, tau_minus):

    gamma_plus, gamma_minus = gamma

    M_plus_tilde = Mtplus(gamma_plus, gamma_minus, tau_plus)
    M_minus_tilde = Mtminus(gamma_plus, gamma_minus, tau_minus)   
    # print('G ', Gval(gamma_plus, gamma_minus)) 
    # print('M_plus_tilde ', M_plus_tilde) 
    # print('M_minus_tilde ', M_minus_tilde)

    return M_plus_tilde, M_minus_tilde


def log_likelihood(gamma, tau, M, M_err):  # calculate log likelihood
    
    M_plus, M_minus = M
    M_plus_err, M_minus_err = M_err
    tau_plus, tau_minus = tau

    chi_sq = 0

    for i in range(len(tau_plus)):

        M_plus_tilde, M_minus_tilde = M_tildes(gamma, tau_plus[i], tau_minus[i])

        chi_plus = (M_plus[i] - M_plus_tilde) / ((2**0.5) * M_plus_err[i])
        chi_minus = (M_minus[i] - M_minus_tilde) / ((2**0.5) * M_minus_err[i])
        chi_sum = (chi_plus**2) + (chi_minus**2)
        # print('chi_sum ', chi_sum)

        chi_sq = chi_sq + (-((chi_plus**2) + (chi_minus**2)))

    return chi_sq


def log_prior(gamma):  # create flat prior for the gammas

    gamma_plus, gamma_minus = gamma
    if 1 < gamma_plus < 10 and 1 < gamma_minus < 10:
        return 0.0

    return -np.inf


def log_posterior(gamma, tau, M, M_err):

    lp = log_prior(gamma)

    if not np.isfinite(lp):
        return -np.inf

    return lp + log_likelihood(gamma, tau, M, M_err)


# Set up initial guesses for gamma and tau for the sampler
init_guess = [1, 1]  # in ms^-1
n_dim = 2  # number of paramaters to guess
n_walkers = 50
n_steps = 1000

# Starting position for walkers around the initial guess
p0 = [init_guess + (1e-4 * np.random.randn(n_dim)) for i in range(n_walkers)]

# Start the Monte Carlo sampler for the first M, M_err and gradually add more M values
gamma_samples_list = []

# for n_data_points in range(1, num_meas + 1):
for n_data_points in range(1, num_meas+1):

    # Slice data to use the first 'n' data points
    tau_plus_slice = tau_plus_list[:n_data_points]     # x1 data
    print('tau_plus_slice ', tau_plus_slice)
    tau_minus_slice = tau_minus_list[:n_data_points]      # x2 data
    M_plus_slice = M_plus[:n_data_points]   # y1 data
    M_minus_slice = M_minus[:n_data_points] # y2 data
    M_plus_err_slice = M_plus_err[:n_data_points]   # y1 data
    M_minus_err_slice = M_minus_err[:n_data_points] # y2 data

    # Run the MCMC sampler
    sampler = emcee.EnsembleSampler(
        n_walkers, n_dim, log_posterior, args = ((tau_plus_slice, tau_minus_slice), (M_plus_slice, M_minus_slice), (M_plus_err_slice, M_minus_err_slice)))
    sampler.run_mcmc(p0, n_steps, progress=True)

    # Store samples after burn-in
    burn_in = 500  # discard first 500 steps
    gamma_samples = sampler.get_chain(discard=burn_in, flat=True)
    gamma_samples_list.append(gamma_samples)

    gamma_plus_samples = gamma_samples[:, 0]
    gamma_minus_samples = gamma_samples[:, 1]

    mean_gamma_plus = np.mean(gamma_plus_samples)
    mean_gamma_minus = np.mean(gamma_minus_samples)

    print('gamma_plus ', mean_gamma_plus)
    print('gamma_minus ', mean_gamma_minus)

    # tau_opt_plus, tau_opt_minus = calculate_tau_opt(mean_gamma_plus, mean_gamma_minus, np.geomspace(0.005, 0.8, 101))
    # print('tau_opt_plus ', tau_opt_plus)
    # print('tau_opt_minus ', tau_opt_minus)

    plt.figure(figsize=(3, 3))
    plt.scatter(gamma_plus_samples, gamma_minus_samples, color="blue", s=0.1, alpha=0.5)
    plt.xlabel("gamma_plus")
    plt.ylabel("gamma_minus")
    plt.xlim(1, 10)
    plt.ylim(1, 10)
    plt.title(f"Parameter space with first {n_data_points} data points")
    plt.show()

    # plt.subplot(121)
    # plt.errorbar(tau_plus_slice, M_plus_slice, yerr = M_plus_err_slice, fmt = '.', capsize = 2)
    # plt.plot(tau_plus_slice, [Mtplus(mean_gamma_plus, mean_gamma_minus, tau_plus_slice) for tau in tau_plus_slice])
    # plt.show()
