import numpy as np
import emcee
import matplotlib.pyplot as plt

# Simulated data for demonstration (replace with your data)
np.random.seed(42)
x_data = np.linspace(0, 10, 10)
true_m = 2.0
true_c = 1.0
sigma = 0.5
y_data = true_m * x_data + true_c + np.random.normal(0, sigma, len(x_data))  # add some noise

# Set up the likelihood, prior, and posterior

def log_likelihood(theta, x, y, sigma=sigma):
    m, c = theta
    model = m * x + c
    chi_sq = -0.5 * np.sum(((y - model) / sigma) ** 2)
    return chi_sq

def log_prior(theta):
    m, c = theta
    if -10.0 < m < 10.0 and -10.0 < c < 10.0:  # Flat priors for m and c
        return 0.0
    return -np.inf

def log_posterior(theta, x, y, sigma=sigma):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, sigma)

# Set up the sampler (initial guesses for m and c)
initial_guess = [0.0, 0.0]
ndim = 2  # Number of parameters (m, c)
nwalkers = 50
nsteps = 1000

# Starting position for the walkers (slightly varied around initial guess)
p0 = [initial_guess + 1e-4 * np.random.randn(ndim) for i in range(nwalkers)]

# Start the emcee sampler for the first data point and gradually add more data points
samples_list = []

for n_data_points in range(1, len(x_data) + 1):
    print(f"Using first {n_data_points} data points")

    # Slice the data to use the first 'n_data_points'
    x_sub = x_data[:n_data_points]
    y_sub = y_data[:n_data_points]
    
    # Run the MCMC sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(x_sub, y_sub, sigma))
    sampler.run_mcmc(p0, nsteps, progress=True)
    
    # Store the samples after the burn-in
    burn_in = 500  # Discard first 200 steps
    samples = sampler.get_chain(discard=burn_in, flat=True)
    samples_list.append(samples)

    # Optionally, plot or analyze the result at each step
    m_samples = samples[:, 0]
    c_samples = samples[:, 1]
    
    plt.figure(figsize=(5, 3))
    plt.scatter(x_data, y_data, label="Data", color="black", alpha=0.5)
    plt.plot(x_data, np.mean(m_samples) * x_data + np.mean(c_samples), color="red", label="Best fit line")
    plt.title(f"Fit with first {n_data_points} data points")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.show()

    plt.figure(figsize=(3, 3))
    plt.scatter(m_samples, c_samples, color="blue", s=0.1, alpha=0.5)
    plt.xlabel("m")
    plt.ylabel("c")
    plt.xlim(-10, 10)
    plt.ylim(-10, 10)
    plt.title(f"Parameter space with first {n_data_points} data points")
    plt.show()


# At the end, you will have a list of samples for each incremental step
# You can analyze them further or plot the evolution of the best-fit line over time