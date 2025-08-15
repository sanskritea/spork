import numpy as np
import emcee
import corner

class MCMCFitter:
    def __init__(self, data, errors):
        """
        Initialize the MCMC fitter
        
        Parameters:
        data (array): Observed data points
        errors (array): Measurement uncertainties
        """
        self.data = data
        self.errors = errors
        
    def model(self, params):
        """
        Your 6-parameter model function.
        Modify this to match your specific model.
        """
        param1, param2, param3, param4, param5, param6 = params
        # Replace this with your actual model calculation
        return param1 * np.exp(-param2 * self.data) + param3 * np.sin(param4 * self.data + param5) + param6
    
    def log_prior(self, params):
        """
        Log of prior probability distribution.
        Define parameter bounds and priors here.
        """
        param1, param2, param3, param4, param5, param6 = params
        
        # Example of uniform priors with bounds
        if (0 < param1 < 10 and 
            0 < param2 < 5 and 
            -5 < param3 < 5 and 
            0 < param4 < 2*np.pi and 
            0 < param5 < 2*np.pi and 
            -10 < param6 < 10):
            return 0.0  # Log of 1
        return -np.inf  # Log of 0
    
    def log_likelihood(self, params):
        """
        Log-likelihood function assuming Gaussian errors
        """
        model_prediction = self.model(params)
        return -0.5 * np.sum(((self.data - model_prediction) / self.errors) ** 2)
    
    def log_probability(self, params):
        """
        Log probability = log prior + log likelihood
        """
        lp = self.log_prior(params)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.log_likelihood(params)
    
    def run_mcmc(self, nwalkers=32, nsteps=10000, initial_guess=None):
        """
        Run the MCMC sampling
        """
        # Number of parameters
        ndim = 6
        
        # Set up initial positions of walkers
        if initial_guess is None:
            initial_guess = [5., 2.5, 0., np.pi, np.pi, 0.]  # Modify based on your problem
            
        # Add small random offsets for each walker
        pos = initial_guess + 1e-4 * np.random.randn(nwalkers, ndim)
        
        # Initialize and run sampler
        sampler = emcee.EnsembleSampler(nwalkers, ndim, self.log_probability)
        sampler.run_mcmc(pos, nsteps, progress=True)
        
        return sampler
    
    def analyze_chains(self, sampler, burnin=1000, params_of_interest=[0, 1]):
        """
        Analyze the MCMC chains and extract parameters of interest
        """
        # Remove burn-in phase
        samples = sampler.get_chain(discard=burnin, flat=True)
        
        # Extract parameters of interest
        selected_params = samples[:, params_of_interest]
        
        # Calculate statistics
        medians = np.median(selected_params, axis=0)
        percentiles = np.percentile(selected_params, [16, 84], axis=0)
        
        return selected_params, medians, percentiles

# Example usage:
if __name__ == "__main__":
    # Generate synthetic data
    x = np.linspace(0, 10, 100)
    true_params = [5.0, 2.0, 1.0, 1.5, 0.5, 2.0]
    y_true = true_params[0] * np.exp(-true_params[1] * x) + \
             true_params[2] * np.sin(true_params[3] * x + true_params[4]) + true_params[5]
    errors = 0.1 * np.ones_like(x)
    y_obs = y_true + errors * np.random.normal(size=len(x))
    
    # Initialize and run MCMC
    fitter = MCMCFitter(y_obs, errors)
    sampler = fitter.run_mcmc()
    
    # Analyze results for parameters 0 and 1
    samples, medians, percentiles = fitter.analyze_chains(sampler, params_of_interest=[0, 1])
    
    # Plot results using corner plot
    import corner
    figure = corner.corner(samples, labels=["param1", "param2"],
                          truths=[true_params[0], true_params[1]])