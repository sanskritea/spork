"""
Experiment to run T1 measuerment on the NV-AFM Spork setup


"""

import time
from itertools import count

import numpy as np
import matplotlib.pyplot as plt
import random
import scipy as scp

from nspyre import DataSource
from nspyre import InstrumentGateway

from rpyc.utils.classic import obtain
from guis.guiElements_general import flexSave

from experiments.NewPulses import Pulses
from experiments.NewportSpatialFeedback import SpatialFeedback
from drivers.ni.nidaq_final import NIDAQ
# from experiments.BayesianFunctions import BF

from scipy.interpolate import interp1d
import sympy
from sympy import Symbol
from sympy import lambdify


class Bayesian_T1_Meas:

    def trap_normalize_2D_pdf(self, pdf, gamma_grid):

        gamma_plus, gamma_minus = gamma_grid[0][0], np.transpose(gamma_grid[1])[0]
        norm = np.trapz(np.trapz(pdf, gamma_minus), gamma_plus)
        # print('norm ', norm)

        return norm


    def calculate_conf_intervals(self, prior_gamma, gamma_grid, delta_gamma):

        gamma_plus = gamma_grid[0][0]
        gamma_minus = np.transpose(gamma_grid[1])[0]

        # print('Normalizing pdfs')
        gamma_plus_distr = np.sum(prior_gamma, 0) * delta_gamma
        norm_gamma_plus_distr = gamma_plus_distr / (np.sum(gamma_plus_distr) * delta_gamma)

        gamma_minus_distr = np.sum(prior_gamma, 1) * delta_gamma
        norm_gamma_minus_distr = gamma_minus_distr / (np.sum(gamma_minus_distr) * delta_gamma)

        # print('Calculating cdfs')
        gamma_plus_cdf = np.cumsum(norm_gamma_plus_distr) * delta_gamma
        gamma_minus_cdf = np.cumsum(norm_gamma_minus_distr) * delta_gamma

        # print('Interpolating')
        gamma_plus_cdf_interp = interp1d(gamma_plus_cdf, gamma_plus) #, fill_value = 'extrapolate')
        gamma_minus_cdf_interp = interp1d(gamma_minus_cdf, gamma_minus) #, fill_value = 'extrapolate')

        gamma_plus_vals = gamma_plus_cdf_interp(np.array([0.5]))
        gamma_minus_vals = gamma_minus_cdf_interp(np.array([0.5]))

        return gamma_plus_vals, gamma_minus_vals


    def nob_calculate_tau_opt(self, tau_grid, repetitions, gamma_plus, gamma_minus, T_overhead):
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

        # print('Time taken to optimize tau ', time.time() - start_tau_optimization)

        return np.array([tp[min_cost_idx], tm[min_cost_idx]])


    def calculate_likelihood(self, M_measured_mean, M_measured_sigma, tau_opt, gamma_grid):
        
        print('tau_opt received ', tau_opt)
        start_likelihood_calculation = time.time()
        M_plus_measured, M_minus_measured = M_measured_mean
        sigma_M_plus_measured, sigma_M_minus_measured = M_measured_sigma

        M_plus_tilde, M_minus_tilde = self.calculate_M_tildes(gamma_grid, tau_opt)

        chi_plus = (M_plus_measured - M_plus_tilde) / ((2**0.5) * (sigma_M_plus_measured))
        chi_minus = (M_minus_measured - M_minus_tilde) / (
            (2**0.5) * (sigma_M_minus_measured)
        )

        chi_sq = (chi_plus**2) + (chi_minus**2)
        # chi_sq_final = chi_sq - np.min(chi_sq)
        # likelihood = np.exp(-(chi_sq_final))
        log_likelihood = - chi_sq

        # print('Time taken to calculate likelihood ', time.time() - start_likelihood_calculation)

        return log_likelihood


    def calculate_M_tildes(self, gamma_grid, tau_opt):

        gamma_plus, gamma_minus = gamma_grid
        tau_plus, tau_minus = tau_opt
        # pp, pm, mm, mp, Mtplus, Mtminus = BF.give_sympy_functions()
        
        M_plus_tilde = Symbols.Mtplus(gamma_plus, gamma_minus, tau_plus)
        M_minus_tilde = Symbols.Mtminus(gamma_plus, gamma_minus, tau_minus)

        return [M_plus_tilde, M_minus_tilde]


    def fit_data_to_gaussian(self, data_list):

        mean, stdev = scipy.stats.norm.fit(data_list)

        return mean, stdev




    def BayesianT1(
        self,
        datasetName: str,
        samplingFreq: float,
        maxIterations: int,
        freq: float,
        rf_power: float,
        laser_power: float,
        num_samples: int,
        clock_time: int,        # DAQ counting trigger pulse (from Swabian) duration, ~ 10ns
        init_time: int,         # NV initialization laser ON duration
        laser_lag: int,         # laser stabilizing time, usually ~100ns
        probe_time: int,        # readout laser ON time
        singlet_decay: int,     # NV singlet state emptying duration
        # x_init_position: float, 
        # y_init_position: float,
        # z_init_position: float,
        bayesian_iterations: int,
        pi_time: int,           # ideally should be pi_x/y for gamma_plus/minus but let's run with this until IQ works        
        # gamma_lower: float,   # for now, let's use the preset gamma bounds
        # gamma_upper: float,
    ):

        # connect to the instrument server
        # connect to the data server and create a data set, or connect to an
        # existing one with the same name if it was created earlier.
        readout_time = probe_time
        with InstrumentGateway() as gw, DataSource(datasetName) as BayesianT1Data:

            # setup a relevant iterator depending on maxIterations
            if maxIterations < 0:
                iters = count()  # infinite iterator
            else:
                iters = range(maxIterations)

            freq_list = [2.914e9, 2.828e9]
            pi_time_list = [30, 60]

            # PLOTTING VARIABLES
            self.iteration_number = np.zeros(0)
            self.gamma_plus_list = np.zeros(0)
            self.gamma_minus_list = np.zeros(0)
            self.M_plus = np.zeros(0)
            self.M_minus = np.zeros(0)
            tau_plus_list = np.zeros(bayesian_iterations)
            tau_minus_list = np.zeros(bayesian_iterations)
            S_0_0_0_list = np.zeros((2, bayesian_iterations))
            S_0_0_tau_list = np.zeros((2, bayesian_iterations))
            S_level_0_0_list = np.zeros((2, bayesian_iterations))
            S_level_0_tau_list = np.zeros((2, bayesian_iterations))

            # SETUP BAYESIAN VARIABLES
            gamma_lower = 1 # in ms^-1 (correspoonding T1: 5ms)
            gamma_upper = 10   # in ms^-1 (correspoonding T1: 100us)
            n_gamma = 1000
            tau_lower = 0.01  # in ms (10us)
            tau_upper = 5  # in ms
            n_tau = 10000
            repetitions = num_samples 

            # decay rates grid
            gamma_plus_arr = np.linspace(gamma_lower, gamma_upper, n_gamma)
            gamma_minus_arr = np.linspace(gamma_lower, gamma_upper, n_gamma)
            gamma_grid = np.meshgrid(gamma_plus_arr, gamma_minus_arr)
            delta_gamma = gamma_plus_arr[1] - gamma_plus_arr[0]

            # decay rates distribution
            gamma_distr = np.ones((n_gamma, n_gamma))
            gamma_distr_norm =  self.trap_normalize_2D_pdf(gamma_distr, gamma_grid)
            gamma_distr = gamma_distr / gamma_distr_norm

            # relaxometry delay tau grid
            tau_plus_arr = np.geomspace(tau_lower, tau_upper, n_tau)
            tau_minus_arr = np.geomspace(tau_lower, tau_upper, n_tau)
            tau_grid = np.meshgrid(tau_plus_arr, tau_minus_arr)

            tau_plus_forced = np.linspace(1000*50, 1000*500*1, bayesian_iterations)
            random.shuffle(tau_plus_forced)
            tau_minus_forced = np.linspace(1000*50, 1000*500*1, bayesian_iterations)
            random.shuffle(tau_minus_forced)
            z_list = []

            # tau_high_slope = np.space(1000*100, 1000*500, bayesian_iterations)
            # tau_high_slope = 1000 * np.random.normal(100, 0.01, bayesian_iterations)
            # print('tau_high_slope in us', 1e-3 * tau_high_slope)

            # begin with a flat prior in gammas
            prior_gamma = gamma_distr.copy()
            tau_plus, tau_minus = 0, 0
            start_time = time.time()

            # Feedback parameters
            feedback_trigger_rate = int(20e3)
            feedback_time_per_point = 0.05
            feedback_num_samples = int(feedback_trigger_rate * feedback_time_per_point)
            x_init_position = 7.4196
            y_init_position = 10.307
            z_init_position = 4.5083
            feedback_timer = time.time()
            feedback_counter = 0

            sum_log_likelihood = np.zeros((n_gamma, n_gamma))

            for num in range(bayesian_iterations):

                # # SPATIAL FEEDBACK EVERY 5 minutes
                # if ((time.time() - feedback_timer) > 300):
                #     feedback_counter = feedback_counter + 1
                #     print('Feedback')
                #     begin_feedback = time.time()
                #     SpatialFeedback.Feedback(x_init_position, y_init_position, z_init_position)
                #     feedback_duration = time.time() - begin_feedback
                #     print('Feedback duration: ', feedback_duration)
                #     print('Feedback counter ', feedback_counter)
                #     feedback_timer = time.time()

                # MAIN EXPERIMENT
                print('Bayesian iteration number ', int(num + 1))
                self.iteration_number = np.append(self.iteration_number, int(num + 1))

                # USE CDF TO CALCULATE MEAN GAMMAS FROM PRIOR DISTRIBUTIONS
                start_mean_gamma_calculation = time.time()
                print('Calculating mean gammas from priors')
                gamma_vals = self.calculate_conf_intervals(prior_gamma, gamma_grid, delta_gamma)
                gamma_plus = gamma_vals[0][0]
                gamma_minus = gamma_vals[1][0]
                # itnum = self.iteration_number[num]
                self.gamma_plus_list = np.append(self.gamma_plus_list, gamma_plus)
                self.gamma_minus_list = np.append(self.gamma_minus_list, gamma_minus)
                print('gamma_plus (ms^-1) ', gamma_plus)
                print('gamma_minus (ms^-1) ', gamma_minus)

                # GET OPTIMIZED TAUS FROM NOB TO USE IN EXPERIMENT
                # print('Optimizing taus for pulse sequences')
                # T_overhead = 0
                # tau_opt = self.nob_calculate_tau_opt(tau_grid, repetitions, gamma_plus, gamma_minus, T_overhead)
                # tau_plus, tau_minus = tau_opt
                # tau_plus_list[num] = tau_plus
                # tau_minus_list[num] = tau_minus
                # print('calculated tau_plus (ms) ', tau_plus)
                # print('calculated tau_minus (ms) ', tau_minus)

                # CREATE PULSE SEQUENCES FOR OPTIMIZED TAUS
                print('Creating pulse sequences')
                seqs = []
                # seqs.append(Pulses(gw).BAYESIAN_T1(
                #             int(1e6 * tau_plus), clock_time, init_time, readout_time, laser_lag, singlet_decay, pi_time_list[0]
                #         ))
                # seqs.append(Pulses(gw).BAYESIAN_T1(
                #             int(1e6 * tau_minus), clock_time, init_time, readout_time, laser_lag, singlet_decay, pi_time_list[1]
                #         ))

                print('tau_plus_forced[num] in ms ', 1e-6 * tau_plus_forced[num])
                print('tau_minus_forced[num] in ms ', 1e-6 * tau_minus_forced[num])
                tau_opt = 1e-6 * tau_plus_forced[num], 1e-6 * tau_minus_forced[num]
                print('tau_opt ', tau_opt)
                seqs.append(Pulses(gw).BAYESIAN_T1(
                            int(tau_plus_forced[num]), clock_time, init_time, readout_time, laser_lag, singlet_decay, pi_time_list[0]
                        ))
                seqs.append(Pulses(gw).BAYESIAN_T1(
                            int(tau_minus_forced[num]), clock_time, init_time, readout_time, laser_lag, singlet_decay, pi_time_list[1]
                        ))

                # print('tau_high_slope[num] in us', 1e-3 * tau_high_slope[num])
                # tau_opt = tau_high_slope[num] * 1e-6, tau_high_slope[num] * 1e-6
                # print('tau_opt in ms ', tau_opt)
                # seqs.append(Pulses(gw).BAYESIAN_T1(
                #             int(tau_high_slope[num]), clock_time, init_time, readout_time, laser_lag, singlet_decay, pi_time_list[0]
                #         ))
                # seqs.append(Pulses(gw).BAYESIAN_T1(
                #             int(tau_high_slope[num]), clock_time, init_time, readout_time, laser_lag, singlet_decay, pi_time_list[1]
                #         ))

                # MEASURE COUNTS
                start_counting = time.time()
                print('Measuring counts')
                M_mean_list = np.zeros(2)
                M_std_list = np.zeros(2)
                with NIDAQ() as mynidaq:

                    # SET LASER POWER
                    mynidaq.laser_power_atten(laser_power)

                    for t in range(
                        2 
                    ):  # measure NV with delay = tau_plus and then delay = tau_minus

                        # SRS ACTIONS
                        gw.sg.set_rf_amplitude(rf_power)    # set ouput power
                        gw.sg.set_frequency(freq_list[t])
                        gw.sg.set_mod_state(False) # make QAM/external after IQ is ready
                        gw.sg.set_rf_state("1")

                        # START READ TASK
                        mynidaq.start_external_read_task(samplingFreq, ((10 * num_samples) + 1))

                        # START PULSESTREAMER
                        gw.swabian.runSequenceInfinitely(seqs[t])

                        # COUNT WITH EXTERNAL TRIGGER
                        counts = obtain(mynidaq.external_read_task(samplingFreq, ((10 * num_samples) + 1)))

                        # SEPARATE COUNTS S_(initial level)_(final level)_delay
                        S_0_0_0_list = 1e9 * counts[0::8] / readout_time
                        S_0_0_tau_list = 1e9 * counts[2::8] / readout_time 
                        S_level_0_0_list = 1e9 * counts[4::8] / readout_time
                        S_level_0_tau_list = 1e9 * counts[6::8] / readout_time

                        S_0_0_0_mean = np.mean(S_0_0_0_list)
                        S_0_0_tau_mean = np.mean(S_0_0_tau_list) 
                        S_level_0_0_mean = np.mean(S_level_0_0_list) 
                        S_level_0_tau_mean = np.mean(S_level_0_tau_list) 

                        S_0_0_0_err = np.std(S_0_0_0_list) / (num_samples ** 0.5)
                        S_0_0_tau_err = np.std(S_0_0_tau_list) / (num_samples ** 0.5)
                        S_level_0_0_err = np.std(S_level_0_0_list) / (num_samples ** 0.5)
                        S_level_0_tau_err = np.std(S_level_0_tau_list) / (num_samples ** 0.5)

                        # FINAL MEASURE AND MOMENTS (Appendix E)
                        top = S_0_0_tau_mean - S_level_0_tau_mean
                        bottom = S_0_0_0_mean - S_level_0_0_mean
                        top_err = (S_0_0_tau_err ** 2 + S_level_0_tau_err ** 2) ** 0.5
                        bottom_err = (S_0_0_0_err ** 2 + S_level_0_0_err ** 2) ** 0.5

                        ratio = bottom/(bottom_err * bottom_err)
                        relative = bottom_err / bottom
                        z0 = 0.25*ratio*( ( 1 + 8*relative*relative )**0.5 - 1 )
                        z02 = z0*z0 
                        eA02 = bottom_err * bottom_err    
                        dL2dz2 = 2*ratio/(z02*z0) + 2/z02 - 3/(eA02*z02*z02) 
                        ez0 = 1/ (-dL2dz2)**0.5
                        inv_bottom = z0
                        inv_bottom_err = ez0

                        # z_list.append(inv_bottom)

                        z = top * inv_bottom
                        z_err = ((top * inv_bottom_err)**2 + (inv_bottom_err * top_err) ** 2) ** 0.5

                        M_mean = z
                        M_std = z_err
                        
                        print('M_mean ', M_mean)
                        print('M_std ', M_std)

                        M_mean_list[t] = M_mean
                        M_std_list[t] = M_std
                        
                        gw.swabian.reset()

                    print('Counting time ', time.time() - start_counting)

                # MEASUREMENT STATISTICS
                self.M_plus = np.append(self.M_plus, M_mean_list[0])
                self.M_minus = np.append(self.M_minus, M_mean_list[1])
                M_measured_mean = [M_mean_list[0], M_mean_list[1]]
                M_measured_sigma = [M_std_list[0], M_std_list[1]]

                # CALCULATE LIKELIHOOD FROM THE MEASUREMENT RESULT
                print('Calculating likelihood')
                log_likelihood = self.calculate_likelihood(M_measured_mean, M_measured_sigma, tau_opt, gamma_grid)
                sum_log_likelihood = sum_log_likelihood + log_likelihood
                plt.pcolormesh(sum_log_likelihood)
                plt.colorbar()
                # plt.show()
                # print('np.max(log_likelihood) ', np.max(log_likelihood))
                # print('np.mim(log_likelihood) ', np.min(log_likelihood))

                # CALCULATE POSTERIOR
                print('Calculating posterior')
                log_prior = np.log(prior_gamma)
                log_posterior = log_prior + log_likelihood
                posterior_gamma_unnorm = np.exp(log_posterior - np.max(log_posterior))

                # UPDATE PRIOR
                print('Updating prior')
                prior_gamma = posterior_gamma_unnorm

                # SAVE CURRENT DATA TO DATA SERVER
                BayesianT1Data.push(
                    {'params': {
                        'datasetName': datasetName,
                        'samplingFreq': samplingFreq,
                        'maxIterations': maxIterations,
                        'freq': freq,
                        'rf_power': rf_power,
                        'laser_power': laser_power,
                        'num_samples' : num_samples,
                        'clock_time' : clock_time,
                        'init_time' : init_time,
                        'laser_lag' : laser_lag,
                        'probe_time' : probe_time,
                        'singlet_decay' : singlet_decay,
                        'bayesian_iterations' : bayesian_iterations,
                        'pi_time' : pi_time
                        },
                        'title': 'Bayesian T1',
                        'xlabel': 'Iteration number',
                        'ylabel': 'Gamma (ms^-1)',
                        'datasets': {
                            'iteration_number': self.iteration_number,
                            'M_plus': self.M_plus,
                            'M_minus': self.M_minus,
                            'GammaPlus': self.gamma_plus_list,
                            'GammaMinus': self.gamma_minus_list,
                            'TauPlus': tau_plus_list,
                            'TauMinus': tau_minus_list,
                            'S000Plus': S_0_0_0_list[0],
                            'S000Minus': S_0_0_0_list[1],
                            'S00tauPlus': S_0_0_tau_list[0],
                            'S00tauMinus': S_0_0_tau_list[1],
                            'Slevel00Plus': S_level_0_0_list[0],
                            'Slevel00Minus': S_level_0_0_list[1],
                            'Slevel0tauPlus': S_level_0_tau_list[0],
                            'Slevel0tauMinus': S_level_0_tau_list[1],
                        }
                    }
                )

            flexSave(datasetName, 'Bayesian T1', 'final')
            print('time taken ', time.time() - start_time)
            print("Experiment finished!")


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



    


