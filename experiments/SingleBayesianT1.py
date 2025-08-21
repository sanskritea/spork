"""
Experiment to run T1 measuerment on the NV-AFM Spork setup
Only runs for the lower peak (where NV resonates with magnons)

"""

import time
from itertools import count

import numpy as np
import matplotlib.pyplot as plt
import random
import scipy as scp
# from scipy.interpolate import interp1d
import sympy
from sympy import Symbol
from sympy import lambdify
import emcee
from datetime import datetime
from os import mkdir

from nspyre import DataSource
from nspyre import InstrumentGateway

from rpyc.utils.classic import obtain
from guis.guiElements_general import flexSave

from experiments.NewPulses import Pulses
from experiments.NewportSpatialFeedback import SpatialFeedback
from drivers.ni.nidaq_final import NIDAQ
# from experiments.BayesianFunctions import BF

gamma_bounds = [0.01,9]



class Bayesian_T1_Meas:

    # def trap_normalize_2D_pdf(self, pdf, gamma_grid):

    #     gamma_plus, gamma_minus = gamma_grid[0][0], np.transpose(gamma_grid[1])[0]
    #     norm = np.trapz(np.trapz(pdf, gamma_minus), gamma_plus)
    #     # print('norm ', norm)

    #     return norm


    # def calculate_conf_intervals(self, prior_gamma, gamma_grid, delta_gamma):

    #     gamma_plus = gamma_grid[0][0]
    #     gamma_minus = np.transpose(gamma_grid[1])[0]

    #     # print('Normalizing pdfs')
    #     gamma_plus_distr = np.sum(prior_gamma, 0) * delta_gamma
    #     norm_gamma_plus_distr = gamma_plus_distr / (np.sum(gamma_plus_distr) * delta_gamma)

    #     gamma_minus_distr = np.sum(prior_gamma, 1) * delta_gamma
    #     norm_gamma_minus_distr = gamma_minus_distr / (np.sum(gamma_minus_distr) * delta_gamma)

    #     # print('Calculating cdfs')
    #     gamma_plus_cdf = np.cumsum(norm_gamma_plus_distr) * delta_gamma
    #     gamma_minus_cdf = np.cumsum(norm_gamma_minus_distr) * delta_gamma

    #     # print('Interpolating')
    #     gamma_plus_cdf_interp = interp1d(gamma_plus_cdf, gamma_plus) #, fill_value = 'extrapolate')
    #     gamma_minus_cdf_interp = interp1d(gamma_minus_cdf, gamma_minus) #, fill_value = 'extrapolate')

    #     gamma_plus_vals = gamma_plus_cdf_interp(np.array([0.5]))
    #     gamma_minus_vals = gamma_minus_cdf_interp(np.array([0.5]))

    #     return gamma_plus_vals, gamma_minus_vals


    def nob_calculate_tau_opt(self, tau, repetitions, gamma, T_overhead):
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

        num = 1 / (2 * Symbols.del_two_M_tilde_by_del_gamma_sq)

        den = gamma

        T = 2 * repetitions * (tau) + T_overhead

        # actual cost function to compare to paper
        cost_function = (num / den) * (t ** 0.5)

        # clean up cost function by changing all nans to an arbitrarily large value
        flag = 0
        for nni in range(len(tau)):
            if np.isnan(cost_function[nni]) or np.isinf(cost_function[nni]):
                    # print('AHA!')
                cost_function[nni] = 1e100
                flag += 1

        t = tau.flatten()
        min_cost_idx = np.argmin(cost_function.flatten())

        # print('Time taken to optimize tau ', time.time() - start_tau_optimization)

        return np.array(t[min_cost_idx])


    def log_likelihood(self, gamma, tau, M, M_err):  # calculate log likelihood
    
        M_meas = M
        M_err = M_err

        chi_sq = 0

        for i in range(len(tau)):

            M_tilde = self.calculate_M_tildes(gamma, tau[i])

            chi = (M_meas[i] - M_tilde) / ((2**0.5) * M_err[i])
            chi_sum = (chi**2)
            chi_sq = chi_sq - chi_sum

        return chi_sq


    def log_prior(self, gamma):  # create flat prior for the gammas


        if gamma_bounds[0] < gamma < gamma_bounds[1]:
            return 0.0

        return -np.inf


    def log_posterior(self, gamma, tau, M, M_err):

        lp = self.log_prior(gamma)

        if not np.isfinite(lp):
            return -np.inf

        return lp + self.log_likelihood(gamma, tau, M, M_err)


    def calculate_M_tildes(self, gamma, tau):
        
        M_tilde = Symbols.Mt(gamma, tau)

        return [M_tilde]


    def BayesianT1(
        self,
        datasetName: str,
        num_samples: int,
        rf_power: float,
        laser_power: float,
        Lower_Freq: float,
        Lower_Pi: float,
        bayesian_iterations: int,
    ):

        # connect to the instrument server
        # connect to the data server and create a data set, or connect to an
        # existing one with the same name if it was created earlier.
        readout_time = probe_time
        with InstrumentGateway() as gw, DataSource(datasetName) as BayesianT1Data:

            # Constants
            freq = [Lower_Freq]
            pi_time_list = [Lower_Pi]
            start_time = time.time()

            # PLOTTING VARIABLES
            self.iteration_number = np.zeros(0)
            self.gamma_list = np.zeros(0)
            self.gamma_err_list = np.zeros(0)
            self.elapsed_time = np.zeros(0)
            self.M = np.zeros(0)
            self.M_err = np.zeros(0)
            self.tau_list = np.zeros(0)
            self.elapsed_time = np.zeros(0)
            self.S_0_0_0_list = []
            self.S_0_0_tau_list = []
            self.S_level_0_0_list = []
            self.S_level_0_tau_list = []

            # SETUP BAYESIAN VARIABLES
            tau_lower = 0.001 # in ms (10us)
            tau_upper = 1 # in ms
            n_tau = 1000
            repetitions = num_samples 
            control_tau_list = np.geomspace(tau_lower, tau_upper, bayesian_iterations)

            # relaxometry delay tau grid
            tau_arr = np.geomspace(tau_lower, tau_upper, n_tau)

            # Feedback parameters
            feedback_trigger_rate = int(20e3)
            feedback_time_per_point = 0.05
            feedback_num_samples = int(feedback_trigger_rate * feedback_time_per_point)
            objective_x_init_position = 4.6196
            objective_y_init_position = 11.3896
            objective_z_init_position = -0.940763
            feedback_timer = time.time()
            feedback_counter = 0

            # SAMPLER PARAMETERS
            n_dim = 1   # gamma_plus, gamma_minus
            n_walkers = 50
            n_steps = 1000
            burn_in = 500

            # STARTING POSITION FOR WALKERS AROUND INITIAL GUESS
            init_guess = [1]
            # init_guess = [np.mean(gamma_bounds), np.mean(gamma_bounds)]
            p0 = [init_guess + (1e-4 * np.random.randn(n_dim)) for i in range(n_walkers)]

            # GAMMA SAMPLES LIST
            gamma_samples_list = []

            # # CREATE GAMMA CLOUD DIRECTORY
            dir_title = 'C://Users/awschlab/Desktop/data/gamma_clouds/YIG300nm_TAU_GUESS_70G_GammaBounds=' + str(gamma_bounds) + '_InitialGuess=' + str(init_guess) + '_Freqs=' + str(freq_list) + '_R=' + str(num_samples) + 'NoFeedback'
            mkdir(dir_title)


            for num in range(bayesian_iterations):

                # MAIN EXPERIMENT
                print('Bayesian iteration number ', num)
                self.iteration_number = np.append(self.iteration_number, num)

                ##############################################################################
                #### STARTING WITH GAMMA GUESS ####
                # # USE INITIAL WALKER POSITION AS GAMMA GUESS
                # if num == 0:
                #     self.gamma_plus_list = np.append(self.gamma_plus_list, init_guess[0])
                #     self.gamma_minus_list = np.append(self.gamma_minus_list, init_guess[1])
                #     self.gamma_plus_err_list = np.append(self.gamma_plus_err_list, init_guess[0] / 10)
                #     self.gamma_minus_err_list = np.append(self.gamma_minus_err_list, init_guess[1] / 10)
                # gamma_plus = self.gamma_plus_list[num]
                # gamma_minus = self.gamma_minus_list[num]

                # # OPTIMIZE TAU WITH THIS GAMMA
                # print('Optimizing taus for pulse sequences')
                # T_overhead = 0
                # tau_opt = self.nob_calculate_tau_opt(tau_grid, repetitions, gamma_plus, gamma_minus, T_overhead)
                # tau_plus, tau_minus = tau_opt
                # self.tau_plus_list = np.append(self.tau_plus_list, tau_plus)
                # self.tau_minus_list = np.append(self.tau_minus_list, tau_minus)
                # print('tau_plus_opt (ms) ', tau_plus)
                # print('tau_minus_opt (ms) ', tau_minus)
                ##############################################################################


                ##############################################################################
                #### STARTING WITH TAU GUESS/CALCULATION ####
                # GET OPTIMIZED TAUS FROM NOB TO USE IN EXPERIMENT
                print('Optimizing taus for pulse sequences')
                T_overhead = 0
                tau_opt = []
                if num == 0:
                    # tau_opt = tau_lower, tau_lower
                    tau_opt = 0.005
                else:
                    tau_opt = self.nob_calculate_tau_opt(tau, repetitions, self.gamma_list[num - 1], T_overhead)
                tau = tau_opt
                self.tau_list = np.append(self.tau_list, tau)
                print('tau_opt (ms) ', tau)
                ##############################################################################


                ##############################################################################
                # #### CONTROL TAUS ####
                # # Use geometrically spaced tau array
                # print('Turned off tau optimization, using geometrically spaced taus')
                # tau = control_tau_list[num], control_tau_list[num]
                # self.tau_list = np.append(self.tau_plus_list, tau)
                # print('tau_opt (ms) ', tau)
                ##############################################################################


                # CREATE PULSE SEQUENCES FOR OPTIMIZED TAUS
                print('Creating pulse sequences')
                seqs = None

                ##############################################################################
                #### MW T1 PULSE SEQUENCE ####
                seqs = Pulses(gw).BAYESIAN_T1(int(1e6 * tau), pi_time_list[0])
                # seqs.append(Pulses(gw).BAYESIAN_T1(int(1e6 * tau_minus), pi_time_list[1]))
                ##############################################################################


                # MEASURE COUNTS
                print('Measuring counts')
                with NIDAQ() as mynidaq:

                    # # SPATIAL FEEDBACK EVERY 5 minutes
                    # if ((time.time() - feedback_timer) > 12000):
                    #     feedback_counter = feedback_counter + 1
                    #     print('Feedback')
                    #     begin_feedback = time.time()
                    #     trigger_rate = int(20e3)
                    #     time_per_point = 0.003
                    #     objective_x_init_position, objective_y_init_position, current_max_counts = SpatialFeedback.Feedback(objective_x_init_position, objective_y_init_position, objective_z_init_position)
                    #     print('Current objective_x_position ', objective_x_init_position)
                    #     print('Current objective_y_position ', objective_y_init_position)
                    #     print('Returned max counts ', current_max_counts)
                        
                    #     gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(int(trigger_rate)))
                    #     num_samples = int(trigger_rate * time_per_point)
                    #     post_feedback_counts = np.mean(obtain(mynidaq.internal_read_task(int(trigger_rate), num_samples))) / (1 / trigger_rate)
                    #     print('Post feedback counts ', post_feedback_counts)
                    #     gw.swabian.reset()
                    #     feedback_duration = time.time() - begin_feedback
                    #     print('Feedback duration: ', feedback_duration)
                    #     print('Feedback counter ', feedback_counter)
                    #     feedback_timer = time.time()

                    # SET LASER POWER
                    mynidaq.laser_power_atten(laser_power)

                    # temporary storage lists
                    M = []
                    M_err = []

                    S_0_0_0 = []
                    S_0_0_tau = []
                    S_level_0_0 = []
                    S_level_0_tau = []

                    # SRS ACTIONS
                    gw.sg.set_rf_amplitude(rf_power)    # set ouput power
                    gw.sg.set_frequency(1e9*freq)
                    gw.sg.set_mod_state(False) # make QAM/external after IQ is ready
                    gw.sg.set_rf_state("1")

                    # check counts before measuring
                    gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(feedback_trigger_rate))
                    current_counts = np.mean(obtain(mynidaq.internal_read_task(feedback_trigger_rate, int(feedback_trigger_rate * 0.01)))) / (1 / feedback_trigger_rate)
                    print('Current counts ', current_counts)
                    gw.swabian.reset()

                    ######################################################################
                    # ##### OPTICAL T1 NORMALIZATION
                    # gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(feedback_trigger_rate))
                    # norm = np.mean(obtain(mynidaq.internal_read_task(feedback_trigger_rate, int(feedback_trigger_rate * 0.01)))) / (1 / feedback_trigger_rate)
                    # print('Normalization ', norm)
                    # gw.swabian.reset()
                    ######################################################################

                    # START READ TASK
                    #### FOR BAYESIAN T1 ####
                    mynidaq.start_external_read_task(20e6, ((8 * num_samples) + 1))

                    # START PULSESTREAMER
                    gw.swabian.runSequenceInfinitely(seqs)

                    ######################################################################
                    # #### OPTICAL BAYESIAN T1 ####
                    # raw_counts = (obtain(mynidaq.external_read_task(20e6, ((2 * num_samples) + 1))))
                    # counts = 1e9 * np.mean(raw_counts[0::2]) / readout_time
                    # print('counts ', counts)
                    # norm_counts = (1e9 * np.mean(raw_counts[0::2]) / readout_time)
                    # norm_err = np.std((1e9 * raw_counts[0::2] / readout_time)) / (num_samples ** 0.5)
                    # M.append(norm_counts)
                    # M_err.append(norm_err)
                    # print('norm_counts ', norm_counts)
                    # print('norm_err ', norm_err)
                    ######################################################################


                    ######################################################################
                    # #### MW BAYESIAN T1 ####
                    # COUNT WITH EXTERNAL TRIGGER
                    counts = obtain(mynidaq.external_read_task(20e6, ((8 * num_samples) + 1)))

                    # SEPARATE COUNTS S_(initial level)_(final level)_delay
                    S_0_0_0.append(1e9 * counts[0::8] / readout_time)
                    S_0_0_tau.append(1e9 * counts[2::8] / readout_time) 
                    S_level_0_0.append(1e9 * counts[4::8] / readout_time)
                    S_level_0_tau.append(1e9 * counts[6::8] / readout_time)

                    S_0_0_0_mean = np.mean(S_0_0_0[t])
                    S_0_0_tau_mean = np.mean(S_0_0_tau[t]) 
                    S_level_0_0_mean = np.mean(S_level_0_0[t]) 
                    S_level_0_tau_mean = np.mean(S_level_0_tau[t]) 

                    S_0_0_0_err = np.std(S_0_0_0[t]) / (num_samples ** 0.5)
                    S_0_0_tau_err = np.std(S_0_0_tau[t]) / (num_samples ** 0.5)
                    S_level_0_0_err = np.std(S_level_0_0[t]) / (num_samples ** 0.5)
                    S_level_0_tau_err = np.std(S_level_0_tau[t]) / (num_samples ** 0.5)

                    # FINAL MEASURE AND MOMENTS (Appendix E)
                    top = S_0_0_tau_mean - S_level_0_tau_mean
                    bottom = S_0_0_0_mean - S_level_0_0_mean
                    top_err = (S_0_0_tau_err ** 2 + S_level_0_tau_err ** 2) ** 0.5
                    bottom_err = (S_0_0_0_err ** 2 + S_level_0_0_err ** 2) ** 0.5
                    if bottom == 0:
                        print('bottom ', bottom)
                        bottom = 1e-10
                    if bottom_err == 0:
                        print('bottom_err ', bottom_err)
                        bottom_err = 1e-10
                    ratio = bottom/(bottom_err * bottom_err)
                    relative = bottom_err / bottom
                    z0 = 0.25*ratio*( ( 1 + 8*relative*relative )**0.5 - 1 )
                    z02 = z0*z0 
                    eA02 = bottom_err * bottom_err    
                    dL2dz2 = 2*ratio/(z02*z0) + 2/z02 - 3/(eA02*z02*z02) 
                    ez0 = 1/ (-dL2dz2)**0.5
                    inv_bottom = z0
                    inv_bottom_err = ez0

                    z = top * inv_bottom
                    z_err = ((top * inv_bottom_err)**2 + (inv_bottom_err * top_err) ** 2) ** 0.5

                    M.append(z)
                    M_err.append(z_err)
                    print('M_mean ', z)
                    print('M_err ', z_err)
                    ######################################################################
                    
                    gw.swabian.reset()

                    # UPDATE LISTS
                    print('Updating all lists')
                    self.S_0_0_0_list.append(S_0_0_0)
                    self.S_0_0_tau_list.append(S_0_0_tau)
                    self.S_level_0_0_list.append(S_level_0_0)
                    self.S_level_0_tau_list.append(S_level_0_tau)
                    self.M = np.append(self.M, M[0])
                    self.M_err = np.append(self.M_err, M_err[0])

                    print('self.gamma_list ', self.gamma_list)
                    print('self.M ', self.M)
                    print('self.M_err ', self.M_err)
                    print('self.iteration_number ', self.iteration_number)

                    


                    # SLICE DATA TO RUN SAMPLER WITH FIRST N SAMPLES
                    print('Slicing data')
                    tau_slice = self.tau_list[:num+1]     # x data
                    M_slice = self.M[:num+1]              # y1 data
                    M_err_slice = self.M_err[:num+1]      # y1 data

                    # RUN MCMC SAMPLER
                    print('Running sampler')
                    sampler = emcee.EnsembleSampler(
                        n_walkers, n_dim, self.log_posterior, args = ((tau_slice), (M_slice), (M_err_slice)))
                    sampler.run_mcmc(p0, n_steps)#, progress=True)

                    # STORE GAMMA SAMPLES AFTER BURN-IN
                    print('Storing gamma_samples')
                    gamma_samples = sampler.get_chain(discard=burn_in, flat=True)
                    gamma_samples_list.append(gamma_samples)

                    gamma_samples = gamma_samples[:, 0]

                    # save gamma scatter plot
                    plt.figure(figsize=(6, 6))
                    plt.hist(gamma_samples)
                    plt.xlabel("gamma")
                    plt.ylabel("Posterior")
                    plt.xlim(gamma_bounds)
                    # plt.ylim(gamma_bounds)
                    # plt.xlim(1, 10)
                    # plt.ylim(1, 10)
                    title = 'Gammas after ' + str(num) + ' iterations'
                    plt.suptitle(title)
                    filetitle = dir_title + '/' + str(datetime.now().strftime("%Y_%m_%d_%H_%M_%S")) + '.png' 
                    
                    plt.savefig(filetitle)

                    mean_gamma = np.mean(gamma_samples)
                    err_gamma = np.std(gamma_samples) / np.sqrt(len(gamma_samples))

                    print('Updating gamma lists')
                    self.gamma_list = np.append(self.gamma_list, mean_gamma)
                    self.gamma_err_list = np.append(self.gamma_err_list, err_gamma)
                    print('gamma ', mean_gamma)

                    
                    self.elapsed_time = np.append(self.elapsed_time, time.time() - start_time)

                    # print('Push to GUI')
                    # SAVE CURRENT DATA TO DATA SERVER
                    BayesianT1Data.push(
                        {'params': {
                            'datasetName': datasetName,
                            'rf_power': rf_power,
                            'laser_power': laser_power,
                            'num_samples' : num_samples,
                            'Lower_Freq': Lower_Freq,
                            'Lower_Pi': Lower_Pi,
                            'bayesian_iterations' : bayesian_iterations,
                            },

                            'title': 'Bayesian T1',
                            'xlabel': 'Iteration number',
                            # 'ylabel': 'Gamma (ms^-1)',
                            'datasets': {
                                'iteration_number': self.iteration_number,
                                'M': self.M,
                                'M_err': self.M_err,
                                'Gamma': self.gamma_list,
                                'GammaErr': self.gamma_err_list,
                                'Tau': self.tau_list,
                                'S000': np.array(self.S_0_0_0_list)[:,0,:],
                                'S00tau': np.array(self.S_0_0_tau_list)[:,0,:],
                                'Slevel00': np.array(self.S_level_0_0_list)[:,0,:],
                                'Slevel0tau': np.array(self.S_level_0_tau_list)[:,0,:],
                                'TotalTime': self.elapsed_time
                            }
                        }
                    )

                    notes = ''
                    flexSave(datasetName, notes, 'Bayesian T1')

            print('time taken ', time.time() - start_time)
            print("Experiment finished!")


class Symbols:

    # tau_plus = Symbol("tau_-")
    # tau_minus = Symbol("tau_+")
    # gamma_plus = Symbol("Gamma_+")
    # gamma_minus = Symbol("Gamma_-")

    # G = sympy.sqrt(gamma_plus**2 + gamma_minus**2 - gamma_plus * gamma_minus)
    # beta_plus = gamma_plus + gamma_minus + G
    # beta_minus = gamma_plus + gamma_minus - G

    # M_tilde_plus = (G + gamma_plus) * sympy.exp(-tau_plus * beta_plus) + (
    #     G - gamma_plus
    # ) * sympy.exp(-tau_plus * beta_minus)
    # M_tilde_minus = (G + gamma_minus) * sympy.exp(-tau_minus * beta_plus) + (
    #     G - gamma_minus
    # ) * sympy.exp(-tau_minus * beta_minus)
    # M_tilde_plus = M_tilde_plus / (2 * G)
    # M_tilde_minus = M_tilde_minus / (2 * G)

    # pp = lambdify(
    #     [gamma_plus, gamma_minus, tau_plus],
    #     sympy.diff(M_tilde_plus, gamma_plus).simplify(),
    # )
    # pm = lambdify(
    #     [gamma_plus, gamma_minus, tau_plus],
    #     sympy.diff(M_tilde_plus, gamma_minus).simplify(),
    # )
    # mm = lambdify(
    #     [gamma_plus, gamma_minus, tau_minus],
    #     sympy.diff(M_tilde_minus, gamma_minus).simplify(),
    # )
    # mp = lambdify(
    #     [gamma_plus, gamma_minus, tau_minus],
    #     sympy.diff(M_tilde_minus, gamma_plus).simplify(),
    # )
    # Mtplus = lambdify([gamma_plus, gamma_minus, tau_plus], M_tilde_plus.simplify())
    # Mtminus = lambdify([gamma_plus, gamma_minus, tau_minus], M_tilde_minus.simplify())

    # for single peak experiment
    tau = Symbol("tau")
    gamma = Symbol("Gamma")

    G = sympy.sqrt(gamma**2)
    beta = (2 * gamma) + G

    M_tilde = (G + gamma) * sympy.exp(-tau * beta)

    M = lambdify([gamma, tau], M_tilde.simplify())
    del_two_M_tilde_by_del_gamma_sq = lambdify(
        [gamma, tau],
        sympy.diff(M_tilde, gamma, 2).simplify(),
    )



    


