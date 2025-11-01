"""
Experiment to run T1 measuerment on the NV-AFM Spork setup


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


    # def calculate_likelihood(self, M_measured_mean, M_measured_sigma, tau_opt, gamma_grid):
        
    #     print('tau_opt received ', tau_opt)
    #     start_likelihood_calculation = time.time()
    #     M_plus_measured, M_minus_measured = M_measured_mean
    #     sigma_M_plus_measured, sigma_M_minus_measured = M_measured_sigma

    #     M_plus_tilde, M_minus_tilde = self.calculate_M_tildes(gamma_grid, tau_opt)

    #     chi_plus = (M_plus_measured - M_plus_tilde) / ((2**0.5) * (sigma_M_plus_measured))
    #     chi_minus = (M_minus_measured - M_minus_tilde) / (
    #         (2**0.5) * (sigma_M_minus_measured)
    #     )

    #     chi_sq = (chi_plus**2) + (chi_minus**2)
    #     # chi_sq_final = chi_sq - np.min(chi_sq)
    #     # likelihood = np.exp(-(chi_sq_final))
    #     log_likelihood = - chi_sq

    #     # print('Time taken to calculate likelihood ', time.time() - start_likelihood_calculation)

    #     return log_likelihood


    def log_likelihood(self, gamma, tau, M, M_err):  # calculate log likelihood
    
        M_plus, M_minus = M
        M_plus_err, M_minus_err = M_err
        tau_plus, tau_minus = tau

        chi_sq = 0

        for i in range(len(tau_plus)):

            M_plus_tilde, M_minus_tilde = self.calculate_M_tildes(gamma, tau_plus[i], tau_minus[i])

            chi_plus = (M_plus[i] - M_plus_tilde) / ((2**0.5) * M_plus_err[i])
            chi_minus = (M_minus[i] - M_minus_tilde) / ((2**0.5) * M_minus_err[i])
            chi_sum = (chi_plus**2) + (chi_minus**2)

            chi_sq = chi_sq - chi_sum

        return chi_sq


    def log_prior(self, gamma):  # create flat prior for the gammas

        gamma_plus, gamma_minus = gamma

        if gamma_bounds[0] < gamma_plus < gamma_bounds[1] and gamma_bounds[0] < gamma_minus < gamma_bounds[1]:
            return 0.0

        return -np.inf


    def log_posterior(self, gamma, tau, M, M_err):

        lp = self.log_prior(gamma)

        if not np.isfinite(lp):
            return -np.inf

        return lp + self.log_likelihood(gamma, tau, M, M_err)


    def calculate_M_tildes(self, gamma_grid, tau_plus, tau_minus):

        gamma_plus, gamma_minus = gamma_grid
        # pp, pm, mm, mp, Mtplus, Mtminus = BF.give_sympy_functions()
        
        M_plus_tilde = Symbols.Mtplus(gamma_plus, gamma_minus, tau_plus)
        M_minus_tilde = Symbols.Mtminus(gamma_plus, gamma_minus, tau_minus)

        return [M_plus_tilde, M_minus_tilde]


    def BayesianT1(
        self,
        datasetName: str,
        num_samples: int,
        rf_power: float,
        laser_power: float,
        Lower_Freq: float,
        Lower_Pi: float,
        Higher_Freq: float,
        Higher_Pi: float,
        bayesian_iterations: int,
    ):

        # connect to the instrument server
        # connect to the data server and create a data set, or connect to an
        # existing one with the same name if it was created earlier.
        with InstrumentGateway() as gw, DataSource(datasetName) as BayesianT1Data:

            # Constants
            freq_list = [Higher_Freq, Lower_Freq]
            pi_time_list = [Higher_Pi, Lower_Pi]
            start_time = time.time()
            readout_time = 500

            # PLOTTING VARIABLES
            self.iteration_number = np.zeros(0)
            self.gamma_plus_list = np.zeros(0)
            self.gamma_minus_list = np.zeros(0)
            self.gamma_plus_err_list = np.zeros(0)
            self.gamma_minus_err_list = np.zeros(0)
            self.elapsed_time = np.zeros(0)
            self.M_plus = np.zeros(0)
            self.M_minus = np.zeros(0)
            self.M_plus_err = np.zeros(0)
            self.M_minus_err = np.zeros(0)
            self.tau_plus_list = np.zeros(0)
            self.tau_minus_list = np.zeros(0)
            self.elapsed_time = np.zeros(0)
            self.S_0_0_0_list = []
            self.S_0_0_tau_list = []
            self.S_level_0_0_list = []
            self.S_level_0_tau_list = []

            # SETUP BAYESIAN VARIABLES
            tau_lower = 0.01 # in ms (10us)
            tau_upper = 10 # in ms
            n_tau = 1000
            repetitions = num_samples 
            control_tau_list = np.geomspace(tau_lower, tau_upper, bayesian_iterations)

            # relaxometry delay tau grid
            tau_plus_arr = np.geomspace(tau_lower, tau_upper, n_tau)
            tau_minus_arr = np.geomspace(tau_lower, tau_upper, n_tau)
            tau_grid = np.meshgrid(tau_plus_arr, tau_minus_arr)

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
            n_dim = 2   # gamma_plus, gamma_minus
            n_walkers = 50
            n_steps = 1000
            burn_in = 500

            # STARTING POSITION FOR WALKERS AROUND INITIAL GUESS
            init_guess = [1,1]
            # init_guess = [np.mean(gamma_bounds), np.mean(gamma_bounds)]
            p0 = [init_guess + (1e-4 * np.random.randn(n_dim)) for i in range(n_walkers)]

            # GAMMA SAMPLES LIST
            gamma_samples_list = []

            # # CREATE GAMMA CLOUD DIRECTORY
            dir_title = 'C://Users/awschlab/Desktop/data/gamma_clouds/CONTROL_10ms_GammaBounds=' + str(gamma_bounds) + '_InitialGuess=' + str(init_guess) + '_Freqs=' + str(freq_list) + '_R=' + str(num_samples) + 'NoFeedback'
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
                    tau_opt = 0.005, 0.005
                else:
                    tau_opt = self.nob_calculate_tau_opt(tau_grid, repetitions, self.gamma_plus_list[num - 1], self.gamma_minus_list[num - 1], T_overhead)
                tau_plus, tau_minus = tau_opt
                self.tau_plus_list = np.append(self.tau_plus_list, tau_plus)
                self.tau_minus_list = np.append(self.tau_minus_list, tau_minus)
                print('tau_plus_opt (ms) ', tau_plus)
                print('tau_minus_opt (ms) ', tau_minus)
                ##############################################################################


                #############################################################################
                # #### CONTROL TAUS ####
                # # Use geometrically spaced tau array
                # print('Turned off tau optimization, using geometrically spaced taus')
                # tau_plus, tau_minus = control_tau_list[num], control_tau_list[num]
                # self.tau_plus_list = np.append(self.tau_plus_list, tau_plus)
                # self.tau_minus_list = np.append(self.tau_minus_list, tau_minus)
                # print('tau_plus_opt (ms) ', tau_plus)
                # print('tau_minus_opt (ms) ', tau_minus)
                ##############################################################################


                # CREATE PULSE SEQUENCES FOR OPTIMIZED TAUS
                print('Creating pulse sequences')
                seqs = []

                ##############################################################################
                #### MW T1 PULSE SEQUENCE ####
                seqs.append(Pulses(gw).BAYESIAN_T1(int(1e6 * tau_plus), pi_time_list[0]))
                seqs.append(Pulses(gw).BAYESIAN_T1(int(1e6 * tau_minus), pi_time_list[1]))
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

                    for t in range(2):  # measure NV with delay = tau_plus and then delay = tau_minus

                        # SRS ACTIONS
                        gw.sg.set_rf_amplitude(rf_power)    # set ouput power
                        gw.sg.set_frequency(1e9*freq_list[t])
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
                        gw.swabian.runSequenceInfinitely(seqs[t])

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
                    self.M_plus = np.append(self.M_plus, M[0])
                    self.M_minus = np.append(self.M_minus, M[1])
                    self.M_plus_err = np.append(self.M_plus_err, M_err[0])
                    self.M_minus_err = np.append(self.M_minus_err, M_err[1])

                    print('self.gamma_plus_list ', self.gamma_plus_list)
                    print('self.M_plus ', self.M_plus)
                    print('self.M_plus_err ', self.M_plus_err)
                    print('self.iteration_number ', self.iteration_number)

                    # SLICE DATA TO RUN SAMPLER WITH FIRST N SAMPLES
                    print('Slicing data')
                    tau_plus_slice = self.tau_plus_list[:num+1]     # x data
                    tau_minus_slice = self.tau_minus_list[:num+1]   # x data
                    M_plus_slice = self.M_plus[:num+1]              # y1 data
                    M_minus_slice = self.M_minus[:num+1]            # y2 data
                    M_plus_err_slice = self.M_plus_err[:num+1]      # y1 data
                    M_minus_err_slice = self.M_minus_err[:num+1]    # y2 data

                    # RUN MCMC SAMPLER
                    print('Running sampler')
                    sampler = emcee.EnsembleSampler(
                        n_walkers, n_dim, self.log_posterior, args = ((tau_plus_slice, tau_minus_slice), (M_plus_slice, M_minus_slice), (M_plus_err_slice, M_minus_err_slice)))
                    sampler.run_mcmc(p0, n_steps)#, progress=True)

                    # STORE GAMMA SAMPLES AFTER BURN-IN
                    print('Storing gamma_samples')
                    gamma_samples = sampler.get_chain(discard=burn_in, flat=True)
                    gamma_samples_list.append(gamma_samples)

                    gamma_plus_samples = gamma_samples[:, 0]
                    gamma_minus_samples = gamma_samples[:, 1]

                    # save gamma scatter plot
                    plt.figure(figsize=(6, 6))
                    plt.scatter(gamma_plus_samples, gamma_minus_samples, color="blue", s=0.1, alpha=0.5)
                    plt.xlabel("gamma_plus")
                    plt.ylabel("gamma_minus")
                    plt.xlim(gamma_bounds)
                    plt.ylim(gamma_bounds)
                    # plt.xlim(1, 10)
                    # plt.ylim(1, 10)
                    title = 'Gammas after ' + str(num) + ' iterations'
                    plt.suptitle(title)
                    filetitle = dir_title + '/' + str(datetime.now().strftime("%Y_%m_%d_%H_%M_%S")) + '.png' 
                    
                    plt.savefig(filetitle)

                    mean_gamma_plus = np.mean(gamma_plus_samples)
                    mean_gamma_minus = np.mean(gamma_minus_samples)
                    err_gamma_plus = np.std(gamma_plus_samples) / np.sqrt(len(gamma_plus_samples))
                    err_gamma_minus = np.std(gamma_minus_samples) / np.sqrt(len(gamma_minus_samples))

                    print('Updating gamma lists')
                    self.gamma_plus_list = np.append(self.gamma_plus_list, mean_gamma_plus)
                    self.gamma_minus_list = np.append(self.gamma_minus_list, mean_gamma_minus)
                    self.gamma_plus_err_list = np.append(self.gamma_plus_err_list, err_gamma_plus)
                    self.gamma_minus_err_list = np.append(self.gamma_minus_err_list, err_gamma_minus)
                    print('gamma_plus ', mean_gamma_plus)
                    print('gamma_minus ', mean_gamma_minus)

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
                            'Higher_Freq': Higher_Freq,
                            'Higher_Pi': Higher_Pi,
                            'bayesian_iterations' : bayesian_iterations,
                            },

                            'title': 'Bayesian T1',
                            'xlabel': 'Iteration number',
                            # 'ylabel': 'Gamma (ms^-1)',
                            'datasets': {
                                'iteration_number': self.iteration_number,
                                'M_plus': self.M_plus,
                                'M_minus': self.M_minus,
                                'M_plus_err': self.M_plus_err,
                                'M_minus_err': self.M_minus_err,
                                'GammaPlus': self.gamma_plus_list,
                                'GammaMinus': self.gamma_minus_list,
                                'GammaPlusErr': self.gamma_plus_err_list,
                                'GammaMinusErr': self.gamma_minus_err_list,
                                'TauPlus': self.tau_plus_list,
                                'TauMinus': self.tau_minus_list,
                                # 'S000Plus': np.array(self.S_0_0_0_list)[:,0,:],
                                # 'S000Minus': np.array(self.S_0_0_0_list)[:,1,:],
                                # 'S00tauPlus': np.array(self.S_0_0_tau_list)[:,0,:],
                                # 'S00tauMinus': np.array(self.S_0_0_tau_list)[:,1,:],
                                # 'Slevel00Plus': np.array(self.S_level_0_0_list)[:,0,:],
                                # 'Slevel00Minus': np.array(self.S_level_0_0_list)[:,1,:],
                                # 'Slevel0tauPlus': np.array(self.S_level_0_tau_list)[:,0,:],
                                # 'Slevel0tauMinus': np.array(self.S_level_0_tau_list)[:,1,:],
                                'TotalTime': self.elapsed_time
                            }
                        }
                    )

                    notes = ''
                    flexSave(datasetName, notes, 'Bayesian T1')

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



    


