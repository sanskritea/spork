"""
This is an application to adaptively measure NV T1 using a two-measurement protocol

Copyright (c) February 2024, Sanskriti Chitransh
All rights reserved.

This work is licensed under the terms of the 3-Clause BSD license.
For a copy, see <https://opensource.org/licenses/BSD-3-Clause>.

"""

import time
from itertools import count

import numpy as np

from nspyre import DataSource
from nspyre import InstrumentGateway

from rpyc.utils.classic import obtain

from experiments.NewPulses import Pulses
from drivers.ni.nidaq_final import NIDAQ

# from PulsePatterns import Pulses


class Adaptive_Two_Meas_T1_Measurement:

    def BayesianT1(
        self,
        datasetName: str,
        samplingFreq: float,
        maxIterations: int,
        mw_freq: float,
        rfPower: float,
        laser_power: float,
        # debug=False):
        num_samples: int,
        clock_time: int,
        init_time: int,
        pi_time: int,
        repetitions: int,
        gamma_lower: float,
        gamma_upper: float,
    ):

        # connect to the instrument server
        # connect to the data server and create a data set, or connect to an
        # existing one with the same name if it was created earlier.
        with InstrumentGateway() as gw, DataSource(datasetName) as BayesianT1Data:

            # SRS actions
            gw.sg.set_mod_state(False)
            gw.sg.set_rf_amplitude(rfPower)  # set ouput power
            gw.sg.set_frequency(freq)

            # setup a relevant iterator depending on maxIterations
            if maxIterations < 0:
                iters = count()  # infinite iterator
            else:
                iters = range(maxIterations)

            # MAIN EXPERIMENT LOOP
            with NIDAQ() as mynidaq:

                # set laser power
                mynidaq.laser_power_atten(laser_power)

                # SETUP BAYESIAN VAIRABLES

                # decay rates grid
                n_gamma = num_samples
                gamma_plus_range = np.linspace(gamma_lower, gamma_upper, n_gamma)
                gamma_minus_range = np.linspace(gamma_lower, gamma_upper, n_gamma)

                # decay rates distribution
                gamma_plus_distr = np.ones(n_gamma) / n_gamma
                gamma_minus_distr = np.ones(n_gamma) / n_gamma

                # Total experimental times
                # Questions: should this be log spaced if the gammas are linearly spaced?
                n_tau = 1000
                tau_plus_arr = np.linspace(3, 5500, n_tau)  # in microseconds
                tau_minus_arr = np.linspace(3, 5500, n_tau)  # in microseconds

                T = np.zeros(n_tau, n_tau)  # total experimental time
                T0 = 0  # overhead delays, assume 0 for now
                for i in range(n_tau):
                    for j in range(n_tau):
                        T[i][j] = (
                            2 * repetitions * (tau_plus_list[i] + tau_minus_list[i])
                            + T0
                        )

                prior_plus = gamma_plus_distr  # starting with flat prior in gamma_plus
                prior_minus = (
                    gamma_minus_distr  # starting with flat prior in gamma_minus
                )

                for i in iters:

                    # calculate cost function and minimize it to find wait times
                    gamma_plus_mean = np.mean(prior_plus)
                    gamma_plus_sigma = np.std(prior_plus)
                    gamma_minus_mean = np.mean(prior_minus)
                    gamma_minus_sigma = np.std(prior_minus)

                    tau_plus_opt, tau_minus_opt = tau_opt_from_cost_function(
                        tau_plus_arr,
                        tau_minus_arr,
                        gamma_plus_mean,
                        gamma_minus_mean,
                        gamma_plus_sigma,
                        gamma_minus_sigma,
                        T0,
                    )
                    tau_vals = [tau_plus_opt, tau_minus_opt]

                    # use wait times in measurement
                    M_measured = np.zeros((2, num_samples))

                    for t in range(
                        2
                    ):  # measure NV with delay = tau_plus and then delay = tau_minus

                        # INCORPORATE REPETITIONS AFTER FIGURING OUT THE LIKELIHOOD FUNCTION
                        # Question to self: is R same as num_samples?
                        # for r in repetitions: # number of times the T1 protocol should repeat for a given tau_optimized (tau_plu or tau_minus)

                        # START READ TASK
                        mynidaq.start_read_task(4 * num_samples + 1)

                        # START PULSESTREAMER
                        gw.swabian.runSequenceInfinitely(
                            Pulses(gw).Two_Meas_MW_T1(
                                init_time,
                                clock_time,
                                singlet_decay,
                                tau_vals[t],
                                pi_time,
                            )
                        )

                        # COUNT WITH EXTERNAL TRIGGER
                        counts = obtain(mynidaq.read_samples(num_samples * 4 + 1))

                        # SEPARATE COUNTS
                        # S_(initial level)_(final level)_delay
                        S_0_0_0 = counts[0::4]
                        S_0_0_tau = counts[1::4]
                        S_level_0_0 = counts[2::4]
                        S_level_0_tau = counts[3::4]

                        # Final measure
                        M_measured[t] = (S_0_0_tau - S_level_0_tau) / (
                            S_0_0_0 - S_level_0_0
                        )

                    # CALCULATE LIKELIHOOD FROM THE MEASURED M_plus and M_minus
                    likelihood = calculate_likelihood(M_measured, tau_vals, gamma_plus_arr, gamma_minus_arr, n_gamma)

                    # CALCULATE POSTERIOR
                    prior_plus_unnorm = likelihood * prior_plus
                    prior_minus_unnorm = likelihood * prior_minus

                    # NORMALIZE POSTERIOR TO MAKE UPDATED PRIOR
                    prior_plus = prior_plus_unnorm / np.sum(prior_plus_unnorm)
                    prior_minus = prior_minus_unnorm / np.sum(prior_minus_unnorm)


    def tau_opt_from_cost_function(
        tau_plus_arr,
        tau_minus_arr,
        gamma_plus_mean,
        gamma_minus_mean,
        gamma_plus_sigma,
        gamma_minus_sigma,
        T0,
    ):

        cost_function = (
            (
                (gamma_plus_sigma / gamma_plus_mean) ** 2
                + (gamma_minus_sigma / gamma_minus_mean) ** 2
            )
            ** 0.5
        ) * T
        idx = np.argmin(cost_function)

        return tau_plus_arr[idx[0]], tau_minus_arr[idx[1]]


    def calculate_likelihood(
        M_measured, 
        tau_vals, 
        gamma_plus_arr, 
        gamma_minus_arr, 
        n_gamma
    ):

        M_plus_measured = np.mean(M_measured[0]) * np.ones((n_gamma, n_gamma))
        sigma_M_plus_measured = np.std(M_measured[0])

        M_minus_measured = np.mean(M_measured[1]) * np.ones((n_gamma, n_gamma))
        sigma_M_minus_measured = np.std(M_measured[1])

        tau_plus, tau_minus = tau_vals

        M_plus_tilde, M_minus_tilde = np.zeros(
            n_gamma, n_gamma
        )

        for gamma_plus in gamma_plus_arr:

            tilde_plus_col, tilde_minus_col = []

            for gamma_minus in gamma_minus_arr:

                g = (
                    (gamma_plus_arr**2)
                    + (gamma_minus_arr**2)
                    - (gamma_plus * gamma_minus)
                ) ** 0.5

                beta_plus = gamma_plus + gamma_minus + g
                beta_minus = gamma_plus + gamma_minus - g

                model_func_plus = (((g + gamma_plus) * np.exp(- beta_plus * tau_plus)) + ((g- gamma_plus) * np.exp(- beta_minus * tau_plus))) / g
                model_func_minus = (((g + gamma_minus) * np.exp(- beta_plus * tau_minus)) + ((g - gamma_minus) * np.exp(- beta_minus * tau_minus))) / g

                tilde_plus_col.append(model_func_plus)
                tilde_minus_col.append(model_func_minus)

            M_plus_tilde.append(tilde_plus_col)
            M_minus_tilde.append(tilde_minus_col)

        chi_plus = (M_plus_measured - M_plus_tilde) / ((2 ** 0.5) * (sigma_M_plus_measured))
        chi_minus = (M_minus_measured - M_plus_tilde) / ((2 ** 0.5) * (sigma_M_plus_measured))

        likelihood = np.exp(-((chi_plus ** 2) + (chi_minus ** 2)))

        return likelihood













                # # SAVE CURRENT DATA TO DATA SERVER
                # cwODMRdata.push({'params': {'datasetName': datasetName,
                #                 'samplingFreq': samplingFreq,
                #                 'maxIterations': maxIterations,
                #                 'startFreq': startFreq, 'endFreq': endFreq,
                #                 'numFreqs': numFreqs, 'rfPower': rfPower,
                #                 # 'preReadoutLaserAndMwTime': preReadoutLaserAndMwTime, 'laserAndMwReadOutTime': laserAndMwReadOutTime,
                #                 # 'extraLaserInitTime': extraLaserInitTime, 'waitTime': waitTime,
                #                 'num_samples' : num_samples,
                #                 'clock_time' : clock_time,
                #                 'probe_time' : probe_time},
                #                 'title': 'CW ODMR',
                #                 'xlabel': 'Freq',
                #                 'ylabel': 'Counts',
                #                 'datasets': {'freqs': self.freqs,
                #                         'mwCountsDict': self.mwCountsDict,
                #                         'noMwCountsDict': self.noMwCountsDict
                #                 }
                # })

                # # RESET SWABIAN OUTPUTS
                # gw.swabian.reset()

            print("Experiment finished!")



