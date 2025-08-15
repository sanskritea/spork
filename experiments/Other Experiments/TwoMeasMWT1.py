"""
Experiment to run T1 measuerment on the NV-AFM Spork setup


"""

import time
from itertools import count

import numpy as np
import math

from nspyre import DataSource
from nspyre import InstrumentGateway
from guis.guiElements_general import flexSave

from rpyc.utils.classic import obtain

from experiments.NewPulses import Pulses
from experiments.NewportSpatialFeedback import SpatialFeedback
from drivers.ni.nidaq_final import NIDAQ


class Two_Meas_MW_T1_Meas:

    def TwoMeasMWT1(
        self,
        datasetName: str,
        # samplingFreq: float,
        # maxIterations: int,
        # freq: float,
        tau_min: int,
        tau_max: int,
        tau_num: int,
        rf_power: float,
        laser_power: float,
        num_samples: int,
        clock_time: int,    # DAQ counting trigger pulse (from Swabian) duration, usually 10ns
        init_time: int,     # NV initialization laser ON duration
        laser_lag: int,     # laser stabilizing time, usually ~100ns
        readout_time: int,    # readout laser ON time
        singlet_decay: int, # NV singlet state emptying duration
        # pi_time: int,
    ):

        # connect to the instrument server
        # connect to the data server and create a data set, or connect to an
        # existing one with the same name if it was created earlier.
        with InstrumentGateway() as gw, DataSource(datasetName) as TwoMeasMWT1Data:

            # # setup a relevant iterator depending on maxIterations
            # if maxIterations < 0:
            #     iters = count()  # infinite iterator
            # else:
            #     iters = range(maxIterations)

            # parameters for the two peaks
            freqs = [2.92e9, 2.85e9]
            pi_times = [45, 55]

            # list of tau wait times 
            # tau_short = np.linspace(tau_min, 1000000, int((tau_num - 5)), endpoint = False)
            # tau_long = np.linspace(1000000, tau_max, 5, endpoint = True )
            # self.tau_list = np.concatenate([tau_short, tau_long])
            self.tau_list = np.geomspace(tau_min, tau_max, tau_num)
            tau_list_size = int(np.size(self.tau_list))

            # storing experiment data
            self.M_plus = np.zeros(tau_list_size)
            self.M_minus = np.zeros(tau_list_size)
            self.M_plus_err = []
            self.M_minus_err = []
            self.S_0_0_0_plus = []
            self.S_0_0_tau_plus = []
            self.S_level_0_0_plus = []
            self.S_level_0_tau_plus = []
            self.S_0_0_0_minus = []
            self.S_0_0_tau_minus = []
            self.S_level_0_0_minus = []
            self.S_level_0_tau_minus = []
            self.elapsed_time = np.zeros(0)

            print('np.size(self.M_plus) ', np.size(self.M_plus))

            # Feedback parameters
            feedback_trigger_rate = int(20e3)
            feedback_time_per_point = 0.05
            feedback_num_samples = int(feedback_trigger_rate * feedback_time_per_point)
            x_init_position = 7.419
            y_init_position = 10.3064
            z_init_position = 4.5083
            feedback_timer = time.time()
            feedback_counter = 0
            start_time = time.time()         


            # MAIN EXPERIMENT LOOP
            for j in range(tau_list_size):

                tau = self.tau_list[j]
                print(j, 'th tau (ms) : ', tau / 1e6)

                print('Counting')
                with NIDAQ() as mynidaq:

                    # set laser power
                    mynidaq.laser_power_atten(laser_power)

                    # temporary storage lists
                    M = []
                    M_err = []
                    S_0_0_0 = []
                    S_0_0_tau = []
                    S_level_0_0 = []
                    S_level_0_tau = []

                    # Measure both peaks one by one
                    for k in range(2):

                        # SRS actions
                        gw.sg.set_rf_amplitude(rf_power)  # set ouput power
                        gw.sg.set_frequency(freqs[k])
                        gw.sg.set_mod_state(False)
                        gw.sg.set_rf_state("1")

                        # START TASK
                        mynidaq.start_external_read_task(20e6, ((8 * num_samples) + 1))

                        # PULSE SEQUENCE
                        gw.swabian.runSequenceInfinitely(Pulses(gw).BAYESIAN_T1(int(tau), clock_time, init_time, readout_time, laser_lag, singlet_decay, pi_times[k]))

                        # COUNT WITH EXTERNAL TRIGGER
                        counts = (obtain(mynidaq.external_read_task(20e6, ((8 * num_samples) + 1))))

                        S_0_0_0.append(1e9 * counts[0::8] / readout_time)
                        S_0_0_tau.append(1e9 * counts[2::8] / readout_time)
                        S_level_0_0.append(1e9 * counts[4::8] / readout_time)
                        S_level_0_tau.append(1e9 * counts[6::8] / readout_time)

                        S_0_0_0_mean = np.mean(S_0_0_0[k])
                        S_0_0_tau_mean = np.mean(S_0_0_tau[k]) 
                        S_level_0_0_mean = np.mean(S_level_0_0[k]) 
                        S_level_0_tau_mean = np.mean(S_level_0_tau[k]) 

                        S_0_0_0_err = np.std(S_0_0_0[k]) / (num_samples ** 0.5)
                        S_0_0_tau_err = np.std(S_0_0_tau[k]) / (num_samples ** 0.5)
                        S_level_0_0_err = np.std(S_level_0_0[k]) / (num_samples ** 0.5)
                        S_level_0_tau_err = np.std(S_level_0_tau[k]) / (num_samples ** 0.5)

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

                        z = top * inv_bottom
                        z_err = ((top * inv_bottom_err)**2 + (inv_bottom_err * top_err) ** 2) ** 0.5

                        M.append(z)
                        M_err.append(z_err)
                        print('M_mean ', z)
                        print('M_err ', z_err)

                        # Reset Swabian
                        gw.swabian.reset()

                        # SPATIAL FEEDBACK EVERY 5mins
                        if ((time.time() - feedback_timer) > 300):
                            feedback_counter = feedback_counter + 1
                            print('Feedback')
                            begin_feedback = time.time()
                            SpatialFeedback.Feedback(x_init_position, y_init_position, z_init_position)
                            feedback_duration = time.time() - begin_feedback
                            print('Feedback duration: ', feedback_duration)
                            print('Feedback counter ', feedback_counter)
                            feedback_timer = time.time()

                    # update final data lists
                    self.M_plus[j] = M[0]
                    self.M_minus[j] = M[1]
                    self.M_plus_err.append(M_err[0])
                    self.M_minus_err.append(M_err[1])
                    self.S_0_0_0_plus.append(S_0_0_0[0])
                    self.S_0_0_tau_plus.append(S_0_0_tau[0])
                    self.S_level_0_0_plus.append(S_level_0_0[0])
                    self.S_level_0_tau_plus.append(S_level_0_tau[0])
                    self.S_0_0_0_minus.append(S_0_0_0[1])
                    self.S_0_0_tau_minus.append(S_0_0_tau[1])
                    self.S_level_0_0_minus.append(S_level_0_0[1])
                    self.S_level_0_tau_minus.append(S_level_0_tau[1]) 

                    # SAVE DATA TO DATASERVER
                    TwoMeasMWT1Data.push(
                        {
                            "params": {
                                "datasetName": datasetName,
                                # "samplingFreq": samplingFreq,
                                # "maxIterations": maxIterations,
                                # "freq": freq,
                                "tau_min": tau_min,
                                "tau_max": tau_max,
                                "tau_num": tau_num,
                                "rf_power": rf_power,
                                "laser_power": laser_power,
                                "num_samples": num_samples,
                                "clock_time": clock_time,    
                                "init_time": init_time,    
                                "laser_lag": laser_lag,    
                                "readout_time": readout_time,    
                                "singlet_decay": singlet_decay, 
                                # "pi_time": pi_time
                            },
                            "title": "MW T1",
                            "xlabel": "Tau Time (ms)",
                            "ylabel": "AU",
                            "datasets": {
                                "taus": self.tau_list,
                                'M_plus': self.M_plus,
                                'M_minus': self.M_minus,
                                'M_plus_err': self.M_plus_err,
                                'M_minus_err': self.M_minus_err,
                                'S000Plus': self.S_0_0_0_plus,
                                'S000Minus': self.S_0_0_0_minus,
                                'S00tauPlus': self.S_0_0_tau_plus,
                                'S00tauMinus': self.S_0_0_tau_minus,
                                'Slevel00Plus': self.S_level_0_0_plus,
                                'Slevel00Minus': self.S_level_0_0_minus,
                                'Slevel0tauPlus': self.S_level_0_tau_plus,
                                'Slevel0tauMinus': self.S_level_0_tau_minus,
                            },
                        }
                    )

                    flexSave(datasetName, 'TwoMeasT1', 'final')

            print('Total time taken ', time.time() - start_time)
            print("Experiment finished!")



