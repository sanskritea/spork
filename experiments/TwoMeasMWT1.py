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
        samplingFreq: float,
        maxIterations: int,
        freq: float,
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
        pi_time: int,
    ):

        # connect to the instrument server
        # connect to the data server and create a data set, or connect to an
        # existing one with the same name if it was created earlier.
        with InstrumentGateway() as gw, DataSource(datasetName) as TwoMeasMWT1Data:

            print('MW T1 PARAMETERS')
            print('clock_time ', clock_time)
            print('init_time ', init_time)
            print('laser_lag ', laser_lag)
            print('readout_time ', readout_time)
            print('singlet_decay ', singlet_decay)


            # setup a relevant iterator depending on maxIterations
            if maxIterations < 0:
                iters = count()  # infinite iterator
            else:
                iters = range(maxIterations)

            # list of tau wait times 
            self.tau_list = np.geomspace(tau_min, tau_max, int(tau_num))
            # self.tau_list = np.linspace(tau_min, tau_max, int(tau_num))

            # storing experiment data
            self.M_measured = dict([[tau, []] for tau in self.tau_list])
            self.M_mean = dict([[tau, []] for tau in self.tau_list])
            self.M_std = dict([[tau, []] for tau in self.tau_list])
            S_0_0_0 = np.zeros(maxIterations)
            S_0_0_tau = np.zeros(maxIterations)
            S_level_0_0 = np.zeros(maxIterations)
            S_level_0_tau = np.zeros(maxIterations)

            # SRS actions
            gw.sg.set_rf_amplitude(rf_power)  # set ouput power
            gw.sg.set_frequency(freq)
            gw.sg.set_mod_state(False)
            gw.sg.set_rf_state("1")

            # generate pulse sequences
            seqs = []
            for tau in self.tau_list:
                seqs.append(Pulses(gw).BAYESIAN_T1(int(tau), clock_time, init_time, readout_time, laser_lag, singlet_decay, pi_time))

            # Feedback parameters
            feedback_trigger_rate = int(20e3)
            feedback_time_per_point = 0.05
            feedback_num_samples = int(feedback_trigger_rate * feedback_time_per_point)
            x_init_position = 7.4198
            y_init_position = 10.3078
            z_init_position = 4.5083
            feedback_timer = time.time()
            feedback_counter = 0
            start_time = time.time()

            for i in iters:

                print('Iteration number ', i)
                
                # MAIN EXPERIMENT LOOP
                for j in range(int(tau_num)):
                    print('Counting')
                    with NIDAQ() as mynidaq:
                        # set laser power
                        mynidaq.laser_power_atten(laser_power)

                        # START TASK
                        tau = self.tau_list[j]
                        mynidaq.start_external_read_task(samplingFreq, ((8 * num_samples) + 1))

                        # START PULSESTREAMER
                        gw.swabian.runSequenceInfinitely(seqs[j])

                        # COUNT WITH EXTERNAL TRIGGER
                        counts = (obtain(mynidaq.external_read_task(samplingFreq, ((8 * num_samples) + 1))))

                        # SEPARATE COUNTS S_(initial level)_(final level)_delay
                        # garbage_counts = np.mean(1e9 * counts[0::10] / readout_time)
                        # print('garbage_counts ', garbage_counts)
                        
                        S_0_0_0_list = 1e9 * counts[0::8] / readout_time
                        S_0_0_tau_list = 1e9 * counts[2::8] / readout_time 
                        S_level_0_0_list = 1e9 * counts[4::8] / readout_time
                        S_level_0_tau_list = 1e9 * counts[6::8] / readout_time

                        S_0_0_0_mean = np.mean(S_0_0_0_list)
                        S_0_0_tau_mean = np.mean(S_0_0_tau_list) 
                        S_level_0_0_mean = np.mean(S_level_0_0_list) 
                        S_level_0_tau_mean = np.mean(S_level_0_tau_list) 

                        S_0_0_0_std = np.sqrt(S_0_0_0_mean / num_samples)
                        S_0_0_tau_std = np.sqrt(S_0_0_tau_mean / num_samples)
                        S_level_0_0_std = np.sqrt(S_level_0_0_mean / num_samples)
                        S_level_0_tau_std = np.sqrt(S_level_0_tau_mean / num_samples)

                        # FINAL MEASURE AND MOMENTS (Appendix E)
                        # mean_delta = S_0_0_0_mean - S_level_0_0_mean
                        # sigma_delta_sq = (S_level_0_0_std ** 2) + (S_0_0_0_std ** 2)
                        # sigma_delta = sigma_delta_sq ** 0.5 

                        # print('mean_delta ', mean_delta)
                        # Z_max = (1 / (4 * sigma_delta_sq)) * ((((mean_delta ** 2) + (8 * sigma_delta_sq)) ** 0.5) - mean_delta)
                        # print('1/Z_max ', 1/Z_max)
                        # if np.abs(1/Z_max) < 1:
                        #     print('Correcting counting error')
                        #     Z_max = 1 / mean_delta
                        #     print('New Z_max ', Z_max)
                        # sigma_Z = ((Z_max ** 2) * sigma_delta) / ((2 - (Z_max * mean_delta)) ** 0.5)

                        # A_n = S_0_0_tau_mean - S_level_0_tau_mean
                        # sigma_A_n_sq = (S_level_0_tau_std ** 2) + (S_0_0_tau_std ** 2)
                        # sigma_A_n = sigma_A_n_sq ** 0.5 

                        # M_mean = A_n * Z_max
                        # M_std = M_mean * np.abs((((sigma_A_n / A_n) ** 2) + ((sigma_Z / Z_max) ** 2)) ** 0.5)
                        
                        # print('A_n ', A_n)
                        # # print('1/Z_max ', 1/Z_max)

                        top = S_0_0_tau_mean - S_level_0_tau_mean
                        bottom = S_0_0_0_mean - S_level_0_0_mean
                        top_err = (S_0_0_tau_std ** 2 + S_level_0_tau_std ** 2) ** 0.5
                        bottom_err = (S_0_0_0_std ** 2 + S_level_0_0_std ** 2) ** 0.5

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

                        M_mean = z
                        M_std = z_err
                        print('M_mean ', M_mean)
                        print('M_std ', M_std)

                        self.M_measured[tau].append(M_mean)
                        self.M_mean[tau].append(M_mean)
                        self.M_std[tau].append(M_std)
                        
                        gw.swabian.reset()

                    # SAVE DATA TO DATASERVER
                    TwoMeasMWT1Data.push(
                        {
                            "params": {
                                "datasetName": datasetName,
                                "samplingFreq": samplingFreq,
                                "maxIterations": maxIterations,
                                "freq": freq,
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
                                "pi_time": pi_time
                            },
                            "title": "MW T1",
                            "xlabel": "Tau Time (s)",
                            "ylabel": "AU",
                            "datasets": {
                                "taus": self.tau_list,
                                "M_measured": self.M_measured,
                                "M_mean": self.M_mean,
                                "M_std": self.M_std,
                                "S_0_0_0": S_0_0_0,
                                "S_0_0_tau": S_0_0_tau,
                                "S_level_0_0": S_level_0_0,
                                "S_level_0_tau": S_level_0_tau
                            },
                        }
                    )

                    

                    # SPATIAL FEEDBACK EVERY 5mins
                    if ((time.time() - feedback_timer) > 120):
                        feedback_counter = feedback_counter + 1
                        print('Feedback')
                        begin_feedback = time.time()
                        SpatialFeedback.Feedback(x_init_position, y_init_position, z_init_position)
                        feedback_duration = time.time() - begin_feedback
                        print('Feedback duration: ', feedback_duration)
                        print('Feedback counter ', feedback_counter)
                        feedback_timer = time.time()

            flexSave(datasetName, 'TwoMeasT1', 'final')
            print('Total time taken ', time.time() - start_time)
            print("Experiment finished!")



