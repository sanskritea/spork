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


class MW_T1_Meas:

    def MWT1(
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
        with InstrumentGateway() as gw, DataSource(datasetName) as MWT1Data:

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
            # tau_short = np.linspace(tau_min, 1000000, int((tau_num - 10)), endpoint = False)
            # # tau_med = np.linspace(tau_min * 100, tau_min * 400, int((tau_num / 3) - 5), endpoint = False)
            # tau_long = np.linspace(1000000, tau_max, 30, endpoint = True )
            # self.tau_list = np.concatenate([tau_short, tau_long])
            self.tau_list = np.geomspace(tau_min, tau_max, tau_num)
            tau_list_size = int(np.size(self.tau_list))

            # print('self.tau_list ', self.tau_list)
            # print('np.size(self.tau_list) ', np.size(self.tau_list))
            # self.tau_list = np.geomspace(tau_min, tau_max, int(tau_num))
            # self.tau_list = np.linspace(tau_min, tau_max, int(tau_num))

            # storing experiment data
            self.BrightCounts = dict([[tau, []] for tau in self.tau_list])
            self.DarkCounts = dict([[tau, []] for tau in self.tau_list])
            self.raw_bright_counts_list = []
            self.raw_dark_counts_list = []

            # SRS actions
            gw.sg.set_rf_amplitude(rf_power)  # set ouput power
            gw.sg.set_frequency(freq)
            gw.sg.set_mod_state(False)
            gw.sg.set_rf_state("1")

            # generate pulse sequences
            seqs = []
            for tau in self.tau_list:
                seqs.append(Pulses(gw).MW_T1(int(tau), clock_time, init_time, readout_time, laser_lag, singlet_decay, tau_max, pi_time))

            # Feedback parameters
            feedback_trigger_rate = int(20e3)
            feedback_time_per_point = 0.05
            feedback_num_samples = int(feedback_trigger_rate * feedback_time_per_point)
            objective_x_init_position = 4.6186
            objective_y_init_position = 11.3904
            objective_z_init_position = -0.935003
            feedback_timer = time.time()
            feedback_counter = 0
            start_time = time.time()

            for i in iters:

                print('Iteration number ', i)

                raw_bright_counts = []
                raw_dark_counts = []
                
                # MAIN EXPERIMENT LOOP
                for j in range(tau_list_size):

                    tau = self.tau_list[j]
                    print('tau: ', tau)

                    with NIDAQ() as mynidaq:

                        # check counts before measuring
                        gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(feedback_trigger_rate))
                        current_counts = np.mean(obtain(mynidaq.internal_read_task(feedback_trigger_rate, int(feedback_trigger_rate * 0.01)))) / (1 / feedback_trigger_rate)
                        print('Current counts ', current_counts)
                        gw.swabian.reset()

                        # set laser power
                        mynidaq.laser_power_atten(laser_power)

                        # START TASK
                        mynidaq.start_external_read_task(samplingFreq, ((4 * num_samples) + 1))

                        # START PULSESTREAMER
                        gw.swabian.runSequenceInfinitely(seqs[j])

                        # COUNT WITH EXTERNAL TRIGGER
                        raw_counts = (obtain(mynidaq.external_read_task(samplingFreq, ((4 * num_samples) + 1))))

                        # SEPARATING COUNTS INTO MW ON / MW OFF
                        raw_bright_counts.append(raw_counts[2::4])
                        raw_dark_counts.append(raw_counts[0::4])
                        bright_counts = 1e9 * np.mean(raw_bright_counts[j]) / readout_time
                        dark_counts = 1e9 * np.mean(raw_dark_counts[j]) / readout_time
                        diff = bright_counts - dark_counts
                        self.BrightCounts[tau].append(bright_counts)
                        self.DarkCounts[tau].append(dark_counts)

                        # print('raw bright counts for iteration number ', i, 'and tau ', tau, ' : ', raw_bright_counts[j])

                        # SAVE DATA TO DATASERVER
                        MWT1Data.push(
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
                                "ylabel": "PL (cts/s)",
                                "datasets": {
                                    "taus": self.tau_list,
                                    "BrightCounts": self.BrightCounts,
                                    "DarkCounts": self.DarkCounts,
                                    "raw_bright_counts": self.raw_bright_counts_list,
                                    "raw_dark_counts": self.raw_dark_counts_list
                                },
                            }
                        )

                        # RESET SWABIAN OUTPUTS
                        gw.swabian.reset()

                        # SPATIAL FEEDBACK EVERY 10mins
                        if ((time.time() - feedback_timer) > 3000000):

                            # Parameters
                            trigger_rate = int(20e3)
                            time_per_point = 0.003
                            begin_feedback = time.time()
                            objective_x_init_position, objective_y_init_position, current_max_counts = SpatialFeedback.Feedback(objective_x_init_position, objective_y_init_position, objective_z_init_position)
                            print('Current objective_x_position ', objective_x_init_position)
                            print('Current objective_y_position ', objective_y_init_position)
                            print('Returned max counts ', current_max_counts)
                            
                            gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(int(trigger_rate)))
                            num_samples = int(trigger_rate * time_per_point)
                            post_feedback_counts = np.mean(obtain(mynidaq.internal_read_task(int(trigger_rate), num_samples))) / (1 / trigger_rate)
                            print('Post feedback counts ', post_feedback_counts)
                            gw.swabian.reset()

                            gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(int(trigger_rate)))
                            feedback_duration = time.time() - begin_feedback
                            feedback_timer = time.time()
                            feedback_counter = feedback_counter + 1
                            print('Feedback duration: ', feedback_duration)
                            print('Feedback counter ', feedback_counter)
                            feedback_timer = time.time()

                # print('raw bright counts for iteration number ', i, ' : ', raw_bright_counts)            

                # update raw counts compendium
                # self.raw_bright_counts_list.append(raw_bright_counts)
                # self.raw_dark_counts_list.append(raw_dark_counts)

                # # print('self.raw_bright_counts_list ', self.raw_bright_counts_list)

                # # SAVE DATA TO DATASERVER
                # MWT1Data.push(
                #     {
                #         "params": {
                #             "datasetName": datasetName,
                #             "samplingFreq": samplingFreq,
                #             "maxIterations": maxIterations,
                #             "freq": freq,
                #             "tau_min": tau_min,
                #             "tau_max": tau_max,
                #             "tau_num": tau_num,
                #             "rf_power": rf_power,
                #             "laser_power": laser_power,
                #             "num_samples": num_samples,
                #             "clock_time": clock_time,    
                #             "init_time": init_time,    
                #             "laser_lag": laser_lag,    
                #             "readout_time": readout_time,    
                #             "singlet_decay": singlet_decay, 
                #             "pi_time": pi_time
                #         },
                #         "title": "MW T1",
                #         "xlabel": "Tau Time (s)",
                #         "ylabel": "PL (cts/s)",
                #         "datasets": {
                #             "taus": self.tau_list,
                #             "BrightCounts": self.BrightCounts,
                #             "DarkCounts": self.DarkCounts,
                #             "raw_bright_counts": self.raw_bright_counts_list,
                #             "raw_dark_counts": self.raw_dark_counts_list
                #         },
                #     }
                # )

            flexSave(datasetName, 'MWT1', 'final')
            print('Total time taken ', time.time() - start_time)
            print("Experiment finished!")



