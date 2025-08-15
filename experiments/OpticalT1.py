"""
Experiment to run T1 measuerment on the NV-AFM Spork setup


"""

import time
from itertools import count

import numpy as np
import math

from nspyre import DataSource
from nspyre import InstrumentGateway

from rpyc.utils.classic import obtain

from experiments.NewPulses import Pulses
from drivers.ni.nidaq_final import NIDAQ
from experiments.NewportSpatialFeedback import SpatialFeedback


class Optical_T1_Meas:

    def OpticalT1(
        self,
        datasetName: str,
        samplingFreq: float,
        maxIterations: int,
        # freq: float,
        tau_min: int,
        tau_max: int,
        tau_num: int,
        laser_power: float,
        num_samples: int,
        clock_time: int,    # DAQ counting trigger pulse (from Swabian) duration, usually 10ns
        init_time: int,     # NV initialization laser ON duration
        laser_lag: int,     # laser stabilizing time, usually ~100ns
        # probe_time: int,    # readout laser ON time
        singlet_decay: int, # NV singlet state emptying duration
    ):

        # connect to the instrument server
        # connect to the data server and create a data set, or connect to an
        # existing one with the same name if it was created earlier.
        with InstrumentGateway() as gw, DataSource(datasetName) as OpticalT1Data:

            # setup a relevant iterator depending on maxIterations
            if maxIterations < 0:
                iters = count()  # infinite iterator
            else:
                iters = range(maxIterations)

            # list of tau wait times 
            self.tau_list = np.geomspace(tau_min, tau_max, int(tau_num))
            # self.tau_list = np.linspace(tau_min, tau_max, int(tau_num))

            # storing experiment data
            self.PLcounts = dict([[tau, []] for tau in self.tau_list])

            # Feedback parameters
            # feedback_trigger_rate = int(20e3)
            # feedback_time_per_point = 0.05
            # feedback_num_samples = int(feedback_trigger_rate * feedback_time_per_point)
            x_init_position = 5.8124
            y_init_position = 18.8962
            z_init_position = -5.47199
            feedback_timer = time.time()
            feedback_counter = 0

            # readout time in ns
            readout_time = 300

            # generate pulse sequences
            # seqs = []
            # for tau in self.tau_list:
            #     seqs.append(Pulses(gw).OPTICAL_T1(tau, clock_time, init_time, readout_time, laser_lag, singlet_decay, tau_max))


            # MAIN EXPERIMENT LOOP
            for i in iters:

                print('Doing ', i, 'th iteration')

                # SPATIAL FEEDBACK EVERY 5 minutes
                if ((time.time() - feedback_timer) > 300):
                    feedback_counter = feedback_counter + 1
                    print('Feedback')
                    begin_feedback = time.time()
                    SpatialFeedback.Feedback(x_init_position, y_init_position, z_init_position)
                    feedback_duration = time.time() - begin_feedback
                    print('Feedback duration: ', feedback_duration)
                    print('Feedback counter ', feedback_counter)
                    feedback_timer = time.time()
            
                with NIDAQ() as mynidaq:

                    # set laser power
                    mynidaq.laser_power_atten(laser_power)
                        
                    # for tau in self.tau_list:  # measure NV with delay = tau
                    for j in range(int(tau_num)):

                        tau = self.tau_list[j]
                        print('tau ', tau / 1e6, 'ms')

                        # START TASK
                        # num_samples = int(readout_time  * 1e-9 * 20e6)
                        mynidaq.start_external_read_task(20e6, ((2 * num_samples) + 1))

                        # START PULSESTREAMER
                        # gw.swabian.runSequenceInfinitely(seqs[j])
                        gw.swabian.runSequenceInfinitely(Pulses(gw).OPTICAL_T1(tau, clock_time, init_time, readout_time, laser_lag, singlet_decay, tau_max))

                        # COUNT WITH EXTERNAL TRIGGER
                        raw_counts = (obtain(mynidaq.external_read_task(20e6, ((2 * num_samples) + 1))))
                        counts = 1e9 * np.mean(raw_counts[0::2]) / readout_time
                        self.PLcounts[tau].append(counts)

                        # RESET SWABIAN OUTPUTS
                        gw.swabian.reset()

                        # # with reference normalization
                        # ref = 1e9 * np.mean(raw_counts[2::4]) / readout_time
                        # # print('ref ', ref)
                        # counts = 1e9 * np.mean(raw_counts[0::4]) / readout_time
                        # # print('counts ', counts)
                        # if ref==0.0:
                        #     continue
                        # else:
                        #     norm_counts = counts / ref
                        #     tau = self.tau_list[i]
                        #     self.PLcounts[tau].append(norm_counts)

                        # SAVE DATA TO DATASERVER
                        OpticalT1Data.push(
                            {
                                "params": {
                                    "datasetName": datasetName,
                                    "samplingFreq": samplingFreq,
                                    "maxIterations": maxIterations,
                                    # "freq": freq,
                                    "tau_min": tau_min,
                                    "tau_max": tau_max,
                                    "tau_num": tau_num,
                                    "laser_power": laser_power,
                                    "num_samples": num_samples,
                                    "clock_time": clock_time,    
                                    "init_time": init_time,    
                                    "laser_lag": laser_lag,    
                                    # "probe_time": probe_time,    
                                    "singlet_decay": singlet_decay, 
                                },
                                "title": "Optical T1",
                                "xlabel": "Tau Time (s)",
                                "ylabel": "PL (cts/s)",
                                "datasets": {
                                    "taus": self.tau_list,
                                    "PLcounts": self.PLcounts,
                                },
                            }
                        )

                    

            # mynidaq.laser_power_atten(1)
            print("Experiment finished!")



