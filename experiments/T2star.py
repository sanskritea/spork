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
from guis.guiElements_general import flexSave

from experiments.NewportSpatialFeedback import SpatialFeedback
from drivers.ni.nidaq_final import NIDAQ


class T2_STAR_Meas:

    def T2star(
        self,
        datasetName: str,
        num_samples: int,
        maxIterations: int,
        rf_power: float,
        laser_power: float,
        freq: float,
        pi_time: int,
        tau_min: int,
        tau_max: int,
        tau_num: int,    
    ):

        # connect to the instrument server
        # connect to the data server and create a data set, or connect to an
        # existing one with the same name if it was created earlier.
        with InstrumentGateway() as gw, DataSource(datasetName) as T2starData:

            # setup a relevant iterator depending on maxIterations
            if maxIterations < 0:
                iters = count()  # infinite iterator
            else:
                iters = range(maxIterations)

            readout_time = 500

            # list of tau wait times 
            # self.tau_list = np.geomspace(tau_min, tau_max, int(tau_num))
            self.tau_list = np.linspace(tau_min, tau_max, int(tau_num))

            # storing experiment data
            self.MW_ON_counts = dict([[tau, []] for tau in self.tau_list])
            self.MW_OFF_counts = dict([[tau, []] for tau in self.tau_list])

            # SRS actions
            gw.sg.set_rf_amplitude(rf_power)  # set ouput power
            gw.sg.set_frequency(freq)
            gw.sg.set_mod_state(False)
            gw.sg.set_rf_state("1")

            # Feedback parameters
            feedback_trigger_rate = int(20e3)

            # generate pulse sequences
            seqs = []
            for tau in self.tau_list:
                seqs.append(Pulses(gw).T2_STAR(pi_time, int(tau), tau_max))


            # MAIN EXPERIMENT LOOP
            with NIDAQ() as mynidaq:

                # set laser power
                mynidaq.laser_power_atten(laser_power)

                # # Feedback parameters
                # feedback_trigger_rate = int(20e3)
                # feedback_time_per_point = 0.05
                # feedback_num_samples = int(feedback_trigger_rate * feedback_time_per_point)
                # x_init_position = 3.4192
                # y_init_position = 3.1298
                # z_init_position = 7.7909
                # self.feedback_timer = time.time()
                # feedback_counter = 0

                for i in iters:
                    print("Starting iteration ", i)

                    # # SPATIAL FEEDBACK EVERY 5 minutes
                    # if ((time.time() - self.feedback_timer) > 300):
                    #     feedback_counter = feedback_counter + 1
                    #     print('Feedback')
                    #     begin_feedback = time.time()
                    #     SpatialFeedback.Feedback(x_init_position, y_init_position, z_init_position)
                    #     feedback_duration = time.time() - begin_feedback
                    #     print('Feedback duration: ', feedback_duration)
                    #     print('Feedback counter ', feedback_counter)
                    #     self.feedback_timer = time.time()

                    # T2_STAR MEASUREMENT
                    for j in range(int(tau_num)):

                        print("Measuring tau ", self.tau_list[j], " ns")

                        # check counts before measuring
                        gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(feedback_trigger_rate))
                        current_counts = np.mean(obtain(mynidaq.internal_read_task(feedback_trigger_rate, int(feedback_trigger_rate * 0.01)))) / (1 / feedback_trigger_rate)
                        print('Current counts ', current_counts)
                        gw.swabian.reset()

                        # START TASK
                        tau = self.tau_list[j]
                        # print('tau ', tau, 'ns')
                        # num_samples = int(readout_time  * 1e-9 * 20e6)
                        mynidaq.start_external_read_task(20e6, ((4 * num_samples) + 1))

                        # START PULSESTREAMER
                        gw.swabian.runSequenceInfinitely(seqs[j])
                        # gw.swabian.runSequenceInfinitely(Pulses(gw).OPTICAL_T1(tau, clock_time, init_time, readout_time, laser_lag, singlet_decay, tau_max))

                        # COUNT WITH EXTERNAL TRIGGER
                        raw_counts = (obtain(mynidaq.external_read_task(20e6, ((4 * num_samples) + 1))))

                        # SEPARATING COUNTS INTO MW ON / MW OFF
                        MW_ON_counts = 1e9 * np.mean(raw_counts[2::4]) / readout_time
                        # bright_normalization = 1e9 * np.mean(raw_counts[4::8]) / readout_time
                        MW_OFF_counts = 1e9 * np.mean(raw_counts[0::4]) / readout_time
                        # dark_normalization = 1e9 * np.mean(raw_counts[0::8]) / readout_time
                        self.MW_ON_counts[tau].append(MW_ON_counts) # / bright_normalization)
                        self.MW_OFF_counts[tau].append(MW_OFF_counts) # / dark_normalization)

                        # RESET SWABIAN OUTPUTS
                        gw.swabian.reset()

                        # SAVE DATA TO DATASERVER
                        T2starData.push(
                            {
                                "params": {
                                    "datasetName": datasetName,
                                    "num_samples": num_samples,
                                    "maxIterations": maxIterations,
                                    "rf_power": rf_power,
                                    "laser_power": laser_power,
                                    "freq": freq,
                                    "pi_time": pi_time,
                                    "tau_min": tau_min,
                                    "tau_max": tau_max,
                                    "tau_num": tau_num,
                                },
                                "title": "T2_star",
                                "xlabel": "Tau Time (s)",
                                "ylabel": "PL (cts/s)",
                                "datasets": {
                                    "taus": self.tau_list,
                                    "MW_ON_counts": self.MW_ON_counts,
                                    "MW_OFF_counts": self.MW_OFF_counts,
                                },
                            }
                        )

                    notes = 'Iters(' + str(i) + ')'
                    flexSave(datasetName, notes, 'T2star')
                
                print("Experiment finished!")



