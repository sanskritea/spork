"""
This is an application to run CW ODMR on Jasper

Copyright (c) April 2023, C. Egerstrom
All rights reserved.

This work is licensed under the terms of the 3-Clause BSD license.
For a copy, see <https://opensource.org/licenses/BSD-3-Clause>.

Modified: PMN July '23
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


class Reinit_Test_Measurement:

    def ReinitTest(
        self,
        datasetName: str,
        samplingFreq: float,
        maxIterations: int,
        tau_min: int,
        # tau_max: int,
        laser_power: float,
        num_samples: int,
        clock_time: int,    # DAQ counting trigger pulse (from Swabian) duration, usually 10ns
        init_time: int,     # NV initialization laser ON duration
        readout_time: int,
    ):

        with InstrumentGateway() as gw, DataSource(datasetName) as ReinitTestData:

            # setup a relevant iterator depending on maxIterations
            if maxIterations < 0:
                iters = count()  # infinite iterator
            else:
                iters = range(maxIterations)

            # list of tau wait times 
            print('tau_max ', init_time - (2 * clock_time) - readout_time)
            self.tau_list = np.arange(tau_min, (init_time - (2 * clock_time) - readout_time), int(readout_time))
            length = init_time - (2 * clock_time) - readout_time
            print('length ', length)

            # storing experiment data
            self.BrightCounts = dict([[tau, []] for tau in self.tau_list])
            self.DarkCounts = dict([[tau, []] for tau in self.tau_list])

            # generate pulse sequences
            seqs = []
            # rise_laser_off = 500
            # fall_laser_off = 10000000
            print('creating tau list')
            seqs.append(Pulses(gw).REINIT_TEST(clock_time, init_time, readout_time, singlet_decay, pi_time, length))


            # MAIN EXPERIMENT LOOP
            with NIDAQ() as mynidaq:

                # set laser power
                mynidaq.laser_power_atten(laser_power)
                for ii in iters:

                    # for tau in self.tau_list:  # measure NV with delay = tau
                    # for i in range(length):

                    #     tau = self.tau_list[i]
                        # print('tau ', tau)

                    # START TASK
                    print('start task')
                    mynidaq.start_external_read_task(samplingFreq, ((2 * length) + 1))

                    # START PULSESTREAMER
                    print('pulsing')
                    gw.swabian.runSequenceInfinitely(seqs[i])

                    # COUNT WITH EXTERNAL TRIGGER
                    print('counting')
                    raw_counts = obtain(mynidaq.external_read_task(samplingFreq, ((2 * length) + 1)))
                    bright_counts = 1e9 * np.mean(raw_counts[0::4]) / readout_time
                    dark_counts = 1e9 * np.mean(raw_counts[2::4]) / readout_time
                    self.BrightCounts[tau].append(bright_counts)
                    self.DarkCounts[tau].append(dark_counts)

                    # RESET SWABIAN OUTPUTS
                    gw.swabian.reset()

                    # SAVE DATA TO DATASERVER
                    ReinitTestData.push(
                        {
                            "params": {
                                "datasetName": datasetName,
                                "samplingFreq": samplingFreq,
                                "maxIterations": maxIterations,
                                "tau_min": tau_min,
                                "laser_power": laser_power,
                                "num_samples": num_samples,
                                "clock_time": clock_time,    
                                "init_time": init_time,    
                                "readout_time": readout_time,    
                            },
                            "title": "Reinit Test",
                            "xlabel": "Tau Time (s)",
                            "ylabel": "PL (cts/s)",
                            "datasets": {
                                "taus": self.tau_list,
                                "BrightCounts": self.BrightCounts,
                                "DarkCounts": self.DarkCounts,
                            },
                        }
                    )
                    
            # mynidaq.laser_power_atten(1)
            print("Experiment finished!")