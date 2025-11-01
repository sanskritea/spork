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


class AOM_Lag_Meas:

    def AOMLag(
        self, 
        datasetName: str,
        samplingFreq: float,
        maxIterations: int,
        tau_min: int,
        tau_max: int,
        tau_num: int,
        laser_power: float,
        num_samples: int,
        clock_time: int,    # DAQ counting trigger pulse (from Swabian) duration, usually 10ns
        init_time: int,     # NV initialization laser ON duration
        readout_time: int,
    ):

        with InstrumentGateway() as gw, DataSource(datasetName) as AOMLagData:

            # setup a relevant iterator depending on maxIterations
            if maxIterations < 0:
                iters = count()  # infinite iterator
            else:
                iters = range(maxIterations)

            # list of tau wait times 
            self.tau_list = np.linspace(tau_min, tau_max, int(tau_num))

            # storing experiment data
            self.PLcounts = dict([[tau, []] for tau in self.tau_list])
            self.OtherPLcounts = dict([[tau, []] for tau in self.tau_list])

            # generate pulse sequences
            seqs = []
            rise_laser_off = 500
            fall_laser_off = 10000000
            print('Making tau list')
            for tau in self.tau_list:
                seqs.append(Pulses(gw).AOM_Lag(tau, clock_time, init_time, readout_time, rise_laser_off, fall_laser_off))


            # MAIN EXPERIMENT LOOP
            with NIDAQ() as mynidaq:

                # set laser power
                mynidaq.laser_power_atten(laser_power)
                start_time = time.time()
                for i in iters:

                    print('doing ', i, 'th iteration')

                    # for tau in self.tau_list:  # measure NV with delay = tau
                    for i in range(int(tau_num)):

                        tau = self.tau_list[i]
                        print('tau ', tau)

                        # START TASK
                        mynidaq.start_external_read_task(samplingFreq, ((2 * num_samples) + 1))

                        # START PULSESTREAMER
                        gw.swabian.runSequenceInfinitely(seqs[i])

                        # COUNT WITH EXTERNAL TRIGGER
                        raw_counts = obtain(mynidaq.external_read_task(samplingFreq, ((2 * num_samples) + 1)))
                        counts = 1e9 * np.mean(raw_counts[0::2]) / readout_time
                        other_counts = 1e9 * np.mean(raw_counts[1::2]) / (rise_laser_off + init_time + fall_laser_off - (2 * clock_time) - readout_time)
                        self.PLcounts[tau].append(counts)
                        self.OtherPLcounts[tau].append(other_counts)

                        # RESET SWABIAN OUTPUTS
                        psreset = time.time()
                        gw.swabian.reset()
                        # print('ps reset time ', time.time() - psreset)

                    # SAVE DATA TO DATASERVER
                    push_start = time.time()
                    AOMLagData.push(
                        {
                            "params": {
                                "datasetName": datasetName,
                                "samplingFreq": samplingFreq,
                                "maxIterations": maxIterations,
                                "tau_min": tau_min,
                                "tau_max": tau_max,
                                "tau_num": tau_num,
                                "laser_power": laser_power,
                                "num_samples": num_samples,
                                "clock_time": clock_time,    
                                "init_time": init_time,    
                                "readout_time": readout_time,    
                            },
                            "title": "AOM Lag",
                            "xlabel": "Tau Time (s)",
                            "ylabel": "PL (cts/s)",
                            "datasets": {
                                "taus": self.tau_list,
                                "PLcounts": self.PLcounts,
                                # "OtherPLcounts": self.OtherPLcounts,
                            },
                        }
                    )

                    

            # mynidaq.laser_power_atten(1)
            print("Experiment finished!")
            print('time taken ', time.time() - start_time)

