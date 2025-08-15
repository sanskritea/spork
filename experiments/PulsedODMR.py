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
from guis.guiElements_general import flexSave

from rpyc.utils.classic import obtain

from experiments.NewPulses import Pulses
from drivers.ni.nidaq_final import NIDAQ

# from PulsePatterns import Pulses


class Pulsed_ODMR_Measurement:

    def PulsedODMR(
        self,
        datasetName: str,
        num_samples: int,
        maxIterations: int,
        rfPower: float,
        laser_power: float,
        startFreq: float,
        endFreq: float,
        numFreqs: int,
        pi_time: int

    ):
        """Run a CW ODMR experiment
        Arguments:  *
        """  # TODO: Fill in arg documentation

        # connect to the instrument server
        # connect to the data server and create a data set, or connect to an
        # existing one with the same name if it was created earlier.

        with InstrumentGateway() as gw, DataSource(datasetName) as PulsedODMRdata:

            # SRS actions
            gw.sg.set_rf_amplitude(rfPower)  # set ouput power
            gw.sg.set_mod_state(False)
            print('probe time ', probe_time)

            # set up frequencies to sweep over
            # self.freqs = np.random.permutation(np.linspace(startFreq, endFreq, numFreqs, endpoint=True)) # permute the order of freqs to minimize time effects
            self.freqs = np.linspace(startFreq, endFreq, numFreqs, endpoint=True)
            # self.freqDict = dict([self.freqs])
            # for storing the experiment data
            self.mwCountsDict = dict(
                [[freq, []] for freq in self.freqs]
            )  # should make dict of {freq: [counts1, counts2, ...], ...} pairs
            self.noMwCountsDict = dict([[freq, []] for freq in self.freqs])

            # setup a relevant iterator depending on maxIterations
            if maxIterations < 0:
                iters = count()  # infinite iterator
            else:
                iters = range(maxIterations)

            # MAIN EXPERIMENT LOOP
            with NIDAQ() as mynidaq:

                # set laser power
                mynidaq.laser_power_atten(laser_power)

                for i in iters:

                    for freq in self.freqs:

                        # SET SRS OUTPUT FREQUENCY
                        gw.sg.set_frequency(freq)
                        print("freq : ", freq)
                        gw.sg.set_rf_state("1")

                        # START TASK
                        # num_samples = int(readout_time  * 1e-9 * 20e6)
                        mynidaq.start_external_read_task(20e6, ((4 * num_samples) + 1))

                        # START PULSESTREAMER
                        gw.swabian.runSequenceInfinitely(
                            Pulses(gw).PULSED_ODMR(pi_time)
                        )

                        # COUNTING WITH EXTERNAL TRIGGER
                        print('counting')
                        raw_counts = (obtain(mynidaq.external_read_task(20e6, ((4 * num_samples) + 1))))

                        # SEPARATING COUNTS INTO MW ON / MW OFF
                        # print('counts : ', counts)
                        signal_counts = 1e9 * np.mean(raw_counts[2::4]) / probe_time
                        bg_counts = 1e9 * np.mean(raw_counts[0::4]) / probe_time
                        self.mwCountsDict[freq].append(signal_counts)
                        self.noMwCountsDict[freq].append(bg_counts)

                        # SAVE CURRENT DATA TO DATA SERVER
                        PulsedODMRdata.push(
                            {
                                "params": {
                                    "datasetName": datasetName,
                                    "num_samples": num_samples,
                                    "maxIterations": maxIterations,
                                    "rf_power": rf_power,
                                    "laser_power": laser_power,
                                    "startFreq": startFreq,
                                    "endFreq": endFreq,
                                    "numFreqs": numFreqs,
                                    "pi_time": pi_time
                                },
                                "title": "Pulsed ODMR",
                                "xlabel": "Freq",
                                "ylabel": "Counts",
                                "datasets": {
                                    "freqs": self.freqs,
                                    "mwCountsDict": self.mwCountsDict,
                                    "noMwCountsDict": self.noMwCountsDict,
                                },
                            }
                        )

                        # RESET SWABIAN OUTPUTS
                        gw.swabian.reset()

                notes = ''
                flexSave(datasetName, notes, 'Pulsed ODMR')

            print("Experiment finished!")


# if __name__ == "__main__":

#     exp = Pulsed_ODMR_Measurement()
#     print(
#         """Running ODMR with 1Hz sampling rate, 10us MW+lasers before counts, 50us of counting, 20us of laser reinit, 
#     30us wait with 1mW of laser power. Saving to CW_ODMR on dataserv"""
#     )
#     exp.cwODMR(
#         datasetName="CW_ODMR",
#         sampleFreq=1,
#         maxIterations=3,
#         startFreq=2.75e9,
#         endFreq=2.95e9,
#         numFreqs=5,
#         rfPower=-17,
#         # preReadoutLaserAndMwTime=10000, laserAndMwReadOutTime=50000,
#         # extraLaserInitTime=20000, waitTime=30000, debug=False
#         num_samples=100000,
#         clock_time=11,
#         probe_time=50000,
#     )
#     print("Completed cwODMR.py")
