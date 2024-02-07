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


class CW_ODMR_Measurement:

    def cwODMR(self, 
                datasetName: str, samplingFreq: float, maxIterations: int,
               startFreq: float, endFreq: float, numFreqs: int, rfPower: float, laser_power: float,
               # preReadoutLaserAndMwTime: int, laserAndMwReadOutTime: int, 
               # extraLaserInitTime: int, waitTime: int,
               # debug=False):
               num_samples: int, clock_time: int, probe_time: int, wait_time: int,
               # ivalon: float, ivaloff: float, qvalon: float, qvaloff: float, frequency: float,
               ):
        """Run a CW ODMR experiment
        Arguments:  *
        """ #TODO: Fill in arg documentation

        # connect to the instrument server
        # connect to the data server and create a data set, or connect to an
        # existing one with the same name if it was created earlier.
        with InstrumentGateway() as gw, DataSource(datasetName) as cwODMRdata:

            # SRS actions
            gw.sg.set_mod_state(False)
            gw.sg.set_rf_amplitude(rfPower) # set ouput power

            # set up frequencies to sweep over
            # self.freqs = np.random.permutation(np.linspace(startFreq, endFreq, numFreqs, endpoint=True)) # permute the order of freqs to minimize time effects
            self.freqs = np.linspace(startFreq, endFreq, numFreqs, endpoint=True)
            # self.freqDict = dict([self.freqs])
            # for storing the experiment data
            self.mwCountsDict = dict( [[freq, [] ] for freq in self.freqs]) #should make dict of {freq: [counts1, counts2, ...], ...} pairs
            self.noMwCountsDict = dict( [[freq, [] ] for freq in self.freqs]) 

            # setup a relevant iterator depending on maxIterations
            if maxIterations < 0:
                iters = count() # infinite iterator
            else:
                iters = range(maxIterations) 

            # MAIN EXPERIMENT LOOP
            with NIDAQ() as mynidaq:

                # set laser power
                mynidaq.laser_power_atten(laser_power)

                for i in iters:

                    for freq in self.freqs:

                        # START READ TASK
                        mynidaq.start_read_task(2 * num_samples + 1)
            
                        # SET SRS OUTPUT FREQUENCY
                        gw.sg.set_frequency(freq)
                        print('freq : ', freq)
                       
                        # START PULSESTREAMER
                        gw.swabian.runSequenceInfinitely(Pulses(gw).CW_ODMR(clock_time, probe_time, wait_time))

                        # COUNTING WITH EXTERNAL TRIGGER
                        counts = obtain(mynidaq.read_samples(num_samples * 2 + 1 ))

                        # SEPARATING COUNTS INTO MW ON / MW OFF
                        # print('counts : ', counts)
                        self.mwCountsDict[freq].append(counts[0::2])
                        self.noMwCountsDict[freq].append(counts[1::2])
                        # print('counts[0::2] : ', counts[0::2])
                        # print('counts[1::2] : ', counts[1::2])

                        # SAVE CURRENT DATA TO DATA SERVER
                        cwODMRdata.push({'params': {'datasetName': datasetName, 
                                        'samplingFreq': samplingFreq, 
                                        'maxIterations': maxIterations,
                                        'startFreq': startFreq, 'endFreq': endFreq,
                                        'numFreqs': numFreqs, 'rfPower': rfPower,
                                        # 'preReadoutLaserAndMwTime': preReadoutLaserAndMwTime, 'laserAndMwReadOutTime': laserAndMwReadOutTime, 
                                        # 'extraLaserInitTime': extraLaserInitTime, 'waitTime': waitTime,
                                        'num_samples' : num_samples,
                                        'clock_time' : clock_time,
                                        'probe_time' : probe_time},
                                        'title': 'CW ODMR',
                                        'xlabel': 'Freq',
                                        'ylabel': 'Counts',
                                        'datasets': {'freqs': self.freqs, 
                                                'mwCountsDict': self.mwCountsDict, 
                                                'noMwCountsDict': self.noMwCountsDict
                                        }
                        })

                        # RESET SWABIAN OUTPUTS
                        gw.swabian.reset()

            print('Experiment finished!')


if __name__ == '__main__':

    exp = CW_ODMR_Measurement()
    print('''Running ODMR with 1Hz sampling rate, 10us MW+lasers before counts, 50us of counting, 20us of laser reinit, 
    30us wait with 1mW of laser power. Saving to CW_ODMR on dataserv''')
    exp.cwODMR(datasetName='CW_ODMR', sampleFreq=1, maxIterations=3,
               startFreq=2.75e9, endFreq=2.95e9, numFreqs=5, rfPower=-17,
               # preReadoutLaserAndMwTime=10000, laserAndMwReadOutTime=50000, 
               # extraLaserInitTime=20000, waitTime=30000, debug=False
               num_samples = 100000, clock_time = 11, probe_time = 50000)
    print('Completed cwODMR.py')