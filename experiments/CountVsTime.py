"""
This is a basic CountVsTime Application

Copyright (c) April 2023, C. Egerstrom
All rights reserved.

This work is licensed under the terms of the 3-Clause BSD license.
For a copy, see <https://opensource.org/licenses/BSD-3-Clause>.
"""

import time
from itertools import count

import numpy as np

from nspyre import DataSource
from nspyre import InstrumentGateway
from nspyre import StreamingList

from rpyc.utils.classic import obtain

from guis.guiElements_general import flexSave
from experiments.NewPulses import Pulses
from drivers.ni.nidaq_final import NIDAQ


class CountVsTimeMeasurement:

    def CountVsTime(self, datasetName: str, sampleFreq: float, laser_power: float, maxIterations: int, autosaveParams=None, debug=False):
        """Run a CountVsTime2 experiment

        Args:
            datasetName: name of the dataset to push data to
            sampleFreq (float): how quickly to read data (in Hz)
            ctrChanNums: Which PFI channels to read. Default is [11,1,4,8]
            autosaveParams: Default: None, but will take a list of [shouldAutosave, autosaveInterval] 
            debug: optional (default False), will run TimeVsTime if true
        """

        # connect to the instrument server
        # connect to the data server and create a data set, or connect to an
        # existing one with the same name if it was created earlier.
        with InstrumentGateway() as gw, DataSource(datasetName) as cvt_data:

            # turn on laser
            # gw.swabian.runSequenceInfinitely(Pulses(gw).laser_on())

            # fake input signal from swabian
            # gw.swabian.runSequenceInfinitely(Pulses(gw).fakeDAQinput())

            # for storing the experiment data
            self.times = np.zeros(0)
            self.counts = np.zeros(0)

            # setup a relevant iterator depending on maxIterations
            if maxIterations < 0:
                iters = count() #infinite iterator
            else:
                iters = range(maxIterations)


            with NIDAQ() as mynidaq:

                # set laser power
                mynidaq.laser_power_atten(laser_power)

                for i in iters:

                    self.startTime = time.time() #take time diff measurements

                    if hasattr(autosaveParams, '__len__'): # Unpack the autosaveParams if they exist
                        shouldAutosave, autosaveInterval = autosaveParams
                    else:
                        shouldAutosave = False

                    # creating streaming list for counts
                    PL_data_streaming_list = StreamingList([])

                    # start Swabian external trigger for counting
                    gw.swabian.runSequenceInfinitely(Pulses(gw).CountVsTime(sampleFreq))

                    # if debug: # time vs time
                    #     self.counts = np.append(self.counts, time.time())
                    #     time.sleep(1 / sampleFreq) #need to set a delay somehow since 'read' is instant

                    # else: # start DAQ counting
                    
                    while True:

                        # START READ TASK
                        mynidaq.start_read_task(2)

                        newData = obtain(mynidaq.read_samples(2))
                        self.counts = [np.append(self.counts, newData)]
                        print('counts : ', self.counts)

                        # Read and save new time
                        self.times = np.append(self.times, time.time() - self.startTime)

                        # Using streaming lists
                        # PL_data_streaming_list.append(dict([(f'PFI3 counts', self.counts[j]) if j<1 else ('times', self.times) for j in range(2)]))

                        # save the current data to the data server.
                        cvt_data.push({'params': {'SampleFreq': sampleFreq},
                                        'title': 'CountVsTime',
                                        'xlabel': 'Time (s)',
                                        'ylabel': 'Counts',
                                        # 'datasets': {'times': self.times, 'counts': self.counts}
                                        'datasets': dict([('CountVsTime', self.counts[j]) if j<1 else ('times', self.times) for j in range(2)]) # HOW TO UNDERSTAND THIS LINE?
                                        # 'datasets': PL_data_streaming_list
                        })

                        if shouldAutosave and (i+1)%autosaveInterval == 0: # Autosave logic, +1 so it doesn't autosave first data point
                            flexSave(datasetName, 'CountVsTime', 'autosave')

            flexSave(datasetName, 'CountVsTime', 'final') # after measurement finishes


if __name__ == '__main__':
    exp = CountVsTimeMeasurement()
    print('Running CountVsTime with 1Hz sampling rate for max 1hr, saving to CountVsTime on dataserv')
    exp.CountVsTime('CountVsTime', 1, 3600)

