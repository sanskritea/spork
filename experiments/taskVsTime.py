"""
This is a basic TaskVsTime Application

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



class TaskVsTimeMeasurement:

    def taskVsTime(self, datasetName: str, sampleFreq: float, maxIterations: int, ctrChanNums=[0,1,2,3], autosaveParams=None, debug=False):
        """Run a TaskVsTime experiment

        Args:
            datasetName: name of the dataset to push data to
            sampleFreq (float): how quickly to read data (in Hz)
            maxIterations: max number of data points to collect. If negative, will go infinitely 
            ctrChanNums: Which PFI channels to read. Default is [11,1,4,8]
            autosaveParams: Default: None, but will take a list of [shouldAutosave, autosaveInterval] 
            debug: optional (default False), will run TimeVsTime if true
        """

        # connect to the instrument server
        # connect to the data server and create a data set, or connect to an
        # existing one with the same name if it was created earlier.
        with InstrumentGateway() as gw, DataSource(datasetName) as tvt_data:

            # for storing the experiment data
            self.times = np.zeros(0)
            self.counts = [np.zeros(0) for i in ctrChanNums]

            self.startTime = time.time() #take time diff measurements

            #setup a relevant iterator depending on maxIterations
            if maxIterations < 0:
                iters = count() #infinite iterator
            else:
                iters = range(maxIterations)

            if hasattr(autosaveParams, '__len__'): #Unpack the autosaveParams if they exist
                shouldAutosave, autosaveInterval = autosaveParams
            else:
                shouldAutosave = False

            # creating streaming list for counts
            PL_data_streaming_list = StreamingList([])

            #MAIN EXPERIMENT LOOP
            for i in iters:

                if debug: #time vs time
                    self.counts = np.append(self.counts, time.time())
                    time.sleep(1/sampleFreq) #need to set a delay somehow since 'read' is instant

                else: #normal operating mode
                    # newData = obtain(gw.nidaq.readCtrs_singleRead_intClk(acqRate=sampleFreq, ctrChanNums=ctrChanNums)) #obtaining to prevent NetRef issues
                    newData = obtain(gw.nidaq.readCtrs_singleRead_intClk(acqRate=sampleFreq, ctrChanNums=ctrChanNums))
                    #print(newData) #DEBUG
                    self.counts = [np.append(self.counts[j], newData[j]) for j in range(len(ctrChanNums))]
                
                #Read and save new time
                self.times = np.append(self.times, time.time()-self.startTime)

                # Using streaming lists
                PL_data_streaming_list.append(dict([(f'PFI{ctrChanNums[j]}counts', self.counts[j]) if j<len(ctrChanNums) else ('times', self.times) for j in range(len(ctrChanNums)+1)]))

                # save the current data to the data server.
                tvt_data.push({'params': {'SampleFreq': sampleFreq, 'MaxIters': maxIterations, 'CtrChanNums': ctrChanNums},
                                'title': 'Task vs Time',
                                'xlabel': 'Time (s)',
                                'ylabel': 'Counts',
                                #'datasets': {'times': self.times, 'counts': self.counts}
                                # 'datasets': dict([(f'PFI{ctrChanNums[j]}counts', self.counts[j]) if j<len(ctrChanNums) else ('times', self.times) for j in range(len(ctrChanNums)+1)]) # HOW TO UNDERSTAND THIS LINE?
                                'datasets': PL_data_streaming_list
                })

                if shouldAutosave and (i+1)%autosaveInterval == 0: #Autosave logic, +1 so it doesn't autosave first data point
                    flexSave(datasetName, 'TvT', 'autosave')

        flexSave(datasetName, 'TvT', 'final') #after measurement finishes



if __name__ == '__main__':
    exp = TaskVsTimeMeasurement()
    print('Running TaskVsTime with 1Hz sampling rate for max 1hr, saving to TaskVsTime on dataserv')
    exp.taskVsTime('TaskVsTime', 1, 3600)

