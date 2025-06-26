"""
This is an application to measure the AOM+switch lag on Jasper

Copyright (c) May 2023, C. Egerstrom
All rights reserved.

This work is licensed under the terms of the 3-Clause BSD license.
For a copy, see <https://opensource.org/licenses/BSD-3-Clause>.
"""

import time
from itertools import count

import numpy as np

from nspyre import DataSource
from nspyre import InstrumentGateway

from rpyc.utils.classic import obtain

from experiments.PulsePatterns import Pulses


class AOMandSwitchLag_Measurement:

    def aomAndSwitchLag(self, datasetName: str, sampleFreq: float, maxIterations: int,
               laserStartTime: int, laserDur: int,
               firstCountStartTime: int, lastCountStartTime: int, numCountsStartTimes: int, 
               countDur: int, countChan:int, 
               totalSeqLength:int, laserPower = 1, 
               turnAomOffAtEnd = True, debug=False):
        """Run an AOM + switch delay calibration experiment
        Arguments:  *
        """ #TODO: Fill in arg documentation

        # connect to the instrument server
        # connect to the data server and create a data set, or connect to an
        # existing one with the same name if it was created earlier.
        with InstrumentGateway() as gw, DataSource(datasetName) as aomAndSwitchLagData:

            #Make sure laser is on first
            #This will hang and annoy the user until the laser interlock is flipped
            if not gw.laser.is_on():
                gw.laser.turn_on()
            while not gw.laser.is_on():
                print('Waiting for laser interlock')
                time.sleep(2)
            gw.laser.set_power(laserPower)

            # for storing the experiment data
            self.countStartTimes = np.linspace(firstCountStartTime, lastCountStartTime, numCountsStartTimes, endpoint=True)
            self.countsAtStartTimes = [np.zeros(0) for i in self.countStartTimes]

            #setup a relevant iterator depending on maxIterations
            if maxIterations < 0:
                iters = count() #infinite iterator
            else:
                iters = range(maxIterations) 
            
            if totalSeqLength%8 != 0:
                print('''Total length of sequence is not a multiple of 8ns. 
                      Could lead to desyncing if trying to clock Swab+DAQ together
                      (but that isn't implemented yet, so who cares)''') #slight warning, but not actually an issue, at least for now

            #MAIN EXPERIMENT LOOP
            for i in iters:
                for j,countStartTime in enumerate(self.countStartTimes):

                    gw.swabian.ps.stream(Pulses().aomPlusSwitchDelayCal(laserStartTime, laserDur, countStartTime, countDur, countChan, totalSeqLength), #TODO: Figure out why you made this a class, you clown
                                         gw.swabian.ps.REPEAT_INFINITELY)
                    countsForStartTime = obtain(gw.nidaq.readCtrs_singleRead_intClk(acqRate=sampleFreq, ctrChanNums=[countChan]))
                    
                    self.countsAtStartTimes[j] = np.append(self.countsAtStartTimes[j], countsForStartTime)

                    # save the current data to the data server.
                    aomAndSwitchLagData.push({'params': {'datasetName': datasetName, 'sampleFreq': sampleFreq, 'maxIterations': maxIterations,
                                                         'laserStartTime': laserStartTime, 'laserDur': laserDur, 
                                                         'countStartRange': (firstCountStartTime, lastCountStartTime), 'numCountsStartTimes': numCountsStartTimes, 
                                                         'countDur': countDur, 'countChan':countChan, 'totalSeqLength': totalSeqLength, 'laserPower': laserPower},
                                    'title': 'AOM and Switch Lag',
                                    'xlabel': 'AOM Wait',
                                    'ylabel': 'Counts',
                                    'datasets': {'countStartTimes': self.countStartTimes, 'countsAtStartTimes': self.countsAtStartTimes}
                    })
            
            if turnAomOffAtEnd:
                gw.laser.turn_off()



if __name__ == '__main__':
    pass
    """exp = AOMandSwitchLag_Measurement()
    print('''Running ODMR with 1Hz sampling rate, 10us MW+lasers before counts, 50us of counting, 20us of laser reinit, 
    30us wait with 1mW of laser power. Saving to CW_ODMR on dataserv''')
    exp.cwODMR(datasetName='CW_ODMR', sampleFreq=1, maxIterations=10,
               startFreq=2.75e9, endFreq=2.95e9, numFreqs=20, rfPower=-17,
               preReadoutLaserAndMwTime=10000, laserAndMwReadOutTime=50000, 
               extraLaserInitTime=20000, waitTime=30000, 
               laserPower = 1, turnAomOffAtEnd = True, debug=False)"""

