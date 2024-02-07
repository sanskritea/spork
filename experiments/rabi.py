# -*- coding: utf-8 -*-



#import the fancy doodads
import experiments.PulsePatterns as pulse
import drivers.swabian.SwabianPS82 #pulse confetti
import time 
import numpy as np
from itertools import count
from nspyre import InstrumentGateway, DataSource

from rpyc.utils.classic import obtain


#laser vars 
init = 100 #ns
aomRise = None
ISC = None
initBuffer = 450 #us


readout = 150 #us
readoutBuffer = 2 #us

#MW vars
power = None

#channelList = [laser: '0', rf: '4']

class Rabi_Measurement:

    def rabi(self, minMwTime: int, maxMwTime: int, mwTime: int, numMwTimes: int, freq: float, delayTime: int, aomLag: int, readoutTime: int, 
                initTime: int, aomFallLagAndISC: int, countChan, rfPower: float,
                lasetTime: int, maxIters: int, laserPower = 1):
        ''' readout = pulse.readoutAndInitPulse(aomRiseLag, readOutTime, extraInitLaserTime, aomFallLagAndISC, countChan)
        mwSeq = pulse.rawMwpulse(mwTime)
        wait = pulse.wait(delayTime)
        return(readout + mwSeq + readout + wait)'''

        with InstrumentGateway() as gw, DataSource(datasetName) as rabiData:

            if not gw.laser.is_on():
                gw.laser.turn_on()
            while not gw.laser.is_on():
                print('Waiting for laser interlock.')
                time.sleep(2)
            gw.laser.set_power(laserPower)

            #turn on sig-gen
            gw.vaunix.setFreq(freq)
            gw.vaunix.setPwrLvl(rfPower)
            gw.vaunix.setIfRFOut(True)



            if maxIters < 0:
                iters = count()
            else:
                iters = range(maxIters)

            self.mwTimes = np.linspace(minMwTime, maxMwTime, numMwTimes, endpoint=True)

            #data storage
            self.mwCountsDict = dict([[mwTime, [] ] for mwTime in self.mwTimes])
            self.noMwCountsDict = dict( [[mwTime, [] ] for mwTime in self.mwTimes]) 

            #MAIN EXPERIMENT
            for i in iters:
                for mwTime in self.mwTimes:
                    readout = pulse.readoutAndInitPulse(aomLag, readoutTime, initTime, aomFallLagAndISC, countChan)
                    mwShort = pulse.rawMwpulse(mwTime)
                    mwLong = pulse.rawMwpulse(maxMwTime - mwTime)
                    wait = pulse.wait(delayTime)
                    background = pulse.wait(maxMwTime/2) #add option to disable/enable


                    gw.swabian.ps.stream(readout + mwShort + wait +
                                        readout + mwLong + wait + 
                                        readout + background + wait)
                
                    mwCounts, noMwCounts = obtain(gw.nidaq.readCtrs_singlRead_intClk(acqRate= None, ctrChanNumns=[(1,4)]))
                    self.noMwCountsDict[mwTime].append(noMwCounts)

                    rabiData.push({'params': {'datasetName' : rabiData,
                                            'counts' : noMwCounts,
                                            'rftime': maxMwTime, 
                                            'readoutTime': readout,
                                            'delayTime': wait,
                                            'mwTime': mwTime
                                            },
                                            'title': 'Rabi',
                                            'xlabel': 'Time (s)',
                                            'ylabel': 'PL (cts/s)',
                                            'datasets': {'mwTimes': self.mwTimes, 'noMwCountsDict': self.noMwCountsDict}


                })
                    
if __name__ == '__main__':
    exp = Rabi_Measurement()
    print('''Running ODMR with 1Hz sampling rate, 10us MW+lasers before counts, 50us of counting, 20us of laser reinit, 
    30us wait with 1mW of laser power. Saving to CW_ODMR on dataserv''')
    exp.cwODMR(datasetName='Rabi', sampleFreq=1, maxIterations=10,
                startFreq=2.75e9, endFreq=2.95e9, numFreqs=20, rfPower=-17,
                preReadoutLaserAndMwTime=10000, laserAndMwReadOutTime=50000, 
                extraLaserInitTime=20000, waitTime=30000, 
                laserPower = 1, turnLaserOffAtEnd = True, debug=False)