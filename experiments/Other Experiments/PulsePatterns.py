

from nspyre import InstrumentGateway
from rpyc.utils.classic import obtain

class Pulses():
    #0 connected to AOM, 4 connected to Vaunix switch, Switch1=H => Switch 2, which will be for main counting. Switch1=L => Switch 3, 
    #where Switch3=H will do normalization when both channels from S1 are in use for diff measurements, and Switch1&3=L => 50ohm dump
    # DIG_CHAN_DICT = {'aomSwitch': 0, 'rfSwitch': 4, 'Switch1': 1, 'Switch2': 2, 'Switch3': 3} 

    #COUNTS CHAN0 = PFI11 (always counting, no need to turn on), CHAN1 = PFI1, CHAN2=PFI4, CAN3=PFI8
    # SWITCHES_FOR_COUNT_CHAN = { 0:[], 1: [DIG_CHAN_DICT['Switch1'], DIG_CHAN_DICT['Switch2']], 2:[DIG_CHAN_DICT['Switch1']], 3:[DIG_CHAN_DICT['Switch3']] }
    #TODO: Measure AOM Lag and Switch Lag times, then incorporate those as class vars for all of these pulses?!? <= Maybe the only reason to wrap these as a class 
    
    # SANS: SWABIAN CHANNEL DICT FOR NO-RF-SWITCH OPERATION
    DIG_CHAN_DICT = {'DAQ': 0, 'AOM': 1, 'SRS': 2} 

    #General sequence is of the form:
    #[(duration in ns, [list of channels to turn on, others are off (empty is all off)], A0_Voltage, A1_Voltage), (duration2, [channel2], A0V_2, A1V_2), ...]

    def __init__(self):
        pass

    def rawMwPulse(self, mwTime):
        '''Returns a Swabian Seq-compatible list that will switch MWs on for mwTime
        Arguments:  *mwTime, time to turn MWs on for (in ns)'''
        return( [(mwTime, [self.DIG_CHAN_DICT['DAQ']], 0,0)] )

    def rawLaserPulse(self, laserTime):
        '''Returns a Swabian Seq-compatible list that will switch laser on for laserTime
        Arguments:  *laserTime, time to turn laser on for (in ns)'''
        return( [(laserTime, [self.DIG_CHAN_DICT['AOM']], 0,0)] )

    def wait(self, delayTime):
        '''Returns a Swabian Seq-compatible list with all channels off for delayTime
        Arguments:  *delayTime, time to do nothing (in ns)'''
        return( [(delayTime, [], 0,0)] )

    def rawReadOutPulse(self, readOutTime, countChan):
        '''Returns a Swabian Seq-compatible list with laser on and specified counter channel on for readOutTime
        Arguments:  *readOutTime, time to have laser on and counter channel open for (in ns)
                    *countChan, one of 0-3 corresponding to PFI11,1,4,8 respectively'''
        return countChan
            

    def readoutAndInitPulse(self, aomRiseLag, readOutTime, extraInitLaserTime, aomFallLagAndISC, countChan):
        '''Returns a Swabian Seq-compatible list with laser on for an aomRiseLag time, then a readout pulse with specified counterChan
        for readOutTime, then just laser on for extraInitLaserTime before waiting for aomFallLagAndISC time with everything off
        Arguments:  *readOutTime, time to have laser on and counter channel open for (in ns)
                    *extraInitLaserTime, time to have laser on to reset NV state (in ns)
                    *aomRiseLag, time to have laser on before counting corresponding to how long the AOM takes to rise (in ns)
                    *aomFallLagAndISC, time to wait for AOM to fall and for ISC to depopulate (in ns)
                    *countChan, one of 0-3 corresponding to PFI11,1,4,8 respectively'''

        aomLagCompensation = self.rawLaserPulse(aomRiseLag)
        readOutWindow = self.rawReadOutPulse(readOutTime, countChan)
        extraLaserForInit = self.rawLaserPulse(extraInitLaserTime)
        aomLagAndIscDelay = self.wait(aomFallLagAndISC)
        return(aomLagCompensation + readOutWindow + extraLaserForInit + aomLagAndIscDelay)
    
    def readOutWithMWs(self, preReadoutLaserAndMwTime, readOutTime, extraInitLaserTime, aomFallLagAndISC, countChan):
        '''Returns a Swabian Seq-compatible list with laser+MWs on for an preReadoutLaserAndMwTime, then a readout pulse with specified counterChan
        for readOutTime with laser+MWs, then just laser on for extraInitLaserTime before waiting for aomFallLagAndISC time with everything off
        Arguments:  *preReadoutLaserAndMwTime, time to have laser+MWs on before counting (should at least be long enough for AOM to rise) (in ns)
                    *readOutTime, time to have laser+MWs on and counter channel open for (in ns)
                    *extraInitLaserTime, time to have only laser on to reset NV state (in ns)
                    *aomFallLagAndISC, time to wait for AOM to fall and for ISC to depopulate (extra delay can be added here) (in ns)
                    *countChan, one of 0-3 corresponding to PFI11,1,4,8 respectively'''
        preReadoutLaserAndMws = [(preReadoutLaserAndMwTime, [self.DIG_CHAN_DICT['aomSwitch'],self.DIG_CHAN_DICT['rfSwitch']], 0,0)]
        readOutWindow = [(readOutTime, [self.DIG_CHAN_DICT['aomSwitch'],self.DIG_CHAN_DICT['rfSwitch']]+self.SWITCHES_FOR_COUNT_CHAN[countChan], 0,0)] 
        extraLaserForInit = self.rawLaserPulse(extraInitLaserTime)
        aomLagAndIscDelay = self.wait(aomFallLagAndISC)
        return(preReadoutLaserAndMws + readOutWindow + extraLaserForInit + aomLagAndIscDelay)

    def balancedCwODMR(self, preReadoutLaserAndMwTime, readOutTime, extraInitLaserTime, aomFallLagAndISC):
        '''Returns a Swabian Seq-compatible list with laser+MWs (no counting) for a preReadout time, then readouts with CW laser+MWs to Chan1
        before running just the laser to re-init the NV state to 0 before waiting some time with everything off for the AOM to fall and ISC to depopulate
        for readOutTime with laser+MWs, then just laser on for extraInitLaserTime before waiting for aomFallLagAndISC time with everything off.
        This sequence with idential timings is repeated with no MWs and counts going into Chan2
        Arguments:  *preReadoutLaserAndMwTime, time to have laser+MWs on before counting (should at least be long enough for AOM to rise) (in ns)
                    *readOutTime, time to have laser+MWs on and counter channel open for (in ns)
                    *extraInitLaserTime, time to have only laser on to reset NV state (in ns)
                    *aomFallLagAndISC, time to wait for AOM to fall and for ISC to depopulate (extra delay can be added here) (in ns)
                    *countChan, one of 0-3 corresponding to PFI11,1,4,8 respectively'''
        withMwsSeq = self.readOutWithMWs(preReadoutLaserAndMwTime, readOutTime, extraInitLaserTime, aomFallLagAndISC, countChan=1) #want to use countChans 1&2 thru same switch for differential measurement
        withoutMwsSeq = self.readoutAndInitPulse(preReadoutLaserAndMwTime, readOutTime, extraInitLaserTime, aomFallLagAndISC, countChan=2)
        return(withMwsSeq + withoutMwsSeq)
    
    def aomPlusSwitchDelayCal(self, laserStartTime, laserDur, countStartTime, countDur, countChan, fullLengthWithBuffer):
        '''Returns a Swabian Seq(!) with the AOM turned on at a specified laserStartTime for laserDur, as well as a specified counter opened
        at specified countStartTime for a countDur, then the seq repeats every fullLengthWithBuffer (which must be a multiple of 8ns for the Swab to loop it optimally)
        Arguments:  *laserStartTime, time to turn on the laser (in ns)
                    *laserDur, how long to turn the laser on for (in ns)
                    *countStartTime, time to turn on the desired counting channel (in ns)
                    *countDur, how long to turn the counting channel on for (in ns)
                    *countChan, which counter channel to use (one of 0-3 corresponding to PFI11,1,4,8 respectively)
                    *fullLengthWithBuffer, total time of sequence (in ns), preferably a multiple of 8ns if trying to rough sync the Swab+DAQ for AllCounts measurements'''
        with InstrumentGateway() as gw: #creating sequence this way rather than with a pattern with all channels to simplify logic
            seq = obtain(gw.swabian.ps.createSequence()) #trying to pre-emptively avoid rpyc netref issues

            seq.setDigital(self.DIG_CHAN_DICT['aomSwitch'], [(laserStartTime, 0), (laserDur, 1), (fullLengthWithBuffer-(laserStartTime+laserDur), 0)] ) #switch the laser on only for the laserDur
            for switch in self.SWITCHES_FOR_COUNT_CHAN[countChan]: #make sure relevant switches for each counts channel go on together (ignores diffs in switch times, but should be fine)
                seq.setDigital(switch, [(countStartTime, 0), (countDur, 1), (fullLengthWithBuffer-(countStartTime+countDur), 0)] )
            


