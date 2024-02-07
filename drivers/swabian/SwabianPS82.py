'''
    Driver for a Swabian PulseStreamer 8/2
    C.Egerstrom - Mar 2023
    De-lantz-ed version of Nazar's old driver. Use as you wish, be evil.

    N.Delegan - Sep 2020.
    Use as you wish, be nice.

    :to-do:
    De-jank init so multiple PSs doesn't break it

    leans heavily into API info at https://www.swabianinstruments.com/static/documentation/PulseStreamer/sections/api-doc.html
'''

from pulsestreamer import PulseStreamer

from rpyc.utils.classic import obtain #to deal with inevitable NetRef issues


class SwabianPulseStreamer82:
    """Driver for a Swabian Pulse Streamer 8/2"""

    def __init__(self): 
        self.voltage_sp_ch0 = 0


    def __enter__(self, address="192.168.1.76"): # "169.254.8.2" is the static fallback address
        '''Connects to pulsestreamer at specified address
        Arguments:  *address: PS Address. Default: 169.254.8.2'''
        self.address = address
        self.ps = PulseStreamer(self.address)
        print("Swabian Connected")
        return(self)

    def __exit__(self, exc_type, exc_value, exc_traceback):
        return

    def idn(self):
        '''Get device info from findPulseStreamers()'s return'''
        ser_num = self.ps.getSerial()
        frm_ver = self.ps.getFirmwareVersion()
        return "Serial #: " + ser_num + " // Firmware: " + frm_ver


    def reset(self):
        '''(In-built) Reset the Pulse Streamer device to the default state. 
        All outputs are set to 0V, and all functional configurations 
        are set to default. The automatic rearm functionality is enabled, 
        and the clock source is the internal clock of the device. 
        No specific trigger functionality is enabled, which means that 
        each sequence is streamed immediately when its upload is completed.'''
        self.ps.reset()


    def reboot(self):
        '''(In-built) Perform a soft reboot of the device without power-cycling.'''
        self.ps.reboot()


    def streaming_state(self, verbose=False):
        '''Get streaming status status
        Arguments:  *verbose. Default: True
        Returns:    *[If PS has a sequence, if it's streaming, if it's finishsed] 
                     as booleans if not Verbose, otherwise those are all in a string'''
        bool_seq = self.ps.hasSequence()
        bool_strm = self.ps.isStreaming()
        bool_fin = self.ps.hasFinished()
        if verbose:
            return('Sequence in memory: ' + str(bool_seq) + ' | Is streaming: ' + str(bool_strm)\
                + ' | Is finished: ' + str(bool_fin))
        return([bool_seq, bool_strm, bool_fin]) #if not verbose


    def reset_streamer(self):
        '''Sets all digital and analog outputs to 0'''
        self.ps.constant() #Calling the method without a parameter will result in the default output state with all output set to 0V.


    def test_sequence(self): #, unconnected=[0,5,6,7,'A1','A2']):
        """This is to run a test sequence that aims to test all of the channels. 
        i.e. infinite loop of 1 second TTL pulses on each of the counter channels.
        """       
        patt_d_ch0 = [(2e9,1),(1e9,1),(2e9,1),(1e9,1),(2e9,1),(1e9,1)] #AOM
        patt_d_ch1 = [(3e9,1),(1e9,0),(2e9,1),(1e9,0),(2e9,0),(4e9,0)] #Switch #1
        patt_d_ch2 = [(3e9,1),(1e9,0),(2e9,0),(1e9,0),(2e9,0),(4e9,0)] #Switch #2
        patt_d_ch3 = [(3e9,0),(1e9,0),(2e9,0),(1e9,0),(2e9,1),(4e9,0)] #Switch #2
        patt_d_ch4 = [(3e9,0),(1e9,0),(2e9,0),(1e9,0),(2e9,0),(4e9,0)] #RF
        #3 seconds H-H (counter 0), 2 seconds H-L (counter 1), 1 second L-H (counter 2)
        #followed by a 4 second break
        # 
        test_sequence = self.ps.createSequence() #create the sequence class
        test_sequence.setDigital(0,patt_d_ch0)
        test_sequence.setDigital(1,patt_d_ch1)
        test_sequence.setDigital(2,patt_d_ch2)
        test_sequence.setDigital(3,patt_d_ch3)
        test_sequence.setDigital(4,patt_d_ch4)
        '''test_sequence = self.ps.createSequence() #instantiate the sequence class
        for digI in range(8): 
            if digI not in unconnected: #Loop on for (1+digI) sec, off for (8-digI) sec
                patt_d_chI = [((digI+1)*1e9,1),((8-digI)*1e9,0)] 
                test_sequence.setDigital(digI,patt_d_chI)
        for anlgI in range(2):
            if ('A'+str(anlgI)) not in unconnected:
                test_sequence.setAnalog(anlgI, [(100, 0.5**(anlgI+1)), (100, 0), (100, -1*0.5**(anlgI+1), (100, 0))]) 
                #Loop 100ns (+0.5 for A0/+0.25 for A1), 100ns 0, 100ns (-0.5 for A0/-0.25 for A1), 100ns 0
'''
        n_runs = self.ps.REPEAT_INFINITELY #inifnite number of runs
        self.ps.stream(test_sequence, n_runs)
    

    def runSequenceInfinitely(self, seq):
        '''Main workhorse function when using Swab thru an InstrumentGateway. Obtains the desired sequence and starts streaming it'''
        self.ps.stream(obtain(seq), self.ps.REPEAT_INFINITELY )


if __name__ == '__main__':
    print("Testing the PS82 pulse streamer driver...\n")
    with SwabianPulseStreamer82() as test_ps: 
        print("Asking the Swabian PS82 to self identify:\n", test_ps.deviceInfo())
        print("\nIdentification successful. Driver seems to work.\n")       
