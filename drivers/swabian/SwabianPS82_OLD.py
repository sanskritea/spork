"""
    lantz.drivers.swabian.SwabianPS82.SwabianPulseStreamer82
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    :versions/changelog: 
    * V0.1 - Sep 2020 -  It is birthed, initiates, and can run a dummy sequence. 
    # TODO: update for timing discrepencies. 
    # TODO: you should probably instantiate the class you defined in the name==main section, b/c you're not right now - CE 11/16/22

    :copyright: None; Nazar wrote this, use it as you wish, be nice.
"""

from lantz import Action, Feat, Driver
from pulsestreamer import PulseStreamer, Sequence

class SwabianPulseStreamer82(Driver):
    """Driver for a Swabian Pulse Streamer 8/2.
    """
    def __init__(self, address= "169.254.8.2"): #TODO: get our IP, not backup
        # super().__init__()
        self.address = address
        self.ps = PulseStreamer(address)
        self.voltage_sp_ch0 = 0

    def finalize(self):
        #Reset the device
        self.ps.reset()
        None

    
    @Feat(read_once=True)
    def idn(self):
        ser_num = self.ps.getSerial()
        frm_ver = self.ps.getFirmwareVersion()
        return "Serial #: " + ser_num + " // Firmware: " + frm_ver

    @Feat(read_once=False)
    def streaming_state(self):
        """Get streaming status status
        """
        bool_seq = str(self.ps.hasSequence())
        bool_strm = str(self.ps.isStreaming())
        bool_fin = str(self.ps.hasFinished())
        return "Sequence in memory ("+bool_seq+") | "+"Is streaming ("+bool_strm+") | "+"Finished ("+bool_fin+")"

    @Action()
    def reset_streamer(self):
        """Reset
        """
        self.voltage_sp_ch0 = 0
        self.ps.constant(([],self.voltage_sp_ch0,0))
        self.ps.reset()

    @Action()
    def cntr_test_sequence(self):
        """This is to run a test sequence that aims to test all of the channels. 
        i.e. infinite loop of 1 second TTL pulses on each of the counter channels.
        """
        patt_d_ch0 = [(3e9,0),(1e9,0),(2e9,0),(1e9,0),(2e9,0),(4e9,0)] #EMPTY
        patt_d_ch1 = [(3e9,1),(1e9,0),(2e9,1),(1e9,0),(2e9,0),(4e9,0)] #Switch #1
        patt_d_ch2 = [(3e9,1),(1e9,0),(2e9,0),(1e9,0),(2e9,0),(4e9,0)] #Switch #2
        patt_d_ch3 = [(3e9,0),(1e9,0),(2e9,0),(1e9,0),(2e9,1),(4e9,0)] #Switch #2
        patt_d_ch4 = [(3e9,0),(1e9,0),(2e9,0),(1e9,0),(2e9,0),(4e9,0)] #RF
        #3 seconds H-H (counter 0), 2 seconds H-L (counter 1), 1 second L-H (counter 2)
        #followed by a 4 second break

        test_sequence = self.ps.createSequence() #create the sequence class
        test_sequence.setDigital(0,patt_d_ch0)
        test_sequence.setDigital(1,patt_d_ch1)
        test_sequence.setDigital(2,patt_d_ch2)
        test_sequence.setDigital(3,patt_d_ch3)
        test_sequence.setDigital(4,patt_d_ch4)

        n_runs = self.ps.REPEAT_INFINITELY #inifnite number of runs
        self.ps.stream(test_sequence, n_runs)

    @Action()
    def aom_test_sequence(self):   
        self.ps.constant(([],self.voltage_sp_ch0,0))

    @Feat(units="V",limits=(0, 1))
    def ach0_voltage_sp(self):
        """To handle output power set point (mW) in constant power mode
        """
        return self.voltage_sp_ch0
        # return 1000 * float(self.query("p?"))
        #Power setter, between 0 and 100 mW
    @ach0_voltage_sp.setter
    def ach0_voltage_sp(self, value):
        self.voltage_sp_ch0 = value
        # self.query("p {:.5f}".format(value / 1000))

if __name__ == '__main__':
    print("Testing the PS82 pulse streamer driver...\n")
    test_ps = PulseStreamer("169.254.8.2")
    print("Asking the Swabian PS82 to self identify:\n Serial #: {0} || Firmware: {1}".format(test_ps.getSerial(),test_ps.getFirmwareVersion()))
    print("\nIdentification successfull. Driver seems to work.\n")        