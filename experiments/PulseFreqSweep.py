import time

import numpy as np
from nspyre import DataSource
from nspyre import InstrumentGateway

SWABIAN_OUTPUT_PORT = 0

class PulseFreqSweepMeasurement:

    def PulseFreqSweep(self, datasetName: str, startFreq: float, 
                        endFreq: float, numFreqs: int, time_per_freq: float):

        with InstrumentGateway() as gw, DataSource(datasetName) as data:

            #Pulse frequencies to sweep over
            self.freqs = np.linspace(startFreq,endFreq,numFreqs, endpoint= True)
        
            #Our data with cps at each freq
            # self.DAQ_cps_dict = dict([[freq,[]] for freq in self.freqs]) 

            #print(f'freqs{self.freqs}')
            sequence = gw.swabian.ps.createSequence()
            for freq in self.freqs:
                #print(f'freq:{freq}')
                #Generates and starts the new pulse sequence
                pulse_pattern = self.generate_pulse_pattern(freq)
                sequence.setDigital(SWABIAN_OUTPUT_PORT, pulse_pattern)
                gw.swabian.runSequenceInfinitely(sequence)

                time.sleep(time_per_freq)

    #TODO : test this
    # Generates a pulse sequence with 50% duty cycle of frequency pulse_freq
    def generate_pulse_pattern(self, pulse_freq:float):
        #Default units of pulse streamer is ns
        period = 1/pulse_freq
        converted_period = (10**9)*period

        #convert to int?
        half_period = converted_period/2
        return [(half_period, 1), (half_period, 0)]

if __name__ == '__main__':
    exp = PulseFreqSweepMeasurement()
    print('''Running test Pulse Freq Sweep''')
    exp.PulseFreqSweep("test",1000,10000,10,30)