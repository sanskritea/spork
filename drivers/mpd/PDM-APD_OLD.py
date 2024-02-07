# -*- coding: utf-8 -*-
'''
    lantz.drivers.thorlabs.LEDD1B
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    :copyright: 
    N.Delegan - Sep 2020.
    Use as you wish, be nice.

    :description:
    Driver for MPD PDM series APD
    
    :versions/changelog: 
    * V0 - Oct 2020 -  Driver started
    * V1 - Nov 2020 - Driver tested and working.

    :dependencies: 
    * Needs a properly wired DAQ to gate the APD. Note that this should not be used for fast gating.
    This is limited to turning coutns on and off.
	
	:usage:
    APD gate is wired into BOX2 AO0 so, should be DAQ2/ao1

    For practical usage, just turn APD ON/OFF
  
'''
# Import stuff
from lantz import Action, Feat, Driver
import numpy as np
import nidaqmx as daq

# Driver Class
class MPD_PDMSeries_APD(Driver):
    '''Driver class for the APD.
    In essence this is a nidaqmx wrapper as the gating for the APD.
    '''
    # Constructor function to define how things are hardwired
    def __init__(self):
        self.APDGate_port = 'DAQ1/ao3'
        self.DAQ_name = 'DAQ1'
        self.gate_voltage = 0 #By default gate to 0

        # Initialize function to get the DAQ_DEVICE OBJECT
        # Identify the local machine
        this_system = daq.system.System.local()
        # We know on our current system that our DAQ is named 'DAQ1'
        self.DAQ_device = this_system.devices[self.DAQ_name]
        None

    # Nothing is needed here
    def finalize(self):
        self.gate_sp = 'Counting OFF'

    @Feat(values={'Counting ON' : 5, 'Counting OFF' : 0})
    def gate_sp(self):
        return self.gate_voltage
    @gate_sp.setter
    def gate_sp(self, value):
        self.gate_voltage = value
        return self.set_output_potential(self.gate_voltage)

    # Main function to be called
    def set_output_potential(self,voltage):
        '''Function that converts provided potential into DAC 0-5V (OFF/ON).'''
        self.gate_voltage = voltage
        with daq.Task() as task:
            task.ao_channels.add_ao_voltage_chan(self.APDGate_port, min_val = 0, max_val = 5) #defautls are -10 to 10
            task.write(self.gate_voltage, auto_start=True)
        None

if __name__ == '__main__':
    print('Testing the LEDD1B driver...\n')
    with MPD_PDMSeries_APD() as APD:
        print('Current APD gate output potential is {0} V.\n'.format(APD.gate_voltage))
        print([ao.name for ao in APD.DAQ_device.ao_physical_chans])
        print('\nAPD-DAQ/ao channels verified. Driver seems to work.\n')