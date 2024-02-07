# -*- coding: utf-8 -*-
'''
    lantz.drivers.thorlabs.LEDD1B
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    :copyright: 
    N.Delegan - Sep 2020.
    Use as you wish, be nice.

    :description:
    Driver for LEDD1B LED driver
    
    :versions/changelog: 
    * V0 - Oct 2020 -  Driver started
    * V1 - Nov 2020 - Driver tested and working.

    :dependencies: 
    * Needs a properly wired DAQ to modulate the WLS.
	
	:usage:
    WLS is wired into BOX2 AO0 so, should be DAQ2/ao0
    1. Create a Task and Virtual Channels
    2. Configure the Timing Parameters
    3. Start the Task
    4. Perform a Read operation from the DAQ
    5. Perform a Write operation to the DAQ
    6. Stop and Clear the Task.

    For practical usage, just set the output potential you want to send your LED.
    golden_light function was created as a 'default' light mode.
  
'''
# Import stuff
from lantz import Action, Feat, Driver
import numpy as np
import nidaqmx as daq

# Driver Class
class WhiteLightSourceLEDD1B(Driver):
    '''Driver class for the LEDD1B.
    In essence this is a nidaqmx wrapper as the modulation will be controlled through the DAC
    '''
    # Constructor function to define how things are hardwired
    def __init__(self):
        self.WLS_port = 'Dev1/ao2'
        self.DAQ_name = 'Dev1'
        self.golden_V = 0.225 #goldylocks power setting for faster action
        self.DAQ_voltage_output = 0

    # Initialize function to get the DAq_DEVICE OBJECT
    def initialize(self):
        # Identify the local machine
        this_system = daq.system.System.local()
        # We know on our current system that our DAQ is named 'DAQ1'
        self.DAQ_device = this_system.devices[self.DAQ_name]
        None

    # Nothing is needed here
    def finalize(self):
        self.power_sp = 0

    @Feat(units='V',limits=(0, 5))
    def power_sp(self):
        return self.DAQ_voltage_output
    @power_sp.setter
    def power_sp(self, value):
        self.DAQ_voltage_output = value
        return self.set_output_potential(self.DAQ_voltage_output)

    # Actions
    @Action()
    def max_light(self):
        '''Max lamp light
        '''
        self.power_sp = 5
        
    @Action()
    def off_light(self):
        '''Max lamp light
        '''
        self.power_sp = 0

    @Action()
    def golden_light(self):
        '''Max lamp light
        '''
        self.power_sp = self.golden_V

    # Main function to be called
    def set_output_potential(self,DAQ_voltage):
        '''Function that converts provided power in % into a DAC 0-5V output on the WLS channel.'''
        self.DAQ_voltage_output = DAQ_voltage
        with daq.Task() as task:
            task.ao_channels.add_ao_voltage_chan(self.WLS_port, min_val = 0, max_val = 5) #defautls are -10 to 10
            task.write(DAQ_voltage, auto_start=True)
        None

if __name__ == '__main__':
    print('Testing the LEDD1B driver...\n')
    with WhiteLightSourceLEDD1B() as WL:
        print('Current output potential is {0} V.\n'.format(WL.DAQ_voltage_output))
        print([ao.name for ao in WL.DAQ_device.ao_physical_chans])
        print('\nWLS-DAQ/ao channels verified. Driver seems to work.\n')