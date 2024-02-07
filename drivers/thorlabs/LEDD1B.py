# -*- coding: utf-8 -*-
'''
    lantz.drivers.thorlabs.LEDD1B
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    :copyright:
    C. Egerstrom - Feb 2023
    Use as you wish, be evil 
    N.Delegan - Sep 2020.
    Use as you wish, be nice.

    :description:
    Driver for LEDD1B LED driver
    
    :versions/changelog: 
    * V0 - Oct 2020 -  Driver started
    * V1 - Nov 2020 - Driver tested and working.
    * V2 - Feb 2023 - De-Lantz-ed.

    :dependencies: 
    * Needs a properly wired DAQ to modulate the WLS.
    * nspyre, imports InstrumentGateway so that once all instruments are connected, 
    flipping the APD or white light on will first shut the other one down
    * time, to wait for APD to turn off
	
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

import nidaqmx as daq
from nspyre import InstrumentGateway
import time


# Driver Class
class WhiteLightSourceLEDD1B():
    '''Driver class for the LEDD1B.
    In essence this is a nidaqmx wrapper as the modulation will be controlled through the DAC
    '''
    # Constructor function to define how things are hardwired
    def __init__(self):
        self.WLS_port = 'Dev1/ao2'
        self.DAQ_name = 'Dev1'
        self.GOLDEN_V = 0.225 #goldylocks power setting for faster action
        self.DAQ_voltage_output = 0


    # Initialize function to get the DAq_DEVICE OBJECT
    def __enter__(self):
        # Identify the local machine
        this_system = daq.system.System.local()
        # We know on our current system that our DAQ is named 'DAQ1'
        self.DAQ_device = this_system.devices[self.DAQ_name]
        print('White Light Connected')
        return(self)

    def __exit__(self, *args):
        '''Turns the lamp off'''
        self.set_output_potential(0)


    # Main function to be called
    def set_output_potential(self, DAQ_voltage):
        '''Sets DAC output to WLS to a specified voltage.
        Argumenets:  *DAQ_voltage, in V, between 0 and 5'''
        if DAQ_voltage < 0 or DAQ_voltage > 5:
            raise ValueError(f'Requested DAQ Voltage of {DAQ_voltage}V is out of range [0V,5V]')

        if DAQ_voltage > 0: #if you're trying to turn on the lamp, first kill the APD
            with InstrumentGateway() as gw: #TODO: Implement some sort of error handling if this isn't in the GW
                gw.apdGate.apdOff()
                time.sleep(0.5) #wait 500ms for APD to turn off
        self.DAQ_voltage_output = DAQ_voltage
        with daq.Task() as task:
            task.ao_channels.add_ao_voltage_chan(self.WLS_port, min_val = 0, max_val = 5) #defautls are -10 to 10
            task.write(DAQ_voltage, auto_start=True)


    def golden_light(self):
        '''Set the lamp to a nice, medium light level'''
        self.set_output_potential(self.GOLDEN_V)


    def max_light(self):
        '''Turn the lamp output to max (5V)'''
        self.set_output_potential(5)


    def off_light(self):
        '''Turn the lamp output to 0V'''
        self.set_output_potential(0)


    def getOutputV(self):
        '''Returns the DAQ output V to the WLS (in V)'''
        return(self.DAQ_voltage_output)



if __name__ == '__main__':
    print('Testing the LEDD1B driver...\n')
    with WhiteLightSourceLEDD1B() as WL:
        print('Current output potential is {0} V.\n'.format(WL.DAQ_voltage_output))
        print([ao.name for ao in WL.DAQ_device.ao_physical_chans])
        print('\nWLS-DAQ/ao channels verified. Driver seems to work.\n')