# -*- coding: utf-8 -*-
'''
    drivers.thorlabs.LEDD1B
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    :copyright:
    C.Egerstrom - Jan 2023.
    Use as you wish, be evil.

    Lantz-free version heavily based on driver by
    N.Delegan - Sep 2020.
    Use as you wish, be nice.

    :description:
    Driver for MPD PDM series APD
    
    :versions/changelog: 
    * V0_OLD - Oct 2020 -  Driver started
    * V1_OLD - Nov 2020 - Driver tested and working.
    * V2 - Jan 2023 - De-lantzed the driver

    :dependencies: 
    * nidaqmx library
    * Needs a properly wired DAQ to gate the APD. Note that this should not be used for fast gating.
    This is limited to turning counts on and off.
    * nspyre, imports InstrumentGateway so that once all instruments are connected, 
    flipping the APD or white light on will first shut the other one down
    * time, to wait for lamp to turn off
	
	:usage:
    APD gate is wired into BOX2 AO0 so, should be DAQ2/ao1

    For practical usage, just turn APD ON/OFF
'''

import nidaqmx
from nspyre import InstrumentGateway
import time

class MPD_PDMSeries_APD():

    def __init__(self):
        self.APDGate_port = 'Dev1/ao3' # We know on our current system that our DAQ is named 'DAQ1'
        self.DAQ_name = 'Dev1'
        self.gate_voltage = 0 #By default gate to 0


    def __enter__(self):
        '''Actually connect to device'''
        # Initialize function to get the DAQ_DEVICE OBJECT
        # Identify the local machine
        this_system = nidaqmx.system.System.local()
        self.DAQ_device = this_system.devices[self.DAQ_name]
        print('APD Connected')
        return(self)


    # Nothing is needed here
    def __exit__(self, *args):
        self.gate_sp = 'Counting OFF'


    def getGateVoltage(self):
        '''Accessor for gate voltage'''
        return self.gate_voltage


    def apdOn(self):
        '''Turn the APD on (set output potential to 5V)'''
        self.set_output_potential(5)


    def apdOff(self):
        '''Turn the APD off (set output potential to 0V)'''
        self.set_output_potential(0)


    def set_output_potential(self,voltage):
        '''Function that converts provided potential into DAC 0-5V (OFF/ON).'''
        if voltage > 0: #if you're trying to turn on the APD, first kill the lamp
            with InstrumentGateway() as gw: #TODO: Implement some sort of error handling if this isn't in the GW
                gw.whiteLED.off_light()
                time.sleep(0.5) #wait 500ms for lamp to turn off
        self.gate_voltage = voltage
        with nidaqmx.Task() as task:
            task.ao_channels.add_ao_voltage_chan(self.APDGate_port, min_val = 0, max_val = 5) #defautls are -10 to 10
            task.write(self.gate_voltage, auto_start=True)



if __name__ == '__main__':
    print('Testing the LEDD1B driver...\n')
    with MPD_PDMSeries_APD() as APD:
        print('Current APD gate output potential is {0} V.\n'.format(APD.gate_voltage))
        print([ao.name for ao in APD.DAQ_device.ao_physical_chans])
        print('\nAPD-DAQ/ao channels verified. Driver seems to work.\n')
