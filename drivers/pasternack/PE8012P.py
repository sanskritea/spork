# -*- coding: utf-8 -*-
'''
    drivers.pasternack.PE8012P
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    TODO: Figure out the AO port for the Schottky, manually measure/calibrate the RF power vs V curve, and maybe handle errors
    if/when the voltage is out of the spline's range 

    :copyright: 
    C.Egerstrom - Mar 2023.
    Use as you wish, be evil.

    * Based off of mpdAPD driver by N.Delegan from Mar 2020

    :description:
    Driver for PE8012P Schottky photodiode 
    
    :versions/changelog: 
    * V0 - Mar 2023 - Driver started

    :dependencies: 
    *Numpy, Scipy.interpolate, nidaqmx
    *Properly connected NI DAQ
	
	:usage:
    TO-DO
  
'''

import numpy as np
import nidaqmx as daq
from scipy.interpolate import CubicSpline


class PE8012P_Schottky():

    def __init__(self):
        self.schottky_port = 'Dev1/ai0' #TO-DO, FIGURE THIS OUT
        self.DAQ_name = 'Dev1'

        VvsParr = np.loadtxt(r"C:\Users\adminnd\code\jasper-repo\drivers\pasternack\SpecTransferCurve.csv", delimiter=",", dtype=float, skiprows=1)
        self.PvsVfxn = CubicSpline(VvsParr[:,1], VvsParr[:,0], extrapolate=False)


    # Initialize function to get the DAQ_DEVICE OBJECT
    def __enter__(self):
        '''Connect to the physical device'''
        # Identify the local machine
        this_system = daq.system.System.local()
        # We know on our current system that our DAQ is named 'DAQ1'
        self.DAQ_device = this_system.devices[self.DAQ_name]
        print('Schottky Connected')
        return(self)


    def __exit__(self, *args):
        #Nothing needed here for now
        pass


    # Main functions to be called:
    def get_input_potential(self):
        '''Returns the voltage (in V) into the ADC from the Schottky diode'''
        with daq.Task() as task:
            task.ai_channels.add_ai_voltage_chan(self.schottky_port, min_val= 0.0, max_val=2.5)#, terminal_config=daq.constants.TerminalConfiguration.DEFAULT) #Seems like a reasonable range based on spec
            return(task.read())

    
    def get_RF_power(self):
        '''Returns the RF power (in dbm) into the circulator based on the manufacturer specified transfer fxn
        (for a different Schottky)''' #TO-DO: Do manual calibration with a spectrum analyzer to get better data
        print(self.get_input_potential(), self.PvsVfxn(1000*self.get_input_potential()))
        return(20+self.PvsVfxn(1000*self.get_input_potential()))



if __name__ == '__main__':
    print('Testing the PE8012P Schottky driver...\n')
    with PE8012P_Schottky() as schottky:
        print('Current output potential is {0} V.\n'.format(schottky.get_input_potential()))
        print([ao.name for ao in schottky.DAQ_device.ai_physical_chans])
        print('\nSchottky-DAQ/ao channels verified. Driver seems to work.\n')