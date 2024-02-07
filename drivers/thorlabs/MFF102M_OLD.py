#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
    lantz.drivers.thorlabs.MFF102M
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    :copyright: 
    M.Solomon
    N.Delegan - Sep 2020.
    Use as you wish, be nice.

    :description:
    Driver for the MFF102M thorlabs pellicle flipper
    
    :versions/changelog: 
    * V0 - Nov 2020 -  Driver started
    * V1 - Nov 2020 - Driver works and tested

    :dependencies: 
    Requires Kinesis drivers installed and the path needs to be added
    to the appropriate variable.
	
	:usage:
    Set position and read position
  
"""

#Import stuff
import sys, time, clr
from lantz import Action, Feat
from lantz import Driver

#Some system wide necessities for interacting with methods trough python.net
clr.AddReference('System')
from System import UInt32, Int32

class MFF102MFlipper(Driver):
    '''Driver class for the MFF102MF.
    In essence this is a Kinesis driver wrapper for the flipper.
    '''
    def __init__(self,serialNo = '37869793'):
        # serialNo: The FilterFlipper MFF101(2) serial number

        #Define path of Kinesis drivers on the system
        kinesis_path = r'C:\Program Files\Thorlabs\Kinesis'
        #Add Thorlabs Kinesis Program directory to PATH
        sys.path.append(kinesis_path)

        #Load the relevant assemblies (flipper specific)
        clr.AddReference('Thorlabs.MotionControl.FilterFlipperCLI')
        clr.AddReference('Thorlabs.MotionControl.DeviceManagerCLI')

        #Import namespaces from Kinesis drivers to this class
        from Thorlabs.MotionControl.DeviceManagerCLI import DeviceManagerCLI as DeviceManagerCLI
        from Thorlabs.MotionControl.FilterFlipperCLI import FilterFlipper as FilterFlipperCLI

        #Start interfacing with the hardware
        try:
            # Tell the device manager to get the list of all devices connected to the computer
            DeviceManagerCLI.BuildDeviceList()
        except Exception as ex:
            # An error occurred - see ex for details
            raise RuntimeError(f'Exception raised by BuildDeviceList: {ex}.') 

        # Get available FilterFlipper(s) connected and check our serial number is correct
        serialNumbers = DeviceManagerCLI.GetDeviceList(FilterFlipperCLI.DevicePrefix)
        if serialNo not in serialNumbers:
            # The requested serial number is not a MFF101 or is not connected
            raise ValueError(f'{serialNo} is not a valid FilterFlipper serial number.')

        # Create the FilterFlipper device
        self.flipper = FilterFlipperCLI.CreateFilterFlipper(serialNo)
        if self.flipper == None:
            # An error occured
            raise ValueError(f'{serialNo} is not a FilterFlipper.')

        # Open a connection to the device.
        try:
            self.flipper.Connect(serialNo)
            print(f'Opening device {serialNo}.')
        except:
            # Connection failed
            raise RuntimeError(f'Failed to open device {serialNo}.')

        # Wait for the device settings to initialize - timeout 5000ms
        if not self.flipper.IsSettingsInitialized():
            try:
                self.flipper.WaitForSettingsInitialized(5000)
            except:
                raise RuntimeError('Settings failed to initialize.')
        
        # Start the device polling
        # The polling loop requests regular status requests to the motor to ensure the program keeps track of the device.
        self.flipper.StartPolling(250)
        # Needs a delay so that the current enabled state can be obtained
        time.sleep(0.5)
        # Enable flipper (otherwise any move is ignored)
        self.flipper.EnableDevice()
        # Needs a delay so that the device can be enabled
        time.sleep(0.5)

        # Get device settings
        self.currentDeviceSettings = self.flipper.GetDeviceConfiguration(serialNo, 1) #1 = UseDeviceSettings
        # Display info about device
        self.deviceInfo = self.flipper.GetDeviceInfo()

        print(f'Initialized {self.deviceInfo.Description} with S/N: {self.currentDeviceSettings.SerialNo}, and name: {self.deviceInfo.Name}\n')
        print(f'Current position is {self.flipper.Position}.')

    def finalize(self):
        #Stop Polling to keep track of the flipper device
        self.flipper.StopPolling()
        self.flipper.ShutDown()

    #Main position feature and feedback on current position
    @Feat(values={"Up": 2, "Down": 1})
    def flipper_position(self):
        return self.flipper.Position #Current position
    @flipper_position.setter
    def flipper_position(self,value):
        self.move_flipper(value)

    #Main movement function
    def move_flipper(self, position):
        if not (position == 1 or position == 2):
            raise ValueError('Position can only be 1 or 2.')
        try:
            self.flipper.SetPosition.Overloads[UInt32, Int32](position, 60000) 
            #Overloads nomenclature necessary from python.net (use help() to figure out inputs)
            print(f'Moving Device to {position}.')
            time.sleep(0.5)
        except Exception:
            raise RuntimeError('Failed to move to position.')

# For driver testing
if __name__ == "__main__":
    print("Testing the MFF102M pellicle flipper driver...\n")
    serialNo = '37869793'
    with MFF102MFlipper(serialNo) as FM:
        # Get device settings
        # FM.flipper.SetPosition.Overloads[UInt32, Int32](2, 60000)
        print(r'If you saw no errors before this message, things seem to be working.')