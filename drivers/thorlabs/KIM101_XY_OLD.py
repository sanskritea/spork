# -*- coding: utf-8 -*-
"""
    lantz.drivers.thorlabs.KIM101_XY
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    :copyright: 
    N.Delegan - Sep 2020.
    Use as you wish, be nice.

    :description:
    Driver for the KIM101 controlled cube with X and Y piezos
    
    :versions/changelog: 
    * V0 - Nov 2020 -  Driver started

    :dependencies: 
    Requires Kinesis drivers installed and the path needs to be added
    to the appropriate variable.
	
	:usage:
"""

#Import stuff
import sys, time, clr
from lantz import Action, Feat, DictFeat
from lantz import Driver
from lantz import Q_

# import inspect
#Some system wide necessities for interacting with methods trough python.net
# clr.AddReference('System')
# from System import UInt32, Int32

#Define the class, then the constructor (__init__), it follows similar structure to the flipper
class KIM101_XY_Controller(Driver):
    '''Driver class for the KIM101 that has two piezos (PIA50) for XY motion control.
    '''

    def __init__(self, serialNo = '97000137'):
        # serialNo: the device serialNo, default is the one we use, can be redefined.

        # Define path of Kinesis drivers on the system
        kinesis_path = r'C:\Program Files\Thorlabs\Kinesis'
        # Add Thorlabs Kinesis Program directory to PATH
        sys.path.append(kinesis_path)

        # Load the relevant assemblies (cube specific)
        clr.AddReference('Thorlabs.MotionControl.DeviceManagerCLI')
        clr.AddReference('Thorlabs.MotionControl.KCube.InertialMotorCLI')

        # Import namespaces from the Kinesis drivers to this class
        from Thorlabs.MotionControl.DeviceManagerCLI import DeviceManagerCLI as DeviceManagerCLI
        from Thorlabs.MotionControl.KCube.InertialMotorCLI import KCubeInertialMotor as KCubeInertialMotor
        from Thorlabs.MotionControl.KCube.InertialMotorCLI import ThorlabsInertialMotorSettings as ThorlabsInertialMotorSettings
        from Thorlabs.MotionControl.KCube.InertialMotorCLI import InertialMotorStatus as InertialMotorStatus

        # Do the Kinesis interfacing dance with the hardware
        try:
            # Tell the device manager to get the list of all devices connected to the computer
            DeviceManagerCLI.BuildDeviceList()
        except Exception as ex:
            # An error occurred - see ex for details
            raise RuntimeError(f'Exception raised by BuildDeviceList: {ex}.')
        
        # Get available KIM101s connected and check our serial number is correct
        serialNumbers = DeviceManagerCLI.GetDeviceList(KCubeInertialMotor.DevicePrefix_KIM101)
        if serialNo not in serialNumbers:
            # The requested serial number is not a MFF101 or is not connected
            raise ValueError(f'{serialNo} is not a valid KCubeIntertialMotor serial number.')

        # Create the KCube device
        self.KCubeXY = KCubeInertialMotor.CreateKCubeInertialMotor(serialNo)
        if self.KCubeXY == None:
            # An error occured
            raise ValueError(f'{serialNo} is not a KCube.')

        # Open a connection to the device
        try:
            self.KCubeXY.Connect(serialNo)
            print(f'Opening device {serialNo}.')
        except:
            # If connection failed
            raise RuntimeError(f'Failed to open device {serialNo}.')

        # Wait for settings to initialize
        if not self.KCubeXY.IsSettingsInitialized():
            try:
                self.KCubeXY.WaitForSettingsInitialized(5000)
            except:
                raise RuntimeError('Settings failed to initialize.')
        
        # Start device polling (with built in delay)
        self.KCubeXY.StartPolling(250)
        time.sleep(0.5)
        # Enable device
        self.KCubeXY.EnableDevice()
        time.sleep(0.5)

        # Display info about device
        self.deviceInfo = self.KCubeXY.GetDeviceInfo()
        print(f'Initialized {self.deviceInfo.Description} with S/N: {self.deviceInfo.SerialNumber}, and name: {self.deviceInfo.Name}\n')

        # Call the motor configuration on the device to initialize the settings
        InertialMotorConfiguration = self.KCubeXY.GetInertialMotorConfiguration(serialNo) #Get class structure
        self.currentDeviceSettings = ThorlabsInertialMotorSettings.GetSettings(InertialMotorConfiguration)

        # Finally, get the inertialmotorstatus class into the device class
        self.InertialMotorStatus = InertialMotorStatus

        #Debugging
        # print(inspect.getmembers(self.currentDeviceSettings.Drive.Channel(getattr(self.InertialMotorStatus.MotorChannels, 'Channel1'))))

        # End of the handshaking

    def finalize(self):
        #Stop Polling to keep track of the flipper device
        self.KCubeXY.StopPolling()
        self.KCubeXY.ShutDown() # Can use .Disconnect(True) as well  

    # Feats
    meters_per_step = 20E-9 #Approximate

# Thorlabs.MotionControl.KCube.InertialMotorCLI
#       DriveParams.Position_Max = 2000000000;
#       DriveParams.Position_Min = -2000000000;
#       DriveParams.StepSize_Def = 1000;
#       DriveParams.StepSize_Max = 100000000;
#       DriveParams.StepSize_Min = -100000000;
#       DriveParams.StepRate_Def = 1000;
#       DriveParams.JoystickStepRate_Max = 1000;
#       DriveParams.StepRate_Max = 2000;
#       DriveParams.StepRate_Min = 1;
#       DriveParams.Acceleration_Def = 10000;
#       DriveParams.Acceleration_Max = 100000;
#       DriveParams.Acceleration_Min = 1;

    # DictFeat that deals with current and move_to position
    # @DictFeat(keys = ['Channel1', 'Channel2'], units='m') #50mm (~limit of piezo screws in each direction, see above for hardware limits) limits=(-50000000*meters_per_step, 50000000*meters_per_step, meters_per_step), 
    # def perchannel_position(self, key):
    #     return self.meters_per_step*self.KCubeXY.GetPosition(getattr(self.InertialMotorStatus.MotorChannels, key))
    # @perchannel_position.setter
    # def perchannel_position(self, key, goto_position):
    #     print('moving XY stage')
    #     self.move_to(getattr(self.InertialMotorStatus.MotorChannels, key), (goto_position.magnitude/self.meters_per_step.magnitude))
    
    # DictFeat to set the stepsize value to take in driver units
    # @DictFeat(keys = ['Channel1', 'Channel2'], limits=(-100000000, 100000000)) #Default is 1000
    # def channel_stepsize(self, key):
    #     return self.currentDeviceSettings.Drive.Channel(getattr(self.InertialMotorStatus.MotorChannels, key)).StepSize
    # @channel_stepsize.setter
    # def channel_stepsize(self, key, step):
    #     try:
    #         self.currentDeviceSettings.Drive.Channel(getattr(self.InertialMotorStatus.MotorChannels, key)).StepSize = step
    #         self.KCubeXY.SetSettings(self.currentDeviceSettings, True, True)
    #     except:
    #         raise RuntimeError(f'Unable to set a new step size for {key}.')
    # Note: Step size is not an existing method. This might not apply to our piezo drivers, or we are doing something wrong.

    # TODO: replace with dictfeat once that has been fixed
    @Feat(units='m', limits=(-50000000*meters_per_step, 50000000*meters_per_step, meters_per_step)) #50mm (~limit of piezo screws in each direction, see above for hardware limits) 
    def position_X(self):
        return self.meters_per_step*self.KCubeXY.GetPosition(getattr(self.InertialMotorStatus.MotorChannels, 'Channel1'))
    @position_X.setter
    def position_X(self, goto_position):
        self.move_to(getattr(self.InertialMotorStatus.MotorChannels, 'Channel1'), (goto_position/self.meters_per_step))

    @Feat(units='m', limits=(-50000000*meters_per_step, 50000000*meters_per_step, meters_per_step)) #50mm (~limit of piezo screws in each direction, see above for hardware limits) 
    def position_Y(self):
        return self.meters_per_step*self.KCubeXY.GetPosition(getattr(self.InertialMotorStatus.MotorChannels, 'Channel2'))
    @position_Y.setter
    def position_Y(self, goto_position):
        self.move_to(getattr(self.InertialMotorStatus.MotorChannels, 'Channel2'), (goto_position/self.meters_per_step))

    # DictFeat to set the steprate value to take in driver units
    @DictFeat(keys = ['Channel1', 'Channel2'], limits=(1, 2000, 1)) #Default is 1000
    def channel_steprate(self, key):
        return self.currentDeviceSettings.Drive.Channel(getattr(self.InertialMotorStatus.MotorChannels, key)).StepRate #if errors try get_ and set_
    @channel_steprate.setter
    def channel_steprate(self, key, step):
        try:
            self.currentDeviceSettings.Drive.Channel(getattr(self.InertialMotorStatus.MotorChannels, key)).StepRate = step
            self.KCubeXY.SetSettings(self.currentDeviceSettings, True, True)
        except:
            raise RuntimeError(f'Unable to set a new step rate for {key}.')
    
    # DictFeat to set the stepaccel value to take in driver units
    @DictFeat(keys = ['Channel1', 'Channel2'], limits=(1, 100000, 1)) #Default is 10000
    def channel_stepaccel(self, key):
        return self.currentDeviceSettings.Drive.Channel(getattr(self.InertialMotorStatus.MotorChannels, key)).StepAcceleration
    @channel_stepaccel.setter
    def channel_stepaccel(self, key, stepaccel):
        try:
            self.currentDeviceSettings.Drive.Channel(getattr(self.InertialMotorStatus.MotorChannels, key)).StepAcceleration = stepaccel
            self.KCubeXY.SetSettings(self.currentDeviceSettings, True, True)
        except:
            raise RuntimeError(f'Unable to set a new step acceleration rate for {key}.')

    @Action()
    def zero_X(self):
        self.KCubeXY.SetPositionAs(self.InertialMotorStatus.MotorChannels.Channel1, 0)
        self.position_X
    
    @Action()
    def zero_Y(self):
        self.KCubeXY.SetPositionAs(self.InertialMotorStatus.MotorChannels.Channel2, 0)
        self.position_Y
       
    # Main movement function.
    def move_to(self, channel, position, timeout=60000): #Timeout given 60s to get there
        try: # Choo choo
            self.KCubeXY.MoveTo(channel, round(position), timeout) 
        except:
            # If does not get there.
            raise RuntimeError(f'Failed to move to position {position}.')

# For driver testing
if __name__ == "__main__":
    print("Testing the KIM101_XY piezo controller...\n")
    serialNo = '97000137'
    with KIM101_XY_Controller(serialNo) as KCube:
        print(r'If you saw no errors before this message, things seem to be working.')
        # print(KCube.perchannel_position['Channel1'])
        # print(KCube.perchannel_position['Channel2'])
        # KCube.perchannel_position['Channel1'] = 11.0500
        # print(KCube.perchannel_position['Channel2'])