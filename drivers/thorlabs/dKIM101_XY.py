# -*- coding: utf-8 -*-

#DEFAULT UNIT IS MM, AND IS ONLY APPROXIMATE

#import the things
import sys, time, clr

class KIM_XY_CTLR():

    def __init__(self, serialNo = '97000137'):
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

    '''
    VARIABLES
    '''

    mm_per_step = 20E-6 #Approximate


    def finalize(self):
        #Stop Polling to keep track of the flipper device
        self.KCubeXY.StopPolling()
        self.KCubeXY.ShutDown() # Can use .Disconnect(True) as well  


    #welcome, Alice!
    def getX(self):
        return self.mm_per_step*self.KCubeXY.GetPosition(getattr(self.InertialMotorStatus.MotorChannels, 'Channel1'))
    
    def setX(self, moveTo):
        self.moveTo(getattr(self.InertialMotorStatus.MotorChannels, 'Channel1'), (moveTo/self.mm_per_step) )

    def getY(self):
        return self.mm_per_step*self.KCubeXY.GetPosition(getattr(self.InertialMotorStatus.MotorChannels, 'Channel2'))

    def setY(self, moveTo):
        self.moveTo(getattr(self.InertialMotorStatus.MotorChannels, 'Channel2'), (moveTo/self.mm_per_step))

    #steprate fns

    def chStepRate(self, key):
        return self.currentDeviceSettings.Drive.Channel(getattr(self.InertialMotorStatus.MotorChannels, key)).StepRate

    def setChStepRate(self, key, step):
        try:
            self.currentDeviceSettings.Drive.Channel(getattr(self.InertialMotorStatus.MotorChannels, key)).StepRate = step
            self.KCubeXY.SetSettings(self.currentDeviceSettings, True, True)
        except:
            raise RuntimeError(f'Unable to set a new step rate for {key}.')

    def chStepAccel(self, key):
        return self.currentDeviceSettings.Drive.Channel(getattr(self.InertialMotorStatus.MotorChannels, key)).StepAcceleration

    def setChStepAccel(self, key, stepaccel):
        try:
            self.currentDeviceSettings.Drive.Channel(getattr(self.InertialMotorStatus.MotorChannels, key)).StepAcceleration = stepaccel
            self.KCubeXY.SetSettings(self.currentDeviceSettings, True, True)
        except:
            raise RuntimeError(f'Unable to set a new step acceleration rate for {key}.')

    #return to 0 X or Y positions, respectively
    def zeroX(self):
        self.KCubeXY.SetPositionAs(self.InertialMotorStatus.MotorChannels.Channel1, 0)
        self.posX

    def zeroY(self):
        self.KCubeXY.SetPositionAs(self.InertialMotorStatus.MotorChannels.Channel2, 0)
        self.posY

    #main movement function
    def moveTo(self, channel, position, timeout=60000):
        try: # Choo choo
                self.KCubeXY.MoveTo(channel, round(position), timeout) 
        except:
            # If does not get there.
            raise RuntimeError(f'Failed to move to position {position}.')


if __name__ == "__main__":
    print("Testing the KIM101_XY piezo controller...\n")
    serialNo = '97000137'
    with KIM101_XY_Controller(serialNo) as KCube:
        print(r'If you saw no errors before this message, things seem to be working.')
        # print(KCube.perchannel_position['Channel1'])
        # print(KCube.perchannel_position['Channel2'])
        # KCube.perchannel_position['Channel1'] = 11.0500
        # print(KCube.perchannel_position['Channel2'])