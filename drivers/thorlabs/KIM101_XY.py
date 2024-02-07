# -*- coding: utf-8 -*-
"""
    drivers.thorlabs.KIM101_XY
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    :copyright: 
    D.Mark, C.Egerstrom - Mar 2023.
    Use as you wish, be evil.

    N.Delegan - Sep 2020.
    Use as you wish, be nice.

    :description:
    Driver for the KIM101 controlled cube with X and Y piezos
    
    :versions/changelog: 
    * V0 - Nov 2020 -  Driver started
    * V1 - Mar 2023 -  Driver de-Lantz-ed, relative moves implemented

    :dependencies: 
    Requires Kinesis drivers installed and the path needs to be added
    to the appropriate variable.
	
	:usage:
"""


#TO-DO: DEAL WITH ENDING/RESTARTING CONNECTION BETTER

#DEFAULT DIST UNIT IS MM, AND IS ONLY APPROXIMATE

#import the things
import sys, time, clr

class KIM_XY_CTLR():

    '''    STATIC VARIABLES    '''
    DEFAULT_UNITS_VOLTAGE = 'V'
    DEFAULT_UNITS_DISTANCE = 'mm'

    MM_PER_STEP = 20E-6 #Approximate, 20nm/step
    MAX_STEP_RATE = 2000 #per s, roughly 40us/s
    MIN_STEP_RATE = 1 #per s, roughly 20nm/s
    DEFAULT_TIMEOUT = 60000 #Assuming ms. TO-DO: CONFIRM THIS


    def __init__(self):
        kinesis_path = r'C:\Program Files\Thorlabs\Kinesis'
        # Add Thorlabs Kinesis Program directory to PATH
        sys.path.append(kinesis_path)

        # Load the relevant assemblies (cube specific)
        clr.AddReference('Thorlabs.MotionControl.DeviceManagerCLI')
        clr.AddReference('Thorlabs.MotionControl.KCube.InertialMotorCLI')        


    def __enter__(self, serialNo = '97000137'):
        '''Connect to KIM
        Arguments:  *serialNo. Default: 97000137'''
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
        #print(serialNumbers) #DEBUG
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
            print(f'Opening KIM device {serialNo}.')
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
        return(self)


    def __exit__(self, *args):
        self.KCubeXY.StopPolling()
        self.KCubeXY.Disconnect(True)
        #self.KCubeXY.ShutDown() # Can use .Disconnect(True) as well  


    #steprate fns
    def getChStepRate(self, channel):
        '''Returns step rate for specified channel in (20ish nm) piezo steps/s
        Arguments:  *channel, e.g. 'x' or 'y' '''
        return self.currentDeviceSettings.Drive.Channel(channel).StepRate


    def setChStepRate(self, channel, stepRate): #in steps/s, in [1, 2k] (???)
        '''Sets step rate for specified channel in specified (20ish nm) piezo steps/s
        Arguments:  *channel, e.g. 'x' or 'y'
                    *stepRate, in steps/s, must be between 1 and 2000 '''
        try:
            self.currentDeviceSettings.Drive.Channel(channel).StepRate = stepRate
            self.KCubeXY.SetSettings(self.currentDeviceSettings, True, True)
        except:
            raise RuntimeError(f'Unable to set a new step rate for {stepRate}.')


    def getChStepAccel(self, channel):
        '''Returns step acceleration rate for specified channel in (20ish nm) piezo steps/s^2
        Arguments:  *channel, e.g. 'x' or 'y'
        Note: Deceleration cannot be specified '''
        return self.currentDeviceSettings.Drive.Channel(channel).StepAcceleration


    def setChStepAccel(self, channel, stepAccel): #in steps/s/s, in [1, 100k] (???)
        '''Sets step rate for specified channel in specified (20ish nm) piezo steps/s^2
        Arguments:  *channel, e.g. 'x' or 'y'
                    *stepAccel, in steps/s, must be between 1 and 100000 
        Note: Deceleration cannot be specified'''
        try:
            self.currentDeviceSettings.Drive.Channel(channel).StepAcceleration = stepAccel
            self.KCubeXY.SetSettings(self.currentDeviceSettings, True, True)
        except:
            raise RuntimeError(f'Unable to set a new step acceleration rate for {stepAccel}.')


    #main movement functions
    def internalRawMoveTo(self, channel, position, timeout=DEFAULT_TIMEOUT):
        '''Old main movement function. Moves a specified channel to a given position (in piezo steps)
        Arguments:  *channel, e.g. 'x' or 'y' 
                    *position, to move to
                    *timeout, in ms. Default: 60000'''
        try: 
            # Choo choo
            self.KCubeXY.MoveTo(channel, round(position), timeout) 
        except:
            # If does not get there.
            raise RuntimeError(f'Failed to move to position {position}.')


    def internalAdaptiveMoveTo(self, channel, position, timeout=DEFAULT_TIMEOUT):
        '''Move a channel to position in piezo steps (approximately 20nm) at an adaptive rate 
        based on the distance to travel above 500nm, and those last 500nm are travelled at 1um/s
        since deceleration can't be controlled
        Arguments:  *channel, e.g. 'x' or 'y' 
                    *position, in piezo steps
                    *timeout, in ms. Default: 60000'''
        try: 
            # Choo choo
            startPos = self.KCubeXY.GetPosition(channel)
            deltaPos = position-startPos
            if deltaPos/self.MAX_STEP_RATE > timeout/1000:
                print("You're probably going to timeout") #TO-DO: RAISE AN EXCEPTION HERE???
            if abs(deltaPos) > 25:
                #will either try to smoothly travel from start to end in 1s at 
                #adaptive speed (or 100nm/s if that's faster), or go at MAX_STEP_RATE for a longer time
                rate = max(min(self.MAX_STEP_RATE, 2*round(abs(deltaPos))), 50)
                self.setChStepRate(channel, rate)
                self.setChStepAccel(channel, 2*rate)
                print(position, round(position+25*(-1)**(deltaPos>0))) #DEBUG
                self.KCubeXY.MoveTo(channel, round(position+25*(-1)**(deltaPos>0)), timeout-100)

            self.setChStepRate(channel, 50) #last bit at 1um/s = 50step/s for up to 0.5s = 500nm = 25steps
            self.KCubeXY.MoveTo(channel, round(position), 10000) #should only take 100ms max, not sure why it needs so long but 1k didn't work. Maybe piezo PID takes a minute to update???
        except:
            # If does not get there.
            raise RuntimeError(f'Failed to move to position {position}.')


    #welcome, Alice!
    def getX(self):
        '''Returns current X position in mm (approximate)'''
        return self.MM_PER_STEP*self.KCubeXY.GetPosition(getattr(self.InertialMotorStatus.MotorChannels, 'Channel1'))


    def getY(self):
        '''Returns current Y position in mm (approximate)'''
        return self.MM_PER_STEP*self.KCubeXY.GetPosition(getattr(self.InertialMotorStatus.MotorChannels, 'Channel2'))
    

    def setX(self, moveToPos, adaptive=True):
        '''Moves to a specified X position in mm (approximate) at an adaptive rate unless otherwise specified
        Arguments:  *moveToPos, absolute position to move to in mm
                    *adaptive, bool. Default is True, otherwise it moves at set rate'''
        if adaptive:
            self.internalAdaptiveMoveTo(getattr(self.InertialMotorStatus.MotorChannels, 'Channel1'), (moveToPos/self.MM_PER_STEP) )
        else:
            self.internalRawMoveTo(getattr(self.InertialMotorStatus.MotorChannels, 'Channel1'), (moveToPos/self.MM_PER_STEP) )


    def setY(self, moveToPos, adaptive=True):
        '''Moves to a specified Y position in mm (approximate) at an adaptive rate unless otherwise specified
        Arguments:  *moveToPos, absolute position to move to in mm
                    *adaptive, bool. Default is True, otherwise it moves at set rate'''
        if adaptive:
            self.internalAdaptiveMoveTo(getattr(self.InertialMotorStatus.MotorChannels, 'Channel2'), (moveToPos/self.MM_PER_STEP) )
        else:
            self.internalRawMoveTo(getattr(self.InertialMotorStatus.MotorChannels, 'Channel2'), (moveToPos/self.MM_PER_STEP))


    def relMoveX(self, relPos, adaptive=True):
        '''Moves to a specified relative X position in mm (approximate) at an adaptive rate unless otherwise specified
        Arguments:  *relPos, relative position to move to (in mm)
                    *adaptive, bool. Default is True, otherwise it moves at set rate'''
        moveToPos = relPos + self.getX()
        self.setX(moveToPos, adaptive)


    def relMoveY(self, relPos, adaptive=True):
        '''Moves to a specified relative Y position in mm (approximate) at an adaptive rate unless otherwise specified
        Arguments:  *relPos, relative position to move to (in mm)
                    *adaptive, bool. Default is True, otherwise it moves at set rate'''
        moveToPos = relPos + self.getY()
        self.setY(moveToPos, adaptive)


    #set X or Y position as 0, respectively
    def setXZero(self):
        '''Sets current X position as 0 position (and returns that 0 position)'''
        self.KCubeXY.SetPositionAs(self.InertialMotorStatus.MotorChannels.Channel1, 0)
        return(self.getX())


    def setYZero(self):
        '''Sets current Y position as 0 position (and returns that 0 position)'''
        self.KCubeXY.SetPositionAs(self.InertialMotorStatus.MotorChannels.Channel2, 0)
        return(self.getY())



if __name__ == "__main__":
    print("Testing the KIM101_XY piezo controller...\n")
    serialNo = '97000137'
    with KIM_XY_CTLR(serialNo) as KCube:
        print(r'If you saw no errors before this message, things seem to be working.')
        # print(KCube.perchannel_position['Channel1'])
        # print(KCube.perchannel_position['Channel2'])
        # KCube.perchannel_position['Channel1'] = 11.0500
        # print(KCube.perchannel_position['Channel2'])