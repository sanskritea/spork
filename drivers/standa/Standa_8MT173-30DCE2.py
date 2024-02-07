# -*- coding: utf-8 -*-
'''
    drivers.standa.8MT173-30DCE2.StandaZFocusDriver 
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    :copyright: 
    D.Mark, C.Egerstrom - Mar 2023
    De-lantz-ed version of Nazar's old driver. Use as you wish, be evil.

    N.Delegan - Sep 2020.
    Use as you wish, be nice.

    :description: 
    It is a modified/wrapped/re-skinned version of the developper package provided by Standa. 
    I moved the .dll and the python wrapper functions into the subfolder 'StandaLibraries' for clarity sake. 
    Also, for the device to be found the base driver needs to be installed in windows 
    (from the dev package, there is a driver folder).
    Final note, being a USB driver, one instance of the driver can be used at a time, no double dipping. 
    The driver is written to accomodate one device of this type, if multiple are used, 
    see 'self.dev_count' definition below to modify accordingly. See development kit to modify driver for multiple instances of this hardware.
    Many of the structure definitions and value assignments are done in pyximc.py, code is best understood with that in-hand.
    
    MOVE RIGHT = AWAY FROM THE SAMPLE IN THE JASPER SETUP
    MOVE LEFT = TOWARDS THE SAMPLE IN THE JASPER SETUP

    LOW COUNTS = MOVE LEFT, TOWARDS
    HIGH COUNTS = MOVE RIGHT, AWAY

    :versions/changelog:
    V1.0 - Oct 2020 - Driver written.
    V1.1 - Dec 2020 - Tested and cleaned.

    :dependencies:
    Base driver needs to be installed for hardware detection, otherwise specify the adress. 
    bindy.dll libjximc.dll libximc.dll libximx.pyximc.lib pyximc.py ximwrapper.dll in StandaLibraries subfolder

'''

# Import various dependencies required for the STANDA driver:
from ctypes import byref, c_uint, POINTER, c_int, create_string_buffer, cast, string_at, c_uint32

class StandaZFocusDriver():
    
    DEFAULT_UNITS_DISTANCE = 'mm'
    # Real unit definitions:
    #CALIB_LENGTH = 0.343000 # m
    CALIB_LENGTH = 343 #mm
    CALIB_STEPS = 5926176 # encoder units
    DIST_PER_STEP = CALIB_LENGTH/CALIB_STEPS # length per calibration unit


    # Border state (limit flags), if encoder is used, have it set to BORDER IS ENCODER, otherwise use 0x08 for hard stops
    BORDER_FLAG_VALUES = {
        'CANNOT GET BORDER STATUS (Likely device not properly initiated)':0x00,
        'BORDER_IS_ENCODER':0x01, 
        'BORDER_STOP_LEFT':0x02,
        'BORDER_STOP_RIGHT':0x04,
        'BORDERS_SWAP_MISSET_DETECTION':0x08,
        'BORDER_STOP_LEFT_AND_RIGHT':0x06,
        'BORDER_STOP_LEFT_AND_RIGHT_AND_ENCODER':0x07
    }


    def __init__(self):
        '''Constructor function, looks and opens the device, prints error if unable.'''
        try:
            from StandaLibraries.pyximc import lib, EnumerateFlags, controller_name_t, device_information_t, get_position_t, edges_settings_t, move_settings_t, engine_settings_t, Result
        except:
            raise RuntimeError(f'Unable to load libraries into namespace, check that pyximc.py and its dependencies are localted in the StandaLibraries subfolder.')

        # Place things into appropriate places to follow typical lib. notation
        self.lib = lib
        self.EnumerateFlags = EnumerateFlags
        self.controller_name_t = controller_name_t
        self.device_information_t = device_information_t
        self.get_position_t = get_position_t
        self.edges_settings_t = edges_settings_t
        self.move_settings_t = move_settings_t
        self.engine_settings_t = engine_settings_t
        self.Result = Result


    def __enter__(self):
        '''Connects to device'''
        # This is device search and enumeration with probing. It gives more information about devices.
        enum_hints = b'addr=COM5'
        probe_flags = self.EnumerateFlags.ENUMERATE_PROBE #+ self.EnumerateFlags.ENUMERATE_ALL_COM # Tell the thing where to look
     
        # enum_hints = b'addr=' # Use this hint string for broadcast enumerate
        try:
            self.devenum = self.lib.enumerate_devices(probe_flags, enum_hints) 
            self.dev_count = self.lib.get_device_count(self.devenum) # If more devices are connected, adjust code, see dev package
        except:
            raise RuntimeError(f'Unable to enumerate devices, check if device is connected.')

        # Get device name, takes name of the first one, see dev kit if you are using multiples.
        open_name = None
        open_name = self.lib.get_device_name(self.devenum, 0)
        
        # Get device_id which will serve as an identifier for all future library calls
        try:
            self.device_id = self.lib.open_device(open_name)
        except:
            raise RuntimeError(f'Unable to open device and assign device ID.')
        
        # Initial settings, odds and ends.
        self.step_size = 5 #in motor steps (rev's?), initial default movement step size #TODO: Confirm these units are right
        self.pos_soft_home = 200000 #in motor steps (rev's?), soft home for sample changes and such #TODO: Is this still right?
        self.pos_engaged = int(10.5/(self.CALIB_LENGTH/self.CALIB_STEPS)) #position for 'engaged' length #TODO: Is this still right?
        print('Standa connected')
        return(self)

    # Finalization function, to be ran when closed:
    def __exit__(self, *args):
        '''Finalize function, closes the device that was opened.'''
        # Close the device one this class closes
        self.lib.command_sstp(self.device_id)  
        self.lib.close_device(byref(cast(self.device_id, POINTER(c_int))))


    def idn(self):
        '''Identification function.
        Returns: model name, device info, and serial number in a string'''
        # Get model name
        controller_name = self.controller_name_t() #See pyximc.py for structure definition
        self.lib.get_enumerate_device_controller_name(self.devenum, 0, byref(controller_name))

        # Get device info
        x_device_information = self.device_information_t() #See pyximc.py
        self.lib.get_device_information(self.device_id, byref(x_device_information))
        # if result == Result.Ok: #For Debugging

        # Get serial
        x_serial = c_uint()
        self.lib.get_serial_number(self.device_id, byref(x_serial))

        # Generate output
        return 'Device version: {0}.{1}.{2} || Model: {3} {4} || Serial #: {5}'.format(repr(x_device_information.Major),repr(x_device_information.Minor),repr(x_device_information.Release), repr(controller_name.ControllerName)[2:-2],repr(string_at(x_device_information.ProductDescription).decode())[1:-1],repr(x_serial.value))


    def minIncrement(self):
        '''Identify minimum step size (in nm) based on spec.'''
        return 'Approx. min. step: {:.1f} nm'.format(self.DIST_PER_STEP*(10**6))


    def getPos(self):
        '''Returns position (in mm)'''
        # Get your f position
        f_pos = self.get_position_t() #See pyximc.py
        self.lib.get_position(self.device_id, byref(f_pos))
        return f_pos.Position*self.DIST_PER_STEP


    def setPos(self, value):
        '''Move to a given z-position (in mm)'''
        self.movement_function(value/self.DIST_PER_STEP)


    def relMove(self, value):
        '''Move a relative amount in z (in mm)'''
        self.movement_function( (self.getPos()+value)/self.DIST_PER_STEP)


    def posLimitFlags(self):
        '''Border encoder flags. i.e. limits are either set by encoder (soft limits) or by hardware switches.'''
        self.f_pos_limits = self.edges_settings_t() #See pyximc.py
        self.lib.get_edges_settings(self.device_id, byref(self.f_pos_limits))
        return self.f_pos_limits.BorderFlags


    #units are in mm
    def getPosLimitLeft(self):
        '''Returns soft limit (in mm) for the LEFT border (SAMPLE SIDE)'''
        self.f_pos_limits = self.edges_settings_t() #See pyximc.py
        self.lib.get_edges_settings(self.device_id, byref(self.f_pos_limits))
        return self.f_pos_limits.LeftBorder*self.DIST_PER_STEP


    def setPosLimitLeft(self, lLimit):
        '''Sets soft limit for the LEFT border (SAMPLE SIDE) at the specified lLimit (in mm)
        Arguments:  *lLimit'''
        if lLimit/self.DIST_PER_STEP > self.f_pos_limits.RightBorder:
            print('ERROR: Trying to set left border (sample side) to a higher value than the right border.')
        else:
            self.f_pos_limits.LeftBorder = c_int(int(lLimit/self.DIST_PER_STEP))
            self.lib.set_edges_settings(self.device_id, byref(self.f_pos_limits)) 


    def getPosLimitRight(self):
        '''Returns soft limit (in mm) for the RIGHT border (AWAY FROM SAMPLE SIDE)'''
        self.f_pos_limits = self.edges_settings_t() #See pyximc.py
        self.lib.get_edges_settings(self.device_id, byref(self.f_pos_limits))
        return self.f_pos_limits.RightBorder*self.DIST_PER_STEP
    

    def setPosLimitRight(self, rLimit):
        '''Sets soft limit for the RIGHT border (AWAY FROM SAMPLE SIDE) at the specified rLimit (in mm)
        Arguments:  *rLimit'''
        if rLimit/self.DIST_PER_STEP < self.f_pos_limits.LeftBorder:
            print('ERROR: Trying to set right border (away-from-sample side) to a lower value than the left border.')
        else:
            self.f_pos_limits.RightBorder = c_int(int(rLimit/self.DIST_PER_STEP))
            self.lib.set_edges_settings(self.device_id, byref(self.f_pos_limits))


    def moveSpeed(self):
        '''Returns move speed (in rev/s)''' #TODO: Confirm these units are right
        self.move_settings = self.move_settings_t()
        # Get current move settings from controller
        self.lib.get_move_settings(self.device_id, byref(self.move_settings))
        return self.move_settings.Speed


    def setMoveSpeed(self, newSpeed):
        '''Sets move speed (in rev/s)
        Arguments:  *newSpeed''' #TODO: Confirm these units are right
        self.move_settings.Speed = c_uint(newSpeed)
        self.lib.set_move_settings(self.device_id, byref(self.move_settings))


    def getMoveAccel(self):
        '''Returns move acceleration (in revs/s^2)''' #TODO: Confirm these units are right
        self.move_settings = self.move_settings_t()
        # Get current move settings from controller
        self.lib.get_move_settings(self.device_id, byref(self.move_settings))
        return self.move_settings.Accel


    def setMoveAccel(self, newAcc):
        '''Set move acceleration (in revs/s^2)
        Arguments:  *newAcc''' #TODO: Confirm these units are right
        self.move_settings.Accel = c_uint(newAcc)
        self.lib.set_move_settings(self.device_id, byref(self.move_settings))


    def getMoveDecel(self):
        '''Returns move deceleration (in revs/s^2)''' #TODO: Confirm these units are right
        self.move_settings = self.move_settings_t()
            # Get current move settings from controller
        self.lib.get_move_settings(self.device_id, byref(self.move_settings))
        return self.move_settings.Decel


    def setMoveDecel(self, newDecel):
        '''Set move acceleration (in revs/s^2)
        Arguments:  *newDecel''' #TODO: Confirm these units are right
        self.move_settings.Decel = c_uint(newDecel)
        self.lib.set_move_settings(self.device_id, byref(self.move_settings))


    def getMoveAntiplaySpd(self):
        '''Returns the move_settings antiplay speed'''
        self.move_settings = self.move_settings_t()
        # Get current move settings from controller
        self.lib.get_move_settings(self.device_id, byref(self.move_settings))
        return self.move_settings.AntiplaySpeed   


    def getEngineAntiplaySteps(self):
        '''Returns the engine_settings antiplay speed'''
        self.engine_settings = self.engine_settings_t()
        # Get current move settings from controller
        self.lib.get_engine_settings(self.device_id, byref(self.engine_settings))
        return self.engine_settings.Antiplay


    def getMoveStepSize(self):
        '''Returns step size (in mm)'''
        return self.step_size*self.DIST_PER_STEP


    def setMoveStepSize(self, stepSize):
        '''Sets step size to a given value (in mm)
        Arguments:  *stepSize'''
        self.step_size = stepSize/self.DIST_PER_STEP


    def softStop(self):
        '''Sends command to Standa to soft stop'''
        self.lib.command_sstp(self.device_id)


    #movement relative to sample in set moveStepSize
    def moveAwayStep(self):
        '''Move sample forward +1 step in current stepsize'''
        self.getPos = self.getPos(self.step_size*self.DIST_PER_STEP)


    def moveTowardStep(self):
        '''Move sample backward -1 step in current stepsize'''
        self.getPos = self.getPos(-self.step_size*self.DIST_PER_STEP)


    def moveHome(self):
        '''Return to soft home position''' #TODO: Figure out if 'soft home' position is still relevant
        self.getPos = self.pos_soft_home*self.DIST_PER_STEP


    def moveEngagedPos(self):
        '''Move motor to engaged position''' #TODO: Figure out if 'engaged' position is still relevant
        self.getPos = self.pos_engaged*self.DIST_PER_STEP


    def resetDefaults(self):
        '''Reset default values (that Nazar determined)'''
        # Speed to 2000 (default)
        self.move_speed = 2000

        # Right border is 300000, left border is 111800
        self.setPosLimitLeft = 111800
        self.setPosLimitRight = 300000
        self.f_pos_limits.BorderFlags = c_uint(0x01)
        self.lib.set_edges_settings(self.device_id, byref(self.f_pos_limits)) 

        # Soft Home is 170000
        self.pos_soft_home = 200000 # not that this should ever change, but wtv...
        self.pos_engaged = int(0.0105/self.DIST_PER_STEP) #Confirm this value
        print('Values reset to default for the standa speed ({0}), soft home ({1}), and edge settings (L:{2}.uL:{3}, R:{4}.uR:{5}).'.format(self.move_settings_t().Speed, self.pos_soft_home, self.f_pos_limits.LeftBorder, self.f_pos_limits.uLeftBorder, self.f_pos_limits.RightBorder, self.f_pos_limits.uRightBorder))


    #movement function called for all movements
    def movement_function(self, move_to_pos):
        '''Main movement function that also prevent moving to out-of-bounds positions
        Arguments:  *move_to_pos: desired position to move to (in steps)'''
        # Main movement script

        # Inquiry towards limits
        self.f_pos_limits = self.edges_settings_t()
        self.lib.get_edges_settings(self.device_id, byref(self.f_pos_limits))
       
        # Inquiry towards current position
        self.f_pos = self.get_position_t() #See pyximc.py
        self.lib.get_position(self.device_id, byref(self.f_pos))

        # Check if position asked to move to is inside the boundary conditions
        if move_to_pos < self.f_pos_limits.RightBorder and move_to_pos > self.f_pos_limits.LeftBorder:
            #Do the movement here
            #Figure out timing and wait
            print('FocusMotor : Moving, please wait...')
            self.lib.command_move(self.device_id, int(move_to_pos))
            while (self.lib.command_wait_for_stop(self.device_id,c_uint32(100)) != self.Result.Ok):
                #10ms is default, right now wait is set to 100ms wait time between inqueries, default
                pass
            else:
                print('FocusMotor : Done, position is {0} um.'.format(move_to_pos*self.DIST_PER_STEP*10**6))
        else:
            #Print out error is directed location is out of soft bounds.
            print('FocusMotor : MOVEMENT ERROR\nTrying to move to position {0} um from current position {1} um. \nDictated position is out of the movement limits: {2} um and {3} um.'.format(move_to_pos*self.DIST_PER_STEP, self.f_pos.Position*self.DIST_PER_STEP, self.f_pos_limits.LeftBorder*self.DIST_PER_STEP, self.f_pos_limits.RightBorder*self.DIST_PER_STEP))
    

# Testing and validation
if __name__ == '__main__':
    with StandaZFocusDriver() as FocusMotor:
        print('Testing the Standa motor driver...\n')
        print(FocusMotor.idn)
        print('\nIdentification successfull. Driver seems to work.\n')
        print('Current position is {0}'.format(FocusMotor.getPos))
    