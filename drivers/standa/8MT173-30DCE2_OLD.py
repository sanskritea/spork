# -*- coding: utf-8 -*-
'''
    lantz.drivers.standa.8MT173-30DCE2.StandaZFocusDriver 
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    :copyright: 
    N.Delegan - Sep 2020.
    Use as you wish, be nice.

    Chris - Dec 2022.
    Use as you wish, be evil.

    :description: 
    It is a modified/wrapped/re-skinned version of the developper package provided by Standa. 
    Nazar moved the .dll and the python wrapper functions into the subfolder 'StandaLibraries' for clarity sake. 
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
    V1.2 - Dec 2022 - Removed Lantz dependencies

    :dependencies:
    Base driver needs to be installed for hardware detection, otherwise specify the adress. 
    bindy.dll libjximc.dll libximc.dll libximx.pyximc.lib pyximc.py ximwrapper.dll in StandaLibraries subfolder

'''

# Import various dependencies required for the STANDA driver:
from ctypes import byref, c_uint, POINTER, c_int, create_string_buffer, cast, string_at, c_uint32

# Import the lantz related dependencies:
from lantz import Action, Feat, Q_

# Finally we can start with the driver itself:
class StandaZFocusDriver():
    
    def __init__(self):
        '''Constructor function, looks and opens the device, prints error if unable.'''
        try:
            from StandaLibraries.pyximc import lib, EnumerateFlags, controller_name_t, device_information_t, get_position_t, edges_settings_t, move_settings_t, engine_settings_t, Result
        except:
            raise RuntimeError(f'Unable to load libraries into namespace, check that pyximc.py and its dependencies are localted in the StandaLibraries subfolder.')

        # Place things into appropriate places to follow typical lib. notation
        self.lib = lib
        self.controller_name_t = controller_name_t
        self.device_information_t = device_information_t
        self.get_position_t = get_position_t
        self.edges_settings_t = edges_settings_t
        self.move_settings_t = move_settings_t
        self.engine_settings_t = engine_settings_t
        self.Result = Result

        # This is device search and enumeration with probing. It gives more information about devices.
        enum_hints = b'addr=COM5'
        probe_flags = EnumerateFlags.ENUMERATE_PROBE #+ self.EnumerateFlags.ENUMERATE_ALL_COM # Tell the thing where to look
     
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
        
        '''Initialization function, sets user defined values such as default step size,
         home position (disengaged), and engaged position.'''
        # Initial settings, odds and ends.
        self.step_size = 5 #initial default step size
        self.pos_soft_home = 200000 #soft home for sample changes and such
        self.pos_engaged = int(0.0105/(0.343/5926176)) #position for 'engaged' length
    
    # Finalization function, to be ran when closed:
    def finalize(self):
        '''Finalize function, closes the device that was opened.'''
        # Close the device one this class closes
        self.lib.command_sstp(self.device_id)  
        self.lib.close_device(byref(cast(self.device_id, POINTER(c_int))))
            
    # Real unit definitions:
    calib_length = 0.343000 # m
    calib_steps = 5926176 # encoder units
    dist_per_step = calib_length/calib_steps # length per calibration unit

    # Identify library version
    @Feat(read_once=True)
    def idn_library(self):
        '''Library identification function.'''
        sbuf = create_string_buffer(64) #Create buffer
        self.lib.ximc_version(sbuf) #Place library version in the buffer
        return 'Library version: ' + sbuf.raw.decode().rstrip('\0') #Decode into string what is in the buffer

    # Identify device
    @Feat(read_once=True)
    def idn(self):
        '''Identification function, outputs model name, device info, and serial number'''
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

    @Feat(read_once=True)
    def idn_min_increment(self):
        '''Identify minimum step size based on spec.'''
        return 'Approx. min. step: {:.1f} nm'.format(self.dist_per_step*(10**9))

    # Position of the motor
    @Feat(units='m', limits=(0,500000*dist_per_step,dist_per_step))
    def position_f(self):
        '''Position function, get position, and set position via the general movement function.'''
        # Get your f position
        f_pos = self.get_position_t() #See pyximc.py
        self.lib.get_position(self.device_id, byref(f_pos))
        return f_pos.Position*self.dist_per_step
    @position_f.setter
    def position_f(self,value):
        self.movement_function(value/self.dist_per_step)
    
    # Border state (limit flags), if encoder is used, have it set to BORDER IS ENCODER, otherwise use 0x08 for hard stops
    border_flag_values = {
        'CANNOT GET BORDER STATUS (Likely device not properly initiated)':0x00,
        'BORDER_IS_ENCODER':0x01, 
        'BORDER_STOP_LEFT':0x02,
        'BORDER_STOP_RIGHT':0x04,
        'BORDERS_SWAP_MISSET_DETECTION':0x08,
        'BORDER_STOP_LEFT_AND_RIGHT':0x06,
        'BORDER_STOP_LEFT_AND_RIGHT_AND_ENCODER':0x07
    }   
    @Feat(values=border_flag_values)
    def position_f_limit_flags(self):
        '''Border encoder flags. i.e. limits are either set by encoder (soft limits) or by hardware switches.'''
        self.f_pos_limits = self.edges_settings_t() #See pyximc.py
        self.lib.get_edges_settings(self.device_id, byref(self.f_pos_limits))
        return self.f_pos_limits.BorderFlags
    # @position_f_limit_flags.setter
    # def position_f_limit_flags(self, value):
    #     self.f_pos_limits.BorderFlags = c_uint(value)
    #     self.lib.set_edges_settings(self.device_id, byref(self.f_pos_limits))
    
    # Position limit setter LEFT
    @Feat(units='m', limits=(0,500000*dist_per_step,dist_per_step)) # Real travel range is a bit more, but thats ok. 
    def position_f_limit_left(self):
        '''Soft limits for the LEFT border (SAMPLE SIDE)'''
        self.f_pos_limits = self.edges_settings_t() #See pyximc.py
        self.lib.get_edges_settings(self.device_id, byref(self.f_pos_limits))
        return self.f_pos_limits.LeftBorder*self.dist_per_step
    @position_f_limit_left.setter
    def position_f_limit_left(self,value):
        if value/self.dist_per_step > self.f_pos_limits.RightBorder:
            print('ERROR: Trying to set left border (sample side) to a higher value than the right border.')
        else:
            self.f_pos_limits.LeftBorder = c_int(value/self.dist_per_step)
            self.lib.set_edges_settings(self.device_id, byref(self.f_pos_limits)) 

    # Position limit setter RIGHT
    @Feat(units='m', limits=(0,500000*dist_per_step,dist_per_step)) # Real travel range is a bit more, but thats ok. 
    def position_f_limit_right(self):
        '''Soft limits for the RIGHT border (AWAY FROM SAMPLE SIDE)'''
        self.f_pos_limits = self.edges_settings_t() #See pyximc.py
        self.lib.get_edges_settings(self.device_id, byref(self.f_pos_limits))
        return self.f_pos_limits.RightBorder*self.dist_per_step
    @position_f_limit_right.setter
    def position_f_limit_right(self,value):
        if value/self.dist_per_step < self.f_pos_limits.LeftBorder:
            print('ERROR: Trying to set right border (away-from-sample side) to a lower value than the left border.')
        else:
            self.f_pos_limits.RightBorder = c_int(value/self.dist_per_step)
            self.lib.set_edges_settings(self.device_id, byref(self.f_pos_limits))
    
    # Move settings below, see if changing these is a necessary functionality, if so, simply use set_move_settings with the move_settings_t() structure
    @Feat(limits=(1,20000,1)) #Range: 0..100000.
    def move_speed(self):
        '''Movement speed setter function.'''
        self.move_settings = self.move_settings_t()
        # Get current move settings from controller
        self.lib.get_move_settings(self.device_id, byref(self.move_settings))
        return self.move_settings.Speed
    @move_speed.setter
    def move_speed(self,value):
        self.move_settings.Speed = c_uint(value)
        self.lib.set_move_settings(self.device_id, byref(self.move_settings))

    # Acceleration
    @Feat(limits=(1,10000,1)) # accel Range: 1..65535.
    def move_accel(self):
        '''Acceleration setter function.'''
        self.move_settings = self.move_settings_t()
        # Get current move settings from controller
        self.lib.get_move_settings(self.device_id, byref(self.move_settings))
        return self.move_settings.Accel
    @move_accel.setter
    def move_accel(self, value):
        self.move_settings.Accel = c_uint(value)
        self.lib.set_move_settings(self.device_id, byref(self.move_settings))
    
    # Deceleration
    @Feat(limits=(1,10000,1))
    def move_decel(self):
        '''Deceleration setter function.'''
        self.move_settings = self.move_settings_t()
        # Get current move settings from controller
        self.lib.get_move_settings(self.device_id, byref(self.move_settings))
        return self.move_settings.Decel
    @move_decel.setter
    def move_decel(self, value):
        self.move_settings.Decel = c_uint(value)
        self.lib.set_move_settings(self.device_id, byref(self.move_settings))
    
    # Antiplay speed
    @Feat()
    def move_antiplay_speed(self):
        '''Antiplay speed function display (should not be messed with)'''
        self.move_settings = self.move_settings_t()
        # Get current move settings from controller
        self.lib.get_move_settings(self.device_id, byref(self.move_settings))
        return self.move_settings.AntiplaySpeed

    # # Antiplay steps to take before changing direction
    @Feat()
    def engine_antiplay_steps(self):
        '''Antiplay steps display function display (should not be messed with)'''
        self.engine_settings = self.engine_settings_t()
        # Get current move settings from controller
        self.lib.get_engine_settings(self.device_id, byref(self.engine_settings))
        return self.engine_settings.Antiplay

    # # Engine flags, need to be unpacked, will indicate various engine settings.
    # engine_flag_values = {
    #     'ENGINE_REVERSE':0x01,
	#     'ENGINE_CURRENT_AS_RMS':0x02,
	#     'ENGINE_MAX_SPEED':0x04,
	#     'ENGINE_ANTIPLAY':0x08,
	#     'ENGINE_ACCEL_ON':0x10,
	#     'ENGINE_LIMIT_VOLT':0x20,
	#     'ENGINE_LIMIT_CURR':0x40,
	#     'ENGINE_LIMIT_RPM':0x80,
    # } 
    # @Feat(values=engine_flag_values)
    # def engine_flags(self):
    #     '''Engine status function to display engine flags.'''
    #     self.engine_settings = self.engine_settings_t()
    #     # Get current move settings from controller
    #     self.lib.get_engine_settings(self.device_id, byref(self.engine_settings))
    #     return self.engine_settings.EngineFlags

    # Soft coded step size for move buttons on the driver side
    @Feat(units='m', limits=(1*dist_per_step,500*dist_per_step,dist_per_step))
    def move_step_size(self):
        '''Function to define a typical movement 'step'.'''
        return self.step_size*self.dist_per_step
    @move_step_size.setter
    def move_step_size(self,value):
        self.step_size = value/self.dist_per_step

    # Actions 

    # Soft stop command button
    @Action()
    def stop_soft(self):
        '''Engine soft stop function call'''
        self.lib.command_sstp(self.device_id)    
    
    # Move away from sample by step size
    @Action()
    def move_away_sample_step(self):
        '''Move away from sample by one step'''
        #Moving away from sample means moving to positive values
        self.position_f = self.position_f + Q_(self.step_size*self.dist_per_step, 'm')

    # Move towards sample by step size
    @Action()
    def move_towards_sample_step(self):
        '''Move towards sample by one step'''
        self.position_f = self.position_f - Q_(self.step_size*self.dist_per_step, 'm')
    
    # Move home (disengage)
    @Action()
    def move_home(self):
        '''Move motor home'''
        self.position_f = Q_(self.pos_soft_home*self.dist_per_step, 'm')

    # Move to 'engaged position'
    @Action()
    def move_engaged_pos(self):
        '''Move motor to engaged position'''
        self.position_f = Q_(self.pos_engaged*self.dist_per_step, 'm')

    # Reset to default settings (speed and borders)
    @Action()
    def reset_defaults(self):
        '''Reset default (values I determined) values that you can toy around with'''
        # Speed to 2000 (default)
        self.move_speed = 2000

        # Right border is 300000, left border is 111800
        self.position_f_limit_left = 111800
        self.position_f_limit_right = 300000
        self.f_pos_limits.BorderFlags = c_uint(0x01)
        self.lib.set_edges_settings(self.device_id, byref(self.f_pos_limits)) 

        # Soft Home is 170000
        self.pos_soft_home = 200000 # not that this should ever change, but wtv...
        self.pos_engaged = int(0.0105/self.dist_per_step) #Confirm this value
        print('Values reset to default for the standa speed ({0}), soft home ({1}), and edge settings (L:{2}.uL:{3}, R:{4}.uR:{5}).'.format(self.move_settings.Speed, self.pos_soft_home, self.f_pos_limits.LeftBorder, self.f_pos_limits.uLeftBorder, self.f_pos_limits.RightBorder, self.f_pos_limits.uRightBorder))

    # Shared functions 

    # Movement function to be called for all movement.
    def movement_function(self, move_to_pos):
        '''Main movement function that also serves to minimize stupid errors.'''
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
                print('FocusMotor : Done, position is {0} um.'.format(move_to_pos*self.dist_per_step*10**6))
        else:
            #Print out error is directed location is out of soft bounds.
            print('FocusMotor : MOVEMENT ERROR\nTrying to move to position {0} um from current position {1} um. \nDictated position is out of the movement limits: {2} um and {3} um.'.format(move_to_pos*self.dist_per_step, self.f_pos.Position*self.dist_per_step, self.f_pos_limits.LeftBorder*self.dist_per_step, self.f_pos_limits.RightBorder*self.dist_per_step))

# Testing and validation
if __name__ == '__main__':
    with StandaZFocusDriver() as FocusMotor:
        print('Testing the Standa motor driver...\n')
        print(FocusMotor.idn)
        print('\nIdentification successfull. Driver seems to work.\n')
        print('Current position is {0}'.format(FocusMotor.position_f))
