# -*- coding: utf-8 -*-
'''
    lantz.drivers.cobolt.cobolt0401
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    :copyright: 
    N.Delegan - Sep 2020.
    Use as you wish, be nice.

    :description:
    Driver for Cobolt 04-01 series laser
    
    :versions/changelog: 
    * V1.0 - Sep 2020 -  Driver written, not exhaustively tested.

    :dependencies: 
    * pyvisa.constants
    * lantz Action, Feat, & MessageBasedDriver
    * time sleep
	
	:usage:
'''

from pyvisa import constants
from lantz import Action, Feat, MessageBasedDriver
from time import sleep

class Cobolt0401(MessageBasedDriver):
    '''Driver class for any Cobolt 04-01 Series laser.
    '''
	
	#def __init__:
	#	pass
	
	#kwargs for MessageBasedDriver
    DEFAULTS = {'ASRL': {'write_termination': '\r',
                         'read_termination': '\r',
                         'baud_rate': 115200,
                         # 'bytesize': 8,
                         'parity': constants.Parity.none,
                         'stop_bits': constants.StopBits.one,
                         'encoding': 'ascii',
                         }}
    
    #LASER OPERATION
    @Feat(read_once=True)
    def idn(self):
        '''Get serial number query output as string.
        *query*: 'sn?' returns a 32-bit unsigned integer
        '''
        return 'Cobolt 04-01 Series Laser S/N: ' + self.query('sn?')[1:]

    @Feat(values={True: '1', False: '0'})
    def directcontrol_enabled(self):
        '''Direct control enable state
        *query*: '@cobasdr?' returns 0 = disabled or 1 = enabled
        *command*: '@cobardr' to set 0 = disabled or 1 = enabled
        **Set capability is for OEM models only**
        '''
        #Get direct control enable state
        ans = self.query('@cobasdr?')
        return ans[1:]
    @directcontrol_enabled.setter
    def directcontrol_enabled(self, value):
        self.query('@cobasdr ' + value)

    @Feat(values={True: '1', False: '0'})
    def enabled(self):
        '''Method for laser status inquery
        and turning the laser ON/OFF. Control requires autostart disabled.
        *query*: 'l?' returns 0 = OFF or 1 = ON
        *command*: 
        *'l1' Laser ON, requires autostart disabled. (OEM models)
        *'l0' Laser OFF. (OEM models)
        '''
        return self.query('l?')[-1]
    @enabled.setter
    def enabled(self, value):
        self.query('l' + value)

    @Feat(values={True: '1', False: '0'})
    def autostart(self):
        '''Autostart handling
        *query*: '@cobas?' returns 0 = disabled or 1 = enabled
        *command*: 
        *'@cobas 1' Autostart ON.
        *'@cobas 0' Autostart OFF. (Diable only on OEM models 
        as this will stop the laser from going through the proper warm up routine)
        '''
        ans = self.query('@cobas?')
        return ans[-1]
    @autostart.setter
    def autostart(self, value):
        self.query('@cobas ' + value)

    # LASER INFORMATION METHODS

    @Feat()
    def operating_hours(self):
        '''Get Laser Head operating hours
        *query*: 'hrs?' returns float
        '''
        return self.query('hrs?')[1:]

    @Feat(values={'Interlock open': '1', 'OK': '0'})
    def interlock(self):
        '''Get interlock state
        *query*: 'ilk?' returns 0 = OK or 1 = remove interlock open
        (See pin diagram behind controller for more info)
        '''
        return self.query('ilk?')[1:]

    # LASER'S SET POINT

    @Feat(units='A',limits=(0, 4.2))
    def current_sp(self):
        '''Get drive current and set drive current
        *query*: 'i?' returns a float in amperes
        *command*: 
        *'slc [amperes in float]'
        '''
        return float(self.query('i?'))
    @current_sp.setter
    def current_sp(self, value):
    	#Set output current, max is 4.2 A
        print('Laser in constant power mode.\n Current setting not recommended.')
        #self.query('slc {:.1f}'.format(value))

    @Feat(units='mW',limits=(0, 100))
    def power_sp(self):
        '''Laser power control (laser receives and outputs in W)
        *query*: 'p?' get set output power as a float (in W)
        *command*: 
        *'p [watts in float]' set output power
        '''
        return 1000 * float(self.query('p?'))
        #Power setter, between 0 and 100 mW
    @power_sp.setter
    def power_sp(self, value):
        self.query('p {:.5f}'.format(value / 1000))
        
    @Feat(units='mW')
    def power(self):
        '''Read actual output power
        *query*: 'pa?' get actual output power as a float (in W)
        '''
        return 1000 * float(self.query('pa?'))

    @Feat(values={'Temperature error': 1, 'No errors': 0,
                  'Interlock error': 3, 'Constant power time out': 4})
    def operating_fault(self):
        '''Get operating fault
        *query*: 'f?' outputs are 
        *0 = no fault
        *1 = temperature error
        *3 = open interlock
        *4 = constant power fault
        '''
        return int(self.query('f?')[-1])

    @Feat()
    def LED_status(self):
        '''Get display LED status (from controller)
        *query*: 'leds?' output is an int [0:15] with bit = 1/0 being that LED is on/off
        *bit 0 = power on
        *bit 1 = laser on
        *bit 2 = laser lock
        *bit 3 = error
        '''
        LED_str = format(int(self.query('leds?')), '04b')
        pwr_LED = 'Power LED ('+LED_str[3]+') | '
        lsr_LED = 'Laser LED ('+LED_str[2]+') | '
        lsrlock_LED = 'Lock LED ('+LED_str[1]+') | ' 
        err_LED = 'Error LED ('+LED_str[0]+')'
        return pwr_LED + lsr_LED + lsrlock_LED + err_LED
    
    #Actions and finalize

    @Action()
    def clear_fault(self):
        '''Clear fault command
        *command*: 'cf'
        '''
        self.query('cf')

    @Action()
    def restart(self):
        '''Forces the laser on without checking if autostart is enabled.
        '''
        self.query('@cob1')

    @Action()
    def laser_start(self):
        '''Start the laser properly
        The command clear faults, ensures autostart is enables, then starts the laser.
        A timed loop gives the laser time to warm up before other commands are inputted.
        '''
        if self.interlock == '0':
            print('Interlock is off, cannot start laser.')
        else:
            #Clear faults
            self.query('cf')
            #Laser ON after interlock check
            self.query('@cob1')
            print('Please wait while laser warms up and starts.\n')
            while self.power_sp < 0.95*self.power: #self.query('l?')[-1] == '0' & 
                #Wait 5 seconds before checking laser status again
                sleep(5)
            else:
                print('Laser is now ON.\n Please wait 1-2 minutes for power to stabilize.')

    @Action()
    def laser_stop(self):
        '''Stop laser command'''
        self.query('l0')
        print('Laser turned off.\n')

    def finalize(self):
        '''Finalize command, clear faults and turn laser off.'''
        #Laser off
        # self.query('l0')
        #Clear faults
        self.query('cf')

if __name__ == '__main__':
    print('Testing the Cobolt0401 laser driver...\n')
    with Cobolt0401('COM3') as laser:
        print(laser.idn)
        print('\nIdentification successfull. Driver seems to work.\n')