# -*- coding: utf-8 -*-
'''
    lantz.drivers.newport.FSM_via_nidaqmx
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    :copyright:
    C.Egerstrom - Dec 2022
    De-lantz-ed Nazar's stuff. Use as you wish, but only for evil.
    N.Delegan - Dec 2020.
    Heavily adapted from Uri and Jacob's code; debugged with Michael. Use as you wish, be nice.

    :description:
    Driver for our newport FSM via the generic DAQ motion controller driver.
    
    :versions/changelog: 
    * V1.0 - Dec 2020 -  Driver created

    :dependencies: 
    * lantz
    * ni_motion_controller, ensure it can be imported properly
    * Q_ (units)
	
	:usage:
'''

# Import lantz dependencies
from lantz import Driver, Feat, DictFeat, Action, Q_

# Import ni_motion_controller; Note: On jasper there is a weird pathing thing, had to hard add the path.
try:
    import sys
    sys.path.append(r'C:\Users\adminnd\code\jasper-repo') # TODO: Remove this once Jasper is fixed and import properly

    # Import the generic ni motion controller driver
    from drivers.ni.ni_motion_controller import NIDAQMotionController, NIDAQAxis # For some reason this way of importing does not work...
except:
    RuntimeWarning('Unable to load the ni_motion_controller.py driver. Check the file import path and syntax.')

# Driver class that is set-up specific
class FSM_2D(Driver):
    # Default values from calibration values in um/V (Calibrated on 18th dec 2020)
    def __init__(self, x_ch='DAQ1/ao0', y_ch='DAQ1/ao1', ctr_ch='DAQ1/ctr0', XperV = 10.1665, YperV = 7.3414):
        # Create axis objects; limits are -10 to 10V
        x_axis = NIDAQAxis(x_ch, 'um', Q_(1/XperV, 'V/um'), limits=(Q_(-XperV*10,'um'), Q_(XperV*10,'um')))
        y_axis = NIDAQAxis(y_ch, 'um', Q_(1/YperV, 'V/um'), limits=(Q_(-YperV*10,'um'), Q_(YperV*10,'um')))

        # Create the daq object # TODO: check units here, can acquisition be changed? smoothing?
        self.daq_controller = NIDAQMotionController(ctr_ch, Q_(20, 'kHz'), {'x': x_axis, 'y': y_axis,}, ao_smooth_steps=Q_(5000, '1/V'))

        # Set increments # TODO: Any reason to ever change this?
        self.increment_X = "Q_({},'nm')".format(XperV*10) # Note the units
        self.increment_Y = "Q_({},'nm')".format(YperV*10) # Note the units

    # Initialization function, call the daq_controller initialize
    def initialize(self):
        return self.daq_controller.initialize()
        
    # Finalize function, call the daq_controller finalize
    def finalize(self):
        return self.daq_controller.finalize()
    
    # Feats
    # Acquisition rate setter
    @Feat(units='Hz')
    def acq_rate(self):
        return self.daq_controller.acq_rate
    @acq_rate.setter
    def acq_rate(self, acq_rate):
        self.daq_controller.acq_rate = 1

    #TODO: Incorporate once dictfeat are fixed.
    # @DictFeat(keys = ['x', 'y'], units='um')
    # def per_channel_position(self, key):
    #     return self.daq_controller.position[key]
    # @per_channel_position.setter
    # def per_channel_position(self, key, pos):
    #     # TODO: Fix this to make it more elegant, i.e. is keys fed to the subfunction?
    #     keys = ['x', 'y']
    #     keys.remove(key)
    #     self.daq_controller.move({key: pos, keys: self.daq_controller.position[keys]})

    # Position setter X
    @Feat(units='um')
    def x(self):
        # print(self.daq_controller.position['x'])
        return self.daq_controller.position['x']
    @x.setter
    def x(self, pos):
        # print('x set {}'.format(pos))
        self.daq_controller.move({'x': pos, 'y': self.y})

    # Position setter Y
    @Feat(units='um')
    def y(self):
        return self.daq_controller.position['y']
    @y.setter
    def y(self, pos):
        self.daq_controller.move({'x': self.x, 'y': pos})

    # Actions
    # Create new counter channel TODO: check if thse are thread safe.
    @Action()
    def new_ctr_task(self, ctr_ch):
        return self.daq_controller.new_ctr_task(ctr_ch)

    # Move to point function
    @Action()
    def move(self, point):
        return self.daq_controller.move(point)

    # Line scan
    @Action()
    def line_scan(self, init_point, final_point, steps, pts_per_step):
        return self.daq_controller.line_scan(init_point, final_point, steps, pts_per_step)

    # @Action()
    # def twoD_scan(self, init_point, final_point, steps, pts_per_step):
    #     # TODO: Create 2D scan FSM scan function?
    #     return None

# For testing
if __name__ == '__main__':
    # Importing wait
    import time

    print('Testing the FSM driver...\n')
    with FSM_2D() as fsm:
        print('\nTesting FSM.\n')

        # Just some random test points to see if things move.
        point0 = {'x': Q_(0, 'um'), 'y': Q_(0, 'um')}
        point1 = {'x': Q_(10, 'um'), 'y': Q_(10, 'um')}
        point2 = {'x': Q_(-10, 'um'), 'y': Q_(-10, 'um')}

        # Some 'unit testing'.
        print(f'Moving to point {point1}.')
        fsm.daq_controller.move(point1)
        print(f'Chilling here for a second.')
        time.sleep(1)
        print(f'Moving to point {point2}.')
        fsm.daq_controller.move(point2)
        print(f'Chilling here for a second.')
        time.sleep(1)
        print(f'Moving back to point {point0}.')
        fsm.daq_controller.move(point0)
        print('\nIf nothing yelled at you, the driver is likely fine. FSM driver works.\n')

        None


