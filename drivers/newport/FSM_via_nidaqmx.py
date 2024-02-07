# -*- coding: utf-8 -*-
'''
    drivers.newport.FSM_via_nidaqmx
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    :copyright:
    C.Egerstrom - Dec 2022
    De-lantz-ed version of Nazar's old driver. Use as you wish, be evil.

    N.Delegan - Dec 2020.
    Heavily adapted from Uri and Jacob's code; debugged with Michael. Use as you wish, be nice.

    :description:
    Driver for our newport FSM via the generic DAQ motion controller driver.
    
    :versions/changelog: 
    * V1.0_OLD - Dec 2020 -  Driver created by Nazar
    * V1.1 - Dec 2022 -  Driver de-Lantz-ed (WIP-needs to be verified and debugged)
                         Uses 'x','y' dicts for points, might just make them tuples/list if they annoy me - Chris
    * V1.2 - Mar 2023 - Drivers previously verified, but added functionality for move fxns to take (x,y) tuples as well as {'x': x,'y':y} dicts
                        TODO: Check everything still works, rel moves work, and 1D and 2D scans

    :dependencies: 
    * ni_motion_controller, ensure it can be imported properly
    * numpy (because Chris is lazy, can get rid of it if needed)
	
	:usage:
'''

# Import np for convenience
import numpy as np

#jsut to deal with some weirdness that comes from passing dicts as rpyc netrefs in new Nspyre
from rpyc.utils.classic import obtain

# Import ni_motion_controller; Note: On jasper there is a weird pathing thing, had to hard add the path.
try:
    import sys
    sys.path.append(r'C:\\Users\\adminnd\\code\\jasper-repo\\drivers\\') # TODO: Remove this once Jasper is fixed and import properly

    # Import the generic ni motion controller driver
    from ni.ni_motion_controller import NIDAQMotionController, NIDAQAxis #Maybe figured this out??? -CE 020123
except:
    RuntimeWarning('Unable to load the ni_motion_controller.py driver. Check the file import path and syntax.')



# Driver class that is set-up specific
class FSM_2D():
    DEFAULT_UNITS_DISTANCE = 'um' #these are never used, just want them to be standard and conveniently accessible
    DEFAULT_UNITS_RATE = 'Hz'


    # Default values from calibration values in um/V (Recalibrated on 8th Mar 2023)
    def __init__(self, x_ch='Dev1/ao0', y_ch='Dev1/ao1', ctr_ch='Dev1/ctr0', XperV = 10.9472683, YperV = 7.37113658): #Was DAQ1, not Dev1, but this allowed a move command to work
        # Create axis objects; limits are -10 to 10V
        self.x_axis = NIDAQAxis(x_ch, 1/XperV, limits=(-XperV*10, XperV*10))
        self.y_axis = NIDAQAxis(y_ch, 1/YperV, limits=(-YperV*10, YperV*10))

        self.ctr_ch = ctr_ch


    # Initialization function, call the daq_controller initialize
    def __enter__(self):
        '''Create the daq object''' # TODO: check units here, can acquisition be changed? smoothing?
        self.daq_controller = NIDAQMotionController(self.ctr_ch, 20e3, {'x': self.x_axis, 'y': self.y_axis,})#, ao_smooth_steps=5000) #already default smooth arg
        print('FSM Connected')
        self.daq_controller.__enter__() #TODO does returning this vs self matter?
        return self


    # Finalize function, call the daq_controller finalize
    def __exit__(self, *args):
        self.daq_controller.__exit__()


    def setAcqRate(self, newAcqRate):
        '''Sets the new acq rate (in Hz) to newAcqRate
        Arguments:  *newAcqRate'''
        self.daq_controller.acq_rate = newAcqRate


    def getAcqRate(self):
        '''Returns the acq rate (in Hz)'''
        return(self.daq_controller.acq_rate)


    def getPosFromAxis(self, axName):
        '''Returns position (in um) for given axis
        Arguments:  *axName, possible values: 'x','y' '''
        if axName not in ('x', 'y'):
            raise KeyError("Axis must be 'x' or 'y'")
        return self.daq_controller.position[axName]


    def getPos(self):
        '''Returns the current (x,y) coord (in um) as tuple'''
        return(self.getPosFromAxis('x'), self.getPosFromAxis('y'))


    def convertPtLikeToPt(self, point):
        '''Convert given point to an \{'x':x,'y':y\}) dict
        Arguments:  *point, a len=2 list, tuple, or {'x': x, 'y': y} dict'''
        point = obtain(point)
        if hasattr(point, '__len__') and len(point)==2:
            if isinstance(point, dict): #Not checking for dict type because rpyc screws up dicts and passes them as netrefs
                return(point)
            elif isinstance(point, list) or isinstance(point, tuple):
                return( {'x': point[0], 'y': point[1]} )
        raise ValueError(f"Given point {point} is not valid (i.e. not a len=2 list, tuple, or 'x':x,'y':y) dict")


    def move(self, point):
        '''Moves to a given (x,y) coord (in um)
        Arguments:  *point, either as len=2 list, tuple, or {'x': x, 'y': y} dict'''
        point = self.convertPtLikeToPt(point)
        return(self.daq_controller.move(point))


    def relMove(self, relPoint):
        '''Moves to a point a given (x,y) dist away (in um)
        Arguments:  *relPoint, either as len=2 list, tuple, or {'x': x, 'y': y} dict'''
        relPoint = self.convertPtLikeToPt(relPoint)
        return(self.daq_controller.move({'x': self.getPosFromAxis('x')+relPoint['x'], 'y': self.getPosFromAxis('y')+relPoint['y']}))


    def relMoveX(self, xDist):
        '''Moves a given xDist in x (in um)
        Arguments:  *xDist'''
        return(self.daq_controller.move({'x': self.getPosFromAxis('x')+xDist, 'y': self.getPosFromAxis('y')}))


    def relMoveY(self, yDist):
        '''Moves a given yDist in y (in um)
        Arguments:  *yDist'''
        return(self.daq_controller.move({'x': self.getPosFromAxis('x'), 'y': self.getPosFromAxis('y')+yDist}))


    # TODO: check if thse are thread safe.
    def new_ctr_task(self, ctr_ch):
        '''Create new counter channel
        Arguments:  *ctr_ch, a counter hardware address like DAQ1/ctr0 (default ctr_ch during init)'''
        return(self.daq_controller.new_ctr_task(ctr_ch))


    def line_scan(self, init_point, final_point, steps, pts_per_step):
        '''Line scans between init_point and final_point with a number of steps and pts_per_step
        Arguments:  *init_point, either as len=2 list, tuple, or {'x': x, 'y': y} dict
                    *final_point, either as len=2 list, tuple, or {'x': x, 'y': y} dict
                    *pts_per_step, TODO: figure out what this is. Number of points to scan???'''
        init_point = self.convertPtLikeToPt(init_point)
        final_point = self.convertPtLikeToPt(final_point)
        return(self.daq_controller.line_scan(init_point, final_point, steps, pts_per_step))


    def twoD_scan(self, init_point, final_point, steps, pts_per_step):
        '''2D scans between init_point and final_point with a number of steps (can be a number or 
        len=2 list for x/y number of steps) and pts_per_step. Returns data in the form:
        [(x_i,y_i) (x_i,y_1) ... (x_i,y_n-2) (x_i,y_f)
         (x_1,y_i) (x_1,y_1) ... (x_1,y_n-2) (x_1,y_f)
         ...
         (x_n-2,y_i)     ... (x_n-2,y_n-2) (x_n-2,y_f)
         (x_f,y_i)       ... (x_f,y_n-2) (x_f,y_f) ]
         '''
        init_point = self.convertPtLikeToPt(init_point)
        final_point = self.convertPtLikeToPt(final_point)

        if hasattr(steps, '__len__'):
            if len(steps)==2:
                if type(steps) is dict:
                    xNumSteps = steps['x']
                    yNumSteps = steps['y']
                elif type(steps) is list or tuple:
                    xNumSteps = steps[0]
                    yNumSteps = steps[1]
            else:
                raise IndexError('If passing a list or dict of steps to twoD_scan, must have two (x,y) values')
        else:
            #print(steps) #DEBUG
            xNumSteps = steps
            yNumSteps = steps

        twoDavgdAndNormed = []
        xSteps = np.linspace(init_point['x'], final_point['x'], xNumSteps)
        yStart = init_point['y']
        yStop = final_point['y']
        for xStep in xSteps:
            twoDavgdAndNormed.append(self.daq_controller.line_scan( {'x':xStep, 'y':yStart}, {'x':xStep, 'y':yStop}, yNumSteps, pts_per_step))
        return(np.array(twoDavgdAndNormed))



# For testing
if __name__ == '__main__':
    # Importing wait
    import time

    print('Testing the FSM driver...\n')
    with FSM_2D() as fsm:
        print('\nTesting FSM.\n')

        advTesting = False
        if advTesting:
            data = fsm.twoD_scan({'x':-5,'y':-5}, {'x':5,'y':5}, 5, int(fsm.getAcqRate()))
            print(data)
        else:
            # Just some random test points (in um) to see if things move.
            point0 = {'x': 0, 'y': 10}
            point1 = {'x': 10, 'y': 10}
            point2 = {'x': -10, 'y': -10}

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


