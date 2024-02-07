"""
C.Egerstrom's first foray into creating a GUI to control XYZ motion with new nspyre.

Copyright (c) Mar. 2023, C.Egerstrom
All rights reserved.

This work is licensed under the terms of the 3-Clause BSD license.
For a copy, see <https://opensource.org/licenses/BSD-3-Clause>.
"""

from nspyre import InstrumentGateway


class JasperMotionController():
    """Allows for XYZ positioning on the Jasper setup.
    Arguments:  *testing, where it doesn't connect to an InstrumentGateway or actually do anything. Default=False """

    def __init__(self, testing=False):
        self.testing=testing
        if not testing:
            with InstrumentGateway() as gw:
                try:
                    gw.kim
                    gw.standa
                    gw.fsm
                except Exception as e:
                    print('Missing a device:', e)


    def relMoveStage(self, axis, dist, verbose=False):
        '''If controller isn't running in testing mode, will move a specified dist along a specified 
        stage axis (using KimXY for XY and Standa DC motor for Z)
        Arguments:  *axis, one of 'x','y','z'
                    *dist, distance to move (in mm)
                    *verbose, whether to print when moving. Default=False'''
        if not self.testing:
            with InstrumentGateway() as gw:
                if axis=='x':
                    if verbose:
                        print(f'Moving X stage {dist} mm')
                    gw.kim.relMoveX(dist)
                elif axis=='y':
                    if verbose:
                        print(f'Moving Y stage {dist} mm')
                    gw.kim.relMoveY(dist)
                elif axis=='z':
                    if verbose:
                        print(f'Moving Z stage {dist} mm')
                    gw.standa.relMove(dist)
                else:
                    raise AssertionError("Passed an axis that wasn't x,y,z")


    def relMoveFSM(self, axis, dist, verbose=False):
        '''If controller isn't running in testing mode, will move the FSM a specified dist
        Arguments:  *axis, one of 'x','y'
                    *dist, distance to move (in um)
                    *verbose, whether to print when moving. Default=False'''
        if not self.testing:
            with InstrumentGateway() as gw:
                if axis=='x':
                    if verbose:
                        print(f'Moving X FSM {dist} um')
                    gw.fsm.relMoveX(dist)
                elif axis=='y':
                    if verbose:
                        print(f'Moving Y FSM {dist} um')
                    gw.fsm.relMoveY(dist)
                else:
                    raise AssertionError("Passed an axis that wasn't x,y,z")


    def moveFSMhome(self, verbose=False):
        '''If controller isn't running in testing mode, will move the FSM to (0,0)
        Arguments:  *verbose, whether to print when moving. Default=False'''
        if not self.testing:
            with InstrumentGateway() as gw:
                if verbose:
                    print(f'Moving FSM to (0,0)')
                gw.fsm.move( (0,0) )


    '''def absMove(self, axis, dist,):
        time.sleep(1) #TEMP
        if axis=='x':
            pass
        elif axis=='y':
            pass
        elif axis=='z':
            pass
        else:
            raise AssertionError("Passed an axis that wasn't x,y,z")'''

    '''def getPos(self, axis):
        time.sleep(1) #TEMP
        if axis=='x':
            pass
        elif axis=='y':
            pass
        elif axis=='z':
            pass
        else:
            raise AssertionError("Passed an axis that wasn't x,y,z")'''

if __name__ == '__main__':
        import time #if running with testing=False, de-indent this line and un-comment next line 
    #with InstrumentGateway() as gw:
        myMtnCtlr = JasperMotionController(testing=True)
        
        print('Moving FSM')
        myMtnCtlr.relMoveFsm('x', 10)
        time.sleep(1)
        myMtnCtlr.relMoveFsm('x', -10)
        time.sleep(1)
        myMtnCtlr.relMoveFsm('y', 10)
        time.sleep(1)
        myMtnCtlr.relMoveFsm('y', -10)
        time.sleep(1)
        print('FSM done, now moving stage')
        myMtnCtlr.relMoveStage('x', 0.01)
        time.sleep(1)
        myMtnCtlr.relMoveStage('x', -0.01)
        time.sleep(1)
        myMtnCtlr.relMoveStage('y', 0.01)
        time.sleep(1)
        myMtnCtlr.relMoveStage('y', -0.01)
        time.sleep(1)
        myMtnCtlr.relMoveStage('z', 0.01)
        time.sleep(1)
        myMtnCtlr.relMoveStage('z', -0.01)
        time.sleep(1)
        print('Done moving stage')