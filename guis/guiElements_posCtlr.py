#BETTER LOCK-OUT TO PREVENT MUTLIPLE MOVES, AND AN INDICATOR FOR IF IT'S STILL MOVING

"""
GUI for a Motion Controller Application

Copyright (c) April 2023, Chris Egerstrom
All rights reserved.

This work is licensed under the terms of the 3-Clause BSD license.
For a copy, see <https://opensource.org/licenses/BSD-3-Clause>.
"""

from pyqtgraph.Qt import QtWidgets
from pyqtgraph import SpinBox
from nspyre import ParamsWidget
from experiments.xyzPosCtlr import JasperMotionController


class JasperMotionControllerGUI(QtWidgets.QWidget):
    
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Jasper Motion Controller')

        self.isMoving = False
        self.motionCtlr = JasperMotionController()

        self.xyStageParams_widget = ParamsWidget(
            {
                'xStageStep': {
                    'display_text': 'X Step',
                    'widget': SpinBox(
                        value = 0.01,
                        suffix = 'mm',
                        siPrefix = False,
                        bounds = (0, 10),
                        dec = True,
                    ),
                },
                'yStageStep': {
                    'display_text': 'Y Step',
                    'widget': SpinBox(
                        value = 0.01,
                        suffix = 'mm',
                        siPrefix = False,
                        bounds = (0, 10),
                        dec = True,
                    ),
                },
            }
        )
        self.zStageParams_widget = ParamsWidget(
            {
                'zStageStep': {
                    'display_text': 'Z Step',
                    'widget': SpinBox(
                        value = 0.01,
                        suffix = 'mm',
                        siPrefix = False,
                        bounds = (0, 1),
                        dec = True,
                    ),
                },
            }
        )
        self.xyFsmParams_widget = ParamsWidget(
            {
                'xFSMstep': {
                    'display_text': 'X Step',
                    'widget': SpinBox(
                        value = 1.0,
                        suffix = 'um',
                        siPrefix = False,
                        bounds = (0, 30),
                        dec = True,
                    ),
                },
                'yFSMstep': {
                    'display_text': 'Y Step',
                    'widget': SpinBox(
                        value = 1.0,
                        suffix = 'um',
                        siPrefix = False,
                        bounds = (0, 30),
                        dec = True,
                    ),
                }
            }
        )

        #Create and connect buttons to argument-less functions that should be run when they're clicked
        upXstageButton = QtWidgets.QPushButton('+X')
        upXstageButton.clicked.connect(lambda a: self.moveButtonClicked('stage','x',-1))
        downXstageButton = QtWidgets.QPushButton('-X')
        downXstageButton.clicked.connect(lambda a: self.moveButtonClicked('stage','x',1))

        upYstageButton = QtWidgets.QPushButton('+Y')
        upYstageButton.clicked.connect(lambda a: self.moveButtonClicked('stage','y',1))
        downYstageButton = QtWidgets.QPushButton('-Y')
        downYstageButton.clicked.connect(lambda a: self.moveButtonClicked('stage','y',-1))

        upZstageButton = QtWidgets.QPushButton('+Z')
        upZstageButton.clicked.connect(lambda a: self.moveButtonClicked('stage','z',1)) #NOT SURE IF Z IS RIGHT, BUT OTHERWISE THESE ARE
        downZstageButton = QtWidgets.QPushButton('-Z')
        downZstageButton.clicked.connect(lambda a: self.moveButtonClicked('stage','z',-1)) #INTUITIVE WITH RESPECT TO THE CAMERA

        upXfsmButton = QtWidgets.QPushButton('+X')
        upXfsmButton.clicked.connect(lambda a: self.moveButtonClicked('fsm','x',-1))
        downXfsmButton = QtWidgets.QPushButton('-X')
        downXfsmButton.clicked.connect(lambda a: self.moveButtonClicked('fsm','x',1))

        upYfsmButton = QtWidgets.QPushButton('+Y')
        upYfsmButton.clicked.connect(lambda a: self.moveButtonClicked('fsm','y',-1))
        downYfsmButton = QtWidgets.QPushButton('-Y')
        downYfsmButton.clicked.connect(lambda a: self.moveButtonClicked('fsm','y',1))

        homeFsmButton = QtWidgets.QPushButton('Home') #Home button to move FSM to (0,0)
        homeFsmButton.clicked.connect(self.fsmHomeButtonClicked)

        #create labels for sets of move buttons
        xyStageLabel = QtWidgets.QLabel()
        xyStageLabel.setText('XY Stage')
        zStageLabel = QtWidgets.QLabel()
        zStageLabel.setText('Z Stage')
        xyFsmLabel = QtWidgets.QLabel()
        xyFsmLabel.setText('FSM')

        #setup a grid and fill in with buttons to get a pretty control panel
        moveLayout = QtWidgets.QGridLayout()
        moveLayout.addWidget(xyStageLabel, 0, 1)
        moveLayout.addWidget(upXstageButton, 2, 2)
        moveLayout.addWidget(downXstageButton, 2, 0)
        moveLayout.addWidget(upYstageButton, 1, 1)
        moveLayout.addWidget(downYstageButton, 3, 1)
        moveLayout.addWidget(self.xyStageParams_widget, 4, 1)

        moveLayout.addWidget(zStageLabel, 0, 4)
        moveLayout.addWidget(upZstageButton, 1, 4)
        moveLayout.addWidget(downZstageButton, 3, 4)
        moveLayout.addWidget(self.zStageParams_widget, 4, 4)
        
        moveLayout.addWidget(xyFsmLabel, 0, 7)
        moveLayout.addWidget(upXfsmButton, 2, 8)
        moveLayout.addWidget(downXfsmButton, 2, 6)
        moveLayout.addWidget(upYfsmButton, 1, 7)
        moveLayout.addWidget(downYfsmButton, 3, 7)
        moveLayout.addWidget(homeFsmButton, 2, 7)
        moveLayout.addWidget(self.xyFsmParams_widget, 4, 7)

        #keep empty columns from collapsing to have some spacing
        moveLayout.setColumnMinimumWidth(3, 50)
        moveLayout.setColumnMinimumWidth(5, 50)

        #make layout active
        self.setLayout(moveLayout)


    def moveButtonClicked(self, devType, axis, dir):
        '''Flexible internal function to call with particular args when a move button is pressed
        Arguments:  *devType: one of 'stage','fsm'
                    *axis: one of 'x','y','z'
                    *dir: one of +1,-1 '''
        if not self.isMoving: #should lock out user out of making multiple moves at once, but doesn't actually work #TODO: Make it
            self.isMoving = True

            if devType == 'stage':
                if axis == 'x': 
                    self.motionCtlr.relMoveStage(axis, dir*self.xyStageParams_widget.xStageStep)
                elif axis == 'y': 
                    self.motionCtlr.relMoveStage(axis, dir*self.xyStageParams_widget.yStageStep)
                elif axis == 'z': 
                    self.motionCtlr.relMoveStage(axis, dir*self.zStageParams_widget.zStageStep)
                else:
                    raise AssertionError('Passed an axis not in x,y,z')
                
            elif devType == 'fsm':
                if axis == 'x': 
                    self.motionCtlr.relMoveFSM(axis, dir*self.xyFsmParams_widget.xFSMstep)
                elif axis == 'y': 
                    self.motionCtlr.relMoveFSM(axis, dir*self.xyFsmParams_widget.yFSMstep)
                else:
                    raise AssertionError('Passed an axis not in x,y')
            else:
                raise AssertionError('Passed an invalid devType')
            self.isMoving = False


    def fsmHomeButtonClicked(self):
        '''Wrapper for motion controller's moveFSMhome method to call when relevant button is clicked '''
        if not self.isMoving:
            self.isMoving = True
            self.motionCtlr.moveFSMhome()
            self.isMoving = False