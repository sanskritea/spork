"""
GUI for a Fsm Scan Application

Copyright (c) May 2023, Chris Egerstrom
All rights reserved.

This work is licensed under the terms of the 3-Clause BSD license.
For a copy, see <https://opensource.org/licenses/BSD-3-Clause>.
"""

#TO-DO, NEEDS SOME DEBUGGING
from functools import partial
from importlib import reload

from pyqtgraph.Qt import QtWidgets
from pyqtgraph import SpinBox
from nspyre import ParamsWidget
from nspyre import ProcessRunner
from nspyre import DataSink
from nspyre.gui.widgets.heatmap import HeatMapWidget

import experiments.fsmScan as fsmScan
from guis.guiElements_general import flexSave


class CustomFsmScanWidget(QtWidgets.QWidget):
    
    def __init__(self):
        super().__init__()

        self.setWindowTitle('FSM Scan')

        self.params_widget = ParamsWidget({
            'xCenter': {
                'display_text': 'X Center Coord',
                'widget': SpinBox(
                    value=0,
                    suffix='m',
                    siPrefix=True,
                    bounds=(-50e-6, 50e-6),
                    dec=True,
                ),
            },
            'yCenter': {
                'display_text': 'Y Center Coord',
                'widget': SpinBox(
                    value=0,
                    suffix='m',
                    siPrefix=True,
                    bounds=(-50e-6, 50e-6),
                    dec=True,
                ),
            },
            'xRange': {
                'display_text': 'X Range (one-sided)',
                'widget': SpinBox(
                    value=5e-6,
                    suffix='m',
                    siPrefix=True,
                    bounds=(-50e-6, 50e-6),
                    dec=True,
                ),
            },
            'yRange': {
                'display_text': 'Y Range (one-sided)',
                'widget': SpinBox(
                    value=5e-6,
                    suffix='m',
                    siPrefix=True,
                    bounds=(-50e-6, 50e-6),
                    dec=True,
                ),
            },
            'xPts': {
                'display_text': 'Number of X Pts',
                'widget': SpinBox(
                    value=10,
                    int=True,
                    dec=True,
                ),
            },
            'yPts': {
                'display_text': 'Number of Y Pts',
                'widget': SpinBox(
                    value=10,
                    int=True,
                    dec=True,
                ),
            },
            'collectsPerPt': {
                'display_text': 'Points per Position',
                'widget': SpinBox(
                    value=100,
                    int=True,
                    dec=True,
                ),
            },
            'datasetName': {
                'display_text': 'Dataset Name',
                'widget': QtWidgets.QLineEdit('FsmScan'),
            },
        })

        #Setup run and stop buttons
        runButton = QtWidgets.QPushButton('Run')
        runButton.clicked.connect(self.runClicked)

        stopButton = QtWidgets.QPushButton('Stop')
        stopButton.clicked.connect(self.stop)

        # stop the process if the widget is closed
        self.destroyed.connect(partial(self.stop))

        # the process running the sweep function
        self.sweepProc = ProcessRunner()

        # Qt layout that arranges the params, checkboxes, and buttons vertically
        params_layout = QtWidgets.QVBoxLayout()
        params_layout.addWidget(self.params_widget)
        params_layout.addStretch()
        params_layout.addWidget(runButton)
        params_layout.addWidget(stopButton)
        self.setLayout(params_layout)


    def runClicked(self):
        """Runs when the 'run' button is pressed."""

        # reload the spin measurements module at runtime in case any changes were made to the code
        reload(fsmScan)

        # create an instance of the ODMR class that implements the experimental logic.
        fsmScanMeas = fsmScan.FSMScanMeasurement()

        # run the sweep function in a new thread
        self.sweepProc.run(
            fsmScanMeas.fsmScan,
            self.params_widget.datasetName,
            1e6*self.params_widget.xCenter, #convert m=>um since FSM default unit is um
            1e6*self.params_widget.yCenter,
            2e6*self.params_widget.xRange, #1-sided to 2-sided ranges
            2e6*self.params_widget.yRange,
            self.params_widget.xPts,
            self.params_widget.yPts,
            self.params_widget.collectsPerPt
        )


    def stop(self):
        """Stop the sweep process."""
       #kill the TaskVsTime sweep process
        self.sweepProc.kill()
        flexSave(self.params_widget.datasetName, 'fsmScan', 'closeout') #if process is closed



class FsmScanPlotWidget(HeatMapWidget):
    def __init__(self):
        title = 'FSM Scan'
        super().__init__(title=title, btm_label='X', lft_label='Y')


    def setup(self):
        self.sink = DataSink('FsmScan')
        self.sink.__enter__()


    def teardown(self):
        self.sink.__exit__()


    def update(self):
        self.sink.pop() #wait for some data to be saved to sink
        self.set_data(self.sink.datasets['xSteps'], self.sink.datasets['ySteps'], self.sink.datasets['ScanCounts'])
