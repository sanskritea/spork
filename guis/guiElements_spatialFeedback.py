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
from nspyre import LinePlotWidget

import experiments.spatialFeedback as spatialFeedback
from guis.guiElements_general import flexSave


class CustomSpatialFeedbackWidget(QtWidgets.QWidget):
    
    def __init__(self):
        super().__init__()

        self.setWindowTitle('Spatial Feedback')

        self.params_widget = ParamsWidget({
            'xGuess': {
                'display_text': 'X Guess',
                'widget': SpinBox(
                    value=1e-6,
                    suffix='m',
                    siPrefix=True,
                    bounds=(-50e-6, 50e-6),
                    dec=True,
                ),
            },
            'yGuess': {
                'display_text': 'Y Guess',
                'widget': SpinBox(
                    value=1e-6,
                    suffix='m',
                    siPrefix=True,
                    bounds=(-50e-6, 50e-6),
                    dec=True,
                ),
            },
            'xRange': {
                'display_text': 'X Range (one-sided)',
                'widget': SpinBox(
                    value=2e-6,
                    suffix='m',
                    siPrefix=True,
                    bounds=(-50e-6, 50e-6),
                    dec=True,
                ),
            },
            'yRange': {
                'display_text': 'Y Range (one-sided)',
                'widget': SpinBox(
                    value=2e-6,
                    suffix='m',
                    siPrefix=True,
                    bounds=(-50e-6, 50e-6),
                    dec=True,
                ),
            },
            'xPts': {
                'display_text': 'Number of X Pts',
                'widget': SpinBox(
                    value=30,
                    int=True,
                    dec=True,
                ),
            },
            'yPts': {
                'display_text': 'Number of Y Pts',
                'widget': SpinBox(
                    value=30,
                    int=True,
                    dec=True,
                ),
            },
            'collectsPerPt': {
                'display_text': 'Points per Position',
                'widget': SpinBox(
                    value=1000,
                    int=True,
                    dec=True,
                ),
            },
            'iters': {
                'display_text': 'Num of Iterations',
                'widget': SpinBox(
                    value=2,
                    int=True,
                    dec=True,
                ),
            },
            'datasetName': {
                'display_text': 'Dataset Name',
                'widget': QtWidgets.QLineEdit('SpatialFeedback'),
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
        reload(spatialFeedback)

        # create an instance of the ODMR class that implements the experimental logic.
        spatialFeedbackMeas = spatialFeedback.SpatialFeedbackMeasurement()

        # run the sweep function in a new thread
        self.sweepProc.run(
            spatialFeedbackMeas.spatialFeedback,
            self.params_widget.datasetName,
            1e6*self.params_widget.xGuess, #convert m=>um since FSM default unit is um
            1e6*self.params_widget.yGuess,
            2e6*self.params_widget.xRange, #1-sided to 2-sided ranges
            2e6*self.params_widget.yRange,
            self.params_widget.xPts,
            self.params_widget.yPts,
            self.params_widget.collectsPerPt,
            self.params_widget.iters
        )


    def stop(self):
        """Stop the sweep process."""
       #kill the TaskVsTime sweep process
        self.sweepProc.kill()
        flexSave(self.params_widget.datasetName, 'spatialFeedback', 'closeout') #if process is closed



class XSpatialFeedbackPlotWidget(LinePlotWidget):
    def __init__(self):
        title = 'Latest X Spatial Feedback'
        super().__init__(title=title, xlabel='X', ylabel='Counts')


    def setup(self):
        self.sink = DataSink('SpatialFeedback')
        self.sink.__enter__()
        self.add_plot('Latest X-Scan')
        self.add_plot('Latest X-Fit')


    def teardown(self):
        self.sink.__exit__()


    def update(self):
        self.sink.pop() #wait for some data to be saved to sink
        try: #because -1 indexing occasionally breaks this plotter and I'm lazy #TODO: Implement better error handling
            self.set_data('Latest X-Scan', self.sink.datasets['xSteps'][-1,:], self.sink.datasets['xSweepData'][-1,:])
            fitFxn = lambda x: spatialFeedback.customGaussian(x, self.sink.datasets['xFits'][-1][0], self.sink.datasets['xFits'][-1][1], self.sink.datasets['xFits'][-1][2], self.sink.datasets['xFits'][-1][3])
            self.set_data('Latest X-Fit', self.sink.datasets['xSteps'][-1,:], fitFxn(self.sink.datasets['xSteps'][-1,:]) )
        except:
            pass



class YSpatialFeedbackPlotWidget(LinePlotWidget):
    def __init__(self):
        title = 'Latest Y Spatial Feedback'
        super().__init__(title=title, xlabel='Y', ylabel='Counts')


    def setup(self):
        self.sink = DataSink('SpatialFeedback')
        self.sink.__enter__()
        self.add_plot('Latest Y-Scan')
        self.add_plot('Latest Y-Fit')


    def teardown(self):
        self.sink.__exit__()


    def update(self):
        self.sink.pop() #wait for some data to be saved to sink
        try: #because -1 indexing occasionally breaks this plotter and I'm lazy #TODO: Implement better error handling
            self.set_data('Latest Y-Scan', self.sink.datasets['ySteps'][-1,:], self.sink.datasets['ySweepData'][-1,:])
            fitFxn = lambda x: spatialFeedback.customGaussian(x, self.sink.datasets['yFits'][-1][0], self.sink.datasets['yFits'][-1][1], self.sink.datasets['yFits'][-1][2], self.sink.datasets['yFits'][-1][3])
            self.set_data('Latest Y-Fit', self.sink.datasets['ySteps'][-1,:], fitFxn(self.sink.datasets['ySteps'][-1,:]) )
        except:
            pass
