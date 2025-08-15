"""
GUI for a Newport Linescan application

Sanskriti Chitransh, 2024-June-5

Adapted from Chris Egerstrom

"""

#TO-DO, NEEDS SOME DEBUGGING
from functools import partial
from importlib import reload

from pyqtgraph.Qt import QtWidgets
from pyqtgraph import SpinBox
from nspyre import ParamsWidget
from nspyre import ProcessRunner
from nspyre import DataSink

# from guis.ian_heatmap_crosshair import HeatMapWidget
from nspyre import LinePlotWidget

import experiments.NewportLineScanning as NewportLineScanning
import numpy as np


class NewportLineScanning_Widget(QtWidgets.QWidget):
    
    def __init__(self):
        super().__init__()

        self.setWindowTitle('Line Scanning')

        self.params_widget = ParamsWidget({
            'datasetName': {
                'display_text': 'Dataset Name',
                'widget': QtWidgets.QLineEdit('NewportLineScanningData'),
            },

            'x_init_position': {
                'display_text': 'X Init Pos mm',
                'widget': SpinBox(
                    value=0,
                    dec=True,
                ),
            },

            'y_init_position': {
                'display_text': 'Y Init Pos mm',
                'widget': SpinBox(
                    value=0 ,
                    dec=True,
                ),
            },

            'z_init_position': {
                'display_text': 'Z Init Pos mm',
                'widget': SpinBox(
                    value=0 ,
                    dec=True,
                ),
            },

            'scanning_axis': {
                'display_text': 'Linescan axis',
                'widget': SpinBox(
                    value=1 ,
                    dec=True,
                    int=True
                ),
            },

            'position_steps': {
                'display_text': 'Number of Steps',
                'widget': SpinBox(
                    value=50,
                    dec=True,
                    int=True,
                ),
            },

            'step_length': {
                'display_text': 'Step Length mm',
                'widget': SpinBox(
                    value=0.0002,
                    dec=True,
                ),
            },

            'time_per_pixel': {
                'display_text': 'Time per pixel ',
                'widget': SpinBox(
                    value=0.005,
                    dec=True,
                ),
            },

            'step_wait': {
                'display_text': 'Wait after move',
                'widget': SpinBox(
                    value=0.03,
                    dec=True,
                ),
            },

        })

        # #Set-up check boxes
        # self.turnLaserOffAtEndButton = QtWidgets.QCheckBox('Turn Laser Off?')
        
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
        # params_layout.addWidget(self.turnLaserOffAtEndButton)
        params_layout.addStretch()
        params_layout.addWidget(runButton)
        params_layout.addWidget(stopButton)

        self.setLayout(params_layout)


    def runClicked(self):
        """Runs when the 'run' button is pressed."""

        # reload the spin measurements module at runtime in case any changes were made to the code
        reload(NewportLineScanning)

        # create an instance of the ODMR class that implements the experimental logic.
        NewportLineScanningMeas = NewportLineScanning.LineScan()

        # self.turnLaserOffAtEndButton.setEnabled(False)

        # run the sweep function in a new thread
        self.sweepProc.run(
            NewportLineScanningMeas.scanning,
            self.params_widget.datasetName,
            self.params_widget.x_init_position,
            self.params_widget.y_init_position,
            self.params_widget.z_init_position,
            self.params_widget.scanning_axis, 
            self.params_widget.position_steps,
            self.params_widget.step_length, 
            self.params_widget.time_per_pixel, 
            self.params_widget.step_wait,
        )


    def stop(self):
        """Stop the sweep process."""
        #Re-enable checkbox
        # self.turnLaserOffAtEndButton.setEnabled(True)

        #kill the TaskVsTime sweep process
        self.sweepProc.kill()



class NewportLineScanningPlotWidget(LinePlotWidget):

    def __init__(self):
        title = 'Line Scan'
        super().__init__(title=title, xlabel='Scan position (mm)', ylabel='Counts')


    def setup(self):
        self.sink = DataSink('NewportLineScanningData')
        self.sink.__enter__() 
        self.add_plot('PL')


    def teardown(self):
        self.sink.stop()
        self.sink.__exit__()


    def update(self):
        self.sink.pop() # wait for some data to be saved to sink
        # update the plot
        self.set_data('PL', self.sink.datasets['position'], self.sink.datasets['counts'])
        # self.plot_item.enableAutoRange(True)


