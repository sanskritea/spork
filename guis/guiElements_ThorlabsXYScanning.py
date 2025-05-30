"""
GUI for a XY Scanning (with Attocube scanner and DAQ analog control) application

Sanskriti Chitransh, 2023-Oct-26

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
from nspyre.gui.widgets.heatmap import HeatMapWidget

import experiments.ThorlabsXYScanning as ThorlabsXYScanning
import numpy as np


class ThorlabsXYScanning_Widget(QtWidgets.QWidget):
    
    def __init__(self):
        super().__init__()

        self.setWindowTitle('XY Scanning')

        self.params_widget = ParamsWidget({
            'datasetName': {
                'display_text': 'Dataset Name',
                'widget': QtWidgets.QLineEdit('ScanningData'),
            },

            'device': {
                'display_text': 'DAQ Device Name',
                'widget': QtWidgets.QLineEdit('Dev4'),
            },

            'x_init_position': {
                'display_text': 'X Initial Position mm',
                'widget': SpinBox(
                    value=10.54,
                    dec=True,
                ),
            },

            'y_init_position': {
                'display_text': 'Y Initial Position mm',
                'widget': SpinBox(
                    value=9.34 ,
                    dec=True,
                ),
            },

            'position_steps': {
                'display_text': 'Number of Steps',
                'widget': SpinBox(
                    value=20,
                    dec=True,
                    int=True,
                ),
            },

            'step_length': {
                'display_text': 'Step Length um',
                'widget': SpinBox(
                    value=5,
                    dec=True,
                ),
            },

            'time_per_pixel': {
                'display_text': 'Time per pixel',
                'widget': SpinBox(
                    value=0.005,
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
        reload(ThorlabsXYScanning)

        # create an instance of the ODMR class that implements the experimental logic.
        ThorlabsXYScanningMeas = ThorlabsXYScanning.XYScan()

        # self.turnLaserOffAtEndButton.setEnabled(False)

        # run the sweep function in a new thread
        self.sweepProc.run(
            ThorlabsXYScanningMeas.scanning,
            self.params_widget.datasetName, 
            self.params_widget.device,
            self.params_widget.x_init_position,
            self.params_widget.y_init_position, 
            self.params_widget.position_steps,
            self.params_widget.step_length, 
            self.params_widget.time_per_pixel, 
        )


    def stop(self):
        """Stop the sweep process."""
        #Re-enable checkbox
        # self.turnLaserOffAtEndButton.setEnabled(True)

        #kill the TaskVsTime sweep process
        self.sweepProc.kill()



class ThorlabsXYScanningPlotWidget(HeatMapWidget):

    def __init__(self):
        title = 'XY Scan'
        super().__init__(title=title, btm_label='X position ', lft_label='Y position') #TODO: Switch this over to contrast


    def setup(self):
        self.sink = DataSink('ScanningData')
        self.sink.__enter__() 


    def teardown(self):
        self.sink.stop()
        self.sink.__exit__()


    def update(self):
        self.sink.pop() # wait for some data to be saved to sink
        # update the plot
        self.set_data(self.sink.datasets['x_position'], self.sink.datasets['y_position'], self.sink.datasets['counts'])
        self.plot_item.enableAutoRange(True)


