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

# from nspyre.gui.widgets.heatmap import HeatMapWidget
from guis.ian_heatmap_crosshair import HeatMapWidget

import experiments.AttocubeXYScanning as AttocubeXYScanning
import numpy as np


class AttocubeXYScanning_Widget(QtWidgets.QWidget):
    
    def __init__(self):
        super().__init__()

        self.setWindowTitle('Attocube XY Scanning')

        self.params_widget = ParamsWidget({
            'datasetName': {
                'display_text': 'Dataset Name',
                'widget': QtWidgets.QLineEdit('AttocubeXYScanningData'),
            },

            'device': {
                'display_text': 'DAQ Device Name',
                'widget': QtWidgets.QLineEdit('Dev4'),
            },

            'x_init_voltage': {
                'display_text': 'DAQ Initial Voltage X axis',
                'widget': SpinBox(
                    value=-0.1,
                    dec=True,
                ),
            },

            'x_final_voltage': {
                'display_text': 'DAQ Final Voltage X axis',
                'widget': SpinBox(
                    value=8,
                    dec=True,
                ),
            },

            'y_init_voltage': {
                'display_text': 'DAQ Initial Voltage Y axis',
                'widget': SpinBox(
                    value=-0.1,
                    dec=True,
                ),
            },

            'y_final_voltage': {
                'display_text': 'DAQ Final Voltage Y axis',
                'widget': SpinBox(
                    value=8,
                    dec=True,
                ),
            },

            'x_voltage_steps': {
                'display_text': 'Number of X Steps',
                'widget': SpinBox(
                    value=1001,
                    dec=True,
                    int=True,
                ),
            },

            'y_voltage_steps': {
                'display_text': 'Number of Y Steps',
                'widget': SpinBox(
                    value=1001,
                    dec=True,
                    int=True,
                ),
            },

            'time_per_pixel': {
                'display_text': 'Time per pixel',
                'widget': SpinBox(
                    value=0.01,
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
        reload(AttocubeXYScanning)

        # create an instance of the ODMR class that implements the experimental logic.
        AttocubeXYScanningMeas = AttocubeXYScanning.XYScan()

        # self.turnLaserOffAtEndButton.setEnabled(False)

        # run the sweep function in a new thread
        self.sweepProc.run(
            AttocubeXYScanningMeas.scanning,
            self.params_widget.datasetName, 
            self.params_widget.device,
            self.params_widget.x_init_voltage, 
            self.params_widget.x_final_voltage,
            self.params_widget.y_init_voltage, 
            self.params_widget.y_final_voltage,
            self.params_widget.x_voltage_steps,
            self.params_widget.y_voltage_steps, 
            self.params_widget.time_per_pixel, 
        )


    def stop(self):
        """Stop the sweep process."""
        #Re-enable checkbox
        # self.turnLaserOffAtEndButton.setEnabled(True)

        #kill the TaskVsTime sweep process
        self.sweepProc.kill()



class AttocubeXYScanningPlotWidget(HeatMapWidget):

    def __init__(self):
        title = 'XY Scan'
        super().__init__(title=title, btm_label='Y voltage/ X objective direction', lft_label='X voltage/ Y objective direction') 


    def setup(self):
        self.sink = DataSink('AttocubeXYScanningData')
        self.sink.__enter__() 


    def teardown(self):
        self.sink.stop()
        self.sink.__exit__()


    def update(self):
        self.sink.pop() # wait for some data to be saved to sink
        # update the plot
        # self.plot_item.invertX()
        # self.plot_item.invertY()
        self.set_data(self.sink.datasets['yvoltage'], self.sink.datasets['xvoltage'], self.sink.datasets['counts'])


