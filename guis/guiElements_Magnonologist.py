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

import experiments.Magnonologist as Magnonologist
import numpy as np


class Magnonologist_Widget(QtWidgets.QWidget):
    
    def __init__(self):
        super().__init__()

        self.setWindowTitle('T1 Scanning')

        self.params_widget = ParamsWidget({
            'datasetName': {
                'display_text': 'Dataset Name',
                'widget': QtWidgets.QLineEdit('MagnonologistData'),
            },

            'x_init_voltage': {
                'display_text': 'Initial DAQ AO2',
                'widget': SpinBox(
                    value=-0.1,
                    dec=True,
                ),
            },

            'x_final_voltage': {
                'display_text': 'Final DAQ AO2',
                'widget': SpinBox(
                    value=8,
                    dec=True,
                ),
            },

            'y_init_voltage': {
                'display_text': 'Initial DAQ AO3',
                'widget': SpinBox(
                    value=-0.1,
                    dec=True,
                ),
            },

            'y_final_voltage': {
                'display_text': 'Final DAQ AO3',
                'widget': SpinBox(
                    value=8,
                    dec=True,
                ),
            },

            'x_voltage_steps': {
                'display_text': 'Number of AO2 Steps',
                'widget': SpinBox(
                    value=1001,
                    dec=True,
                    int=True,
                ),
            },

            'y_voltage_steps': {
                'display_text': 'Number of AO3 Steps',
                'widget': SpinBox(
                    value=1001,
                    dec=True,
                    int=True,
                ),
            },

            'rf_power': {
                'display_text': 'RF Power',
                'widget': SpinBox(
                    value=-6,
                    suffix='dBm',
                    siPrefix=False,
                    bounds=(-42, 6),
                    dec=True,
                ),
            },

            'laser_power': {
                'display_text': 'Laser power',
                'widget': SpinBox(
                    value=2,
                    bounds=(0,5),
                    dec=True,
                ),
            },

            'NV_freq': {
                'display_text': 'NV Freq',
                'widget': SpinBox(
                    value=2.82e9,
                    suffix='Hz',
                    siPrefix=True, 
                    bounds=(2e9, 4e9),
                    dec=True,
                ),
            },

            'pi_time': {
                'display_text': 'Pi pulse',
                'widget': SpinBox(
                    value=100e-9,
                    suffix='s',
                    siPrefix=True, 
                    dec=True,
                ),
            },

            'Bayesian_R': {
                'display_text': 'Samples per tau',
                'widget': SpinBox(
                    value=100,
                    dec=True,
                    int=True
                ),
            },

            'objective_x_init_position': {
                'display_text': 'Objective X Pos mm',
                'widget': SpinBox(
                    value=0,
                    dec=True,
                ),
            },

            'objective_y_init_position': {
                'display_text': 'Objective Y Pos mm',
                'widget': SpinBox(
                    value=0 ,
                    dec=True,
                ),
            },

            'objective_z_init_position': {
                'display_text': 'Objective Z Pos mm',
                'widget': SpinBox(
                    value=0 ,
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
        reload(Magnonologist)

        # create an instance of the ODMR class that implements the experimental logic.
        MagnonologistMeas = Magnonologist.Magnonologist_Measurement()

        # self.turnLaserOffAtEndButton.setEnabled(False)

        # run the sweep function in a new thread
        self.sweepProc.run(
            MagnonologistMeas.MagnonologistScan,
            self.params_widget.datasetName, 
            self.params_widget.x_init_voltage, 
            self.params_widget.x_final_voltage,
            self.params_widget.y_init_voltage, 
            self.params_widget.y_final_voltage,
            self.params_widget.x_voltage_steps,
            self.params_widget.y_voltage_steps, 
            self.params_widget.rf_power,
            self.params_widget.laser_power,
            self.params_widget.NV_freq,
            int(1e9*self.params_widget.pi_time),
            self.params_widget.Bayesian_R, 
            self.params_widget.objective_x_init_position,
            self.params_widget.objective_y_init_position,
            self.params_widget.objective_z_init_position,
        )


    def stop(self):
        """Stop the sweep process."""
        #Re-enable checkbox
        # self.turnLaserOffAtEndButton.setEnabled(True)

        #kill the TaskVsTime sweep process
        self.sweepProc.kill()



class MagnonologistPlotWidget(HeatMapWidget):

    def __init__(self):
        title = 'T1 Scan'
        super().__init__(title=title, btm_label='Y voltage/ X objective direction', lft_label='X voltage/ Y objective direction') 


    def setup(self):
        self.sink = DataSink('MagnonologistData')
        self.sink.__enter__() 


    def teardown(self):
        self.sink.stop()
        self.sink.__exit__()


    def update(self):
        self.sink.pop() # wait for some data to be saved to sink
        # update the plot
        # self.plot_item.invertX()
        # self.plot_item.invertY()
        self.set_data(self.sink.datasets['yvoltage'], self.sink.datasets['xvoltage'], self.sink.datasets['t1'])


