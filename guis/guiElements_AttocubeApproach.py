'''
NSpyre v0.6.1 GUI for Attocube Approach application for NV AFM

Sanskriti Chitransh, 2023-Oct-24

'''

#TO-DO, NEEDS SOME DEBUGGING
from functools import partial
from importlib import reload

from pyqtgraph.Qt import QtWidgets
from pyqtgraph import SpinBox
from nspyre import LinePlotWidget
from nspyre import ParamsWidget
from nspyre import ProcessRunner
from nspyre import DataSink

import experiments.AttocubeApproach as AttocubeApproach
from guis.guiElements_general import AutoSaveWidget
from guis.guiElements_general import flexSave
import numpy as np


class Attocube_Approach_Widget(QtWidgets.QWidget):

	def __init__(self):
		super().__init__()

		self.setWindowTitle('Attocube Approach')

		self.params_widget = ParamsWidget({
            'datasetName': {
                'display_text': 'Dataset Name',
                'widget': QtWidgets.QLineEdit('AttocubeApproach'),
            },

            'device': { # use the 6343 (Dev4)
                'display_text': 'DAQ Device Name',
                'widget': QtWidgets.QLineEdit('Dev4'),
            },

            'AOchannel': {
                'display_text': 'DAQ Analog Output Channel Number',
                'widget': QtWidgets.QLineEdit('AO1'),
            },

            'step_wait': {
                'display_text': 'Wait time between Stepping and PID read',
                'widget': SpinBox(
                    value=0.1,
                    bounds=(0.1, 15),
                    dec=True,
                ),
            },

            'stage_min': {
                'display_text': 'Min DAQ AO Voltage to Scanner Z',
                'widget': SpinBox(
                    value=-0.0069,
                    dec=True,
                ),
            },

            'stage_max': {
                'display_text': 'Max DAQ AO Voltage to Scanner Z',
                'widget': SpinBox(
                    value=0,
                    bounds=(-0.0069, 3.9515),
                    dec=True,
                ),
            },

            'stage_voltage': {
                'display_text': 'Number of Voltage Steps',
                'widget': SpinBox(
                    value=1001,
                    dec=True,
                    int=True
                ),
            },

            'threshold': {
            	'display_text': 'PID Setpoint Threshold',
            	'widget': SpinBox(
            		value=0.995,
            		bounds=(0.5, 0.999),
            		dec=True,
            		)
            },

            'A_init': {
            	'display_text': 'Initial AO Voltage',
            	'widget': SpinBox(
            		value=0.001,
            		dec=True,
            		)
            }
        })

        # Setup run and stop buttons
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
	    reload(AttocubeApproach)

	    # create an instance of the ODMR class that implements the experimental logic.
	    AttocubeApproachMeas = AttocubeApproach.Attocube_Approach_Measurement()

	    # run the sweep function in a new thread
	    self.sweepProc.run(
	        AttocubeApproachMeas.attocubeapproach,
	        self.params_widget.datasetName, 
	        self.params_widget.device, 
	        self.params_widget.AOchannel,
	        self.params_widget.step_wait, # TIME VALUE CHECK AGAIN
	        self.params_widget.stage_min,
	        self.params_widget.stage_max,
	        self.params_widget.stage_voltage,
	        self.params_widget.threshold,
	        self.params_widget.A_init
	    )


	def stop(self):
	    """Stop the sweep process."""

	    # kill the sweep process
	    self.sweepProc.kill()


class AttocubeApproachPlotWidget(LinePlotWidget):

    def __init__(self):
        title = 'Attocube Approach'
        super().__init__(title=title, xlabel='Scanner Voltage (V)', ylabel='Fork Amplitude (V)') #TODO: Switch this over to contrast


    def setup(self):
        self.sink = DataSink('AttocubeApproach')
        self.sink.__enter__() 
        self.add_plot('Fork Amplitude')


    def teardown(self):
        self.sink.stop()
        self.sink.__exit__()


    def update(self):
        self.sink.pop() # wait for some data to be saved to sink
        # update the plot
        self.set_data('Fork Amplitude', self.sink.datasets['voltage'], self.sink.datasets['amplitude'])



