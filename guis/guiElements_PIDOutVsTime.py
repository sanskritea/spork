'''
NSpyre v0.6.1 GUI for AUX Out VS Time application for NV AFM

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

import experiments.PIDOutVsTime as PIDOutVsTime
from guis.guiElements_general import AutoSaveWidget
from guis.guiElements_general import flexSave

class PIDOutVsTime_Widget(QtWidgets.QWidget):

	def __init__(self):
		super().__init__()

		self.setWindowTitle('PIDOutVsTime')

		self.params_widget = ParamsWidget({
            'datasetname': {
                'display_text': 'Dataset Name',
                'widget': QtWidgets.QLineEdit('PIDOutVsTime'),
            },

            'time_per_point': {
                'display_text': 'Time per point',
                'widget': SpinBox(
                    value=0.1,
                    dec=True,
                ),
            },

            'aux_chan': {
                'display_chan': 'AUX Output Channel',
                'widget': SpinBox(
                    value=0,
                    dec=True,
                    int=True,
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
	    reload(PIDOutVsTime)

	    # create an instance of the ODMR class that implements the experimental logic.
	    PIDOutVsTimeMeas = PIDOutVsTime.PIDOutVsTime_Measurement()

	    # run the sweep function in a new thread
	    self.sweepProc.run(
	        PIDOutVsTimeMeas.auxoutvstime,
	        self.params_widget.datasetname, 
	        self.params_widget.time_per_point, 
	        self.params_widget.aux_chan,
	    )


	def stop(self):
	    """Stop the sweep process."""

	    # kill the sweep process
	    self.sweepProc.kill()


class PIDOutVsTimePlotWidget(LinePlotWidget):

    def __init__(self):
        title = 'PID Out Vs Time'
        super().__init__(title=title, xlabel='time (s)', ylabel='AUXOUT value (V)') #TODO: Switch this over to contrast


    def setup(self):
        self.sink = DataSink('PIDOutVsTime')
        self.sink.__enter__() 
        self.add_plot('PID output')


    def teardown(self):
        self.sink.stop()
        self.sink.__exit__()


    def update(self):
        self.sink.pop() # wait for some data to be saved to sink
        # update the plot
        self.set_data('PID output', self.sink.datasets['time'], self.sink.datasets['AUX'])



