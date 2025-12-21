
from functools import partial
from importlib import reload

from pyqtgraph.Qt import QtWidgets
from pyqtgraph import SpinBox
from nspyre import LinePlotWidget
from nspyre import ParamsWidget
from nspyre import ProcessRunner
from nspyre import DataSink

import experiments.DEER as DEER
import numpy as np


class DEER_Widget(QtWidgets.QWidget):
    
    def __init__(self):
        super().__init__()

        self.setWindowTitle('DEER')

        self.params_widget = ParamsWidget({

            'datasetName': {
                'display_text': 'Dataset Name',
                'widget': QtWidgets.QLineEdit('DEERData'),
            },

            'num_samples': {
                'display_text': 'Samples per tau',
                'widget': SpinBox(
                    value=10000,
                    dec=True,
                    int=True
                ),
            },

            'maxIterations': {
                'display_text': 'Iterations',
                'widget': SpinBox(
                    value=20,
                    dec=True,
                    int=True,
                ),
            },

            'laser_power': {
                'display_text': 'Laser Power',
                'widget': SpinBox(
                    value=2,
                    bounds=(0,5),
                    dec=True,
                ),
            },

            'NV_freq': {
                'display_text': 'NV Frequency',
                'widget': SpinBox(
                    value=2.87e9,
                    suffix='Hz',
                    siPrefix=True,
                    dec=True,
                ),
            },

            'NV_rf_power': {
                'display_text': 'NV RF Power',
                'widget': SpinBox(
                    value=-6,
                    suffix='dBm',
                    siPrefix=False,
                    bounds=(-42, 6),
                    dec=True,
                ),
            },

            'NV_pi_time': {
                'display_text': 'NV Pi Pulse',
                'widget': SpinBox(
                    value=150e-9, #300ns ish
                    suffix='s',
                    siPrefix=True, 
                    bounds=(0, 1),
                    dec=True,
                ),
            },

            'P1_freq': {
                'display_text': 'P1 Frequency',
                'widget': SpinBox(
                    value=2.87e9,
                    suffix='Hz',
                    siPrefix=True,
                    dec=True,
                ),
            },

            'P1_rf_power': {
                'display_text': 'P1 RF Power',
                'widget': SpinBox(
                    value=0,
                    suffix='dBm',
                    siPrefix=False,
                    bounds=(-42, 6),
                    dec=True,
                ),
            },

            'P1_pi_time': {
                'display_text': 'P1 Pi Pulse',
                'widget': SpinBox(
                    value=150e-9, #300ns ish
                    suffix='s',
                    siPrefix=True, 
                    bounds=(0, 1),
                    dec=True,
                ),
            },

            'tau_min': {
                'display_text': 'Min Tau delay',
                'widget': SpinBox(
                    value=10e-9,
                    suffix='s',
                    siPrefix=True,
                    #bounds=(100e3, 10e9),
                    dec=True,
                ),
            },

            'tau_max': {
                'display_text': 'Max Tau delay',
                'widget': SpinBox(
                    value=1e-6,
                    suffix='s',
                    siPrefix=True,
                    #bounds=(100e3, 10e9),
                    dec=True,
                ),
            },

            'tau_num': {
                'display_text': 'Number of Taus',
                'widget': SpinBox(
                    value=21,
                    #suffix='us',
                    siPrefix=False,
                    #bounds=(100e3, 10e9),
                    dec=True,
                    int=True
                ),
            },

        })

        #Set-up check boxes for each channel
        
        #Set-up check boxes
        self.turnLaserOffAtEndButton = QtWidgets.QCheckBox('Turn Laser Off?')
        
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
        reload(DEER)

        # create an instance of the ODMR class that implements the experimental logic.
        rabiMeas = DEER.DEER_Measurement()

        # run the sweep function in a new thread
        self.sweepProc.run(
            rabiMeas.DEER,
            self.params_widget.datasetName,
            self.params_widget.num_samples,
            self.params_widget.maxIterations,
            self.params_widget.laser_power,
            self.params_widget.NV_freq,
            self.params_widget.NV_rf_power,
            int(1e9*self.self.params_widget.NV_pi_time),
            self.params_widget.P1_freq,
            self.params_widget.P1_rf_power,
            int(1e9*self.self.params_widget.P1_pi_time),
            int(1e9*self.params_widget.tau_min),
            int(1e9*self.params_widget.tau_max),
            self.params_widget.tau_num,    
        )


    def stop(self):
        """Stop the sweep process."""
        
        #kill the TaskVsTime sweep process
        self.sweepProc.kill()



class DEERPlotWidget(LinePlotWidget):

    def __init__(self):
        title = 'DEER'
        super().__init__(title=title, xlabel='Times (ns)', ylabel='PL')


    def setup(self):
        self.sink = DataSink('DEERData')
        self.sink.__enter__() 
        self.add_plot('MW_ON')
        self.add_plot('MW_OFF')
        self.plot_widget.setYRange(-100, 5100)


    def teardown(self):
        self.sink.stop()
        self.sink.__exit__()

    def update(self):
        self.sink.pop() #wait for some data to be saved to sink
        # update the plot
        tau_list = np.sort(self.sink.datasets['tau_times'])
        self.set_data('MW_ON', tau_list, [np.mean(self.sink.datasets['mwCountsDict'][tau_time]) for tau_time in tau_list])
        self.set_data('MW_OFF', tau_list, [np.mean(self.sink.datasets['noMwCountsDict'][tau_time]) for tau_time in tau_list])

