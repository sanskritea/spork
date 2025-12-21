
from functools import partial
from importlib import reload

from pyqtgraph.Qt import QtWidgets
from pyqtgraph import SpinBox
from nspyre import LinePlotWidget
from nspyre import ParamsWidget
from nspyre import ProcessRunner
from nspyre import DataSink

import experiments.P1Rabi as P1Rabi
import numpy as np


class P1Rabi_Widget(QtWidgets.QWidget):
    
    def __init__(self):
        super().__init__()

        self.setWindowTitle('P1Rabi')

        self.params_widget = ParamsWidget({

            'datasetName': {
                'display_text': 'Dataset Name',
                'widget': QtWidgets.QLineEdit('P1RabiData'),
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

            'P1_min_MW_time': {
                'display_text': 'Min P1 MW duration',
                'widget': SpinBox(
                    value=10e-9,
                    suffix='s',
                    siPrefix=True,
                    #bounds=(100e3, 10e9),
                    dec=True,
                ),
            },

            'P1_max_MW_time': {
                'display_text': 'Max P1 MW duration',
                'widget': SpinBox(
                    value=1e-6,
                    suffix='s',
                    siPrefix=True,
                    #bounds=(100e3, 10e9),
                    dec=True,
                ),
            },

            'P1_num_MW_times': {
                'display_text': 'Number of P1 Taus',
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
        reload(P1Rabi)

        # create an instance of the ODMR class that implements the experimental logic.
        rabiMeas = P1Rabi.P1Rabi_Measurement()

        # run the sweep function in a new thread
        self.sweepProc.run(
            rabiMeas.P1Rabi,
            self.params_widget.datasetName,
            self.params_widget.num_samples,
            self.params_widget.maxIterations,
            self.params_widget.laser_power,
            self.params_widget.NV_freq,
            self.params_widget.NV_rf_power,
            int(1e9*self.self.params_widget.NV_pi_time),
            self.params_widget.P1_freq,
            self.params_widget.P1_rf_power,
            int(1e9*self.params_widget.P1_min_MW_time),
            int(1e9*self.params_widget.P1_max_MW_time),
            self.params_widget.P1_num_MW_times,    
        )


    def stop(self):
        """Stop the sweep process."""
        
        #kill the TaskVsTime sweep process
        self.sweepProc.kill()



class P1RabiPlotWidget(LinePlotWidget):

    def __init__(self):
        title = 'P1Rabi'
        super().__init__(title=title, xlabel='Times (ns)', ylabel='PL')


    def setup(self):
        self.sink = DataSink('P1RabiData')
        self.sink.__enter__() 
        self.add_plot('NV_MW_ON')
        self.add_plot('NV_MW_OFF')
        self.plot_widget.setYRange(-100, 5100)


    def teardown(self):
        self.sink.stop()
        self.sink.__exit__()

    def update(self):
        self.sink.pop() #wait for some data to be saved to sink
        # update the plot
        P1_mw_time_list = np.sort(self.sink.datasets['P1_mw_times'])
        self.set_data('NV_MW_ON', P1_mw_time_list, [np.mean(self.sink.datasets['mwCountsDict'][P1_mw_time]) for P1_mw_time in P1_mw_time_list])
        self.set_data('NV_MW_OFF', P1_mw_time_list, [np.mean(self.sink.datasets['noMwCountsDict'][P1_mw_time]) for P1_mw_time in P1_mw_time_list])

