
from functools import partial
from importlib import reload

from pyqtgraph.Qt import QtWidgets
from pyqtgraph import SpinBox
from nspyre import LinePlotWidget
from nspyre import ParamsWidget
from nspyre import ProcessRunner
from nspyre import DataSink

import experiments.Rabi as Rabi
import numpy as np


class Rabi_Widget(QtWidgets.QWidget):
    
    def __init__(self):
        super().__init__()

        self.setWindowTitle('Rabi')

        self.params_widget = ParamsWidget({

            'datasetName': {
                'display_text': 'Dataset Name',
                'widget': QtWidgets.QLineEdit('RabiData'),
            },

            'samplingFreq': {
                'display_text': 'Sampling Freq',
                'widget': SpinBox(
                    value=20e6,
                    suffix='Hz',
                    siPrefix=True,
                    dec=True,
                ),
            },

            'maxIterations': {
                'display_text': 'Iterations',
                'widget': SpinBox(
                    value=100,
                    dec=True,
                    int=True,
                ),
            },

            'freq': {
                'display_text': 'Frequency',
                'widget': SpinBox(
                    value=2.87e9,
                    suffix='Hz',
                    siPrefix=True,
                    dec=True,
                ),
            },

            'min_MW_time': {
                'display_text': 'Min MW pulse duration',
                'widget': SpinBox(
                    value=10e-9,
                    suffix='s',
                    siPrefix=True,
                    #bounds=(100e3, 10e9),
                    dec=True,
                ),
            },

            'max_MW_time': {
                'display_text': 'Max MW pulse duration',
                'widget': SpinBox(
                    value=1e-6,
                    suffix='s',
                    siPrefix=True,
                    #bounds=(100e3, 10e9),
                    dec=True,
                ),
            },

            'num_MW_times': {
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

            'rf_power': {
                'display_text': 'RF Power',
                'widget': SpinBox(
                    value=6,
                    suffix='dBm',
                    siPrefix=False,
                    bounds=(-42, 6),
                    dec=True,
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

            'num_samples': {
                'display_text': 'Samples per MW tau',
                'widget': SpinBox(
                    value=10000,
                    dec=True,
                    int=True
                ),
            },

            'clock_time': {
                'display_text': 'Clock pulse duration',
                'widget': SpinBox(
                    value=11e-9,
                    suffix='s',
                    siPrefix=True,
                    bounds=(1e-9, 1),
                    dec=True,
                ),
            },

            'init_time': {
                'display_text': 'NV Init Pulse',
                'widget': SpinBox(
                    value=2e-6,
                    suffix='s',
                    siPrefix=True, #kiloseconds are cursed
                    #bounds=(100e3, 10e9),
                    dec=True,
                ),
            },

            'laser_lag': {
                'display_text': 'Laser Stabilization Lag',
                'widget': SpinBox(
                    value=300e-9,
                    suffix='s',
                    siPrefix=True, #kiloseconds are cursed
                    #bounds=(100e3, 10e9),
                    dec=True,
                ),
            },

            'probe_time': {
                'display_text': 'Readout Duration',
                'widget': SpinBox(
                    value=300e-9,
                    suffix='s',
                    siPrefix=True, 
                    #bounds=(100e3, 10e9),
                    dec=True,
                ),
            },

            'singlet_decay': {
                'display_text': 'NV Singlet Decay',
                'widget': SpinBox(
                    value=1e-6,
                    suffix='s',
                    siPrefix=True, 
                    #bounds=(100e3, 10e9),
                    dec=True,
                ),
            },

            # 'initial_counts': {
            #     'display_text': 'Starting PL from CountVSTime',
            #     'widget': SpinBox(
            #         value=270000,
            #         dec=True,
            #     ),
            # },

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

            # 'threshold': {
            #     'display_text': 'Feedback threshold',
            #     'widget': SpinBox(
            #         value=1,
            #         dec=True,
            #     ),
            # },

        
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
        reload(Rabi)

        # create an instance of the ODMR class that implements the experimental logic.
        rabiMeas = Rabi.Rabi_Measurement()

        # run the sweep function in a new thread
        self.sweepProc.run(
            rabiMeas.Rabi,
            self.params_widget.datasetName,
            self.params_widget.samplingFreq,
            self.params_widget.maxIterations,
            self.params_widget.freq,
            int(1e9*self.params_widget.min_MW_time),
            int(1e9*self.params_widget.max_MW_time),
            self.params_widget.num_MW_times,
            self.params_widget.rf_power,
            self.params_widget.laser_power,
            self.params_widget.num_samples,
            int(1e9*self.params_widget.clock_time),
            int(1e9*self.params_widget.init_time),
            int(1e9*self.params_widget.laser_lag),
            int(1e9*self.params_widget.probe_time),
            int(1e9*self.params_widget.singlet_decay),
            # self.params_widget.initial_counts,
            self.params_widget.x_init_position,
            self.params_widget.y_init_position,
            self.params_widget.z_init_position,
            # self.params_widget.threshold,
        )


    def stop(self):
        """Stop the sweep process."""
        
        #kill the TaskVsTime sweep process
        self.sweepProc.kill()



class RabiPlotWidget(LinePlotWidget):

    def __init__(self):
        title = 'Rabi'
        super().__init__(title=title, xlabel='Times (ns)', ylabel='PL')


    def setup(self):
        self.sink = DataSink('RabiData')
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
        mw_time_list = np.sort(self.sink.datasets['mw_times'])
        self.set_data('MW_ON', mw_time_list, [np.mean(self.sink.datasets['mwCountsDict'][mw_time]) for mw_time in mw_time_list])
        self.set_data('MW_OFF', mw_time_list, [np.mean(self.sink.datasets['noMwCountsDict'][mw_time]) for mw_time in mw_time_list])

