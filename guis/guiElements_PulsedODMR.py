"""
GUI for a TaskVsTime Application

Copyright (c) April 2023, Chris Egerstrom
All rights reserved.

This work is licensed under the terms of the 3-Clause BSD license.
For a copy, see <https://opensource.org/licenses/BSD-3-Clause>.
"""

#TO-DO, NEEDS SOME DEBUGGING
from functools import partial
from importlib import reload

from pyqtgraph.Qt import QtWidgets
from pyqtgraph import SpinBox
from nspyre import LinePlotWidget
from nspyre import ParamsWidget
from nspyre import ProcessRunner
from nspyre import DataSink

import experiments.PulsedODMR as PulsedODMR
import numpy as np


class Pulsed_ODMR_Widget(QtWidgets.QWidget):
    
    def __init__(self):
        super().__init__()

        self.setWindowTitle('Pulsed ODMR')

        self.params_widget = ParamsWidget({
            'datasetName': {
                'display_text': 'Dataset Name',
                'widget': QtWidgets.QLineEdit('Pulsed_ODMR'),
            },
            'samplingFreq': {
                'display_text': 'Sampling Freq',
                'widget': SpinBox(
                    value=2e7,
                    suffix='Hz',
                    siPrefix=True,
                    dec=True,
                ),
            },
            'maxIterations': {
                'display_text': 'Iterations',
                'widget': SpinBox(
                    value=10,
                    dec=True,
                    int=True,
                ),
            },
            'startFreq': {
                'display_text': 'Starting Freq',
                'widget': SpinBox(
                    value=2.82e9,
                    suffix='Hz',
                    siPrefix=True, 
                    bounds=(2e9, 4e9),
                    dec=True,
                ),
            },
            'endFreq': {
                'display_text': 'Ending Freq',
                'widget': SpinBox(
                    value=2.92e9,
                    suffix='Hz',
                    siPrefix=True,
                    bounds=(2e9, 4e9),
                    dec=True,
                ),
            },
            'numFreqs': {
                'display_text': 'Number of Freqs',
                'widget': SpinBox(
                    value=11,
                    dec=True,
                    int=True
                ),
            },
            'rfPower': {
                'display_text': 'RF Power',
                'widget': SpinBox(
                    value=0,
                    suffix='dBm',
                    siPrefix=False,
                    bounds=(-42, 6),
                    dec=True,
                ),
            },
            'laser_power': {
                'display_text': 'Laser power/ attenuator voltage',
                'widget': SpinBox(
                    value=2,
                    bounds=(0,5),
                    dec=True,
                ),
            },
            # 'preReadoutLaserAndMwTime': {
            #     'display_text': 'Pre-Meas MW+Laser Time',
            #     'widget': SpinBox(
            #         value=20e-6,
            #         suffix='s',
            #         siPrefix=True,
            #         bounds=(1e-9, 1),
            #         dec=True,
            #     ),
            # },
            # 'laserAndMwReadOutTime': {
            #     'display_text': 'Readout Time',
            #     'widget': SpinBox(
            #         value=50e-6,
            #         suffix='s',
            #         siPrefix=True,
            #         bounds=(1e-9, 1),
            #         dec=True,
            #     ),
            # },
            # 'extraLaserInitTime': {
            #     'display_text': 'Laser Re-init Time',
            #     'widget': SpinBox(
            #         value=20e-6,
            #         suffix='s',
            #         siPrefix=True,
            #         bounds=(1e-9, 1),
            #         dec=True,
            #     ),
            # },
            # 'waitTime': {
            #     'display_text': 'Buffer b/w Meas',
            #     'widget': SpinBox(
            #         value=10e-6,
            #         suffix='s',
            #         siPrefix=True,
            #         bounds=(1e-9, 1),
            #         dec=True,
            #     ),
            # },
            'num_samples': {
                'display_text': 'Samples per frequency point',
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

            'pi_time': {
                'display_text': 'Pi pulse duration',
                'widget': SpinBox(
                    value=100e-9,
                    suffix='s',
                    siPrefix=True, 
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
        reload(PulsedODMR)

        # create an instance of the ODMR class that implements the experimental logic.
        PulsedODMRmeas = PulsedODMR.Pulsed_ODMR_Measurement()

        # self.turnLaserOffAtEndButton.setEnabled(False)

        # run the sweep function in a new thread
        self.sweepProc.run(
            PulsedODMRmeas.PulsedODMR,
            self.params_widget.datasetName, 
            self.params_widget.samplingFreq, 
            self.params_widget.maxIterations,
            self.params_widget.startFreq,
            self.params_widget.endFreq,
            self.params_widget.numFreqs,
            self.params_widget.rfPower,
            self.params_widget.laser_power,
            self.params_widget.num_samples,
            int(1e9*self.params_widget.clock_time),
            int(1e9*self.params_widget.init_time),
            int(1e9*self.params_widget.laser_lag),
            int(1e9*self.params_widget.probe_time),
            int(1e9*self.params_widget.singlet_decay),
            int(1e9*self.params_widget.pi_time),
        )


    def stop(self):
        """Stop the sweep process."""
        #Re-enable checkbox
        # self.turnLaserOffAtEndButton.setEnabled(True)

        #kill the TaskVsTime sweep process
        self.sweepProc.kill()



class PulsedODMRplotWidget(LinePlotWidget):

    def __init__(self):
        title = 'Pulsed ODMR'
        super().__init__(title=title, xlabel='Freq (GHz)', ylabel='Counts') #TODO: Switch this over to contrast


    def setup(self):
        self.sink = DataSink('Pulsed_ODMR')
        self.sink.__enter__() 
        self.add_plot('MW_ON')
        self.add_plot('MW_OFF')
        self.plot_widget.setYRange(-100, 5100)


    def teardown(self):
        self.sink.stop()
        self.sink.__exit__()


    def update(self):
        self.sink.pop() # wait for some data to be saved to sink
        # update the plot
        freq_list = np.sort(self.sink.datasets['freqs'])
        self.set_data('MW_ON', freq_list / 1e9, [np.mean(self.sink.datasets['mwCountsDict'][freq]) for freq in freq_list])
        self.set_data('MW_OFF', freq_list / 1e9, [np.mean(self.sink.datasets['noMwCountsDict'][freq]) for freq in freq_list])


