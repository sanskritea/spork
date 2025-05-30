"""
GUI for a CountVsTime Application

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

import experiments.SRSTesting as SRSTesting
from guis.guiElements_general import AutoSaveWidget
from guis.guiElements_general import flexSave


class SRS_Testing_Widget(QtWidgets.QWidget):
    
    def __init__(self):
        super().__init__()

        self.setWindowTitle('SRS Testing')

        self.params_widget = ParamsWidget({

            'freq': {
                'display_text': 'SRS Freq',
                'widget': SpinBox(
                    value=2.5e9,
                    suffix='Hz',
                    siPrefix=True,
                    dec=True,
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

            'I_On': {
                'display_text': 'I_On',
                'widget': SpinBox(
                    value=0.5,
                    bounds=(-1, 1),
                    dec=True,
                ),
            },

            'I_Off': {
                'display_text': 'I_Off',
                'widget': SpinBox(
                    value=0.5,
                    bounds=(-1, 1),
                    dec=True,
                ),
            },

            'Q_On': {
                'display_text': 'Q_On',
                'widget': SpinBox(
                    value=0.5,
                    bounds=(-1, 1),
                    dec=True,
                ),
            },

            'Q_Off': {
                'display_text': 'Q_Off',
                'widget': SpinBox(
                    value=0.5,
                    bounds=(-1, 1),
                    dec=True,
                ),
            },

            'period': {
                'display_text': 'MW period',
                'widget': SpinBox(
                    value=0.5,
                    suffix='s',
                    siPrefix=True,
                    dec=True,
                ),
            },

        })

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
        params_layout.addStretch()
        params_layout.addWidget(runButton)
        params_layout.addWidget(stopButton)

        self.setLayout(params_layout)


    def runClicked(self):
        """Runs when the 'run' button is pressed."""

        # reload the spin measurements module at runtime in case any changes were made to the code
        reload(SRSTesting)

        # create an instance of the ODMR class that implements the experimental logic.
        srstestingMeas = SRSTesting.SRS_Testing_Measurement()

        # run the sweep function in a new thread
        self.sweepProc.run(
            srstestingMeas.SRSTesting,
            self.params_widget.freq,
            self.params_widget.rf_power,
            self.params_widget.I_On,
            self.params_widget.I_Off,
            self.params_widget.Q_On,
            self.params_widget.Q_Off,
            int(1e9*self.params_widget.period),
        )


    def stop(self):
        """Stop the sweep process."""
        #Re-enable all checkboxes
        # for checkBox in self.checkBoxes.values():
        #     checkBox.setEnabled(True)

        #kill the CountVsTime sweep process
        self.sweepProc.kill()



