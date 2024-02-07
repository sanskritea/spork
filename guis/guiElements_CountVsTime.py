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

import experiments.CountVsTime as CountVsTime
from guis.guiElements_general import AutoSaveWidget
from guis.guiElements_general import flexSave


class CustomCountVsTimeWidget(QtWidgets.QWidget):
    
    def __init__(self):
        super().__init__()

        self.setWindowTitle('CountVsTime')

        self.params_widget = ParamsWidget({
            'samplingFreq': {
                'display_text': 'Sampling Freq',
                'widget': SpinBox(
                    value=1,
                    suffix='Hz',
                    siPrefix=True,
                    #bounds=(100e3, 10e9),
                    dec=True,
                ),
            },
            'maxTimeout': {
                'display_text': 'Timeout',
                'widget': SpinBox(
                    value=3600,
                    suffix='s',
                    siPrefix=False, #kiloseconds are cursed
                    #bounds=(100e3, 10e9),
                    dec=True,
                ),
            },
            'datasetName': {
                'display_text': 'Dataset Name',
                'widget': QtWidgets.QLineEdit('CountVsTime'),
            },
            'laser_power': {
                'display_text': 'Laser power / attenuator voltage',
                'widget': SpinBox(
                    value=0.5,
                    bounds=(0, 5),
                    dec=True,
                ),
            },
            'maxIterations': {
                'display_text': 'some iterator',
                'widget': SpinBox(
                    value=1,
                    int=True,
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
        # for checkBox in self.checkBoxes.values():
        #     params_layout.addWidget(checkBox)
        params_layout.addStretch()
        params_layout.addWidget(runButton)
        params_layout.addWidget(stopButton)

        self.autosaveWidget = AutoSaveWidget(False, 30*self.params_widget.samplingFreq) #update once every 30-ish s
        params_layout.addWidget(self.autosaveWidget)

        self.setLayout(params_layout)


    def runClicked(self):
        """Runs when the 'run' button is pressed."""

        # reload the spin measurements module at runtime in case any changes were made to the code
        reload(CountVsTime)

        # create an instance of the ODMR class that implements the experimental logic.
        cvt_meas = CountVsTime.CountVsTimeMeasurement()

        #figure out which channels to care about
        # ctrChanNums = [0]
        # for pfiChanNum, checkBox in self.checkBoxes.items():
        #     if checkBox.isChecked():
        #         ctrChanNums.append(pfiChanNum)
        #     checkBox.setEnabled(False) #disable changing all checkboxes while running

        # run the sweep function in a new thread
        self.sweepProc.run(
            cvt_meas.CountVsTime,
            self.params_widget.datasetName,
            self.params_widget.samplingFreq,
            self.params_widget.laser_power,
            self.params_widget.maxIterations,
            int(self.params_widget.maxTimeout/self.params_widget.samplingFreq),
            # ctrChanNums,
            [self.autosaveWidget.shouldAutosave(), self.autosaveWidget.getAutosaveInterval()], #or None #Autosave params
            #True #runs debug mode. True => TimeVsTime
        )


    def stop(self):
        """Stop the sweep process."""
        #Re-enable all checkboxes
        # for checkBox in self.checkBoxes.values():
        #     checkBox.setEnabled(True)

        #kill the CountVsTime sweep process
        self.sweepProc.kill()
        flexSave(self.params_widget.datasetName, 'CountVsTime', 'closeout') #if process is closed


class CountVsTimePlotWidget(LinePlotWidget):

    def __init__(self):
        title = 'CountVsTime'
        super().__init__(title=title, xlabel='Times', ylabel='Counts')


    def setup(self):
        self.sink = DataSink('CountVsTime')
        self.sink.__enter__() #NVM, Jacob says this is fine #TODO: Dejank this jank
        #not dealing with data here, wait for something to pop in the update file
        self.oldDatasets = []
        self.plot_widget.setYRange(-100, 5100)


    def teardown(self):
        self.sink.__exit__()


    def update(self):
        self.sink.pop() #wait for some data to be saved to sink
        # update the plot
        if self.oldDatasets != self.sink.datasets.keys(): #if different chans are being collected than previous run, re-setup the widget
            self.restartPlot()

        for datasetName in self.sink.datasets:
            if datasetName != 'times':
                self.set_data(f'CountVsTime', self.sink.datasets['times'], self.sink.datasets[datasetName])


    def restartPlot(self):
        #Still not happy wtih some sort of timer threading issues when going from no data => data, but fine after first 'restart' 
        for datasetName in self.oldDatasets: #remove the old plots
            if datasetName != 'times':
                self.remove_plot(f'CountVsTime')
        for datasetName in self.sink.datasets: #add the new plots
            if datasetName != 'times':
                self.add_plot(f'CountVsTime')
        self.oldDatasets = self.sink.datasets.keys() #update this tracking var


