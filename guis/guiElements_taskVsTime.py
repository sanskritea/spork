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

import experiments.taskVsTime as taskVsTime
from guis.guiElements_general import AutoSaveWidget
from guis.guiElements_general import flexSave


class CustomTaskVsTimeWidget(QtWidgets.QWidget):
    
    def __init__(self):
        super().__init__()

        self.setWindowTitle('Task Vs Time')

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
                'widget': QtWidgets.QLineEdit('TaskVsTime'),
            },
        })

        #Set-up check boxes for each channel
        # ifCtr0Check = QtWidgets.QCheckBox('Counts Chan0')
        # ifCtr0Check.setChecked(True)
        # ifCtr1Check = QtWidgets.QCheckBox('Counts Chan1')
        # ifCtr1Check.setChecked(False)
        # ifCtr2Check = QtWidgets.QCheckBox('Counts Chan2')
        # ifCtr2Check.setChecked(False)
        # ifCtr3Check = QtWidgets.QCheckBox('Counts Chan3')
        # ifCtr3Check.setChecked(False)
        # #wrap the check boxes together nicely, indexed by PFI channel number
        # self.checkBoxes = {0: ifCtr0Check, 1: ifCtr1Check, 2: ifCtr2Check, 3: ifCtr3Check}
        
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
        reload(taskVsTime)

        # create an instance of the ODMR class that implements the experimental logic.
        tvt_meas = taskVsTime.TaskVsTimeMeasurement()

        #figure out which channels to care about
        ctrChanNums = [0]
        # for pfiChanNum, checkBox in self.checkBoxes.items():
        #     if checkBox.isChecked():
        #         ctrChanNums.append(pfiChanNum)
        #     checkBox.setEnabled(False) #disable changing all checkboxes while running

        # run the sweep function in a new thread
        self.sweepProc.run(
            tvt_meas.taskVsTime,
            self.params_widget.datasetName,
            self.params_widget.samplingFreq,
            int(self.params_widget.maxTimeout/self.params_widget.samplingFreq),
            ctrChanNums,
            [self.autosaveWidget.shouldAutosave(), self.autosaveWidget.getAutosaveInterval()], #or None #Autosave params
            #True #runs debug mode. True => TimeVsTime
        )


    def stop(self):
        """Stop the sweep process."""
        #Re-enable all checkboxes
        # for checkBox in self.checkBoxes.values():
        #     checkBox.setEnabled(True)

        #kill the TaskVsTime sweep process
        self.sweepProc.kill()
        flexSave(self.params_widget.datasetName, 'TvT', 'closeout') #if process is closed



class TaskVsTimePlotWidget(LinePlotWidget):
    def __init__(self):
        title = 'Task vs Time'
        super().__init__(title=title, xlabel='Times', ylabel='Counts')


    def setup(self):
        self.sink = DataSink('TaskVsTime')
        self.sink.__enter__() #NVM, Jacob says this is fine #TODO: Dejank this jank
        #not dealing with data here, wait for something to pop in the update file
        self.oldDatasets = []
        self.plot_widget.setYRange(-100, 5100)


    def teardown(self):
        self.sink.__exit__()


    def update(self):
        self.sink.pop() #wait for some data to be saved to sink
        # update the plot
        # if self.oldDatasets != self.sink.datasets.keys(): #if different chans are being collected than previous run, re-setup the widget
            
        #     self.restartPlot()
        for datasetName in self.sink.datasets:
            # if datasetName != 'times':
            #     self.set_data(f'task_vs_time-{datasetName}', self.sink.datasets['times'], self.sink.datasets[datasetName])
            if datasetName != 'times':
                self.set_data(f'task_vs_time-{datasetName}', self.sink.datasets['times'], self.sink.datasets[datasetName])


    def restartPlot(self):
        #Still not happy wtih some sort of timer threading issues when going from no data => data, but fine after first 'restart' 
        for datasetName in self.oldDatasets: #remove the old plots
            if datasetName != 'times':
                self.remove_plot(f'task_vs_time-{datasetName}')
        for datasetName in self.sink.datasets: #add the new plots
            if datasetName != 'times':
                self.add_plot(f'task_vs_time-{datasetName}')
        self.oldDatasets = self.sink.datasets.keys() #update this tracking var


