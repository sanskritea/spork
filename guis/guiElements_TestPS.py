from functools import partial
from importlib import reload

from pyqtgraph.Qt import QtWidgets
from pyqtgraph import SpinBox
from nspyre import LinePlotWidget
from nspyre import ParamsWidget
from nspyre import ProcessRunner
from nspyre import DataSink

import experiments.TestPS as TestPS
from guis.guiElements_general import AutoSaveWidget
from guis.guiElements_general import flexSave


class CustomTestPSWidget(QtWidgets.QWidget):
    
    def __init__(self):
        super().__init__()

        self.setWindowTitle('Testing PS')

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
        # params_layout.addWidget(self.params_widget)
        # # for checkBox in self.checkBoxes.values():
        # #     params_layout.addWidget(checkBox)
        # params_layout.addStretch()
        params_layout.addWidget(runButton)
        params_layout.addWidget(stopButton)

        # self.autosaveWidget = AutoSaveWidget(False, 30*self.params_widget.samplingFreq) #update once every 30-ish s
        # params_layout.addWidget(self.autosaveWidget)

        self.setLayout(params_layout)


    def runClicked(self):
        """Runs when the 'run' button is pressed."""

        # reload the spin measurements module at runtime in case any changes were made to the code
        reload(TestPS)

        # create an instance of the ODMR class that implements the experimental logic.
        cvt_meas = TestPS.TestPSMeasurement()

        #figure out which channels to care about
        # ctrChanNums = [0]
        # for pfiChanNum, checkBox in self.checkBoxes.items():
        #     if checkBox.isChecked():
        #         ctrChanNums.append(pfiChanNum)
        #     checkBox.setEnabled(False) #disable changing all checkboxes while running

        # run the sweep function in a new thread
        self.sweepProc.run(
            cvt_meas.TestPS,
            # self.params_widget.datasetName,
            # self.params_widget.samplingFreq,
            # int(self.params_widget.maxTimeout/self.params_widget.samplingFreq),
            # ctrChanNums,
            # [self.autosaveWidget.shouldAutosave(), self.autosaveWidget.getAutosaveInterval()], #or None #Autosave params
            #True #runs debug mode. True => TimeVsTime
        )


    def stop(self):
        """Stop the sweep process."""
        #Re-enable all checkboxes
        # for checkBox in self.checkBoxes.values():
        #     checkBox.setEnabled(True)

        #kill the TestPS sweep process
        self.sweepProc.kill()
        # flexSave(self.params_widget.datasetName, 'CvT', 'closeout') #if process is closed

