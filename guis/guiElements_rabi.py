
from functools import partial
from importlib import reload

from pyqtgraph.Qt import QtWidgets
from pyqtgraph import SpinBox
from nspyre import LinePlotWidget
from nspyre import ParamsWidget
from nspyre import ProcessRunner
from nspyre import DataSink

import experiments.rabi as rabi


class RabiWidget(QtWidgets.QWidget):
    
    def __init__(self):
        super().__init__()

        self.setWindowTitle('Rabi')

        self.params_widget = ParamsWidget({

            'datasetName': {
                'display_text': 'Dataset Name',
                'widget': QtWidgets.QLineEdit('Rabi'),
            },

           """ 'mwTime': {
                'display_text': 'Tau',
                'widget': SpinBox(
                    value=1,
                    suffix='us',
                    siPrefix=True,
                    #bounds=(100e3, 10e9),
                    dec=True,
                ),
            },"""

            'minMwTime': {
                'display_text': 'Min MW Time',
                'widget': SpinBox(
                    value=0e-9,
                    suffix='s',
                    siPrefix=True,
                    #bounds=(100e3, 10e9),
                    dec=True,
                ),
            },

            'maxMwTime': {
                'display_text': 'Max MW Time',
                'widget': SpinBox(
                    value=1e-6,
                    suffix='s',
                    siPrefix=True,
                    #bounds=(100e3, 10e9),
                    dec=True,
                ),
            },

            'mwNum': {
                'display_text': 'Taus',
                'widget': SpinBox(
                    value=21,
                    #suffix='us',
                    siPrefix=False,
                    #bounds=(100e3, 10e9),
                    dec=True,
                ),
            },

            'delayTime': {
                'display_text': 'delayTime',
                'widget': SpinBox(
                    value=3600,
                    suffix='s',
                    siPrefix=False, #kiloseconds are cursed
                    #bounds=(100e3, 10e9),
                    dec=True,
                ),
            },

            'aomLag': {
                'display_text': 'AOM Lag',
                'widget': SpinBox(
                    value=1e-6,
                    suffix='s',
                    siPrefix=True, 
                    #bounds=(100e3, 10e9),
                    dec=True,
                ),
            },

            'readoutTime': {
                'display_text': 'Readout Time',
                'widget': SpinBox(
                    value=1,
                    suffix='s',
                    siPrefix=True, 
                    #bounds=(100e3, 10e9),
                    dec=True,
                ),
            },

            'initTime': {
                'display_text': 'Readout Time',
                'widget': SpinBox(
                    value=1,
                    suffix='s',
                    siPrefix=True, 
                    #bounds=(100e3, 10e9),
                    dec=True,
                ),
            },

            'samplingRate': {
                'display_text': 'Sampling Rate',
                'widget': SpinBox(
                    value=1,
                    suffix='Hz',
                    siPrefix=True,
                    dec=True,
                ),
            },

            'timePerPoint': {
                'display_text': 'Time/Pt',
                'widget': SpinBox(
                    value=1,
                    suffix='s',
                    siPrefix=True,
                    dec=True,
                ),
            },

            'frequency': {
                'display_text': 'Frequency',
                'widget': SpinBox(
                    value=2.87e9,
                    suffix='Hz',
                    siPrefix=True,
                    dec=True,
                ),
            },

            'laserPower': {
                'display_text': 'Laser Power',
                'widget': SpinBox(
                    value=1e-3,
                    suffix='W',
                    siPrefix=True,
                    bounds=(1e-3, 1e-1),
                    dec=True,
                ),
            },

            'rfPower': {
                'display_text': 'RF Power',
                'widget': SpinBox(
                    value=-17,
                    suffix='dBm',
                    siPrefix=False,
                    bounds=(-42, 0),
                    dec=True,
                ),
            },

            """'mwTimes': {
                'display_text': 'MW Times',
                'widget': SpinBox(
                    value=0e-9,
                    suffix='s',
                    siPrefix=True,
                    bounds=(1e-9, 1e-6),
                    dec=True,
                ),
            },   """ 

            'initTime': {
                'display_text': 'MW Time',
                'widget': SpinBox(
                    value=0e-9,
                    suffix='s',
                    siPrefix=True,
                    bounds=(1e-9, 1e-6),
                    dec=True,
                ),
            },   

            'sweeps': {
                'display_text': 'Sweeps',
                'widget': SpinBox(
                    value=100,
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
        params_layout.addWidget(self.turnLaserOffAtEndButton)
        params_layout.addStretch()
        params_layout.addWidget(runButton)
        params_layout.addWidget(stopButton)

        self.setLayout(params_layout)


    def runClicked(self):
        """Runs when the 'run' button is pressed."""

        # reload the spin measurements module at runtime in case any changes were made to the code
        reload(rabi)

        # create an instance of the ODMR class that implements the experimental logic.
        rabiMeas = rabi.Rabi_Measurement()

        # run the sweep function in a new thread
        self.sweepProc.run(
            rabiMeas.rabi,
            self.params_widget.datasetName,
            self.params_widget.aomLag,
            self.params_widget.delayTime,
            self.params_widget.mwTime,
            self.params_widget.maxMwTime,
            self.params_widget.minMwTime,
            self.params_widget.numMwTimes,
            self.params_widget.readoutTime,
            self.params_widget.initTime,
            self.params_widget.samplingRate,
            self.params_widget.frequency,
            self.params_widget.laserPower,
            self.params_widget.timePerPoint,
            self.params_widget.sweeps
        )


    def stop(self):
        """Stop the sweep process."""
        #Re-enable all checkboxes
        for checkBox in self.checkBoxes.values():
            checkBox.setEnabled(True)

        #kill the TaskVsTime sweep process
        self.sweepProc.kill()



class rabiPlotWidget(LinePlotWidget):
    def __init__(self):
        title = 'Rabi'
        super().__init__(title=title, xlabel='Times', ylabel='PL')


    def setup(self):
        self.sink = DataSink('Rabi')
        self.sink.__enter__() 
        self.add_plot('PL')
        self.plot_widget.setYRange(-100, 5100)



    def teardown(self):
        self.sink.stop()
        self.sink.__exit__()


    def update(self):
        self.sink.pop() #wait for some data to be saved to sink
        # update the plot
        mwTime = np.sort(self.sink.datasets['mwTime'])
        self.set_data('rabi-signal', mwTime, [np.mean(self.sink.datasets[mwTime]['noMwCountsDict']) for freq in freqs])
       # self.set_data('cwODMR-background', freqs, [np.mean(self.sink.datasets['noMwCountsDict'][freq]) for freq in freqs])

