
from functools import partial
from importlib import reload
import numpy as np

from pyqtgraph.Qt import QtWidgets
from pyqtgraph import SpinBox
from nspyre import LinePlotWidget
from nspyre import ParamsWidget
from nspyre import ProcessRunner
from nspyre import DataSink

import experiments.ReinitTest as ReinitTest


class Reinit_Test_Widget(QtWidgets.QWidget):
    
    def __init__(self):
        super().__init__()

        self.setWindowTitle('Optical T1')

        self.params_widget = ParamsWidget({

            'datasetName': {
                'display_text': 'Dataset Name',
                'widget': QtWidgets.QLineEdit('Reinit Test'),
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
                    value=10,
                    dec=True,
                    int=True,
                ),
            },

            'tau_min': {
                'display_text': 'Min tau delay',
                'widget': SpinBox(
                    value=0e-9,
                    suffix='s',
                    siPrefix=True,
                    #bounds=(100e3, 10e9),
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
                'display_text': 'Samples per tau',
                'widget': SpinBox(
                    value=10000,
                    dec=True,
                    int=True
                ),
            },

            'clock_time': {
                'display_text': 'Clock Pulse ',
                'widget': SpinBox(
                    value=10e-9,
                    suffix='s',
                    siPrefix=True,
                    bounds=(1e-9, 1),
                    dec=True,
                ),
            },

            'init_time': {
                'display_text': 'NV Init Pulse',
                'widget': SpinBox(
                    value=10e-6, # 1-2us at saturation, can go longer
                    suffix='s',
                    siPrefix=True, #kiloseconds are cursed
                    bounds=(0, 1),
                    dec=True,
                ),
            },

            'readout_time': {
                'display_text': 'Readout Pulse',
                'widget': SpinBox(
                    value=50e-9, #300ns ish
                    suffix='s',
                    siPrefix=True, 
                    bounds=(0, 1),
                    dec=True,
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
        reload(ReinitTest)

        # create an instance of the ODMR class that implements the experimental logic.
        ReinitTestMeas = ReinitTest.Reinit_Test_Measurement()

        # run the sweep function in a new thread
        self.sweepProc.run(
            ReinitTestMeas.ReinitTest,
            self.params_widget.datasetName,
            self.params_widget.samplingFreq,
            self.params_widget.maxIterations,
            int(1e9*self.params_widget.tau_min),
            self.params_widget.laser_power,
            self.params_widget.num_samples,
            int(1e9*self.params_widget.clock_time),
            int(1e9*self.params_widget.init_time),
            int(1e9*self.params_widget.readout_time),
        )


    def stop(self):
        """Stop the sweep process."""
        
        #kill the TaskVsTime sweep process
        self.sweepProc.kill()



class ReinitTestPlotWidget(LinePlotWidget):

    def __init__(self):
        title = 'Reinit Test'
        super().__init__(title=title, xlabel='Time (us)', ylabel='PL')


    def setup(self):
        self.sink = DataSink('Reinit Test')
        self.sink.__enter__() 
        self.add_plot('Bright Counts')
        self.add_plot('Dark Counts')
        self.plot_widget.setYRange(-100, 5100)


    def teardown(self):
        self.sink.stop()
        self.sink.__exit__()


    def update(self):
        self.sink.pop() #wait for some data to be saved to sink
        # update the plot
        tau_list = np.sort(self.sink.datasets['taus']) 
        self.set_data('Bright Counts', tau_list / 1000, [np.mean(self.sink.datasets['BrightCounts'][tau]) for tau in tau_list])
        self.set_data('Dark Counts', tau_list / 1000, [np.mean(self.sink.datasets['DarkCounts'][tau]) for tau in tau_list])



