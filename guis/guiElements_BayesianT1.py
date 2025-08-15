
from functools import partial
from importlib import reload
import numpy as np

from pyqtgraph.Qt import QtWidgets
from pyqtgraph import SpinBox
from nspyre import LinePlotWidget
from nspyre import ParamsWidget
from nspyre import ProcessRunner
from nspyre import DataSink

import experiments.BayesianT1 as BayesianT1


class BayesianT1_Widget(QtWidgets.QWidget):
    
    def __init__(self):
        super().__init__()

        self.setWindowTitle('BayesianT1')

        self.params_widget = ParamsWidget({

            'datasetName': {
                'display_text': 'Dataset Name',
                'widget': QtWidgets.QLineEdit('BayesianT1'),
            },


            'num_samples': {
                'display_text': 'Samples per tau',
                'widget': SpinBox(
                    value=100,
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

            'Lower_Freq': {
                'display_text': 'Lower NV Freq',
                'widget': SpinBox(
                    value=2.82e9,
                    suffix='Hz',
                    siPrefix=True, 
                    bounds=(2e9, 4e9),
                    dec=True,
                ),
            },

            'Lower_Pi': {
                'display_text': 'Lower Pi pulse',
                'widget': SpinBox(
                    value=100e-9,
                    suffix='s',
                    siPrefix=True, 
                    dec=True,
                ),
            },

            'Higher_Freq': {
                'display_text': 'Higher NV Freq',
                'widget': SpinBox(
                    value=2.82e9,
                    suffix='Hz',
                    siPrefix=True, 
                    bounds=(2e9, 4e9),
                    dec=True,
                ),
            },

            'Higher_Pi': {
                'display_text': 'Higher Pi pulse',
                'widget': SpinBox(
                    value=100e-9,
                    suffix='s',
                    siPrefix=True, 
                    dec=True,
                ),
            },

            'bayesian_iterations': {
                'display_text': 'Bayesian cycles',
                'widget': SpinBox(
                    value=100,
                    dec=True,
                    int=True
                ),
            },

            })

        #Set-up check boxes for each channel
        
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
        reload(BayesianT1)

        # create an instance of the ODMR class that implements the experimental logic.
        BayesianT1Meas = BayesianT1.Bayesian_T1_Meas()

        # run the sweep function in a new thread
        self.sweepProc.run(
            BayesianT1Meas.BayesianT1,
            self.params_widget.datasetName,
            self.params_widget.num_samples,
            self.params_widget.rf_power,
            self.params_widget.laser_power,
            self.params_widget.Lower_Freq,
            int(1e9*self.params_widget.Lower_Pi),
            self.params_widget.Higher_Freq,
            int(1e9*self.params_widget.Higher_Pi),
            self.params_widget.bayesian_iterations,
        )


    def stop(self):
        """Stop the sweep process."""
        
        #kill the TaskVsTime sweep process
        self.sweepProc.kill()



class BayesianT1PlotWidget(LinePlotWidget):

    def __init__(self):
        title = 'Bayesian T1'
        super().__init__(title=title, xlabel='Iteration number')#, ylabel='AU')


    def setup(self):
        self.sink = DataSink('BayesianT1')
        self.sink.__enter__() 
        self.add_plot('M_plus')
        self.add_plot('M_minus')
        self.add_plot('GammaPlus')
        self.add_plot('GammaMinus')


    def teardown(self):
        self.sink.stop()
        self.sink.__exit__()


    def update(self):
        self.sink.pop() #wait for some data to be saved to sink
        # update the plot
        self.set_data('M_plus', self.sink.datasets['iteration_number'], self.sink.datasets['M_plus'])
        self.set_data('M_minus', self.sink.datasets['iteration_number'], self.sink.datasets['M_minus'])
        self.set_data('GammaPlus', self.sink.datasets['iteration_number'], self.sink.datasets['GammaPlus'])
        self.set_data('GammaMinus', self.sink.datasets['iteration_number'], self.sink.datasets['GammaMinus'])

