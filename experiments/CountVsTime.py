"""
This is a basic CountVsTime Application

Copyright (c) April 2023, C. Egerstrom
All rights reserved.

This work is licensed under the terms of the 3-Clause BSD license.
For a copy, see <https://opensource.org/licenses/BSD-3-Clause>.
"""

import time
from itertools import count

import numpy as np

from nspyre import DataSource
from nspyre import InstrumentGateway
from nspyre import StreamingList

from rpyc.utils.classic import obtain

from guis.guiElements_general import flexSave
from experiments.NewPulses import Pulses
from experiments.NewportSpatialFeedback import SpatialFeedback
from drivers.ni.nidaq_final import NIDAQ


class CountVsTimeMeasurement:

    def CountVsTime(self, datasetName: str, samplingFreq: float, trigger_rate: float, laser_power: float, time_per_point: float, x_init_position: float, 
        y_init_position: float,
        z_init_position: float): #,   
        # threshold: float):
        # , autosaveParams=None, debug=False):
        """Run a CountVsTime2 experiment

        Args:
            datasetName: name of the dataset to push data to
            sampleFreq (float): how quickly to read data (in Hz)
            ctrChanNums: Which PFI channels to read. Default is [11,1,4,8]
            autosaveParams: Default: None, but will take a list of [shouldAutosave, autosaveInterval] 
            debug: optional (default False), will run TimeVsTime if true
        """

        # connect to the instrument server
        # connect to the data server and create a data set, or connect to an
        # existing one with the same name if it was created earlier.
        with InstrumentGateway() as gw, DataSource(datasetName) as cvt_data:

            # turn on laser
            gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(int(trigger_rate)))

            # for storing the experiment data
            self.times = np.zeros(0)
            self.counts = np.zeros(0)

            with NIDAQ() as mynidaq:

                # set laser power
                mynidaq.laser_power_atten(laser_power)

                # start_time
                self.start_time = time.time() #take time diff measurements
                self.outer_start_time = self.start_time

                while True:

                    # start DAQ counting
                    num_samples = int(trigger_rate * time_per_point)
                    raw_counts = obtain(mynidaq.internal_read_task(int(trigger_rate), num_samples))

                    # Read and save new time
                    elapsed_time = time.time() - self.start_time
                    self.times = np.append(self.times, elapsed_time)
                    trigger_period = 1 / trigger_rate
                    self.counts = [np.append(self.counts, (np.mean(raw_counts / trigger_period)))]
                    # print('counts std ', np.std(raw_counts / trigger_period))

                    # save the current data to the data server.
                    cvt_data.push({'params': {},
                                    'title': 'CountVsTime',
                                    'xlabel': 'Time (s)',
                                    'ylabel': 'Counts',
                                    # 'datasets': {'times': self.times, 'counts': self.counts}
                                    'datasets': dict([('CountVsTime', self.counts[j]) if j<1 else ('times', self.times) for j in range(2)])
                    })

                    # SPATIAL FEEDBACK (almost) EVERY 10-15mins (IDEALLY)
                    # FORGET THE THRESHOLD, JUST SCAN EVERY FEW MINUTES/SECONDS/HOURS AND COME TO MAX LOCATION 
                    feedback_time = int(time.time() - self.outer_start_time)
                    # print('feedback_time ', feedback_time)
                    # if (feedback_time >= 300):
                    
                    #     # Parameters
                    #     print('elapsed_time ', elapsed_time)
                    #     begin_feedback = time.time()
                    #     # mean_initial_counts = np.mean(self.counts[0][0:100])
                    #     # std_initial_counts = np.std(self.counts[0][0:100])
                    #     # threshold = 100 * std_initial_counts / mean_initial_counts

                    #     # Perform feedback
                    #     # print('Feedback')
                    #     SpatialFeedback.Feedback(x_init_position, y_init_position, z_init_position)

                    #     # Measure counts right after feedback
                    #     gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(int(trigger_rate)))
                    #     num_samples = int(trigger_rate * time_per_point)
                    #     post_feedback_counts = np.mean(obtain(mynidaq.internal_read_task(int(trigger_rate), num_samples))) / (1 / trigger_rate)
                    #     print('Post feedback counts ', post_feedback_counts)
                    #     gw.swabian.reset()

                    #     # Feedback closeout
                    #     self.outer_start_time = time.time()
                    #     gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(int(trigger_rate)))
                    #     feedback_duration = time.time() - begin_feedback
                    #     # print('Feedback duration: ', feedback_duration)
                    #     self.start_time = self.start_time + feedback_duration

            #         if shouldAutosave and (i+1)%autosaveInterval == 0: # Autosave logic, +1 so it doesn't autosave first data point
            #             flexSave(datasetName, 'CountVsTime', 'autosave')

            # flexSave(datasetName, 'CountVsTime', 'final') # after measurement finishes


if __name__ == '__main__':
    exp = CountVsTimeMeasurement()
    print('Running CountVsTime with 1Hz sampling rate for max 1hr, saving to CountVsTime on dataserv')
    exp.CountVsTime('CountVsTime', 1, 3600)

