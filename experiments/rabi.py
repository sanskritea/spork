
import time
from itertools import count

import numpy as np

from nspyre import DataSource
from nspyre import InstrumentGateway

from rpyc.utils.classic import obtain

from experiments.NewPulses import Pulses
from experiments.NewportSpatialFeedback import SpatialFeedback
from drivers.ni.nidaq_final import NIDAQ


class Rabi_Measurement:

    def Rabi(
        self,
        datasetName: str,
        samplingFreq: float,
        maxIterations: int,
        freq: float,
        min_MW_time: int,
        max_MW_time: int,
        num_MW_times: int,
        rf_power: float,
        laser_power: float,
        num_samples: int,
        clock_time: int,    # DAQ counting trigger pulse (from Swabian) duration, usually 10ns
        init_time: int,     # NV initialization laser ON duration
        laser_lag: int,     # laser stabilizing time, usually ~100ns
        probe_time: int,    # readout laser ON time
        singlet_decay: int, # NV singlet state emptying duration
        # initial_counts: float,
        x_init_position: float, 
        y_init_position: float,
        z_init_position: float,   
        # threshold: float, 
    ):
        """Run a Rabi experiment
        Arguments:  *
        """  # TODO: Fill in arg documentation

        # connect to the instrument server
        # connect to the data server and create a data set, or connect to an
        # existing one with the same name if it was created earlier.

        with InstrumentGateway() as gw, DataSource(datasetName) as RabiData:

            print('RABI PARAMETERS')
            print('clock_time ', clock_time)
            print('init_time ', init_time)
            print('laser_lag ', laser_lag)
            print('probe_time ', probe_time)
            print('singlet_decay ', singlet_decay)

            # set up MW duration list to sweep over
            self.mw_times = np.linspace(min_MW_time, max_MW_time, num_MW_times, endpoint=True)
            # for storing the experiment data
            self.mwCountsDict = dict(
                [[mw_time, []] for mw_time in self.mw_times]
            )
            self.noMwCountsDict = dict([[mw_time, []] for mw_time in self.mw_times])

            # setup a relevant iterator depending on maxIterations
            if maxIterations < 0:
                iters = count()  # infinite iterator
            else:
                iters = range(maxIterations)

            # make sequences beforehand
            seqs = []
            for mw_time in self.mw_times:
                seqs.append(Pulses(gw).RABI(int(mw_time), clock_time, init_time, laser_lag, probe_time, singlet_decay, max_MW_time))

            # SRS actions
            gw.sg.set_rf_amplitude(rf_power)  # set ouput power
            gw.sg.set_frequency(freq)
            gw.sg.set_mod_state(False)
            gw.sg.set_rf_state("1")

            # FEEDBACK PARAMETERS
            self.start_time = time.time()

            # MAIN EXPERIMENT LOOP
            with NIDAQ() as mynidaq:

                # set laser power
                mynidaq.laser_power_atten(laser_power)

                for i in iters:

                    # # SPATIAL FEEDBACK (almost) EVERY 2 HOURS
                    # if (int(time.time() - self.start_time) >= 2):
                    
                    #     # Measure counts after 2 hours
                    #     # turn on laser
                    #     trigger_rate = int(20e3)
                    #     time_per_point = 0.01
                    #     num_samples = int(trigger_rate * time_per_point)
                    #     print('Measuring counts')
                    #     gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(int(trigger_rate)))
                    #     current_counts = np.mean(obtain(mynidaq.internal_read_task(int(trigger_rate), num_samples))) / (1 / trigger_rate)
                    #     gw.swabian.reset()

                    #     # If counts dropped, do spatial feedback
                    #     # if (np.abs(initial_counts - current_counts) / initial_counts)  > threshold:
                    #     print('Feedback')
                    #     x_final_position, y_final_position, z_final_position = SpatialFeedback.Feedback(x_init_position, y_init_position, z_init_position, initial_counts, current_counts, 10)
                    #     print('New locations')
                    #     print('X ', x_final_position, 'mm')
                    #     print('Y ', y_final_position, 'mm')
                    #     print('Z ', z_final_position, 'mm')

                    #     self.start_time = time.time()

                    for j in np.arange(num_MW_times):

                        # START TASK
                        mw_time = self.mw_times[j]
                        # print('MW duration ', mw_time, ' ns')
                        # time.sleep(10)
                        mynidaq.start_external_read_task(samplingFreq, ((4 * num_samples) + 1))

                        # START PULSESTREAMER
                        gw.swabian.runSequenceInfinitely(seqs[j])

                        # COUNTING WITH EXTERNAL TRIGGER
                        # print('counting')
                        raw_counts = (obtain(mynidaq.external_read_task(20e6, ((4 * num_samples) + 1))))

                        # SEPARATING COUNTS INTO MW ON / MW OFF
                        # print('counts ', raw_counts)
                        signal_counts = 1e9 * np.mean(raw_counts[2::4]) / probe_time
                        bg_counts = 1e9 * np.mean(raw_counts[0::4]) / probe_time
                        print(raw_counts[2::4])
                        self.mwCountsDict[mw_time].append(signal_counts)
                        self.noMwCountsDict[mw_time].append(bg_counts)

                        # RESET SWABIAN OUTPUTS
                        gw.swabian.reset()

                        # SAVE CURRENT DATA TO DATA SERVER
                        RabiData.push(
                            {
                                "params": {
                                    "datasetName": datasetName,
                                    "samplingFreq": samplingFreq,
                                    "maxIterations": maxIterations,
                                    "freq": freq,
                                    "min_MW_time": min_MW_time,
                                    "max_MW_time": max_MW_time,
                                    "num_MW_times": num_MW_times,
                                    "rf_power": rf_power,
                                    "laser_power": laser_power,
                                    "num_samples": num_samples,
                                    "clock_time": clock_time,
                                    "init_time": init_time,
                                    "laser_lag": laser_lag,
                                    "probe_time": probe_time,
                                    "singlet_decay": singlet_decay,
                                },
                                "title": "RABI DATA",
                                "xlabel": "MW Times",
                                "ylabel": "Counts",
                                "datasets": {
                                    "mw_times": self.mw_times,
                                    "mwCountsDict": self.mwCountsDict,
                                    "noMwCountsDict": self.noMwCountsDict,
                                },
                            }
                        )

                        
            print("Experiment finished!")

           