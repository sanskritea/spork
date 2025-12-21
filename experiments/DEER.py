
import time
from itertools import count

import numpy as np

from nspyre import DataSource
from nspyre import InstrumentGateway

from rpyc.utils.classic import obtain
from guis.guiElements_general import flexSave

from experiments.NewPulses import Pulses
# from experiments.fitter import DEER
from experiments.NewportSpatialFeedback import SpatialFeedback
from drivers.ni.nidaq_final import NIDAQ


class DEER_Measurement:

    def DEER(
        self,
        datasetName: str,
        num_samples: int,
        maxIterations: int,
        laser_power: float,
        NV_freq: float,
        NV_rf_power: float,
        NV_pi_time: float,
        P1_freq: float,
        P1_rf_power: float,
        P1_pi_time: float,
        tau_min: int,
        tau_max: int,
        tau_num: int,
    ):
        """Run a DEER experiment
        Arguments:  *
        """  # TODO: Fill in arg documentation

        # connect to the instrument server
        # connect to the data server and create a data set, or connect to an
        # existing one with the same name if it was created earlier.

        with InstrumentGateway() as gw, DataSource(datasetName) as DEERData:

            # set up MW duration list to sweep over
            self.tau_times = np.linspace(tau_min, tau_max, tau_num, endpoint=True)
            # for storing the experiment data
            self.mwCountsDict = dict(
                [[tau_time, []] for tau_time in self.tau_times]
            )
            self.noMwCountsDict = dict([[tau_time, []] for tau_time in self.tau_times])
            probe_time = 500

            # setup a relevant iterator depending on maxIterations
            if maxIterations < 0:
                iters = count()  # infinite iterator
            else:
                iters = range(maxIterations)

            # make sequences beforehand
            seqs = []
            for tau_time in self.tau_times:
                seqs.append(Pulses(gw).DEER(NV_pi_time, P1_pi_time, tau_time, tau_max))

            # SRS actions
            gw.sg.set_rf_amplitude(NV_rf_power)  # set ouput power
            gw.sg.set_frequency(NV_freq)
            gw.sg.set_mod_state(False)
            gw.sg.set_rf_state("1")

            gw.sg2.set_rf_amplitude(P1_rf_power)  # set ouput power
            gw.sg2.set_frequency(P1_freq)
            gw.sg2.set_mod_state(False)
            gw.sg2.set_rf_state("1")

            # FEEDBACK PARAMETERS
            self.start_time = time.time()

            # MAIN EXPERIMENT LOOP
            with NIDAQ() as mynidaq:

                # set laser power
                mynidaq.laser_power_atten(laser_power)

                for i in iters:

                    for j in np.arange(tau_min):

                        # START TASK
                        tau_time = self.tau_times[j]
                        # print('MW duration ', mw_time, ' ns')
                        # time.sleep(10)
                        mynidaq.start_external_read_task(20e6, ((4 * num_samples) + 1))

                        # START PULSESTREAMER
                        gw.swabian.runSequenceInfinitely(seqs[j])

                        # COUNTING WITH EXTERNAL TRIGGER
                        # print('counting')
                        raw_counts = (obtain(mynidaq.external_read_task(20e6, ((4 * num_samples) + 1))))

                        # SEPARATING COUNTS INTO MW ON / MW OFF
                        # print('counts ', raw_counts)
                        signal_counts = 1e9 * np.mean(raw_counts[2::4]) / probe_time
                        bg_counts = 1e9 * np.mean(raw_counts[0::4]) / probe_time
                        self.mwCountsDict[tau_time].append(signal_counts)
                        self.noMwCountsDict[tau_time].append(bg_counts)

                        # RESET SWABIAN OUTPUTS
                        gw.swabian.reset()

                        # SAVE CURRENT DATA TO DATA SERVER
                        DEERData.push(
                            {
                                "params": {
                                    "datasetName": datasetName,
                                    "num_samples": num_samples,
                                    "maxIterations": maxIterations,
                                    "laser_power": laser_power,
                                    "NV_freq": NV_freq,
                                    "NV_rf_power": NV_rf_power,
                                    "NW_pi_time": NV_pi_time,
                                    "P1_freq": P1_freq,
                                    "P1_rf_power": P1_rf_power,
                                    "P1_pi_time": P1_pi_time,
                                    "tau_min": tau_min,
                                    "tau_max": tau_max,
                                    "tau_num": tau_num,
                                },
                                "title": "DEER DATA",
                                "xlabel": "Tau Times",
                                "ylabel": "Counts",
                                "datasets": {
                                    "tau_times": self.tau_times,
                                    "mwCountsDict": self.mwCountsDict,
                                    "noMwCountsDict": self.noMwCountsDict,
                                },
                            }
                        )

                notes = ''
                flexSave(datasetName, notes, 'DEER')

                # # Run fitter to find peaks
                # pi_time = rabi(self.mw_times, self.mwCountsDict, self.noMwCountsDict)
                # # print('Fitted pi-time ', pi_time)
                # return pi_time

            print("Experiment finished!")

           