
import time
from itertools import count

import numpy as np

from nspyre import DataSource
from nspyre import InstrumentGateway

from rpyc.utils.classic import obtain

from experiments.NewPulses import Pulses
from drivers.ni.nidaq_final import NIDAQ


class Ramsey_Measurement:

    def Ramsey(
        self,
        datasetName: str,
        samplingFreq: float,
        maxIterations: int,
        freq: float,
        min_tau_time: int,
        max_tau_time: int,
        num_tau_times: int,
        rf_power: float,
        laser_power: float,
        num_samples: int,
        clock_time: int,    # DAQ counting trigger pulse (from Swabian) duration, usually 10ns
        init_time: int,     # NV initialization laser ON duration
        laser_lag: int,     # laser stabilizing time, usually ~100ns
        probe_time: int,    # readout laser ON time
        singlet_decay: int, # NV singlet state emptying duration
        pi_time: int,       # pi pulse duration from Rabi
    ):
        """Run a Ramsey experiment
        Arguments:  *
        """  # TODO: Fill in arg documentation

        # connect to the instrument server
        # connect to the data server and create a data set, or connect to an
        # existing one with the same name if it was created earlier.

        with InstrumentGateway() as gw, DataSource(datasetName) as RamseyData:

            print('RAMSEY PARAMETERS')
            print('clock_time ', clock_time)
            print('init_time ', init_time)
            print('laser_lag ', laser_lag)
            print('probe_time ', probe_time)
            print('singlet_decay ', singlet_decay)

            # set up MW duration list to sweep over
            self.tau_times = np.linspace(min_MW_time, max_MW_time, num_MW_times, endpoint=True)
            # for storing the experiment data
            self.mwCountsDict = dict(
                [[tau_time, []] for tau_time in self.tau_times]
            )
            self.noMwCountsDict = dict([[tau_time, []] for tau_time in self.tau_times])

            # setup a relevant iterator depending on maxIterations
            if maxIterations < 0:
                iters = count()  # infinite iterator
            else:
                iters = range(maxIterations)

            # make sequences beforehand
            seqs = []
            print('clock_time ', clock_time)
            print('init_time ', init_time)
            pi_half_time = int(pi_time / 2)
            for tau_time in self.tau_times:
                seqs.append(Pulses(gw).RAMSEY(clock_time, init_time, laser_lag, probe_time, singlet_decay, pi_half_time, tau_time, max_tau_time))

            # SRS actions
            gw.sg.set_rf_amplitude(rf_power)  # set ouput power
            gw.sg.set_mod_state(True)
            gw.sg.set_mod_type("7")
            gw.sg.set_mod_func("5")
            gw.sg.set_frequency(freq)
            gw.sg.set_rf_state("1")

            # print('test streaming')
            # gw.swabian.runSequenceInfinitely(seqs[0])

            # MAIN EXPERIMENT LOOP
            with NIDAQ() as mynidaq:

                # set laser power
                mynidaq.laser_power_atten(laser_power)

                for i in iters:

                    for j in np.arange(num_MW_times):

                        # START TASK
                        tau_time = self.tau_times[j]
                        print('MW duration ', tau_time, ' ns')
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
                        self.mwCountsDict[tau_time].append(signal_counts)
                        self.noMwCountsDict[tau_time].append(bg_counts)

                        # RESET SWABIAN OUTPUTS
                        gw.swabian.reset()

                        # SAVE CURRENT DATA TO DATA SERVER
                        RamseyData.push(
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
                                "title": "RAMSEY DATA",
                                "xlabel": "MW Times",
                                "ylabel": "Counts",
                                "datasets": {
                                    "tau_times": self.tau_times,
                                    "mwCountsDict": self.mwCountsDict,
                                    "noMwCountsDict": self.noMwCountsDict,
                                },
                            }
                        )

                        

            print("Experiment finished!")

           