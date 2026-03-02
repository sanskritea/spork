
import time
from itertools import count

import numpy as np

from nspyre import DataSource
from nspyre import InstrumentGateway

from rpyc.utils.classic import obtain
from guis.guiElements_general import flexSave

from experiments.NewPulses import Pulses
from experiments.fitter import rabi_fit
from experiments.NewportSpatialFeedback import SpatialFeedback
from drivers.ni.nidaq_final import NIDAQ


class Rabi_Measurement:

    def safe_set_rf_state(self, gw, state, max_retries=3):
        '''Set RF state with retry on I/O error.'''
        for attempt in range(max_retries):
            try:
                gw.sg.set_rf_state(state)
                return  # Success!
            except Exception as e:
                if "VI_ERROR_IO" in str(e) and attempt < max_retries - 1:
                    print(f"Warning: RF state set failed (attempt {attempt+1}/{max_retries}), retrying...")
                    time.sleep(0.5)  # Wait before retry
                    # Try to reset the connection
                    try:
                        gw.sg.close()
                        time.sleep(0.2)
                        gw.sg.initialize()
                    except:
                        pass
                else:
                    raise  # Re-raise if final attempt or different error


    def Rabi(
        self,
        datasetName: str,
        num_samples: int,
        maxIterations: int,
        rf_power: float,
        laser_power: float,
        freq: float,
        min_MW_time: int,
        max_MW_time: int,
        num_MW_times: int,
    ):
        """Run a Rabi experiment
        Arguments:  *
        """  # TODO: Fill in arg documentation

        # connect to the instrument server
        # connect to the data server and create a data set, or connect to an
        # existing one with the same name if it was created earlier.

        with InstrumentGateway() as gw, DataSource(datasetName) as RabiData:

            # NV freq
            print('NV freq : ', freq)

            # set up MW duration list to sweep over
            self.mw_times = np.linspace(min_MW_time, max_MW_time, num_MW_times, endpoint=True)
            # for storing the experiment data
            self.mwCountsDict = dict(
                [[mw_time, []] for mw_time in self.mw_times]
            )
            self.noMwCountsDict = dict([[mw_time, []] for mw_time in self.mw_times])
            probe_time = 500

            # setup a relevant iterator depending on maxIterations
            if maxIterations < 0:
                iters = count()  # infinite iterator
            else:
                iters = range(maxIterations)

            # make sequences beforehand
            seqs = []
            for mw_time in self.mw_times:
                seqs.append(Pulses(gw).RABI(int(mw_time), max_MW_time))

            # SRS actions
            gw.sg.set_rf_amplitude(rf_power)  # set ouput power
            gw.sg.set_frequency(freq)
            gw.sg.set_mod_state(False)
            gw.sg.set_rf_state('1')  # turn on RF, will be turned off at end of each iteration in case of error to prevent overheating
            # self.safe_set_rf_state(gw, "1")

            # FEEDBACK PARAMETERS
            self.start_time = time.time()

            # MAIN EXPERIMENT LOOP
            with NIDAQ() as mynidaq:

                # set laser power
                mynidaq.laser_power_atten(laser_power)

                for i in iters:
                    
                    gw.sg.set_rf_state('1')  # turn on RF, will be turned off at end of each iteration in case of error to prevent overheating


                    for j in np.arange(num_MW_times):

                        # START TASK
                        mw_time = self.mw_times[j]
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
                        self.mwCountsDict[mw_time].append(signal_counts)
                        self.noMwCountsDict[mw_time].append(bg_counts)

                        # RESET SWABIAN OUTPUTS
                        gw.swabian.reset()

                        # SAVE CURRENT DATA TO DATA SERVER
                        RabiData.push(
                            {
                                "params": {
                                    "datasetName": datasetName,
                                    "num_samples": num_samples,
                                    "maxIterations": maxIterations,
                                    "rf_power": rf_power,
                                    "laser_power": laser_power,
                                    "freq": freq,
                                    "min_MW_time": min_MW_time,
                                    "max_MW_time": max_MW_time,
                                    "num_MW_times": num_MW_times,
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
                    
                    # TURN OFF SRS
                    gw.sg.set_rf_state('0')  # turn off RF
                    time.sleep(0.1)  # wait for RF to turn off before next iteration to prevent overheating

                notes = ''
                flexSave(datasetName, notes, 'Rabi')

                # Run fitter to find peaks
                pi_time, fit, tau = rabi_fit(self.mw_times, self.mwCountsDict, self.noMwCountsDict)
                print('Fitted pi-time ', np.abs(pi_time))
                return np.abs(pi_time)

            print("Experiment finished!")

           