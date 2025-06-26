"""

XY Scanning through DAQ controlled ANM300 Attocube Scanner

Sanskriti Chitransh, 2023-Oct-25

Note: The ANC300/ANM200 Scanner has DC input control enabled from the inserv file and only the DAQ's analog output drives the scanner. 

"""

import numpy as np
import time as time

from nspyre import DataSource
from nspyre import InstrumentGateway
from nspyre import StreamingList

from guis.guiElements_general import flexSave
from rpyc.utils.classic import obtain

from drivers.ni.nidaq_final import NIDAQ
from experiments.NewPulses import Pulses


class XZScan:

    def scanning(
        self,
        datasetname: str,
        trigger_rate: float,
        device: str,
        x_init_position: float,
        z_init_position: float,
        position_steps: int,
        step_length: float,
        time_per_pixel: float,
        step_wait: float,
    ):
        """

        x_init_position: inital X actuator position
        x_final_position: last X actuator position to move to
        z_init_position: inital Y actuator position
        z_final_position: last Y actuator position to move to
        x_position_steps: steps in X direction
        z_position_steps: steps in Y direction
        time_per_pixel: photon clicks for each (X,Y) location

        """

        with InstrumentGateway() as gw, DataSource(datasetname) as NewportXZScanningData: 

            pos_len = int(2 * position_steps + 1)

            # position lists
            x_position_list = np.linspace(
                x_init_position - (position_steps * step_length),
                x_init_position + (position_steps * step_length),
                pos_len,
            )
            z_position_list = np.linspace(
                z_init_position - (position_steps * step_length),
                z_init_position + (position_steps * step_length),
                pos_len,
            )
            # to counter the backlash, let's create a position list that overshoots the z_position_list (the inner loop) where the stepper can idle and then move to the actual starting y position

            print("x_position_list ", x_position_list)
            print("z_position_list ", z_position_list)

            counts = np.ones((pos_len, pos_len))

            with NIDAQ() as mynidaq:

                # turn on laser
                # trigger_rate = 20e3
                gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(int(trigger_rate)))

                for n in range(pos_len):

                    ## MOVE X
                    # print('moving to ', x_position_list[n])
                    gw.esp.espZ.move_to(z_position_list[n])
                    time.sleep(step_wait)
                    start_linescan_time = time.time()

                    # compensate for X backlash
                    # gw.esp.espX.move_to(x_position_list[0] - 0.025)
                    gw.esp.espY.move_to(x_position_list[0] - 0.025)
                    time.sleep(step_wait)
                    x_backlash = np.arange(x_position_list[0] - 0.025, x_position_list[0] - 0.001, 0.001)
                    for xx in x_backlash:
                        # gw.esp.espX.move_to(xx)
                        gw.esp.espY.move_to(xx)
                        time.sleep(step_wait)


                    for nn in range(pos_len):

                        ## MOVE Y
                        # print('moving to ', z_position_list[nn])
                        # gw.esp.espX.move_to(x_position_list[nn])
                        gw.esp.espY.move_to(x_position_list[nn])
                        time.sleep(step_wait)

                        ## COUNTING
                        ## start Swabian external trigger for counting
                        daq_sampling_rate = 20e6

                        ## start DAQ counting tasks
                        # sampling_rate =  trigger_rate
                        # num_samples = int(time_per_pixel * daq_sampling_rate)
                        num_samples = int(time_per_pixel * trigger_rate)
                        reading_period = 1 / trigger_rate

                        counts[n][nn] = np.mean(obtain(mynidaq.internal_read_task(trigger_rate, num_samples))) / reading_period

                        ## PUSH DATA
                        # print("start push")
                        NewportXZScanningData.push(
                            {
                                "params": {
                                    "datasetname": datasetname,
                                    "trigger_rate": trigger_rate,
                                    "device": device,
                                    "x_init_position": x_init_position,
                                    "z_init_position": z_init_position,
                                    "position_steps": position_steps,
                                    "step_length": step_length,
                                    "time_per_pixel": time_per_pixel,
                                },
                                "title": "XZScanning",
                                "xlabel": "X position",
                                "ylabel": "Z position",
                                "datasets": {
                                    "x_position": x_position_list,
                                    "z_position": z_position_list,
                                    "counts": counts,
                                },
                            }
                        )

            # Save data 
            flexSave(datasetname, 'YZ', 'final')

            # go back to starting point
            print('Moving to start location')
            # compensate for X backlash
            # gw.esp.espX.move_to(x_init_position - 0.025)
            gw.esp.espY.move_to(x_init_position - 0.025)
            time.sleep(step_wait)
            x_backlash = np.arange(x_init_position - 0.025, x_init_position + 0.001, 0.001)
            for xx in x_backlash:
                # gw.esp.espX.move_to(xx)
                gw.esp.espY.move_to(xx)
                time.sleep(step_wait)
            gw.esp.espZ.move_to(z_init_position)

            print("XZ-scan finished")
