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


class XYScan:

    def scanning(
        self,
        datasetname: str,
        trigger_rate: float,
        device: str,
        x_init_position: float,
        y_init_position: float,
        position_steps: int,
        step_length: float,
        time_per_pixel: float,
        step_wait: float,
    ):
        """

        x_init_position: inital X actuator position
        x_final_position: last X actuator position to move to
        y_init_position: inital Y actuator position
        y_final_position: last Y actuator position to move to
        x_position_steps: steps in X direction
        y_position_steps: steps in Y direction
        time_per_pixel: photon clicks for each (X,Y) location

        """

        with InstrumentGateway() as gw, DataSource(datasetname) as NewportXYScanningData: 

            pos_len = int(2 * position_steps + 1)

            # position lists
            x_position_list = np.linspace(
                x_init_position - (position_steps * step_length),
                x_init_position + (position_steps * step_length),
                pos_len,
            )
            y_position_list = np.linspace(
                y_init_position - (position_steps * step_length),
                y_init_position + (position_steps * step_length),
                pos_len,
            )

            print("x_position_list ", x_position_list)
            print("y_position_list ", y_position_list)

            counts = np.ones((pos_len, pos_len))

            # compensate for Y backlash
            gw.esp.espY.move_to(y_position_list[0] - 0.025)
            time.sleep(step_wait)
            y_backlash = np.arange(y_position_list[0] - 0.025, y_position_list[0] + 0.001, 0.001)
            for yy in y_backlash:
                gw.esp.espY.move_to(yy)
                time.sleep(step_wait)

            with NIDAQ() as mynidaq:

                start_time = time.time()

                # start DAQ counting tasks
                num_samples = int(time_per_pixel * trigger_rate)
                # mynidaq.start_external_read_task(trigger_rate, num_samples)

                # start Swabian external trigger for counting
                gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(int(trigger_rate)))

                for n in range(pos_len):

                    # MOVE Y
                    gw.esp.espY.move_to(y_position_list[n])
                    time.sleep(step_wait)
                    start_linescan_time = time.time()

                    # compensate for X backlash
                    gw.esp.espX.move_to(x_position_list[0] - 0.025)
                    time.sleep(step_wait)
                    x_backlash = np.arange(x_position_list[0] - 0.025, x_position_list[0] + 0.001, 0.001)
                    for xx in x_backlash:
                        gw.esp.espX.move_to(xx)
                        time.sleep(step_wait)

                    for nn in range(pos_len):

                        ## MOVE X
                        # print('moving to ', y_position_list[nn])
                        gw.esp.espX.move_to(x_position_list[nn])
                        # gw.esp.espY.wait()
                        time.sleep(step_wait)

                        ## COUNTING
                        reading_period = 1 / trigger_rate
                        counts[n][nn] = np.mean(obtain(mynidaq.internal_read_task(trigger_rate, num_samples))) / (reading_period) #

                        ## PUSH DATA
                        # print("start push")
                        NewportXYScanningData.push(
                            {
                                "params": {
                                    "datasetname": datasetname,
                                    "trigger_rate": trigger_rate,
                                    "device": device,
                                    "x_init_position": x_init_position,
                                    "y_init_position": y_init_position,
                                    "position_steps": position_steps,
                                    "step_length": step_length,
                                    "time_per_pixel": time_per_pixel,
                                },
                                "title": "XYScanning",
                                "xlabel": "X position",
                                "ylabel": "Y position",
                                "datasets": {
                                    "x_position": x_position_list,
                                    "y_position": y_position_list,
                                    "counts": counts,
                                },
                            }
                        )

            # total time
            print('Snac time: ', time.time() - start_time)

            # go back to starting point while accounting for the backlash
            print('Moving to start location')

            # compensate for X backlash
            gw.esp.espX.move_to(x_init_position - 0.025)
            time.sleep(step_wait)
            x_backlash = np.arange(x_init_position - 0.025, x_init_position + 0.001, 0.001)
            for xx in x_backlash:
                gw.esp.espX.move_to(xx)
                time.sleep(step_wait)

            # compensate for Y backlash
            gw.esp.espY.move_to(y_init_position - 0.025)
            time.sleep(step_wait)
            y_backlash = np.arange(y_init_position - 0.025, y_init_position + 0.001, 0.001)
            for yy in y_backlash:
                gw.esp.espY.move_to(yy)
                time.sleep(step_wait)

            print("Plane-scan finished")
