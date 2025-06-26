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
from drivers.thorlabs.KDC101 import KDC
from experiments.NewPulses import Pulses
from pylablib.devices import Thorlabs


class XYScan:

    def scanning(
        self,
        datasetname: str,
        device: str,
        x_init_position: float,
        y_init_position: float,
        position_steps: int,
        step_length: float,
        time_per_pixel: float,
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

        with InstrumentGateway() as gw, DataSource(datasetname) as ScanningData: 

            print('Beginning X location (mm): ', gw.kdc.kdcX.get_position() / 34555)
            print('Beginning Y location (mm): ', gw.kdc.kdcY.get_position() / 34555)

            pos_len = int(2 * position_steps + 1)
            x_init = int(x_init_position * 34555)  # scale: 34.5 encoder units = 1 um
            y_init = int(y_init_position * 34555)

            # position lists
            x_position_list = np.linspace(
                x_init - (position_steps * step_length * 34.5),
                x_init + (position_steps * step_length * 34.5),
                pos_len,
            )
            y_position_list = np.linspace(
                y_init - (position_steps * step_length * 34.5),
                y_init + (position_steps * step_length * 34.5),
                pos_len,
            )
            # print("x_position_list ", x_position_list)
            # print("y_position_list ", y_position_list)

            counts = np.ones((pos_len, pos_len))

            with NIDAQ() as mynidaq:

                ## start Swabian external trigger for counting
                trigger_rate = 20e3
                gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(int(trigger_rate)))

                for n in range(pos_len):

                    # MOVE X
                    gw.kdc.kdcX.move_to(int(x_position_list[n]))
                    gw.kdc.kdcX.wait_move()
                    start_linescan_time = time.time()

                    for nn in range(pos_len):

                        # MOVE Y
                        gw.kdc.kdcY.move_to(int(y_position_list[nn]))
                        gw.kdc.kdcY.wait_move()

                        # COUNTING
                        start_count_time = time.time()
                        reading_period = 1 / trigger_rate
                        num_samples = int(time_per_pixel * trigger_rate)
                        counts[n][nn] = np.mean(obtain(mynidaq.internal_read_task(trigger_rate, num_samples))) / (reading_period)
                        counting_time = time.time() - start_count_time
                        # if n == 0 and nn == 0:
                        # 	print('Swabian trigger_rate (Hz) ', trigger_rate)
                        # 	print('DAQ sampling_rate (Hz) ', sampling_rate)
                        # 	print('counting time (s) ', counting_time)

                        ## PUSH DATA
                        # print("start push")
                        x_axis = x_position_list / 34555
                        y_axis = y_position_list / 34555
                        ScanningData.push(
                            {
                                "params": {
                                    "datasetname": datasetname,
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
                                    "x_position": x_axis,
                                    "y_position": y_axis,
                                    "counts": counts,
                                },
                            }
                        )

                    if n == 0 and nn == pos_len - 1:
                        linescan_time = time.time() - start_linescan_time
                        print('linescan time ', linescan_time)
                        print('total scan time (s) ~ ', linescan_time * pos_len)
                    print("Line-scan finished")

            # instead of jumping back to (x_init, y_init) (which causes backlash offset error), stepping slowly back to inital location with a 1s delay
            # x_pos_rev = np.flip(x_position_list)
            # y_pos_rev = np.flip(y_position_list)

            # for n in range(pos_len):
            #     gw.kdc.kdcX.move_to(x_pos_rev[n])
            #     gw.kdc.kdcX.wait_move()
            #     gw.kdc.kdcY.move_to(y_pos_rev[n])
            #     gw.kdc.kdcY.wait_move()
            #     time.sleep(1)

            # # trying to figure out where the problem comes from
            # print('Ending X location (mm): ', gw.kdc.kdcX.get_position() / 34555)
            # print('Ending Y location (mm): ', gw.kdc.kdcY.get_position() / 34555)

            print("Plane-scan finished")
