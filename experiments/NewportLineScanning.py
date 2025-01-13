"""

Linecut scanning with Newport stages

Sanskriti Chitransh, 2024-June-5


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


class LineScan:

    def scanning(
        self,
        datasetname: str,
        trigger_rate: float,
        device: str,
        x_init_position: float,
        y_init_position: float,
        z_init_position: float,
        scanning_axis: int,
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

        with InstrumentGateway() as gw, DataSource(datasetname) as NewportLineScanningData: 

            pos_len = int(2 * position_steps + 1)
            position_list = np.ones(pos_len)
            counts = np.ones(pos_len)

            with NIDAQ() as mynidaq:

                # COUNTING
                # start DAQ counting tasks
                num_samples = int(time_per_pixel * trigger_rate)
                reading_period = 1 / trigger_rate

                # start Swabian external trigger for counting
                gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(int(trigger_rate)))
                start_time = time.time()
                daq_sampling_rate = 20e6

                if scanning_axis == 1:
                    position_list = np.linspace(
                        z_init_position - (position_steps * step_length),
                        z_init_position + (position_steps * step_length),
                        pos_len,
                    )
                    # MOVE Z
                    for n in range(pos_len):
                        gw.esp.espZ.move_to(position_list[n])
                        time.sleep(step_wait)
                        counts[n]= np.mean(obtain(mynidaq.internal_read_task(trigger_rate, num_samples))) / (reading_period)
                        # PUSH DATA
                        NewportLineScanningData.push(
                            {
                                "params": {},
                                "title": "NewportLineScanning",
                                "xlabel": "Position (mm)",
                                "ylabel": "Counts",
                                "datasets": {
                                    "position": position_list,
                                    "counts": counts,
                                },
                            }
                        ) 

                if scanning_axis == 2:
                    position_list = np.linspace(
                        y_init_position - (position_steps * step_length),
                        y_init_position + (position_steps * step_length),
                        pos_len,
                    )
                    # MOVE Y
                    for n in range(pos_len):
                        gw.esp.espY.move_to(position_list[n])
                        time.sleep(step_wait)
                        counts[n]= np.mean(obtain(mynidaq.internal_read_task(trigger_rate, num_samples))) / (reading_period) 
                        # PUSH DATA
                        NewportLineScanningData.push(
                            {
                                "params": {},
                                "title": "NewportLineScanning",
                                "xlabel": "Position (mm)",
                                "ylabel": "Counts",
                                "datasets": {
                                    "position": position_list,
                                    "counts": counts,
                                },
                            }
                        )
            
                if scanning_axis == 3:
                    position_list = np.linspace(
                        x_init_position - (position_steps * step_length),
                        x_init_position + (position_steps * step_length),
                        pos_len,
                    )
                    # MOVE X
                    for n in range(pos_len):
                        gw.esp.espX.move_to(position_list[n])
                        time.sleep(step_wait)
                        counts[n]= np.mean(obtain(mynidaq.internal_read_task(trigger_rate, num_samples))) / (reading_period) 
                        # PUSH DATA
                        NewportLineScanningData.push(
                            {
                                "params": {},
                                "title": "Line Scan",
                                "xlabel": "Position (mm)",
                                "ylabel": "Counts",
                                "datasets": {
                                    "position": position_list,
                                    "counts": counts,
                                },
                            }
                        )
                    

            # total time
            print('Scan time: ', time.time() - start_time)

            # go back to starting point while accounting for the backlash
            print('Moving to start location')

            if scanning_axis == 1: 
                # compensate for Z backlash
                gw.esp.espZ.move_to(z_init_position - 0.025)
                time.sleep(step_wait)
                z_backlash = np.arange(z_init_position - 0.025, z_init_position + 0.001, 0.001)
                for zz in z_backlash:
                    gw.esp.espZ.move_to(zz)
                    time.sleep(step_wait)

            if scanning_axis == 2: 
                # compensate for Y backlash
                gw.esp.espY.move_to(y_init_position - 0.025)
                time.sleep(step_wait)
                y_backlash = np.arange(y_init_position - 0.025, y_init_position + 0.001, 0.001)
                for yy in y_backlash:
                    gw.esp.espY.move_to(yy)
                    time.sleep(step_wait)

            if scanning_axis == 3: 
                # compensate for X backlash
                gw.esp.espX.move_to(x_init_position - 0.025)
                time.sleep(step_wait)
                x_backlash = np.arange(x_init_position - 0.025, x_init_position + 0.001, 0.001)
                for xx in x_backlash:
                    gw.esp.espX.move_to(xx)
                    time.sleep(step_wait)

            print("Line-scan finished")
