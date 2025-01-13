"""

Spatial feedback with Newport stages

Sanskriti Chitransh, 2024-August-20


"""

import numpy as np
import time as time

from nspyre import InstrumentGateway
from rpyc.utils.classic import obtain
from drivers.ni.nidaq_final import NIDAQ
from experiments.NewPulses import Pulses


class SpatialFeedback():


    def Feedback(x_init_position, y_init_position, z_init_position):

        with InstrumentGateway() as gw:

            # paramters
            position_steps = 10
            step_length = 0.0002
            step_wait = 0.05
            pos_len = int(2 * position_steps + 1)
            trigger_rate = int(20e3)
            reading_period = 1 / trigger_rate
            time_per_pixel = 0.005
            num_samples = int(time_per_pixel * trigger_rate)
            position_list = np.ones(pos_len)
            counts = np.ones((pos_len, pos_len))
            x_max_counts = 0
            y_max_counts = 0
            max_counts = 0

            with NIDAQ() as mynidaq:

                # set laser power
                mynidaq.laser_power_atten(2)
                gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(trigger_rate))

                # position lists
                x_max_counts = x_init_position
                x_position_list = np.linspace(
                    x_init_position - (position_steps * step_length),
                    x_init_position + (position_steps * step_length),
                    pos_len,
                )
                y_max_counts = y_init_position
                y_position_list = np.linspace(
                    y_init_position - (position_steps * step_length),
                    y_init_position + (position_steps * step_length),
                    pos_len,
                )
                # print('x_position_list ', x_position_list)
                # print('y_position_list ', y_position_list)

                print('XY Scan')
                # compensate for Y backlash
                gw.esp.espY.move_to(y_position_list[0] - 0.01)
                time.sleep(step_wait)
                gw.esp.espY.move_to(y_position_list[0] - 0.005)
                time.sleep(step_wait)
                y_backlash = np.arange(y_position_list[0] - 0.005, y_position_list[0] + 0.0002, 0.0002)
                for yy in y_backlash:
                    gw.esp.espY.move_to(yy)
                    time.sleep(step_wait)

                # start scanning
                for n in range(pos_len):

                    # MOVE Y
                    gw.esp.espY.move_to(y_position_list[n])
                    time.sleep(step_wait)

                    # compensate for X backlash
                    gw.esp.espX.move_to(x_position_list[0] - 0.01)
                    time.sleep(step_wait)
                    gw.esp.espX.move_to(x_position_list[0] - 0.005)
                    time.sleep(step_wait)
                    x_backlash = np.arange(x_position_list[0] - 0.005, x_position_list[0] + 0.0002, 0.0002)
                    for xx in x_backlash:
                        gw.esp.espX.move_to(xx)
                        time.sleep(step_wait)

                    for nn in range(pos_len):

                        # MOVE X
                        gw.esp.espX.move_to(x_position_list[nn])
                        time.sleep(step_wait)

                        # COUNT
                        counts[n][nn] = np.mean(obtain(mynidaq.internal_read_task(trigger_rate, num_samples))) / reading_period

                        # FIND MAX COUNT LOCATION
                        if counts[n][nn] > max_counts:
                            max_counts = counts[n][nn]
                            x_max_counts = x_position_list[nn]
                            y_max_counts = y_position_list[n]

                # RESET SWABIAN
                gw.swabian.reset()
                time.sleep(5)

                # MOVE TO LOCATION
                print('Move to max counts location')
                print('X location ', x_max_counts, 'mm')
                print('Y location ', y_max_counts, 'mm')
                
                # compensate for X backlash
                gw.esp.espX.move_to(x_max_counts - 0.01)
                time.sleep(step_wait)
                gw.esp.espX.move_to(x_max_counts - 0.005)
                time.sleep(step_wait)
                gw.esp.espX.move_to(x_max_counts - 0.002)
                time.sleep(step_wait)
                x_backlash = np.arange(x_max_counts - 0.002, x_max_counts, 0.0002)
                for xx in x_backlash:
                    gw.esp.espX.move_to(xx)
                    time.sleep(step_wait)

                # compensate for Y backlash
                gw.esp.espY.move_to(y_max_counts - 0.01)
                time.sleep(step_wait)
                gw.esp.espY.move_to(y_max_counts - 0.005)
                time.sleep(step_wait)
                gw.esp.espY.move_to(y_max_counts - 0.002)
                time.sleep(step_wait)
                y_backlash = np.arange(y_max_counts - 0.002, y_max_counts, 0.0002)
                for yy in y_backlash:
                    gw.esp.espY.move_to(yy)
                    time.sleep(step_wait)

                # Check counts after moving
                gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(trigger_rate))
                current_counts = np.mean(obtain(mynidaq.internal_read_task(trigger_rate, int(trigger_rate * 0.01)))) / (1 / trigger_rate)
                gw.swabian.reset()

                print('Max counts ', max_counts)
                print('Counts after ', current_counts)

                return x_max_counts, y_max_counts, max_counts



    # def XYScanning(x_init_position, y_init_position, initial_counts, position_steps, step_length, step_wait, pos_len, trigger_rate, reading_period, num_samples): 

    #     print('XY Feedback')

    #     with InstrumentGateway() as gw: 

    #         position_list = np.ones(pos_len)
    #         counts = np.ones((pos_len, pos_len))

    #         with NIDAQ() as mynidaq:

    #             # DAQ parameters
    #             gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(trigger_rate))

    #             # position lists
    #             x_max_counts = x_init_position
    #             x_position_list = np.linspace(
    #                 x_init_position - (position_steps * step_length),
    #                 x_init_position + (position_steps * step_length),
    #                 pos_len,
    #             )
    #             y_max_counts = y_init_position
    #             y_position_list = np.linspace(
    #                 y_init_position - (position_steps * step_length),
    #                 y_init_position + (position_steps * step_length),
    #                 pos_len,
    #             )
    #             # print('x_position_list ', x_position_list)
    #             # print('y_position_list ', y_position_list)

    #             # compensate for Y backlash
    #             gw.esp.espY.move_to(y_position_list[0] - 0.025)
    #             time.sleep(step_wait)
    #             y_backlash = np.arange(y_position_list[0] - 0.025, y_position_list[0] + 0.001, 0.001)
    #             for yy in y_backlash:
    #                 gw.esp.espY.move_to(yy)
    #                 time.sleep(step_wait)

    #             # start scanning
    #             for n in range(pos_len):

    #                 # MOVE Y
    #                 gw.esp.espY.move_to(y_position_list[n])
    #                 time.sleep(step_wait)

    #                 # compensate for X backlash
    #                 gw.esp.espX.move_to(x_position_list[0] - 0.025)
    #                 time.sleep(step_wait)
    #                 x_backlash = np.arange(x_position_list[0] - 0.025, x_position_list[0] + 0.001, 0.001)
    #                 for xx in x_backlash:
    #                     gw.esp.espX.move_to(xx)
    #                     time.sleep(step_wait)

    #                 for nn in range(pos_len):

    #                     ## MOVE X
    #                     # print('moving to ', y_position_list[nn])
    #                     gw.esp.espX.move_to(x_position_list[nn])
    #                     # gw.esp.espY.wait()
    #                     time.sleep(step_wait)

    #                     # COUNT
    #                     counts[n][nn] = np.mean(obtain(mynidaq.internal_read_task(trigger_rate, num_samples))) / reading_period

    #             # RESET SWABIAN
    #             gw.swabian.reset()
                
    #             # FIND MAX COUNT LOCATION, MOVE
    #             max_counts = np.max(counts)
    #             XY_flag = 100 * np.abs(max_counts - initial_counts) / initial_counts 
    #             x_loc, y_loc = np.unravel_index(np.argmax(counts), np.shape(counts))
    #             x_max_counts, y_max_counts = x_position_list[x_loc], y_position_list[y_loc]

    #             # print('Moving after XY scan')
    #             # # compensate for Y backlash
    #             # gw.esp.espY.move_to(y_max_counts - 0.01)
    #             # time.sleep(step_wait)
    #             # gw.esp.espY.move_to(y_max_counts - 0.01 + 0.005)
    #             # time.sleep(step_wait)
    #             # # final move
    #             # y_backlash = np.arange(y_max_counts - 0.005, y_max_counts + 0.0002, 0.0002)
    #             # for yy in y_backlash:
    #             #     gw.esp.espY.move_to(yy)
    #             #     time.sleep(step_wait)

    #             # # compensate for X backlash
    #             # gw.esp.espX.move_to(x_max_counts - 0.01)
    #             # time.sleep(step_wait)
    #             # gw.esp.espX.move_to(x_max_counts - 0.01 + 0.005)
    #             # time.sleep(step_wait)
    #             # # final move
    #             # x_backlash = np.arange(x_max_counts - 0.005, x_max_counts + 0.0002, 0.0002)
    #             # for xx in x_backlash:
    #             #     gw.esp.espX.move_to(xx)
    #             #     time.sleep(step_wait)

    #             # print('XY_flag ', XY_flag)
    #             return XY_flag, x_max_counts, y_max_counts


    # def ZScanning(z_init_position, initial_counts, position_steps, step_length, step_wait, pos_len, trigger_rate, reading_period, num_samples):

    #     print('Z Feedback')

    #     with InstrumentGateway() as gw: 

    #         position_list = np.ones(pos_len)
    #         counts = np.ones(pos_len)

    #         with NIDAQ() as mynidaq:

    #             # DAQ parameters
    #             gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(trigger_rate))

    #             # position list
    #             z_max_counts = z_init_position
    #             z_position_list = np.linspace(
    #                     z_init_position - (position_steps * step_length),
    #                     z_init_position + (position_steps * step_length),
    #                     pos_len,
    #                 )

    #             # MOVE Z
    #             for n in range(pos_len):
    #                 gw.esp.espZ.move_to(z_position_list[n])
    #                 time.sleep(step_wait)
    #                 counts[n]= np.mean(obtain(mynidaq.internal_read_task(trigger_rate, num_samples))) / reading_period

    #             # RESET SWABIAN
    #             gw.swabian.reset()

    #             # FIND MAX COUNT LOCATION, MOVE
    #             max_counts = np.max(counts) 
    #             Z_flag = 100 * np.abs(max_counts - initial_counts) / initial_counts
    #             z_loc = np.unravel_index(np.argmax(counts), np.shape(counts))[0]
    #             z_max_counts = z_position_list[z_loc]

    #             print('Moving after Z scan')
    #             # compensate for backlash
    #             gw.esp.espZ.move_to(z_max_counts - 0.01)
    #             time.sleep(step_wait)
    #             gw.esp.espZ.move_to(z_max_counts - 0.01 + 0.005)
    #             time.sleep(step_wait)
    #             # final move
    #             z_backlash = np.arange(z_max_counts - 0.005, z_max_counts + 0.0002, 0.0002)
    #             for zz in z_backlash:
    #                 gw.esp.espZ.move_to(zz)
    #                 time.sleep(step_wait)

    #             print('Z_flag ', Z_flag)
    #             return Z_flag, z_max_counts


               

            


            

        