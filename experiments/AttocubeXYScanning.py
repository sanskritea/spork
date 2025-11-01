'''

XY Scanning through DAQ controlled ANM300 Attocube Scanner

Sanskriti Chitransh, 2023-Oct-25

Note: The ANC300/ANM200 Scanner has DC input control enabled from the inserv file and only the DAQ's analog output drives the scanner. 

'''

import numpy as np
import time

from nspyre import DataSource
from nspyre import InstrumentGateway

from guis.guiElements_general import flexSave
from rpyc.utils.classic import obtain

from drivers.ni.nidaq_final import NIDAQ
from experiments.NewPulses import Pulses
from experiments.NewportSpatialFeedback import SpatialFeedback


class XYScan:

	def scanning(
		self, 
		datasetname: str, 
		x_init_voltage: float, 
		x_final_voltage: float, 
		x_voltage_steps: int,
		y_init_voltage: float, 
		y_final_voltage: float, 
		y_voltage_steps: int, 
		time_per_pixel: float):


		with InstrumentGateway() as gw, DataSource(datasetname) as AttocubeXYScanningData:

			# VOLTAGE LISTS
			x_voltage_list = np.linspace(x_init_voltage, x_final_voltage, x_voltage_steps)
			y_voltage_list = np.linspace(y_init_voltage, y_final_voltage, y_voltage_steps)
			counts = np.zeros((x_voltage_steps, y_voltage_steps))
			diff = x_voltage_list[1] - x_voltage_list[0]

			print('x_voltage_list ', x_voltage_list)
			print('y_voltage_list ', y_voltage_list)

			objective_x_init_position = 5.3646
			objective_y_init_position = 18.7688
			objective_z_init_position = -0.1089

			sleep_time = 0.1

			# start_time
			self.start_time = time.time()

			with NIDAQ() as mynidaq:

				# current AO2/AO3 voltages
				current_AO2, current_AO3 = mynidaq.analog_read_task()
				print('current_AO2 ', current_AO2)
				print('current_AO3 ', current_AO3)


				# ANALOG DAQ TASKS
				# print('Creating DAQ tasks')
				self.ao_x_task = mynidaq.create_task()
				self.ao_x_task.ao_channels.add_ao_voltage_chan('DEV4/AO2')
				self.ao_y_task = mynidaq.create_task()
				self.ao_y_task.ao_channels.add_ao_voltage_chan('DEV4/AO3')

				# move slowly to initial position (assume scanner is at (-0.1, -0.1)
				# print('x_init_voltage ', x_init_voltage)
				# print('current_AO2 ', current_AO2)
				print('Moving to initial voltage')
				if x_init_voltage != current_AO2:
					# print('Move x to initial voltage')
					steps = np.abs(int((current_AO2 - x_init_voltage) / 0.002))
					x_list = np.linspace(current_AO2, x_init_voltage, steps)
					# print('x_list ', x_list)
					for x in x_list:
						# print('Applying x voltage ', x)
						self.ao_x_task.write(x)
						time.sleep(sleep_time)

				if y_init_voltage != current_AO3:
					steps = np.abs(int((current_AO3 - y_init_voltage) / 0.002))
					y_list = np.linspace(current_AO3, y_init_voltage, steps)
					for y in y_list:
						self.ao_y_task.write(y)
						time.sleep(sleep_time)

				# setup
				# print('Setting up Swabian')
				trigger_rate = 20e3
				num_samples = int(time_per_pixel * trigger_rate)
				gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(int(trigger_rate)))

				for i in range(x_voltage_steps):

					vx = x_voltage_list[i]
					print('Applying x voltage ', vx)
					self.ao_x_task.write(vx)
					time.sleep(sleep_time)

					for j in range(y_voltage_steps):

						vy = y_voltage_list[j]
						print('Applying y voltage ', vy)
						self.ao_y_task.write(vy)
						time.sleep(sleep_time)

						# print('start counting')
						reading_period = 1 / trigger_rate
						counts[i][j] = np.mean(obtain(mynidaq.internal_read_task(trigger_rate, num_samples))) / (reading_period)

						# PUSH DATA
						AttocubeXYScanningData.push(
							{'params': {
								'datasetname': datasetname, 
								'x_init_voltage': x_init_voltage, 
								'x_final_voltage': x_final_voltage, 
								'x_voltage_steps': x_voltage_steps, 
								'y_init_voltage': y_init_voltage, 
								'y_final_voltage': y_final_voltage, 
								'y_voltage_steps': y_voltage_steps, 
								'time_per_pixel': time_per_pixel, 
							},
							'title': 'XYScanning',
							'xlabel': 'X Voltage',
							'ylabel': 'Y voltage',
							'datasets': {'xvoltage': x_voltage_list,
										'yvoltage': y_voltage_list,
										'counts': counts
							}
						})

					# Bring inner dimension to starting voltage
					print('Reset inner dimension')
					y_list = np.flip(np.arange(y_init_voltage, y_final_voltage + diff, diff))
					for y in y_list:
						self.ao_y_task.write(y)
						time.sleep(sleep_time)

				# Bring outer dimension to starting voltage
				print('Reset outer dimension')
				x_list = np.flip(np.arange(x_init_voltage, x_final_voltage + diff, diff))
				for x in x_list:
					self.ao_x_task.write(x)
					time.sleep(sleep_time)

				# Save XY Scan data
				notes = ''
				flexSave(datasetname, notes, 'AttocubeXYScan')

				# CLOSE ALL TASKS
				self.ao_x_task.stop()
				self.ao_x_task.close()
				self.ao_x_task = None

				self.ao_y_task.stop()
				self.ao_y_task.close()
				self.ao_y_task = None
				
				
			print('Scan finished')





