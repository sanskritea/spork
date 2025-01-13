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


class XYScan:

	def scanning(
		self, 
		datasetname: str, 
		device: str, 
		x_init_voltage: float, 
		x_final_voltage: float, 
		y_init_voltage: float, 
		y_final_voltage: float, 
		x_voltage_steps: int, 
		y_voltage_steps: int, 
		time_per_pixel: float):

		'''
		
		x_init_voltage: starting DAQ analog output voltage to Scanner X
		x_final_voltage: last DAQ analog output voltage to Scanner X
		y_init_voltage: starting DAQ analog output voltage to Scanner Y
		y_final_voltage: last DAQ analog output voltage to Scanner Y
		x_voltage_steps: steps in X direction
		y_voltage_steps: steps in Y direction
		counts_per_pixel: photon clicks for each (X,Y) location

		'''

		with InstrumentGateway() as gw, DataSource(datasetname) as AttocubeXYScanningData:

			# VOLTAGE LISTS
			x_voltage_list = np.linspace(x_init_voltage, x_final_voltage, x_voltage_steps)
			y_voltage_list = np.linspace(y_init_voltage, y_final_voltage, y_voltage_steps)
			counts = np.zeros((y_voltage_steps, x_voltage_steps))

			print('x_voltage_list ', x_voltage_list)
			print('y_voltage_list ', y_voltage_list)

			with NIDAQ() as mynidaq:

				# setup
				trigger_rate = 20e3
				num_samples = int(time_per_pixel * trigger_rate)
				gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(int(trigger_rate)))

				# ANALOG DAQ TASKS
				self.ao_x_task = mynidaq.create_task()
				self.ao_x_task.ao_channels.add_ao_voltage_chan(device + '/' + 'AO2')
				self.ao_y_task = mynidaq.create_task()
				self.ao_y_task.ao_channels.add_ao_voltage_chan(device + '/' + 'AO3')

				for i in range(x_voltage_steps):

					vx = x_voltage_list[i]
					# print('Applying x voltage')
					self.ao_x_task.write(vx)
					time.sleep(0.03)

					for j in range(y_voltage_steps):

						vy = y_voltage_list[j]
						# print('Applying y voltage')
						self.ao_y_task.write(vy)
						time.sleep(0.03)

						# print('start counting')
						reading_period = 1 / trigger_rate
						counts[i][j] = np.mean(obtain(mynidaq.internal_read_task(trigger_rate, num_samples))) / (reading_period)

						# print('Push data')
						# PUSH DATA
						AttocubeXYScanningData.push(
							{'params': {
								'datasetname': datasetname, 
								'device': device, 
								'x_init_voltage': x_init_voltage, 
								'x_final_voltage': x_final_voltage, 
								'y_init_voltage': y_init_voltage, 
								'y_final_voltage': y_final_voltage, 
								'x_voltage_steps': x_voltage_steps, 
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

				print('check counts')
				print('First location and counts: X voltage ', x_voltage_list[0], ' Y voltage ', y_voltage_list[0], 'counts ',  counts[0][0])
				print('Last location and counts: X voltage ', x_voltage_list[x_voltage_steps - 1], ' Y voltage ', y_voltage_list[y_voltage_steps - 1], 'counts ',  counts[x_voltage_steps - 1][y_voltage_steps - 1])

				# return to initial location
				self.ao_y_task.write(y_init_voltage)
				self.ao_x_task.write(x_init_voltage)
				
			print('Scan finished')





