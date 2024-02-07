'''

XY Scanning through DAQ controlled ANM300 Attocube Scanner

Sanskriti Chitransh, 2023-Oct-25

Note: The ANC300/ANM200 Scanner has DC input control enabled from the inserv file and only the DAQ's analog output drives the scanner. 

'''

import numpy as np

from nspyre import DataSource
from nspyre import InstrumentGateway

from guis.guiElements_general import flexSave
from rpyc.utils.classic import obtain

from drivers.ni.nidaq_final import NIDAQ


class XYScan:

	def scanning(self, datasetname: str, device: str, x_init_voltage: float, x_final_voltage: float, y_init_voltage: float, y_final_voltage: float, x_voltage_steps: int, y_voltage_steps: int, counts_per_pixel: float):

		'''
		
		x_init_voltage: starting DAQ analog output voltage to Scanner X
		x_final_voltage: last DAQ analog output voltage to Scanner X
		y_init_voltage: starting DAQ analog output voltage to Scanner Y
		y_final_voltage: last DAQ analog output voltage to Scanner Y
		x_voltage_steps: steps in X direction
		y_voltage_steps: steps in Y direction
		counts_per_pixel: photon clicks for each (X,Y) location

		'''

		with DataSource(datasetname) as ScanningData:

			# VOLTAGE LISTS
			x_voltage_list = np.linspace(x_init_voltage, x_final_voltage, x_voltage_steps)
			y_voltage_list = np.linspace(y_init_voltage, y_final_voltage, y_voltage_steps)
			counts = np.zeros((y_voltage_steps, x_voltage_steps))

			with NIDAQ() as mynidaq:

				# ANALOG DAQ TASKS
				self.ao_x_task = mynidaq.create_task()
				self.ao_x_task.ao_channels.add_ao_voltage_chan(device + '/' + 'AO2')
				self.ao_y_task = mynidaq.create_task()
				self.ao_y_task.ao_channels.add_ao_voltage_chan(device + '/' + 'AO3')

				for vx in x_voltage_list:

					print('Applying x voltage')
					self.ao_x_task.write(vx)

					for vy in y_voltage_list:

						print('Applying y voltage')
						self.ao_y_task.write(vy)

						print('start reading task')
						mynidaq.start_read_task(counts_per_pixel)

						print('start counting')
						raw_counts = mynidaq.read_samples(counts_per_pixel)
						
						print('raw counts : ', raw_counts) 
						np.append(counts, np.mean(raw_counts))

						print('Push data')
						# PUSH DATA
						ScanningData.push({'params': {'datasetname': datasetname, 'device': device, 'x_init_voltage': x_init_voltage, 'x_final_voltage': x_final_voltage, 'y_init_voltage': y_init_voltage, 'y_final_voltage': y_final_voltage, 'x_voltage_steps': x_voltage_steps, 'y_voltage_steps': y_voltage_steps, 'counts_per_pixel': counts_per_pixel
							},
							'title': 'XYScanning',
							'xlabel': 'X Voltage',
							'ylabel': 'Y voltage',
							'datasets': {'xvoltage': x_voltage_list,
										'yvoltage': y_voltage_list,
										'counts': counts
							}
						})

				print('Scan finished')





