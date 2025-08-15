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
		z_init_voltage: float, 
		z_final_voltage: float,
		z_voltage_steps: int, 		# number of AO voltages to step
		step_wait: float, 		 	
		threshhold: float, 			# compare PID error to this value
		A_init: float
		):

		'''
		
		x_init_voltage: starting DAQ analog output voltage to Scanner X
		x_final_voltage: last DAQ analog output voltage to Scanner X
		y_init_voltage: starting DAQ analog output voltage to Scanner Y
		y_final_voltage: last DAQ analog output voltage to Scanner Y
		x_voltage_steps: steps in X direction
		y_voltage_steps: steps in Y direction
		counts_per_pixel: photon clicks for each (X,Y) location

		'''

		with InstrumentGateway() as gw, DataSource(datasetname) as AttocubeApproachXYScanningData:

			# VOLTAGE LISTS
			x_voltage_list = np.linspace(x_init_voltage, x_final_voltage, x_voltage_steps)
			y_voltage_list = np.linspace(y_init_voltage, y_final_voltage, y_voltage_steps)
			z_voltage_list = np.linspace(z_init_voltage, z_final_voltage, z_voltage_steps)

			engagement_voltage = np.zeros((y_voltage_steps, x_voltage_steps))

			print('x_voltage_list ', x_voltage_list)
			print('y_voltage_list ', y_voltage_list)
			# print('z_voltage_list ', z_voltage_list)

			with NIDAQ() as mynidaq:

				# ANALOG DAQ TASKS
				self.ao_z_task = mynidaq.create_task()
				self.ao_z_task.ao_channels.add_ao_voltage_chan(device + '/' + 'AO1')
				self.ao_x_task = mynidaq.create_task()
				self.ao_x_task.ao_channels.add_ao_voltage_chan(device + '/' + 'AO2')
				self.ao_y_task = mynidaq.create_task()
				self.ao_y_task.ao_channels.add_ao_voltage_chan(device + '/' + 'AO3')

				# GET PID VALUE AND ENGAGE MODULATION
				self.pid_setpoint = gw.mfli.get_PID_setpoint
				print('self.pid_setpoint : ', self.pid_setpoint)
				self.engaged = False
				time.sleep(5)
				av_num = 10
				print('Averaging ', av_num, ' times')
				
				for i in range(x_voltage_steps):

					# GO TO Y LOCATION (X VOLTAGE)
					vx = x_voltage_list[i]
					print('Applying x voltage ', vx)
					self.ao_x_task.write(vx)
					time.sleep(0.03)

					for j in range(y_voltage_steps):

						# GO TO X LOCATION (Y VOLTAGE)
						vy = y_voltage_list[j]
						print('Applying y voltage ', vy)
						self.ao_y_task.write(vy)
						time.sleep(0.03)

						# APPROACH
						print('X location ', vy, ' Y location ', vx)
						print('Starting approach at this XY location')

						# Check fork amplitude
						amplitude_array = 0
						for k in range(av_num):
							amplitude_read = gw.mfli.AUXOUT_read(0)     
							amplitude_array = amplitude_array + amplitude_read
							time.sleep(.1)

						# print('summed amplitude_array:',amplitude_array)
						amplitude = amplitude_array / av_num
						print('Pre-coarse approach Fork amplitude : ', amplitude)

						# Coarse stepping
						if z_init_voltage != 0.0:
							print('Coarse approach')
							coarse_approach_voltage_list = np.linspace(0, z_init_voltage, 11)
							for cav in coarse_approach_voltage_list:
								self.ao_z_task.write(cav) 
								time.sleep(step_wait)

						# Check fork amplitude
						amplitude_array = 0
						for k in range(av_num):
							amplitude_read = gw.mfli.AUXOUT_read(0)     
							amplitude_array = amplitude_array + amplitude_read
							time.sleep(.1)

						# print('summed amplitude_array:',amplitude_array)
						amplitude = amplitude_array / av_num
						print('Pre-fine approach Fork amplitude : ', amplitude)
						A_init = amplitude

						# Fine stepping
						print('Fine approach')
						for voltage in z_voltage_list:

							print(self.engaged)
							if ((not self.engaged) and (voltage >= z_init_voltage) and (voltage <= z_final_voltage)):

								print('DAQ analog output voltage : ', voltage)
								self.ao_z_task.write(voltage)

								time.sleep(step_wait)
								amplitude_array = 0

								for k in range(av_num):
									amplitude_read = gw.mfli.AUXOUT_read(0)     
									amplitude_array = amplitude_array + amplitude_read
									time.sleep(.1)

								# print('summed amplitude_array:',amplitude_array)
								amplitude = amplitude_array / av_num
								print('Fork amplitude : ', amplitude)

								if ((amplitude / A_init < threshhold) or (amplitude / A_init > (2 - threshhold))): # decide threshhold
									engagement_voltage[i][j] = voltage
									print('Engagement voltage ', voltage)
									self.engaged = True
									print('engaged')

							else:
								print('Either engaged or voltage out of range')

						print('Engagement voltage matrix ', engagement_voltage)

						# PUSH DATA
						AttocubeApproachXYScanningData.push(
							{'params': {
								'datasetname': datasetname, 
								'device': device, 
								'x_init_voltage': x_init_voltage, 
								'x_final_voltage': x_final_voltage, 
								'y_init_voltage': y_init_voltage, 
								'y_final_voltage': y_final_voltage, 
								'x_voltage_steps': x_voltage_steps, 
								'y_voltage_steps': y_voltage_steps, 
								'z_init_voltage': z_init_voltage, 
								'z_final_voltage': z_final_voltage, 
								'z_voltage_steps': z_voltage_steps,
								'step_wait': step_wait, 
								'threshhold': threshhold, 
								'A_init': A_init
							},
							'title': 'XYScanning',
							'xlabel': 'X Voltage',
							'ylabel': 'Y voltage',
							'datasets': {'xvoltage': x_voltage_list,
										'yvoltage': y_voltage_list,
										'engagement_voltage': engagement_voltage,
							}
						})

						print('Reset engagement flag')
						self.engaged = False

						print('Move Z to 0 before moving XY')
						self.ao_z_task.write(0)

						print('Wait for tip to stabilize before moving')
						time.sleep(200)

				# return to initial XY location
				self.ao_y_task.write(y_init_voltage)
				self.ao_x_task.write(x_init_voltage)
				
			print('Scan finished')





