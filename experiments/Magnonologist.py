'''
Sanskriti Chitransh, 20205-August-14

Parent code to run scanning experiments to obtain an X-Y-T1 plot for a given height
	1. Begin with a given NV frequency, pi time at some engagement distance from the YIG surface (<1um)
	2. Measure Bayesian T1
	3. Move sample, measure NV frequency (measure of stray field) and apropriate pi-pulse again for next T1 measurement

Modifications for later:
	1. Automatically approach to find engagement voltage and move sample to required separation
	2. Setup NSpyre control of the magnet to avoid manual magnet adjustment (?)

Questions:
	1. Can I call the PIDOutVSTime from here and monitor in parallel? 
	2. How to get PIDOutVSTime to stop everything anytime it retracts? Can the retraction be recorded as 'saved' dat point?

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
from experiments.BayesianT1 import BayesianT1
from experiments.Rabi import Rabi
from experiments.PulsedODMR import PulsedODMR


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
		NV_freq_L: float,
		pi_time_L: float,
		NV_freq_H: float,
		pi_time_H: float,
		Bayesian_R: int, # 500000, can be hardcoded into the Bayesian T1 experiment as well
        objective_x_init_position: float, 
        objective_y_init_position: float,
        objective_z_init_position: float,  
		)


		with InstrumentGateway() as gw, DataSource(datasetname) as MagnonologyData:

			# ATTOCUBE/DAQ VOLTAGE LISTS
			x_voltage_list = np.linspace(x_init_voltage, x_final_voltage, x_voltage_steps)
			y_voltage_list = np.linspace(y_init_voltage, y_final_voltage, y_voltage_steps)
			counts = np.zeros((x_voltage_steps, y_voltage_steps))
			diff = x_voltage_list[1] - x_voltage_list[0]

			# OBJECTIVE INITAL POSITION (for spatial feedback)
			objective_x_init_position = 5.3646
			objective_y_init_position = 18.7688
			objective_z_init_position = -0.1089

			sleep_time = 0.2

			# start_time
			self.start_time = time.time()

			with NIDAQ() as mynidaq:

				# GET CURRENT AO2, AO3
				current_AO2, current_AO3 = mynidaq.analog_read_task()

				# ANALOG DAQ TASKS
				self.ao_x_task = mynidaq.create_task()
				self.ao_x_task.ao_channels.add_ao_voltage_chan(device + '/' + 'AO2')
				self.ao_y_task = mynidaq.create_task()
				self.ao_y_task.ao_channels.add_ao_voltage_chan(device + '/' + 'AO3')

				# MOVE SLOWLY TO INITIAL VOLTAGE 
				if x_init_voltage != current_AO2:
					steps = int((current_AO2 - x_init_voltage) / 0.005)
					x_list = np.linspace(current_AO2, x_init_voltage, steps)
					for x in x_list:
						self.ao_x_task.write(x)
						time.sleep(sleep_time)

				if y_init_voltage != current_AO3:
					steps = int((current_AO3 - y_init_voltage) / 0.005)
					y_list = np.linspace(current_AO3, y_init_voltage, steps)
					for y in y_list:
						self.ao_y_task.write(y)
						time.sleep(sleep_time)

				# SAMPLE SCAN 
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

						# SWABIAN SETUP
						trigger_rate = 20e3
						gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(int(trigger_rate)))

						# BAYESIAN T1
						BayesianT1(datasetName, rf_power, laser_power, R_Bayesian, 1





						counts[i][j] = 

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

						# # SPATIAL FEEDBACK
						# feedback_time = int(time.time() - self.start_time)

						# if (feedback_time >= 180000):

						# 	print('Objective spatial feedback')

						# 	# Parameters
						# 	trigger_rate = int(20e3)
						# 	time_per_point = 0.003
						# 	begin_feedback = time.time()
						# 	objective_x_init_position, objective_y_init_position, current_max_counts = SpatialFeedback.Feedback(objective_x_init_position, objective_y_init_position, objective_z_init_position)
						# 	print('Current objective_x_position ', objective_x_init_position)
						# 	print('Current objective_y_position ', objective_y_init_position)
						# 	print('Returned max counts ', current_max_counts)
							
						# 	gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(int(trigger_rate)))
						# 	num_samples = int(trigger_rate * time_per_point)
						# 	post_feedback_counts = np.mean(obtain(mynidaq.internal_read_task(int(trigger_rate), num_samples))) / (1 / trigger_rate)
						# 	print('Post feedback counts ', post_feedback_counts)
						# 	gw.swabian.reset()

						# 	gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(int(trigger_rate)))
						# 	feedback_duration = time.time() - begin_feedback
						# 	print('Feedback duration: ', feedback_duration)
						# 	self.start_time = time.time()

					# Bring inner dimension to starting voltage
					print('Reset inner dimension')
					y_list = np.flip(np.arange(y_init_voltage, y_final_voltage + diff, diff))
					for y in y_list:
						self.ao_y_task.write(y)
						time.sleep(sleep_time)


				# # return to initial location
				# print('Moving to zero')
				# if x_final_voltage != -0.1:
				# 	print('Move x to zero')
				# 	x_list = np.flip(np.arange(-0.1, x_final_voltage + 0.001, 0.001))
				# 	for x in x_list:
				# 		# print('Applying x voltage ', x)
				# 		self.ao_x_task.write(x)
				# 		time.sleep(0.2)

				# if y_init_voltage != -0.1:
				# 	y_list = np.flip(np.arange(-0.1, y_init_voltage + 0.001, 0.001))
				# 	for y in y_list:
				# 		self.ao_y_task.write(y)
				# 		time.sleep(0.2)

				# CLOSE ALL TASKS
				self.ao_x_task.stop()
				self.ao_x_task.close()
				self.ao_x_task = None

				self.ao_y_task.stop()
				self.ao_y_task.close()
				self.ao_y_task = None
				
				
			print('Scan finished')
