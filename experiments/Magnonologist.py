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
	2. How to get PIDOutVSTime to stop everything anytime it retracts? Can the retraction be recorded as 'saved' data point?

'''

import numpy as np
import time

from nspyre import DataSource
from nspyre import InstrumentGateway
from datetime import datetime

from guis.guiElements_general import flexSave
from rpyc.utils.classic import obtain

from drivers.ni.nidaq_final import NIDAQ
from guis.guiElements_general import flexSave

from experiments.NewPulses import Pulses
# from experiments.NewportSpatialFeedback import SpatialFeedback
from experiments.SingleBayesianT1 import Single_Bayesian_T1_Meas
from experiments.Rabi import Rabi_Measurement
from experiments.PulsedODMR import Pulsed_ODMR_Measurement


class Magnonologist_Measurement:

	def MagnonologistScan(
		self, 
		datasetname: str, 
		x_init_voltage: float, 
		x_final_voltage: float, 
		y_init_voltage: float, 
		y_final_voltage: float, 
		x_measurement_points: int,  # Number of points to MEASURE at
		y_measurement_points: int,  # Number of points to MEASURE at
		step_size_volts: float,     # Step size for movement (e.g., 0.002 V)
		rf_power: float,
        # laser_power: float,
		NV_freq: float,
		pi_time: float,
		# Bayesian_R: int, # 500000, can be hardcoded into the Bayesian T1 experiment as well
        objective_x_init_position: float, 
        objective_y_init_position: float,
        objective_z_init_position: float,  
		):

		with InstrumentGateway() as gw, DataSource(datasetname) as MagnonologistData:

			print('NV_freq ', NV_freq)
			print('pi_time ', pi_time)
			laser_power = 0 # hardcoded for now since not being used, can add as parameter later if needed

			# # ATTOCUBE/DAQ VOLTAGE LISTS
			# x_voltage_list = np.linspace(x_init_voltage, x_final_voltage, x_voltage_steps)
			# y_voltage_list = np.linspace(y_init_voltage, y_final_voltage, y_voltage_steps)
			# counts = np.zeros((y_voltage_steps, x_voltage_steps))
			# diff = x_voltage_list[1] - x_voltage_list[0]
			# print('x_voltage_list ', x_voltage_list)
			# print('y_voltage_list ', y_voltage_list)

			# MEASUREMENT POINTS - where we actually take data
			x_measurement_voltages = np.linspace(x_init_voltage, x_final_voltage, x_measurement_points)
			y_measurement_voltages = np.linspace(y_init_voltage, y_final_voltage, y_measurement_points)

			print('x_measurement_voltages ', x_measurement_voltages)
			print('y_measurement_voltages ', y_measurement_voltages)
			print('step_size_volts ', step_size_volts)

			# # Measure T1 at only these points
			# # measure_x_voltage_list = np.linspace(x_init_voltage, x_final_voltage, x_voltage_steps)
			# # meau_voltage_list = np.linspace(y_init_voltage, y_final_voltage, y_voltage_steps)

			# STORAGE 
			t1_table = np.zeros((x_measurement_points, y_measurement_points))
			current_counts_table = np.zeros((x_measurement_points, y_measurement_points))
			current_time_table = np.zeros((x_measurement_points, y_measurement_points), dtype=object)
			current_pi_time_table = np.zeros((x_measurement_points, y_measurement_points), dtype=object)
			current_podmr_freq_table = np.zeros((x_measurement_points, y_measurement_points), dtype=object)
			# start_time
			self.start_time = time.time()

			# other constants
			sleep_time = 0.2

			with NIDAQ() as mynidaq:

				# current AO2/AO3 voltages
				current_AO2, current_AO3 = mynidaq.analog_read_task()
				print('current_AO2 ', current_AO2)
				print('current_AO3 ', current_AO3)

				# ANALOG DAQ TASKS
				# print('Creating DAQ tasks')
				self.ao_x_task = mynidaq.create_task()
				self.ao_x_task.ao_channels.add_ao_voltage_chan('DEV4/AO2') # INNER LOOP
				self.ao_y_task = mynidaq.create_task()
				self.ao_y_task.ao_channels.add_ao_voltage_chan('DEV4/AO3') # OUTER LOOP

				# # MOVE SLOWLY TO INITIAL VOLTAGE 
				# print('Move to initial voltage')
				# # move both x and y together from current position to starting position for scan, note that the number of steps for x and y might be different
				# x_steps = np.abs(int((current_AO2 - x_init_voltage) / 0.002))
				# y_steps = np.abs(int((current_AO3 - y_init_voltage) / 0.002))
				# max_steps = max(x_steps, y_steps)		
				# x_list = np.linspace(current_AO2, x_init_voltage, max_steps)
				# y_list = np.linspace(current_AO3, y_init_voltage, max_steps)	
				# for x, y in zip(x_list, y_list):
				# 	self.ao_x_task.write(x)
				# 	self.ao_y_task.write(y)
				# 	time.sleep(sleep_time)

				# MOVE SLOWLY TO INITIAL VOLTAGE 
				print('Move to initial voltage')
				self._move_to_voltage(
					current_AO2, current_AO3, 
					x_init_voltage, y_init_voltage, 
					step_size_volts, sleep_time
				)

				# SAMPLE SCAN 
				for i in range(y_measurement_points):

					vy = y_measurement_voltages[i]
					print(f'\n=== Y MEASUREMENT POINT {i+1}/{y_measurement_points} ===')
					print('Target y voltage ', vy)
					
					# Move to this y position in steps
					current_AO2, current_AO3 = mynidaq.analog_read_task()
					current_y = current_AO3 if i > 0 else y_init_voltage
					self._move_axis(self.ao_y_task, current_y, vy, step_size_volts, sleep_time, 'Y')

					# vy = y_voltage_list[i]
					# print('Applying y voltage ', vy)
					# self.ao_y_task.write(vy)
					# time.sleep(sleep_time)

					for j in range(x_measurement_points):
						vx = x_measurement_voltages[j]
						print(f'\n--- X MEASUREMENT POINT {j+1}/{x_measurement_points} ---')
						print('Target x voltage ', vx)
						
						# Move to this x position in steps
						current_AO2, current_AO3 = mynidaq.analog_read_task()
						current_x = current_AO2 if j > 0 else x_init_voltage
						self._move_axis(self.ao_x_task, current_x, vx, step_size_volts, sleep_time, 'X')

						# save current date and time
						current_time_table[i][j] = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
						gw.swabian.reset()

					# for j in range(x_voltage_steps):
					# 	vx = x_voltage_list[j]
					# 	print('Applying x voltage ', vx)
					# 	self.ao_x_task.write(vx)
					# 	time.sleep(sleep_time)

						# === TAKE MEASUREMENTS AT THIS POINT ===

						# SWABIAN SETUP
						feedback_trigger_rate = 20e3
						gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(int(feedback_trigger_rate)))
						current_counts = np.mean(obtain(mynidaq.internal_read_task(feedback_trigger_rate, int(feedback_trigger_rate * 0.01)))) / (1 / feedback_trigger_rate)
						current_counts_table[i][j] = current_counts
						gw.swabian.reset()

						

						# Re-run Rabi at current frequency and get pi time
						print('Updating NV freq and pi-time')
						print('Check current pi time')
						intermediate_pi_time = Rabi_Measurement().Rabi('RabiData', 10000, 20, rf_power, laser_power, NV_freq, 0, 1500, 31)
						print('Intermediate pi_time (ns) ', intermediate_pi_time)

						if intermediate_pi_time < 100:
							print('Fitting error, using previous pi time')
							intermediate_pi_time = pi_time

						# Use updated pi time to find new NV resonance peak
						print('Check current NV freq')
						intermediate_freq = Pulsed_ODMR_Measurement().PulsedODMR('PulsedODMRdata', 10000, 5, rf_power, laser_power, (NV_freq - 0.010e9), (NV_freq + 0.010e9), 41, int(intermediate_pi_time)) 
						print('Intermediate NV_freq (GHz) ', intermediate_freq)
						intermediate_freq = intermediate_freq*1e9
						
						# Find updated pi time for this peak
						print('Confirm pi time for this frequency')
						pi_time = Rabi_Measurement().Rabi('RabiData', 10000, 20, rf_power, laser_power, intermediate_freq, 0, 1500, 31)
						current_pi_time_table[i][j] = pi_time
						print('updated_pi_time (ns)', pi_time)

						# Final peak location
						print('Confirm NV freq with final pi pulse')
						NV_freq = Pulsed_ODMR_Measurement().PulsedODMR('PulsedODMRdata', 10000, 5, rf_power, laser_power, (intermediate_freq - 0.010e9), (intermediate_freq + 0.010e9), 41, int(pi_time))
						current_podmr_freq_table[i][j] = NV_freq
						print('updated_freq (GHz) ', NV_freq)
						NV_freq = NV_freq*1e9

						# BAYESIAN T1
						print('Measuring T1')
						gamma_local = Single_Bayesian_T1_Meas().SingleBayesianT1(
							'SingleBayesianT1Data', 
							200000, # Bayesian R 
							rf_power, 
							laser_power, 
							NV_freq, 
							pi_time, 
							10, # 3 Bayesian iterations should be enough for testing, increase later if needed
							'Scanning_X=' + str(vx) + '_Y=' + str(vy)
							)
						t1_local = np.round(1000 / gamma_local, 2) # T1 in us since gamma is in ms^-1
						print('T1_local (μs) ', t1_local)
						t1_table[i][j] = t1_local

						# PUSH DATA
						MagnonologistData.push(
							{'params': {
								'datasetname': datasetname,  
								'x_init_voltage': x_init_voltage, 
								'x_final_voltage': x_final_voltage, 
								'y_init_voltage': y_init_voltage, 
								'y_final_voltage': y_final_voltage,
								'x_measurement_points': x_measurement_points, 
								'y_measurement_points': y_measurement_points,
								'step_size_volts': step_size_volts, 
								# 'x_voltage_steps': x_voltage_steps, 
								# 'y_voltage_steps': y_voltage_steps,
								'rf_power': rf_power,
								# 'laser_power': laser_power,
								'NV_freq': NV_freq,
								'pi_time': pi_time,
								# 'Bayesian_R': Bayesian_R,
								'objective_x_init_position': objective_x_init_position,
								'objective_y_init_position': objective_y_init_position,
								'objective_z_init_position': objective_z_init_position  
							},
							'title': 'T1 Scan',
							'xlabel': 'X Voltage',
							'ylabel': 'Y voltage',
							'datasets': {'xvoltage': x_measurement_voltages,
										'yvoltage': y_measurement_voltages,
										# 'xvoltage': x_voltage_list,
										# 'yvoltage': y_voltage_list,
										't1': t1_table,
										'current_counts': current_counts_table,
										'current_time': current_time_table,
										'current_pi_time': current_pi_time_table,
										'current_podmr_freq': current_podmr_freq_table
							}
						})

						# Reset Swabian before next point
						gw.swabian.reset()

					# save each linescan
					notes = ''
					flexSave(datasetname, notes, 'T1AttocubeScan') # Intermediate save after each linescan

					# After completing X scan at this Y position, return X to initial
					print('\nResetting X to initial position')
					current_AO2, current_AO3 = mynidaq.analog_read_task()
					current_x = current_AO2
					self._move_axis(self.ao_x_task, current_x, x_init_voltage, step_size_volts, sleep_time, 'X')

				# final save 
				notes = 'Final'
				flexSave(datasetname, notes, 'T1AttocubeScan') # Final save after total scan

				# After completing scan, move both X and Y to (3.9V, 3.9V)
				print('\nResetting both axes to position (3.9V, 3.9V)')
				current_x, current_y = mynidaq.analog_read_task()
				self._move_to_voltage(current_x, current_y, 3.9, 3.9, step_size_volts, sleep_time)

				# CLOSE ALL TASKS
				self.ao_x_task.stop()
				self.ao_x_task.close()
				self.ao_x_task = None

				self.ao_y_task.stop()
				self.ao_y_task.close()
				self.ao_y_task = None
				
				
			print('Scan finished')

	
	def _move_axis(self, task, current_voltage, target_voltage, step_size, sleep_time, axis_name):
		"""
		Move a single axis from current to target voltage in fixed steps.
		
		Args:
			task: DAQ task for this axis
			current_voltage: Starting voltage
			target_voltage: Ending voltage
			step_size: Step size in volts
			sleep_time: Time to wait after each step
			axis_name: 'X' or 'Y' for logging
		"""
		if abs(target_voltage - current_voltage) < step_size/10:
			# Already at target
			return
			
		# Generate steps
		num_steps = int(abs(target_voltage - current_voltage) / step_size)
		voltage_steps = np.linspace(current_voltage, target_voltage, num_steps + 1)
		
		print(f'  Moving {axis_name} from {current_voltage:.4f}V to {target_voltage:.4f}V in {num_steps} steps')
		
		for v in voltage_steps[1:]:  # Skip first point (current position)
			task.write(v)
			time.sleep(sleep_time)
		
		# Ensure we're exactly at target
		task.write(target_voltage)


	def _move_to_voltage(self, x_current, y_current, x_target, y_target, step_size, sleep_time):
		"""
		Move both X and Y axes from current to target positions in synchronized steps.
		
		Args:
			x_current, y_current: Current positions
			x_target, y_target: Target positions
			step_size: Step size in volts
			sleep_time: Time to wait after each step
		"""
		# Calculate number of steps needed for each axis
		x_steps = int(abs(x_target - x_current) / step_size)
		y_steps = int(abs(y_target - y_current) / step_size)
		max_steps = max(x_steps, y_steps)
		
		if max_steps == 0:
			return  # Already at target
		
		# Generate synchronized step lists
		x_list = np.linspace(x_current, x_target, max_steps + 1)
		y_list = np.linspace(y_current, y_target, max_steps + 1)
		
		print(f'Moving to initial position in {max_steps} steps')
		
		for x, y in zip(x_list[1:], y_list[1:]):  # Skip first point (current position)
			self.ao_x_task.write(x)
			self.ao_y_task.write(y)
			time.sleep(sleep_time)
