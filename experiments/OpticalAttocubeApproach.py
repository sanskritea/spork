'''
NSpyre v0.6.1 compatible application to control Attocube Scanner approach towards the AFM tip with PID feedback

Sanskriti Chitransh, 2023-Oct-24

Adapted from Jonathan Marcks, Michael Solomon

'''

import time
from itertools import count

import numpy as np

from nspyre import DataSource
from nspyre import InstrumentGateway

from rpyc.utils.classic import obtain
from experiments.NewPulses import Pulses
from guis.guiElements_general import flexSave
from drivers.zurich.mfli import MFLI
from drivers.ni.nidaq_final import NIDAQ

from experiments.SingleBayesianT1 import Single_Bayesian_T1_Meas
from experiments.Rabi import Rabi_Measurement
from experiments.PulsedODMR import Pulsed_ODMR_Measurement


class Optical_Attocube_Approach_Measurement:

	def opticalattocubeapproach(
		self, 
		datasetname: str,  
		step_wait: float, 
		stage_min_voltage: float, 	# min AO voltage to output to scanner Z from DAQ
		stage_max_voltage: float, 	# max AO voltage to output to scanner Z from DAQ
		stage_voltage_steps: int, 	# number of AO voltages to step
		spin_measurement_steps: int, # take spin data but not at every step moved
		threshhold: float, 			# compare PID error to this value
		A_init: float,
		rf_power: float,
		NV_freq: float,
		pi_time: float,
		):

		'''
		This steps the attoube Z stage via DAQ AO control until
	    the MFLI PID error reaches a set threshhold value, at which point it stops
	    '''

		with InstrumentGateway() as gw, DataSource(datasetname) as OpticalAttocubeApproach:

			# for laser counting
			notes = ''
			trigger_rate = 20e3
			num_samples = int(trigger_rate * step_wait)
			trigger_period = 1 / trigger_rate
			laser_power = 0
			

			with NIDAQ() as mynidaq:

				# CREATING ANALOG DAQ TASKS
				self.ao_task = mynidaq.create_task()
				self.ao_task.ao_channels.add_ao_voltage_chan('Dev4/AO1')
				self.ao_task.write(0)

				# GET PID VALUE AND ENGAGE MODULATION
				self.pid_setpoint = gw.mfli.get_PID_setpoint
				print('self.pid_setpoint : ', self.pid_setpoint)
				self.engaged = False
				time.sleep(5)

				# DATATSETS
				self.voltages = np.zeros(0)
				self.amplitudes = np.zeros(0)
				self.counts = np.zeros(0)
				self.pi_times = np.zeros(0)
				self.freqs = np.zeros(0)
				self.t1s = np.zeros(0)	

				# COARSE APPROACH: MOVE TO LOWEST VOLTAGE FROM 0 
				# step up to lowest position from zero before scanning to help hysteresis
				if stage_min_voltage != 0.0:
					print('Coarse stepping')
					coarse_approach_voltage_list = np.linspace(0, stage_min_voltage, 11)
					for cav in coarse_approach_voltage_list:
						self.ao_task.write(cav) 
						time.sleep(step_wait)

				av_num = 10
				print('Averaging ', av_num, ' times')

				start_time = time.time()

				# FINE APPROACH: MOVE SCANNER THROUGH SPECIFIED VOLTAGES
				for voltage in np.linspace(stage_min_voltage, stage_max_voltage, stage_voltage_steps):

					print(self.engaged)
					if ((not self.engaged) and (voltage >= stage_min_voltage) and (voltage <= stage_max_voltage)):

						print('DAQ analog output voltage : ', voltage)
						self.ao_task.write(voltage)
						self.voltages = np.append(self.voltages, voltage)

						time.sleep(step_wait)
						amplitude_array = 0
						counts_array = 0
						for i in range(av_num):
							amplitude_read = gw.mfli.AUXOUT_read(0)     
							amplitude_array = amplitude_array + amplitude_read
							time.sleep(.15)
						print('summed amplitude_array:', amplitude_array)
						amplitude = amplitude_array / av_num
						print('Fork amplitude : ', amplitude)
						self.amplitudes = np.append(self.amplitudes, amplitude)

						gw.swabian.runSequenceInfinitely(Pulses(gw).counting_trigger(int(trigger_rate)))
						raw_counts = np.mean(obtain(mynidaq.internal_read_task(int(trigger_rate), num_samples)) / trigger_period)
						print('counts ', raw_counts)
						gw.swabian.reset()
						self.counts	= np.append(self.counts, raw_counts)	

						# # Re-run Rabi at current frequency and get pi time
						# print('Updating NV freq and pi-time')
						# print('Check current pi time')
						# intermediate_pi_time = Rabi_Measurement().Rabi('RabiData', 10000, 10, rf_power, laser_power, NV_freq, 0, 1500, 31)
						# print('Intermediate pi_time (ns) ', intermediate_pi_time)

						# # Use updated pi time to find new NV resonance peak
						# print('Check current NV freq')
						# intermediate_freq = Pulsed_ODMR_Measurement().PulsedODMR('PulsedODMRdata', 10000, 5, rf_power, laser_power, (NV_freq - 0.010e9), (NV_freq + 0.010e9), 41, int(intermediate_pi_time)) 
						# print('Intermediate NV_freq (GHz) ', intermediate_freq)
						# intermediate_freq = intermediate_freq*1e9
						
						# # Find updated pi time for this peak
						# print('Confirm pi time for this frequency')
						# pi_time = Rabi_Measurement().Rabi('RabiData', 10000, 10, rf_power, laser_power, intermediate_freq, 0, 1500, 31)
						# self.pi_times = np.append(self.pi_times, pi_time)
						# print('updated_pi_time (ns)', pi_time)

						# # Final peak location
						# print('Confirm NV freq with final pi pulse')
						# NV_freq = Pulsed_ODMR_Measurement().PulsedODMR('PulsedODMRdata', 10000, 5, rf_power, laser_power, (intermediate_freq - 0.010e9), (intermediate_freq + 0.010e9), 41, int(pi_time))
						# self.freqs = np.append(self.freqs, NV_freq)
						# print('updated_freq (GHz) ', NV_freq)
						# NV_freq = NV_freq*1e9

						# # BAYESIAN T1
						# print('Measuring T1')
						# gamma_local = Single_Bayesian_T1_Meas().SingleBayesianT1(
						# 	'SingleBayesianT1Data', 
						# 	100000, # Bayesian R 
						# 	rf_power, 
						# 	laser_power, 
						# 	NV_freq, 
						# 	pi_time, 
						# 	5, # 3 Bayesian iterations should be enough for testing, increase later if needed
						# 	'Scanning_Z=' + str(voltage)
						# 	)
						# t1_local = np.round(1000 / gamma_local, 2) # T1 in us since gamma is in ms^-1
						# print('T1_local (μs) ', t1_local)
						# self.t1s = np.append(self.t1s, t1_local)

				        # SAVE CURRENT DATA TO DATA SERVER     
						OpticalAttocubeApproach.push(
								{'params':{
								'datasetname': datasetname, 
								'step_wait': step_wait, 
								'stage_min_voltage': stage_min_voltage, 
								'stage_max_voltage': stage_max_voltage, 
								'stage_voltage_steps': stage_voltage_steps, 
								'threshhold': threshhold, 
								'A_init': A_init,
							},

							'title': 'Optical Attocube Approach',
							'xlabel': 'Scanner Voltage (V)',
							'ylabel': 'Readouts',
							'datasets': {'voltage': self.voltages,
										'amplitude': self.amplitudes,
										'counts': self.counts,
										'pi_times': self.pi_times,
										'freqs': self.freqs,
										't1s': self.t1s
							}
							}
						)

						if ((amplitude / A_init < threshhold) or (amplitude / A_init > (2 - threshhold))): #must decide how to threshhold
							self.engaged = True
							print('engaged')
							# Bring scanner to 0 for safety
							self.ao_task.write(0)

						else:
							print('either engaged or voltage out of range')
							
						flexSave(datasetname, notes, 'OpticalAttocubeApproach') # intermediate save

				notes = 'Final'
				flexSave(datasetname, notes, 'OpticalAttocubeApproach') # final save

				# CLOSE DAQ TASKS
				self.ao_task.stop()
				self.ao_task.close()
				self.ao_task = None
				print('Closed DAQ tasks')

				print('Total time taken: ', time.time() - start_time, ' s')
				return


