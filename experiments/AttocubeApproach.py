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

from guis.guiElements_general import flexSave
from drivers.zurich.mfli import MFLI
from drivers.ni.nidaq_final import NIDAQ


class Attocube_Approach_Measurement:

	def attocubeapproach(
		self, 
		datasetname: str, 
		step_wait: float, 
		stage_min_voltage: float, 	# min AO voltage to output to scanner Z from DAQ
		stage_max_voltage: float, 	# max AO voltage to output to scanner Z from DAQ
		stage_voltage_steps: int, 	# number of AO voltages to step
		threshhold: float, 			# compare PID error to this value
		A_init: float):

		'''
		This steps the attoube Z stage via DAQ AO control until
	    the MFLI PID error reaches a set threshhold value, at which point it stops
	    '''

		with InstrumentGateway() as gw, DataSource(datasetname) as AttocubeApproach:

			with NIDAQ() as mynidaq:

				# CREATING ANALOG DAQ TASKS
				self.ao_task = mynidaq.create_task()
				self.ao_task.ao_channels.add_ao_voltage_chan('DEV4/AO1')
				self.ao_task.write(0)

				# DATATSETS
				self.voltages = np.zeros(0)
				self.amplitudes = np.zeros(0)

				# COARSE APPROACH: MOVE TO LOWEST VOLTAGE FROM 0 
				# step up to lowest position from zero before scanning to help hysteresis
				if stage_min_voltage != -0.1:
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
						time.sleep(step_wait)
						amplitude_array = 0
						for i in range(av_num):
							amplitude_read = gw.mfli.AUXOUT_read(0)     
							amplitude_array = amplitude_array + amplitude_read
							time.sleep(.15)
						print('summed amplitude_array:',amplitude_array)
						amplitude = amplitude_array / av_num
						print('Fork amplitude : ', amplitude)

						self.voltages = np.append(self.voltages, voltage)
						self.amplitudes = np.append(self.amplitudes, amplitude)		

				        # SAVE CURRENT DATA TO DATA SERVER     
						AttocubeApproach.push(
							{'params':{
								'datasetname': datasetname, 
								'step_wait': step_wait, 
								'stage_min_voltage': stage_min_voltage, 
								'stage_max_voltage': stage_max_voltage, 
								'stage_voltage_steps': stage_voltage_steps, 
								'threshhold': threshhold, 
								'A_init': A_init
							},
							
							'title': 'Attocube Approach',
							'xlabel': 'Scanner Voltage (V)',
							'ylabel': 'Fork Amplitude (V)',
							'datasets': {'voltage': self.voltages,
										'amplitude': self.amplitudes
							}
						})

						if ((amplitude / A_init < threshhold) or (amplitude / A_init > (2 - threshhold))): #must decide how to threshhold
							self.engaged = True
							self.ao_task.write(-0.1)
							print('engaged, activating PID pullback')
							# gw.mfli.set_PID_state(True)


					else:
						print('either engaged or voltage out of range')

				# Bring scanner to 0 for safety, turn off PID, return offset to 0
				# self.ao_task.write(-0.1)

				# CLOSE DAQ TASKS
				self.ao_task.stop()
				self.ao_task.close()
				self.ao_task = None
				print('Closed DAQ tasks')

				print('Total time taken: ', time.time() - start_time, ' s')
				return

				flexSave(datasetName, 'AttocubeApproach', 'final') # after measurement finishes

