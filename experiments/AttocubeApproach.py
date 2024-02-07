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

	def attocubeapproach(self, datasetname: str, device: str, AOchannel: str, 
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
				self.ao_task.ao_channels.add_ao_voltage_chan(device + '/' + AOchannel)
				self.ao_task.write(stage_min_voltage) # go back to lowest position before scanning to help hysteresis

				# GET PID VALUE AND ENGAGE MODULATION
				self.pid_setpoint = gw.mfli.get_PID_setpoint
				print('self.pid_setpoint : ', self.pid_setpoint)
				self.engaged = False
				time.sleep(5)

				# DATATSETS
				self.voltages = np.zeros(0)
				self.amplitudes = np.zeros(0)

				# START SCANNER Z MOTION
				for voltage in np.linspace(stage_min_voltage, stage_max_voltage, stage_voltage_steps):

					print(self.engaged)
					if ((not self.engaged) and (voltage >= stage_min_voltage) and (voltage <= stage_max_voltage)):

						print('DAQ analog output voltage : ', voltage)
						self.ao_task.write(voltage)
						time.sleep(step_wait)
						amplitude_array = 0
						for i in range(100):
							amplitude_read = gw.mfli.AUXOUT_read(0)     
							amplitude_array = amplitude_array + amplitude_read
							time.sleep(.01)
						print('summed amplitude_array:',amplitude_array)
						amplitude = amplitude_array / 100
						print('Fork amplitude : ', amplitude)

						self.voltages = np.append(self.voltages, voltage)
						self.amplitudes = np.append(self.amplitudes, amplitude)

						

				        # SAVE CURRENT DATA TO DATA SERVER     
						AttocubeApproach.push({'params':{'datasetname': datasetname, 'device': device, 'AOchannel': AOchannel, 'step_wait': step_wait, 'stage_min_voltage': stage_min_voltage, 'stage_max_voltage': stage_max_voltage, 'stage_voltage_steps': stage_voltage_steps, 'threshhold': threshhold, 'A_init': A_init},
							'title': 'Attocube Approach',
							'xlabel': 'Scanner Voltage (V)',
							'ylabel': 'Fork Amplitude (V)',
							'datasets': {'voltage': self.voltages,
										'amplitude': self.amplitudes
							}
						})

						if ((amplitude / A_init < threshhold) or (amplitude / A_init > (2 - threshhold))): #must decide how to threshhold
							self.engaged = True
							print('engaged')

					else:
						print('either engaged or voltage out of range')

				# Bring scanner to 0 for safety
				self.ao_task.write(0.0)


				# CLOSE DAQ TASKS
				self.ao_task.stop()
				self.ao_task.close()
				self.ao_task = None
				print('done')
				return

				flexSave(datasetName, 'AttocubeApproach', 'final') # after measurement finishes

