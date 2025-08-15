'''
NSpyre v0.6.1 compatible application to readout MLFI Threshold error signal

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


class ThresholdVsTime_Measurement:

	def auxoutvstime(
		self, 
		datasetname: str
		):

		with InstrumentGateway() as gw, DataSource(datasetname) as ThresholdVsTime:

			with NIDAQ() as mynidaq:

				# CREATING ANALOG DAQ TASKS
				self.ao_z_task = mynidaq.create_task()
				self.ao_z_task.ao_channels.add_ao_voltage_chan('Dev4/AO1')
				self.ao_x_task = mynidaq.create_task()
				self.ao_x_task.ao_channels.add_ao_voltage_chan('Dev4/AO2')
				self.ao_y_task = mynidaq.create_task()
				self.ao_y_task.ao_channels.add_ao_voltage_chan('Dev4/AO3')

				self.start_time = time.time()
				iterator = count()

				self.times = np.zeros(0)
				self.auxop = np.zeros(0)

				for i in iterator:
					value = gw.mfli.AUXOUT_read(3)
					if value != 0:
						self.ao_z_task.write(-0.1)
						# self.ao_x_task.write(-0.1)
						# self.ao_y_task.write(-0.1)
						print('RETRACTED!')
						time.sleep(1)

					self.auxop = np.append(self.auxop, value)
					self.times = np.append(self.times, time.time() - self.start_time)

					ThresholdVsTime.push({
				    	'params': {
							'datasetname': datasetname
						},
				    	'title': 'Threshold VS Time',
				    	'xlabel': 'time (s)',
				    	'ylabel': 'Threshold (V)',
				    	'datasets': {'time': self.times,
				    				'AUX': self.auxop,
				    	}
				    })

					time.sleep(time_per_point)

			flexSave(datasetName, 'ThresholdVsTime', 'final') # after measurement finishes


