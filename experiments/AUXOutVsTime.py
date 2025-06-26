'''
NSpyre v0.6.1 compatible application to readout MLFI PID error signal

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


class AUXOutVsTime_Measurement:

	def auxoutvstime(self, datasetname: str, time_per_point: float, aux_chan: int):

		with InstrumentGateway() as gw, DataSource(datasetname) as AUXOutVsTime:
			self.start_time = time.time()
			iterator = count()

			self.times = np.zeros(0)
			self.auxop = np.zeros(0)

			for i in iterator:
				self.auxop = np.append(self.auxop, gw.mfli.AUXOUT_read(aux_chan))
				self.times = np.append(self.times, time.time() - self.start_time)

				AUXOutVsTime.push({
			    						'params': {
			    									'datasetname': datasetname, 
			    									'time_per_point': time_per_point, 
			    									'aux_chan':aux_chan
			    								},
			    	'title': 'AOUX Out VS Time',
			    	'xlabel': 'time (s)',
			    	'ylabel': 'AUX OP (V)',
			    	'datasets': {'time': self.times,
			    				'AUX': self.auxop,
			    	}
			    })

				time.sleep(time_per_point)

			flexSave(datasetName, 'AUXOutVsTime', 'final') # after measurement finishes


