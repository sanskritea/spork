'''

ANM200/ANC300 driver for Attocube Scanner.  

Sanskriti Chitransh, 2023-Oct-25

Note: This uses a serial connection to the controller and directly writes the commands to it.

'''

import numpy as np
import time
import math
from itertools import cycle
import logging
import scipy as sp
from scipy import signal
import datetime as Dt


import serial
from serial.tools import list_ports
from serial import SerialException

# from rpyc.utils.classic import obtain 

import time
# from ctypes import POINTER, c_double, c_uint, c_void_p, pointerr 


class scanner():

	def __init__(self, port='COM4'):
        
		self.port = port

		if self.port != None:
			try:
				device = serial.Serial(self.port, timeout=1)
			except Exception as err:
				self.address = None
				raise SerialException(f"{self.port} not accesible.") from err

		# ENABLE DC INPUT 
		# device.write(b'getdci 1\r\n')
		# device.write(b'getdci 2\r\n')
		# device.write(b'getdci 3\r\n')

		# i = 0
		# while i == 0:
		# 	line = device.readline()
		# 	if line == b'':
		# 		i = 1
		# 	else:
		# 		print(line)

		device.write(b'setdci 1 on\r\n')
		device.write(b'setdci 2 on\r\n')
		device.write(b'setdci 3 on\r\n')

		# i = 0
		# while i == 0:
		# 	line = device.readline()
		# 	if line == b'':
		# 		i = 1
		# 	else:
		# 		print(line)

		return


    # def initialize(self, devNo=None):

    #     if not devNo is None: self.devNo = devNo
    #     device = c_void_p()
    #     self.check_error(self.lib.connect(self.dev_no, pointer(device)))
    #     self.device = device


    # def finalize(self):

    #     self.check_error(self.lib.disconnect(self.device))
    #     self.device = None


    # def check_error(self, err):

    #     if err != 0:
    #         raise Exception("Driver Error {}: {}".format(err, self.RETURN_STATUS[err]))
    #     return


    # def get_frequency(self, axis):
    #     ret_freq = c_double()
    #     self.check_error(self.lib.getFrequency(self.device, axis, pointer(ret_freq)))
    #     return ret_freq.value


    # def set_frequency(self, axis, freq):
    #     self.check_error(self.lib.setFrequency(self.device, axis, freq))
    #     return err


    # def get_position(self, axis):
    #     ret_pos = c_double()
    #     self.check_error(self.lib.getPosition(self.device, axis, pointer(ret_pos)))
    #     return ret_pos.value * 1e6


    # def set_position(self, axis, pos):
    #     return self.absolute_move(axis, pos * 1e-6)


    # def capacitance(self, axis):
    #     ret_c = c_double()
    #     self.check_error(self.lib.measureCapacitance(self.device, axis, pointer(ret_c)))
    #     return ret_c.value


    # def get_status(self, axis):
    #     status_names = [
    #         'connected',
    #         'enabled',
    #         'moving',
    #         'target',
    #         'eot_fwd',
    #         'eot_bwd',
    #         'error',
    #     ]
    #     status_flags = [c_uint() for _ in range(7)]
    #     status_flags_p = [pointer(flag) for flag in status_flags]
    #     self.check_error(self.lib.getAxisStatus(self.device, axis, *status_flags_p))

    #     ret = dict()
    #     for status_name, status_flag in zip(status_names, status_flags):
    #         ret[status_name] = True if status_flag.value else False
    #     return ret


    # def dc_bias(self, axis, voltage):
    #     self.check_error(self.lib.setDcVoltage(self.device, axis, voltage))
    #     return








