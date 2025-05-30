'''

Python driver for Newport ESP302 actuator STAGE
Adapted from serial driver for ESP300 and Ethernet driver for XPS 

Sanskriti Chitransh, 2024-May-10

'''

import time
import os
import sys
from ctypes import *

import NewportESP302

from lantz.messagebased import MessageBasedDriver
import serial 

class NESP:
    
    def __init__(self):

        self.esp = NewportESP302.ESP()  # open communication with controller
        self.espZ = self.esp.axis(1)    # open Z axis 
        self.espY = self.esp.axis(2)    # open Y axis    
        self.espX = self.esp.axis(3)    # open X axis


       