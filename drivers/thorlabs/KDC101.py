'''

Python driver for Thorlabs KDC101 actuator controlling cubes
Adapted from Thorlabs, Daniel Mark, Chris Egerstrom, Nazar Delegan

Sanskriti Chitransh, 2024-Feb-15

'''

import time
import os
import sys
from ctypes import *
from pylablib.devices import Thorlabs


class KDC():


    def __init__(self):

        # Connect devices
        self.kdcX = Thorlabs.KinesisMotor("27267826", is_rack_system=False)
        self.kdcY = Thorlabs.KinesisMotor("27267845", is_rack_system=False)

        # Open devices
        self.kdcX.open()
        self.kdcY.open()


    def set_x_position(self, new_pos): # position in encoder units

        self.kdcX.move_to(new_pos)


    def set_y_position(self, new_pos): # position in encoder units

        self.kdcY.move_to(new_pos)


    def get_x_position(self, new_pos): # position in encoder units

        self.kdcX.get_position()


    def get_y_position(self, new_pos): # position in encoder units

        self.kdcY.get_position()


    def x_wait_move(self):

        self.kdcX.wait_move()


    def y_wait_move(self):

        self.kdcY.wait_move()


    def __exit__(self, *args):

        # Close devices
        self.kdcX.close()
        self.kdcY.close()

# 10.33296
# 9.67505