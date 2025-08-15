"""
This is an application to run CW ODMR on Jasper

Copyright (c) April 2023, C. Egerstrom
All rights reserved.

This work is licensed under the terms of the 3-Clause BSD license.
For a copy, see <https://opensource.org/licenses/BSD-3-Clause>.

Modified: PMN July '23
"""


import time
from itertools import count

import numpy as np

from nspyre import DataSource
from nspyre import InstrumentGateway

from rpyc.utils.classic import obtain

from experiments.NewPulses import Pulses
# from PulsePatterns import Pulses


class ATTENTesting:

    def ATTENTest(self, atten_voltage: float):

        with InstrumentGateway() as gw:

            # test ATTENUATOR laser intensity modulation
            gw.nidaq.laser_power_atten(atten_voltage)

