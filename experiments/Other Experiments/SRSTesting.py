"""
Test IQ values to control SRS output

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
from drivers.ni.nidaq_final import NIDAQ


class SRS_Testing_Measurement:

    def SRSTesting(
        self, 
        freq: float,
        rf_power: float,
        I_On: float,
        I_Off: float,
        Q_On: float,
        Q_Off: float,
        period: int,
    ):

        print('starting')

        with InstrumentGateway() as gw:

            # SRS SETTINGS
            gw.sg.set_rf_amplitude(rf_power)    # set ouput power
            gw.sg.set_frequency(freq)           # set output frequency
            gw.sg.set_mod_state(True)           # enable modulation
            gw.sg.set_mod_type("7")             # QAM modulation
            gw.sg.set_mod_func("5")             # External modulation (from Swabian)
            gw.sg.set_rf_state("1")             # start output

            # START PULSESTREAMER
            print('STREAMING')
            gw.swabian.runSequenceInfinitely(Pulses(gw).SRS_TESTING(I_On, I_Off, Q_On, Q_Off, period))
            # for f in np.linspace(2.77e9, 2.97e9, 11):
            #     gw.sg.set_frequency(f)
            #     time.sleep(2)




