import time
from itertools import count

import numpy as np

from nspyre import DataSource
from nspyre import InstrumentGateway

class TestPSMeasurement:

    def TestPS(self):
        
        with InstrumentGateway() as gw:

            # Swabian trigger for counting
            count_seq = [(20, 1), (1e9, 0)]

            # Swabian fake input to DAQ
            input_seq = [(20, 1), (1e6, 0)]

            # putting them together
            total_seq = gw.swabian.ps.createSequence()
            total_seq.setDigital(0, count_seq)
            total_seq.setDigital(7, input_seq)
            gw.swabian.runSequenceInfinitely(total_seq)
