# -*- coding: utf-8 -*-
"""
    lantz.drivers.stanford.sg396
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Implementation of SG396 signal generator

    Author: Kevin Miao & Berk Diler

    Date: 12/15/2015 & 8/21/17

    Modified GS,LW 12/6/22
    Modified PMN 7/11/23
"""

import numpy as np
# from lantz import Action, Feat, DictFeat, ureg
from lantz.core import MessageBasedDriver
# from lantz.messagebased import MessageBasedDriver
from collections import OrderedDict

class SG396(MessageBasedDriver):

        DEFAULTS = {
            'COMMON': {
                'write_termination': '\r\n',
                'read_termination': '\r\n',
            }
        }

        MODULATION_TYPE = OrderedDict([
            ('AM', 0),
            ('FM', 1),
            ('Phase',2),
            ('Sweep',3),
            ('Pulse',4),
            ('Blank',5),
            ('QAM',7),
            ('CPM',8),
            ('VSB',9)
        ])

        MODULATION_FUNCTION = OrderedDict([
            ('sine', 0),
            ('ramp', 1),
            ('triangle', 2),
            ('square', 3),
            ('noise', 4),
            ('external', 5)
        ])

        # Signal synthesis commands

        def get_lf_amplitude(self):
            """
            low frequency amplitude (BNC output)
            """
            return float(self.query('AMPL?'))


        def set_lf_amplitude(self, value):
            self.write('AMPL{:.2f}'.format(value))


        def get_rf_amplitude(self):
            """
            RF amplitude (Type N output)
            """
            return float(self.query('AMPR?'))
        

        def set_rf_amplitude(self, value):
            self.write('AMPR{:.2f}'.format(value))


        def get_lf_state(self):
            """
            low frequency output state
            """
            return self.query('ENBL?')


         # note: can only be on if freq in range
        def set_lf_state(self, value):
            self.write('ENBL{:s}'.format(value))
        

        def get_rf_state(self):
            """
            RF output state
            """
            return self.query('ENBR?')


        # note: can only be on if freq in range
        def set_rf_state(self, value):
            self.write('ENBR{:s}'.format(value))
  

        def get_frequency(self):
            """
            signal frequency
            """
            return float(self.query('FREQ?'))

    
        def set_frequency(self, value):
            self.write('FREQ{:.2f}'.format(value))


        def get_lf_offset(self):
            """
            low frequency offset voltage
            """
            return self.query('OFSL?')


        def set_lf_offset(self, value):
            self.write('OFSL{:.2f}'.format(value))


        def get_phase(self):
            """
            carrier phase
            """
            return self.query('PHAS?')


        def set_phase(self, value):
            self.write('PHAS{:.2f}'.format(value))


        def set_rel_phase(self):
            """
            sets carrier phase to 0 degrees
            """
            self.write('RPHS')


        def get_mod_state(self):
            """
            Modulation State
            """
            return int(self.query('MODL?'))


        def set_mod_state(self, value):
            self.write('MODL {}'.format(value))


        def get_mod_type(self):
            """
            Modulation State
            """
            return int(self.query('TYPE?'))


        def set_mod_type(self, value):
            self.write('TYPE {}'.format(value))


        def get_mod_func(self):
            """
            Modulation Function
            """
            return int(self.query('MFNC?'))


        def set_mod_func(self, value):
            self.write('MFNC {}'.format(value))


        def get_mod_rate(self):
            """
            Modulation Rate
            """
            return float(self.query('RATE?'))


        def set_mod_rate(self, val):
            self.write('RATE {}'.format(val))


        def get_AM_mod_depth(self):
            """
            AM Modulation Depth
            """
            return float(self.query('ADEP?'))


        def set_AM_mod_depth(self, val):
            self.write('ADEP {}'.format(val))


        def get_FM_mod_dev(self):
            """
            FM Modulation Deviation
            """
            return float(self.query("FDEV?"))


        #limits=(0.1, 8.e6)) Hz
        def set_FM_mod_dev(self, val):
            self.write('FDEV {}'.format(val))