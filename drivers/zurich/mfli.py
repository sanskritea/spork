'''

Zurich Instruments MFLI Driver for NSpyre v0.6.1

Sanskriti Chitransh, 2023-Oct-24

Adapted from Jonathan Marcks

'''


# std lib
# from collections import OrderedDict
import time

# 3rd party
import numpy as np
import zhinst.ziPython as ziPython # Zurich python

# lantz
# from lantz import Driver
# from lantz.core.feat import Feat, DictFeat
# from lantz.core.action import Action

# nspyre
# from lantz import Q_


class MFLI:

    # AUXOUT_TYPE = OrderedDict([
    #     ('manual', -1),
    #     ('demod_x', 0),
    #     ('demod_y', 1),
    #     ('demod_r', 2),
    #     ('demod_theta', 3),
    #     ('pid', 5), # PID output
    #     ('pid_shift', 9),
    #     ('pid_error', 10),
    #     ('tu_filtered', 11),
    #     ('tu_output', 13),
    # ])


    def __init__(self, ip="192.168.1.58", dev="dev5302"):

        self.mfli = ziPython.ziDAQServer(ip,8004,6)
        self.dev = dev
    

    def get_PID_setpoint(self): 	
        # get current PID loop setpoint value
        return float(self.mfli.getDouble("/%s/PIDS/0/setpoint" % (self.dev)))
    

    def set_PID_setpoint(self, value): 	 
        # set PID loop setpoint value
        self.mfli.set("/%s/PIDS/0/setpoint" % (self.dev), value)
            
    
    def PID_error(self):
    	return float(self.mfli.getDouble("/%s/PIDS/0/error" % (self.dev)))


    def get_PID_state(self):
    	return self.mfli.get("/%s/PIDS/0/enable" % (self.dev))


    def set_PID_state(self,value):
    	self.mfli.set("/%s/PIDS/0/enable" % (self.dev),value)


    def get_PID_values(self):
    	return float(self.mfli.getDouble("/%s/PIDS/0/value" % (self.dev)))

    
    def AUXOUT_read(self, channel):
    	return float(self.mfli.getDouble("/%s/AUXOUTS/%d/value" % (self.dev, channel)))
    

    def zero_OFFSET(self):
        self.mfli.set("/%s/AUXOUTS/2/offset" % (self.dev), 0)

    # @Feat(values=AUXOUT_TYPE)
    # def AUXOUT_select(self,channel):
    #     """
    #     AUXOUT output selection
    #     """
    #     return float(self.mfli.get("/%s/AUXOUTS/%d/outputselect" % (self.dev, channel)))
