# -*- coding: utf-8 -*-
"""
    lantz.drivers.newport.xpsq8
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Implementation of XPS Q8 controller
    NOTE: XPS_Q8_drivers.py must be placed within the same directory
    as this script
    Author: Kevin Miao, Cyrus Zeledon, Elena Glen, Jacob Feder
    Date: 07/06/2020
"""

from lantz.driver import Driver
from lantz import DictFeat

from collections import OrderedDict
from time import sleep

from lantz.core import Action, Feat, Driver

from lantz.drivers.newport import XPS_Q8_drivers


class XPSQ8(Driver):


    def __init__(self, address, channels=['GROUP1', 'GROUP2', 'GROUP3'], port=5001, timeout=20.0):
        super(XPSQ8, self).__init__()
        self._xps = XPS_Q8_drivers.XPS()
        self._socket_id = self._xps.TCP_ConnectToServer(address, port, timeout)
        self.channels = channels
        for i in range(len(channels)):
            self._xps.GroupInitialize(0,self.channels[i])
            self._xps.GroupHomeSearch(0,self.channels[i])


    @Action()
    def reboot(self):
        self._xps.Reboot(self._socket_id)
        return

    @Action()
    def reinitialize(self):
        for i in range(len(self.channels)):
                self._xps.GroupMotionEnable(0,self.channels[i])

   @DictFeat()
   def travel_limits(self, channel):
       retval = self._xps.PositionerUserTravelLimitsGet(self._socket_id, channel)

    @Action()
    def home(self, channel):
        retval = self._xps.GroupHomeSearch(self._socket_id, channel)

    @DictFeat(units='mm')
    def abs_position(self, channel, position):
        retval = self._xps.GroupMoveAbsolute(self._socket_id, channel, [position])
        if retval == '[-42,'']':
            reinitialize()
            print('Jog value out of range. Try again')
        return
    
        if retval == '[-22,'']':
            reinitialize()
            print('Not allowed action. Try again')
        return
    
        if retval == '[-17,'']':
            print('Parameter out of range. Try again')
        return

    @DictFeat(units = 'mm', keys=['GROUP1', 'GROUP2','GROUP3'])
    def abs_position(self, channel):
        return self._xps.GroupPositionCurrentGet(self._socket_id, channel, 1)[1]

    @abs_position.setter
    def abs_position(self, channel, val):
        retval = self._xps.GroupMoveAbsolute(self._socket_id, channel, [val])
        if retval == '[-42,'']':
            reinitialize()
            print('Jog value out of range. Try again')
        return

        if retval == '[-22,'']':
            reinitialize()
            print('Not allowed action. Try again')
        return

        if retval == '[-17,'']':
            print('Parameter out of range. Try again')
        return


    @Action()
    def rel_position(self, channel, dposition):
        retval = self._xps.GroupMoveRelative(self._socket_id, channel, [dposition])
        if retval == '[-42,'']':
            reinitialize()
            print('Jog value out of range. Try again')
        return

        if retval == '[-22,'']':
            reinitialize()
            print('Not allowed action. Try again')
        return

        if retval == '[-17,'']':
            print('Parameter out of range. Try again')
        return

    @Action()
    def jog(self, channel, velocity, acceleration):
        retval = self._xps.GroupJogParametersSet(self._socket_id, channel, [velocity], [acceleration])
        if retval == '[-42,'']':
            reinitialize()
            print('Jog value out of range. Try again')
        return

        if retval == '[-22,'']':
            reinitialize()
            print('Not allowed action. Try again')
        return

        if retval == '[-17,'']':
            print('Parameter out of range. Try again')
        return

    @Feat(units = 'mm' )
    def get_abs_position(self, channel):
        retval = self._xps.GroupPositionCurrentGet (self._socket_id, channel, 1)

def main():
    import logging
    import sys
    from lantz.log import log_to_screen
    import numpy as np
    log_to_screen(logging.CRITICAL)
    res_name = sys.argv[1]
    with XPSQ8(res_name) as inst:
        value = inst._xps.GroupJogParametersGet(inst._socket_id, 'Group1.Pos', 1)
        ret = inst._xps.GroupJogParametersSet(inst._socket_id, 'Group2.Pos', [-0.0,], [1.0,])
        print(ret)
        value = inst._xps.GroupJogParametersGet(inst._socket_id, 'Group2.Pos', 1)
        print(value)
        return
        positions = np.linspace(-12.5, 12.5, 20)
        for val in positions:
            print(val)
            inst.abs_position['Group1.Pos'] = val
            print(inst.abs_position['Group1.Pos'])

if __name__ == '__main__':
    main()
