
'''
Pulse sequence class with all sequences used for Nanoscale NMR experiments for new nspyre

Edited & rewritten by Evan Villafranca - September 10, 2022
'''

from time import time
import numpy as np
# import pandas as pd
from math import sin, cos, radians
import typing as t

from rpyc.utils.classic import obtain

from pulsestreamer import Sequence
from drivers.ni.nidaq_final import NIDAQ

class Pulses():
    '''
    ALL UNITS: [ns]
    '''

    clock_time = 10
    init_time = 2000
    readout_time = 500
    laser_lag = 300
    singlet_decay = 1000


    def __init__(self, gateway):

        self.channel_dict = {"DAQ_CLOCK": 0, "AOM": 1, "SWITCH": 4,"FAKE_DAQ_INPUT": 7}  
        
        # self.laser_time = laser_time
        # self.initial_delay = initial_delay
        # self.singlet_decay = singlet_decay
        # self.readout_time = readout_time
        # self.MW_buffer_time = MW_buffer_time
        # self.probe_time = probe_time
        # self.clock_time = clock_time
        # self.sampling_time = sampling_time
        # self.trig_spot = trig_spot
        # self.awg_trig_time = awg_trig_time
        # self.awg_pulse_delay = awg_pulse_delay
        # self.rest_time_btw_seqs = rest_time_btw_seqs

        self.Pulser = gateway.swabian.ps
        self.sequence = Sequence()

        # self.latest_streamed = pd.DataFrame({})
        self.total_time = 0 #update when a pulse sequence is streamed

        ## Analog voltage levels to set for sig gen I and Q. If no crosstalk should be 0.5 and 0.
        self.IQ0 = [-0.007, -0.007]
        self.IQ = self.IQ0

        self.IQpx = [0.5, 0]
        # self.IQnx = [-0.4862,-0.0030]

        # self.IQpy = [-0.0038, 0.4870]
        # self.IQny = [-0.0038, -0.4870]
        
        # self.IQboth = [0.4862, 0.4870]
        # self.IQtest = [.95, 0]


    def PiHalf(self, axis, pi_half_time):
        iq_on = pi_half_time
        if axis == 'x':
            mw_I_on = (iq_on, self.IQpx[0])
            mw_Q_on = (iq_on, self.IQpx[1])
        elif axis == '-x':
            mw_I_on = (iq_on, self.IQnx[0])
            mw_Q_on = (iq_on, self.IQnx[1])
        elif axis == 'y':
            mw_I_on = (iq_on, self.IQpy[0])
            mw_Q_on = (iq_on, self.IQpy[1])
        elif axis == '-y':
            mw_I_on = (iq_on, self.IQny[0])
            mw_Q_on = (iq_on, self.IQny[1])
        
        return mw_I_on, mw_Q_on
        
    def Pi(self, axis, pi_time):
        iq_on = pi_time
        if axis == 'x':
            mw_I_on = (iq_on, self.IQpx[0])
            mw_Q_on = (iq_on, self.IQpx[1])
        elif axis == '-x':
            mw_I_on = (iq_on, self.IQnx[0])
            mw_Q_on = (iq_on, self.IQnx[1])
        elif axis == 'y':
            mw_I_on = (iq_on, self.IQpy[0])
            mw_Q_on = (iq_on, self.IQpy[1])
        elif axis == '-y':
            mw_I_on = (iq_on, self.IQny[0])
            mw_Q_on = (iq_on, self.IQny[1])
        
        return mw_I_on, mw_Q_on


    def PSDAQtesting(self):
        # Swabian trigger for counting
        count_seq = [(20, 1), (20, 0)]

        # Swabian fake input to DAQ
        input_seq = [(20, 1), (20, 0)]

        # putting them together
        total_seq = self.Pulser.createSequence()
        total_seq.setDigital(0, count_seq)
        total_seq.setDigital(7, input_seq)

        return total_seq


    def fakeDAQinput(self):
        # fake input signal for DAQ
        fake_daq_input_seq = [(20, 1), (1e3, 0)]

        seq = self.Pulser.createSequence()
        seq.setDigital(7, fake_daq_input_seq)
        return seq


    def laser_on(self):

        seq = self.Pulser.createSequence()
        laser_seq = [(-1, 1)]
        seq.setDigital(1, laser_seq)

        return seq


    def counting_trigger(self, buffer_rate):

        seq = self.Pulser.createSequence()

        laser_seq = [(-1, 1)]

        daq_seq = [(10, 1), (int(1e9 / buffer_rate), 0)]

        seq.setDigital(0, daq_seq)
        seq.setDigital(1, laser_seq)

        return seq


    def SRS_TESTING(self, ivalon, ivaloff, qvalon, qvaloff, period):

        seq = self.Pulser.createSequence()
        mw_on = period
        mw_off = period
        mw_I_seq = [(mw_off, ivaloff), (mw_on, ivalon)]
        mw_Q_seq = [(mw_off, qvaloff), (mw_on, qvalon)]
        seq.setAnalog(0, mw_I_seq) # mw_I
        seq.setAnalog(1, mw_Q_seq) # mw_Q

        return seq


    def AOM_Lag(self, tau_time, clock_time, init_time, readout_time, rise_laser_off, fall_laser_off):
        '''
        Optical T1 sequence (without any MW), all sequences have buffer/padding to have same dute cycle for the AOM
        '''
        ## Initialize with green, then measure the signal
        ## counts from NV after variable tau delay

        # create sequence object
        seq = self.Pulser.createSequence()
        tau_time = int(tau_time)

        # Laser
        laser_off_1 = rise_laser_off
        laser_on_1 = init_time
        laser_off_2 = fall_laser_off
        laser_seq = [(laser_off_1, 0), (laser_on_1, 1), (laser_off_2, 0)]
        print('laser duration ', (laser_off_1 + laser_on_1 + laser_off_2) )

        # DAQ
        daq_off_1 = tau_time # tau_time
        daq_on_1 = clock_time
        daq_off_2 = readout_time
        daq_on_2 = clock_time
        daq_off_3 = rise_laser_off + init_time + fall_laser_off - tau_time - (2 * clock_time) - readout_time
        daq_clock_seq = [(daq_off_1, 0), (daq_on_1, 1), (daq_off_2, 0), (daq_on_2, 1), (daq_off_3, 0)]
        print('daq duration ', (daq_off_1 + daq_on_1 + daq_off_2 + daq_on_2 + daq_off_3))

        # assign sequences to respective channels
        seq.setDigital(0, daq_clock_seq) # DAQ clock -- record data
        seq.setDigital(1, laser_seq)

        return seq


    def CW_ODMR(self, cw_odmr_probe_time):
    
        # create sequence object
        seq = self.Pulser.createSequence()

        # daq counting pulses
        daq_clock_seq = [(clock_time, 1), (cw_odmr_probe_time - clock_time, 0), (clock_time, 1), (cw_odmr_probe_time - clock_time, 0)]

        # microwave switch 
        switch_ttl_seq = [(cw_odmr_probe_time, 1), (cw_odmr_probe_time, 0)]
        # I_seq = [(probe_time, 0.5), (probe_time, -0.007)]
        # Q_seq = [(probe_time, 0), (probe_time, -0.007)]

        # turn on laser throughout the measurement
        laser_seq = [(-1, 1)]

        # assign sequences to respective channels
        seq.setDigital(0, daq_clock_seq) # DAQ clock -- record data
        seq.setDigital(1, laser_seq)
        seq.setDigital(4, switch_ttl_seq)
        # seq.setAnalog(0, I_seq)
        # seq.setAnalog(1, Q_seq)


        return seq


    def RABI(self, mw_time, max_MW_time):
        '''
        Rabi sequence
        '''
        ## Run a MW pulse of varying duration, then measure the signal
        ## and reference counts from NV.

        # create sequence object
        seq = self.Pulser.createSequence()

        # define laser sequence
        laser_on_1 = init_time
        laser_off_1 = singlet_decay + max_MW_time
        laser_on_2 = init_time
        laser_off_2 = singlet_decay + max_MW_time
        laser_seq = [(laser_on_1, 1), (laser_off_1, 0), (laser_on_2, 1), (laser_off_2, 0)] 

        # define DAQ counting sequence triggers
        daq_off_0 = laser_lag 
        daq_on_1 = clock_time 
        daq_off_1 = probe_time
        daq_on_2 = clock_time
        daq_off_2 = init_time - (2 * clock_time) - probe_time - laser_lag + singlet_decay + max_MW_time + laser_lag
        daq_on_3 = clock_time
        daq_off_3 = probe_time
        daq_on_4 = clock_time
        daq_off_4 = init_time - (2 * clock_time) - probe_time - laser_lag + singlet_decay + max_MW_time 
        daq_clock_seq = [(daq_off_0, 0), (daq_on_1, 1), (daq_off_1, 0), (daq_on_2, 1), (daq_off_2, 0), (daq_on_3, 1), (daq_off_3, 0), (daq_on_4, 1), (daq_off_4, 0)]

        # define sequence for MW (switch)
        switch_off_0 = init_time + singlet_decay 
        switch_on_1 = mw_time
        switch_off_1 =  max_MW_time - mw_time + init_time + singlet_decay + max_MW_time 
        switch_ttl_seq = [(switch_off_0, 0), (switch_on_1, 1), (switch_off_1, 0)]
        I_seq = [(switch_off_0, 0), (switch_on_1, 0.5), (switch_off_1, 0)]
        Q_seq = [(switch_off_0, 0), (switch_on_1, 0), (switch_off_1, 0)]

        # assign sequences to respective channels
        seq.setDigital(0, daq_clock_seq) # DAQ clock -- record data
        seq.setDigital(1, laser_seq)
        seq.setDigital(4, switch_ttl_seq)
        # seq.setAnalog(0, I_seq)
        # seq.setAnalog(1, Q_seq)

        return seq


    def PULSED_ODMR(self, pi_time):

        # create sequence object
        seq = self.Pulser.createSequence()

        # define laser sequence
        laser_on_1 = init_time
        laser_off_1 = singlet_decay + pi_time
        laser_on_2 = init_time
        laser_off_2 = singlet_decay + pi_time
        laser_seq = [(laser_on_1, 1), (laser_off_1, 0), (laser_on_2, 1), (laser_off_2, 0)]

        # define DAQ counting sequence triggers
        daq_off_0 = laser_lag 
        daq_on_1 = clock_time 
        daq_off_1 = probe_time
        daq_on_2 = clock_time
        daq_off_2 = init_time - (2 * clock_time) - probe_time - laser_lag + singlet_decay + pi_time + laser_lag
        daq_on_3 = clock_time
        daq_off_3 = probe_time
        daq_on_4 = clock_time
        daq_off_4 = init_time - (2 * clock_time) - probe_time - laser_lag + singlet_decay + pi_time 
        daq_clock_seq = [(daq_off_0, 0), (daq_on_1, 1), (daq_off_1, 0), (daq_on_2, 1), (daq_off_2, 0), (daq_on_3, 1), (daq_off_3, 0), (daq_on_4, 1), (daq_off_4, 0)]

        # define sequence for MW (switch)
        switch_off_0 = init_time + singlet_decay
        switch_on_1 = pi_time
        switch_off_1 = init_time + singlet_decay + pi_time 
        switch_ttl_seq = [(switch_off_0, 0), (switch_on_1, 1), (switch_off_1, 0)]
        I_seq = [(switch_off_0, 0), (switch_on_1, 0.5), (switch_off_1, 0)]
        Q_seq = [(switch_off_0, 0), (switch_on_1, 0.5), (switch_off_1, 0)]

        # assign sequences to respective channels
        seq.setDigital(0, daq_clock_seq) # DAQ clock -- record data
        seq.setDigital(1, laser_seq)
        seq.setDigital(4, switch_ttl_seq)
        # seq.setAnalog(0, I_seq)
        # seq.setAnalog(1, Q_seq)

        return seq


    def REINIT_TEST(self, clock_time, init_time, readout_time, singlet_decay, pi_time, num_readouts):

        # create sequence object
        seq = self.Pulser.createSequence()

        # define laser sequence
        laser_on_1 = 50000
        laser_off_1 = singlet_decay
        laser_on_2 = init_time
        laser_off_2 = singlet_decay
        laser_on_3 = laser_on_1
        laser_off_3 = singlet_decay
        laser_on_4 = init_time
        laser_off_4 = singlet_decay
        laser_seq = [(laser_on_1, 1), (laser_off_1, 0)]

        # define DAQ sequence
        daq_off_0 = laser_on_1 + singlet_decay
        daq_clock_seq = [(daq_off_0, 0)]
        readout_seq = [(clock_time, 1), (readout_time, 0), (clock_time, 1)]
        for i in range(num_readouts):
            daq_clock_seq.append(readout_seq)
        daq_clock_seq.append([(laser_off_2 + laser_on_3 + laser_off_3, 0)])
        for i in range(num_readouts):
            daq_clock_seq.append(readout_seq)
        daq_clock_seq.append([(laser_off_4, 0)])
        
        # define sequence for MW (switch)
        switch_ttl_seq = [(laser_on_1 + laser_off_1 + laser_on_2+ laser_off_2 + laser_on_3 + singlet_decay - pi_time, 0), (pi_time, 1), (init_time + singlet_decay, 0)]


        # assign sequences to respective channels
        seq.setDigital(0, daq_clock_seq) # DAQ clock -- record data
        seq.setDigital(1, laser_seq)
        seq.setDigital(4, switch_ttl_seq)

        return seq


    def RAMSEY(self, clock_time, init_time, laser_lag, probe_time, singlet_decay, pi_half_time, tau_time, max_tau_time):

        # create sequence object
        seq = self.Pulser.createSequence()
        padding = 1000 + int(max_tau_time - tau_time) 

        # define laser sequence
        laser_off_0 = padding
        laser_on_1 = init_time
        laser_off_1 = singlet_decay + (2 * pi_half_time) + tau_time
        laser_on_2 = init_time
        laser_seq = [(laser_off_0, 0), (laser_on_1, 1), (laser_off_1, 0), (laser_on_2, 1)]

        # define sequence for MW (switch)
        switch_off_0 = padding + init_time + singlet_decay
        switch_on_1 = pi_half_time
        switch_off_1 = tau_time 
        switch_on_2 = pi_half_time
        switch_off_2 = init_time

        switch_ttl_seq = [(switch_off_0, 0), (switch_on_1, 1), (switch_off_1, 0), (switch_on_2, 1), (switch_off_2, 0)]

        # define DAQ counting sequence triggers
        daq_off_0 = padding + laser_lag 
        daq_on_1 = clock_time 
        daq_off_1 = probe_time
        daq_on_2 = clock_time
        daq_off_2 = init_time - (2 * clock_time) - probe_time - laser_lag + singlet_decay + ( 2 * pi_half_time) + tau_time + laser_lag
        daq_on_3 = clock_time
        daq_off_3 = probe_time
        daq_on_4 = clock_time
        daq_off_4 = init_time - (2 * clock_time) - probe_time - laser_lag
        daq_clock_seq = [(daq_off_0, 0), (daq_on_1, 1), (daq_off_1, 0), (daq_on_2, 1), (daq_off_2, 0), (daq_on_3, 1), (daq_off_3, 0), (daq_on_4, 1), (daq_off_4, 0)]

        # assign sequences to respective channels
        seq.setDigital(0, daq_clock_seq) # DAQ clock -- record data
        seq.setDigital(1, laser_seq)
        seq.setDigital(4, switch_ttl_seq)

        return seq    


    def OPTICAL_T1(self, tau_time, clock_time, init_time, readout_time, laser_lag, singlet_decay, tau_max):
        '''
        Optical T1 sequence (without any MW), all sequences have buffer/padding to have same dute cycle for the AOM
        '''
        ## Initialize with green, then measure the signal
        ## counts from NV after variable tau delay

        print('Checking wait time from Bayesian T1 experiment ', tau_time)

        # create sequence object
        seq = self.Pulser.createSequence()
        tau_time = int(tau_time)
        if tau_max > 0:
            padding_time = 1000 + tau_max - tau_time    # sequence padding, extra laser off time of 1us to prevent first point jumps, for normal T1
        else:
            padding_time = 1000   # for Bayesian T1

        # Laser
        laser_off_1 = padding_time 
        laser_on_1 = init_time
        laser_off_2 = singlet_decay + tau_time
        laser_on_2 = init_time
        laser_seq = [(laser_off_1, 0), (laser_on_1, 1), (laser_off_2, 0), (laser_on_2, 1)] 

        # DAQ
        daq_off_1 = padding_time + init_time + singlet_decay + tau_time + laser_lag 
        daq_on_1 = clock_time
        daq_off_2 = readout_time
        daq_on_2 = clock_time
        daq_off_3 = init_time - (2 * clock_time) - readout_time - laser_lag 
        daq_clock_seq = [(daq_off_1, 0), (daq_on_1, 1), (daq_off_2, 0), (daq_on_2, 1), (daq_off_3, 0)]

        # assign sequences to respective channels
        seq.setDigital(0, daq_clock_seq) # DAQ clock -- record data
        seq.setDigital(1, laser_seq)

        return seq


    def MW_T1(self, tau_time, clock_time, init_time, probe_time, laser_lag, singlet_decay, tau_max, pi_time):

        # create sequence object
        seq = self.Pulser.createSequence()

        # define laser sequence
        laser_off_0 = 2000 + int(tau_max - tau_time)
        laser_on_1 = init_time
        laser_off_1 = singlet_decay + pi_time + tau_time
        laser_on_2 = init_time
        laser_off_2 = 2000 + int(tau_max - tau_time) 
        laser_on_3 = init_time
        laser_off_3 = singlet_decay + tau_time
        laser_on_4 = init_time
        laser_seq = [(laser_off_0, 0), (laser_on_1, 1), (laser_off_1, 0), (laser_on_2, 1), (laser_off_2, 0), (laser_on_3, 1), (laser_off_3, 0), (laser_on_4, 1)] 

        # define DAQ counting sequence triggers
        daq_off_0 = 2000 + tau_max + init_time + singlet_decay + pi_time + laser_lag
        daq_on_1 = clock_time 
        daq_off_1 = probe_time
        daq_on_2 = clock_time
        daq_off_2 = init_time - (2 * clock_time) - probe_time - laser_lag + 2000 + tau_max + init_time + singlet_decay + laser_lag
        daq_on_3 = clock_time
        daq_off_3 = probe_time
        daq_on_4 = clock_time
        daq_off_4 = init_time - (2 * clock_time) - probe_time - laser_lag 
        daq_clock_seq = [(daq_off_0, 0), (daq_on_1, 1), (daq_off_1, 0), (daq_on_2, 1), (daq_off_2, 0), (daq_on_3, 1), (daq_off_3, 0), (daq_on_4, 1), (daq_off_4, 0)]

        # define sequence for MW (switch)
        switch_off_0 = 2000 + int(tau_max - tau_time) + init_time + singlet_decay 
        switch_on_1 = pi_time
        switch_off_1 = tau_time + init_time + 2000 + tau_max + init_time + singlet_decay + init_time
        switch_ttl_seq = [(switch_off_0, 0), (switch_on_1, 1), (switch_off_1, 0)]

        # assign sequences to respective channels
        seq.setDigital(0, daq_clock_seq) # DAQ clock -- record data
        seq.setDigital(1, laser_seq)
        seq.setDigital(4, switch_ttl_seq)

        # # define laser sequence
        # laser_off_0 = 2000  
        # laser_on_1 = init_time
        # laser_off_1 = singlet_decay + pi_time + tau_time
        # laser_on_2 = init_time
        # laser_off_2 = 2000 
        # laser_on_3 = init_time
        # laser_off_3 = singlet_decay + pi_time + int(tau_max - tau_time)
        # laser_on_4 = init_time
        # laser_seq = [(laser_off_0, 0), (laser_on_1, 1), (laser_off_1, 0), (laser_on_2, 1), (laser_off_2, 0), (laser_on_3, 1), (laser_off_3, 0), (laser_on_4, 1)] 

        # # define DAQ counting sequence triggers
        # daq_off_0 = 2000 + init_time + singlet_decay + pi_time + tau_time + laser_lag
        # daq_on_1 = clock_time 
        # daq_off_1 = probe_time
        # daq_on_2 = clock_time
        # daq_off_2 = init_time - (2 * clock_time) - probe_time - laser_lag + 2000 + init_time + singlet_decay + pi_time + int(tau_max - tau_time) + laser_lag
        # daq_on_3 = clock_time
        # daq_off_3 = probe_time
        # daq_on_4 = clock_time
        # daq_off_4 = init_time - (2 * clock_time) - probe_time - laser_lag 
        # daq_clock_seq = [(daq_off_0, 0), (daq_on_1, 1), (daq_off_1, 0), (daq_on_2, 1), (daq_off_2, 0), (daq_on_3, 1), (daq_off_3, 0), (daq_on_4, 1), (daq_off_4, 0)]

        # # define sequence for MW (switch)
        # switch_off_0 = 2000 + init_time + singlet_decay 
        # switch_on_1 = pi_time
        # switch_off_1 = tau_time + init_time + 2000 + init_time + singlet_decay + pi_time + int(tau_max - tau_time) + init_time
        # switch_ttl_seq = [(switch_off_0, 0), (switch_on_1, 1), (switch_off_1, 0)]

        # # assign sequences to respective channels
        # seq.setDigital(0, daq_clock_seq) # DAQ clock -- record data
        # seq.setDigital(1, laser_seq)
        # seq.setDigital(4, switch_ttl_seq)

        return seq


    def T2_STAR(self, tau_time, clock_time, init_time, probe_time, laser_lag, singlet_decay, pi_half_time, tau_max_time):
        '''
        MW (differential) T1 sequence for the two longitudinal relaxation rates measured with different wait times
        '''
        # create sequence objects for MW on and off blocks
        seq = self.Pulser.createSequence()
        pi_time = int(2 * pi_half_time)
        padding_time = int(tau_max_time - tau_time) + 2000

        # laser sequence
        laser_seq = [
            (padding_time, 0), (init_time, 1), (singlet_decay + pi_half_time + tau_time + pi_half_time, 0), (init_time, 1),
            (padding_time, 0), (init_time, 1), (singlet_decay + pi_half_time + tau_time + int(3 * pi_half_time), 0), (init_time, 1)
            ]

        # define sequence structure for DAQ counting triggers
        daq_clock_seq = [
            (padding_time + laser_lag, 0), 
            (clock_time, 1), (probe_time, 0), (clock_time, 1),
            (init_time - (2 * clock_time) - probe_time - laser_lag 
                + singlet_decay + pi_half_time + tau_time + pi_half_time + laser_lag, 0), 
            (clock_time, 1), (probe_time, 0), (clock_time, 1), 
            (init_time - (2 * clock_time) - probe_time - laser_lag 
                + padding_time + laser_lag, 0), 
            (clock_time, 1), (probe_time, 0), (clock_time, 1),
            (init_time - (2 * clock_time) - probe_time - laser_lag 
                + singlet_decay + pi_half_time + tau_time + int(3 * pi_half_time) + laser_lag, 0), 
            (clock_time, 1), (probe_time, 0), (clock_time, 1), 
            (init_time - (2 * clock_time) - probe_time - laser_lag, 0)
            ]

        # define sequence for MW (switch)
        switch_ttl_seq = [(padding_time + init_time + singlet_decay , 0), (pi_half_time, 1), (tau_time, 0), (pi_half_time, 1), (init_time + padding_time + init_time + singlet_decay, 0), (pi_half_time, 1), (tau_time, 0), (int(3 * pi_half_time), 1), (init_time, 0)]

        # assign sequences to respective channels for seq_on
        seq.setDigital(0, daq_clock_seq) 
        seq.setDigital(1, laser_seq) 
        seq.setDigital(4, switch_ttl_seq)

        return seq


    def HAHN(self, tau_time, clock_time, init_time, probe_time, laser_lag, singlet_decay, pi_half_time, tau_max_time):
        '''
        MW (differential) T1 sequence for the two longitudinal relaxation rates measured with different wait times
        '''
        # create sequence objects for MW on and off blocks
        seq = self.Pulser.createSequence()
        pi_time = int(2 * pi_half_time)
        tau_half_time = int(tau_time / 2)
        padding_time = int(tau_max_time - tau_time) + 2000

        # laser sequence
        laser_seq = [
            (padding_time, 0), (init_time, 1), (singlet_decay + pi_half_time + tau_half_time + pi_time + tau_half_time + pi_half_time, 0), (init_time, 1),
            (padding_time, 0), (init_time, 1), (singlet_decay + pi_half_time + tau_half_time + pi_time + tau_half_time + int(3 * pi_half_time), 0), (init_time, 1)
            ]

        # define sequence structure for DAQ counting triggers
        daq_clock_seq = [
            (padding_time + laser_lag, 0), 
            (clock_time, 1), (probe_time, 0), (clock_time, 1),
            (init_time - (2 * clock_time) - probe_time - laser_lag 
                + singlet_decay + pi_half_time + tau_half_time + pi_time + tau_half_time + pi_half_time + laser_lag, 0), 
            (clock_time, 1), (probe_time, 0), (clock_time, 1), 
            (init_time - (2 * clock_time) - probe_time - laser_lag 
                + padding_time + laser_lag, 0), 
            (clock_time, 1), (probe_time, 0), (clock_time, 1),
            (init_time - (2 * clock_time) - probe_time - laser_lag 
                + singlet_decay + pi_half_time + tau_half_time + pi_time + tau_half_time + int(3 * pi_half_time) + laser_lag, 0), 
            (clock_time, 1), (probe_time, 0), (clock_time, 1), 
            (init_time - (2 * clock_time) - probe_time - laser_lag, 0)
            ]

        # define sequence for MW (switch)
        switch_ttl_seq = [
            (padding_time + init_time + singlet_decay , 0), 
            (pi_half_time, 1), 
            (tau_half_time, 0), (pi_time, 1), (tau_half_time, 0), 
            (pi_half_time, 1), 
            (init_time + padding_time + init_time + singlet_decay, 0), 
            (pi_half_time, 1), 
            (tau_half_time, 0), (pi_time, 1), (tau_half_time, 0), 
            (int(3 * pi_half_time), 1), 
            (init_time, 0)]

        # assign sequences to respective channels for seq_on
        seq.setDigital(0, daq_clock_seq) 
        seq.setDigital(1, laser_seq) 
        seq.setDigital(4, switch_ttl_seq)

        return seq



    def BAYESIAN_T1(self, tau_time, pi_time):
        '''
        MW (differential) T1 sequence for the two longitudinal relaxation rates measured with different wait times
        '''
        # create sequence objects for MW on and off blocks
        seq = self.Pulser.createSequence()

        # need to create a laser buffer to avoid weirdly high counts when the laser first turns on
        laser_buffer_time = 5000

        # laser sequence
        laser_seq = [
            (init_time, 1), (singlet_decay, 0), 
            (init_time, 1), (singlet_decay + tau_time, 0), 
            (init_time, 1), (singlet_decay + pi_time, 0), 
            (init_time, 1), (singlet_decay + pi_time + tau_time, 0), 
            (init_time, 1), (laser_buffer_time, 0)]

        # define sequence structure for DAQ counting triggers
        daq_clock_seq = [
            (init_time + singlet_decay + laser_lag, 0), 
            (clock_time, 1), (probe_time, 0), (clock_time, 1), 
            (init_time - (2 * clock_time) - probe_time - laser_lag + singlet_decay + tau_time + laser_lag, 0), 
            (clock_time, 1), (probe_time, 0), (clock_time, 1), 
            (init_time - (2 * clock_time) - probe_time - laser_lag + singlet_decay + pi_time + laser_lag, 0), 
            (clock_time, 1), (probe_time, 0), (clock_time, 1), 
            (init_time - (2 * clock_time) - probe_time - laser_lag + singlet_decay + pi_time + tau_time + laser_lag, 0),
            (clock_time, 1), (probe_time, 0), (clock_time, 1), 
            (init_time - (2 * clock_time) - probe_time - laser_lag, 0),
            (laser_buffer_time, 0)            
            ]

        # define sequence for MW (switch)
        # NEED TO USE IQ FOR +/- BUT LET'S LET IT BE FOR NOW
        switch_ttl_seq = [
            (init_time + singlet_decay + init_time + singlet_decay + tau_time + init_time + singlet_decay , 0), 
            (pi_time, 1), 
            (init_time + singlet_decay, 0), 
            (pi_time, 1), 
            (tau_time + init_time, 0), (laser_buffer_time, 0)]

        # # define sequence structure for MW I and Q when MW = ON
        # mw_I_on_seq = [(iq_off1, self.IQ0[0]), self.Pi(pi_xy, pi_time)[0], (iq_off2, self.IQ0[0])]
        # mw_Q_on_seq = [(iq_off1, self.IQ0[1]), self.Pi(pi_xy, pi_time)[1], (iq_off2, self.IQ0[1])]
        # # when MW = OFF
        # mw_I_off_seq = [(iq_off1, self.IQ0[0]), (pi_time, self.IQ0[0]), (iq_off2, self.IQ0[0])]
        # mw_Q_off_seq = [(iq_off1, self.IQ0[1]), (pi_time, self.IQ0[1]), (iq_off2, self.IQ0[1])]

        # assign sequences to respective channels for seq_on
        seq.setDigital(0, daq_clock_seq) 
        seq.setDigital(1, laser_seq) 
        seq.setDigital(4, switch_ttl_seq)
        # seq.setAnalog(0, mw_I_on_seq) # mw_I
        # seq.setAnalog(1, mw_Q_on_seq) # mw_Q

        return seq


#################################################################################################################

    

    def Pulsed_ODMR(self, pi_xy, pi_time):
        '''
        Pulsed ODMR sequence with integrator
        '''
        ## Run a pi pulse, then measure the signal
        ## and reference counts from NV.
        pi_time = self.convert_type(round(pi_time), float)
       
        ## we can measure the pi time on x and on y.
        ## they should be the same, but they technically
        ## have different offsets on our pulse streamer.
        def Pi(axis):
            iq_on = pi_time 
            
            if axis == 'x':
                mw_I_on = (iq_on, self.IQpx[0])
                mw_Q_on = (iq_on, self.IQpx[1])
            else:
                mw_I_on = (iq_on, self.IQpy[0])
                mw_Q_on = (iq_on, self.IQpy[1])
            
            return mw_I_on, mw_Q_on

        def SinglePulsed_ODMR():
            '''
            CREATE SINGLE RABI SEQUENCE TO REPEAT THROUGHOUT EXPERIMENT
            '''

            '''
            DEFINE SPECIAL TIME INTERVALS FOR EXPERIMENT

            pad_time: padding time to equalize duration of every run (for different vsg_on durations)
            '''
            pad_time = 50000 - self.initial_delay - self.laser_time - self.singlet_decay - pi_time - self.MW_buffer_time - self.readout_time 

            '''
            DEFINE RELEVANT ON, OFF TIMES FOR DEVICES
            '''
            laser_off1 = self.initial_delay
            laser_off2 = self.singlet_decay + pi_time + self.MW_buffer_time
            laser_off3 = pad_time

            # integrator trigger windows     
            int_trig_off1 = laser_off1 + self.laser_time + (laser_off2 - self.trig_delay)
            int_trig_off2 = (self.trig_delay - self.clock_time) + self.readout_time + pad_time            

            # mw I & Q off windows
            iq_off1 = laser_off1 + self.laser_time + self.singlet_decay
            iq_off2 = self.MW_buffer_time + self.readout_time + pad_time

            '''
            CONSTRUCT PULSE SEQUENCE
            '''
            # create sequence objects for MW on and off blocks
            seq_on = self.Pulser.createSequence()
            seq_off = self.Pulser.createSequence()

            # define sequence structure for laser
            laser_seq = [(laser_off1, 0), (self.laser_time, 1), (laser_off2, 0), (self.readout_time, 1), (laser_off3, 0)]

            # define sequence structure for integrator trigger
            int_trig_seq = [(int_trig_off1, 0), (self.clock_time, 1), (int_trig_off2, 0)]
            
            # define sequence structure for MW I and Q when MW = ON
            mw_I_on_seq = [(iq_off1, self.IQ0[0]), Pi(pi_xy)[0], (iq_off2, self.IQ0[0])]
            mw_Q_on_seq = [(iq_off1, self.IQ0[1]), Pi(pi_xy)[1], (iq_off2, self.IQ0[1])]
            # when MW = OFF
            mw_I_off_seq = [(iq_off1, self.IQ0[0]), (pi_time, self.IQ0[0]), (iq_off2, self.IQ0[0])]
            mw_Q_off_seq = [(iq_off1, self.IQ0[1]), (pi_time, self.IQ0[1]), (iq_off2, self.IQ0[1])]

            # assign sequences to respective channels for seq_on
            seq_on.setDigital(3, laser_seq) # laser
            # seq_on.setDigital(4, int_trig_seq) # integrator trigger
            seq_on.setAnalog(0, mw_I_on_seq) # mw_I
            seq_on.setAnalog(1, mw_Q_on_seq) # mw_Q

            # assign sequences to respective channels for seq_off
            seq_off.setDigital(3, laser_seq) # laser
            # seq_off.setDigital(4, int_trig_seq) # integrator trigger
            seq_off.setAnalog(0, mw_I_off_seq) # mw_I
            seq_off.setAnalog(1, mw_Q_off_seq) # mw_Q

            return seq_on + seq_off

        seqs = self.Pulser.createSequence()

        for i in range(self.runs):
            seqs += SinglePulsed_ODMR()

        return seqs



    def Rabi_AWG(self, params, pi_xy):
        '''
        Rabi sequence
        '''
        ## Run a MW pulse of varying duration, then measure the signal
        ## and reference counts from NV.
        # self.total_time = 0
        longest_time = self.convert_type(round(params[-1]), float)
        ## we can measure the pi time on x and on y.
        ## they should be the same, but they technically
        ## have different offsets on our pulse streamer.
        if pi_xy == 'x':
            self.IQ_ON = self.IQpx
        else:
            self.IQ_ON = self.IQpy

        def SingleRabi(iq_on):
            '''
            CREATE SINGLE RABI SEQUENCE TO REPEAT THROUGHOUT EXPERIMENT
            '''

            iq_on = float(round(iq_on)) # convert to proper data type to avoid undesired rpyc netref data type

            '''
            DEFINE SPECIAL TIME INTERVALS FOR EXPERIMENT
            '''
            # padding time to equalize duration of every run (for different vsg_on durations)
            # pad_time = 50000 - self.initial_delay - self.laser_time - self.singlet_decay - iq_on - self.MW_buffer_time - self.readout_time 
            pad_time = longest_time - iq_on

            '''
            DEFINE RELEVANT ON, OFF TIMES FOR DEVICES
            '''

            laser_off1 = self.initial_delay 
            laser_off2 = self.singlet_decay + iq_on + self.MW_buffer_time
            laser_off3 = 100 + pad_time 
            # laser_off3 = pad_time + self.rest_time_btw_seqs
            # laser_off4 = laser_off2
            # laser_off5 = self.rest_time_btw_seqs

            # mw I & Q off windows
            iq_off1 = laser_off1 + self.laser_time + self.singlet_decay
            iq_off2 = self.MW_buffer_time + 1*self.readout_time + laser_off3 # + self.laser_time # + laser_off4 + laser_off5

            # DAQ trigger windows
            clock_off1 = laser_off1 + self.laser_time + laser_off2 + self.readout_time - self.trig_spot - self.clock_time
            clock_off2 = self.trig_spot + laser_off3
            
            print("SEQ TOTAL TIME = ", iq_off1 + iq_on + iq_off2)
            print("TIME AFTER MW PULSE = ", self.MW_buffer_time + laser_off3)
            '''
            CONSTRUCT PULSE SEQUENCE
            '''
            # create sequence objects for MW on and off blocks
            seq_on = self.Pulser.createSequence()
            seq_off = self.Pulser.createSequence()

            # define sequence structure for laser            
            laser_seq = [(laser_off1, 0), (self.laser_time, 1), (laser_off2, 0), (self.readout_time, 1), (laser_off3, 0)]
                        #  (laser_off3, 0), (self.laser_time, 1), (laser_off4, 0), (self.readout_time, 1), (laser_off5, 0)]
        
            # define sequence structure for DAQ trigger
            daq_clock_seq = [(clock_off1, 0), (self.clock_time, 1), (clock_off2, 0)]

            # define sequence structure for MW I and Q when MW = ON
            mw_I_on_seq = [(iq_off1, self.IQ0[0]), (iq_on, self.IQ_ON[0]), (iq_off2, self.IQ0[0])]
            mw_Q_on_seq = [(iq_off1, self.IQ0[1]), (iq_on, self.IQ_ON[1]), (iq_off2, self.IQ0[1])]
            
            # when MW = OFF
            mw_I_off_seq = [(iq_off1, self.IQ0[0]), (iq_on, self.IQ0[0]), (iq_off2, self.IQ0[0])]
            mw_Q_off_seq = [(iq_off1, self.IQ0[1]), (iq_on, self.IQ0[1]), (iq_off2, self.IQ0[1])]

            awg_seq = [(iq_off1, 0), (self.awg_trig_time, 1), (iq_off2 + iq_on - self.awg_trig_time, 0)]
            awg_ref_seq = [(iq_off1, 0), (self.awg_trig_time, 0), (iq_off2 + iq_on - self.awg_trig_time, 0)]
            # switch_on_seq = [(iq_off1 - 20, 0), (iq_on + 40, 1), (iq_off2 - 20, 0)]
            # switch_off_seq = [(iq_off1 - 20, 0), (iq_on + 40, 0), (iq_off2 - 20, 0)]

            # assign sequences to respective channels for seq_on
            seq_on.setDigital(3, laser_seq) # laser 
            seq_on.setDigital(0, daq_clock_seq) # integrator trigger
            seq_on.setDigital(4, awg_seq)
            # seq_on.setDigital(1, switch_on_seq) # RF control switch
            seq_on.setAnalog(0, mw_I_on_seq) # mw_I
            seq_on.setAnalog(1, mw_Q_on_seq) # mw_Q
            # seq_on.plot()
            print("LASER SEQ = ", laser_seq)
            print("DAQ SEQ = ", daq_clock_seq)
            print("AWG SEQ = ", awg_seq)

            # assign sequences to respective channels for seq_off
            seq_off.setDigital(3, laser_seq) # laser
            seq_off.setDigital(0, daq_clock_seq) # integrator trigger
            seq_off.setDigital(4, awg_ref_seq)
            # seq_off.setDigital(1, switch_off_seq) # RF control switch
            seq_off.setAnalog(0, mw_I_off_seq) # mw_I
            seq_off.setAnalog(1, mw_Q_off_seq) # mw_Q
            return seq_on + seq_off

        seqs = self.Pulser.createSequence()

        for mw_time in params:
            seqs += SingleRabi(mw_time)

        return seqs
    
    # def Rabi_AWG(self, params, pi_xy):
    #     '''
    #     Rabi sequence
    #     '''
    #     ## Run a MW pulse of varying duration, then measure the signal
    #     ## and reference counts from NV.
    #     # self.total_time = 0
    #     longest_time = self.convert_type(round(params[-1]), float)
    #     ## we can measure the pi time on x and on y.
    #     ## they should be the same, but they technically
    #     ## have different offsets on our pulse streamer.
    #     if pi_xy == 'x':
    #         self.IQ_ON = self.IQpx
    #     else:
    #         self.IQ_ON = self.IQpy

    #     def SingleRabi(iq_on):
    #         '''
    #         CREATE SINGLE RABI SEQUENCE TO REPEAT THROUGHOUT EXPERIMENT
    #         '''

    #         iq_on = float(round(iq_on)) # convert to proper data type to avoid undesired rpyc netref data type

    #         '''
    #         DEFINE SPECIAL TIME INTERVALS FOR EXPERIMENT
    #         '''
    #         # padding time to equalize duration of every run (for different vsg_on durations)
    #         # pad_time = 50000 - self.initial_delay - self.laser_time - self.singlet_decay - iq_on - self.MW_buffer_time - self.readout_time 
    #         pad_time = longest_time - iq_on

    #         '''
    #         DEFINE RELEVANT ON, OFF TIMES FOR DEVICES
    #         '''

    #         laser_off1 = self.initial_delay 
    #         laser_off2 = self.singlet_decay + iq_on + self.MW_buffer_time
    #         laser_off3 = 100 + pad_time 
    #         # laser_off3 = pad_time + self.rest_time_btw_seqs
    #         # laser_off4 = laser_off2
    #         # laser_off5 = self.rest_time_btw_seqs

    #         # mw I & Q off windows
    #         iq_off1 = laser_off1 + self.laser_time + self.singlet_decay
    #         iq_off2 = self.MW_buffer_time + 1*self.readout_time + laser_off3 # + self.laser_time # + laser_off4 + laser_off5

    #         # DAQ trigger windows
    #         clock_off1 = laser_off1 + self.laser_time + laser_off2 + self.readout_time - self.trig_spot - self.clock_time
    #         clock_off2 = self.trig_spot + laser_off3
            
    #         print("SEQ TOTAL TIME = ", iq_off1 + iq_on + iq_off2)
    #         print("TIME AFTER MW PULSE = ", self.MW_buffer_time + laser_off3)
    #         '''
    #         CONSTRUCT PULSE SEQUENCE
    #         '''
    #         # create sequence objects for MW on and off blocks
    #         seq_on = self.Pulser.createSequence()
    #         seq_off = self.Pulser.createSequence()

    #         # define sequence structure for laser            
    #         laser_seq = [(laser_off1, 0), (self.laser_time, 1), (laser_off2, 0), (self.readout_time, 1), (laser_off3, 0)]
    #                     #  (laser_off3, 0), (self.laser_time, 1), (laser_off4, 0), (self.readout_time, 1), (laser_off5, 0)]
        
    #         # define sequence structure for DAQ trigger
    #         daq_clock_seq = [(clock_off1, 0), (self.clock_time, 1), (clock_off2, 0)]

    #         # define sequence structure for MW I and Q when MW = ON
    #         mw_I_on_seq = [(iq_off1, self.IQ0[0]), (iq_on, self.IQ_ON[0]), (iq_off2, self.IQ0[0])]
    #         mw_Q_on_seq = [(iq_off1, self.IQ0[1]), (iq_on, self.IQ_ON[1]), (iq_off2, self.IQ0[1])]
            
    #         # when MW = OFF
    #         mw_I_off_seq = [(iq_off1, self.IQ0[0]), (iq_on, self.IQ0[0]), (iq_off2, self.IQ0[0])]
    #         mw_Q_off_seq = [(iq_off1, self.IQ0[1]), (iq_on, self.IQ0[1]), (iq_off2, self.IQ0[1])]

    #         awg_seq = [(iq_off1, 0), (self.awg_trig_time, 1), (iq_off2 + iq_on - self.awg_trig_time, 0)]
    #         awg_ref_seq = [(iq_off1, 0), (self.awg_trig_time, 0), (iq_off2 + iq_on - self.awg_trig_time, 0)]
    #         # switch_on_seq = [(iq_off1 - 20, 0), (iq_on + 40, 1), (iq_off2 - 20, 0)]
    #         # switch_off_seq = [(iq_off1 - 20, 0), (iq_on + 40, 0), (iq_off2 - 20, 0)]

    #         # assign sequences to respective channels for seq_on
    #         seq_on.setDigital(3, laser_seq) # laser 
    #         seq_on.setDigital(0, daq_clock_seq) # integrator trigger
    #         seq_on.setDigital(4, awg_seq)
    #         # seq_on.setDigital(1, switch_on_seq) # RF control switch
    #         seq_on.setAnalog(0, mw_I_on_seq) # mw_I
    #         seq_on.setAnalog(1, mw_Q_on_seq) # mw_Q
    #         # seq_on.plot()
    #         print("LASER SEQ = ", laser_seq)
    #         print("DAQ SEQ = ", daq_clock_seq)
    #         print("AWG SEQ = ", awg_seq)

    #         # assign sequences to respective channels for seq_off
    #         seq_off.setDigital(3, laser_seq) # laser
    #         seq_off.setDigital(0, daq_clock_seq) # integrator trigger
    #         seq_off.setDigital(4, awg_ref_seq)
    #         # seq_off.setDigital(1, switch_off_seq) # RF control switch
    #         seq_off.setAnalog(0, mw_I_off_seq) # mw_I
    #         seq_off.setAnalog(1, mw_Q_off_seq) # mw_Q
    #         return seq_on + seq_off

    #     seqs = self.Pulser.createSequence()

    #     for mw_time in params:
    #         seqs += SingleRabi(mw_time)

    #     return seqs


    def Diff_T1(self, params, pi_xy, pi_time):
        '''
        MW (differential) T1 sequence with integrator
        '''
        ## Run a pi pulse, then measure the signal
        ## and reference counts from NV.
        longest_time = self.convert_type(round(params[-1]), float)
        pi_time = self.convert_type(round(pi_time), float)

        ## we can measure the pi time on x and on y.
        ## they should be the same, but they technically
        ## have different offsets on our pulse streamer.

        def SingleDiff_T1(tau_time):
            '''
            CREATE SINGLE T1 SEQUENCE TO REPEAT THROUGHOUT EXPERIMENT
            '''

            tau_time = self.convert_type(round(tau_time), float) # convert to proper data type to avoid undesired rpyc netref data type

            '''
            DEFINE SPECIAL TIME INTERVALS FOR EXPERIMENT
            '''
            # padding time to equalize duration of every run
            pad_time = longest_time - tau_time 

            '''
            DEFINE RELEVANT ON, OFF TIMES FOR DEVICES
            '''
            laser_off1 = self.initial_delay
            laser_off2 = self.singlet_decay + pi_time + tau_time
            # laser_off3 = pad_time + self.rest_time_btw_seqs
            laser_off3 = 100 + pad_time

            # DAQ trigger windows
            clock_off1 = laser_off1 + self.laser_time + laser_off2 + self.readout_time - self.trig_spot - self.clock_time
            clock_off2 = self.trig_spot + laser_off3

            # mw I & Q off windows
            iq_off1 = laser_off1 + self.laser_time + self.singlet_decay
            iq_off2 = tau_time + self.readout_time + laser_off3

            '''
            CONSTRUCT PULSE SEQUENCE
            '''
            # create sequence objects for MW on and off blocks
            seq_on = self.Pulser.createSequence()
            seq_off = self.Pulser.createSequence()

            # define sequence structure for laser
            laser_seq = [(laser_off1, 0), (self.laser_time, 1), (laser_off2, 0), (self.readout_time, 1), (laser_off3, 0)]

            # define sequence structure for DAQ trigger
            daq_clock_seq = [(clock_off1, 0), (self.clock_time, 1), (clock_off2, 0)]

            # define sequence structure for MW I and Q when MW = ON
            mw_I_on_seq = [(iq_off1, self.IQ0[0]), self.Pi(pi_xy, pi_time)[0], (iq_off2, self.IQ0[0])]
            mw_Q_on_seq = [(iq_off1, self.IQ0[1]), self.Pi(pi_xy, pi_time)[1], (iq_off2, self.IQ0[1])]
            # when MW = OFF
            mw_I_off_seq = [(iq_off1, self.IQ0[0]), (pi_time, self.IQ0[0]), (iq_off2, self.IQ0[0])]
            mw_Q_off_seq = [(iq_off1, self.IQ0[1]), (pi_time, self.IQ0[1]), (iq_off2, self.IQ0[1])]

            # assign sequences to respective channels for seq_on
            seq_on.setDigital(3, laser_seq) # laser
            seq_on.setDigital(0, daq_clock_seq) # integrator trigger
            seq_on.setAnalog(0, mw_I_on_seq) # mw_I
            seq_on.setAnalog(1, mw_Q_on_seq) # mw_Q

            # assign sequences to respective channels for seq_off
            seq_off.setDigital(3, laser_seq) # laser
            seq_off.setDigital(0, daq_clock_seq) # integrator trigger
            seq_off.setAnalog(0, mw_I_off_seq) # mw_I
            seq_off.setAnalog(1, mw_Q_off_seq) # mw_Q

            return seq_on + seq_off

        seqs = self.Pulser.createSequence()

        for tau in params:
            seqs += SingleDiff_T1(tau)

        return seqs
    

    def Calibrate_LaserLag(self, params, buffer_time):
        
        longest_time = self.convert_type(round(params[-1]), float)
        buffer_time = self.convert_type(round(buffer_time), float)

        def SingleLag(read_time):
            '''
            CREATE SINGLE T1 SEQUENCE TO REPEAT THROUGHOUT EXPERIMENT
            '''

            read_time = int(round(read_time)) # convert to proper data type to avoid undesired rpyc netref data type

            '''
            DEFINE SPECIAL TIME INTERVALS FOR EXPERIMENT
            '''
            # padding time to equalize duration of every run
            pad_time = longest_time - read_time 

            '''
            DEFINE RELEVANT ON, OFF TIMES FOR DEVICES
            '''
            laser_off1 = self.initial_delay + buffer_time + self.trig_delay
            laser_off2 = buffer_time + self.rest_time_btw_seqs

            # integrator trigger windows     
            int_trig_off1 = self.initial_delay + read_time
            
            int_trig_off2 = (self.trig_delay - self.clock_time) + pad_time + self.rest_time_btw_seqs           

            '''
            CONSTRUCT PULSE SEQUENCE
            '''
            # create sequence objects for MW on and off blocks
            seq1 = self.Pulser.createSequence()
            seq2 = self.Pulser.createSequence()

            # define sequence structure for laser
            laser_seq = [(laser_off1, 0), (self.laser_time, 1), (laser_off2, 0)]

            # define sequence structure for integrator trigger
            int_trig_seq = [(int_trig_off1, 0), (self.clock_time, 1), (int_trig_off2, 0)]

            # assign sequences to respective channels for seq_on
            seq1.setDigital(7, laser_seq) # laser
            seq1.setDigital(4, int_trig_seq) # integrator trigger
            seq2.setDigital(7, laser_seq) # laser
            seq2.setDigital(4, int_trig_seq) # integrator trigger

            return seq1 + seq2 + seq2 + seq1

        seqs = self.Pulser.createSequence()

        for read in params:
            seqs += SingleLag(read)

        return seqs


    def Calibrate_Initialize(self, params, init_pulse_length):
        
        longest_time = self.convert_type(round(params[-1]), float)
        init_pulse_length = self.convert_type(round(init_pulse_length), float)
        self.initialize = init_pulse_length

        def SingleInitialize(read_time):
            '''
            CREATE SINGLE T1 SEQUENCE TO REPEAT THROUGHOUT EXPERIMENT
            '''

            read_time = int(round(read_time)) # convert to proper data type to avoid undesired rpyc netref data type

            '''
            DEFINE SPECIAL TIME INTERVALS FOR EXPERIMENT
            '''
            # padding time to equalize duration of every run
            pad_time = longest_time - read_time 

            '''
            DEFINE RELEVANT ON, OFF TIMES FOR DEVICES
            '''
            laser_off1 = self.initial_delay
            laser_off2 = pad_time + self.rest_time_btw_seqs

            # integrator trigger windows     
            int_trig_off1 = laser_off1 + (read_time - self.trig_delay)
            
            int_trig_off2 = (self.trig_delay - self.clock_time) + self.initialize + laser_off2            

            '''
            CONSTRUCT PULSE SEQUENCE
            '''
            # create sequence objects for MW on and off blocks
            seq1 = self.Pulser.createSequence()
            seq2 = self.Pulser.createSequence()

            # define sequence structure for laser
            laser_seq = [(laser_off1, 0), (self.initialize, 1), (laser_off2, 0)]

            # define sequence structure for integrator trigger
            int_trig_seq = [(int_trig_off1, 0), (self.clock_time, 1), (int_trig_off2, 0)]

            print("LASER SEQ: ", laser_seq)
            print("TRIG SEQ: ", int_trig_seq)

            # assign sequences to respective channels for seq_on
            seq1.setDigital(7, laser_seq) # laser
            seq1.setDigital(4, int_trig_seq) # integrator trigger
            seq2.setDigital(7, laser_seq) # laser
            seq2.setDigital(4, int_trig_seq) # integrator trigger

            return seq1 + seq2 + seq2 + seq1

        seqs = self.Pulser.createSequence()

        for read in params:
            seqs += SingleInitialize(read)

        return seqs

    def Calibrate_SingletDecay(self, params):
        longest_time = self.convert_type(round(params[-1]), float)

        def SingleDecay(decay_time):
            '''
            CREATE SINGLE T1 SEQUENCE TO REPEAT THROUGHOUT EXPERIMENT
            '''

            decay_time = int(round(decay_time)) # convert to proper data type to avoid undesired rpyc netref data type

            '''
            DEFINE SPECIAL TIME INTERVALS FOR EXPERIMENT
            '''
            # padding time to equalize duration of every run
            pad_time = longest_time - decay_time 

            '''
            DEFINE RELEVANT ON, OFF TIMES FOR DEVICES
            '''
            laser_off1 = self.initial_delay
            laser_off2 = decay_time
            laser_off3 = pad_time + self.rest_time_btw_seqs

            # integrator trigger windows     
            int_trig_off1 = laser_off1 + self.laser_time + (laser_off2 - self.trig_delay)
            
            int_trig_off2 = (self.trig_delay - self.clock_time) + self.readout_time + laser_off3            

            '''
            CONSTRUCT PULSE SEQUENCE
            '''
            # create sequence objects for MW on and off blocks
            seq1 = self.Pulser.createSequence()
            seq2 = self.Pulser.createSequence()

            # define sequence structure for laser
            laser_seq = [(laser_off1, 0), (self.laser_time, 1), (laser_off2, 0), (self.readout_time, 1), (laser_off3, 0)]
            
            # define sequence structure for integrator trigger
            int_trig_seq = [(int_trig_off1, 0), (self.clock_time, 1), (int_trig_off2, 0)]

            # assign sequences to respective channels for seq_on
            seq1.setDigital(7, laser_seq) # laser
            seq1.setDigital(4, int_trig_seq) # integrator trigger
            seq2.setDigital(7, laser_seq) # laser
            seq2.setDigital(4, int_trig_seq) # integrator trigger

            return seq1 + seq2 + seq2 + seq1

        seqs = self.Pulser.createSequence()

        for decay in params:
            seqs += SingleDecay(decay)

        return seqs

    def Calibrate_IntSNR(self, pi, pi_xy):
        '''
        Calibrate integration window duration for experiments w/ integrator.
        Assess SNR for these integration windows.
        
        Run a pi pulse, then measure the signal and reference counts from NV.
        '''

        pi_ns = pi.to("ns").magnitude
        iq_buffer = 15
        iq_on = 2*iq_buffer + pi_ns

        def Pi(axis):
            vsg_on = (pi_ns, 1)

            if axis == 'x':
                mw_I_on = (iq_on, self.IQpx[0])
                mw_Q_on = (iq_on, self.IQpx[1])
            elif axis == '-x':
                mw_I_on = (iq_on, self.IQnx[0])
                mw_Q_on = (iq_on, self.IQnx[1])
            elif axis == 'y':
                mw_I_on = (iq_on, self.IQpy[0])
                mw_Q_on = (iq_on, self.IQpy[1])
            elif axis == '-y':
                mw_I_on = (iq_on, self.IQny[0])
                mw_Q_on = (iq_on, self.IQny[1])
            
            return vsg_on, mw_I_on, mw_Q_on

        def SingleSNR():
            '''
            CREATE SINGLE SNR (RABI-LIKE) SEQUENCE TO REPEAT THROUGHOUT EXPERIMENT
            '''

            '''
            DEFINE SPECIAL TIME INTERVALS FOR EXPERIMENT
            '''
            # padding time to equalize duration of every run (for different vsg_on durations)
            pad_time = 50000 - self.laser_time - self.MW_buffer_time - 2*iq_buffer - pi_ns

            '''
            DEFINE RELEVANT ON, OFF TIMES FOR DEVICES
            '''
            # mw I & Q off windows (on slightly longer than VSG to ensure it's set)
            iq_off1 = self.laser_time + self.singlet_decay
            iq_off2 = self.MW_buffer_time + self.laser_time + pad_time

            # VSG disable windows
            vsg_off1 = self.laser_time + self.singlet_decay + iq_buffer
            vsg_off2 = iq_buffer + iq_off2

            laser_off1 = self.singlet_decay + iq_on + self.MW_buffer_time
            laser_off2 = pad_time

            # integrator trigger windows
            int_trig_off1 = self.laser_time + (laser_off1 - self.trig_delay)
            int_trig_off2 = (self.trig_delay - self.clock_time) + self.laser_time + pad_time            


            '''
            CONSTRUCT PULSE SEQUENCE
            '''
            # create sequence objects for MW on and off blocks
            seq_on = self.Pulser.createSequence()
            seq_off = self.Pulser.createSequence()

            # define sequence structure for laser
            laser_seq = [(self.laser_time, 1), (laser_off1, 0), (self.laser_time, 1), (laser_off2, 0)]
            
            # define sequence structure for integrator trigger
            int_trig_seq = [(int_trig_off1, 0), (self.clock_time, 1), (int_trig_off2, 0)]
            
            # define sequence structure for MW I and Q when MW = ON
            mw_I_on_seq = [(vsg_off1, self.IQ0[0]), (pi_ns, self.IQpx[0]), (vsg_off2, self.IQ0[0])]
            mw_Q_on_seq = [(vsg_off1, self.IQ0[1]), (pi_ns, self.IQpx[1]), (vsg_off2, self.IQ0[1])]
            
            # when MW = OFF
            mw_I_off_seq = [(vsg_off1, self.IQ0[0]), (pi_ns, self.IQ0[0]), (vsg_off2, self.IQ0[0])]
            mw_Q_off_seq = [(vsg_off1, self.IQ0[1]), (pi_ns, self.IQ0[1]), (vsg_off2, self.IQ0[1])]

            # assign sequences to respective channels for seq_on
            seq_on.setDigital(7, laser_seq) # laser
            seq_on.setDigital(4, int_trig_seq) # integrator trigger
            seq_on.setAnalog(0, mw_I_on_seq) # mw_I
            seq_on.setAnalog(1, mw_Q_on_seq) # mw_Q
            
            # assign sequences to respective channels for seq_off
            seq_off.setDigital(7, laser_seq) # laser
            seq_off.setDigital(4, int_trig_seq) # integrator trigger
            seq_off.setAnalog(0, mw_I_off_seq) # mw_I
            seq_off.setAnalog(1, mw_Q_off_seq) # mw_Q

            return seq_on + seq_off + seq_off + seq_on # sequence order to deal with int. offset

        seqs = SingleSNR()

        return seqs


    def Calibrate_Switch_Echo(self, tau_set, pihalf_x, pihalf_y, pi_x, pi_y):
        '''
        Spin Echo pulse sequence for calibrating RF switch.
        MW sequence: pi/2(x) - tau - pi(y) - tau - pi/2(x)

        '''
        tau_time = self.convert_type(round(tau_set), float)
        pihalf_x = self.convert_type(round(pihalf_x), float)
        pihalf_y = self.convert_type(round(pihalf_y), float)
        pi_x = self.convert_type(round(pi_x), float)
        pi_y = self.convert_type(round(pi_y), float)

        def SingleCalEcho():
            '''
            CREATE SINGLE HAHN-ECHO SEQUENCE TO REPEAT THROUGHOUT EXPERIMENT
            '''

            '''
            DEFINE RELEVANT ON, OFF TIMES FOR DEVICES
            '''            

            # mw I & Q off windows (on slightly longer than VSG to ensure it's set)
            iq_off1 = 200
            iq_off2 = tau_time # /2
            iq_off3 = tau_time # /2
            iq_off4 = 200

            '''
            CONSTRUCT PULSE SEQUENCE
            '''
            # create sequence objects for MW on and off blocks
            seq = self.Pulser.createSequence()
            seq_ref = self.Pulser.createSequence()
            
            # sequence structure for I & Q MW channels 
            mw_I_seq = [(iq_off1, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0], (iq_off2, self.IQ0[0]), self.Pi('y', pi_y)[0], (iq_off3, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0], (iq_off4, self.IQ0[0])]
            mw_Q_seq = [(iq_off1, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1], (iq_off2, self.IQ0[1]), self.Pi('y', pi_y)[1], (iq_off3, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1], (iq_off4, self.IQ0[1])]

            # sequence structure for I & Q MW channels (MW off)
            mw_I_ref_seq = [(iq_off1, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0], (iq_off2, self.IQ0[0]), self.Pi('y', pi_y)[0], (iq_off3, self.IQ0[0]), self.PiHalf('-x', pihalf_x)[0], (iq_off4, self.IQ0[0])]
            mw_Q_ref_seq = [(iq_off1, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1], (iq_off2, self.IQ0[1]), self.Pi('y', pi_y)[1], (iq_off3, self.IQ0[1]), self.PiHalf('-x', pihalf_x)[1], (iq_off4, self.IQ0[1])]

            # for testing switch, keep MW I and Q constantly on
            mw_I_TEST_seq = [(iq_off1, self.IQpx[0]), self.PiHalf('x', pihalf_x)[0], (iq_off2, self.IQpx[0]), self.Pi('y', pi_y)[0], (iq_off3, self.IQpx[0]), self.PiHalf('x', pihalf_x)[0], (iq_off4, self.IQpx[0])]
            mw_Q_TEST_seq = [(iq_off1, self.IQpx[1]), self.PiHalf('x', pihalf_x)[1], (iq_off2, self.IQpx[1]), self.Pi('y', pi_y)[1], (iq_off3, self.IQpx[1]), self.PiHalf('x', pihalf_x)[1], (iq_off4, self.IQpx[1])]

            switch_on_seq = [(iq_off1, 0), (pihalf_x, 1), (iq_off2, 0), (pi_y, 1), (iq_off3, 0), (pihalf_x, 1), (iq_off4, 0)]
            switch_off_seq = [(iq_off1, 0), (pihalf_x, 0), (iq_off2, 0), (pi_y, 0), (iq_off3, 0), (pihalf_x, 0), (iq_off4, 0)]

            # switch_on_seq = [(iq_off1 - 20, 0), (pihalf_x + 40, 1), (iq_off2 - 40, 0), (pi_y + 40, 1), (iq_off3 - 40, 0), (pihalf_x + 40, 1), (iq_off4 - 20, 0)]
            # switch_off_seq = [(iq_off1 - 20, 0), (pihalf_x + 40, 0), (iq_off2 - 40, 0), (pi_y + 40, 0), (iq_off3 - 40, 0), (pihalf_x + 40, 0), (iq_off4 - 20, 0)]

            # assign sequences to respective channels for seq_on
            # seq.setDigital(3, laser_seq) # laser
            # seq.setDigital(0, daq_clock_seq) # integrator trigger
            seq.setDigital(5, switch_on_seq) # RF switch 
            seq.setDigital(2, switch_on_seq) # RF switch to oscilloscope
            # seq.setAnalog(0, mw_I_seq) # mw_I
            # seq.setAnalog(1, mw_Q_seq) # mw_Q
            seq.setAnalog(0, mw_I_TEST_seq) # mw_I
            seq.setAnalog(1, mw_Q_TEST_seq) # mw_Q
            
            # assign sequences to respective channels for seq_off
            # seq_ref.setDigital(3, laser_seq) # laser
            # seq_ref.setDigital(0, daq_clock_seq) # integrator trigger
            seq_ref.setDigital(5, switch_off_seq) # RF switch 
            seq_ref.setDigital(2, switch_off_seq) # RF switch to oscilloscope
            # seq_ref.setAnalog(0, mw_I_ref_seq) # mw_I
            # seq_ref.setAnalog(1, mw_Q_ref_seq) # mw_Q
            seq_ref.setAnalog(0, mw_I_TEST_seq) # mw_I
            seq_ref.setAnalog(1, mw_Q_TEST_seq) # mw_Q

            return seq # + seq_ref 

        seqs = SingleCalEcho()

        return seqs
    

    def Calibrate_DEER_offset(self):

        def SingleCalDEER():
            '''
            CREATE SINGLE DEER SEQUENCE TO REPEAT THROUGHOUT EXPERIMENT
            '''
            # mw I & Q off windows (on slightly longer than VSG to ensure it's set)
            iq_off1 = 0
            iq_on = 200
            
            awg_off = 1000
            awg_trig = 10
            '''
            CONSTRUCT PULSE SEQUENCE
            '''
            # create sequence objects for MW on and off blocks
            seq = self.Pulser.createSequence()
            
            # sequence structure for I & Q MW channels on SRS SG396
            # srs_I_seq = [(iq_off1, self.IQpx[0]), (iq_on, self.IQpx[0]), (1000, self.IQpx[0])]
            # srs_Q_seq = [(iq_off1, self.IQpx[1]), (iq_on, self.IQpx[1]), (1000, self.IQpx[1])]

            srs_I_seq = [(iq_off1, self.IQ0[0]), (iq_on, 1), (10, self.IQ0[0])]
            srs_Q_seq = [(iq_off1, self.IQ0[1]), (iq_on, 0), (10, self.IQ0[1])]

            # sequence structure for I & Q MW channels (MW off)
            awg_seq = [(0, 0), (awg_trig, 1), (10, 0)]

            # print("SRS I: ", srs_I_seq)
            print("AWG: ", awg_seq)
            # assign sequences to respective channels for seq_on
            seq.setDigital(4, awg_seq)
            seq.setAnalog(0, srs_I_seq) # mw_I
            seq.setAnalog(1, srs_Q_seq) # mw_Q

            return seq
        
        # concatenate single ODMR sequence "runs" number of times

        seqs = SingleCalDEER()

        return seqs
 
    


    def Ramsey(self, params, pihalf_x, pihalf_y):
        
        '''
        Ramsey pulse sequence.
        MW sequence: pi/2(x) - tau - pi/2(x)

        '''
        longest_time = self.convert_type(round(params[-1]), int)
        pihalf_x = self.convert_type(round(pihalf_x), int)
        pihalf_y = self.convert_type(round(pihalf_y), int)

        def SingleRamsey(tau_time):
            '''
            CREATE SINGLE RAMSEY SEQUENCE TO REPEAT THROUGHOUT EXPERIMENT
            '''

            tau_time = self.convert_type(round(tau_time), int) # convert to proper data type to avoid undesired rpyc netref data type

            '''
            DEFINE SPECIAL TIME INTERVALS FOR EXPERIMENT
            '''
            # padding time to equalize duration of every run (for different tau durations)
            pad_time = longest_time - tau_time 

            '''
            DEFINE RELEVANT ON, OFF TIMES FOR DEVICES
            '''

            laser_off1 = self.initial_delay
            laser_off2 = self.singlet_decay + pihalf_x + tau_time + pihalf_x + self.MW_buffer_time
            laser_off3 = 100 + pad_time 
            
            # DAQ trigger windows
            clock_off1 = laser_off1 + self.laser_time + laser_off2 + self.readout_time - self.trig_spot - self.clock_time
            clock_off2 = self.trig_spot + laser_off3          

            # mw I & Q off windows (on slightly longer than VSG to ensure it's set)
            iq_off1 = laser_off1 + self.laser_time + self.singlet_decay
            iq_off2 = tau_time
            iq_off3 = self.MW_buffer_time + self.readout_time + laser_off3

            '''
            CONSTRUCT PULSE SEQUENCE
            '''
            # create sequence objects for MW on and off blocks
            seq = self.Pulser.createSequence()
            seq_ref = self.Pulser.createSequence()

            # define sequence structure for laser
            laser_seq = [(laser_off1, 0), (self.laser_time, 1), (laser_off2, 0), (self.readout_time, 1), (laser_off3, 0)]
            
            # define sequence structure for integrator trigger
            daq_clock_seq = [(clock_off1, 0), (self.clock_time, 1), (clock_off2, 0)]
            
            # sequence structure for I & Q MW channels 
            mw_I_seq = [(iq_off1, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0], (iq_off2, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0], (iq_off3, self.IQ0[0])]
            mw_Q_seq = [(iq_off1, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1], (iq_off2, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1], (iq_off3, self.IQ0[1])]

            # sequence structure for I & Q MW channels (MW off)
            mw_I_ref_seq = [(iq_off1, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0], (iq_off2, self.IQ0[0]), self.PiHalf('-x', pihalf_x)[0], (iq_off3, self.IQ0[0])]
            mw_Q_ref_seq = [(iq_off1, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1], (iq_off2, self.IQ0[1]), self.PiHalf('-x', pihalf_x)[1], (iq_off3, self.IQ0[1])]

            # assign sequences to respective channels for seq_on
            seq.setDigital(3, laser_seq) # laser
            seq.setDigital(0, daq_clock_seq) # integrator trigger
            seq.setAnalog(0, mw_I_seq) # mw_I
            seq.setAnalog(1, mw_Q_seq) # mw_Q
            
            # assign sequences to respective channels for seq_off
            seq_ref.setDigital(3, laser_seq) # laser
            seq_ref.setDigital(0, daq_clock_seq) # integrator trigger
            seq_ref.setAnalog(0, mw_I_ref_seq) # mw_I
            seq_ref.setAnalog(1, mw_Q_ref_seq) # mw_Q

            return seq + seq_ref

        # concatenate single ODMR sequence "runs" number of times
        seqs = self.Pulser.createSequence()
        
        for tau in params:
            seqs += SingleRamsey(tau)
        
        return seqs


    def Echo(self, params, pihalf_x, pihalf_y, pi_x, pi_y):
        '''
        Spin Echo pulse sequence.
        MW sequence: pi/2(x) - tau - pi(y) - tau - pi/2(x)

        '''
        longest_time = self.convert_type(round(params[-1]), float)
        pihalf_x = self.convert_type(round(pihalf_x), float)
        pihalf_y = self.convert_type(round(pihalf_y), float)
        pi_x = self.convert_type(round(pi_x), float)
        pi_y = self.convert_type(round(pi_y), float)

        def SingleEcho(tau_time):
            '''
            CREATE SINGLE HAHN-ECHO SEQUENCE TO REPEAT THROUGHOUT EXPERIMENT
            '''

            tau_time = self.convert_type(round(tau_time), float) # convert to proper data type to avoid undesired rpyc netref data type
            
            '''
            DEFINE SPECIAL TIME INTERVALS FOR EXPERIMENT
            '''
            # padding time to equalize duration of every run (for different tau durations)
            pad_time = longest_time - tau_time

            '''
            DEFINE RELEVANT ON, OFF TIMES FOR DEVICES
            '''            
            laser_off1 = self.initial_delay
            # laser_off2 = self.singlet_decay + pihalf_x + tau_time/2 + pi_y + tau_time/2 + pihalf_x + self.MW_buffer_time
            laser_off2 = self.singlet_decay + pihalf_x + tau_time + pi_y + tau_time + pihalf_x + self.MW_buffer_time
            laser_off3 = 100 + pad_time 

            # DAQ trigger windows
            clock_off1 = laser_off1 + self.laser_time + laser_off2 + self.readout_time - self.trig_spot - self.clock_time
            clock_off2 = self.trig_spot + laser_off3

            # mw I & Q off windows (on slightly longer than VSG to ensure it's set)
            iq_off1 = laser_off1 + self.laser_time + self.singlet_decay
            iq_off2 = tau_time # /2
            iq_off3 = tau_time # /2
            iq_off4 = self.MW_buffer_time + self.readout_time + laser_off3

            '''
            CONSTRUCT PULSE SEQUENCE
            '''
            # create sequence objects for MW on and off blocks
            seq = self.Pulser.createSequence()
            seq_ref = self.Pulser.createSequence()

            # define sequence structure for laser
            laser_seq = [(laser_off1, 0), (self.laser_time, 1), (laser_off2, 0), (self.readout_time, 1), (laser_off3, 0)]
            
            # define sequence structure for integrator trigger
            daq_clock_seq = [(clock_off1, 0), (self.clock_time, 1), (clock_off2, 0)]
            
            # sequence structure for I & Q MW channels 
            mw_I_seq = [(iq_off1, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0], (iq_off2, self.IQ0[0]), self.Pi('y', pi_y)[0], (iq_off3, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0], (iq_off4, self.IQ0[0])]
            mw_Q_seq = [(iq_off1, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1], (iq_off2, self.IQ0[1]), self.Pi('y', pi_y)[1], (iq_off3, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1], (iq_off4, self.IQ0[1])]

            # sequence structure for I & Q MW channels (MW off)
            mw_I_ref_seq = [(iq_off1, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0], (iq_off2, self.IQ0[0]), self.Pi('y', pi_y)[0], (iq_off3, self.IQ0[0]), self.PiHalf('-x', pihalf_x)[0], (iq_off4, self.IQ0[0])]
            mw_Q_ref_seq = [(iq_off1, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1], (iq_off2, self.IQ0[1]), self.Pi('y', pi_y)[1], (iq_off3, self.IQ0[1]), self.PiHalf('-x', pihalf_x)[1], (iq_off4, self.IQ0[1])]

            # for testing switch, keep MW I and Q constantly on
            # mw_I_TEST_seq = [(iq_off1, self.IQpx[0]), self.PiHalf('x', pihalf_x)[0], (iq_off2, self.IQpx[0]), self.Pi('y', pi_y)[0], (iq_off3, self.IQpx[0]), self.PiHalf('x', pihalf_x)[0], (iq_off4, self.IQpx[0])]
            # mw_Q_TEST_seq = [(iq_off1, self.IQpx[1]), self.PiHalf('x', pihalf_x)[1], (iq_off2, self.IQpx[1]), self.Pi('y', pi_y)[1], (iq_off3, self.IQpx[1]), self.PiHalf('x', pihalf_x)[1], (iq_off4, self.IQpx[1])]

            # switch_on_seq = [(iq_off1, 0), (pihalf_x, 1), (iq_off2, 0), (pi_y, 1), (iq_off3, 0), (pihalf_x, 1), (iq_off4, 0)]
            # switch_off_seq = [(iq_off1, 0), (pihalf_x, 0), (iq_off2, 0), (pi_y, 0), (iq_off3, 0), (pihalf_x, 0), (iq_off4, 0)]

            # switch_on_seq = [(iq_off1 - 20, 0), (pihalf_x + 40, 1), (iq_off2 - 40, 0), (pi_y + 40, 1), (iq_off3 - 40, 0), (pihalf_x + 40, 1), (iq_off4 - 20, 0)]
            # switch_off_seq = [(iq_off1 - 20, 0), (pihalf_x + 40, 0), (iq_off2 - 40, 0), (pi_y + 40, 0), (iq_off3 - 40, 0), (pihalf_x + 40, 0), (iq_off4 - 20, 0)]

            # assign sequences to respective channels for seq_on
            seq.setDigital(3, laser_seq) # laser
            seq.setDigital(0, daq_clock_seq) # integrator trigger
            # seq.setDigital(5, switch_on_seq) # RF switch 
            # seq.setDigital(2, switch_on_seq) # RF switch to oscilloscope
            seq.setAnalog(0, mw_I_seq) # mw_I
            seq.setAnalog(1, mw_Q_seq) # mw_Q
            # seq.setAnalog(0, mw_I_TEST_seq) # mw_I
            # seq.setAnalog(1, mw_Q_TEST_seq) # mw_Q
            
            # assign sequences to respective channels for seq_off
            seq_ref.setDigital(3, laser_seq) # laser
            seq_ref.setDigital(0, daq_clock_seq) # integrator trigger
            # seq_ref.setDigital(5, switch_off_seq) # RF switch 
            # seq_ref.setDigital(2, switch_off_seq) # RF switch to oscilloscope
            seq_ref.setAnalog(0, mw_I_ref_seq) # mw_I
            seq_ref.setAnalog(1, mw_Q_ref_seq) # mw_Q
            # seq_ref.setAnalog(0, mw_I_TEST_seq) # mw_I
            # seq_ref.setAnalog(1, mw_Q_TEST_seq) # mw_Q

            return seq + seq_ref 
        
        # concatenate single ODMR sequence "runs" number of times
        seqs = self.Pulser.createSequence()
        
        for tau in params:
            seqs += SingleEcho(tau)

        return seqs

    def WAHUHA(self, params, pihalf_x, pihalf_y, pi_x, pi_y):
        '''
        WAHUHA pulse sequence applied to NV centers.
        MW sequence: pi/2(x) - tau/2 - (pi/2(x) - tau - pi/2(-y) - tau - pi/2(y) - tau - pi/2(-x))^N - tau/2 - pi/2(x)

        '''
        longest_time = self.convert_type(round(params[-1]), float)
        pihalf_x = self.convert_type(round(pihalf_x), float)
        pihalf_y = self.convert_type(round(pihalf_y), float)
        
        def PiPulsesN(axes, tau, N):  
            wahuha_I_seq = [self.PiHalf('x', pihalf_x)[0], (iq_off2, self.IQ0[0]), self.PiHalf('-y', pihalf_y)[0], (iq_off3, self.IQ0[0]), 
                        self.PiHalf('y', pihalf_y)[0], (iq_off4, self.IQ0[0]), self.PiHalf('-x', pihalf_x)[0]]
            wahuha_Q_seq = [self.PiHalf('x', pihalf_x)[1], (iq_off2, self.IQ0[1]), self.PiHalf('-y', pihalf_y)[1], (iq_off3, self.IQ0[1]), 
                        self.PiHalf('y', pihalf_y)[1], (iq_off4, self.IQ0[1]), self.PiHalf('-x', pihalf_x)[1]]
            
            tau_half_I_seq = [((tau/2), self.IQ0[0])]
            tau_half_Q_seq = [((tau/2), self.IQ0[1])]
            tau_I_seq = [(tau, self.IQ0[0])]
            tau_Q_seq = [(tau, self.IQ0[1])]

            xy4_I_seq = [self.Pi('x', pi_x)[0], (tau, self.IQ0[0]), self.Pi('y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('x', pi_x)[0], (tau, self.IQ0[0]), self.Pi('y', pi_y)[0]]
            xy4_Q_seq = [self.Pi('x', pi_x)[1], (tau, self.IQ0[1]), self.Pi('y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('x', pi_x)[1], (tau, self.IQ0[1]), self.Pi('y', pi_y)[1]]
            mw_I = (tau_half_I_seq + xy4_I_seq + tau_I_seq + list(reversed(xy4_I_seq)) + tau_half_I_seq)*N
            mw_Q = (tau_half_Q_seq + xy4_Q_seq + tau_Q_seq + list(reversed(xy4_Q_seq)) + tau_half_Q_seq)*N

            return mw_I, mw_Q
        

        def SingleWAHUHA(tau_time):
            '''
            CREATE SINGLE WAHUHA SEQUENCE TO REPEAT THROUGHOUT EXPERIMENT
            '''

            tau_time = self.convert_type(round(tau_time), float) # convert to proper data type to avoid undesired rpyc netref data type
            
            '''
            DEFINE SPECIAL TIME INTERVALS FOR EXPERIMENT
            '''
            # padding time to equalize duration of every run (for different tau durations)
            pad_time = longest_time - tau_time

            '''
            DEFINE RELEVANT ON, OFF TIMES FOR DEVICES
            '''            
            laser_off1 = self.initial_delay
            # laser_off2 = self.singlet_decay + pihalf_x + tau_time/2 + pi_y + tau_time/2 + pihalf_x + self.MW_buffer_time
            laser_off2 = self.singlet_decay + pihalf_x + tau_time + pihalf_y + tau_time + pihalf_y + tau_time + pihalf_x + self.MW_buffer_time
            laser_off3 = 100 + pad_time
            
            # DAQ trigger windows
            clock_off1 = laser_off1 + self.laser_time + laser_off2 + self.readout_time - self.trig_spot - self.clock_time
            clock_off2 = self.trig_spot + laser_off3

            # mw I & Q off windows (on slightly longer than VSG to ensure it's set)
            iq_off1 = laser_off1 + self.laser_time + self.singlet_decay
            iq_off2 = tau_time # /2
            iq_off3 = 2*tau_time # /2
            iq_off4 = tau_time # /2
            iq_off5 = self.MW_buffer_time + self.readout_time + laser_off3

            '''
            CONSTRUCT PULSE SEQUENCE
            '''
            # create sequence objects for MW on and off blocks
            seq = self.Pulser.createSequence()
            seq_ref = self.Pulser.createSequence()

            # define sequence structure for laser
            laser_seq = [(laser_off1, 0), (self.laser_time, 1), (laser_off2, 0), (self.readout_time, 1), (laser_off3, 0)]
            
            # define sequence structure for integrator trigger
            daq_clock_seq = [(clock_off1, 0), (self.clock_time, 1), (clock_off2, 0)]
            
            # sequence structure for I & Q MW channels 
            mw_I_seq = [(iq_off1, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0], self.PiHalf('x', pihalf_x)[0], (iq_off5, self.IQ0[0])]
            mw_Q_seq = [(iq_off1, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1], (iq_off2, self.IQ0[1]), self.PiHalf('-y', pihalf_y)[1], (iq_off3, self.IQ0[1]), 
                        self.PiHalf('y', pihalf_y)[1], (iq_off4, self.IQ0[1]), self.PiHalf('-x', pihalf_x)[1], (iq_off5, self.IQ0[1])]

            # sequence structure for I & Q MW channels (MW off)
            mw_I_ref_seq = [(iq_off1, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0], (iq_off2, self.IQ0[0]), self.Pi('y', pi_y)[0], (iq_off3, self.IQ0[0]), self.PiHalf('-x', pihalf_x)[0], (iq_off4, self.IQ0[0])]
            mw_Q_ref_seq = [(iq_off1, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1], (iq_off2, self.IQ0[1]), self.Pi('y', pi_y)[1], (iq_off3, self.IQ0[1]), self.PiHalf('-x', pihalf_x)[1], (iq_off4, self.IQ0[1])]

            # switch_on_seq = [(iq_off1, 0), (pihalf_x, 1), (iq_off2, 0), (pi_y, 1), (iq_off3, 0), (pihalf_x, 1), (iq_off4, 0)]
            # switch_off_seq = [(iq_off1, 0), (pihalf_x, 0), (iq_off2, 0), (pi_y, 0), (iq_off3, 0), (pihalf_x, 0), (iq_off4, 0)]

            # switch_on_seq = [(iq_off1 - 20, 0), (pihalf_x + 40, 1), (iq_off2 - 40, 0), (pi_y + 40, 1), (iq_off3 - 40, 0), (pihalf_x + 40, 1), (iq_off4 - 20, 0)]
            # switch_off_seq = [(iq_off1 - 20, 0), (pihalf_x + 40, 0), (iq_off2 - 40, 0), (pi_y + 40, 0), (iq_off3 - 40, 0), (pihalf_x + 40, 0), (iq_off4 - 20, 0)]

            # assign sequences to respective channels for seq_on
            seq.setDigital(3, laser_seq) # laser
            seq.setDigital(0, daq_clock_seq) # integrator trigger
            # seq1.setDigital(1, switch_on_seq) # RF switch 
            seq.setAnalog(0, mw_I_seq) # mw_I
            seq.setAnalog(1, mw_Q_seq) # mw_Q
            
            # assign sequences to respective channels for seq_off
            seq_ref.setDigital(3, laser_seq) # laser
            seq_ref.setDigital(0, daq_clock_seq) # integrator trigger
            # seq2.setDigital(1, switch_off_seq) # RF switch 
            seq_ref.setAnalog(0, mw_I_ref_seq) # mw_I
            seq_ref.setAnalog(1, mw_Q_ref_seq) # mw_Q

            return seq + seq_ref

        # concatenate single ODMR sequence "runs" number of times
        seqs = self.Pulser.createSequence()
        
        for tau in params:
            seqs += SingleWAHUHA(tau)

        return seqs 
        
    def Echo_WAHUHA():
        pass

    def XY4_N(self, params, pulse_axes, pihalf_x, pihalf_y, pi_x, pi_y, n):
        '''
        XY4-N pulse sequence.
        MW sequence: pi/2(x) - tau/2 - (pi(x) - tau - pi(y) - tau - pi(x) - tau - pi(y))^N - tau/2 - pi/2(x, -x)
        '''
        longest_time = self.convert_type(round(params[-1]), float)
        pihalf_x = self.convert_type(round(pihalf_x), float)
        pihalf_y = self.convert_type(round(pihalf_y), float)
        pi_x = self.convert_type(round(pi_x), float)
        pi_y = self.convert_type(round(pi_y), float)
        n = self.convert_type(round(n), int)

        def PiPulsesN(axes, tau, N):            
            if axes == 'xy':
                xy4_I_seq = [((tau/2), self.IQ0[0]), self.Pi('x', pi_x)[0], (tau, self.IQ0[0]), self.Pi('y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('x', pi_x)[0], (tau, self.IQ0[0]), self.Pi('y', pi_y)[0], ((tau/2), self.IQ0[0])]
                xy4_Q_seq = [((tau/2), self.IQ0[1]), self.Pi('x', pi_x)[1], (tau, self.IQ0[1]), self.Pi('y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('x', pi_x)[1], (tau, self.IQ0[1]), self.Pi('y', pi_y)[1], ((tau/2), self.IQ0[1])]
                mw_I = (xy4_I_seq)*N
                mw_Q = (xy4_Q_seq)*N
                
            elif axes == 'yy':
                yy4_I_seq = [((tau/2), self.IQ0[0]), self.Pi('y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('y', pi_y)[0], ((tau/2), self.IQ0[0])]
                yy4_Q_seq = [((tau/2), self.IQ0[1]), self.Pi('y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('y', pi_y)[1], ((tau/2), self.IQ0[1])]
                mw_I = (yy4_I_seq)*N
                mw_Q = (yy4_Q_seq)*N

            return mw_I, mw_Q

        def SingleXY4(tau):
            '''
            CREATE SINGLE HAHN-ECHO SEQUENCE TO REPEAT THROUGHOUT EXPERIMENT
            '''

            tau = self.convert_type(round(tau), float) # convert to proper data type to avoid undesired rpyc netref data type

            '''
            DEFINE SPECIAL TIME INTERVALS FOR EXPERIMENT

            pad_time = padding time to equalize duration of every run (for different tau durations)
            '''
            pad_time = longest_time - tau 
            # NOTICE: change if using PiHalf['y'] to pihalf_y
            # xy4_time = 2*pihalf_x + (2*(tau/2)/(4*n) + 2*pi_x + 2*pi_y + 3*tau/(4*n))*n
            xy4_time = 2*pihalf_x + (2*(tau/2) + 2*pi_x + 2*pi_y + 3*tau)*n
            yy4_time = 2*pihalf_x + (2*(tau/2) + 4*pi_y + 3*tau)*n

            '''
            DEFINE RELEVANT ON, OFF TIMES FOR DEVICES
            '''            
            laser_off1 = self.initial_delay

            if pulse_axes == 'xy':
                laser_off2 = self.singlet_decay + xy4_time + self.MW_buffer_time
            else:
                laser_off2 = self.singlet_decay + yy4_time + self.MW_buffer_time

            laser_off3 = 100 + pad_time 
            
            # DAQ trigger windows
            clock_off1 = laser_off1 + self.laser_time + laser_off2 + self.readout_time - self.trig_spot - self.clock_time
            clock_off2 = self.trig_spot + laser_off3          

            # mw I & Q off windows (on slightly longer than VSG to ensure it's set)
            iq_off_start = laser_off1 + self.laser_time + self.singlet_decay
            iq_off_end = self.MW_buffer_time + self.readout_time + laser_off3
            
            '''
            CONSTRUCT PULSE SEQUENCE
            '''
            # create sequence objects for MW on and off blocks
            seq = self.Pulser.createSequence()
            seq_ref = self.Pulser.createSequence()

            # define sequence structure for laser
            laser_seq = [(laser_off1, 0), (self.laser_time, 1), (laser_off2, 0), (self.readout_time, 1), (laser_off3, 0)]
            
            # define sequence structure for integrator trigger
            daq_clock_seq = [(clock_off1, 0), (self.clock_time, 1), (clock_off2, 0)]
            
            # sequence structure for I & Q MW channels 
            mw_I_seq = [(iq_off_start, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0]] + PiPulsesN(pulse_axes, tau, n)[0] + [self.PiHalf('x', pihalf_x)[0], (iq_off_end, self.IQ0[0])]
            mw_Q_seq = [(iq_off_start, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1]] + PiPulsesN(pulse_axes, tau, n)[1] + [self.PiHalf('x', pihalf_x)[1], (iq_off_end, self.IQ0[1])]
            
            # sequence structure for I & Q MW channels (MW off)
            mw_I_ref_seq = [(iq_off_start, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0]] + PiPulsesN(pulse_axes, tau, n)[0] + [self.PiHalf('-x', pihalf_x)[0], (iq_off_end, self.IQ0[0])]
            mw_Q_ref_seq = [(iq_off_start, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1]] + PiPulsesN(pulse_axes, tau, n)[1] + [self.PiHalf('-x', pihalf_x)[1], (iq_off_end, self.IQ0[1])]

            # assign sequences to respective channels for seq_on
            seq.setDigital(3, laser_seq) # laser
            seq.setDigital(0, daq_clock_seq) # integrator trigger
            seq.setAnalog(0, mw_I_seq) # mw_I
            seq.setAnalog(1, mw_Q_seq) # mw_Q
            
            # assign sequences to respective channels for seq_off
            seq_ref.setDigital(3, laser_seq) # laser
            seq_ref.setDigital(0, daq_clock_seq) # integrator trigger
            seq_ref.setAnalog(0, mw_I_ref_seq) # mw_I
            seq_ref.setAnalog(1, mw_Q_ref_seq) # mw_Q

            return seq + seq_ref

        # concatenate single ODMR sequence "runs" number of times
        seqs = self.Pulser.createSequence()
        
        for tau in params:
            seqs += SingleXY4(tau)
        
        return seqs

    def XY8_N(self, params, pulse_axes, pihalf_x, pihalf_y, pi_x, pi_y, n):
        '''
        XY8-N pulse sequence.
        MW sequence: pi/2(x) - (tau/2 - pi(x) - tau - pi(y) - tau - pi(x) - tau... - pi(y) - tau/2)^N - pi/2(x, -x)

        '''
        longest_time = self.convert_type(round(params[-1]), float)
        pihalf_x = self.convert_type(round(pihalf_x), float)
        pihalf_y = self.convert_type(round(pihalf_y), float)
        pi_x = self.convert_type(round(pi_x), float)
        pi_y = self.convert_type(round(pi_y), float)
        n = self.convert_type(round(n), int)
        
        def PiPulsesN(axes, tau, N):
            tau_half_I_seq = [((tau/2), self.IQ0[0])]
            tau_half_Q_seq = [((tau/2), self.IQ0[1])]
            tau_I_seq = [(tau, self.IQ0[0])]
            tau_Q_seq = [(tau, self.IQ0[1])]

            if axes == 'xy':
                xy4_I_seq = [self.Pi('x', pi_x)[0], (tau, self.IQ0[0]), self.Pi('y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('x', pi_x)[0], (tau, self.IQ0[0]), self.Pi('y', pi_y)[0]]
                xy4_Q_seq = [self.Pi('x', pi_x)[1], (tau, self.IQ0[1]), self.Pi('y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('x', pi_x)[1], (tau, self.IQ0[1]), self.Pi('y', pi_y)[1]]
                mw_I = (tau_half_I_seq + xy4_I_seq + tau_I_seq + list(reversed(xy4_I_seq)) + tau_half_I_seq)*N
                mw_Q = (tau_half_Q_seq + xy4_Q_seq + tau_Q_seq + list(reversed(xy4_Q_seq)) + tau_half_Q_seq)*N
                
            elif axes == 'yy':
                yy4_I_seq_1 = [self.Pi('-y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('-y', pi_y)[0]]
                yy4_Q_seq_1 = [self.Pi('-y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('-y', pi_y)[1]]
                yy4_I_seq_2 = [self.Pi('-y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('-y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('y', pi_y)[0]]
                yy4_Q_seq_2 = [self.Pi('-y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('-y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('y', pi_y)[1]]
                mw_I = (tau_half_I_seq + yy4_I_seq_1 + tau_I_seq + yy4_I_seq_2 + tau_half_I_seq)*N
                mw_Q = (tau_half_Q_seq + yy4_Q_seq_1 + tau_Q_seq + yy4_Q_seq_2 + tau_half_Q_seq)*N
            
            return mw_I, mw_Q

        def SingleXY8(tau):
            '''
            CREATE SINGLE HAHN-ECHO SEQUENCE TO REPEAT THROUGHOUT EXPERIMENT
            '''

            tau = self.convert_type(round(tau), float) # convert to proper data type to avoid undesired rpyc netref data type

            '''
            DEFINE SPECIAL TIME INTERVALS FOR EXPERIMENT

            pad_time = padding time to equalize duration of every run (for different tau durations)
            '''
            pad_time = longest_time - tau 

            # NOTICE: change if using PiHalf['y'] to pihalf_y
            # xy8_time = 2*pihalf_x + ((tau/2)/(8*n) + 4*pi_x + 4*pi_y + 7*tau/(8*n) + (tau/2)/(8*n))*n
            # yy8_time = 2*pihalf_x + ((tau/2)/(8*n) + 8*pi_y + 7*tau/(8*n) + (tau/2)/(8*n))*n
            xy8_time = 2*pihalf_x + ((tau/2) + 4*pi_x + 4*pi_y + 7*tau + (tau/2))*n
            yy8_time = 2*pihalf_x + ((tau/2) + 8*pi_y + 7*tau + (tau/2))*n

            '''
            DEFINE RELEVANT ON, OFF TIMES FOR DEVICES
            '''            
            laser_off1 = self.initial_delay
            if pulse_axes == 'xy':
                laser_off2 = self.singlet_decay + xy8_time + self.MW_buffer_time
            else:
                laser_off2 = self.singlet_decay + yy8_time + self.MW_buffer_time
            laser_off3 = 100 + pad_time

            # DAQ trigger windows
            clock_off1 = laser_off1 + self.laser_time + laser_off2 + self.readout_time - self.trig_spot - self.clock_time
            clock_off2 = self.trig_spot + laser_off3           

            # mw I & Q off windows
            iq_off_start = laser_off1 + self.laser_time + self.singlet_decay
            iq_off_end = self.MW_buffer_time + self.readout_time + laser_off3
            
            '''
            CONSTRUCT PULSE SEQUENCE
            '''
            # create sequence objects for MW on and off blocks
            seq = self.Pulser.createSequence()
            seq_ref = self.Pulser.createSequence()

            # define sequence structure for laser
            laser_seq = [(laser_off1, 0), (self.laser_time, 1), (laser_off2, 0), (self.readout_time, 1), (laser_off3, 0)]
            
            # define sequence structure for integrator trigger
            daq_clock_seq = [(clock_off1, 0), (self.clock_time, 1), (clock_off2, 0)]
            
            # sequence structure for I & Q MW channels 
            mw_I_seq = [(iq_off_start, self.IQ0[0]), self.PiHalf('y', pihalf_y)[0]] + PiPulsesN(pulse_axes, tau, n)[0] + [self.PiHalf('y', pihalf_y)[0], (iq_off_end, self.IQ0[0])]
            mw_Q_seq = [(iq_off_start, self.IQ0[1]), self.PiHalf('y', pihalf_y)[1]] + PiPulsesN(pulse_axes, tau, n)[1] + [self.PiHalf('y', pihalf_y)[1], (iq_off_end, self.IQ0[1])]
            
            # sequence structure for I & Q MW channels (MW off)
            mw_I_ref_seq = [(iq_off_start, self.IQ0[0]), self.PiHalf('y', pihalf_y)[0]] + PiPulsesN(pulse_axes, tau, n)[0] + [self.PiHalf('-y', pihalf_y)[0], (iq_off_end, self.IQ0[0])]
            mw_Q_ref_seq = [(iq_off_start, self.IQ0[1]), self.PiHalf('y', pihalf_y)[1]] + PiPulsesN(pulse_axes, tau, n)[1] + [self.PiHalf('-y', pihalf_y)[1], (iq_off_end, self.IQ0[1])]

            # assign sequences to respective channels for seq_on
            seq.setDigital(3, laser_seq) # laser
            seq.setDigital(0, daq_clock_seq) # integrator trigger
            seq.setAnalog(0, mw_I_seq) # mw_I
            seq.setAnalog(1, mw_Q_seq) # mw_Q
            
            # assign sequences to respective channels for seq_off
            seq_ref.setDigital(3, laser_seq) # laser
            seq_ref.setDigital(0, daq_clock_seq) # integrator trigger
            seq_ref.setAnalog(0, mw_I_ref_seq) # mw_I
            seq_ref.setAnalog(1, mw_Q_ref_seq) # mw_Q

            return seq + seq_ref

        # concatenate single XY8 sequence "runs" number of times
        seqs = self.Pulser.createSequence()
        
        for tau in params:
            seqs += SingleXY8(tau)

        return seqs

    def CPMG_N(self, params, pulse_axis, pihalf_x, pihalf_y, pi_x, pi_y, n):
        '''
        CPMG-N pulse sequence.
        MW sequence: pi/2(x) - tau/2N - (pi(y) - tau/N - ...)^N - tau/2N - pi/2(x, -x)

        '''
        longest_time = self.convert_type(round(params[-1]), float)
        pihalf_x = self.convert_type(round(pihalf_x), float)
        pihalf_y = self.convert_type(round(pihalf_y), float)
        pi_x = self.convert_type(round(pi_x), float)
        pi_y = self.convert_type(round(pi_y), float)
        n = self.convert_type(round(n), int)
        
        print("long time = ", longest_time)
        print("pi/2 x = ", pihalf_x)
        print("n = ", n)

        def PiPulsesN(axis, tau, N):
            if axis == 'X':
                CPMG_I_seq = [self.Pi('x', pi_x)[0], (tau, self.IQ0[0])]
                CPMG_Q_seq = [self.Pi('x', pi_x)[1], (tau, self.IQ0[1])]
                mw_I = CPMG_I_seq*(N-1) + [self.Pi('x', pi_x)[0]]
                mw_Q = CPMG_Q_seq*(N-1) + [self.Pi('x', pi_x)[1]]

            elif axis == 'Y':
                CPMG_I_seq = [self.Pi('y', pi_y)[0], (tau, self.IQ0[0])]
                CPMG_Q_seq = [self.Pi('y', pi_y)[1], (tau, self.IQ0[1])]
                mw_I = CPMG_I_seq*(N-1) + [self.Pi('y', pi_y)[0]]
                mw_Q = CPMG_Q_seq*(N-1) + [self.Pi('y', pi_y)[1]]

            return mw_I, mw_Q

        def SingleCPMG(tau):
            '''
            CREATE SINGLE HAHN-ECHO SEQUENCE TO REPEAT THROUGHOUT EXPERIMENT
            '''
            print("Orig Tau = ", tau)
            tau = self.convert_type(round(tau), float) # convert to proper data type to avoid undesired rpyc netref data type
            print("Converted Tau = ", tau)
            '''
            DEFINE SPECIAL TIME INTERVALS FOR EXPERIMENT

            pad_time = padding time to equalize duration of every run (for different tau durations)
            '''
            pad_time = longest_time - tau 
            
            # if pulse_axis == 'X':
            #     cpmg_time = pihalf_x + tau/(2*n) + (pi_x + tau/n)*(n-1) + pi_x + tau/(2*n) + pihalf_x
            # elif pulse_axis == 'Y':
            #     cpmg_time = pihalf_x + tau/(2*n) + (pi_y + tau/n)*(n-1) + pi_y + tau/(2*n) + pihalf_x

            if pulse_axis == 'X':
                cpmg_time = pihalf_x + tau/(2) + (pi_x + tau)*(n-1) + pi_x + tau/(2) + pihalf_x
            elif pulse_axis == 'Y':
                cpmg_time = pihalf_x + tau/(2) + (pi_y + tau)*(n-1) + pi_y + tau/(2) + pihalf_x
            '''
            DEFINE RELEVANT ON, OFF TIMES FOR DEVICES
            '''       
            laser_off1 = self.initial_delay     
            laser_off2 = self.singlet_decay + cpmg_time + self.MW_buffer_time
            laser_off3 = pad_time + self.rest_time_btw_seqs
            
            # DAQ trigger windows
            clock_off1 = laser_off1 + self.laser_time + laser_off2 + self.readout_time - self.trig_spot - self.clock_time
            clock_off2 = self.trig_spot + laser_off3

            # mw I & Q off windows
            iq_off1 = self.laser_time + self.singlet_decay
            iq_off2 = self.MW_buffer_time + self.readout_time + laser_off3
            
            '''
            CONSTRUCT PULSE SEQUENCE
            '''
            # create sequence objects for MW on and off blocks
            seq = self.Pulser.createSequence()
            seq_ref = self.Pulser.createSequence()

            # define sequence structure for laser
            laser_seq = [(laser_off1, 0), (self.laser_time, 1), (laser_off2, 0), (self.readout_time, 1), (laser_off3, 0)]
            
            # define sequence structure for integrator trigger
            daq_clock_seq = [(clock_off1, 0), (self.clock_time, 1), (clock_off2, 0)]
            
            # sequence structure for I & Q MW channels 
            mw_I_seq = [(iq_off1, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0], (tau/(2), self.IQ0[0])] + PiPulsesN(pulse_axis, tau, n)[0] + [(tau/(2), self.IQ0[0]), self.PiHalf('x', pihalf_x)[0], (iq_off2, self.IQ0[0])]
            mw_Q_seq = [(iq_off1, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1], (tau/(2), self.IQ0[1])] + PiPulsesN(pulse_axis, tau, n)[1] + [(tau/(2), self.IQ0[1]), self.PiHalf('x', pihalf_x)[1], (iq_off2, self.IQ0[1])]
            
            # sequence structure for I & Q MW channels (MW off)
            mw_I_ref_seq = [(iq_off1, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0], (tau/(2), self.IQ0[0])] + PiPulsesN(pulse_axis, tau, n)[0] + [(tau/(2), self.IQ0[0]), self.PiHalf('-x', pihalf_x)[0], (iq_off2, self.IQ0[0])]
            mw_Q_ref_seq = [(iq_off1, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1], (tau/(2), self.IQ0[1])] + PiPulsesN(pulse_axis, tau, n)[1] + [(tau/(2), self.IQ0[1]), self.PiHalf('-x', pihalf_x)[1], (iq_off2, self.IQ0[1])]

            # assign sequences to respective channels for seq_on
            seq.setDigital(3, laser_seq) # laser
            seq.setDigital(0, daq_clock_seq) # integrator trigger
            seq.setAnalog(0, mw_I_seq) # mw_I
            seq.setAnalog(1, mw_Q_seq) # mw_Q
            
            # assign sequences to respective channels for seq_off
            seq_ref.setDigital(3, laser_seq) # laser
            seq_ref.setDigital(0, daq_clock_seq) # integrator trigger
            seq_ref.setAnalog(1, mw_Q_ref_seq) # mw_Q
            seq_ref.setAnalog(0, mw_I_ref_seq) # mw_I

            return seq + seq_ref

        # concatenate single CPMG sequence "runs" number of times
        seqs = self.Pulser.createSequence()
        
        for tau in params:
            seqs += SingleCPMG(tau)
        
        return seqs

    def DEER(self, pihalf_x, pihalf_y, pi_x, pi_y, tau, num_freqs=0):
        '''
        DEER pulse sequence.
        MW sequence: pi/2(x) - tau - pi(y) - tau - pi/2(x)

        '''
        pihalf_x = self.convert_type(round(pihalf_x), float)
        pihalf_y = self.convert_type(round(pihalf_y), float)
        pi_x = self.convert_type(round(pi_x), float)
        pi_y = self.convert_type(round(pi_y), float)
        tau = self.convert_type(round(tau), float)

        def SingleDEER():
            '''
            CREATE SINGLE DEER SEQUENCE TO REPEAT THROUGHOUT EXPERIMENT
            '''

            '''
            DEFINE RELEVANT ON, OFF TIMES FOR DEVICES
            '''            
            # laser_off = self.initial_delay + pihalf_x + tau + pi_y + tau + pihalf_x + self.MW_buffer_time
            # laser_on = self.laser_time
            laser_off1 = self.initial_delay
            laser_off2 = self.singlet_decay + pihalf_x + tau + pi_y + tau + pihalf_x + self.MW_buffer_time
            laser_off3 = self.initial_delay + 1000

            # DAQ trigger windows
            # clock_off1 = laser_off + self.readout_time - self.clock_time
            # clock_off2 = laser_on - self.readout_time
            clock_off1 = laser_off1 + self.laser_time + laser_off2 + self.readout_time - self.trig_spot - self.clock_time
            clock_off2 = self.trig_spot + laser_off3

            # mw I & Q off windows (on slightly longer than VSG to ensure it's set)
            # iq_off1 = self.initial_delay
            # iq_off2 = tau
            # iq_off3 = tau
            # iq_off4 = self.MW_buffer_time + laser_on
            iq_off1 = laser_off1 + self.laser_time + self.singlet_decay
            iq_off2 = tau 
            iq_off3 = tau
            iq_off4 = self.MW_buffer_time + self.readout_time + laser_off3

            awg_off1 = - self.initial_delay + iq_off1 + pihalf_x + self.awg_pulse_delay # additional initial delay at beginning to offset entire AWG pulse seq
            awg_off2 = (tau - self.awg_pulse_delay - self.awg_trig_time) + pi_y + self.awg_pulse_delay
            awg_off3 = (tau - self.awg_pulse_delay - self.awg_trig_time) + pihalf_x + iq_off4 + self.initial_delay

            '''
            CONSTRUCT PULSE SEQUENCE
            '''
            # create sequence objects for MW on and off blocks
            dark_seq = self.Pulser.createSequence()
            dark_seq_ref = self.Pulser.createSequence()
            echo_seq = self.Pulser.createSequence()
            echo_seq_ref = self.Pulser.createSequence()

            # define sequence structure for laser
            # laser_seq = [(laser_off, 0), (laser_on, 1)]
            laser_seq = [(laser_off1, 0), (self.laser_time, 1), (laser_off2, 0), (self.readout_time, 1), (laser_off3, 0)]
            # define sequence structure for integrator trigger
            daq_clock_seq = [(clock_off1, 0), (self.clock_time, 1), (clock_off2, 0)]
            
            # sequence structure for I & Q MW channels on SRS SG396
            srs_I_seq = [(iq_off1, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0], (iq_off2, self.IQ0[0]), self.Pi('y', pi_y)[0], (iq_off3, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0], (iq_off4, self.IQ0[0])]
            srs_Q_seq = [(iq_off1, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1], (iq_off2, self.IQ0[1]), self.Pi('y', pi_y)[1], (iq_off3, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1], (iq_off4, self.IQ0[1])]

            # sequence structure for I & Q MW channels (MW off)
            srs_I_ref_seq = [(iq_off1, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0], (iq_off2, self.IQ0[0]), self.Pi('y', pi_y)[0], (iq_off3, self.IQ0[0]), self.PiHalf('-x', pihalf_x)[0], (iq_off4, self.IQ0[0])]
            srs_Q_ref_seq = [(iq_off1, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1], (iq_off2, self.IQ0[1]), self.Pi('y', pi_y)[1], (iq_off3, self.IQ0[1]), self.PiHalf('-x', pihalf_x)[1], (iq_off4, self.IQ0[1])]

            # sequence structure for I & Q MW channels (MW off)
            awg_seq = [(awg_off1, 0), (self.awg_trig_time, 0), (awg_off2, 0), (self.awg_trig_time, 1), (awg_off3, 0)]
            awg_ref_seq = [(laser_off1 + self.laser_time + laser_off2 + self.readout_time + laser_off3, 0)] # off the entire time

            print("LASER SEQ: ", laser_seq)
            print("DAQ TRIG SEQ: ", daq_clock_seq)
            print("SRS SEQ: ", srs_I_seq)
            print("AWG_seq: ", awg_seq)

            # assign sequences to respective channels for seq_on
            dark_seq.setDigital(3, laser_seq) # laser
            dark_seq.setDigital(0, daq_clock_seq) # integrator trigger
            dark_seq.setDigital(4, awg_seq)
            dark_seq.setAnalog(0, srs_I_seq) # mw_I
            dark_seq.setAnalog(1, srs_Q_seq) # mw_Q
            
            # assign sequences to respective channels for seq_off
            dark_seq_ref.setDigital(3, laser_seq) # laser
            dark_seq_ref.setDigital(0, daq_clock_seq) # integrator trigger
            dark_seq_ref.setDigital(4, awg_seq)
            dark_seq_ref.setAnalog(0, srs_I_ref_seq) # mw_I
            dark_seq_ref.setAnalog(1, srs_Q_ref_seq) # mw_Q

            # assign sequences to respective channels for seq_on
            echo_seq.setDigital(3, laser_seq) # laser
            echo_seq.setDigital(0, daq_clock_seq) # integrator trigger
            echo_seq.setDigital(4, awg_ref_seq)
            echo_seq.setAnalog(0, srs_I_seq) # mw_I
            echo_seq.setAnalog(1, srs_Q_seq) # mw_Q
            
            # assign sequences to respective channels for seq_off
            echo_seq_ref.setDigital(3, laser_seq) # laser
            echo_seq_ref.setDigital(0, daq_clock_seq) # integrator trigger
            echo_seq_ref.setDigital(4, awg_ref_seq)
            echo_seq_ref.setAnalog(0, srs_I_ref_seq) # mw_I
            echo_seq_ref.setAnalog(1, srs_Q_ref_seq) # mw_Q

            return echo_seq + echo_seq_ref # dark_seq + dark_seq_ref + echo_seq + echo_seq_ref
        
        # seqs = self.Pulser.createSequence()

        # for i in range(num_freqs):
        #     seqs += SingleDEER()

        return SingleDEER()
        # return seqs 
    
    def DEER_FID(self, params, pihalf_x, pihalf_y, pi_x, pi_y):
        '''
        DEER FID pulse sequence.
        '''
        longest_time = self.convert_type(round(params[-1]), float)
        pihalf_x = self.convert_type(round(pihalf_x), float)
        pihalf_y = self.convert_type(round(pihalf_y), float)
        pi_x = self.convert_type(round(pi_x), float)
        pi_y = self.convert_type(round(pi_y), float)
        
        print("RUNNING DEER FID")
        def SingleDEERFID(tau):
            '''
            CREATE SINGLE DEER SEQUENCE TO REPEAT THROUGHOUT EXPERIMENT
            '''
            tau = self.convert_type(round(tau), float) # convert to proper data type to avoid undesired rpyc netref data type
            pad_time = longest_time - tau

            '''
            DEFINE RELEVANT ON, OFF TIMES FOR DEVICES
            '''            
            # laser_off = self.initial_delay + pihalf_x + tau + pi_y + tau + pihalf_x + self.MW_buffer_time
            # laser_on = self.laser_time
            laser_off1 = self.initial_delay
            laser_off2 = self.singlet_decay + pihalf_x + tau + pi_y + tau + pihalf_x + self.MW_buffer_time
            laser_off3 = self.initial_delay + pad_time

            # DAQ trigger windows
            # clock_off1 = laser_off + self.readout_time - self.clock_time
            # clock_off2 = laser_on - self.readout_time
            clock_off1 = laser_off1 + self.laser_time + laser_off2 + self.readout_time - self.trig_spot - self.clock_time
            clock_off2 = self.trig_spot + laser_off3

            # mw I & Q off windows (on slightly longer than VSG to ensure it's set)
            # iq_off1 = self.initial_delay
            # iq_off2 = tau
            # iq_off3 = tau
            # iq_off4 = self.MW_buffer_time + laser_on
            iq_off1 = laser_off1 + self.laser_time + self.singlet_decay
            iq_off2 = tau 
            iq_off3 = tau
            iq_off4 = self.MW_buffer_time + self.readout_time + laser_off3

            awg_off1 = iq_off1 + pihalf_x + iq_off2 + pi_y - self.awg_pulse_delay
            awg_off2 = self.awg_pulse_delay + iq_off3 + pihalf_x + iq_off4

            '''
            CONSTRUCT PULSE SEQUENCE
            '''
            # create sequence objects for MW on and off blocks
            dark_seq = self.Pulser.createSequence()
            dark_seq_ref = self.Pulser.createSequence()
            echo_seq = self.Pulser.createSequence()
            echo_seq_ref = self.Pulser.createSequence()

            # define sequence structure for laser
            # laser_seq = [(laser_off, 0), (laser_on, 1)]
            laser_seq = [(laser_off1, 0), (self.laser_time, 1), (laser_off2, 0), (self.readout_time, 1), (laser_off3, 0)]
            # define sequence structure for integrator trigger
            daq_clock_seq = [(clock_off1, 0), (self.clock_time, 1), (clock_off2, 0)]
            
            # sequence structure for I & Q MW channels on SRS SG396
            srs_I_seq = [(iq_off1, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0], (iq_off2, self.IQ0[0]), self.Pi('y', pi_y)[0], (iq_off3, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0], (iq_off4, self.IQ0[0])]
            srs_Q_seq = [(iq_off1, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1], (iq_off2, self.IQ0[1]), self.Pi('y', pi_y)[1], (iq_off3, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1], (iq_off4, self.IQ0[1])]

            # sequence structure for I & Q MW channels (MW off)
            srs_I_ref_seq = [(iq_off1, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0], (iq_off2, self.IQ0[0]), self.Pi('y', pi_y)[0], (iq_off3, self.IQ0[0]), self.PiHalf('-x', pihalf_x)[0], (iq_off4, self.IQ0[0])]
            srs_Q_ref_seq = [(iq_off1, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1], (iq_off2, self.IQ0[1]), self.Pi('y', pi_y)[1], (iq_off3, self.IQ0[1]), self.PiHalf('-x', pihalf_x)[1], (iq_off4, self.IQ0[1])]

            # sequence structure for I & Q MW channels (MW off)
            awg_seq = [(awg_off1, 0), (self.awg_trig_time, 1), (awg_off2, 0)]
            awg_ref_seq = [(awg_off1, 0), (self.awg_trig_time, 0), (awg_off2, 0)] # off the entire time

            print("LASER SEQ: ", laser_seq)
            print("DAQ TRIG SEQ: ", daq_clock_seq)
            print("SRS SEQ: ", srs_I_seq)
            print("AWG_seq: ", awg_seq)

            # assign sequences to respective channels for seq_on
            dark_seq.setDigital(3, laser_seq) # laser
            dark_seq.setDigital(0, daq_clock_seq) # integrator trigger
            dark_seq.setDigital(4, awg_seq)
            dark_seq.setAnalog(0, srs_I_seq) # mw_I
            dark_seq.setAnalog(1, srs_Q_seq) # mw_Q
            
            # assign sequences to respective channels for seq_off
            dark_seq_ref.setDigital(3, laser_seq) # laser
            dark_seq_ref.setDigital(0, daq_clock_seq) # integrator trigger
            dark_seq_ref.setDigital(4, awg_seq)
            dark_seq_ref.setAnalog(0, srs_I_ref_seq) # mw_I
            dark_seq_ref.setAnalog(1, srs_Q_ref_seq) # mw_Q

            # assign sequences to respective channels for seq_on
            echo_seq.setDigital(3, laser_seq) # laser
            echo_seq.setDigital(0, daq_clock_seq) # integrator trigger
            echo_seq.setDigital(4, awg_ref_seq)
            echo_seq.setAnalog(0, srs_I_seq) # mw_I
            echo_seq.setAnalog(1, srs_Q_seq) # mw_Q
            
            # assign sequences to respective channels for seq_off
            echo_seq_ref.setDigital(3, laser_seq) # laser
            echo_seq_ref.setDigital(0, daq_clock_seq) # integrator trigger
            echo_seq_ref.setDigital(4, awg_ref_seq)
            echo_seq_ref.setAnalog(0, srs_I_ref_seq) # mw_I
            echo_seq_ref.setAnalog(1, srs_Q_ref_seq) # mw_Q

            return dark_seq + dark_seq_ref + echo_seq + echo_seq_ref
        
        seqs = self.Pulser.createSequence()

        for tau in params:
            seqs += SingleDEERFID(tau)

        return seqs 

    def Corr_Spectroscopy(self, params, tau, pihalf_x, pihalf_y, pi_x, pi_y, n):
        '''
        Correlation Spectroscopy sequence using YY8-N.
        MW sequence: pi/2(x) - tau/2 - (pi(x) - tau - pi(y) - tau - pi(x) - tau...)^N - tau/2 - pi/2(x, -x)
        '''
        longest_time = self.convert_type(round(params[-1]), float)
        tau = self.convert_type(round(tau), float)
        pihalf_x = self.convert_type(round(pihalf_x), float)
        pihalf_y = self.convert_type(round(pihalf_y), float)
        pi_x = self.convert_type(round(pi_x), float)
        pi_y = self.convert_type(round(pi_y), float)
        n = self.convert_type(round(n), int)
        
        def PiPulsesN(tau, N):
            tau_half_I_seq = [((tau/2), self.IQ0[0])]
            tau_half_Q_seq = [((tau/2), self.IQ0[1])]
            tau_I_seq = [(tau, self.IQ0[0])]
            tau_Q_seq = [(tau, self.IQ0[1])]

            # xy4_I_seq = [Pi('x')[0], (tau, self.IQ0[0]), Pi('y')[0], (tau, self.IQ0[0]), Pi('x')[0], (tau, self.IQ0[0]), Pi('y')[0]]
            # xy4_Q_seq = [Pi('x')[1], (tau, self.IQ0[1]), Pi('y')[1], (tau, self.IQ0[1]), Pi('x')[1], (tau, self.IQ0[1]), Pi('y')[1]]
            # mw_I = (tau_half_I_seq + xy4_I_seq + tau_I_seq + list(reversed(xy4_I_seq)) + tau_half_I_seq)*N
            # mw_Q = (tau_half_Q_seq + xy4_Q_seq + tau_Q_seq + list(reversed(xy4_Q_seq)) + tau_half_Q_seq)*N

            yy4_I_seq_1 = [self.Pi('-y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('-y', pi_y)[0]]
            yy4_Q_seq_1 = [self.Pi('-y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('-y', pi_y)[1]]
            yy4_I_seq_2 = [self.Pi('-y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('-y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('y', pi_y)[0]]
            yy4_Q_seq_2 = [self.Pi('-y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('-y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('y', pi_y)[1]]
            mw_I = (tau_half_I_seq + yy4_I_seq_1 + tau_I_seq + yy4_I_seq_2 + tau_half_I_seq)*N
            mw_Q = (tau_half_Q_seq + yy4_Q_seq_1 + tau_Q_seq + yy4_Q_seq_2 + tau_half_Q_seq)*N

            return mw_I, mw_Q

        def SingleCorrSpec(t_corr):
            '''
            CREATE SINGLE HAHN-ECHO SEQUENCE TO REPEAT THROUGHOUT EXPERIMENT
            '''

            t_corr = self.convert_type(round(t_corr), float) # convert to proper data type to avoid undesired rpyc netref data type

            '''
            DEFINE SPECIAL TIME INTERVALS FOR EXPERIMENT

            pad_time = padding time to equalize duration of every run (for different tau durations)
            '''
            pad_time = longest_time - t_corr 

            # total time for correlation spectroscopy MW pulse sequence
            # corr_spec_time = pihalf_x + ((tau/2)/(8*n) + 4*pi_x + 4*pi_y + 7*tau/(8*n) + (tau/2)/(8*n))*n + pihalf_y + t_corr + \
            #                  pihalf_x + ((tau/2)/(8*n) + 4*pi_x + 4*pi_y + 7*tau/(8*n) + (tau/2)/(8*n))*n + pihalf_y
            corr_spec_time = pihalf_x + ((tau/2) + 0*pi_x + 8*pi_y + 7*tau + (tau/2))*n + pihalf_y + t_corr + \
                             pihalf_x + ((tau/2) + 0*pi_x + 8*pi_y + 7*tau + (tau/2))*n + pihalf_y
            '''
            DEFINE RELEVANT ON, OFF TIMES FOR DEVICES
            '''            
            laser_off1 = self.initial_delay
            laser_off2 = self.singlet_decay + corr_spec_time + self.MW_buffer_time
            laser_off3 = 100 + pad_time 

            # DAQ trigger windows
            clock_off1 = laser_off1 + self.laser_time + laser_off2 + self.readout_time - self.trig_spot - self.clock_time
            clock_off2 = self.trig_spot + laser_off3     

            # mw I & Q off windows (on slightly longer than VSG to ensure it's set)
            iq_off_start = laser_off1 + self.laser_time + self.singlet_decay
            iq_off_end = self.MW_buffer_time + self.readout_time + laser_off3
            
            '''
            CONSTRUCT PULSE SEQUENCE
            '''
            # create sequence objects for MW on and off blocks
            seq = self.Pulser.createSequence()
            seq_ref = self.Pulser.createSequence()

            # define sequence structure for laser
            laser_seq = [(laser_off1, 0), (self.laser_time, 1), (laser_off2, 0), (self.readout_time, 1), (laser_off3, 0)]
            
            # define sequence structure for integrator trigger
            daq_clock_seq = [(clock_off1, 0), (self.clock_time, 1), (clock_off2, 0)]
            
            # sequence structure for I & Q MW channels 
            mw_I_seq = [(iq_off_start, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0]] + PiPulsesN(tau, n)[0] + [self.PiHalf('y', pihalf_y)[0], (t_corr, self.IQ0[0]), self.PiHalf('x', pihalf_y)[0]] + PiPulsesN(tau, n)[0] + [self.PiHalf('y', pihalf_x)[0], (iq_off_end, self.IQ0[0])]
            mw_Q_seq = [(iq_off_start, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1]] + PiPulsesN(tau, n)[1] + [self.PiHalf('y', pihalf_y)[1], (t_corr, self.IQ0[1]), self.PiHalf('x', pihalf_y)[1]] + PiPulsesN(tau, n)[1] + [self.PiHalf('y', pihalf_x)[1], (iq_off_end, self.IQ0[1])]
            
            # sequence structure for I & Q MW channels (MW off)
            mw_I_ref_seq = [(iq_off_start, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0]] + PiPulsesN(tau, n)[0] + [self.PiHalf('y', pihalf_y)[0], (t_corr, self.IQ0[0]), self.PiHalf('x', pihalf_y)[0]] + PiPulsesN(tau, n)[0] + [self.PiHalf('-y', pihalf_x)[0], (iq_off_end, self.IQ0[0])]
            mw_Q_ref_seq = [(iq_off_start, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1]] + PiPulsesN(tau, n)[1] + [self.PiHalf('y', pihalf_y)[1], (t_corr, self.IQ0[1]), self.PiHalf('x', pihalf_y)[1]] + PiPulsesN(tau, n)[1] + [self.PiHalf('-y', pihalf_x)[1], (iq_off_end, self.IQ0[1])]
            
            # assign sequences to respective channels for seq_on
            seq.setDigital(3, laser_seq) # laser
            seq.setDigital(0, daq_clock_seq) # integrator trigger
            seq.setAnalog(0, mw_I_seq) # mw_I
            seq.setAnalog(1, mw_Q_seq) # mw_Q
        
            # assign sequences to respective channels for seq_off
            seq_ref.setDigital(3, laser_seq) # laser
            seq_ref.setDigital(0, daq_clock_seq) # integrator trigger
            seq_ref.setAnalog(0, mw_I_ref_seq) # mw_I
            seq_ref.setAnalog(1, mw_Q_ref_seq) # mw_Q

            return seq + seq_ref

        # concatenate single correlation spectroscopy sequence "runs" number of times
        seqs = self.Pulser.createSequence()
        
        for t_corr in params:
            seqs += SingleCorrSpec(t_corr)

        return seqs

    def CASR(self, params, tau, pihalf_x, pihalf_y, pi_x, pi_y, n):
        '''
        Coherently Averaged Synchronized Readout sequence using YY8-N.
        MW sequence: pi/2(x) - tau/2 - (pi(x) - tau - pi(y) - tau - pi(x) - tau...)^N - tau/2 - pi/2(x, -x)
        '''
        longest_time = self.convert_type(round(params[-1]), float)
        tau = self.convert_type(round(tau), float)
        pihalf_x = self.convert_type(round(pihalf_x), float)
        pihalf_y = self.convert_type(round(pihalf_y), float)
        pi_x = self.convert_type(round(pi_x), float)
        pi_y = self.convert_type(round(pi_y), float)
        n = self.convert_type(round(n), int)

        def PiPulsesN(tau, N):
            tau_half_I_seq = [((tau/2), self.IQ0[0])]
            tau_half_Q_seq = [((tau/2), self.IQ0[1])]
            tau_I_seq = [(tau, self.IQ0[0])]
            tau_Q_seq = [(tau, self.IQ0[1])]

            # xy4_I_seq = [Pi('x')[0], (tau, self.IQ0[0]), Pi('y')[0], (tau, self.IQ0[0]), Pi('x')[0], (tau, self.IQ0[0]), Pi('y')[0]]
            # xy4_Q_seq = [Pi('x')[1], (tau, self.IQ0[1]), Pi('y')[1], (tau, self.IQ0[1]), Pi('x')[1], (tau, self.IQ0[1]), Pi('y')[1]]
            # mw_I = (tau_half_I_seq + xy4_I_seq + tau_I_seq + list(reversed(xy4_I_seq)) + tau_half_I_seq)*N
            # mw_Q = (tau_half_Q_seq + xy4_Q_seq + tau_Q_seq + list(reversed(xy4_Q_seq)) + tau_half_Q_seq)*N

            yy4_I_seq_1 = [self.Pi('-y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('-y', pi_y)[0]]
            yy4_Q_seq_1 = [self.Pi('-y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('-y', pi_y)[1]]
            yy4_I_seq_2 = [self.Pi('-y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('-y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('y', pi_y)[0], (tau, self.IQ0[0]), self.Pi('y', pi_y)[0]]
            yy4_Q_seq_2 = [self.Pi('-y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('-y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('y', pi_y)[1], (tau, self.IQ0[1]), self.Pi('y', pi_y)[1]]
            mw_I = (tau_half_I_seq + yy4_I_seq_1 + tau_I_seq + yy4_I_seq_2 + tau_half_I_seq)*N
            mw_Q = (tau_half_Q_seq + yy4_Q_seq_1 + tau_Q_seq + yy4_Q_seq_2 + tau_half_Q_seq)*N

            return mw_I, mw_Q

        def SingleCASR(t_corr):
            '''
            CREATE SINGLE HAHN-ECHO SEQUENCE TO REPEAT THROUGHOUT EXPERIMENT
            '''

            t_corr = self.convert_type(round(t_corr), float) # convert to proper data type to avoid undesired rpyc netref data type

            '''
            DEFINE SPECIAL TIME INTERVALS FOR EXPERIMENT

            pad_time = padding time to equalize duration of every run (for different tau durations)
            '''
            pad_time = longest_time - t_corr 

            casr_time = pihalf_x + ((tau/2) + 0*pi_x + 8*pi_y + 7*tau + (tau/2))*n + pihalf_y + t_corr + \
                             pihalf_x + ((tau/2) + 0*pi_x + 8*pi_y + 7*tau + (tau/2))*n + pihalf_y
            # TODO: update CASR time

            '''
            DEFINE RELEVANT ON, OFF TIMES FOR DEVICES
            '''            
            laser_off1 = self.initial_delay
            laser_off2 = self.singlet_decay + casr_time + self.MW_buffer_time
            laser_off3 = 100 + pad_time

            # DAQ trigger windows
            clock_off1 = laser_off1 + self.laser_time + laser_off2 + self.readout_time - self.trig_spot - self.clock_time
            clock_off2 = self.trig_spot + laser_off3           

            # mw I & Q off windows (on slightly longer than VSG to ensure it's set)
            iq_off_start = laser_off1 + self.laser_time + self.singlet_decay
            iq_off_end = self.MW_buffer_time + self.readout_time + laser_off3
            
            '''
            CONSTRUCT PULSE SEQUENCE
            '''
            # create sequence objects for MW on and off blocks
            seq = self.Pulser.createSequence()
            seq_ref = self.Pulser.createSequence()

            # define sequence structure for laser
            laser_seq = [(laser_off1, 0), (self.laser_time, 1), (laser_off2, 0), (self.readout_time, 1), (laser_off3, 0)]
            
            # define sequence structure for integrator trigger
            daq_clock_seq = [(clock_off1, 0), (self.clock_time, 1), (clock_off2, 0)]
            
            # sequence structure for I & Q MW channels 
            mw_I_seq = [(iq_off_start, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0]] + PiPulsesN(tau, n)[0] + [self.PiHalf('y', pihalf_y)[0], (t_corr, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0]] + PiPulsesN(tau, n)[0] + [self.PiHalf('y', pihalf_y)[0], (iq_off_end, self.IQ0[0])]
            mw_Q_seq = [(iq_off_start, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1]] + PiPulsesN(tau, n)[1] + [self.PiHalf('y', pihalf_y)[1], (t_corr, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1]] + PiPulsesN(tau, n)[1] + [self.PiHalf('y', pihalf_y)[1], (iq_off_end, self.IQ0[1])]
            
            # sequence structure for I & Q MW channels (MW off)
            mw_I_ref_seq = [(iq_off_start, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0]] + PiPulsesN(tau, n)[0] + [self.PiHalf('y', pihalf_y)[0], (t_corr, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0]] + PiPulsesN(tau, n)[0] + [self.PiHalf('-y', pihalf_y)[0], (iq_off_end, self.IQ0[0])]
            mw_Q_ref_seq = [(iq_off_start, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1]] + PiPulsesN(tau, n)[1] + [self.PiHalf('y', pihalf_y)[1], (t_corr, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1]] + PiPulsesN(tau, n)[1] + [self.PiHalf('-y', pihalf_y)[1], (iq_off_end, self.IQ0[1])]

            print(mw_I_seq)
            
            # assign sequences to respective channels for seq_on
            seq.setDigital(3, laser_seq) # laser
            seq.setDigital(0, daq_clock_seq) # integrator trigger
            seq.setAnalog(0, mw_I_seq) # mw_I
            seq.setAnalog(1, mw_Q_seq) # mw_Q
            
            # print("SEQ 1: ", seq1)

            # assign sequences to respective channels for seq_off
            seq_ref.setDigital(3, laser_seq) # laser
            seq_ref.setDigital(0, daq_clock_seq) # integrator trigger
            seq_ref.setAnalog(0, mw_I_ref_seq) # mw_I
            seq_ref.setAnalog(1, mw_Q_ref_seq) # mw_Q

            # print("SEQ 2: ", seq2)

            return seq + seq_ref 

        # concatenate single correlation spectroscopy sequence "runs" number of times
        seqs = self.Pulser.createSequence()
        
        for t_corr in params:
            seqs += SingleCASR(t_corr)

        return seqs

    def AERIS(self, params, pulse_axes, pihalf_x, pihalf_y, pi_x, pi_y, n):
        '''
        Amplitude-Encoded Radio Induced Signal (AERIS) pulse sequence.
        '''
        longest_time = self.convert_type(round(params[-1]), float)
        pihalf_x = self.convert_type(round(pihalf_x), float)
        pihalf_y = self.convert_type(round(pihalf_y), float)
        pi_x = self.convert_type(round(pi_x), float)
        pi_y = self.convert_type(round(pi_y), float)
        n = self.convert_type(round(n), int)
        
        def PiPulsesN(axes, tau, N):
            if axes == 'xy':
                xy4_I_seq = [((tau/2)/(4*N), self.IQ0[0]), self.Pi('x', pi_x)[0], (tau/(4*N), self.IQ0[0]), self.Pi('y', pi_y)[0], (tau/(4*N), self.IQ0[0]), self.Pi('x', pi_x)[0], (tau/(4*N), self.IQ0[0]), self.Pi('y', pi_y)[0], ((tau/2)/(4*N), self.IQ0[0])]
                xy4_Q_seq = [((tau/2)/(4*N), self.IQ0[1]), self.Pi('x', pi_x)[1], (tau/(4*N), self.IQ0[1]), self.Pi('y', pi_y)[1], (tau/(4*N), self.IQ0[1]), self.Pi('x', pi_x)[1], (tau/(4*N), self.IQ0[1]), self.Pi('y', pi_y)[1], ((tau/2)/(4*N), self.IQ0[1])]
                mw_I = (xy4_I_seq)*N
                mw_Q = (xy4_Q_seq)*N
            elif axes == 'yy':
                yy4_I_seq = [((tau/2)/(4*N), self.IQ0[0]), self.Pi('y', pi_y)[0], (tau/(4*N), self.IQ0[0]), self.Pi('y', pi_y)[0], (tau/(4*N), self.IQ0[0]), self.Pi('y', pi_y)[0], (tau/(4*N), self.IQ0[0]), self.Pi('y', pi_y)[0], ((tau/2)/(4*N), self.IQ0[0])]
                yy4_Q_seq = [((tau/2)/(4*N), self.IQ0[1]), self.Pi('y', pi_y)[1], (tau/(4*N), self.IQ0[1]), self.Pi('y', pi_y)[1], (tau/(4*N), self.IQ0[1]), self.Pi('y', pi_y)[1], (tau/(4*N), self.IQ0[1]), self.Pi('y', pi_y)[1], ((tau/2)/(4*N), self.IQ0[1])]
                mw_I = (yy4_I_seq)*N
                mw_Q = (yy4_Q_seq)*N

            return mw_I, mw_Q

        def SingleAERIS(tau):
            '''
            CREATE SINGLE HAHN-ECHO SEQUENCE TO REPEAT THROUGHOUT EXPERIMENT
            '''

            tau = self.convert_type(round(tau), float) # convert to proper data type to avoid undesired rpyc netref data type

            '''
            DEFINE SPECIAL TIME INTERVALS FOR EXPERIMENT

            pad_time = padding time to equalize duration of every run (for different tau durations)
            '''
            pad_time = longest_time - tau 
            # NOTICE: change if using PiHalf['y'] to pihalf_y
            xy4_time = 2*pihalf_x + (2*(tau/2)/(4*n) + 2*pi_x + 2*pi_y + 3*tau/(4*n))*n
            
            '''
            DEFINE RELEVANT ON, OFF TIMES FOR DEVICES
            '''            
            laser_off1 = self.initial_delay
            laser_off2 = self.singlet_decay + xy4_time + self.MW_buffer_time
            laser_off3 = 100 + pad_time
            
            # DAQ trigger windows
            clock_off1 = laser_off1 + self.laser_time + laser_off2 + self.readout_time - self.trig_spot - self.clock_time
            clock_off2 = self.trig_spot + laser_off3         

            # mw I & Q off windows (on slightly longer than VSG to ensure it's set)
            iq_off_start = laser_off1 + self.laser_time + self.singlet_decay
            iq_off_end = self.MW_buffer_time + self.readout_time + laser_off3
            
            '''
            CONSTRUCT PULSE SEQUENCE
            '''
            # create sequence objects for MW on and off blocks
            seq = self.Pulser.createSequence()
            seq_ref = self.Pulser.createSequence()
            
            # define sequence structure for laser
            laser_seq = [(laser_off1, 0), (self.laser_time, 1), (laser_off2, 0), (self.readout_time, 1), (laser_off3, 0)]
            
            # define sequence structure for integrator trigger
            daq_clock_seq = [(clock_off1, 0), (self.clock_time, 1), (clock_off2, 0)]
            
            # sequence structure for I & Q MW channels 
            mw_I_seq = [(iq_off_start, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0]] + PiPulsesN(pulse_axes, tau, n)[0] + [self.PiHalf('x', pihalf_x)[0], (iq_off_end, self.IQ0[0])]
            mw_Q_seq = [(iq_off_start, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1]] + PiPulsesN(pulse_axes, tau, n)[1] + [self.PiHalf('x', pihalf_x)[1], (iq_off_end, self.IQ0[1])]
            
            # sequence structure for I & Q MW channels (MW off)
            mw_I_ref_seq = [(iq_off_start, self.IQ0[0]), self.PiHalf('x', pihalf_x)[0]] + PiPulsesN(pulse_axes, tau, n)[0] + [self.PiHalf('-x', pihalf_x)[0], (iq_off_end, self.IQ0[0])]
            mw_Q_ref_seq = [(iq_off_start, self.IQ0[1]), self.PiHalf('x', pihalf_x)[1]] + PiPulsesN(pulse_axes, tau, n)[1] + [self.PiHalf('-x', pihalf_x)[1], (iq_off_end, self.IQ0[1])]

            # sequence structure for nuclear spin RF generator 
            rf_seq = []
            rf_ref_seq = []

            # assign sequences to respective channels for seq_on
            seq.setDigital(3, laser_seq) # laser
            seq.setDigital(0, daq_clock_seq) # integrator trigger
            seq.setDigital(1, rf_seq) # VSG switch to enable MW
            seq.setAnalog(0, mw_I_seq) # mw_I
            seq.setAnalog(1, mw_Q_seq) # mw_Q
            
            # assign sequences to respective channels for seq_off
            seq_ref.setDigital(3, laser_seq) # laser
            seq_ref.setDigital(0, daq_clock_seq) # integrator trigger
            seq_ref.setDigital(1, rf_ref_seq) # VSG switch to enable MW
            seq_ref.setAnalog(0, mw_I_ref_seq) # mw_I
            seq_ref.setAnalog(1, mw_Q_ref_seq) # mw_Q

            return seq + seq_ref

        # concatenate single ODMR sequence "runs" number of times
        seqs = self.Pulser.createSequence()
        
        for tau in params:
            seqs += SingleAERIS(tau)
        
        return seqs