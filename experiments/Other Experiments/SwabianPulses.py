
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

# from drivers.swabian.pulsestreamer.lib.pulse_streamer_grpc import PulseStreamer
# from drivers.swabian.pulsestreamer.lib.Sequence import Sequence
from nspyre import InstrumentGateway as gw

class Pulses():
    '''
    ALL UNITS: [ns]
    '''
    # def __init__(self, 
    #             # laser_time = 5e3, initial_delay = 100, singlet_decay = 300, readout_time = 300, 
    #              # MW_buffer_time = 100, probe_time = 50e3, clock_time = 11, sampling_time = 50000, 
    #              # trig_spot = 50, awg_trig_time = 10, awg_pulse_delay = 0,
    #              # rest_time_btw_seqs = 100e3, 
    #              ip="192.168.1.76"):

    def __init__(self):
        pass

        '''
        :param channel_dict: Dictionary of which channels correspond to which instr controls
        :param readout_time: Laser+gate readout time in ns
        :param laser_time: Laser time to reinit post readout
        :param initial_delay: Delay in laser turning on
        :param MW_buffer_time: Buffer after MW turns off
        :param IQ: IQ modulation/analog channels
        '''
        self.channel_dict = {"DAQ": 0, "AOM": 1, "SRS": 4} 
        # SANS: CHECK ABOUT ANALOG DIRECTORY FOR SRS IQ MODULATION

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

        self.Pulser = PulseStreamer(ip)
        # print('creating sequence')
        self.sequence = Sequence()
        # print('done creating sequence')
        # self.latest_streamed = pd.DataFrame({})
        self.total_time = 0 #update when a pulse sequence is streamed

        ## Analog voltage levels to set for sig gen I and Q. If no crosstalk should be 0.5 and 0.
        self.IQ0 = [-0.0038, -0.0030]
        self.IQ = self.IQ0

        self.IQpx = [0.4862, -.0030]
        self.IQnx = [-0.4862,-0.0030]

        self.IQpy = [-0.0038, 0.4870]
        self.IQny = [-0.0038, -0.4870]
        
        self.IQboth = [0.4862, 0.4870]
        self.IQtest = [.95, 0]

    # def has_sequence(self):
    #     """
    #     Has Sequence
    #     """
    #     return self.Pulser.hasSequence()
    
    # def has_finished(self):
    #     """
    #     Has Finished
    #     """
    #     return self.Pulser.hasFinished()
    
    # def laser_on(self):
    #     return self.Pulser.constant((1, [7], 0.0, 0.0))

    # def stream(self,seq,n_runs):
    #     seq = obtain(seq)
    #     # print(type(seq))
    #     # print(seq)
    #     self.Pulser.stream(seq,n_runs)

    # def clocksource(self,clk_src): # SANS: what is this clock source?
    #     self.Pulser.selectClock(clk_src)

    # def _normalize_IQ(self, IQ):
    #     self.IQ = IQ/(2.5*np.linalg.norm(IQ))


    # _T = t.TypeVar('_T')

    # def convert_type(self, arg: t.Any, converter: _T) -> _T:
    #     return converter(arg)

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


    '''
    PULSE SEQUENCES FOR NanoNMR EXPERIMENTS
    '''

    def CountVsTime(self):
        # seq = self.Pulser.createSequence()
        
        # trig_off = self.sampling_time - self.clock_time
        # daq_clock_seq = [(trig_off, 0), (self.clock_time, 1)]

        # seq.setDigital(0, daq_clock_seq) # integrator trigger

        # return seq

        # constant trigger input to DAQ to count
        return self.Pulser.constant((0, [0], 0.0, 0.0))


    # def CW_ODMR(self, clock_time, probe_time):
        
    #     '''
    #     CW ODMR Sequence
    #     Laser on for entire sequence. 
    #     MW on for probe_time.
    #     MW off for probe_time.
    #     User sets how many voltage samples (num_clocks) to take during each MW on/off window.
    #     '''
    
    #     def SingleCW_ODMR(clock_time, probe_time):
            
    #         # create sequence object
    #         seq = self.Pulser.createSequence()
            
    #         # set DAQ trigger off time based on optimal readout window during MW on/off window
    #         # clock_on = self.clock_time
    #         # clock_off1 = 50000 - 2*self.clock_time # self.probe_time - 10*clock_on
    #         # clock_off2 = self.clock_time + clock_off1 # 9*clock_on + clock_off1
    #         # clock_off3 = self.clock_time # 9*clock_on

    #         # iq_on = self.probe_time 
    #         # iq_off = self.probe_time

    #         clock_on = clock_time
    #         clock_off1 = 50000 - 2*clock_time # self.probe_time - 10*clock_on
    #         clock_off2 = clock_time + clock_off1 # 9*clock_on + clock_off1
    #         clock_off3 = clock_time # 9*clock_on

    #         iq_on = probe_time 
    #         iq_off = probe_time

    #         # define sequence structure for clock and MW I/Q channels
    #         daq_clock_seq = [(clock_off1, 0), (clock_on, 1), (clock_off2, 0), (clock_on, 1), (clock_off3, 0)]
    #         mw_I_seq = [(iq_on, self.IQpx[0]), (iq_off, self.IQ0[0])]
    #         mw_Q_seq = [(iq_on, self.IQpx[1]), (iq_off, self.IQ0[1])]

    #         print("DAQ CLOCK: ", daq_clock_seq)
    #         print("MW I: ", mw_I_seq)
    #         # switch_seq = [(iq_on, 1), (iq_off, 0), (iq_off, 0), (iq_on, 1)]

    #         # assign sequences to respective channels
    #         seq.setDigital(0, daq_clock_seq) # DAQ clock -- record data
    #         # seq.setDigital(1, switch_seq) # RF switch
    #         seq.setAnalog(0, mw_I_seq) # mw_I
    #         seq.setAnalog(1, mw_Q_seq) # mw_Q

    #         return seq

    #     seqs = SingleCW_ODMR(clock_time, probe_time)

    #     return seqs


 #    def Pulsed_ODMR(self, pi_xy, pi_time):
 #        '''
 #        Pulsed ODMR sequence with integrator
 #        '''
 #        ## Run a pi pulse, then measure the signal
 #        ## and reference counts from NV.
 #        pi_time = self.convert_type(round(pi_time), float)
       
 #        ## we can measure the pi time on x and on y.
 #        ## they should be the same, but they technically
 #        ## have different offsets on our pulse streamer.
 #        def Pi(axis):
 #            iq_on = pi_time 
            
 #            if axis == 'x':
 #                mw_I_on = (iq_on, self.IQpx[0])
 #                mw_Q_on = (iq_on, self.IQpx[1])
 #            else:
 #                mw_I_on = (iq_on, self.IQpy[0])
 #                mw_Q_on = (iq_on, self.IQpy[1])
            
 #            return mw_I_on, mw_Q_on

 #        def SinglePulsed_ODMR():
 #            '''
 #            CREATE SINGLE RABI SEQUENCE TO REPEAT THROUGHOUT EXPERIMENT
 #            '''

 #            '''
 #            DEFINE SPECIAL TIME INTERVALS FOR EXPERIMENT

 #            pad_time: padding time to equalize duration of every run (for different vsg_on durations)
 #            '''
 #            pad_time = 50000 - self.initial_delay - self.laser_time - self.singlet_decay - pi_time - self.MW_buffer_time - self.readout_time 

 #            '''
 #            DEFINE RELEVANT ON, OFF TIMES FOR DEVICES
 #            '''
 #            laser_off1 = self.initial_delay
 #            laser_off2 = self.singlet_decay + pi_time + self.MW_buffer_time
 #            laser_off3 = pad_time

 #            # integrator trigger windows     
 #            int_trig_off1 = laser_off1 + self.laser_time + (laser_off2 - self.trig_delay)
 #            int_trig_off2 = (self.trig_delay - self.clock_time) + self.readout_time + pad_time            

 #            # mw I & Q off windows
 #            iq_off1 = laser_off1 + self.laser_time + self.singlet_decay
 #            iq_off2 = self.MW_buffer_time + self.readout_time + pad_time

 #            '''
 #            CONSTRUCT PULSE SEQUENCE
 #            '''
 #            # create sequence objects for MW on and off blocks
 #            seq_on = self.Pulser.createSequence()
 #            seq_off = self.Pulser.createSequence()

 #            # define sequence structure for laser
 #            laser_seq = [(laser_off1, 0), (self.laser_time, 1), (laser_off2, 0), (self.readout_time, 1), (laser_off3, 0)]

 #            # define sequence structure for integrator trigger
 #            int_trig_seq = [(int_trig_off1, 0), (self.clock_time, 1), (int_trig_off2, 0)]
            
 #            # define sequence structure for MW I and Q when MW = ON
 #            mw_I_on_seq = [(iq_off1, self.IQ0[0]), Pi(pi_xy)[0], (iq_off2, self.IQ0[0])]
 #            mw_Q_on_seq = [(iq_off1, self.IQ0[1]), Pi(pi_xy)[1], (iq_off2, self.IQ0[1])]
 #            # when MW = OFF
 #            mw_I_off_seq = [(iq_off1, self.IQ0[0]), (pi_time, self.IQ0[0]), (iq_off2, self.IQ0[0])]
 #            mw_Q_off_seq = [(iq_off1, self.IQ0[1]), (pi_time, self.IQ0[1]), (iq_off2, self.IQ0[1])]

 #            # assign sequences to respective channels for seq_on
 #            seq_on.setDigital(3, laser_seq) # laser
 #            # seq_on.setDigital(4, int_trig_seq) # integrator trigger
 #            seq_on.setAnalog(0, mw_I_on_seq) # mw_I
 #            seq_on.setAnalog(1, mw_Q_on_seq) # mw_Q

 #            # assign sequences to respective channels for seq_off
 #            seq_off.setDigital(3, laser_seq) # laser
 #            # seq_off.setDigital(4, int_trig_seq) # integrator trigger
 #            seq_off.setAnalog(0, mw_I_off_seq) # mw_I
 #            seq_off.setAnalog(1, mw_Q_off_seq) # mw_Q

 #            return seq_on + seq_off

 #        seqs = self.Pulser.createSequence()

 #        for i in range(self.runs):
 #            seqs += SinglePulsed_ODMR()

 #        return seqs


 # 