''' 

NI-DAQ 6343 USB driver

Sanskriti Chitransh, 2023-Sept-09

Adapted from Evan Villafranca, Jacob Feder

Note: this DAQ doesn't live on the instrument server. The driver gets called whenever needed purely as a python file to communicate with the DAQ


'''


import numpy as np
import time
import math
from itertools import cycle
import logging
import scipy as sp
from scipy import signal
import datetime as Dt

from rpyc.utils.classic import obtain 

# nidaqmx STAYS SAME
import nidaqmx

# CHECK FOR ADDITIONAL CONSTANTS
from nidaqmx.constants import Edge, READ_ALL_AVAILABLE, TaskMode, TriggerType, AcquisitionType, WAIT_INFINITELY
from nidaqmx.stream_readers import CounterReader
from nidaqmx.stream_writers import AnalogSingleChannelWriter 


class NIDAQ():

    def __init__(self):

        self.read_task = None
        self.clock_task = None
        self.reader_stream = None
        
        self.clk_channel = '/Dev4/PFI3'     # external clock from the Swabian
        self.dev_channel = '/Dev4/PFI0'     # data TTL clicks from the APD


    def __enter__(self):

        self.read_task = None
        self.clock_task = None
        self.reader_stream = None

        return self


    def start_external_read_task(self, sampling_rate, num_samples):
    # def start_read_task(self, num_samples):     # create read task, set up counter, source clock and reader stream.

        # creating DAQ task  
        self.read_task = nidaqmx.Task()
        
        # adding digital input channel as counter (APD)
        self.read_task.ci_channels.add_ci_count_edges_chan(f'/Dev4/ctr1')

        # connecting physical input to virtual counter
        self.read_task.ci_channels.all.ci_count_edges_term = self.dev_channel

        # setting up timing clock (external)
        self.read_task.timing.cfg_samp_clk_timing(
                                sampling_rate,   # max DAQ sampling rate for convenience
                                source = self.clk_channel, # Swabian clock ticks
                                active_edge = nidaqmx.constants.Edge.RISING,
                                sample_mode = nidaqmx.constants.AcquisitionType.FINITE,
                                samps_per_chan = num_samples # number of ticks to be collected
        )



        # starting counting task
        self.read_task.start()

        return self.read_task


    def external_read_task(self, sampling_rate, num_samples):

        # creating counter stream object for counting
        self.reader_stream = CounterReader(self.read_task.in_stream)

        # reading out the counter
        raw_counts = np.zeros(num_samples, dtype=np.uint32) 
        self.reader_stream.read_many_sample_uint32(
                                raw_counts,
                                number_of_samples_per_channel = nidaqmx.constants.READ_ALL_AVAILABLE,
                                timeout = nidaqmx.constants.WAIT_INFINITELY
        )

        # stop read task after each read 
        self.read_task.stop()

        # close read task after each read
        self.read_task.close()
        self.read_task = None
        self.reader_stream = None

        # calculate difference in counts between consecutive clock ticks
        counts = np.diff(raw_counts)
        return counts


    def internal_read_task(self, clock_rate, num_samples):

        # setup fake internal clock task
        self.clock_task = nidaqmx.Task()
        self.clock_task.di_channels.add_di_chan(f'/Dev4/port0')
        self.clock_task.timing.cfg_samp_clk_timing(
                                clock_rate,   # max DAQ clock rate
                                sample_mode = nidaqmx.constants.AcquisitionType.CONTINUOUS,
        )
        self.clock_task.control(TaskMode.TASK_COMMIT)

        # creating DAQ reading task  
        self.read_task = nidaqmx.Task()
        self.read_task.ci_channels.add_ci_count_edges_chan(f'/Dev4/ctr1')
        self.read_task.ci_channels.all.ci_count_edges_term = self.dev_channel
        self.read_task.timing.cfg_samp_clk_timing(
                                clock_rate,   # max DAQ sampling rate for convenience
                                source = f'/Dev4/di/SampleClock', 
                                active_edge = nidaqmx.constants.Edge.RISING,
                                sample_mode = nidaqmx.constants.AcquisitionType.CONTINUOUS,
        )
        
        # set buffer size (how much data you can fit): 
        self.read_task.in_stream.input_buf_size = num_samples

        # Triggers the read off of the sample clock: 
        self.read_task.triggers.arm_start_trigger.trig_type = TriggerType.DIGITAL_EDGE
        self.read_task.triggers.arm_start_trigger.dig_edge_edge = Edge.RISING
        self.read_task.triggers.arm_start_trigger.dig_edge_src = f'/Dev4/di/SampleClock'

        # setup reaser stream
        self.reader_stream = CounterReader(self.read_task.in_stream)

        # start tasks
        self.clock_task.start()
        self.read_task.start()

        # reading out the counter
        raw_counts = np.zeros(num_samples, dtype=np.uint32) 
        self.reader_stream.read_many_sample_uint32(
                                raw_counts,
                                number_of_samples_per_channel = num_samples,
                                timeout = nidaqmx.constants.WAIT_INFINITELY
        )

        # stop tasks after each read 
        self.read_task.stop()
        self.clock_task.stop()

        # close tasks after each read
        self.read_task.close()
        self.clock_task.close()
        self.read_task = None
        self.clock_task = None
        self.reader_stream = None

        # calculate difference in counts between consecutive clock ticks
        counts = np.diff(raw_counts)

        return counts


    def __exit__(self, *args):  # execute when the DAQ object gets killed unexpectedly, for example, when the STOP button is clicked.

        print('in exit')
        return
        # if self.read_task != None:  # in case the DAQ object was killed before the reading was over, close and destroy all read tasks
        #     self.read_task.stop()
        #     self.read_task.close()
        #     self.read_task = None
        #     self.reader_stream = None
        # else:
        #     return


    def laser_power_atten(self, atten_voltage):
        # use the DAQ analog output to control attenuator voltage, therefore, laser power. DAQ AO between 0-5V.

        # creating voltage output task
        self.write_task = nidaqmx.Task()

        # adding analog output channel
        self.write_task.ao_channels.add_ao_voltage_chan(f'/Dev4/ao0', min_val = 0, max_val = 5)

        # creating analog writer object for voltage output
        self.writer_stream = AnalogSingleChannelWriter(self.write_task.out_stream)

        # output attenuator voltage
        self.writer_stream.write_one_sample(np.array(atten_voltage))

        # closeout
        self.write_task.stop()
        self.write_task.close()
        self.write_task = None


    def create_task(self):

        # create NIDAQ task
        self.read_task = nidaqmx.Task()
        return self.read_task


    def start_task(self):
        
        # starting counting task
        self.read_task.start()


    def stop_task(self):

        # stop_task
        self.read_task.stop()


    def close_task(self):

        # close_task
        self.read_task.close()
        

    

#########################################################################################################################################################################
#########################################################################################################################################################################
   

        


    
