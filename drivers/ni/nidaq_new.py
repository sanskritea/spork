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

#CHECK FOR ADDITIONAL CONSTANTS
from nidaqmx.constants import Edge, READ_ALL_AVAILABLE, TaskMode, TriggerType, AcquisitionType
from nidaqmx.stream_readers import CounterReader
from nidaqmx.stream_readers import AnalogSingleChannelReader 

from contextlib import ExitStack

class NIDAQ():

    def __init__(self):
        pass
        
        # USB 6343
        # self.sampling_rate = 50e3
        # self.period = 1 / self.sampling_rate
        # self.read_task = nidaqmx.Task()
        # self.reader_stream = []


    def int_clk_read_many(self, acq_rate, num_samples):

        # empty counts and stream
        counts = []
        reader_stream = []

        # applying user's acquisition rate
        sampling_rate = acq_rate
        period = 1 / sampling_rate

        # creating DAQ tasks
        with nidaqmx.Task() as dummy, ExitStack() as stack:
            # creating dummy signal for internal clock
            dummy = nidaqmx.Task()
            
            dummy.di_channels.add_di_chan('Dev4/port0')
            dummy.timing.cfg_samp_clk_timing(
                                    acq_rate, 
                                    sample_mode=AcquisitionType.CONTINUOUS
            )
            dummy.control(TaskMode.TASK_COMMIT)

            # creating counter task
            read_task = nidaqmx.Task()
    
            # adding digital input channel as counter (APD)
            read_task.ci_channels.add_ci_count_edges_chan(f'Dev4/ctr0')

            # connecting physical input to virtual counter
            read_task.ci_channels.all.ci_count_edges_term = '/Dev4/PFI0'

            # setting up internal timing clock 
            read_task.timing.cfg_samp_clk_timing(
                                    sampling_rate,
                                    source = '/Dev4/di/SampleClock',
                                    active_edge = nidaqmx.constants.Edge.RISING,
                                    sample_mode = nidaqmx.constants.AcquisitionType.FINITE,
                                    samps_per_chan = num_samples
            )

            # creating counter stream object for counting
            reader_stream.append(CounterReader(read_task.in_stream))

            # loading up tasks
            read_task.control(TaskMode.TASK_COMMIT)

        # starting counting and dummy tasks
        dummy.start()
        read_task.start()

        buffer = np.zeros(num_samples, dtype=np.uint32) + 1
        reader_stream[0].read_many_sample_uint32(
                                    buffer,
                                    number_of_samples_per_channel = nidaqmx.constants.READ_ALL_AVAILABLE,
                                    timeout = num_samples * period + 1
            )

        # calculate difference in counts between each period
        counts.append(np.diff(buffer))

        # return counts
        return np.array(counts)

        # reading counter out starting with the buffer
        # try:
        #     buffer = np.zeros(num_samples, dtype=np.uint32) + 1
        #     # print(type(obtain(buffer)))
        #     # print(buffer.dtype)
        #     self.reader_stream.read_many_sample_uint32(
        #                             buffer,
        #                             number_of_samples_per_channel = nidaqmx.constants.READ_ALL_AVAILABLE,
        #                             timeout = num_samples * self.period + 1
        #     )
            
        #     # calculate difference in counts between each period
        #     counts.append(np.diff(buffer))

        #     # return counts
        #     return np.array(counts)
            
        # except TypeError:
        #     self.stop_task()
        #     self.close_task()


    def ext_clk_read_many(self, acq_rate, num_samples):
        
        self.clk_channel = 'Dev4/PFI3' # ctr 1
        self.dev_channel = 'Dev4/PFI5' # ctr 3

        # applying user's acquisition rate
        self.sampling_rate = acq_rate
        self.period = 1 / self.sampling_rate

        # empty counts
        counts = []

        # creating DAQ task  
        stack = ExitStack()
        self.read_task = stack.enter_context(nidaqmx.Task())
        
        # adding diginal input channel (APD)
        self.read_task.ci_channels.add_ci_count_edges_chan(self.dev_channel)

        # setting up internal timing clock 
        self.read_task.timing.cfg_samp_clk_timing(
                                self.sampling_rate,
                                source = self.clk_channel,
                                active_edge = nidaqmx.constants.Edge.RISING,
                                sample_mode = nidaqmx.constants.AcquisitionType.FINITE,
                                samps_per_chan = num_samples
        )

        # creating counter stream object for counting
        self.reader.append(CounterReader(self.read_task.in_stream))
        
        # starting counting task
        self.read_task.start()

        # reading counter out starting with the buffer
        try:
            buffer = obtain(buffer)
            print("DAQ received empty buffer. Reading samples...")
            # print(type(obtain(buffer)))
            # print(buffer.dtype)
            self.reader.read_many_sample_uint32(
                                    buffer,
                                    number_of_samples_per_channel = nidaqmx.constants.READ_ALL_AVAILABLE,
                                    timeout = num_samples * self.period + 1
            )
            
            # calculate difference in counts between each period
            counts.append(np.diff(buffer))

            # return counts
            return np.array(counts)
            
        except TypeError:
            self.stop_task()
            self.close_task()


    def single_read(self, acqRate):
        print('Entering read_ctrs_single')
        return(self.int_clk_read_many(acqRate, 2))


    def stop_task(self):  
        print(f"{self.read_task} TASK STOPPED")
        self.read_task.stop()


    def close_task(self):
        print(f"{self.read_task} TASK CLOSED")
        self.read_task.close()
        self.read_task = None
        self.reader = None
    
