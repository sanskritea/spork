#import time

import numpy as np
import nidaqmx
from nidaqmx.stream_readers import CounterReader
from nidaqmx.constants import Edge, TriggerType, TaskMode, AcquisitionType, READ_ALL_AVAILABLE
from contextlib import ExitStack


class NIDAQ_PhotonCounter():
    def __init__(self):
        self.rawAPDchannel = '/Dev4/PFI3' #ctr1
        self.trigChannel = '/Dev4/PFI0' #ctr2
        self.clockChannel = '/Dev4/PFI5'#ctr3
        self.apdChannelsDict = {0: self.rawAPDchannel, 
                                1: self.clockChannel,
                                2: self.trigChannel,  
                            } 

    def __enter__(self):
        print('nidaq connected')
        return(self)

    def __exit__(self, exc_type, exc_value, exc_traceback):
        return

    def read_ctrs_many(self, acq_rate, numSamples:int, ctrChanNums=[0,1,2]):
        '''Reads specified counter channels for a given period, a designated number of times (based off internal timing)'''
        # print(f'ctrChanNums: '{ctrChanNums})
        print('Entering read_ctrs_many')
        print('ctrChanNums : ', ctrChanNums)
        ctrChannels = [self.apdChannelsDict[ctrChanNum] for ctrChanNum in ctrChanNums]
        print('ctrChannels : ', ctrChannels)
        print('acq_rate : ', acq_rate)
        period = 1 / acq_rate
        numSamples += 1 #so np.diff works later. Also, a DAQ buffer size of 1 isn't supported
        
        # list of np arrays containing raw counts for each counter
        # e.g. 3x ctr_channels with num_samples = 5:
        # [ np.array([5, 10, 20, 30, 35]), np.array([7, 7, 10, 11, 12]), np.array([0, 0, 2, 4, 5]) ]
        all_counts = []

        #create DAQ tasks
        with nidaqmx.Task() as dummyClkTask, ExitStack() as stack:
            # create a digital input dummy task to start the di/SampleClock@acqRate for clocking the counter input task

            dummyClkTask.di_channels.add_di_chan('Dev4/port0')
            dummyClkTask.timing.cfg_samp_clk_timing(
                acq_rate, 
                sample_mode=AcquisitionType.CONTINUOUS, 
                samps_per_chan=numSamples
            )
            dummyClkTask.control(TaskMode.TASK_COMMIT)

            # array containing counter reader streams
            readerStreams = []
            ctrTasks = []
            for i, ctrChan in enumerate(ctrChannels):
                # this automatically calls the __enter__ and __exit__
                # methods for each task as if they were used in a 'with' statement
                ctrTask = stack.enter_context(nidaqmx.Task())

                # create a counter input task
                ctrTask.ci_channels.add_ci_count_edges_chan(f'Dev4/ctr{i}') #connect to a ctr
                ctrTask.ci_channels.all.ci_count_edges_term = ctrChan #connect counter to relevant PFI channel

                # configure counter input task
                ctrTask.timing.cfg_samp_clk_timing(
                    acq_rate,
                    source='/Dev4/di/SampleClock',
                    active_edge=nidaqmx.constants.Edge.RISING,
                    sample_mode=nidaqmx.constants.AcquisitionType.FINITE,
                    samps_per_chan=numSamples
                )
                # create counter input stream object for later
                readerStreams.append(CounterReader(ctrTask.in_stream))

            # Load up tasks to be run quickly together later. Could also just ctrTask.start(), 
            # but worried about more substantial desyncing issues than just starting everything in a python for loop
            ctrTask.control(TaskMode.TASK_COMMIT) # TODO: Get a 2nd opinion about this from Michael
            ctrTasks.append(ctrTask)

            # start the counter tasks
            for ctrTask in ctrTasks:
                ctrTask.start()
            # start the internal clock
            dummyClkTask.start()

            # read counters out
            for readerStream in readerStreams:
                thisCtrRawCts = np.zeros(numSamples, dtype=np.uint32) + 1

                # read the counts out of the buffer
                readerStream.read_many_sample_uint32(thisCtrRawCts,
                                                number_of_samples_per_channel=nidaqmx.constants.READ_ALL_AVAILABLE,
                                                timeout=numSamples*period + 1) #1s overhead should be fine
                
                # calculate the difference in counts between each period
                all_counts.append(np.diff(thisCtrRawCts))
        #print(all_counts) #DEBUG
        return np.array(all_counts)


    def read_ctrs_single(self, acqRate, ctrChanNums):
        print('Entering read_ctrs_single')
        return(self.read_ctrs_many(acqRate, 1, ctrChanNums))


 ######################################################################################### SANSKRITI: FIX AFTER INTERNAL CLOCK COUNTING WORKS #############################################################################################   

    def read_ctrs_ext_clk(self, ctrChanNums=[0], num_samples, acq_rate):
        """Repeatedly read the number of pulses received in a given sampling period using a hardware-timed clock
        ctr_channels: iterable of counter channels to use, e.g. or ['Dev1/ctr1'] or ['Dev1/ctr1', 'Dev1/ctr2']
        period: time to read for each sample before moving onto the next sample (as a quantity object)
        num_samples: number of samples
        return: list of (np arrays of counts corresponding to each sample) for each counter channel
        """

        ctr_channels = [self.apdChannelsDict[ctrChanNum] for ctrChanNum in ctrChanNums]

        device_name = ctr_channels[0].split('/')[0] # e.g. "Dev1"
        # clk_channel = '/{0}/PFI0'.format(device_name)
        # using external trigger from Swabian but internal clock
        trig_channel = '/{0}/PFI5'.format(device_name)

        # list of np arrays containing raw counts for each counter
        # e.g. 3x ctr_channels with num_samples = 5:
        # [ np.array([5, 10, 20, 30, 35]), np.array([7, 7, 10, 11, 12]), np.array([0, 0, 2, 4, 5]) ]
        all_counts = []
        # create DAQ tasks
        with ExitStack() as stack:
            # array containing counter reader streams
            reader_streams = []
            # iterate through all the counter channels and perform setup
            for ctr_channel in ctr_channels:
                # this should automatically call the __enter__ and __exit__
                # methods for each task as if they were used in a 'with' statement
                counter_task = stack.enter_context(nidaqmx.Task())

                counter_task.ci_channels.add_ci_count_edges_chan(ctr_channel)
                # configure counter input task
                counter_task.timing.cfg_samp_clk_timing(
                    acq_rate,
                    source='/Dev4/di/SampleClock',
                    active_edge=nidaqmx.constants.Edge.RISING,
                    sample_mode=nidaqmx.constants.AcquisitionType.FINITE,
                    samps_per_chan=num_samples
                )
                
                #create trigger for reading to begin at the beginning of the sequence
                counter_task.triggers.arm_start_trigger.trig_type = TriggerType.DIGITAL_EDGE
                counter_task.triggers.arm_start_trigger.dig_edge_edge = Edge.RISING
                counter_task.triggers.arm_start_trigger.dig_edge_src = trig_channel

                # create counter input stream object
                sample_reader_stream = CounterReader(counter_task.in_stream)
                reader_streams.append(sample_reader_stream)
                # start the task
                counter_task.start()

            for sample_reader_stream in reader_streams:
                this_counter_raw_counts = np.zeros(num_samples, dtype=np.uint32)
                # read the counts out of the buffer
                # TODO check for timeout
                sample_reader_stream.read_many_sample_uint32(this_counter_raw_counts,
                                                number_of_samples_per_channel=nidaqmx.constants.READ_ALL_AVAILABLE,
                                                timeout=num_samples*(1/acq_rate).to('s').m + 1)
                # calculate the difference in counts between each period
                all_counts.append(this_counter_raw_counts)

        return all_counts




############################################################################################################################################################



###
### SANSKRITI: MODIFY THE PREVIOUS CODE TO ALIGN WITH THIS ONE BELOW ######

    # def readCtrs_singleRead_intClk(self, acqRate, ctrChanNums=[0,1,2,3]):
    #     return(self.readCtrs_multiRead_intClk(acqRate, 1, ctrChanNums,))
    
    """

    def readCtrs_multiRead_extTrigAndClk_basedOnBen(self, ctrChanNums=[1,4,8,11], numSamples=5, acqRate = 100e6):
        '''Reads specified counter channels a designated number of times based off external timing'''
        
        ctrChans = [self.apdChannelsDict[ctrChanNum] for ctrChanNum in ctrChanNums]

        # list of np arrays containing raw counts for each counter
        # e.g. 3x ctr_channels with num_samples = 5:
        # [ np.array([5, 10, 20, 30, 35]), np.array([7, 7, 10, 11, 12]), np.array([0, 0, 2, 4, 5]) ]
        all_counts = []

        # create DAQ tasks
        readerStreams = []
        for ctrChan in ctrChans:


            with nidaqmx.Task() as ctrTask:
                #create a counter input task
                ctrTask.ci_channels.add_ci_count_edges_chan(ctrChan)

                #configure the counter input task
                ctrTask.timing.cfg_samp_clk_timing(acqRate, samps_per_chan = numSamples)
                #source='' (for default onboard clock), active_edge=nidaqmx.constants.Edge.RISING, sample_mode=nidaqmx.constants.AcquisitionType.FINITE, samps_per_chan=numSamples)
                #TODO: source should be changed to clk_channel once we have one
                
                #create trigger for reading to begin at the beginning of the sequence
                ctrTask.triggers.arm_start_trigger.trig_type = TriggerType.DIGITAL_EDGE
                ctrTask.triggers.arm_start_trigger.dig_edge_edge = Edge.RISING
                ctrTask.triggers.arm_start_trigger.dig_edge_src = self.trigChannel
                #TODO: Look into setting read windows up as watchdog tasks that can be killed with some sort of expiration trigger

                #create counter input stream object for later
                readerStreams.append(CounterReader(ctrTask.in_stream))
                
                #start the counter task
                ctrTask.start()
            

            for readerStream in readerStreams:
                thisCtrRawCts = np.zeros(numSamples, dtype=np.uint32)

                # read the counts out of the buffer
                readerStream.read_many_sample_uint32(thisCtrRawCts,
                                                number_of_samples_per_channel=nidaqmx.constants.READ_ALL_AVAILABLE,
                                                timeout=numSamples*(1/acqRate) + 1) #1s overhead should be fine
                
                # calculate the difference in counts between each period
                all_counts.append(np.append(np.zeros(1), thisCtrRawCts)) #add a leading 0 so diff works right

        return all_counts"""


##################### SANSKRITI: IGNORE THIS PART ####################

    """def readCtrs_multiRead_extTrig_MoreMine(self, ctrChanNums=[1,4,8,11], numSamples=5, acqRate = 100e6):
        '''Reads specified counter channels a designated number of times based off external timing'''
        
        ctrChans = [self.apdChannelsDict[ctrChanNum] for ctrChanNum in ctrChanNums]

        # list of np arrays containing raw counts for each counter
        # e.g. 3x ctr_channels with num_samples = 5:
        # [ np.array([5, 10, 20, 30, 35]), np.array([7, 7, 10, 11, 12]), np.array([0, 0, 2, 4, 5]) ]
        all_counts = []

        # create DAQ tasks
        readerStreams = []
        for ctrChan in ctrChans:
            with nidaqmx.Task() as ctrTask:
                #create a counter input task
                ctrTask.ci_channels.add_ci_count_edges_chan(ctrChan)

                #configure the counter input task
                ctrTask.timing.cfg_samp_clk_timing(acqRate, samps_per_chan = numSamples)
                #source='' (for default onboard clock), active_edge=nidaqmx.constants.Edge.RISING, sample_mode=nidaqmx.constants.AcquisitionType.FINITE, samps_per_chan=numSamples)
                
                #create trigger for reading to begin at the beginning of the sequence
                ctrTask.triggers.arm_start_trigger.trig_type = TriggerType.DIGITAL_EDGE
                ctrTask.triggers.arm_start_trigger.dig_edge_edge = Edge.RISING
                ctrTask.triggers.arm_start_trigger.dig_edge_src = self.trigChannel 
                #TODO: Look into setting read windows up as watchdog tasks that can be killed with some sort of expiration trigger

                #create counter input stream object for later
                readerStreams.append(CounterReader(ctrTask.in_stream))
                
                #start the counter task
                ctrTask.start()"""
    
if __name__ == '__main__':
    daq = NIDAQ_PhotonCounter()
    print(daq.readCtrs_singleRead_intClk(acqRate=1)) #1s-long readout
