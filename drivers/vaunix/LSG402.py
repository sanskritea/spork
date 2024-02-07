# -*- coding: utf-8 -*-
"""
    drivers.vaunix.LSG402
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	:copyright: 
	D.Mark - Feb 2023., C. Egerstrom April 2023
    Use as you wish, be nice.

    N.Delegan - Sep 2020.
    Use as you wish, be nice.

    :description:
    Driver for Vaunix LSG402
    
    :versions/changelog: 
    * V1.0 - Sep 2020 -  Driver tested and deemed not yet functionnal. 
	Driver limited to one Vaunix device at a time right now.
	* V2.0 - Feb 2023 -  Driver de-Lantz-ed 

    :dependencies: 
    * ctypes c_int, byref, create_string_buffer, WinDLL
"""

from ctypes import c_int, byref, create_string_buffer, WinDLL


class SignalGeneratorLSG402:
	'''Driver class for any Vaunix LSG402 signal generator.
    '''
    # Definition of constants and markers
	FREQ_STEP = 100000  # Hz
	POWER_STEP = 0.25  # dBm
	sweepModeMarker = 1
	sweepDirMarker = 1

	# Definition of boundaries (these are then dynamically determined in the initialize function)
	HARDWARE_MIN_FREQ = 10000 * FREQ_STEP  # Hz
	HARDWARE_MAX_FREQ = 40000 * FREQ_STEP  # Hz
	HARDWARE_MIN_POWER = -42 * POWER_STEP  # dBm
	HARDWARE_MAX_POWER = 10 * POWER_STEP  # dBm
	SOFTWARE_MAX_POWER = 0 * POWER_STEP #dBm, rough estimate from link budget

	def __init__(self, runTestMode = False):
		'''Load relevant DLLs for controlling Vaunix Sig-gen.
		Arguments:  *runTestMode, sets whether it will connect try to connect to a real
					sig-gen or create a virtual device. Default = False (i.e. real sig-gen)'''
		# Declare where the library path is.
		library_path = r'C:\Users\adminnd\code\jasper-repo\drivers\vaunix\Libraries\vnx_fsynth.dll'
		try:
			self.lib = WinDLL(library_path)
		except:
			print('\nUnable to load Vaunix library at {0}.\n'.format(library_path))

		# Proceed with driver
		# Firstly, define test or operation mode
		# Set to true for testing, will create a virtual device, not connec to real one.
		self.lib.fnLSG_SetTestMode(runTestMode)


	def __enter__(self):
		'''Initialization function. Connects to first sig-gen in list (so will break if there are >1 sig-gen).
		 Also updates hardware limits (in hardware units) for freq and power '''
		# Find devices and count them
		# This is important, finds the devices.
		dev_number = self.lib.fnLSG_GetNumDevices()
		print('{0} devices were recognized by the Vaunix API.'.format(dev_number))

		# Create an IDArray for the device identifiers (currently set to 5 devices)
		DeviceIDArray = c_int * 5  # Create a c_int array for device IDs
		Devices = DeviceIDArray()  # Declare deviced in a format that the DLL can understand
		# Get list of active device handles and place them into Devices
		self.lib.fnLSG_GetDevInfo(Devices)

		# Get device ID (first device in the list)
		# TODO: For more than one vaunix, rewrite the next bit of code, right now it only considers the first one.
		# First device on list is kept (code only written for one right now)
		self.devID = Devices[0]

		# Device initialization
		self.lib.fnLSG_InitDevice(self.devID)

		# Update limits
		self.HARDWARE_MIN_FREQ = self.lib.fnLSG_GetMinFreq(
		    self.devID)*self.FREQ_STEP  # in Hz
		self.HARDWARE_MAX_FREQ = self.lib.fnLSG_GetMaxFreq(
		    self.devID)*self.FREQ_STEP  # in Hz
		self.HARDWARE_MIN_POWER = self.lib.fnLSG_GetMinPwr(
		    self.devID)*self.POWER_STEP  # in dBm
		self.HARDWARE_MAX_POWER = self.lib.fnLSG_GetMaxPwr(
		    self.devID)*self.POWER_STEP  # in dBm
		
		return self


	def __exit__(self, *args):
		'''Finalize function used to disconnect device.'''
		# Device disconnect
		# TODO: add checking and exception handling.
		self.lib.fnLSG_CloseDevice(self.devID)


	def idn(self):
		'''Returns a verbose string about the connected Vaunix sig-gen'''
		sn = self.lib.fnLSG_GetSerialNumber(self.devID)
		dll_ver = self.lib.fnLSG_GetDLLVersion()
		model_name = create_string_buffer(12)
		self.lib.fnLSG_GetModelNameA(self.devID, byref(model_name))
		return 'Vaunix ' + model_name.value.decode('utf-8') + ', serial number: ' + str(sn) + '. DLL version: ' + str(dll_ver)


	def getStatus(self):
		'''Returns the status of the Vaunix sig-gen'''
		return self.lib.fnLSG_GetSeviceStatus(self.devID)


	#FREQUENCY COMMANDS

	#Frequency return and set functions
	def getFreq(self):
		'''Returns the current frequency of the sig-gen (in Hz)'''
		return self.lib.fnLSG_GetFrequency(self.devID)*self.FREQ_STEP
	

	def setFreq(self, val):
		'''Set the frequency of the sig-gen 
		Arguments:  *val, new frequency of the sig-gen (in Hz)'''
		self.lib.fnLSG_SetFrequency(self.devID, int(val/self.FREQ_STEP))


	#obtain internal reference frequency
	def getIntRef(self):
		'''Returns the internal reference frequency of the sig-gen (in Hz)'''
		return self.lib.fnLSG_GetUseInternalRef(self.devID)


	#set internal reference frequency
	def setIntRef(self, val):
		'''Set the internal reference frequency of the sig-gen 
		Arguments:  *val, new reference frequency of the sig-gen (in Hz)'''
		self.lib.fnLSG_SetUseInternalRef(self.devID,bool(val))


	#TODO: Comment getters and setters below once you read the docs and figure out what they are
	#get starting frequency
	def getStartFreq(self):
		return self.lib.fnLSG_GetStartFrequency(self.devID)*self.FREQ_STEP

	#get end frequency
	def getEndFreq(self):
		return self.lib.fnLSG_GetEndFrequency(self.devID)*self.FREQ_STEP

	#set starting frequency
	def setStartF(self, val):
		self.lib.fnLSG_SetStartFrequency(self.devID, int(val/self.FREQ_STEP))

	#set ending frequency
	def setEndFreq(self, val):
		self.lib.fnLSG_SetEndFrequency(self.devID, int(val/self.FREQ_STEP))
	
	def getFreqStep(self):
		return self.lib.fnLSG_GetFrequencyStep(self.devID)*self.FREQ_STEP

	def setFreqStep(self, val):
		self.lib.fnLSG_SetFrequencyStep(self.devID, int(val/self.FREQ_STEP))

	def getFreqDwell(self):
		return self.lib.fnLSG_GetDwellTime(self.devID)*1000
	
	def setFreqDwell(self, val):
		self.lib.fnLSG_SetDwellTime(self.devID, int(val*1000))	


#Sweep Mode functions. Note a return/input of 1 = 'Repeating Sweep' and 0 = 'Single Sweep'
	def getSweepMode(self):
		'''Returns the sweep mode. A return/input of 1 = "Repeating Sweep" and 0 = "Single Sweep" '''
		return self.sweepModeMarker

	def setSweepMode(self, value):
		'''Sets the sweep mode
		Arguments:  *value, input of 1 = "Repeating Sweep" and 0 = "Single Sweep" '''
		self.lib.fnLSG_SetSweepMode(self.devID,value)
		self.sweepModeMarker = value
	

	#TODO: What exactly do these mean? Write better help strings when you figure it out
	#sweep up = 1, sweep down = 0 
	def getSweepDir(self):
		'''Gets the sweep direction. Sweep up = 1, Sweep down = 0'''
		return self.sweepDirMarker

	def setSweepDir(self, val):
		'''Sets the sweep direction
		Arguments:  *val, input of 1=Sweep up, 0=Sweep down'''
		self.lib.fnLSG_SetSweepDirection(self.devID, val)


	#RF Power Setting Methods 
	def getPwrLvl(self):
		'''Returns the RF power setpoint of the sig-gen (in dBm)'''
		return self.lib.fnLSG_GetPowerLevelAbs(self.devID)*self.POWER_STEP

	def setPwrLvl(self, val):
		'''Sets the RF power setpoint of the sig-gen
		Arguments:  *val, RF power to set sig-gen to (in dBm)'''
		if val <= self.SOFTWARE_MAX_POWER:
			self.lib.fnLSG_SetPowerLevel(self.devID, int(val/self.POWER_STEP))
		else:
			raise ValueError('Warning: This RF Power out of the Vaunix could damage RF components post-amp')


	def getIfRfOut(self):
		'''Returns if the sig-gen is actually outputting RF power
		'RF Output ON'=> True, 'RF Output OFF'=> False, 'Device not ready'=> 0x80030000'''
		return self.lib.fnLSG_GetRF_On(self.devID)

	def setIfRFOut(self, val):
		'''Sets if the sig-gen should output RF power
		True =>"RF Output ON", False => "RF Output OFF"'''
		self.lib.fnLSG_SetRFOn(self.devID, val)


	#Other methods
	def sweepStart(self):
		self.lib.fnLSG_StarSweep(self.devID, True)

	def saveSettings(self):
		self.lib.fnLSG_SaveSettings(self.devID)



#TESTING PURPOSES
if __name__ == '__main__':
	'''Solo call testing function '''
	print('Testing the Vaunix LSG402 signal generator driver... \n')
	with SignalGeneratorLSG402() as SigGen:
		print('Asking the LSG402 to self-identify:\n')
		print(SigGen.idn())
		print('\nIdentification successfull. Driver seems to work.\n')
