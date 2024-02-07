# -*- coding: utf-8 -*-
'''
	Driver for Vaunix Sig Gen LSF402
	D.Mark - Mar 2023
	De-lantzed version of Nazar's previous driver. Use as you wish, be grey.

	N.Delegan - Sep 2020.
    Use as you wish, be nice.

'''

from ctypes import c_int, byref, create_string_buffer, WinDLL


class SignalGeneratorLSG402:
	'''Driver class for any Vaunix LSG402 signal generator.
    '''
    # Definition of constants and markers
	freqStep = 100000  # Hz
	pwrStep = 0.25  # dBm
	sweepModeMarker = 1
	sweepDirMarker = 1

	# Definition of boundaries (these are then dynamically determined in the initialize function)
	#note this accounts for the 3 dB loss , so upper cap is -0.5 dBm which is safe
	#Switch Model: ZASWA-2-50DR+
	minFreq = 10000 * freqStep  # Hz
	maxFreq = 40000 * freqStep  # Hz
	minPwr = -42 * pwrStep  # dBm
	maxPwr = 10 * pwrStep  # dBm

	dev_status_values = {
		'Invalid devID' : 0x80000000,
		'Dev connected' : 0x00000001,
		'Dev opened'    : 0x00000002,
		'Dev is sweeping': 0x00000004,
		'Dev is sweeping up in freq.':0x00000008,
		'Dev in continuous sweep mode':0x00000010,
		'Dev in bi-directional sweep mode':0x00000020,
		'Both PLLs are locked':0x00000040,
		"UNKNOWN #1 (Device Initialized?)": 16395,
		"UNKNOWN #2 (Device Closed?)": 9,
		"UNKNOWN #3 (Device in test mode?)" : 195,
		"UNKNOWN #4 (Device in test mode but closed?)" : 193,
		"UNKNOWN #5 (Device RF-ON NO SWEEP)" : 16459,
		"UNKNOWN #6" : 0,
	}

	def __init__(self):
		# Declare where the library path is.
		library_path = r'C:\Users\adminnd\code\jasper-repo\drivers\vaunix\Libraries\vnx_fsynth.dll'
		try:
			self.lib = WinDLL(library_path)
		except:
			print('\nUnable to load Vaunix library at {0}.\n'.format(library_path))

		# Proceed with driver
		# Firstly, define test or operation mode
		# Set to true for testing, will create a virtual device, not connec to real one.
		TestMode = False
		self.lib.fnLSG_SetTestMode(TestMode)

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

	def __enter__(self):
		'''Initialization function'''
		# Update limits
		self.minFreq = self.lib.fnLSG_GetMinFreq(
		    self.devID)*self.freqStep  # in Hz
		self.maxFreq = self.lib.fnLSG_GetMaxFreq(
		    self.devID)*self.freqStep  # in Hz
		self.minPwr = self.lib.fnLSG_GetMinPwr(
		    self.devID)*self.pwrStep  # in dBm
		self.maxPwr = self.lib.fnLSG_GetMaxPwr(
		    self.devID)*self.pwrStep  # in dBm
		return self

	def __exit__(self, *args):
		'''Finalization function used to disconnect device.'''
		# Device disconnect
		# TODO: add checking and exception handling.
		self.lib.fnLSG_CloseDevice(self.devID)
		return

'''
Self-identification returns a string with identifying features of the hardware

Calls and combined output from the following library functions:
*GetDLLVersion() with output being the .dll version (hex)
*GetSerialNumber() with output being the serial number (integer)
*GetModelName() with output being the model name placed in a byref entity
'''
	def id(self):
		sn = self.lib.fnLSG_GetSerialNumber(self.devID)
		dll_ver = self.lib.fnLSG_GetDLLVersion()
		model_name = create_string_buffer(12)
		self.lib.fnLSG_GetModelNameA(self.devID, byref(model_name))
		return 'Vaunix ' + model_name.value.decode('utf-8') + ', serial number: ' + str(sn) + '. DLL version: ' + str(dll_ver)

	def status(self):
		return self.lib.fnLSG_GetSeviceStatus(self.devID)

#FREQUENCY COMMANDS

	#Frequency return and set functions
	'''get and set frequency functions in Hz'''
	def getfreq(self):
		return self.lib.fnLSG_GetFrequency(self.devID)*self.freqStep
	
	def setFreq(self, val):
		self.lib.fnLSG_SetFrequency(self.devID, int(val/self.freqStep))
	
	#Reference frequency
	ref_frequency_values ={
		'Internal ref. freq.': 1,
		'External ref. freq.': 0,
		'Device not ready': 0x80030000,
	}

	'''get/set internal or external frequency reference clock'''
	def getIntRef(self):
		return self.lib.fnLSG_GetUseInternalRef(self.devID)


	def setIntRef(self, val):
		self.lib.fnLSG_SetUseInternalRef(self.devID,bool(val))

	'''return/define start frequency for sweeps in Hz'''
	def getStartF(self):
		return self.lib.fnLSG_GetStartFrequency(self.devID)*self.freqStep

	def setStartF(self, val):
		self.lib.fnLSG_SetStartFrequency(self.devID, int(val/self.freqStep))

	'''return/define stop frequency for sweeps in Hz'''
	def getEndF(self):
		return self.lib.fnLSG_GetEndFrequency(self.devID)*self.freqStep

	def setEndFreq(self, val):
		self.lib.fnLSG_SetEndFrequency(self.devID, int(val/self.freqStep))
	
	'''get/set frequency sweep step in Hz'''
	def getFStep(self):
		return self.lib.fnLSG_GetFrequencyStep(self.devID)*self.freqStep

	def setFStep(self, val):
		self.lib.fnLSG_SetFrequencyStep(self.devID, int(val/self.freqStep))

	'''get/set dwell time (s) at each step'''
	def getFDwell(self):
		return self.lib.fnLSG_GetDwellTime(self.devID)*1000
	
	def setFDwell(self, val):
		self.lib.fnLSG_SetDwellTime(self.devID, int(val*1000))	



	'''Sweep Mode functions. 
		Note: a return/input of 1 = 'Repeating Sweep' and 0 = 'Single Sweep' '''
	def sweepMode(self):
		return self.sweepModeMarker

	def setSweepMode(self, value):
		self.lib.fnLSG_SetSweepMode(self.devID,value)
		self.sweepModeMarker = value
	

	'''
	sweep up = 1, sweep down = 0'''
	def sweepDir(self):
		return self.sweepDirMarker

	def setSweepDir(self, val):
		self.lib.fnLSG_SetSweepDirection(self.devID, val)

	'''get/set output power in dBm'''
	def getPwrLvl(self):
		return self.lib.fnLSG_GetPowerLevelAbs(self.devID)*self.pwrStep

	def setPwrLvl(self, val):
		self.lib.fnLSG_SetPowerLevel(self.devID, int(val/self.pwrStep))

	''' 'RF Output ON': True, 'RF Output OFF': False, 'Device not ready': 0x80030000 '''
	def rfOut(self):
		return self.lib.fnLSG_GetRF_On(self.devID)

	def setRFOut(self, val):
		self.lib.fnLSG_SetRFOn(self.devID, val)

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
		print(SigGen.id())
		print('\nIdentification successfull. Driver seems to work.\n')
