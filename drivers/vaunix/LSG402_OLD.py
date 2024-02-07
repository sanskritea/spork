# -*- coding: utf-8 -*-
"""
    lantz.drivers.vaunix.LSG402
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	:copyright: 
    N.Delegan - Sep 2020.
    Use as you wish, be nice.

    :description:
    Driver for Vaunix LSG402
    
    :versions/changelog: 
    * V1.0 - Sep 2020 -  Driver tested and deemed not yet functionnal. 
	Driver limited to one Vaunix device at a time right now.

    :dependencies: 
    * ctypes c_int, byref, create_string_buffer
    * lantz Action, Feat,  LibraryDriver
	* vnx_fsynth.dll
"""
from lantz import Action, Feat
from lantz.core.foreign import LibraryDriver
from ctypes import c_int, byref, create_string_buffer

#TODO: get the VNX_fsynth.h file for the header info on the error codes
class SignalGeneratorLSG402(LibraryDriver):
	"""Driver class for any Vaunix LSG402 signal generator.
    """
	LIBRARY_NAME = "vnx_fsynth.dll"
	LIBRARY_PREFIX = "fnLSG_"

	def __init__(self):
		super(SignalGeneratorLSG402, self).__init__()
		TestMode=False
		self.lib.SetTestMode(TestMode) #Set to True for testing
		dev_number = self.lib.GetNumDevices() #This is important, finds the devices.
		print(dev_number + " devices were recognized by the Vaunix API.")

		DeviceIDArray = c_int * 5 #Create a c_int array for device IDs
		Devices = DeviceIDArray() #Declare deviced in a format that the DLL can understand
		self.lib.GetDevInfo(Devices) #Get list of active device handles and place them into Devices

		#TODO: For more than one vaunix, rewrite the next bit of code, right now it only considers the first one.
		self.dev_ID = Devices[0] #First device on list is kept (code only written for one right now)

		#Device inialization
		self.lib.InitDevice(self.dev_ID)
		return

	#Definition of constants and markers
	freq_step = 100000 #Hz
	power_step = 0.25 #dBm
	sweep_mode_marker = 1
	sweep_direction_marker = 1

	#Definition of boundaries (these are then dynamically determined in the initialize function)
	min_freq = 10000 * freq_step #Hz
	max_freq = 40000 * freq_step #Hz
	min_power = -42 * power_step #dBm
	max_power = 10 * power_step #dBm

	def initialize(self):
		"""Initialization function"""
		#Update limits
		self.min_freq = self.lib.GetMinFreq(self.dev_ID)*self.freq_step #in Hz
		self.max_freq = self.lib.GetMaxFreq(self.dev_ID)*self.freq_step #in Hz
		self.min_power = self.lib.GetMinPwr(self.dev_ID)*self.power_step #in dBm
		self.max_power = self.lib.GetMaxPwr(self.dev_ID)*self.power_step #in dBm

		#TODO: check what following stuff is, Michael's sophistry
		#MISSING = _NamedObject('MISSING')
		# print('aoidsfhoasidhfoiashfoidhsafoi')
		# print(MISSING.name)
		# print(self.frequency)
		#print(self.frequency.modifiers) #['limits'] = (self.min_freq, self.max_freq, self.freq_step)
		return

	def finalize(self):
		"""Finalization function used to disconnect device."""
		#Device disconnect
		#TODO: add checking and exception handling.
		self.lib.CloseDevice(self.dev_ID)
		return

	#FEATURES
	@Feat()
	def idn(self):
		"""Self-identification returns a string with identifying features of the hardware

		Calls and combined output from the following library functions:
		*GetDLLVersion() with output being the .dll version (hex)
		*GetSerialNumber() with output being the serial number (integer)
		*GetModelName() with output being the model name placed in a byref entity
		"""
		#Get serial number, DLL version, and model name
		sn = self.lib.GetSerialNumber(self.dev_ID)

		dll_ver = self.lib.GetDLLVersion()
		dll_ver = int(dll_ver, 16)

		model_name = create_string_buffer(12)
		self.lib.GetModelNameA(self.dev_ID,byref(model_name))
		return "Vaunix " + model_name.value.decode('utf-8') + ", serial number: " + str(sn) + ". DLL version: " + str(dll_ver)

	#TODO figure out what 16395,16459,16387, 195, and 9 are (email the company) 
	@Feat(values={"Invalid DeviceID": -2147483648, 
					"Device Connected": 1,
					"Device Opened": 2,
					"Device is sweeping": 4,
					"Device is sweeping up in freq.": 8,
					"Device in continuous sweep mode": 16,
					"Device in bi-directional sweep mode": 32,
					"PLL lock status is True": 64,
					"UNKNOWN #1 (Device Initialized?)": 16395,
					"UNKNOWN #2 (Device Closed?)": 9,
					"UNKNOWN #3 (Device in test mode?)" : 195,
					"UNKNOWN #4 (Device in test mode but closed?)" : 193,
					"UNKNOWN #5 (Device RF-ON NO SWEEP)" : 16459,
					"UNKNOWN #6" : 0}) 
	def device_status(self):
		"""Get device status function"""
		return self.lib.GetDeviceStatus(self.dev_ID) 

	#Device frequency
	@Feat(units='Hz', limits=(min_freq,max_freq,freq_step))
	def frequency(self):
		"""Frequency get and set function, limits and units (Hz) defined on the Feat level"""
		return self.lib.GetFrequency(self.dev_ID)*self.freq_step
	@frequency.setter
	def frequency(self, value):
		self.lib.SetFrequency(self.dev_ID, int(value/self.freq_step))

	#Reference frequency
	@Feat(values={"Internal ref. freq.": 1,
				 "External ref. freq.": 0,
				 "Device not ready": -2147287040})
	def RF_internal_ref(self):
		return self.lib.GetUseInternalRef(self.dev_ID)

	@RF_internal_ref.setter
	def RF_internal_ref(self, value):
		self.lib.SetUseInternalRef(self.dev_ID,bool(value))

	#Frequency start and end
	@Feat(units="Hz", limits=(min_freq,max_freq,freq_step))
	def frequency_start(self):
		#Frequency readout with limits defined by the hardware.
		return self.lib.GetStartFrequency(self.dev_ID)*self.freq_step

	@frequency_start.setter
	def frequency_start(self, value):
		self.lib.SetStartFrequency(self.dev_ID, int(value/self.freq_step))

	@Feat(units="Hz", limits=(min_freq,max_freq,freq_step))
	def frequency_stop(self):
		return self.lib.GetEndFrequency(self.dev_ID)*self.freq_step

	@frequency_stop.setter
	def frequency_stop(self, value):
		self.lib.SetEndFrequency(self.dev_ID, int(value/self.freq_step))

	#Frequency step and dwell time
	@Feat(units="Hz", limits=(freq_step,max_freq-min_freq,freq_step))
	def frequency_step(self):
		return self.lib.GetFrequencyStep(self.dev_ID)*self.freq_step

	@frequency_step.setter
	def frequency_step(self, value):
		self.lib.SetFrequencyStep(self.dev_ID, int(value/self.freq_step))

	@Feat(units="s", limits=(0.01,24*3600)) #From minimum to 24 hours
	def frequency_dwell(self):
		return self.lib.GetDwellTime(self.dev_ID)*1000

	@frequency_dwell.setter
	def frequency_dwell(self, value):
		self.lib.SetDwellTime(self.dev_ID, int(value*1000))	

	#Sweep modes and direction
	@Feat(values={"Repeating Sweep": 1, "Single Sweep": 0})
	def sweep_mode(self):
		return self.sweep_mode_marker

	@sweep_mode.setter
	def sweep_mode(self, value):
		self.lib.SetSweepMode(self.dev_ID,value)
		self.sweep_mode_marker = value

	@Feat(values={"Sweep Up": 1, "Sweep Down": 0})
	def sweep_direction(self):
		return self.sweep_direction_marker

	@sweep_direction.setter
	def sweep_direction(self, value):
		self.lib.SetSweepDirection(self.dev_ID,value)
		self.sweep_direction_marker = value

	#Power level setter
	@Feat(limits=(min_power,max_power,power_step)) #In dBm
	def power_level(self):
		#Frequency readout with limits defined by the hardware.
		return self.lib.GetPowerLevelAbs(self.dev_ID)*self.power_step

	@power_level.setter
	def power_level(self, value):
		self.lib.SetPowerLevel(self.dev_ID, int(value/self.power_step))

	#RF on setter
	@Feat(values={"RF Output ON": True, "RF Output OFF": False, "Device not ready": -2147287040})
	def RF_output(self):
		return self.lib.GetRF_On(self.dev_ID)

	@RF_output.setter
	def RF_output(self, value):
		self.lib.SetRFOn(self.dev_ID,value)

	#Actions-------------------------------------------------------
	@Action()
	def sweep_start(self,value):
		self.lib.StartSweep(self.dev_ID, value) #TODO: Fix this, check if value is necessary here

	@Action()
	def settings_save(self,value):
		self.lib.SaveSettings(self.dev_ID)

#TESTING PURPOSES
if __name__ == "__main__":
	"""Solo call testing function """
	print("Testing the Vaunix LSG402 signal generator driver... \n")
	with SignalGeneratorLSG402(library_name=r"vnx_fsynth.dll", library_folder=r"Users\CMEgroup\Desktop\ArgonneDrivers\vaunix\Libraries") as SigGen:
		print("Asking the LSG402 to self-identify:\n")
		print(SigGen.idn)
		print("\nIdentification successfull. Driver seems to work.\n")