"""
    lantz.drivers.thorlabs.DCC1240C
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    :copyright: 
    N.Delegan - Dec 2020.
    Use as you wish, be nice.

    :description:
    Driver for the DCC1240C Thorlabs camera.
    
    :versions/changelog: 
    * V0 - Dec 2020 -  Driver started

    :dependencies: 
        	
	:usage:
"""

# Import stuff
from lantz import Action, Feat, Driver
import sys, time, clr

clr.AddReference('System') # For some C# (.NET) native stuff via clr import
from System import Int32, Byte, Array

import numpy as np

class ThorlabsDCC1240CCamera(Driver):
    '''Driver class for the DCC1240C Thorlabs camera.
    '''
    def __init__(self, serialNo = '4103124864'):
        # serialNo: the device serialNo, default is the one we use, can be redefined.

        # Import the library location
        try:
            # Define path of Scientific Imaging drivers on the system
            thorcam_path = r'C:\Program Files\Thorlabs\Scientific Imaging\DCx Camera Support\Develop\DotNet\signed'
            # Add Thorlabs Scientific Imaging Program directory to PATH
            sys.path.append(thorcam_path)

            # Add reference to library (using clr)
            clr.AddReference('uc480DotNet')
            
        except:
            raise RuntimeError(f'Unable to load library, check location and reference.')

        # Import namespaces
        try:
            import uc480
            import uc480.Info as uc480Info
            import uc480.Defines as uc480Defines
            import uc480.Configuration as uc480Configs          
        except:
            raise RuntimeError(f'Cannot import the uc480 and other related namespaces.')
        
        # Create camera object
        cameraID = 0 # First camera found will be initialized
        self.cam = uc480.Camera(cameraID) # Internally calls the .Init() method unless no cameraID is specified.

        # # Initialize the camera (Only needed when no parameter is fed to .Camera())
        # call_result = self.cam.Init()
        # if call_result != uc480Defines.Status.Success:
        #     raise RuntimeError(f'Unable to initiate camera.')

        # Get display mode
        dummy_var = uc480Defines.DisplayMode.DiB
        call_result, DisplayMode = self.cam.Display.Mode.Get(dummy_var)

        # Dictionary for convenience
        DisplayModeDictionary = {
            1:'DiB',
            4:'Direct3D',
            8:'OpenGL',
            0x00000800:'Mono',
            0x00001000:'Bayer'}
        if call_result != uc480Defines.Status.Success:
            raise RuntimeError(f'Unable to get display mode, camera likely not able to connect.')
        else:
            # print(f'Camera is in display mode: {DisplayModeDictionary[DisplayMode]}')
            None

        # Setting the display mode
        call_result = self.cam.Display.Mode.Set(uc480Defines.DisplayMode.DiB)
        if call_result != uc480Defines.Status.Success:
            raise RuntimeError(f'Unable to set the desired display mode with error:{call_result}. See the uc480.Defines.Status definitions for more info.')

        # Allocate memory
        dummy_s32MemID = Int32(10) # Dummy int value
        call_result, s32MemID = self.cam.Memory.Allocate(dummy_s32MemID, True) # True to set as active memory
        if call_result != uc480Defines.Status.Success:
            raise RuntimeError(f'Unable to allocate memory for camera.')

        # Acquire image
        call_result = self.cam.Acquisition.Freeze(uc480Defines.DeviceParameter.Wait) # Wait for image to be acquired before returning.

        dummy_Actv_s32MemID = Int32(10) # Dummy int value
        call_result, Actv_s32MemID = self.cam.Memory.GetActive(dummy_Actv_s32MemID)
        if call_result != uc480Defines.Status.Success:
            raise RuntimeError(f'Unable to get active memory ID, not sure what can cause this.')

        if Actv_s32MemID != s32MemID:
            raise RuntimeError(f'Allocated memory and active memory do not seem to match. Really not sure how this happened.')

        # Lock memory temporarily
        call_result = self.cam.Memory.Lock(Actv_s32MemID) # Note sure if this is necessary, will debug later.
        if call_result != uc480Defines.Status.Success:
            raise RuntimeError(f'Unable to lock memory after acquisition.')

        # dummy_ImgBitmap = Drawing.Bitmap(2,2) # Create a dummy bitmap
        # call_result, ImgBitmap = self.cam.Memory.CopyToBitmap(Actv_s32MemID, dummy_ImgBitmap) # Should output a 1280x1024 WxH bitmap
        # # result, bitmap = cam.Memory.ToBitmap(mem_id, dummy_var) #What is the difference with copytobitmap?
        # if call_result != uc480Defines.Status.Success:
        #     raise RuntimeError(f'Unable to copy memory image to bitmap.')

        # Get the image (Copy to Array)
        dummy_ImgArray = Array[int]([0])
        # dummy_ImgArray = bytearray(5)
        call_result, ImgArray = self.cam.Memory.CopyToArray(Actv_s32MemID, uc480Defines.ColorMode.RGB8Packed, dummy_ImgArray) #
        if call_result != uc480Defines.Status.Success:
            raise RuntimeError(f'Unable to copy memory image to array.')

        # Unlock the memory that was locked
        call_result = self.cam.Memory.Unlock(Actv_s32MemID) # Note sure if this is necessary, will debug later.
        if call_result != uc480Defines.Status.Success:
            raise RuntimeError(f'Unable to unlock memory after copying image.')

        # Create np array
        npImgArray = np.fromiter(ImgArray, int, ImgArray.LongLength) # 3 by 1280 by 1024
        ImgWidth = 1280
        ImgHeight = 1024
        # Reshape into image matrix
        npImgArray_shaped = npImgArray.reshape((ImgHeight, ImgWidth, 3)) # Encoding is RGB8 packed 0...255 each unsigned hex

        # output_arr = numpy.take_along_axis(data_arr, index_arr[:, :, numpy.newaxis], axis=2)
        # output_arr = output_arr[:,:,0]  # Since take_along_axis keeps the same number of dimensions
        # demos_arr = np.take_along_axis(img_arr, bayer_arr.reshape([img_height, img_width, 1]), axis=2).reshape([img_height, img_width])

        # FOR DEBUGGING: now display the image from the raw NumPy array:
        from matplotlib import pyplot as PLT
        print(npImgArray_shaped.shape)
        PLT.imshow(npImgArray_shaped)
        PLT.show()

        # # Set the display mode to bitmap
        # self.cam.Display.Mode.Set(uc480Defines.DisplayMode.DiB)

        # # Set the color mode
        # self.cam.PixelFormat.Set(uc480Defines.ColorMode.RGBA8Packed)

        # # Trigger mode
        # self.cam.Trigger.Set(uc480Defines.TriggerMode.Software)

        # # some_int_generic = Dictionary[str, int]()

        # self.cam.Acquisition.Freeze(uc480Defines.DeviceParameter.Wait)
        # # print(self.cam.Memory.ImageBuffer.__dir__())
        # # self.cam.Memory.ImageBuffer.SetEnable(True)
        # # self.cam.Memory.ImageBuffer.TransferImage(1, -1)
        # # print(self.cam.Memory.ImageBuffer.GetEnable)
        call_result = self.cam.Exit()
        if call_result != uc480Defines.Status.Success:
            raise RuntimeError(f'Unable to exit camera object.')
        return None



    def finalize(self):
        None

    #TODO: Crete rest of the driver once the __init__ unit testing works


# For driver testing
if __name__ == "__main__":
    with ThorlabsDCC1240CCamera(serialNo = '4103124864') as thorcam:
        print("Testing the camera driver...\n")