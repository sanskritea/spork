#TO-DO, NEEDS SOME DEBUGGING
from pyqtgraph.Qt import QtWidgets
from pyqtgraph import SpinBox
from nspyre import ParamsWidget
from nspyre import DataSink

from nspyre.gui.widgets.save import save_json
from os import path, mkdir
from datetime import datetime



class AutoSaveWidget(QtWidgets.QWidget):

    def __init__(self, defaultShouldAutosave: bool, defaultInterval: int):
        super().__init__() #create a widget to get autosave interval. This is definitely not the easiest way, but this should work
        self.paramsWidget = ParamsWidget({
            'autosaveInterval': {
                'display_text': 'Autosave Interval',
                'widget': SpinBox(
                    value=defaultInterval,
                    int=True,
                    dec=True,
                ),
            },
        })

        self.autosaveCheck = QtWidgets.QCheckBox('Autosave?') #create a checkbox to determine if autosave should be default 
        if defaultShouldAutosave:
            self.autosaveCheck.setChecked(True)

        autosaveLayout = QtWidgets.QVBoxLayout() #set them up in a box together
        autosaveLayout.addWidget(self.autosaveCheck)
        autosaveLayout.addWidget(self.paramsWidget)

        self.setLayout(autosaveLayout)


    def shouldAutosave(self):
        '''Accessor method for if autosaving box is checked'''
        return( self.autosaveCheck.isChecked() )
    

    def getAutosaveInterval(self):
        '''Accessor method for autosaving interval'''
        return( self.paramsWidget.autosaveInterval )
    

# def flexSave(datasetName:str, expType:str, saveType:str, dirs:list = ['C:/Users/awschlab/Data/']):
def flexSave(datasetName:str, expType:str, saveType:str, dirs:list = ['C:\\Users\\awschlab\\Desktop\\Data']): #TODO: Make it so this hangs up data acquisition. Will need to use ProcessRunner, multithread acq and saving, and then make them talk to each other nicely. Gross
    '''Creates a save of the data a specified directory(ies) in a Dir\\DATE(YYMMDD)\\EXPERIMENT_TYPE\\EXP_TYPE TIME(HHMMSS) SAVE_TYPE.json structure
    Arguments:  *datasetName:str, name of data to be saved from dataserv
                *expType:str, name of the experiment (or name of folder to save in under the date)
                *saveType:str, typically something like auto, closeout, final
                *dirs:list, dirs ending in \\ to save data to. Default is Jasper's Data driver'''

    if not len(dirs) > 0:
        raise ValueError('No directories specified for custom autosaver')
    
    # testing
    print('dirs : ', dirs)
    now = datetime.now()
    with DataSink(datasetName) as dataSink:
        try:
            dataSink.pop(1)
            for dir in dirs:
                # datePath = datePath = dir + now.strftime('%y%m%d')+ '/'
                datePath = datePath = dir + now.strftime('%y%m%d')+ '\\' #dir for the date
                # expPath = datePath + expType + '/'
                expPath = datePath + expType + '\\' #dir for the exp on date
                # print('expPath : ', expPath)
                filePath = expPath + expType + now.strftime('%H%M%S') + saveType + '.json' #filename for the json, assumes autosaving is at <1Hz to avoid name collisions
                # print('filePath : ', filePath)
                if not path.isdir(datePath): #create the date dir if it doesn't already exist
                    mkdir(datePath)
                if not path.isdir(expPath): #create the exp dir inside of date dir if it doesn't already exist
                    mkdir(expPath)
                save_json(filePath, dataSink.data) #actually save the data as a json there
        except TimeoutError:
            raise ValueError(f'No data with name \'{datasetName}\' to save (or timed out)')
