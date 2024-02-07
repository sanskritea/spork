
#TO-DO, NEEDS SOME DEBUGGING
from functools import partial
from importlib import reload

from pyqtgraph.Qt import QtWidgets
from pyqtgraph import SpinBox
from nspyre import LinePlotWidget
from nspyre import ParamsWidget
from nspyre import ProcessRunner
from nspyre import DataSink

import experiments.taskVsTime as taskVsTime
from guis.guiElements_general import AutoSaveWidget
from guis.guiElements_general import flexSave
