from PyQt6.QtWidgets import QApplication, QWidget, QHBoxLayout, QVBoxLayout, QPushButton, QComboBox, QToolBox, QLineEdit
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore
from pyqtgraph.Qt import QtGui
from pyqtgraph.Qt import QtWidgets
from pyqtgraph import SpinBox 

import numpy as np
from pyqtgraph import functions as fn
import os 

from nspyre import InstrumentGateway
from nspyre import DataSource
from nspyre import StreamingList
from nspyre import DataSink



class Crosshair(QtCore.QObject):
    sigPositionChanged = QtCore.pyqtSignal(object)
    sigPositionChangeFinished = QtCore.pyqtSignal(object)

    def __init__(self, plot_item, pos, **kwargs):
        super().__init__()
        self.kwargs = kwargs
        self.moving = False
        self.hovering = False

        #pen = fn.mkPen((255, 255, 0, 127))
        pen = fn.mkPen((0, 255, 0, 127))
        self.vLine = pg.InfiniteLine(angle=90, pen=pen, movable=False)
        self.hLine = pg.InfiniteLine(angle=0,  pen=pen, movable=False)
        self.vLine.hoverEvent, self.hLine.hoverEvent = self.hoverEvent, self.hoverEvent
        self.vLine.mouseDragEvent, self.hLine.mouseDragEvent = self.mouseDragEvent, self.mouseDragEvent

        #self.center_dot = pg.ScatterPlotItem(pos=[pos], pen=fn.mkPen((255,0,0, 127)), brush=(255,0,0), symbol='o', size=3)
        self.center_dot = pg.ScatterPlotItem(pos=[pos], pen=fn.mkPen((0,255,255, 127)), brush=(255,0,0), symbol='o', size=3)

        plot_item.addItem(self.vLine, ignoreBounds=True)
        plot_item.addItem(self.hLine, ignoreBounds=True)
        plot_item.addItem(self.center_dot, ignoreBounds=True)
        self.plot_item = plot_item
        self.set_pos(pos)

    def set_pos(self, pos, emit_sig=True):
        if isinstance(pos, QtCore.QPointF):
            self.pos = [pos.x(), pos.y()]
        else:
            self.pos = list(pos)
        self.vLine.setPos(self.pos[0])
        self.hLine.setPos(self.pos[1])
        self.center_dot.setData(pos=[pos])
        if emit_sig: self.sigPositionChanged.emit(self.get_pos())

    def mouseDragEvent(self, ev):
        if ev.button() == QtCore.Qt.MouseButton.LeftButton:
            if ev.isStart():
                self.moving = True
            ev.accept()

            if not self.moving:
                return
            self.set_pos(self.plot_item.vb.mapSceneToView(ev.scenePos()))
            if ev.isFinish():
                self.moving = False
                self.sigPositionChangeFinished.emit(self.get_pos())

    def hoverEvent(self, ev):
        if (not ev.isExit()) and ev.acceptDrags(QtCore.Qt.MouseButton.LeftButton):
            self.hovering = True
            for line in [self.vLine, self.hLine]: line.currentPen = fn.mkPen(255, 0,0)
        else:
            self.hovering = False
            for line in [self.vLine, self.hLine]: line.currentPen = line.pen
        for line in [self.vLine, self.hLine]:
            line.update()

    def get_pos(self):
        return self.pos

    def delete(self):
        self.hLine.deleteLater()
        self.vLine.deleteLater()
        self.center_dot.deleteLater()


class CrosshairAddon(QtWidgets.QWidget):
    sigCrosshairAdded = QtCore.pyqtSignal(object)
    sigCrosshairRemoved = QtCore.pyqtSignal(object)
    def __init__(self, plot_item, **kwargs):
        super().__init__(**kwargs)
        self.plot_item = plot_item
        self.cross_list = list()
        self._spinbox_decimals = 4
        
        self.streaming_cross = [] 

        self.build_ui()


    @property
    def spinbox_decimals(self):
        return self._spinbox_decimals

    @spinbox_decimals.setter
    def spinbox_decimals(self, val):
        if self._spinbox_decimals != val:
            for r in range(self.table.rowCount()):
                self.table.cellWidget(r,0).setDecimals(val)
                self.table.cellWidget(r,1).setDecimals(val)
            self._spinbox_decimals = val

    def build_ui(self):
        main_layout = QVBoxLayout()


        #FSM channel: 
        ao_layout = QHBoxLayout()
        self.dropdownX = QComboBox()
        self.dropdownY = QComboBox() 
        self.dropdownX.addItems(['ao12', 'ao0', 'ao1', 'ao2', 'ao3', 'ao4'])
        self.dropdownY.addItems(['ao13', 'ao0', 'ao1', 'ao2', 'ao3', 'ao4'])

        ao_layout.addWidget(self.dropdownX)
        ao_layout.addWidget(self.dropdownY)

        main_layout.addLayout(ao_layout)


        ###go to 0,0
        zero_layout = QHBoxLayout()
        self.buttonZ = QPushButton('Go to Zero')
        self.buttonZ.clicked.connect(self.on_buttonZ_click)
        zero_layout.addWidget(self.buttonZ)

        main_layout.addLayout(zero_layout)


        ##move xhair
        top_layout = QHBoxLayout()
        self.buttonMove = QPushButton('Move to Crosshair #')
        self.buttonMove.clicked.connect(self.on_buttonMove_click)
        top_layout.addWidget(self.buttonMove)

        self.dropdown = QComboBox() 
        top_layout.addWidget(self.dropdown)

        main_layout.addLayout(top_layout)


        ## input crosshairs
        input_layout = QHBoxLayout()
        self.buttonIn = QPushButton('Import Crosshairs: ')
        self.buttonIn.clicked.connect(self.on_buttonImport_click)
        input_layout.addWidget(self.buttonIn)

        self.inputText = QLineEdit(self)
        input_layout.addWidget(self.inputText)

        main_layout.addLayout(input_layout)


        ## Export crosshairs: 

        ex_layout = QHBoxLayout()
        self.buttonOut = QPushButton('Export Crosshairs: ')
        self.buttonOut.clicked.connect(self.on_buttonExport_click)
        ex_layout.addWidget(self.buttonOut)

        self.exText = QLineEdit(self)
        ex_layout.addWidget(self.exText)

        main_layout.addLayout(ex_layout)


        ###DATASERV BUTTONS: 

        ds_layout = QHBoxLayout() 

        self.buttonPush = QPushButton('Push to Dataserv')
        self.buttonPush.clicked.connect(self.pushdataserv)
        ds_layout.addWidget(self.buttonPush) 
        self.pushText = QLineEdit(self)
        self.pushText.setText("xhair0")
        ds_layout.addWidget(self.pushText)

        self.buttonPull = QPushButton('Pull from Dataserv')
        self.buttonPull.clicked.connect(self.pulldataserv)
        ds_layout.addWidget(self.buttonPull) 
        self.pullText = QLineEdit(self)
        self.pullText.setText("xhair0")
        ds_layout.addWidget(self.pullText)

        main_layout.addLayout(ds_layout)



        ####Crosshair table
        self.table = QtWidgets.QTableWidget()
        self.table.setColumnCount(3)
        self.table.setHorizontalHeaderLabels(['x','y','Delete'])

        #Add control (for now just an add button)
        add_btn = QtWidgets.QPushButton('+ Add')
        def add():
            r = np.array(self.plot_item.getViewBox().viewRange())
            self.add_crosshair(r.mean(axis=1))
        add_btn.clicked.connect(lambda: add())

        #Add a decimal precision box
        decimal_input = SpinBox(value=self.spinbox_decimals, minStep=1, dec=False, int=True, bounds=(0, None), step=1)
        
        decimal_input.valueChanged.connect(lambda x: setattr(self, 'spinbox_decimals', decimal_input.value()))

        ctrl_layout = QtWidgets.QFormLayout()
        ctrl_layout.addRow('Add Crosshair', add_btn)
        ctrl_layout.addRow('Floating point precision', decimal_input)
        ctrl_widget = QtWidgets.QWidget()
        ctrl_widget.setLayout(ctrl_layout)


        #layout = QtWidgets.QGridLayout()
        main_layout.addWidget(ctrl_widget)
        main_layout.addWidget(self.table)
        main_layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(main_layout)


    def pulldataserv(self): 
        
        with DataSink(self.pullText.text()) as sink :
            #print(sink)
           # print(sink.__dir__())
            try:
                sink.pop(timeout=3)
                
                for currentcross in reversed(self.cross_list): 

                    self.remove_crosshair(currentcross)

                if len(sink.datasets['crosslist']) > 0:
                    for cross in sink.datasets['crosslist']:
                        self.add_crosshair(np.array(cross))
                    print("sinked from: " + self.pullText.text())
            except TimeoutError:
                print("Sink does not exist, or timed out after 3 sec.")

        
         
    def pushdataserv(self): 
        cross_pos = []
        for cross in self.cross_list: 
            cross_pos.append(cross.get_pos())
        with DataSource(self.pushText.text()) as crossdata: 
            crossL = StreamingList() 
            crossdata.push({'datasets': {'crosslist': cross_pos}})
        print("pushed to: " + self.pushText.text())
        
    def on_buttonZ_click(self): 
        print("zero")
    def on_buttonMove_click(self): 
        selected_row = self.dropdown.currentIndex() 
        if selected_row >= 0: 
            cross = self.cross_list[selected_row]
            cords=cross.get_pos()
            print(f'Moving: ({cords[0]}, {cords[1]})')
            #print(len(cords))
    def on_buttonImport_click(self): 
        
        if os.path.exists(self.inputText.text()):
            data = np.genfromtxt(self.inputText.text(), delimiter=',', skip_header=1)
            for cross in data: 
                self.add_crosshair(np.array(cross))


            print("Imported from: " + self.inputText.text())

        if not os.path.exists(self.inputText.text()):
            print("ERROR: FILE DOES NOT EXIST")


    def on_buttonExport_click(self):
        cross_cords = []
        for cross in self.cross_list:
            cross_cords.append(cross.get_pos())

        data_array = np.array(cross_cords)

        file_path = self.exText.text() + '.csv'


        if not os.path.exists(file_path):
            with open(file_path, 'w') as file:
                pass 


        np.savetxt(file_path, data_array, delimiter=',', header='X,Y', comments='')
        print("Exported to: " + self.exText.text()+ ".csv")

    def update_dropdown(self): 
        self.dropdown.clear()
        for crossNumb in range(len(self.cross_list)):
            self.dropdown.addItem(str(crossNumb +1 ))

    def add_crosshair(self, pos, **kwargs):
        # Add the Crosshair to the list
        cross = Crosshair(self.plot_item, pos, **kwargs)
        self.cross_list.append(cross)

        # Add the table entry
        row = len(self.cross_list)-1
        self.table.insertRow(row)

        # Add the x,y widgets
        def update_pos(axis, value):
            if cross.moving:
                return
            cur = cross.get_pos()
            if axis==0:
                cross.set_pos([value, cur[1]], emit_sig=False)
            elif axis==1:
                cross.set_pos([cur[0], value], emit_sig=False)
        def lambda_gen(axis):
            return lambda obj: update_pos(axis, obj.value())

        for i in range(2):
            w = SpinBox(value = pos[i], dec=True, decimals=self.spinbox_decimals)

            self.table.setCellWidget(row, i, w)
            w.sigValueChanged.connect(lambda_gen(i))
            w.setMinimumHeight(30)
        self.update_dropdown()

        # Add a remove button
        btn = QtWidgets.QPushButton('X')
        btn.clicked.connect(lambda: self.remove_crosshair(cross))
        self.table.setCellWidget(row, 2, btn)

        # Link the position of the cross to the numbers in the table
        cross.sigPositionChanged.connect(lambda: self.update_table_entry(cross))
        self.sigCrosshairAdded.emit(cross)


    def _find_index(self, cross):
        for i in range(len(self.cross_list)):
            if self.cross_list[i] == cross:
                return i

    def update_table_entry(self, cross):
        row = self._find_index(cross)
        self.table.cellWidget(row, 0).setValue(cross.get_pos()[0])
        self.table.cellWidget(row, 1).setValue(cross.get_pos()[1])

    def remove_crosshair(self, cross):
        index = self._find_index(cross)
        self.table.removeRow(index)
        cross.delete()
        self.cross_list.pop(index)
        self.sigCrosshairRemoved.emit(index)

        self.update_dropdown()


    def __getitem__(self, k):
        return self.cross_list[k].get_pos()

    def __iter__(self):
        for cross in self.cross_list:
            yield cross

    def __len__(self):
        return len(self.cross_list)