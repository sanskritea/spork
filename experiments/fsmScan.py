"""
This is a basic Fsm Scan Application

Copyright (c) May 2023, C. Egerstrom
All rights reserved.

This work is licensed under the terms of the 3-Clause BSD license.
For a copy, see <https://opensource.org/licenses/BSD-3-Clause>.
"""

import numpy as np

from nspyre import DataSource
from nspyre import InstrumentGateway

from guis.guiElements_general import flexSave
from rpyc.utils.classic import obtain


class FSMScanMeasurement:

    def fsmScan(self, datasetName:str, centerX, centerY, xSweepDist, ySweepDist, xPoints:int, yPoints:int, collectsPerPt:int=100):
        """Run an FSM scan (@20kHz sampling freq)

        Arguments:  *datasetName: name of the dataset to push data to
                    *centerX: X position of center of scan (in um) 
                    *centerY: Y position of center of scan (in um)
                    *xSweepDist: how far to sweep scan in X (from top to bottom, in um)
                    *ySweepDist: how far to sweep scan in Y (from left to right, in um)
                    *xPoints: number of points to cover the xSweepDist
                    *yPoints: number of points to cover the ySweepDist
                    *collectsPerPt: How many reads to do at 20kHz at each (x,y) point. Default: 100
        """

        # connect to the instrument server
        # connect to the data server and create a data set, or connect to an
        # existing one with the same name if it was created earlier.
        with InstrumentGateway() as gw, DataSource(datasetName) as fsmScanData:

            gw.swabian.ps.constant( ([0],0,0) ) #Make sure the Swab is off (except for the AOM) from any previous measurements
            gw.apdGate.apdOn() #make sure APD is on (which should should kill the LED if it's on)

            xSteps = np.array([]) #will fill this in as data comes in
            ySteps =  np.linspace(centerY-ySweepDist/2, centerY+ySweepDist/2, yPoints, endpoint=True) #not needed for sweeping, but useful for plotting later
            fsmScanCounts = np.zeros( (0,yPoints) )

            for xStep in np.linspace(centerX-xSweepDist/2, centerX+xSweepDist/2, xPoints, endpoint=True):
                lineScanData = obtain(gw.fsm.line_scan( {'x':xStep, 'y':centerY-ySweepDist/2}, {'x':xStep, 'y':centerY+ySweepDist/2},
                                                 yPoints, collectsPerPt))
                fsmScanCounts = np.append(fsmScanCounts, lineScanData.reshape( (1,yPoints) ), axis=0)
                xSteps = np.append(xSteps, xStep)

                fsmScanData.push({'params': {'CenterOfScan': (centerX, centerY), 'sweepRanges': (xSweepDist, ySweepDist), 'scanPointsPerAxis': (xPoints, yPoints), 'collectsPerPt': collectsPerPt},
                                'title': 'FsmScan',
                                'xLabel': 'X Position',
                                'yLabel': 'Y Position',
                                'zLabel': 'Counts (avg\'d and normed)',
                                'datasets': {'xSteps':xSteps, 'ySteps':ySteps, 'ScanCounts': fsmScanCounts}
                })


            flexSave(datasetName, 'fsmScan', 'final') #after measurement finishes
            gw.apdGate.apdOff()
            gw.fsm.move( (0,0) ) #move to (0,0) after scanning
            print('Moved to (0,0) after scanning') #I'll forget this happens without the print, but at some point I'll align with an offset FSM if I don't do this



if __name__ == '__main__':
    pass

