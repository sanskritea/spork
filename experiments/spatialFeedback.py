"""
This is a basic spatial feedback application

Copyright (c) May 2023, C. Egerstrom
All rights reserved.

This work is licensed under the terms of the 3-Clause BSD license.
For a copy, see <https://opensource.org/licenses/BSD-3-Clause>.
"""

import numpy as np
from scipy.optimize import curve_fit

from nspyre import DataSource
from nspyre import InstrumentGateway

from guis.guiElements_general import flexSave
from rpyc.utils.classic import obtain


def customGaussian(x, A, x0, sigma, K):
        '''Defines a Gaussian function of x with parameters A (vertical stretch factor), x0 (x of peak), sigma (std dev), and K (y-offset)'''
        return K + A*np.exp(-(x-x0)**2/(2*sigma**2))


class SpatialFeedbackMeasurement:

    def guessAndFit(self, lineScanData, steps):
        '''Generates a guess of gaussian paramters based on some lineScanData and its associated steps, then returns the optimal fit parameters
        (with optimal position being arg1)'''
        #Guess height of Gauss is max-min, center is where signal is max, std dev is 500nm, y-offset is min #TODO: Handle std dev smarter
        roughGuessParams = ( np.max(lineScanData)-np.min(lineScanData), steps[np.argmax(lineScanData)], 0.5, np.min(lineScanData) )

        #might want to try-catch this in case optimization fails, but hopefully everything breaks for now. #TODO: Implement this
        optParams, pcov = curve_fit(customGaussian, steps, lineScanData, roughGuessParams )

        return(optParams) #only really care about the position of max (=arg 1)


    def spatialFeedback(self, datasetName:str, guessX, guessY, xSweepDist, ySweepDist, xPoints:int, yPoints:int, collectsPerPt:int=100, iters:int=2):
        """Run spatial feedback (with 20kHz sampling freq) to maximize counts by repeatedly fitting Gaussians 

        Arguments:  *datasetName: name of the dataset to push data to
                    *guessX: starting X position to start fitting (in um) 
                    *guessY: starting Y position to start fitting (in um)
                    *xSweepDist: how far to sweep scan in X (from top to bottom, in um)
                    *ySweepDist: how far to sweep scan in Y (from left to right, in um)
                    *xPoints: number of points to cover the xSweepDist
                    *yPoints: number of points to cover the ySweepDist
                    *collectsPerPt: How many reads to do at 20kHz at each (x,y) point. Default: 100
                    *iters: How many to repeat fitting processing. Default: 2
        """

        # connect to the instrument server
        # connect to the data server and create a data set, or connect to an
        # existing one with the same name if it was created earlier.
        with InstrumentGateway() as gw, DataSource(datasetName) as spatialFeedbackData:

            gw.swabian.ps.constant( ([0],0,0) ) #Make sure the Swab is off (except for the AOM) from any previous measurements
            gw.apdGate.apdOn() #make sure APD is on (which should should kill the LED if it's on)
            
            #need to track sweep positions, sweep data, and fit params as the feedback runs 
            xAllSteps = np.zeros( (0, xPoints) )
            yAllSteps = np.zeros( (0, yPoints) )
            xSweepsData = np.zeros( (0, xPoints) )
            ySweepsData = np.zeros( (0, yPoints) )
            xFits = []
            yFits = []

            for i in range(iters):
                if i==0:
                    latestYguess = guessY
                    latestXguess = guessX

                ySteps =  np.linspace(latestYguess-ySweepDist/2, latestYguess+ySweepDist/2, yPoints, endpoint=True) #not needed for sweeping, but useful for plotting later
                yLineScanData = obtain(gw.fsm.line_scan( {'x':latestXguess, 'y':latestYguess-ySweepDist/2}, {'x':latestXguess, 'y':latestYguess+ySweepDist/2},
                                                yPoints, collectsPerPt) )
                
                #push details of this sweep to tracking, then push those to dataserv
                yAllSteps = np.append(yAllSteps, ySteps.reshape( (1,yPoints) ), axis=0)
                ySweepsData = np.append(ySweepsData, yLineScanData.reshape( (1,yPoints) ), axis=0)
                yFits.append(self.guessAndFit(yLineScanData, ySteps))
                spatialFeedbackData.push({'params': {'sweepRanges': (xSweepDist, ySweepDist), 'scanPointsPerAxis': (xPoints, yPoints), 'collectsPerPt': collectsPerPt, 'NumIters': iters},
                                'title': 'SpatialFeedback',
                                'datasets': {'xSteps':xAllSteps, 'ySteps':yAllSteps, 'xSweepData': xSweepsData, 'ySweepData': ySweepsData,
                                             'xFits': xFits, 'yFits': yFits}
                })
                #then update y guess
                latestYguess = yFits[i][1]

                #And now for X
                xSteps =  np.linspace(latestXguess-xSweepDist/2, latestXguess+xSweepDist/2, xPoints, endpoint=True) #not needed for sweeping, but useful for plotting later
                xLineScanData = obtain( gw.fsm.line_scan( {'x':latestXguess-xSweepDist/2, 'y':latestYguess}, {'x':latestXguess+xSweepDist/2, 'y':latestYguess},
                                                xPoints, collectsPerPt) )
                
                #push details of this sweep to tracking, then push those to dataserv
                xAllSteps = np.append(xAllSteps, xSteps.reshape( (1,xPoints) ), axis=0)
                xSweepsData = np.append(xSweepsData, xLineScanData.reshape( (1,xPoints) ), axis=0)
                xFits.append(self.guessAndFit(xLineScanData, xSteps))
                spatialFeedbackData.push({'params': {'sweepRanges': (xSweepDist, ySweepDist), 'scanPointsPerAxis': (xPoints, yPoints), 'collectsPerPt': collectsPerPt, 'NumIters': iters},
                                'title': 'SpatialFeedback',
                                'datasets': {'xSteps':xAllSteps, 'ySteps':yAllSteps, 'xSweepData': xSweepsData, 'ySweepData': ySweepsData,
                                             'xFits': xFits, 'yFits': yFits}
                })
                #now update x guess
                latestXguess = xFits[i][1] 


            flexSave(datasetName, 'spatialFeedback', 'final') #after measurement finishes
            gw.apdGate.apdOff()
            gw.fsm.move( (latestXguess, latestYguess) )
            print(f'Moved to best guess ({latestXguess},{latestYguess}) after fitting {iters} times')



if __name__ == '__main__':
    pass

