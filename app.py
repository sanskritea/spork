import logging
from pathlib import Path

from nspyre import nspyre_init_logger

import nspyre.gui.widgets.save
import nspyre.gui.widgets.load
#import nspyre.gui.widgets.flex_line_plot_widget
from nspyre import MainWidget
from nspyre import MainWidgetItem
from nspyre import nspyreApp
#import nspyre.gui.widgets.snake

import guis.guiElements_CountVsTime as guiElements_CountVsTime
import guis.guiElements_cwODMR as guiElements_cwODMR
import guis.guiElements_Rabi as guiElements_Rabi
import guis.guiElements_PulsedODMR as guiElements_PulsedODMR
import guis.guiElements_BayesianT1 as guiElements_BayesianT1

import guis.guiElements_NewportLineScanning as guiElements_NewportLineScanning
import guis.guiElements_NewportXYScanning as guiElements_NewportXYScanning
import guis.guiElements_NewportXZScanning as guiElements_NewportXZScanning
import guis.guiElements_ThorlabsXYScanning as guiElements_ThorlabsXYScanning

import guis.guiElements_AUXOutVsTime as guiElements_AUXOutVsTime
import guis.guiElements_AttocubeApproach as guiElements_AttocubeApproach
import guis.guiElements_OpticalAttocubeApproach as guiElements_OpticalAttocubeApproach
import guis.guiElements_AttocubeXYScanning as guiElements_AttocubeXYScanning
import guis.guiElements_ThresholdVsTime as guiElements_ThresholdVsTime


# import guis.guiElements_TestPS as guiElements_TestPS
# import guis.guiElements_AttocubeApproachXYScanning as guiElements_AttocubeApproachXYScanning
# import guis.guiElements_Ramsey as guiElements_Ramsey
# import guis.guiElements_OpticalT1 as guiElements_OpticalT1
# import guis.guiElements_MWT1 as guiElements_MWT1
# import guis.guiElements_TwoMeasMWT1 as guiElements_TwoMeasMWT1
# import guis.guiElements_T2star as guiElements_T2star
# import guis.guiElements_Hahn as guiElements_Hahn
# import guis.guiElements_AOMLag as guiElements_AOMLag
# import guis.guiElements_SRSTesting as guiElements_SRSTesting
# import guis.guiElements_ReinitTest as guiElements_ReinitTest


HERE = Path(__file__).parent


if __name__ == '__main__':
    # Log to the console as well as a file inside the logs folder.
    nspyre_init_logger(
        log_level=logging.INFO,
        log_path=HERE / 'logs',
        log_path_level=logging.DEBUG,
        prefix='mango_app',
        file_size=10_000_000,
    )

    # Create Qt application and apply nspyre visual settings.
    app = nspyreApp()

    # Create the GUI.
    main_widget = MainWidget(
        {
            'CountVsTime': MainWidgetItem(guiElements_CountVsTime, 'CustomCountVsTimeWidget', stretch=(1, 1)),

            'ThorlabsXYScanning': MainWidgetItem(guiElements_ThorlabsXYScanning, 'ThorlabsXYScanning_Widget', stretch=(1, 1)),

            'NewportLineScanning': MainWidgetItem(guiElements_NewportLineScanning, 'NewportLineScanning_Widget', stretch=(1, 1)),

            'NewportXYScanning': MainWidgetItem(guiElements_NewportXYScanning, 'NewportXYScanning_Widget', stretch=(1, 1)),

            'NewportXZScanning': MainWidgetItem(guiElements_NewportXZScanning, 'NewportXZScanning_Widget', stretch=(1, 1)),

            'CW ODMR': MainWidgetItem(guiElements_cwODMR, 'CW_ODMR_Widget', stretch=(1, 1)),

            'Pulsed ODMR': MainWidgetItem(guiElements_PulsedODMR, 'Pulsed_ODMR_Widget', stretch=(1, 1)),

            'Rabi': MainWidgetItem(guiElements_Rabi, 'Rabi_Widget', stretch=(1, 1)),

            # 'SRS Testing': MainWidgetItem(guiElements_SRSTesting, 'SRS_Testing_Widget', stretch=(1, 1)),

            # 'Reinit Test': MainWidgetItem(guiElements_ReinitTest, 'Reinit_Test_Widget', stretch=(1, 1)),

            # 'Ramsey': MainWidgetItem(guiElements_Ramsey, 'Ramsey_Widget', stretch=(1, 1)),

            'AttocubeApproach': MainWidgetItem(guiElements_AttocubeApproach, 'Attocube_Approach_Widget', stretch=(1, 1)),

            'OpticalAttocubeApproach': MainWidgetItem(guiElements_OpticalAttocubeApproach, 'Optical_Attocube_Approach_Widget', stretch=(1, 1)),

            'AUXOutVsTime': MainWidgetItem(guiElements_AUXOutVsTime, 'AUXOutVsTime_Widget', stretch=(1, 1)),

            'AttocubeXYScanning': MainWidgetItem(guiElements_AttocubeXYScanning, 'AttocubeXYScanning_Widget', stretch=(1, 1)),

            # 'AttocubeApproachXYScanning': MainWidgetItem(guiElements_AttocubeApproachXYScanning, 'AttocubeApproachXYScanning_Widget', stretch=(1, 1)),

            'ThresholdVsTime': MainWidgetItem(guiElements_ThresholdVsTime, 'ThresholdVsTime_Widget', stretch=(1, 1)),

            # 'OpticalT1': MainWidgetItem(guiElements_OpticalT1, 'OpticalT1_Widget', stretch=(1, 1)),

            # 'MWT1': MainWidgetItem(guiElements_MWT1, 'MWT1_Widget', stretch=(1, 1)),

            # 'TwoMeasMWT1': MainWidgetItem(guiElements_TwoMeasMWT1, 'TwoMeasMWT1_Widget', stretch=(1, 1)),

            'BayesianT1': MainWidgetItem(guiElements_BayesianT1, 'BayesianT1_Widget', stretch=(1, 1)),

            # 'T2star': MainWidgetItem(guiElements_T2star, 'T2star_Widget', stretch=(1, 1)),

            # 'Hahn': MainWidgetItem(guiElements_Hahn, 'Hahn_Widget', stretch=(1, 1)),

            # 'AOMLag': MainWidgetItem(guiElements_AOMLag, 'AOMLag_Widget', stretch=(1, 1)),

            'Plots': {
                'CountVsTime': MainWidgetItem(
                    guiElements_CountVsTime,
                    'CountVsTimePlotWidget',
                    stretch=(100, 100),
                ),

                'Sample XY Scanning': MainWidgetItem(
                    guiElements_ThorlabsXYScanning,
                    'ThorlabsXYScanningPlotWidget',
                    stretch=(100, 100),
                ),

                'Objective Line Scanning': MainWidgetItem(
                    guiElements_NewportLineScanning,
                    'NewportLineScanningPlotWidget',
                    stretch=(100, 100),
                ),

                'Objective XY Scanning': MainWidgetItem(
                    guiElements_NewportXYScanning,
                    'NewportXYScanningPlotWidget',
                    stretch=(100, 100),
                ),

                'Objective XZ Scanning': MainWidgetItem(
                    guiElements_NewportXZScanning,
                    'NewportXZScanningPlotWidget',
                    stretch=(100, 100),
                ),

                'CW_ODMR': MainWidgetItem(
                    guiElements_cwODMR,
                    'cwODMRplotWidget',
                    stretch=(100, 100),
                ),

                'Pulsed_ODMR': MainWidgetItem(
                    guiElements_PulsedODMR,
                    'PulsedODMRplotWidget',
                    stretch=(100, 100),
                ),

                'Rabi': MainWidgetItem(
                    guiElements_Rabi,
                    'RabiPlotWidget',
                    stretch=(100, 100),
                ),

                # 'Ramsey': MainWidgetItem(
                #     guiElements_Ramsey,
                #     'RamseyPlotWidget',
                #     stretch=(100, 100),
                # ),

                # 'AOMLag': MainWidgetItem(
                #     guiElements_AOMLag,
                #     'AOMLagPlotWidget',
                #     stretch=(100, 100),
                # ),

                # 'ReinitTest': MainWidgetItem(
                #     guiElements_ReinitTest,
                #     'ReinitTestPlotWidget',
                #     stretch=(100, 100),
                # ),

                'Attocube Approach': MainWidgetItem(
                    guiElements_AttocubeApproach,
                    'AttocubeApproachPlotWidget',
                    stretch=(100, 100),
                ),

                'Optical Attocube Approach': MainWidgetItem(
                    guiElements_OpticalAttocubeApproach,
                    'OpticalAttocubeApproachPlotWidget',
                    stretch=(100, 100),
                ),

                'AUX Out VS Time': MainWidgetItem(
                    guiElements_AUXOutVsTime,
                    'AUXOutVsTimePlotWidget',
                    stretch=(100, 100),
                ),

                'Attocube XY Scanning': MainWidgetItem(
                    guiElements_AttocubeXYScanning,
                    'AttocubeXYScanningPlotWidget',
                    stretch=(100, 100),
                ),

                # 'Attocube Aproach XY Scanning': MainWidgetItem(
                #     guiElements_AttocubeApproachXYScanning,
                #     'AttocubeApproachXYScanningPlotWidget',
                #     stretch=(100, 100),
                # ),

                'PID Out VS Time': MainWidgetItem(
                    guiElements_ThresholdVsTime,
                    'ThresholdVsTimePlotWidget',
                    stretch=(100, 100),
                ),

                # 'Optical T1': MainWidgetItem(
                #     guiElements_OpticalT1,
                #     'OpticalT1PlotWidget',
                #     stretch=(100, 100),
                # ),

                # 'MW T1': MainWidgetItem(
                #     guiElements_MWT1,
                #     'MWT1PlotWidget',
                #     stretch=(100, 100),
                # ),

                # 'Two Meas MW T1': MainWidgetItem(
                #     guiElements_TwoMeasMWT1,
                #     'TwoMeasMWT1PlotWidget',
                #     stretch=(100, 100),
                # ),

                'Bayesian T1': MainWidgetItem(
                    guiElements_BayesianT1,
                    'BayesianT1PlotWidget',
                    stretch=(100, 100),
                ),

                # 'T2star': MainWidgetItem(
                #     guiElements_T2star,
                #     'T2starPlotWidget',
                #     stretch=(100, 100),
                # ),

                # 'Hahn': MainWidgetItem(
                #     guiElements_Hahn,
                #     'HahnPlotWidget',
                #     stretch=(100, 100),
                # ),
                
            },
            'Save': MainWidgetItem(nspyre.gui.widgets.save, 'SaveWidget', stretch=(1, 1)),
            'Load': MainWidgetItem(nspyre.gui.widgets.load, 'LoadWidget', stretch=(1, 1)),
            #'FlexPlotting': MainWidgetItem(guiElements_taskVsTime, 'FlexLinePlotWidgetWithSomeDefaults', stretch=(1, 1)),
            #'Snake': MainWidgetItem(nspyre.gui.widgets.snake, 'sssss'),
        }
    )
    main_widget.show()
    # Run the GUI event loop.
    app.exec()
        