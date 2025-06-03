import logging
from pathlib import Path

from nspyre import nspyre_init_logger
from nspyre import serve_instrument_server_cli
from nspyre import InstrumentServer


HERE = Path(__file__).parent

# log to the console as well as a file inside the logs folder
nspyre_init_logger(
    logging.INFO,
    log_path=HERE / 'logs',
    log_path_level=logging.DEBUG,
    prefix='mango_inserv',
    file_size=10_000_000,
)


# create a new instrument server
with InstrumentServer() as inserv:

   inserv.add('swabian', HERE / 'drivers' / 'swabian' / 'SwabianPS82.py', 'SwabianPulseStreamer82')

   inserv.add('sg', 'drivers.stanford.sg396', 'SG396', ['TCPIP0::192.168.1.79::inst0::INSTR'], import_or_file='import')

   inserv.add('mfli', HERE / 'drivers' / 'zurich' / 'mfli.py', 'MFLI')

   # inserv.add('kdc', HERE / 'drivers' / 'thorlabs' / 'KDC101.py', 'KDC')

   inserv.add('esp', HERE / 'drivers' / 'newport' / 'ESP302.py', 'NESP')





   # inserv.add('nidaq', HERE / 'drivers' / 'ni' / 'nidaq_final.py', 'NIDAQ')
   # inserv.add('nidaq', HERE / 'drivers' / 'ni' / 'ni_photonCounting.py', 'NIDAQ_PhotonCounter')
   # inserv.add('scanner', HERE / 'drivers' / 'attocube' / 'scanner.py', 'scanner')
   # inserv.add('laser', HERE / 'drivers' / 'cobolt' / 'cobolt.py', 'CoboltLaser')
   # inserv.add('vaunix', HERE / 'drivers' / 'vaunix' / 'LSG402.py', 'SignalGeneratorLSG402')
   # inserv.add('fsm', HERE / 'drivers' / 'newport' / 'FSM_via_nidaqmx.py', 'FSM_2D')
   # inserv.add('whiteLED', HERE / 'drivers' / 'thorlabs' / 'LEDD1B.py', 'WhiteLightSourceLEDD1B')
   # inserv.add('flipper', HERE / 'drivers' / 'thorlabs' / 'MFF102M.py', 'MFF102MFlipper')
   # inserv.add('kim', HERE / 'drivers' / 'thorlabs' / 'KIM101_XY.py', 'KIM_XY_CTLR')
   # inserv.add('standa', HERE / 'drivers' / 'standa' / 'Standa_8MT173-30DCE2.py', 'StandaZFocusDriver')
   # inserv.add('schottky', HERE / 'drivers' / 'pasternack' / 'PE8012P.py', 'PE8012P_Schottky')
   # inserv.add('apdGate', HERE / 'drivers' / 'mpd' / 'mpdAPD.py', 'MPD_PDMSeries_APD')
   

   # run a CLI (command-line interface) that allows the user to enter
   # commands to control the server
   serve_instrument_server_cli(inserv)