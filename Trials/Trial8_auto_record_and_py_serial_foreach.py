# from lib import Test2_baseband_rx_file_sink_Mod
# Doesn't work. Python won't reload the gnu radio script for each execution.

import serial, time, os

strsToSend = ['FL500000']*5
#strsToSend = ['Hello ...', 'My name is Python ...', 'What is your name?', '...']

port = "COM5"
baud = 19200
ser = serial.Serial(port, baud, timeout=1)

for str in strsToSend:
    ser.write(str.encode('ascii'))
    time.sleep(5) # In second.

    # Test2_baseband_rx_file_sink_Mod.main()
    # execfile('./ModifiedPyFiles/Test2_baseband_rx_file_sink_Mod.py')
    os.system('py -2.7 ./ModifiedPyFiles/Test2_baseband_rx_file_sink_Mod.py')