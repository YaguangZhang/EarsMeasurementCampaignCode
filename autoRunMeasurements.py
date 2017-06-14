# AUTORUNMEASUREMENTS.PY
#
# Load the commands for controlling servos for the antenna, move the antenna
# accordingly, and carry out a signal measurement whenever there is a "Q".
#
# Yaguang Zhang, Purdue University, 2017-06-13

import serial, time, os, sys
from lib import safelyOpenSerial, stripNewlines

'''
Custom variables.
'''
BAUD = 9600
PORT_X = 'COM5'
CMDS_FILE_NAME = 'MIMO_Measurements_Code.txt'

'''
The Script.
'''
print('Summary for autoRunMeasuremnts.py')
print('    BAUD: '+str(BAUD))
print('    PORT_X: '+PORT_X)
#TODO: PORT_Z
print('')

print('Summary for USRP')
sys.stdout.flush()
# Make sure GNURadio is available.
os.system('uhd_find_devices')
print('')

# Read in the commands.
print('Loading commands from .txt file: ' + CMDS_FILE_NAME +'...')
with open(CMDS_FILE_NAME) as f:
    cmds = f.readlines()
# Remove \n and \r (if there is any)
cmds = stripNewlines.stripStrArray(cmds)
print('Commands loaded!')
print('')

print('Opening RS-232 ports...')
# Open the RS-232 ports.
serX = safelyOpenSerial.openPort(PORT_X, BAUD, 1)
# TODO: serZ
print('All ports opened!')
print('\r\n')

print('Starting auto-measurements...')
# Send the commands for the servos and run a signal measurement whenever there
# is a "Q". Note that a servo command will start with the name for the axis
# (i.e. X and Z) to which that command should be sent.
for idx, cmd in enumerate(cmds):
    if cmd[0]=='X':
        serX.write((cmd[1:]+'\r\n').encode('ascii'))
        print(' - Sent to X: ' + cmd[1:])
        # For debugging.
        # time.sleep(1) # In second.
    elif cmd=='Q':
        # Make sure motor for axis X is done before the measurement.
        print('------------------------------------')
        print('Checking whether X is done moving...')
        serX.write('SSDONE\r\n'.encode('ascii'))
        try:
            respX = stripNewlines.stripStr(serX.readline())
            print('    Initial read trial received: "' + respX + '"' )

            numReadsX = 1;
            while(respX != 'DONE'):
                respX = stripNewlines.stripStr(serX.readline())
                numReadsX += 1
                print('    #'+str(numReadsX)+' read trial received: "' \
                    + respX + '"' )

            print('X done moving!')
            print('------------------------------------')
        except Exception as e:
            print('Error: Not able to read from serX!')
            print(e)

            # TODO: serZ

        if stripNewlines.stripStr(respX) =='DONE':
            print('==============================')
            print('Initiating new measurement...')
            print('')
            print('    Time Stamp: ' + str(int(time.time())))
            print('')
            sys.stdout.flush()
            os.system('py -2.7 ./measureSignal.py')
            print('')
            print('Measurement done!')
            print('==============================')
        else:
            print("Servo for axis X is not ready as expected!")
            print("    Expecting DONE but received: " + respX)
# EOF