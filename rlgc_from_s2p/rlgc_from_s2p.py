# extract rlgc line model, single frequency exctration
# data reader and extraction based on scikit-rf functionality

import skrf as rf
import numpy as np
import argparse
from matplotlib import pyplot as plt
import os


# create a log that we can dump to terminal and log file
log = []
def append_log (txt):
    log.append(txt + '\n')


# evaluate commandline
parser = argparse.ArgumentParser()
parser.add_argument("s2p",  help="S2P input filename (Touchstone format)")
parser.add_argument("f_ghz", help="extraction frequency in GHz", type=float)
parser.add_argument("l_um", help="line length in micron", type=float)
args = parser.parse_args()


# input data, must be 2-port S2P data
sub = rf.Network(args.s2p)

# physical length must be supplied by user
length = args.l_um*1e-6

# target frequency for pi model extraction
f_target = args.f_ghz*1e9


# frequency class, see https://github.com/scikit-rf/scikit-rf/blob/master/skrf/frequency.py
freq = sub.frequency
append_log('Extract RLGC transmission line model from S2P S-parameter file')
append_log(f'S2P frequency range is {freq.start/1e9} to {freq.stop/1e9} GHz')
append_log(f'Extraction frequency {args.f_ghz}  GHz')
append_log(f'Physical line length: {args.l_um} micron')
assert f_target < freq.stop

# get index for exctraction
f = freq.f
ftarget_index = rf.find_nearest_index(freq.f, f_target)
omega = 2*np.pi*f[ftarget_index]

z11=sub.z[0::,0,0]
y11=sub.y[0::,0,0]
Zline = np.sqrt(z11/y11)

gamma0 = 1/length * np.arctanh(1/(Zline*y11))
gamma_wideband =  gamma0.real + 1j*np.unwrap (gamma0.imag*length, period=np.pi/2)/length


gamma_ftarget = gamma_wideband[ftarget_index]
Zline_ftarget = Zline[ftarget_index]


R = (gamma_ftarget*Zline_ftarget).real
L = (gamma_ftarget*Zline_ftarget).imag / omega
G = (gamma_ftarget/Zline_ftarget).real
C = (gamma_ftarget/Zline_ftarget).imag / omega


append_log('_________________________________________________')
append_log('RLGC line parameters')
append_log(f"Your input for physical line length: {length:.3e} m")
append_log(f"Extraction frequency {f[ftarget_index]/1e9:.3f} GHz")
append_log(f"R   [Ohm/m]: {R:.3e}") 
append_log(f"L'  [H/m]  : {L:.3e}")  
append_log(f"G   [S/m]  : {G:.3e}") 
append_log(f"C'  [H/m]  : {C:.3e}")  
append_log(f"Zline [Ohm]: {Zline_ftarget.real:.3f}")  
append_log('')

log_text = "".join(log)

# print log to terminal
print(log_text)

# output log message to file also, same basename as *.s2p input file but file extension .txt
log_filename = os.path.splitext(args.s2p)[0] + '.txt'
with open(log_filename, "w", encoding="utf-8") as file:
    file.write(log_text)


