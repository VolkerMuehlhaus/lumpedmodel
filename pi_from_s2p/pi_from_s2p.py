# extract simple pi model, single frequency exctration, not accurate for wideband use
# data reader and extraction based on scikit-rf functionality

import skrf as rf
import math
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
args = parser.parse_args()


# input data, must be 2-port S2P data
sub = rf.Network(args.s2p)

# target frequency for pi model extraction
f_target = args.f_ghz*1e9


# frequency class, see https://github.com/scikit-rf/scikit-rf/blob/master/skrf/frequency.py
append_log('Extract simple inductor pi model from S2P S-parameter file')
append_log(f'S2P frequency range is {sub.frequency.start/1e9} to {sub.frequency.stop/1e9} GHz')
append_log(f'Extraction frequency: {args.f_ghz} GHz')
assert f_target < sub.frequency.stop

# if the input data has DC point, remove that because it will throw warnungs later
if sub.frequency.start == 0:
    # resample to start at 1 GHz (or closest value)
    newrange = '1-' + str(sub.frequency.stop/1e9) + 'ghz'
    sub = sub[newrange]



freq = sub.frequency
f = freq.f

z11=sub.z[0::,0,0]
z21=sub.z[0::,1,0]
z12=sub.z[0::,0,1]
z22=sub.z[0::,1,1]


# 2-port to 1-port conversion
Zdiff = z11-z12-z21+z22
freq = sub.f
omega = sub.f*2*math.pi
Qdiff = Zdiff.imag/Zdiff.real
Ldiff = Zdiff.imag/omega
Rdiff = Zdiff.real

# find frequency of maximum Q factor 
Qmax = max(Qdiff)
Qmax_index = rf.find_nearest_index(Qdiff, Qmax)
f_Qmax = freq[Qmax_index]

# calculate inductor circuit model

# calculate pi model 
# Zser = series element
# Zshunt1 = left shunt element
# Zshunt2 = right shunt element

y11=sub.y[0::,0,0]
y21=sub.y[0::,1,0]
y12=sub.y[0::,0,1]
y22=sub.y[0::,1,1]
ymn = (y12+y21)/2

Zshunt1 =  1 / (y11 + ymn)
Zshunt2 =  1 / (y22 + ymn)
Zseries = -1 / (ymn)

# values over frequency
Rseries = Zseries.real
Lseries = Zseries.imag/omega
Cshunt1 = -1 / (omega*Zshunt1.imag)
Cshunt2 = -1 / (omega*Zshunt2.imag)
Rshunt1 = (1 / (y11+ymn)).real
Rshunt2 = (1 / (y22+ymn)).real
Rshunt = (Rshunt1+Rshunt2)/2

ftarget_index = rf.find_nearest_index(freq, f_target)
Rseries_ftarget = Rseries[ftarget_index]
Lseries_ftarget = Lseries[ftarget_index]
Cshunt1_ftarget = Cshunt1[ftarget_index]
Cshunt2_ftarget = Cshunt2[ftarget_index]
Rshunt1_ftarget = Rshunt1[ftarget_index]
Rshunt2_ftarget = Rshunt2[ftarget_index]

Qdiff_ftarget = Qdiff[ftarget_index]

append_log('\nDifferential inductor parameters')
append_log(f"Effective series L  [nH] : {Ldiff[ftarget_index]*1e9:.3f} @ {f[ftarget_index]/1e9:.3f} GHz")  
append_log(f"Effective series R  [Ohm]: {Rdiff[ftarget_index]:.3f} @ {f[ftarget_index]/1e9:.3f} GHz") 
append_log(f"Differential Q factor    : {Qdiff[ftarget_index]:.2f} @ {f[ftarget_index]/1e9:.3f} GHz")  
append_log('----------------------')
append_log(f"L_DC      [nH] : {Ldiff[1]*1e9:.3f}") 
append_log(f"R_DC      [Ohm]: {Rdiff[0]:.3f}")  
append_log(f"Peak Q         : {max(Qdiff):.2f} @ {f_Qmax/1e9:.3f} GHz") 
append_log('')
append_log(f"\nPi model extraction (narrowband) at {f[ftarget_index]/1e9:.3f} GHz")
append_log(f"Series L  [nH] : {Lseries_ftarget*1e9:.3f}")  
append_log(f"Series R  [Ohm]: {Rseries_ftarget:.3f}") 
append_log(f"Shunt C @ port 1 [fF] : {Cshunt1_ftarget*1e15:.3f}")  
append_log(f"Shunt R @ port 1 [Ohm]: {Rshunt1_ftarget:.3f}")  
append_log(f"Shunt C @ port 2 [fF] : {Cshunt2_ftarget*1e15:.3f}")  
append_log(f"Shunt R @ port 2 [Ohm]: {Rshunt2_ftarget:.3f}")  
append_log('')

log_text = "".join(log)

# append_log log to terminal
print(log_text)

# output log message to file also, same basename as *.s2p input file but file extension .txt
log_filename = os.path.splitext(args.s2p)[0] + '.txt'
with open(log_filename, "w", encoding="utf-8") as f:
    f.write(log_text)


plt.figure()

plt.subplot(121)
plt.plot(freq/1e9, Ldiff*1e9,label='Lseries [nH]')
plt.xlabel('f (GHz)')
plt.legend()

plt.subplot(122)
plt.plot(freq/1e9, Qdiff,label='Q factor')
plt.xlabel('f (GHz)')
plt.legend()

plt.figure()

plt.subplot(221)
plt.plot(freq/1e9, Lseries*1e9,label='Lseries [nH]')
plt.plot(f_target/1e9, Lseries_ftarget*1e9, 'ro', label='fit')
plt.xlabel('f (GHz)')
plt.ylim(0, 2*Lseries_ftarget*1e9)
plt.legend()

plt.subplot(222)
plt.plot(freq/1e9, Rseries, label='Rseries [Ohm]')
plt.plot(f_target/1e9, Rseries_ftarget,'ro', label='fit')
plt.xlabel('f (GHz)')
plt.ylim(0, 2*Rseries_ftarget)
plt.legend()

plt.subplot(223)
plt.plot(freq/1e9, Cshunt1, 'r', label='Cshunt1 [fF]')
plt.plot(freq/1e9, Cshunt2, 'g', label='Cshunt2 [fF]')
plt.plot(f_target/1e9, Cshunt1_ftarget,'ro', label='fit' )
plt.plot(f_target/1e9, Cshunt2_ftarget,'go', label='fit' )
plt.xlabel('f (GHz)')
plt.ylim(0, 5*Cshunt1_ftarget)
plt.legend()

plt.subplot(224)
plt.plot(freq/1e9, Rshunt1, 'r', label='Rshunt1 [Ohm]')
plt.plot(freq/1e9, Rshunt2, 'g', label='Rshunt2 [Ohm]')
plt.plot(f_target/1e9, Rshunt1_ftarget,'ro', label='fit' )
plt.plot(f_target/1e9, Rshunt2_ftarget,'go', label='fit' )
plt.xlabel('Frequency (GHz)')
plt.ylim(0, 5*Rshunt1_ftarget)
plt.legend()

plt.show()

