# extract simple pi model, single frequency exctration, not accurate for wideband use
# data reader and extraction based on scikit-rf functionality

import skrf as rf
import math
import argparse
from matplotlib import pyplot as plt

print('Extract simple inductor pi model from S2P S-parameter file')

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
print('S2P frequency range is ',sub.frequency.start/1e9, ' to ', sub.frequency.stop/1e9, ' GHz')
print('Extraction frequency: ', args.f_ghz, ' GHz')
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

print('\nDifferential inductor parameters')
print(f"Effective series L  [nH] : {Ldiff[ftarget_index]*1e9:.3f} @ {f[ftarget_index]/1e9:.3f} GHz")  
print(f"Effective series R  [Ohm]: {Rdiff[ftarget_index]:.3f} @ {f[ftarget_index]/1e9:.3f} GHz") 
print(f"Differential Q factor    : {Qdiff[ftarget_index]:.2f} @ {f[ftarget_index]/1e9:.3f} GHz")  
print('----------------------')
print(f"L_DC      [nH] : {Ldiff[1]*1e9:.3f}") 
print(f"R_DC      [Ohm]: {Rdiff[0]:.3f}")  
print(f"Peak Q         : {max(Qdiff):.2f} @ {f_Qmax/1e9:.3f} GHz") 
print('')
print(f"\nPi model extraction (narrowband) at {f[ftarget_index]/1e9:.3f} GHz")
print(f"Series L  [nH] : {Lseries_ftarget*1e9:.3f}")  
print(f"Series R  [Ohm]: {Rseries_ftarget:.3f}") 
print(f"Shunt C @ port 1 [fF] : {Cshunt1_ftarget*1e15:.3f}")  
print(f"Shunt R @ port 1 [Ohm]: {Rshunt1_ftarget:.3f}")  
print(f"Shunt C @ port 2 [fF] : {Cshunt2_ftarget*1e15:.3f}")  
print(f"Shunt R @ port 2 [Ohm]: {Rshunt2_ftarget:.3f}")  
print('')


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


# ===================== LOGGING WRAPPER =====================
# This block mirrors all stdout to a text file whose basename matches
# the input S2P file. It does NOT modify any existing logic above.
# Safe to remove without affecting functionality.

import sys
from pathlib import Path
from contextlib import redirect_stdout

# Prevent infinite recursion
if "__LOGGING_ACTIVE__" not in globals():

    __LOGGING_ACTIVE__ = True

    if len(sys.argv) > 1 and sys.argv[1].lower().endswith(".s2p"):
        log_path = Path(sys.argv[1]).with_suffix(".txt")

        class _Tee:
            def __init__(self, *streams):
                self.streams = streams
            def write(self, data):
                for s in self.streams:
                    s.write(data)
            def flush(self):
                for s in self.streams:
                    s.flush()

        with open(log_path, "w") as _logfile:
            _tee = _Tee(sys.stdout, _logfile)
            with redirect_stdout(_tee):
                exec(open(__file__).read(), globals(), globals())
        sys.exit(0)



# ===================== AUTO-OPEN LOG FILE =====================
# Open the generated log file in gedit after script execution.
# This runs only in the logging wrapper execution path.

import subprocess
try:
    subprocess.Popen(["gedit", str(log_path)])
except Exception:
    pass
# =================== END AUTO-OPEN ===================

# =================== END LOGGING WRAPPER ===================
