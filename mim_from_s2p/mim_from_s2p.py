# extract simple pi model, single frequency exctration, not accurate for wideband use
# data reader and extraction based on scikit-rf functionality

import skrf as rf
import math
import argparse
from matplotlib import pyplot as plt

print('Extract MIM capacitor model from S2P S-parameter file')

# evaluate commandline
parser = argparse.ArgumentParser()
parser.add_argument("s2p",  help="S2P input filename (Touchstone format)")
parser.add_argument("f_ghz", help="extraction frequency in GHz", type=float)
args = parser.parse_args()


# input data, must be 2-port S2P data
sub = rf.Network(args.s2p)

# if the input data has DC point, remove that because it will throw warnungs later
if sub.frequency.start == 0:
    # resample to start at 1 GHz (or closest value)
    newrange = '1-' + str(sub.frequency.stop/1e9) + 'ghz'
    sub = sub[newrange]


# target frequency for pi model extraction
f_target = args.f_ghz*1e9


# frequency class, see https://github.com/scikit-rf/scikit-rf/blob/master/skrf/frequency.py
freq = sub.frequency
print('S2P frequency range is ',freq.start/1e9, ' to ', freq.stop/1e9, ' GHz')
assert f_target < freq.stop

f = freq.f

# get a target frequency and low frequency for extraction
f_low = max(f_target/20, 1e9)
ftarget_index = rf.find_nearest_index(f, f_target)
flow_index = rf.find_nearest_index(f, f_low)

omega = f*2*math.pi

# to evaluate target freq, we separate series and shunt path values

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

Cser = 1/(-Zseries.imag*omega)
Cser_target = Cser[ftarget_index]  # target frequency value
Cser_low =  Cser[flow_index]       # low frequency value, where Lseries has very small effect on X

print(f"Choosing  {f[flow_index]/1e9:.1f} GHz for low frequency data")  
print('_________________________________________________________')
print(f"Series C extracted at low frequency: {Cser_low*1e15:.3f} fF")  
# print(f"Series C extracted at {f[ftarget_index]/1e9:.1f} GHz: : {Cser_target*1e15:.3f} fF")  

# calculate series L
Lser = (Zseries.imag + 1 / (omega*Cser_low))/omega
Lser_target = Lser[ftarget_index]
print(f"Series L extracted at {f[ftarget_index]/1e9:.1f} GHz: {Lser_target*1e12:.3f} pH")  

# calculate series R
Rser_target = Zseries.real[ftarget_index]
print(f"Series R extracted at {f[ftarget_index]/1e9:.1f} GHz: {Rser_target:.3f} Ohm")  


# check  effective series impedance from calculated L and C
X = omega*Lser_target - 1 / (omega * Cser_low)
C_from_X = 1/(-X*omega)

# shunt capacitance
Cshunt1 = -1 / (omega*Zshunt1.imag)
Cshunt2 = -1 / (omega*Zshunt2.imag)

# distribute shunt C equally to both sides, because we can't separate precisely
# the real MIM will have different Cshunt on the two ports, but we don't care 
# because they are connected tightly by the large MIM value anyway
Cshunt_target = (Cshunt1[ftarget_index] + Cshunt2[ftarget_index])/2

print(f"Shunt C distributed to both ports:  {Cshunt_target*1e15:.3f} fF each side")  
print('_________________________________________________________')


plt.figure()

plt.plot(f/1e9, Cser*1e15,label='Cseries [fF]')
plt.plot(f_low/1e9, Cser_low*1e15,'go', label = f"{Cser_low*1e15:.1f} fF @ {f_low/1e9:.1f} GHz")
plt.plot(f_target/1e9, Cser_target*1e15,'ro',label = f"{Cser_target*1e15:.1f} fF @ {f_target/1e9:.1f} GHz")
plt.plot(f/1e9, C_from_X*1e15,'r:', label='Model fit')
plt.xlim([0, max(f/1e9)])
plt.ylim([0, max(2*Cser_low*1e15,1.5*Cser_target*1e15)])
plt.xlabel('f (GHz)')
plt.ylabel('Cseries (fF)')
plt.title('Effective series C, input data vs. model')
plt.grid()
plt.legend()


plt.show()
