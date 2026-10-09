# Copyright (C) 2026 Volker Muehlhaus <volker@muehlhaus.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# extract simple pi model, single frequency exctration, not accurate for wideband use
# data reader and extraction based on scikit-rf functionality

import skrf as rf
from skrf.util import find_nearest_index
import math
import argparse
from matplotlib import pyplot as plt
import os
import numpy as np
from scipy.optimize import least_squares, fsolve

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
Qmax_index = find_nearest_index(Qdiff, Qmax)
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

ftarget_index = find_nearest_index(freq, f_target)
Rseries_ftarget = Rseries[ftarget_index]
Lseries_ftarget = Lseries[ftarget_index]
Cshunt1_ftarget = Cshunt1[ftarget_index]
Cshunt2_ftarget = Cshunt2[ftarget_index]
Rshunt1_ftarget = Rshunt1[ftarget_index]
Rshunt2_ftarget = Rshunt2[ftarget_index]


def y_series_model(params, omega):
    R, L, C = params
    return 1.0 / (R + 1j * omega * L) + 1j * omega * C

def y_series_residuals(params, omega, y_data):
    diff = y_series_model(params, omega) - y_data
    return np.concatenate([diff.real, diff.imag])

def z_skin_branch(Rskin, Lskin, omega):
    # Parallel Rskin||Lskin: Zskin -> 0 at DC (Lskin shorts), Zskin -> Rskin at high frequency
    # (Lskin opens), adding an AC resistance rise on top of the DC series R.
    return (Rskin * 1j * omega * Lskin) / (Rskin + 1j * omega * Lskin)

def y_series_model_skin(params, omega):
    R, L, Rskin, Lskin, C = params
    Zbranch = R + 1j * omega * L + z_skin_branch(Rskin, Lskin, omega)
    return 1.0 / Zbranch + 1j * omega * C

def y_series_residuals_skin(params, omega, y_data):
    diff = y_series_model_skin(params, omega) - y_data
    return np.concatenate([diff.real, diff.imag])

# ---- Wideband series-branch fit ----
Yseries_data = -ymn  # = 1/Zseries

# A single R-L(+skin)||Csrf branch can only represent one resonance. For a very wideband input,
# the raw data may extend well past that self resonance into behaviour this topology can't
# represent, and trying to fit that region too just pollutes the fit everywhere else. Find where
# the RAW (unmodeled) Lseries data itself first goes from inductive to capacitive, and restrict
# the fit to data at or below that frequency - the same self resonance crossing the narrowband
# section already computes implicitly, just located directly from data instead of from a model.
raw_crossings = np.where(np.diff(np.sign(Lseries)) < 0)[0]
if len(raw_crossings) > 0:
    i0 = raw_crossings[0]
    f_lo, f_hi = f[i0], f[i0 + 1]
    L_lo, L_hi = Lseries[i0], Lseries[i0 + 1]
    f_srf_raw = f_lo + (f_hi - f_lo) * L_lo / (L_lo - L_hi)
    fit_mask = f <= f_srf_raw
else:
    f_srf_raw = None
    fit_mask = np.ones_like(f, dtype=bool)

# Also cap the fit band at 2x the extraction frequency, on top of (or instead of, if no self
# resonance was found) the self-resonance-based cap above - keeps the fit focused on the region
# the user actually cares about rather than being pulled around by far-off-band behaviour.
f_hi_cutoff = 2 * f_target
fit_mask_hi = f <= f_hi_cutoff
if np.count_nonzero(fit_mask & fit_mask_hi) >= 5:
    fit_mask = fit_mask & fit_mask_hi

# Also drop data well below the extraction frequency: very low frequency points mainly pin down
# the DC resistance/inductance, which is not what a "wideband fit near f_target" should optimize
# for, and including them just dilutes the fit's focus on the region the user actually cares about.
f_lo_cutoff = 0.1 * f_target
fit_mask_lo = f >= f_lo_cutoff
if np.count_nonzero(fit_mask & fit_mask_lo) >= 5:
    fit_mask = fit_mask & fit_mask_lo
else:
    # not enough points left (e.g. f_target very close to self resonance) - skip the low cutoff
    # rather than fail with too few points to fit
    f_lo_cutoff = None

omega_fit = omega[fit_mask]
Yseries_data_fit = Yseries_data[fit_mask]

# -- Stage 1: R-L||Csrf fit (no skin effect) - used only to seed stage 2, not logged/plotted --
R0 = Rseries_ftarget
L0 = Lseries_ftarget
C0 = 1e-15  # generic small seed; dY/dC = j*omega is nonzero everywhere since DC point is trimmed

fit_result_stage1 = least_squares(
    y_series_residuals, x0=[R0, L0, C0], args=(omega_fit, Yseries_data_fit),
    bounds=([1e-9, 1e-15, 0.0], [np.inf, np.inf, np.inf]),
    method='trf', x_scale='jac',  # auto-scale: R~Ohm, L~1e-9H, C~1e-15F differ by orders of magnitude
)
Rseries_fit_stage1, Lseries_fit_stage1, Csrf_fit_stage1 = fit_result_stage1.x

# -- Stage 2: add a skin-effect section (Rskin||Lskin), warm-started from stage 1.
#    Rskin0 assumes an AC resistance rise roughly comparable to the DC resistance as a generic
#    starting point. Lskin0 places the skin section's corner frequency (Rskin/(2*pi*Lskin)) at
#    f_target - the one "this frequency matters" datum the script has from the user. --
Rskin0 = Rseries_fit_stage1
Lskin0 = Rskin0 / (2 * math.pi * f_target)

fit_result_skin = least_squares(
    y_series_residuals_skin,
    x0=[Rseries_fit_stage1, Lseries_fit_stage1, Rskin0, Lskin0, Csrf_fit_stage1],
    args=(omega_fit, Yseries_data_fit),
    # Lskin is capped at 1uH (>1000x a typical Lseries) so a skin corner below the measured band
    # converges to a bounded, sane Lskin instead of running away to an unphysical mH-scale value
    bounds=([1e-9, 1e-15, 1e-9, 1e-15, 0.0], [np.inf, np.inf, np.inf, 1e-6, np.inf]),
    method='trf', x_scale='jac',
)
Rseries_fit, Lseries_fit, Rskin_fit, Lskin_fit, Csrf_fit = fit_result_skin.x

# -- Calibration: adjust Rseries and Lseries (DC resistance and inductance) together so the model
#    matches both closely at the requested extraction frequency, which is what the user actually
#    cares about matching - while leaving Rskin/Lskin/Csrf exactly as found by the wideband fit above,
#    so the rest of the curve's shape (and self resonance) is undisturbed. This is a well-posed
#    2-unknown/2-equation problem (match Re(Z) and Im(Z) at one frequency), solved jointly rather
#    than as two independent 1D solves, since adjusting either parameter shifts both Re(Z) and
#    Im(Z) at once. fsolve (Newton-based root finding for exactly-determined systems) is used
#    instead of least_squares, since least_squares's finite-difference Jacobian estimate is poorly
#    scaled here (R ~ Ohms, L ~ 1e-9 H) and was found to barely move from the starting point. --
Ztarget = Rseries_ftarget + 1j * omega[ftarget_index] * Lseries_ftarget

def _z_model_residual_at_ftarget(params):
    Rseries_val, Lseries_val = params
    full_params = (Rseries_val, Lseries_val, Rskin_fit, Lskin_fit, Csrf_fit)
    Zm = 1.0 / y_series_model_skin(full_params, omega[ftarget_index])
    return [Zm.real - Ztarget.real, Zm.imag - Ztarget.imag]

Rseries_fit, Lseries_fit = fsolve(_z_model_residual_at_ftarget, x0=[Rseries_fit, Lseries_fit])
Rseries_fit = max(Rseries_fit, 1e-9)   # keep physically non-negative
Lseries_fit = max(Lseries_fit, 1e-15)

fitted_params = (Rseries_fit, Lseries_fit, Rskin_fit, Lskin_fit, Csrf_fit)

# RMS error is computed from the calibrated parameter set, over the same sub-band that was
# actually fit (at/below the data's self resonance) - including data the fit was deliberately
# not asked to explain would just make this number meaningless
raw_diff = y_series_model_skin(fitted_params, omega_fit) - Yseries_data_fit
rms_error_Y = np.sqrt(np.mean(raw_diff.real ** 2 + raw_diff.imag ** 2))
rms_error_Y_rel = rms_error_Y / np.sqrt(np.mean(np.abs(Yseries_data_fit) ** 2)) * 100

Rseries_hf = Rseries_fit + Rskin_fit           # asymptotic high-frequency (AC) resistance
f_corner = Rskin_fit / (2 * math.pi * Lskin_fit) # skin-effect onset ("corner") frequency
# if the corner sits outside the measured band, the skin section's transition never happens
# in-band, so the Rdc/Rskin split is unconstrained by the data (many (Rdc,Rskin,Lskin) combinations
# give the same in-band fit) - only their sum Rseries_hf is then a meaningful number
skin_corner_in_range = sub.frequency.start <= f_corner <= sub.frequency.stop

# ---- Self resonant frequency (SRF): numeric zero-crossing of Im(Y_model) ----
# The skin section makes Z_branch(omega) a more complex rational function of omega than the
# simple R-L||C case, so an exact closed-form solve is no longer attractive. Evaluate the fitted
# model on a dense grid extending beyond the measured band (to allow for extrapolated SRF) and
# find the first upward zero crossing of Im(Y_model) - the same physical crossing the old
# closed-form solution located - with linear interpolation between grid points for accuracy.
f_grid = np.linspace(sub.frequency.start, sub.frequency.stop * 1.5, 20000)
im_Y_grid = y_series_model_skin(fitted_params, f_grid * 2 * math.pi).imag
crossings = np.where(np.diff(np.sign(im_Y_grid)) > 0)[0]
if len(crossings) > 0:
    i0 = crossings[0]
    f_lo, f_hi = f_grid[i0], f_grid[i0 + 1]
    im_lo, im_hi = im_Y_grid[i0], im_Y_grid[i0 + 1]
    SRF = f_lo + (f_hi - f_lo) * (-im_lo) / (im_hi - im_lo)
else:
    SRF = float('nan')
srf_in_range = (not math.isnan(SRF)) and (sub.frequency.start <= SRF <= sub.frequency.stop)

Zseries_fit_curve = 1.0 / y_series_model_skin(fitted_params, omega)
Lseries_fit_curve = Zseries_fit_curve.imag / omega
Rseries_fit_curve = Zseries_fit_curve.real

# how well the wideband model (including Csrf loading) actually matches Rseries/Lseries right at
# the requested extraction frequency - should equal the narrowband targets closely after the
# calibration above
Rseries_fit_at_ftarget = Rseries_fit_curve[ftarget_index]
Lseries_fit_at_ftarget = Lseries_fit_curve[ftarget_index]

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

fit_band_start = f[fit_mask][0]
fit_band_stop = f[fit_mask][-1]
append_log(f"\nWideband series branch fit (R + skin Rskin||Lskin + Csrf, fit band {fit_band_start/1e9:.3f}-{fit_band_stop/1e9:.3f} GHz, capped at 2x extraction freq = {f_hi_cutoff/1e9:.3f} GHz)")
append_log(f"Series L (DC)              [nH] : {Lseries_fit*1e9:.3f}")
append_log(f"Series R (DC)              [Ohm]: {Rseries_fit:.3f}")
append_log(f"Skin effect R (Rskin)      [Ohm]: {Rskin_fit:.3f}")
append_log(f"Skin effect L (Lskin)      [nH] : {Lskin_fit*1e9:.3f}")
append_log(f"Csrf (parallel)            [fF] : {Csrf_fit*1e15:.3f}")
append_log('')
append_log(f"Self resonant frequency (SRF) [GHz]: {SRF/1e9:.3f}")
append_log(f"Effective series L at {args.f_ghz:.3f} GHz (model)  [nH] : {Lseries_fit_at_ftarget*1e9:.3f} (narrowband target: {Lseries_ftarget*1e9:.3f})")
append_log(f"Effective series R at {args.f_ghz:.3f} GHz (model)  [Ohm]: {Rseries_fit_at_ftarget:.3f} (narrowband target: {Rseries_ftarget:.3f})")

if not srf_in_range:
    append_log(f"WARNING: fitted SRF ({SRF/1e9:.3f} GHz) is outside the measured frequency range "
                f"({sub.frequency.start/1e9:.3f}-{sub.frequency.stop/1e9:.3f} GHz) - extrapolated, less trustworthy")
if not skin_corner_in_range:
    append_log(f"WARNING: skin effect corner frequency ({f_corner/1e9:.3f} GHz) is outside the measured "
                f"frequency range ({sub.frequency.start/1e9:.3f}-{sub.frequency.stop/1e9:.3f} GHz) - the Rskin/Lskin "
                f"split is underdetermined for this dataset, only their sum (effective high-frequency series R = "
                f"{Rseries_hf:.3f} Ohm) is meaningful")
if f_srf_raw is not None and f_target > f_srf_raw:
    append_log(f"WARNING: extraction frequency ({args.f_ghz:.3f} GHz) is above the data's self resonance "
                f"({f_srf_raw/1e9:.3f} GHz) - the wideband series branch model was fit below self resonance only, "
                f"so this region is extrapolated")
append_log('')

log_text = "".join(log)

# append_log log to terminal
print(log_text)

# output log message to file also, same basename as *.s2p input file but file extension .txt
log_filename = os.path.splitext(args.s2p)[0] + '.txt'
with open(log_filename, "w", encoding="utf-8") as f:
    f.write(log_text)


plt.figure(figsize=(12.8, 9.6))  # 2x the matplotlib default (6.4, 4.8)

plt.subplot(121)
plt.plot(freq/1e9, Ldiff*1e9,label='Lseries [nH]')
plt.xlabel('f (GHz)')
plt.legend()
plt.grid(color='lightgray')

plt.subplot(122)
plt.plot(freq/1e9, Qdiff,label='Q factor')
plt.xlabel('f (GHz)')
plt.legend()
plt.grid(color='lightgray')

plt.figure(figsize=(12.8, 9.6))  # 2x the matplotlib default (6.4, 4.8)

plt.subplot(221)
plt.plot(freq/1e9, Lseries*1e9,label='Lseries [nH]')
plt.plot(freq/1e9, Lseries_fit_curve*1e9, 'b--', label='wideband fit')
plt.plot(f_target/1e9, Lseries_ftarget*1e9, 'ro', label='narrowband fit')
if not math.isnan(SRF):
    plt.axvline(SRF/1e9, color='m', linestyle=':', label=f'SRF = {SRF/1e9:.2f} GHz')
plt.xlabel('f (GHz)')
plt.ylim(0, 2*Lseries_ftarget*1e9)
plt.legend()
plt.grid(color='lightgray')

plt.subplot(222)
plt.plot(freq/1e9, Rseries, label='Rseries [Ohm]')
plt.plot(freq/1e9, Rseries_fit_curve, 'b--', label='wideband fit')
plt.plot(f_target/1e9, Rseries_ftarget,'ro', label='narrowband fit')
if not math.isnan(SRF):
    plt.axvline(SRF/1e9, color='m', linestyle=':', label=f'SRF = {SRF/1e9:.2f} GHz')
plt.xlabel('f (GHz)')
plt.ylim(0, max(2*Rseries_ftarget, 1.2*Rseries_hf))
plt.legend()
plt.grid(color='lightgray')

plt.subplot(223)
plt.plot(freq/1e9, Cshunt1, 'r', label='Cshunt1 [fF]')
plt.plot(freq/1e9, Cshunt2, 'g', label='Cshunt2 [fF]')
plt.plot(f_target/1e9, Cshunt1_ftarget,'ro', label='Cshunt1 narrowband fit' )
plt.plot(f_target/1e9, Cshunt2_ftarget,'go', label='Cshunt2 narrowband fit' )
plt.xlabel('f (GHz)')
plt.ylim(0, 5*Cshunt1_ftarget)
plt.legend()
plt.grid(color='lightgray')

plt.subplot(224)
plt.plot(freq/1e9, Rshunt1, 'r', label='Rshunt1 [Ohm]')
plt.plot(freq/1e9, Rshunt2, 'g', label='Rshunt2 [Ohm]')
plt.plot(f_target/1e9, Rshunt1_ftarget,'ro', label='Rshunt1 narrowband fit' )
plt.plot(f_target/1e9, Rshunt2_ftarget,'go', label='Rshunt2 narrowband fit' )
plt.xlabel('Frequency (GHz)')
plt.ylim(0, 5*Rshunt1_ftarget)
plt.legend()
plt.grid(color='lightgray')

plt.show()

