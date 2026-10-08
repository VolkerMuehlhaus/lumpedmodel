# Wideband RFIC inductor model fit from S2P data, including the substrate network.
# Fully automatic: analytic seed values from the pi decomposition of the Y-parameters,
# followed by a global least-squares fit of the complete 2-port model over the band up to
# (and a bit beyond) the self resonance.
#
# Model topology (up to 15 elements, schematic in doc/inductor_fit_model.png), 1 coil segment:
#
#   p1 --+-- Rs -- Ls -- (Rskin1 || Lskin1) -- (Rskin2 || Lskin2) --+-- p2
#        |                                                          |
#        +---------------------------- Cs --------------------------+   interwinding/underpass C
#
#   p1 -- Cox1 -- s1                s2 -- Cox2 -- p2  oxide capacitance
#                 |                  |
#           Rsi1 || Csi1       Rsi2 || Csi2           bulk silicon, to ground
#                 |                  |
#                 s1 -- Rsub12 || Csub12 -- s2        substrate coupling (only if needed)
#
# Two skin sections (Rskin||Lskin "rungs") with different corner frequencies model the
# frequency dependent R and L from skin and proximity effect / current redistribution.
# For electrically long coils (mm-wave inductors measured up to and beyond SRF), the coil is
# split into several segments with the substrate network distributed along the coil, which
# approximates the distributed (transmission line) behaviour. The element values are always
# reported as totals; the number of segments is chosen automatically.
# The substrate network is kept physically meaningful and as simple as possible: symmetric if
# the data is symmetric, substrate coupling only if it is physically plausible and needed.
# Option --basic always uses 1 segment, 1 skin section and no substrate coupling, for a
# consistent topology across many fits, e.g. when generating data tables for machine learning.

import skrf as rf
import math
import argparse
import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import least_squares

# weights of the error function terms used in the global fit
W_Y = 1.0        # complex Y11, Y22 and series branch -Y12 (relative error)
W_SHUNT = 1.0    # complex shunt branches Y11+Y12 and Y22+Y12 (relative error)
W_SHUNT_RE = 0.5 # real part of the shunt branches = substrate loss (relative error)
W_L = 0.5        # single-ended L11/L22 and differential L (relative error)
W_Q = 0.5        # single-ended Q11/Q22 and differential Q (relative error)

SRF_FIT_FACTOR = 1.2  # automatic fit band stops at this factor times the measured SRF
LQ_SRF_FACTOR = 0.8   # L and Q goals are only used below this factor times SRF (L, Q diverge at SRF)
NUM_RANDOM_STARTS = 2 # additional global fit runs from randomly perturbed seed values
MAX_NFEV = 200        # cap per global fit run - the last iterations crawl along flat valleys without real gain
SEGMENT_CHOICES = (1, 2, 3)   # number of coil segments tried in automatic mode
SEGMENT_COST_TOLERANCE = 1.1  # use the fewest segments with cost within this factor of the best

# Substrate network variants (symmetric = port 2 values tied to port 1, coupling = Rsub12||Csub12),
# ordered from simple to complex. In automatic mode, symmetry is decided from the data. With
# coupling, the substrate network is not always unique: a "shared substrate" solution can fit as
# well, but is not physically meaningful. Coupling is only used if it passes a physical
# plausibility check and improves the fit cost by more than this factor.
SUBSTRATE_VARIANTS = {'symmetric': (True, False), 'asymmetric': (False, False),
                      'symmetric+coupling': (True, True), 'asymmetric+coupling': (False, True)}
SUBSTRATE_COST_TOLERANCE = 1.1
SYMMETRY_TOLERANCE = 0.05  # symmetric substrate network if the shunt branches differ by less (RMS)
COX_SPAN = 0.5  # Cox bound in decades around the low frequency shunt capacitance (seed value)

PARAM_NAMES = ['Rs', 'Ls', 'Rskin1', 'Lskin1', 'Rskin2', 'Lskin2', 'Cs',
               'Cox1', 'Rsi1', 'Csi1', 'Cox2', 'Rsi2', 'Csi2',
               'Rsub12', 'Csub12']

# create a log that we can dump to terminal and log file
log = []
def append_log (txt):
    log.append(txt + '\n')


# ---------------- model equations ----------------

def z_skin_branch(Rskin, Lskin, omega):
    # Parallel Rskin||Lskin: Zskin -> jwLskin at low frequency (adds internal inductance) and
    # Zskin -> Rskin at high frequency (adds AC resistance) - models skin/proximity effect
    return (Rskin * 1j * omega * Lskin) / (Rskin + 1j * omega * Lskin)

def z_coil(Rs, Ls, Rskin1, Lskin1, Rskin2, Lskin2, omega):
    return (Rs + 1j * omega * Ls + z_skin_branch(Rskin1, Lskin1, omega)
            + z_skin_branch(Rskin2, Lskin2, omega))

def y_shunt_branch(Cox, Rsi, Csi, omega):
    # Cox in series with Rsi||Csi
    yox = 1j * omega * Cox
    ysi = 1.0 / Rsi + 1j * omega * Csi
    return yox * ysi / (yox + ysi)

def shunt_weights(segments):
    # Share of the port 1 (a) and port 2 (b) shunt elements at each of the segments+1 coil nodes:
    # trapezoidal weights along the coil, linearly blended from the port 1 to the port 2 values.
    # segments=1 puts the port 1 shunt network fully at p1 and the port 2 network fully at p2.
    k = np.arange(segments + 1)
    trapezoid = np.where((k == 0) | (k == segments), 0.5, 1.0) / segments
    return 2 * trapezoid * (1 - k / segments), 2 * trapezoid * k / segments

def segment_shunt_values(p, segments):
    # Cox, Rsi, Csi at each coil node, from the total port 1 / port 2 values
    Cox1, Rsi1, Csi1, Cox2, Rsi2, Csi2 = p[7:13]
    a, b = shunt_weights(segments)
    return a * Cox1 + b * Cox2, 1.0 / (a / Rsi1 + b / Rsi2), a * Csi1 + b * Csi2

def y_model(p, omega, segments):
    # 2-port Y matrix (shape N,2,2) of the full model. The coil (Rs, Ls, skin sections) is split
    # into 'segments' equal cells, with the Cox-(Rsi||Csi) substrate network distributed over the
    # coil nodes (see shunt_weights). Cs connects p1-p2, Rsub12||Csub12 connects the substrate
    # nodes under p1 and p2. Nodal admittance matrix with coil nodes 0..segments (0 = p1,
    # segments = p2) and one substrate node per coil node, all internal nodes eliminated by Kron
    # reduction: Y = Ypp - Ypi * inv(Yii) * Yip. Also valid at DC (Cox open, inductors short).
    Cs, Rsub12, Csub12 = p[6], p[13], p[14]
    num_coil = segments + 1
    num_nodes = 2 * num_coil
    M = np.zeros((len(omega), num_nodes, num_nodes), dtype=complex)

    def stamp(n1, n2, y):
        M[:, n1, n1] += y
        if n2 is not None:  # None = ground
            M[:, n2, n2] += y
            M[:, n1, n2] -= y
            M[:, n2, n1] -= y

    y_cell = segments / z_coil(*p[:6], omega)
    for k in range(segments):
        stamp(k, k + 1, y_cell)
    stamp(0, segments, 1j * omega * Cs)

    Cox, Rsi, Csi = segment_shunt_values(p, segments)
    for k in range(num_coil):
        stamp(k, num_coil + k, 1j * omega * Cox[k])
        stamp(num_coil + k, None, 1.0 / Rsi[k] + 1j * omega * Csi[k])
    stamp(num_coil, num_coil + segments, 1.0 / Rsub12 + 1j * omega * Csub12)

    ports = [0, segments]
    internal = [n for n in range(num_nodes) if n not in ports]
    Mpp = M[:, ports][:, :, ports]
    Mpi = M[:, ports][:, :, internal]
    Mip = M[:, internal][:, :, ports]
    Mii = M[:, internal][:, :, internal]
    return Mpp - Mpi @ np.linalg.solve(Mii, Mip)


# ---------------- derived quantities (same for measured data and model) ----------------

def pi_branches(Y):
    # exact pi decomposition of any 2-port: shunt branch at port 1, shunt branch at port 2, series branch
    y12 = (Y[:, 0, 1] + Y[:, 1, 0]) / 2
    return Y[:, 0, 0] + y12, Y[:, 1, 1] + y12, -y12

def y_diff(Y):
    # differential admittance between port 1 and port 2 with floating ground,
    # equal to 1/(z11 - z12 - z21 + z22) but without inverting a near-singular Y at low frequency
    det = Y[:, 0, 0] * Y[:, 1, 1] - Y[:, 0, 1] * Y[:, 1, 0]
    return det / (Y[:, 0, 0] + Y[:, 1, 1] + Y[:, 0, 1] + Y[:, 1, 0])

def l_and_q(Yin, omega):
    Zin = 1.0 / Yin
    return Zin.imag / omega, Zin.imag / Zin.real

def find_srf(f, Yin):
    # first frequency where Im(Yin) goes from inductive (negative) to capacitive (positive),
    # linear interpolation between data points
    crossings = np.where(np.diff(np.sign(Yin.imag)) > 0)[0]
    if len(crossings) == 0:
        return None
    i0 = crossings[0]
    im_lo, im_hi = Yin.imag[i0], Yin.imag[i0 + 1]
    return f[i0] + (f[i0 + 1] - f[i0]) * (-im_lo) / (im_hi - im_lo)


# ---------------- error function helpers ----------------

def floored_norm(x):
    # per-frequency normalization |x|, floored at 10% of the RMS value, so that relative errors
    # at points where the data passes through (or near) zero don't dominate the fit
    mag = np.abs(x)
    return np.maximum(mag, 0.1 * np.sqrt(np.mean(mag ** 2)))

def complex_residual(model, meas, norm, weight):
    diff = (model - meas) / norm
    return weight * np.concatenate([diff.real, diff.imag]) / math.sqrt(len(meas))

def real_residual(model, meas, norm, weight):
    return weight * (model - meas) / norm / math.sqrt(len(meas))


# evaluate commandline
parser = argparse.ArgumentParser(description='Wideband RFIC inductor model fit (series L/R/skin/Cs + Cox-Rsi||Csi substrate network)')
parser.add_argument("s2p",  help="S2P input filename (Touchstone format)")
parser.add_argument("--fmin", help="lower end of fit band in GHz (default: lowest non-DC data point)", type=float)
parser.add_argument("--fmax", help=f"upper end of fit band in GHz (default: {SRF_FIT_FACTOR} x measured SRF)", type=float)
parser.add_argument("--segments", help=f"number of coil segments (default: automatic choice from {SEGMENT_CHOICES})", type=int)
parser.add_argument("--substrate", help="substrate network (default: auto = simplest physical variant that fits)",
                    choices=['auto'] + list(SUBSTRATE_VARIANTS), default='auto')
parser.add_argument("--basic", help="basic model with fixed topology: 1 coil segment, 1 skin section, no substrate coupling", action='store_true')
parser.add_argument("--noplot", help="don't show plots", action='store_true')
args = parser.parse_args()
if args.basic:
    assert args.segments in (None, 1), '--basic model always uses 1 coil segment'
    assert not SUBSTRATE_VARIANTS.get(args.substrate, (False, False))[1], '--basic model has no substrate coupling'
    args.segments = 1

# In the basic model, skin section 2 is fixed at a negligible value (Lskin2 -> 0 shorts it) and
# excluded from fitting. All fits work on log10 of the parameter values.
num_skin_sections = 1 if args.basic else 2
series_free = np.array([num_skin_sections == 2 or name not in ('Rskin2', 'Lskin2') for name in PARAM_NAMES[:7]])
X_SKIN2_FIXED = (0.0, -30.0)  # log10 of Rskin2 = 1 Ohm, Lskin2 = 1e-30 H


# input data, must be 2-port S2P data
sub_full = rf.Network(args.s2p)
assert sub_full.nports == 2, 'input data must be 2-port data'
z0 = sub_full.z0[0, 0].real

append_log('Wideband inductor model fit from S2P S-parameter file')
append_log(f'Input file: {args.s2p}')
append_log(f'S2P frequency range is {sub_full.frequency.start/1e9} to {sub_full.frequency.stop/1e9} GHz')

# Drop only a DC point, where the oxide capacitance makes the shunt branches singular.
# Unlike the narrowband tools, all other low frequency data is kept, because it determines
# Rdc, Ls and Cox.
sub = sub_full[sub_full.f > 0] if sub_full.frequency.start == 0 else sub_full

f = sub.f
omega = 2 * math.pi * f
Ymeas = sub.y

Ysh1_meas, Ysh2_meas, Yser_meas = pi_branches(Ymeas)
Ydiff_meas = y_diff(Ymeas)
L11_meas, Q11_meas = l_and_q(Ymeas[:, 0, 0], omega)
L22_meas, Q22_meas = l_and_q(Ymeas[:, 1, 1], omega)
Ldiff_meas, Qdiff_meas = l_and_q(Ydiff_meas, omega)

# ---- measured self resonance and fit band ----
SRF1_meas = find_srf(f, Ymeas[:, 0, 0])
SRF2_meas = find_srf(f, Ymeas[:, 1, 1])
SRFdiff_meas = find_srf(f, Ydiff_meas)
srf_candidates = [x for x in (SRF1_meas, SRF2_meas) if x is not None]
SRF_meas = min(srf_candidates) if srf_candidates else None

fmin = args.fmin * 1e9 if args.fmin is not None else f[0]
if args.fmax is not None:
    fmax = args.fmax * 1e9
elif SRF_meas is not None:
    fmax = min(SRF_FIT_FACTOR * SRF_meas, f[-1])
else:
    fmax = f[-1]

fit_mask = (f >= fmin) & (f <= fmax)
assert np.count_nonzero(fit_mask) >= 10, 'not enough data points in fit band'
f_fit = f[fit_mask]
omega_fit = omega[fit_mask]

# L and Q goals only well below SRF, where they are well defined
lq_mask = fit_mask.copy()
if SRF_meas is not None:
    lq_mask &= f <= LQ_SRF_FACTOR * SRF_meas
lqdiff_mask = fit_mask.copy()
if SRFdiff_meas is not None:
    lqdiff_mask &= f <= LQ_SRF_FACTOR * SRFdiff_meas
# use these as masks into the fit band arrays
lq_fit = lq_mask[fit_mask]
lqdiff_fit = lqdiff_mask[fit_mask]

append_log(f'Fit band: {f_fit[0]/1e9:.3f} to {f_fit[-1]/1e9:.3f} GHz ({len(f_fit)} points)')


# ---------------- analytic seed values ----------------

# -- series branch: staged fit of R-L(+skin)||Cs to the series branch -Y12, similar to the
#    wideband series branch fit in pi_from_s2p, but on log parameters for robust scaling --
Zser_meas = 1.0 / Yser_meas
Rs0 = max(Zser_meas.real[0], 1e-3)
Ls0 = max(np.median(Zser_meas.imag[:3] / omega[:3]), 1e-12)

def y_series_branch(p_series, omega):
    # Rs, Ls, Rskin1, Lskin1, Rskin2, Lskin2, Cs
    return 1.0 / z_coil(*p_series[:6], omega) + 1j * omega * p_series[6]

def series_residuals_rlc(logp, omega, y_data):
    R, L, C = 10 ** logp
    model = y_series_branch([R, L, 1.0, 1e-20, 1.0, 1e-20, C], omega)  # skin sections shorted
    return complex_residual(model, y_data, floored_norm(y_data), 1.0)

def series_residuals_skin(x_free, x_template, omega, y_data):
    x = x_template.copy()
    x[series_free] = x_free
    model = y_series_branch(10 ** x, omega)
    return complex_residual(model, y_data, floored_norm(y_data), 1.0)

Yser_fit = Yser_meas[fit_mask]
Cs0 = 1e-15
res = least_squares(series_residuals_rlc, x0=np.log10([Rs0, Ls0, Cs0]), args=(omega_fit, Yser_fit),
                    bounds=(np.log10([1e-4, 1e-13, 1e-19]), np.log10([1e4, 1e-6, 1e-11])))
Rs1, Ls1, Cs1 = 10 ** res.x

# The skin effect corner frequencies are not known in advance - try (pairs of) seed corners
# spread over the fit band and keep the best. Lskin is capped at 1uH, so a corner below the
# band stays bounded.
corner_seeds = np.geomspace(max(f_fit[0], f_fit[-1] / 1000), f_fit[-1] / 2, 5)
if num_skin_sections == 2:
    corner_pairs = [(fc1, fc2) for i, fc1 in enumerate(corner_seeds) for fc2 in corner_seeds[i + 1:]]
else:
    corner_pairs = [(fc1, None) for fc1 in corner_seeds]
series_lower = np.log10([1e-4, 1e-13, 1e-4, 1e-16, 1e-4, 1e-16, 1e-19])[series_free]
series_upper = np.log10([1e4, 1e-6, 1e5, 1e-6, 1e5, 1e-6, 1e-11])[series_free]
best = None
for f_corner1, f_corner2 in corner_pairs:
    Rskin0 = max(Rs1, 0.1) / num_skin_sections
    skin2 = (Rskin0, Rskin0 / (2 * math.pi * f_corner2)) if f_corner2 else 10 ** np.array(X_SKIN2_FIXED)
    x_template = np.log10([Rs1, Ls1, Rskin0, Rskin0 / (2 * math.pi * f_corner1), *skin2, Cs1])
    res = least_squares(series_residuals_skin, x0=x_template[series_free],
                        args=(x_template, omega_fit, Yser_fit), bounds=(series_lower, series_upper))
    if best is None or res.cost < best.cost:
        best = res
        x_series = x_template.copy()
        x_series[series_free] = res.x
series_seed = 10 ** x_series

# -- shunt branches: Cox from the low frequency capacitance, Csi and Rsi from the high frequency
#    asymptotes Cp -> Cox*Csi/(Cox+Csi) and Rp -> Rsi*(Cox+Csi)^2/Cox^2, then refined by a
#    small least squares fit of Cox-(Rsi||Csi) to the measured shunt branch --
def shunt_residuals(logp, omega, y_data):
    model = y_shunt_branch(*(10 ** logp), omega)
    return np.concatenate([
        complex_residual(model, y_data, floored_norm(y_data), 1.0),
        real_residual(model.real, y_data.real, floored_norm(y_data.real), W_SHUNT_RE)])

def seed_shunt(Ysh):
    Cp = Ysh.imag / omega
    Cp_valid = Cp[Cp > 0]
    Cox = np.median(Cp_valid[:3])
    Cp_band = Cp[fit_mask & (Cp > 0)]
    Cinf = np.min(Cp_band) if len(Cp_band) else Cox
    Csi = Cox * Cinf / (Cox - Cinf) if Cinf < 0.95 * Cox else 20 * Cox
    Rp_band = 1.0 / Ysh.real[fit_mask & (Ysh.real > 0)]
    Rsi = np.min(Rp_band) * Cox ** 2 / (Cox + Csi) ** 2 if len(Rp_band) else 100.0
    Rsi = min(max(Rsi, 1.0), 1e5)
    x0 = np.log10([Cox, Rsi, Csi])
    res = least_squares(shunt_residuals, x0=x0, args=(omega_fit, Ysh[fit_mask]),
                        bounds=(x0 - [1, 3, 3], x0 + [1, 3, 3]))
    return 10 ** res.x

shunt1_seed = seed_shunt(Ysh1_meas)
shunt2_seed = seed_shunt(Ysh2_meas)

# -- substrate coupling: start "weak" and let the global fit decide whether it is needed --
Rsub12_seed = 100 * (shunt1_seed[1] + shunt2_seed[1]) / 2
Csub12_seed = 0.01 * (shunt1_seed[2] + shunt2_seed[2]) / 2

seed = np.concatenate([series_seed, shunt1_seed, shunt2_seed, [Rsub12_seed, Csub12_seed]])
seed[PARAM_NAMES.index('Cs')] = max(seed[PARAM_NAMES.index('Cs')], 1e-16)


# ---------------- global fit of the complete 2-port model ----------------

def model_quantities(p, omega, segments):
    Y = y_model(p, omega, segments)
    Ysh1, Ysh2, Yser = pi_branches(Y)
    L11, Q11 = l_and_q(Y[:, 0, 0], omega)
    L22, Q22 = l_and_q(Y[:, 1, 1], omega)
    Ldiff, Qdiff = l_and_q(y_diff(Y), omega)
    return dict(Y11=Y[:, 0, 0], Y22=Y[:, 1, 1], Yser=Yser, Ysh1=Ysh1, Ysh2=Ysh2,
                ReYsh1=Ysh1.real, ReYsh2=Ysh2.real,
                L11=L11, L22=L22, Ldiff=Ldiff, Q11=Q11, Q22=Q22, Qdiff=Qdiff, Y=Y)

# measured target data in the fit band, normalizations are fixed in advance from these
meas_fit = dict(Y11=Ymeas[fit_mask, 0, 0], Y22=Ymeas[fit_mask, 1, 1], Yser=Yser_meas[fit_mask],
                Ysh1=Ysh1_meas[fit_mask], Ysh2=Ysh2_meas[fit_mask],
                ReYsh1=Ysh1_meas[fit_mask].real, ReYsh2=Ysh2_meas[fit_mask].real,
                L11=L11_meas[fit_mask], L22=L22_meas[fit_mask], Ldiff=Ldiff_meas[fit_mask],
                Q11=Q11_meas[fit_mask], Q22=Q22_meas[fit_mask], Qdiff=Qdiff_meas[fit_mask])

all_fit = np.ones(len(f_fit), dtype=bool)
# fit goals: (quantity, weight, mask into fit band); complex quantities fit real and imaginary part
goals = [('Y11', W_Y, all_fit), ('Y22', W_Y, all_fit), ('Yser', W_Y, all_fit),
         ('Ysh1', W_SHUNT, all_fit), ('Ysh2', W_SHUNT, all_fit),
         ('ReYsh1', W_SHUNT_RE, all_fit), ('ReYsh2', W_SHUNT_RE, all_fit),
         ('L11', W_L, lq_fit), ('L22', W_L, lq_fit), ('Ldiff', W_L, lqdiff_fit),
         ('Q11', W_Q, lq_fit), ('Q22', W_Q, lq_fit), ('Qdiff', W_Q, lqdiff_fit)]
# skip L/Q goals with too few points (SRF close to the lower end of the fit band)
goals = [(key, w, mask, floored_norm(meas_fit[key][mask])) for key, w, mask in goals
         if np.count_nonzero(mask) >= 3]

def residuals(x, segments):
    model = model_quantities(10 ** x, omega_fit, segments)
    res = []
    for key, w, mask, norm in goals:
        if np.iscomplexobj(meas_fit[key]):
            res.append(complex_residual(model[key][mask], meas_fit[key][mask], norm, w))
        else:
            res.append(real_residual(model[key][mask], meas_fit[key][mask], norm, w))
    return np.concatenate(res)

# bounds in decades around the seed values; the skin and coupling elements get a wider range,
# Cox is well determined by the low frequency shunt capacitance and gets a narrow range
span = np.array([2, 2, 3, 3, 3, 3, 3,  COX_SPAN, 3, 3, COX_SPAN, 3, 3,  4, 4])
x_base_seed = np.log10(seed)
PORT2_TIED = {'Cox2': 'Cox1', 'Rsi2': 'Rsi1', 'Csi2': 'Csi1'}
X_COUPLING_OFF = {'Rsub12': 15.0, 'Csub12': -30.0}  # log10: 1e15 Ohm (open), 1e-30 F

def make_variant(name):
    # Free parameters, tied parameters (port 2 = port 1 for symmetric) and fixed parameters
    # (skin section 2 in the basic model, coupling elements if no coupling), with seed and bounds.
    symmetric, coupling = SUBSTRATE_VARIANTS[name]
    x_seed = x_base_seed.copy()
    fixed = {}
    if num_skin_sections == 1:
        fixed.update(Rskin2=X_SKIN2_FIXED[0], Lskin2=X_SKIN2_FIXED[1])
    if not coupling:
        fixed.update(X_COUPLING_OFF)
    ties = PORT2_TIED if symmetric else {}
    for dst, src in ties.items():  # symmetric seed: geometric mean of the port 1 and port 2 seeds
        x_seed[[PARAM_NAMES.index(src), PARAM_NAMES.index(dst)]] = (
            x_base_seed[PARAM_NAMES.index(src)] + x_base_seed[PARAM_NAMES.index(dst)]) / 2
    lower, upper = x_seed - span, x_seed + span
    for lskin in ('Lskin1', 'Lskin2'):
        upper[PARAM_NAMES.index(lskin)] = min(upper[PARAM_NAMES.index(lskin)], -6)  # Lskin <= 1uH
    for fixed_name, value in fixed.items():
        x_seed[PARAM_NAMES.index(fixed_name)] = value
    free = np.array([n not in fixed and n not in ties for n in PARAM_NAMES])
    return dict(name=name, symmetric=symmetric, coupling=coupling, fixed=fixed, ties=ties,
                free=free, x_seed=x_seed, lower=lower, upper=upper)

def expand(variant, x_free):
    # full log10 parameter vector from the free parameters of a variant
    x = variant['x_seed'].copy()
    x[variant['free']] = x_free
    for dst, src in variant['ties'].items():
        x[PARAM_NAMES.index(dst)] = x[PARAM_NAMES.index(src)]
    return x

def fit_variant(variant, segments, extra_starts):
    # best fit over all starts, returns (cost, full log10 parameter vector);
    # extra_starts are full parameter vectors from previous fits
    free = variant['free']
    lo, up = variant['lower'][free], variant['upper'][free]
    clip = lambda x: np.clip(x, lo + 1e-6, up - 1e-6)
    rng = np.random.default_rng(0)  # fixed seed: results are reproducible
    x0_seed = clip(variant['x_seed'][free])
    starts = [x0_seed] + [clip(x[free]) for x in extra_starts] + \
             [clip(x0_seed + rng.normal(0, 0.3, len(x0_seed))) for _ in range(NUM_RANDOM_STARTS)]
    best = None
    for x0 in starts:
        res = least_squares(lambda x_free: residuals(expand(variant, x_free), segments), x0=x0,
                            bounds=(lo, up), method='trf', x_scale=1.0, max_nfev=MAX_NFEV, ftol=1e-6)
        if best is None or res.cost < best.cost:
            best = res
    return best.cost, expand(variant, best.x)

def unphysical_reasons(variant, x):
    # Physical plausibility of the substrate network: Cox must not run away from the measured
    # low frequency shunt capacitance, and the substrate coupling must be a secondary path,
    # weaker than the path to ground of each port's own substrate network (otherwise one port
    # reaches ground only through the other port's network - a "shared substrate" solution).
    p = dict(zip(PARAM_NAMES, 10 ** x))
    reasons = []
    for name in ('Cox1', 'Cox2'):
        i = PARAM_NAMES.index(name)
        if variant['free'][i] and min(x[i] - variant['lower'][i], variant['upper'][i] - x[i]) < 0.005:
            reasons.append(f'{name} at bound')
    if variant['coupling']:
        ysub = np.abs(1 / p['Rsub12'] + 1j * omega_fit * p['Csub12'])
        for port in ('1', '2'):
            ysi = np.abs(1 / p['Rsi' + port] + 1j * omega_fit * p['Csi' + port])
            if np.any(ysub > ysi):
                reasons.append(f'coupling dominates port {port} substrate path')
    return reasons

# Symmetry is decided from the data, not from the fit cost: the asymmetric variant always fits
# a little better, by using its extra freedom on weakly determined substrate elements.
shunt_asymmetry = math.sqrt(np.mean(np.abs(Ysh1_meas[fit_mask] - Ysh2_meas[fit_mask]) ** 2) /
                            np.mean(np.abs((Ysh1_meas[fit_mask] + Ysh2_meas[fit_mask]) / 2) ** 2))
data_symmetric = shunt_asymmetry < SYMMETRY_TOLERANCE
if args.substrate != 'auto':
    substrate_choices = [args.substrate]
else:
    sym = 'symmetric' if data_symmetric else 'asymmetric'
    substrate_choices = [sym] if args.basic else [sym, sym + '+coupling']
variants = {name: make_variant(name) for name in substrate_choices}

# Step 1: number of coil segments, fitted with the substrate network without coupling (unique,
# no degenerate solutions). Parameters are totals, so the result for fewer segments is a good
# additional starting point for more segments.
segment_variant = substrate_choices[0]
segment_choices = [args.segments] if args.segments is not None else SEGMENT_CHOICES
segment_results = {}
for n in segment_choices:
    extra_starts = [segment_results[max(segment_results)][1]] if segment_results else []
    segment_results[n] = fit_variant(variants[segment_variant], n, extra_starts)
min_cost = min(cost for cost, _ in segment_results.values())
segments = min(n for n, (cost, _) in segment_results.items() if cost <= SEGMENT_COST_TOLERANCE * min_cost)

# Step 2: substrate network variant for this number of segments - the simplest physically
# meaningful variant with cost within tolerance of the best physical one
substrate_results = {segment_variant: segment_results[segments]}
for name in substrate_choices:
    if name not in substrate_results:
        substrate_results[name] = fit_variant(variants[name], segments, [segment_results[segments][1]])
substrate_reasons = {name: unphysical_reasons(variants[name], x) for name, (_, x) in substrate_results.items()}
physical = [name for name in substrate_choices if not substrate_reasons[name]]
candidates = physical if physical else substrate_choices
min_cost = min(substrate_results[name][0] for name in candidates)
substrate = next(name for name in candidates if substrate_results[name][0] <= SUBSTRATE_COST_TOLERANCE * min_cost)
variant = variants[substrate]
best_cost, x_fit = substrate_results[substrate]
lower, upper = variant['lower'].copy(), variant['upper'].copy()
x_seed = variant['x_seed'].copy()
# names of the elements in the selected model (fixed = removed, tied = shown with port 1 value)
active_names = [name for name in PARAM_NAMES if name not in variant['fixed']]

# sort the two skin sections by corner frequency (section 1 = lower corner), together with bounds
i1, i2 = PARAM_NAMES.index('Rskin1'), PARAM_NAMES.index('Rskin2')
corner = lambda x, i: x[i] - x[i + 1]  # log10(Rskin/Lskin), monotonic in corner frequency
if num_skin_sections == 2 and corner(x_fit, i1) > corner(x_fit, i2):
    for arr in (x_fit, lower, upper, x_seed):
        arr[[i1, i1 + 1, i2, i2 + 1]] = arr[[i2, i2 + 1, i1, i1 + 1]]

p_fit = 10 ** x_fit
fit = dict(zip(PARAM_NAMES, p_fit))
seed_values = dict(zip(PARAM_NAMES, 10 ** x_seed))


# ---------------- post processing ----------------

model = model_quantities(p_fit, omega, segments)
Y_model, Ysh1_model, Ysh2_model, Yser_model = model['Y'], model['Ysh1'], model['Ysh2'], model['Yser']
L11_model, L22_model, Ldiff_model = model['L11'], model['L22'], model['Ldiff']
Q11_model, Q22_model, Qdiff_model = model['Q11'], model['Q22'], model['Qdiff']

# model SRF on a dense grid, extended beyond the data range to allow for an extrapolated SRF
f_grid = np.linspace(f[0], 1.5 * f[-1], 20000)
Y_grid = y_model(p_fit, 2 * math.pi * f_grid, segments)
SRF1_model = find_srf(f_grid, Y_grid[:, 0, 0])
SRF2_model = find_srf(f_grid, Y_grid[:, 1, 1])
SRFdiff_model = find_srf(f_grid, y_diff(Y_grid))

def rms_rel_error(model, meas, mask):
    m, t = model[mask], meas[mask]
    return 100 * np.sqrt(np.mean(np.abs(m - t) ** 2) / np.mean(np.abs(t) ** 2))

def peak_q(Q, mask):
    i = np.argmax(np.where(mask, Q, -np.inf))
    return Q[i], f[i]

below_srf = f <= (SRF_meas if SRF_meas is not None else f[-1])
below_srfdiff = f <= (SRFdiff_meas if SRFdiff_meas is not None else f[-1])
Qpk11_meas, fQpk11_meas = peak_q(Q11_meas, below_srf)
Qpk11_model, fQpk11_model = peak_q(Q11_model, below_srf)
Qpkdiff_meas, fQpkdiff_meas = peak_q(Qdiff_meas, below_srfdiff)
Qpkdiff_model, fQpkdiff_model = peak_q(Qdiff_model, below_srfdiff)

f_skin_corner1 = fit['Rskin1'] / (2 * math.pi * fit['Lskin1'])
f_skin_corner2 = fit['Rskin2'] / (2 * math.pi * fit['Lskin2'])

def fmt_ghz(x):
    return f'{x/1e9:.3f}' if x is not None else 'n/a'

def fmt_value(name, value):
    if name[0] == 'R':
        return f'{value:12.4g} Ohm'
    if name[0] == 'L':
        return f'{value*1e9:12.4g} nH'
    return f'{value*1e15:12.4g} fF'

append_log('')
if args.basic:
    append_log('Basic model: 1 coil segment, 1 skin section, no substrate coupling (option --basic)')
if len(segment_results) > 1:
    append_log(f'Fit cost for number of coil segments ({segment_variant} substrate network): ' +
               ', '.join(f'{n}: {cost:.4g}' for n, (cost, _) in segment_results.items()))
append_log(f'Number of coil segments: {segments}' +
           (f' (fewest segments with cost within {SEGMENT_COST_TOLERANCE} x best)' if len(segment_results) > 1 else ''))
append_log(f'Port symmetry: shunt branches differ by {100*shunt_asymmetry:.2f} % RMS ' +
           ('(symmetric)' if data_symmetric else '(asymmetric)') +
           (f', threshold {100*SYMMETRY_TOLERANCE:.0f} %' if args.substrate == 'auto' else ''))
if len(substrate_results) > 1:
    append_log(f'Fit cost for substrate network variants:')
    for name in substrate_choices:
        reasons = substrate_reasons[name]
        append_log(f'  {name:20s}: {substrate_results[name][0]:.4g}' +
                   (f'   not physical: {", ".join(reasons)}' if reasons else ''))
append_log(f'Substrate network: {substrate}' +
           (f' (simplest physical variant with cost within {SUBSTRATE_COST_TOLERANCE} x best)' if len(substrate_results) > 1 and physical else ''))
if not physical:
    append_log('WARNING: no substrate network variant passed the physical plausibility check - '
               'substrate element values may not be physically meaningful')
append_log('')
append_log('Fitted model element values, totals over all segments (seed values in brackets)')
shunt2_title = 'Shunt branch port 2' + (' (symmetric: same as port 1)' if variant['symmetric'] else '')
for title, names in (('Series branch', PARAM_NAMES[:7]), ('Shunt branch port 1', PARAM_NAMES[7:10]),
                     (shunt2_title, PARAM_NAMES[10:13]), ('Substrate coupling port 1 to port 2', PARAM_NAMES[13:])):
    names = [name for name in names if name in active_names]
    if names:
        append_log(title)
    for name in names:
        append_log(f'  {name:7s}: {fmt_value(name, fit[name])}   ({fmt_value(name, seed_values[name]).strip()})')

append_log('')
append_log('Derived values')
if num_skin_sections == 2:
    append_log(f'  L at DC (Ls + Lskin1 + Lskin2)          [nH] : {(fit["Ls"] + fit["Lskin1"] + fit["Lskin2"])*1e9:.4g}')
    append_log(f'  R at high frequency (Rs+Rskin1+Rskin2) [Ohm]: {fit["Rs"] + fit["Rskin1"] + fit["Rskin2"]:.4g}')
    append_log(f'  Skin section 1 corner frequency        [GHz]: {f_skin_corner1/1e9:.4g}')
    append_log(f'  Skin section 2 corner frequency        [GHz]: {f_skin_corner2/1e9:.4g}')
else:
    append_log(f'  L at DC (Ls + Lskin1)           [nH] : {(fit["Ls"] + fit["Lskin1"])*1e9:.4g}')
    append_log(f'  R at high frequency (Rs+Rskin1) [Ohm]: {fit["Rs"] + fit["Rskin1"]:.4g}')
    append_log(f'  Skin section corner frequency  [GHz]: {f_skin_corner1/1e9:.4g}')

append_log('')
append_log('                                  measured     model')
append_log(f'  SRF port 1 (port 2 shorted) [GHz]: {fmt_ghz(SRF1_meas):>9s} {fmt_ghz(SRF1_model):>9s}')
append_log(f'  SRF port 2 (port 1 shorted) [GHz]: {fmt_ghz(SRF2_meas):>9s} {fmt_ghz(SRF2_model):>9s}')
append_log(f'  SRF differential            [GHz]: {fmt_ghz(SRFdiff_meas):>9s} {fmt_ghz(SRFdiff_model):>9s}')
append_log(f'  Peak Q11                         : {Qpk11_meas:9.2f} {Qpk11_model:9.2f}')
append_log(f'  Frequency of peak Q11       [GHz]: {fQpk11_meas/1e9:9.3f} {fQpk11_model/1e9:9.3f}')
append_log(f'  Peak Q differential              : {Qpkdiff_meas:9.2f} {Qpkdiff_model:9.2f}')
append_log(f'  Frequency of peak Q diff    [GHz]: {fQpkdiff_meas/1e9:9.3f} {fQpkdiff_model/1e9:9.3f}')
append_log(f'  L11 at {f[0]/1e9:.3f} GHz         [nH] : {L11_meas[0]*1e9:9.4f} {L11_model[0]*1e9:9.4f}')
append_log(f'  R11 at {f[0]/1e9:.3f} GHz        [Ohm] : {(1/Ymeas[0, 0, 0]).real:9.4f} {(1/Y_model[0, 0, 0]).real:9.4f}')

append_log('')
append_log(f'Fit cost (weighted error function): {best_cost:.4g}')
append_log(f'RMS relative error over fit band [%]: Y11 {rms_rel_error(Y_model[:, 0, 0], Ymeas[:, 0, 0], fit_mask):.2f}, '
           f'Y22 {rms_rel_error(Y_model[:, 1, 1], Ymeas[:, 1, 1], fit_mask):.2f}, '
           f'Y21 {rms_rel_error(Y_model[:, 1, 0], Ymeas[:, 1, 0], fit_mask):.2f}, '
           f'S21 {rms_rel_error(rf.network.y2s(Y_model, z0)[:, 1, 0], sub.s[:, 1, 0], fit_mask):.2f}')
append_log(f'RMS relative error below {LQ_SRF_FACTOR} x SRF [%]   : L11 {rms_rel_error(L11_model, L11_meas, lq_mask):.2f}, '
           f'Q11 {rms_rel_error(Q11_model, Q11_meas, lq_mask):.2f}, '
           f'L22 {rms_rel_error(L22_model, L22_meas, lq_mask):.2f}, '
           f'Q22 {rms_rel_error(Q22_model, Q22_meas, lq_mask):.2f}')

# ---- warnings ----
append_log('')
at_lower = x_fit - lower < 0.005  # within ~1% of the bound
at_upper = upper - x_fit < 0.005
# Elements that ran to the bound where they no longer have an effect (parallel C -> 0, parallel
# R -> open, series C -> short, skin section shorted) are simply not needed for this data. Any
# other element at a bound means the fit wanted to go further than the allowed range.
NEGLIGIBLE_AT = {'Cs': 'lower', 'Csi1': 'lower', 'Csi2': 'lower', 'Csub12': 'lower',
                 'Rsi1': 'upper', 'Rsi2': 'upper', 'Rsub12': 'upper',
                 'Rskin1': 'lower', 'Lskin1': 'lower', 'Rskin2': 'lower', 'Lskin2': 'lower'}
# If Rsi is negligible against Cox over the whole fit band, Cox effectively connects to ground
# (e.g. ground shield or low ohmic substrate contact) and Rsi/Csi are not determined by the data.
shorted_substrate = []
for port in ('1', '2'):
    if fit['Rsi' + port] * omega_fit[-1] * fit['Cox' + port] < 0.05 and not (variant['symmetric'] and port == '2'):
        shorted_substrate += ['Rsi' + port, 'Csi' + port]
        append_log(f'NOTE: Rsi{port} is negligible against Cox{port} in the fit band - Cox{port} effectively connects '
                   f'to ground, Rsi{port}/Csi{port} values are not determined by the data')
skin_negligible = []
for i, name in enumerate(PARAM_NAMES):
    if not variant['free'][i] or not (at_lower[i] or at_upper[i]) or name in shorted_substrate:
        continue
    bound = 'lower' if at_lower[i] else 'upper'
    if NEGLIGIBLE_AT.get(name) != bound:
        append_log(f'WARNING: {name} ended at its {bound} bound - check the fit, this element may not be determined by the data')
    elif name.startswith(('Rskin', 'Lskin')):
        skin_negligible.append(name[-1])
    else:
        append_log(f'NOTE: {name} ran to its {bound} bound - this element has no effect for this data')
for section in sorted(set(skin_negligible)):
    append_log(f'NOTE: skin section {section} is shorted (Rskin{section} or Lskin{section} at lower bound) - '
               f'it has no effect for this data')
for section, f_corner in (('1', f_skin_corner1), ('2', f_skin_corner2))[:num_skin_sections]:
    if section in skin_negligible:
        continue
    if f_corner < f_fit[0]:
        append_log(f'NOTE: skin section {section} corner frequency ({f_corner/1e9:.3g} GHz) is below the fit band - '
                   f'in band it acts as an additional series resistance Rskin{section}')
    elif f_corner > f_fit[-1]:
        append_log(f'NOTE: skin section {section} corner frequency ({f_corner/1e9:.3g} GHz) is above the fit band - '
                   f'in band it acts as an additional series inductance Lskin{section}')
if SRF_meas is None:
    append_log(f'WARNING: no self resonance found in the measured data - model SRF is extrapolated')
elif SRF1_model is None or SRF1_model > f[-1]:
    append_log(f'WARNING: model SRF is outside the data range - extrapolated')
if fmax < f[-1]:
    append_log(f'Data above {fmax/1e9:.3f} GHz was not used for fitting (beyond {SRF_FIT_FACTOR} x SRF)')
append_log('')


# ---------------- output files ----------------

base_filename = os.path.splitext(args.s2p)[0]

# SPICE netlist with explicit elements per coil segment. Coil nodes p1, c1 ... c<N-1>, p2 with
# substrate node s<k> under coil node k. Series elements per segment are total/N.
coil_node = lambda k: 'p1' if k == 0 else 'p2' if k == segments else f'c{k}'
netlist = [
    '* WIDEBAND RFIC INDUCTOR MODEL WITH SUBSTRATE NETWORK',
    f'* Created using inductor_fit.py from {os.path.basename(args.s2p)}',
    f'* Fit band {f_fit[0]/1e9:.3f} to {f_fit[-1]/1e9:.3f} GHz, {segments} coil segment(s)',
    '* Total element values:',
] + [f'*   {name} = {fit[name]:.6g}' for name in active_names] + [
    '*',
    '.SUBCKT inductor_model p1 p2',
]
for k in range(segments):
    sfx = f'_{k+1}' if segments > 1 else ''
    n_in, n_out = coil_node(k), coil_node(k + 1)
    netlist += [
        f'* coil segment {k+1}' if segments > 1 else '* series branch',
        f'Rs{sfx} {n_in} a{k+1} {fit["Rs"]/segments:.6g}',
        f'Ls{sfx} a{k+1} b{k+1} {fit["Ls"]/segments:.6g}',
    ]
    n_skin1 = f'd{k+1}' if num_skin_sections == 2 else n_out
    netlist += [
        f'Rskin1{sfx} b{k+1} {n_skin1} {fit["Rskin1"]/segments:.6g}',
        f'Lskin1{sfx} b{k+1} {n_skin1} {fit["Lskin1"]/segments:.6g}',
    ]
    if num_skin_sections == 2:
        netlist += [
            f'Rskin2{sfx} d{k+1} {n_out} {fit["Rskin2"]/segments:.6g}',
            f'Lskin2{sfx} d{k+1} {n_out} {fit["Lskin2"]/segments:.6g}',
        ]
netlist.append(f'Cs p1 p2 {fit["Cs"]:.6g}')
Cox_k, Rsi_k, Csi_k = segment_shunt_values(p_fit, segments)
for k in range(segments + 1):
    sfx = str(k + 1) if segments == 1 else f'_n{k}'
    netlist += [
        f'* substrate network at {coil_node(k)}',
        f'Cox{sfx} {coil_node(k)} s{k} {Cox_k[k]:.6g}',
        f'Rsi{sfx} s{k} 0 {Rsi_k[k]:.6g}',
        f'Csi{sfx} s{k} 0 {Csi_k[k]:.6g}',
    ]
if variant['coupling']:
    netlist += [
        '* substrate coupling',
        f'Rsub12 s0 s{segments} {fit["Rsub12"]:.6g}',
        f'Csub12 s0 s{segments} {fit["Csub12"]:.6g}',
    ]
netlist.append('.ENDS')
netlist_filename = base_filename + '_model.sp'
with open(netlist_filename, 'w', encoding='utf-8') as netlist_file:
    netlist_file.write('\n'.join(netlist) + '\n')
append_log(f'SPICE netlist written to {netlist_filename}')

# model S-parameters on the original frequency grid (including a DC point, if present)
model_filename = base_filename + '_model.s2p'
s_model_full = rf.network.y2s(y_model(p_fit, 2 * math.pi * sub_full.f, segments), z0)
with open(model_filename, 'w', encoding='utf-8') as snp_file:
    snp_file.write(f'! Inductor model fitted by inductor_fit.py from {os.path.basename(args.s2p)}\n')
    snp_file.write(f'#   Hz   S  RI   R   {z0:g}\n')
    for fi, s in zip(sub_full.f, s_model_full):
        values = [s[0, 0], s[1, 0], s[0, 1], s[1, 1]]  # Touchstone 2-port order: S11 S21 S12 S22
        snp_file.write(f'{fi:.6e} ' + ' '.join(f'{v.real:.9e} {v.imag:.9e}' for v in values) + '\n')
append_log(f'Model S-parameters written to {model_filename}')

log_text = "".join(log)

# log to terminal
print(log_text)

# output log message to file also, same basename as *.s2p input file but file extension .txt
log_filename = base_filename + '.txt'
with open(log_filename, "w", encoding="utf-8") as log_file:
    log_file.write(log_text)


# ---------------- plots ----------------

if args.noplot:
    raise SystemExit

fghz = f / 1e9

def mark_band(ax, srf):
    if srf is not None:
        ax.axvline(srf / 1e9, color='m', linestyle=':', label=f'SRF = {srf/1e9:.2f} GHz')
    if fmax < f[-1]:
        ax.axvspan(fmax / 1e9, f[-1] / 1e9, color='lightgray', alpha=0.4, label='not fitted')
    ax.set_xlabel('f (GHz)')
    ax.grid(color='lightgray')
    ax.legend()

def plot_pair(ax, meas, model, label, scale=1.0, srf=SRF_meas, ylim=None):
    ax.plot(fghz, meas * scale, 'b', label=f'{label} measured')
    ax.plot(fghz, model * scale, 'r--', label=f'{label} model')
    if ylim is not None:
        ax.set_ylim(*ylim)
    mark_band(ax, srf)

def l_ylim(L):
    L0 = L[0] * 1e9
    return (min(0, L0), 3 * abs(L0))

def q_ylim(Q, mask):
    Qmax = np.max(Q[mask])
    return (-0.5 * Qmax, 1.3 * Qmax)

title = f'{os.path.basename(args.s2p)}, {segments} coil segment(s)'

fig, axs = plt.subplots(2, 2, figsize=(12.8, 9.6))
fig.suptitle(f'{title}: inductance and Q factor')
plot_pair(axs[0, 0], L11_meas, L11_model, 'L11 [nH]', 1e9, ylim=l_ylim(L11_meas))
plot_pair(axs[0, 1], Q11_meas, Q11_model, 'Q11', ylim=q_ylim(Q11_meas, below_srf))
plot_pair(axs[1, 0], Ldiff_meas, Ldiff_model, 'Ldiff [nH]', 1e9, srf=SRFdiff_meas, ylim=l_ylim(Ldiff_meas))
plot_pair(axs[1, 1], Qdiff_meas, Qdiff_model, 'Qdiff', srf=SRFdiff_meas, ylim=q_ylim(Qdiff_meas, below_srfdiff))
fig.tight_layout()

fig, axs = plt.subplots(2, 2, figsize=(12.8, 9.6))
fig.suptitle(f'{title}: pi model branches')
for ysh_m, ysh_f, port in ((Ysh1_meas, Ysh1_model, 1), (Ysh2_meas, Ysh2_model, 2)):
    axs[0, 0].plot(fghz, ysh_m.imag / omega * 1e15, label=f'Cshunt{port} [fF] measured')
    axs[0, 0].plot(fghz, ysh_f.imag / omega * 1e15, '--', label=f'Cshunt{port} [fF] model')
    axs[0, 1].plot(fghz, 1 / ysh_m.real, label=f'Rshunt{port} [Ohm] measured')
    axs[0, 1].plot(fghz, 1 / ysh_f.real, '--', label=f'Rshunt{port} [Ohm] model')
Cp_max = 1.5 * max(np.max(Ysh1_meas.imag[fit_mask] / omega_fit), np.max(Ysh2_meas.imag[fit_mask] / omega_fit))
axs[0, 0].set_ylim(0, Cp_max * 1e15)
mark_band(axs[0, 0], None)
axs[0, 1].set_ylim(0, max(3 * np.median(1 / Ysh1_meas.real[fit_mask]), 1.0))
mark_band(axs[0, 1], None)
Zser_model = 1 / Yser_model
plot_pair(axs[1, 0], Zser_meas.imag / omega, Zser_model.imag / omega, 'Lseries [nH]', 1e9,
          srf=None, ylim=l_ylim(Zser_meas.imag / omega))
plot_pair(axs[1, 1], Zser_meas.real, Zser_model.real, 'Rseries [Ohm]', srf=None,
          ylim=(0, 3 * np.max(Zser_meas.real[lq_mask])))
fig.tight_layout()

fig, axs = plt.subplots(2, 2, figsize=(12.8, 9.6))
fig.suptitle(f'{title}: S-parameters')
s_model = rf.network.y2s(Y_model, z0)
for col, (i, j) in enumerate(((0, 0), (1, 0))):
    plot_pair(axs[0, col], 20 * np.log10(np.abs(sub.s[:, i, j])), 20 * np.log10(np.abs(s_model[:, i, j])),
              f'S{i+1}{j+1} [dB]')
    plot_pair(axs[1, col], np.angle(sub.s[:, i, j], deg=True), np.angle(s_model[:, i, j], deg=True),
              f'S{i+1}{j+1} [deg]')
fig.tight_layout()

plt.show()
