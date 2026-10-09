# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

A collection of independent Python scripts that extract lumped circuit models from S-parameter
(Touchstone `.s2p`/`.snp`) data measured or simulated for RFIC passives. Each subfolder is a
standalone, single-script tool with its own README — there is no shared package, build system,
or test suite. Everything runs via `python <script>.py <args>`.

## Prerequisites

Python 3 with `scikit-rf` (`skrf`) and `matplotlib` (all tools), `numpy` (all except
`mim_from_s2p`), `scipy` (`inductor_fit`, `pi_from_s2p`) and `packaging` (`vector_fit`, for a
scikit-rf version check). `inductor_fit/doc/draw_model.py` also needs `schemdraw`. There is no
`requirements.txt` — install manually:

```
pip install scikit-rf numpy scipy matplotlib packaging
```

The root `README.md` has the requirements table with the dependency licenses. The repo is GPL-3.0-or-later;
all dependencies are permissive (BSD, MIT, Apache-2.0, Matplotlib/PSF-based) and only imported,
not bundled. When adding a dependency, check its license for GPL-3.0 compatibility and add it to
that table and to the tool's own README.

## Running the tools

Each tool is invoked directly with argparse-style positional arguments; run with `-h` for help.

- `pi_from_s2p/pi_from_s2p.py <s2p> <f_ghz>` — narrowband pi model for a 2-port RFIC inductor
  (series L/R + shunt C/R per port) at a single extraction frequency.
- `mim_from_s2p/mim_from_s2p.py <s2p> <f_ghz>` — wideband-ish MIM capacitor model (series L/R/C
  + shunt C split across ports), using both the target frequency and an auto-picked low
  frequency (`max(f_target/20, 1 GHz)`) to separate C from parasitic L.
- `rlgc_from_s2p/rlgc_from_s2p.py <s2p> <f_ghz> <l_um> <z0_ohm>` — per-meter RLGC transmission
  line parameters from a 2-port line measurement; physical length and port impedance must be
  supplied explicitly.
- `inductor_fit/inductor_fit.py <s2p|s3p> [--basic] [--segments N] [--substrate VARIANT] [--ct-port N] [--fmin/--fmax GHz] [--peak-q-weight [A]] [--shunt-accuracy [K]] [--spectre] [--noplot]` —
  fully automatic wideband RFIC inductor model with physical substrate network (series Rs/Ls +
  two Rskin||Lskin sections + Cs, per-port Cox-(Rsi||Csi), optional Rsub12||Csub12 coupling),
  global least-squares fit up to 1.2x SRF. Auto-selects 1-3 coil segments for distributed
  behaviour and the simplest physically meaningful substrate variant; `--basic` fixes 1 segment /
  1 skin section / no coupling (consistent, unique topology, e.g. for ML data tables).
  3-port input switches to the center-tapped inductor model (two coupled half coils, effective
  k, Rct lead, substrate at p1/p2/ct; `--ct-port` names the tap port).
  Writes `.txt`, `<base>_model.sp` (SPICE subckt) and `<base>_model.s2p`/`.s3p`; needs `scipy`.
  `doc/draw_model.py` regenerates both schematic images (needs `schemdraw`).
- `vector_fit/vectorfit_sparam.py <snp> [numpoles]` — full n-port wideband vector fit (via
  `skrf.vectorFitting.VectorFitting`), auto-determining pole count if `numpoles` is omitted;
  emits a SPICE subcircuit netlist and a `_predicted.s2p`/`.snp` file for validation.

Example:
```
python pi_from_s2p/pi_from_s2p.py pi_from_s2p/sample_inductor.s2p 20
```

## Common script structure

The three `*_from_s2p.py` extraction scripts (`pi_from_s2p`, `mim_from_s2p`, `rlgc_from_s2p`)
follow the same pattern — when modifying one, check whether the same fix applies to the others:

1. Parse args with `argparse` (input file + extraction frequency in GHz, plus tool-specific
   params like line length or port impedance).
2. Load the Touchstone file via `skrf.Network(...)`; if the data includes a DC point, resample
   to start at 1 GHz to avoid warnings/singularities (`pi_from_s2p`, `mim_from_s2p` only).
3. Convert to Z/Y parameters at each frequency point and derive the lumped-element values
   analytically (no optimizer — closed-form extraction from S/Y/Z at the target frequency).
4. Build up output via a local `log`/`append_log()` list, print it, and also write it to
   `<input_basename>.txt` next to the input file.
5. Plot results with matplotlib (`plt.show()`), marking the extraction frequency on the curves.

`inductor_fit/inductor_fit.py` is also structurally different: analytic seed values from the pi
(2-port: `Y11+Y12`, `Y22+Y12`, `-Y12`) or delta (3-port) decomposition, then a global
`least_squares` fit of the complete model (nodal matrix + Kron reduction: `y_model` for 2-port,
`y_model_ct` with `stamp_coupled` for the coupled half coils of the 3-port model) on log-scaled
parameters. Shared engine (`make_variant`/`expand`/`fit_variant`/`fit_segments`, seeding, bound
notes, plot helpers) at the top; `run_two_port()` and `run_center_tap()` hold the topology
specific flow, dispatched on the number of ports. Parameters are always totals;
`shunt_weights()` distributes the substrate network over coil segments. The SPICE netlist
writers must stay consistent with `y_model`/`y_model_ct` - verify by solving the written netlist
(including `K` elements) independently and comparing against `_model.s2p`/`.s3p`. `--spectre`
output (`.scs`) is not a separate writer: `spice_to_spectre()` translates the generated SPICE lines,
so a new element type in the SPICE writers also needs a mapping there. Refactoring
must keep 2-port results byte-identical (compare `.txt`/`.sp`/`.s2p` of the samples before/after).

`vector_fit/vectorfit_sparam.py` is structurally different: it fits a rational-function model
across the whole band using scikit-rf's `VectorFitting`, checks/enforces passivity, and writes
a SPICE subcircuit (`<basename>.sp`) plus a predicted S-parameter file
(`<basename>_predicted<ext>`) instead of doing per-frequency analytic extraction. It branches on
`skrf.__version__` (>= 1.11.0 vs older) because the `enforce_dc` kwarg to `vector_fit()` was
added in 1.11 and changes fit behavior — preserve this version branch when touching the fitting
calls.

## Model accuracy notes (relevant when changing extraction logic)

- The three analytic extraction tools are **exact only at the chosen extraction frequency**;
  they are not wideband models. `pi_from_s2p` in particular can show non-intuitive frequency
  trends (e.g. series R decreasing with frequency) that are correct for a single-frequency-exact
  model — don't "fix" this as if it were a bug.
- `mim_from_s2p` deliberately splits shunt capacitance equally across both ports rather than
  extracting it per-port, because the two ports are tightly coupled through the large series MIM
  value.
- Only `vector_fit/vectorfit_sparam.py` (black-box rational fit) and `inductor_fit/inductor_fit.py`
  (physical equivalent circuit) produce genuinely wideband models; check accuracy by comparing
  the `_predicted`/`_model` S2P output against the original input data.
- In `inductor_fit`, the substrate network is not always unique from 2-port data; the design
  rule is "physically meaningful first, symmetric if the data is, simpler when in doubt".
  Symmetry is decided from the data (shunt branch difference < 5%), not from fit cost (the
  asymmetric variant always fits slightly better by abusing weakly determined elements).
  Rsub12||Csub12 coupling is only kept if it passes `unphysical_reasons()` (Cox near the
  low-frequency shunt C; coupling weaker than each port's own path to ground - rejects "shared
  substrate" solutions) and improves cost by >10%. Small mm-wave inductors (e.g.
  `inductor_74p.s2p`) need 2-3 coil segments to follow the data.
- Center-tapped model: the center tap lead inductance and the half coil mutual inductance are
  exactly degenerate in terminal data (T-equivalent), so there is deliberately no `Lct` and `k`
  is an *effective* coupling - don't add `Lct` back as a free parameter. Symmetry is decided per
  group (half coil L, half coil R, p1/p2 substrate): real inductor3 data had symmetric L but 6%
  different half coil R, and a single all-or-nothing tie let the fit make L asymmetric
  (unphysical). No substrate coupling
  element in the 3-port model (simplest physical choice). Verified with synthetic S3P data
  (exact recovery) and Palace EM data of an IHP SG13G2 inductor3 (2 turns).
- `inductor_fit` low-frequency behaviour:
  - Skin section corners are softly constrained to the fit band (`corner_penalty`, plus the
    same clip on the series seed, because the global bounds are centered on the seeds).
  - Without that constraint, FDTD data gave corners at 0.1-1 MHz. The model then had a DC
    inductance of 330-985 nH instead of about 5 nH, and too low Rdc.
  - Don't add an option that drops or zero-weights low-frequency data. The series branch there
    is the only anchor for Rdc and L: tested with `--fmin` 1-2 GHz and weight ramps, Rdc
    collapsed to 0.01-0.5 Ohm.
  - The low-frequency weakness of measured and FDTD data is in the shunt branch (1e-4 to 1e-3 of
    the series branch). `--shunt-accuracy` addresses it by flooring the shunt goal
    normalization at K x |series|.
  - A truncated FDTD run (openEMS energy limit -40 dB) gives a smooth, plausible-looking R offset
    that no fit option can detect or fix.
- `--peak-q-weight` is a trade-off, so it is opt-in with default A = 2. A = 5 doubled the 3-port
  Y21 error on flat-Q data. Fit costs are only comparable between runs with the same weight
  options (segment and substrate selection compare costs within one run, which stays valid).
- `make_goals` weights can be scalars or per-point arrays. With the options off, the 2-port and
  3-port results are byte-identical to the code before the options were added; only the corner
  constraint changes results.

## Input/output data files

Sample `.s2p`/`.snp` and generated `.txt`/`_predicted.s2p`/`.sp` files live alongside each
script and are checked into the repo as usage examples — treat existing samples as fixtures
when testing changes, and note that running any tool regenerates the `.txt` log and (for
`vector_fit`) `.sp`/`_predicted` files, (for `inductor_fit`) `_model.sp`/`_model.s2p`/`.s3p` files
next to the input.
`inductor_fit` samples live in `inductor_fit/examples/` (inputs and generated outputs together).
