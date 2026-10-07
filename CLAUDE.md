# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

A collection of independent Python scripts that extract lumped circuit models from S-parameter
(Touchstone `.s2p`/`.snp`) data measured or simulated for RFIC passives. Each subfolder is a
standalone, single-script tool with its own README — there is no shared package, build system,
or test suite. Everything runs via `python <script>.py <args>`.

## Prerequisites

Python 3 with `scikit-rf` (`skrf`). `vector_fit/vectorfit_sparam.py` additionally needs
`packaging` (for a scikit-rf version check) and uses `numpy` directly; the other scripts only
need `skrf`, `matplotlib`, and stdlib. There is no `requirements.txt` — install manually:

```
pip install scikit-rf matplotlib numpy packaging
```

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
- Only `vector_fit/vectorfit_sparam.py` produces a genuinely wideband model; its accuracy should
  be checked by comparing the `_predicted` S2P output against the original input data.

## Input/output data files

Sample `.s2p`/`.snp` and generated `.txt`/`_predicted.s2p`/`.sp` files live alongside each
script and are checked into the repo as usage examples — treat existing samples as fixtures
when testing changes, and note that running any tool regenerates the `.txt` log and (for
`vector_fit`) `.sp`/`_predicted` files next to the input.
