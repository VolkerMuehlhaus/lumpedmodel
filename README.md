# Create lumped models from S-parameters

This is a collection of tools to convert S-parameters from RFIC EM simulation 
into a lumped model that can be used in all sorts of circuit simulation. 

The code for all model extractions requires Python3 with the skitkit-rf library.

## Inductor (two port, no center tap) from S2P data
[pi_from_s2p](./pi_from_s2p)

[<img src="./doc/inductor_model.png" width="500" />](./doc/inductor_model.png)

## Inductor (two port) and center-tapped inductor (three port), wideband model with substrate network
[inductor_fit](./inductor_fit)

Fully automatic wideband fit up to and beyond self resonance, with oxide capacitance,
bulk silicon network and skin effect. Exports a SPICE subcircuit.

[<img src="./inductor_fit/doc/inductor_fit_model.png" width="500" />](./inductor_fit/doc/inductor_fit_model.png)

## MIM capacitor from S2P data
[mim_from_s2p](./mim_from_s2p)

[<img src="./doc/mim_model.png" width="500" />](./doc/mim_model.png)

## Transmission line from S2P data
[rlgc_from_s2p](./rlgc_from_s2p)

[<img src="./doc/rlgc_segments.png" width="500" />](./doc/rlgc_segments.png)


