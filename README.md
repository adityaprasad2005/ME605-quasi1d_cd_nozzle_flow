# CFD Simulation of Quasi-1D CD Nozzle Flow using Conserved form of Euler equations

This repository contains the source code and results for a project simulating quasi-1D compressible flow through a convergent-divergent (CD) nozzle using the MacCormack predictor-corrector scheme. The discretized forms of the conserved Euler equations are employed to model the flow behavior.

## Project Objectives

* Implement the MacCormack predictor-corrector technique for solving the conserved variables in a quasi-1D compressible flow simulation.
* Analyze the flow characteristics within a CD nozzle for subsonic and sub-supersonic cases.
* Generate pressure, temperature, and Mach number profiles along the nozzle length for both flow regimes.

## Files

* `subsonic_flow.m`: MATLAB script for simulating purely subsonic flow in the CD nozzle.
* `subsonic_supersonic_flow.m`: MATLAB script for simulating the sub-supersonic flow case.
