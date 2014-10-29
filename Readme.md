# Turbulence Monte Carlo

Simulates propagation of second and fourth-order optical beams through atmospheric turbulence. Built upon routines available in scientific literature with added utility methods that provide feedback to the user and, more importantly, allows for the propagation of correlation beams (two-photon optics).

## Getting started
To run a simulation, write a JSON file with the simulation parameters and run
```matlab
runSimulation
```
in the MATLAB or Octave shell.

## Simulation parameters
Parameters for the simulations are provided via a JSON file. This allows the user to have different files, each for a specific simulation scenario. 

## Simulation types
Different simulation types have been implemented, each performing different operations and providing different outputs. Currently these are:
* PSIPARITY
* PSICHI2
* IRRADIANCE
* MODEMATCHING
* MODEPARITY
