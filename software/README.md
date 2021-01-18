# Subsurface Energy Storage Software
The goal of this software is to provide a simple simulation tool to model the flow of single and two phase synthetic fuels in a depleted oil or gas reservoir, mainly for the calculation of energy loss during injection and production of the energy storage medium.

# Dependencies
The PDE's are solved in FVTool and the physical propoerties are calculated using CoolProp. The whole package is written in Matlab, although it is possible to run things in Octave with some changes to the CoolProp call. You will also need a Python and CoolProp installation, and also need to introduce it to Matlab by
```
pyversion /home/ali/miniconda3/bin/python
```
It is possible to do the same in Windows and probably mac too.
Make sure that Matlab recognizes your Python environment by running `pyenv` command in Matlab. The output shoul NOT be an empty string. 

# Input files
TBD

## Input time series
It is possible to give a variable input power and the software will calculate the efficiency of fuel production, storage, and recuperation.