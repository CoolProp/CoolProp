This document discusses the `PropsSI` interface for MATLAB, based on the MATLAB⇒Python interface introduced in R2014b.

##Setup

Several things are required:
 1. A MATLAB release that is R2014b or newer.
 2. A python release with the [`numpy`](https://pypi.python.org/pypi/numpy) and [`CoolProp`](http://www.coolprop.org/coolprop/wrappers/Python/index.html) packages installed.
 It is not required to follow the ["supported python versions"](https://www.mathworks.com/help/matlab/matlab_external/system-requirements-for-matlab-engine-for-python.html#buijfe8) (this tool was tested with R2015a & R2017a + python 3.6).

If you have several python versions installed, make sure that your MATLAB recognizes the correct python executable by running `pyversion`.

To confirm that your system is configured correctly, please execute the provided tests file using `run(CoolPropsWrapperTests);` and verify that it completes with no errors.

##Usage
This tool is compatible with the 2- and 6-parameter syntax of `PropsSI`.

 - 2-input syntax (for ["trivial" inputs](http://www.coolprop.org/coolprop/HighLevelAPI.html#trivial-inputs "i.e. properties that do not depend on the thermodynamic state")):
  `PropsSI('Tcrit','Water')`. 
 - 6-input syntax ([standard](http://www.coolprop.org/coolprop/HighLevelAPI.html#sample-code "Many examples available...")): `PropsSI('D','T',298.15,'P',101325,'Air')`.

Whenever more than one output is requested, a 3D array will be returned.
  
##Troubleshooting

Please open an issue in the repository and tag @Dev-iL if you have problems.

###License
Copyright (C) 2017 Iliya Romm, under the MIT license.