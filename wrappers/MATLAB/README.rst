This document discusses the ``PropsSI`` interface for MATLAB, based on
the MATLAB⇒Python interface introduced in R2014b.

Setup
=====

Several things are required: 

1. A MATLAB release that is R2014b or newer. 
2. A python release with the |numpy|_ and |CoolProp|_ packages installed. It is not required to follow the `"supported python
versions" <https://www.mathworks.com/help/matlab/matlab_external/system-requirements-for-matlab-engine-for-python.html#buijfe8>`__
(this tool was tested with R2015a & R2017a + python 3.6).

.. |numpy| replace:: ``numpy``
.. _numpy: https://pypi.python.org/pypi/numpy

.. |CoolProp| replace:: ``CoolProp``
.. _CoolProp: http://www.coolprop.org/coolprop/wrappers/Python/index.html

If you have several python versions installed, make sure that your
MATLAB recognizes the correct python executable by running
``pyversion``.

To confirm that your system is configured correctly, please execute the
provided tests file using ``run(CoolPropsWrapperTests);`` and verify
that it completes with no errors.

Usage
=====

This tool is compatible with the 2- and 6-parameter syntax of
``PropsSI``.

-  2-input syntax (for `"trivial"
   inputs <http://www.coolprop.org/coolprop/HighLevelAPI.html#trivial-inputs>`__):
   ``PropsSI('Tcrit','Water')``.
-  6-input syntax
   (`standard <http://www.coolprop.org/coolprop/HighLevelAPI.html#sample-code>`__):
   ``PropsSI('D','T',298.15,'P',101325,'Air')``.

Whenever more than one output is requested, a 3D array will be returned.

Troubleshooting
===============

If you have problems please open an issue in the repository and tag @Dev-iL.

License
-------

Copyright (C) 2017 Iliya Romm, under the MIT license.
