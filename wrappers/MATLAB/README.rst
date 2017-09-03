This document discusses the ``PropsSI`` interface for MATLAB, based on
the MATLAB⇒Python interface introduced in R2014b.

Setup
=====

Several things are required: 

1. A MATLAB release that is R2014b or newer.
2. A python release with the |numpy|_ and |CoolProp|_ packages installed. If you have several python environments/versions available, make sure that your MATLAB recognizes the desired executable using |pyversion|_.

.. |numpy| replace:: ``numpy``
.. _numpy: https://pypi.python.org/pypi/numpy

.. |CoolProp| replace:: ``CoolProp``
.. _CoolProp: http://www.coolprop.org/coolprop/wrappers/Python/index.html

.. |pyversion| replace:: ``pyversion``
.. _pyversion: https://www.mathworks.com/help/matlab/ref/pyversion.html

It is not required to follow the MATLAB's `supported python
versions <https://www.mathworks.com/help/matlab/matlab_external/system-requirements-for-matlab-engine-for-python.html#buijfe8>`__ (this tool was tested with R2015a & R2017a + python 3.6). 

To confirm that your system is configured correctly, please execute the provided tests file using ``run(CoolPropWrapperTests);`` and verify that it completes with no errors.

Usage
=====

This tool is compatible with the 2- and 6-parameter syntax of ``PropsSI``.

-  2-input short syntax (for `"trivial"
   inputs <http://www.coolprop.org/coolprop/HighLevelAPI.html#trivial-inputs>`__):
   ``PropsSI('Tcrit','Water')``.
-  6-input `standard <http://www.coolprop.org/coolprop/HighLevelAPI.html#sample-code>`__ syntax:
   ``PropsSI('D','T',298.15,'P',101325,'Air')``.
-  The first input to ``PropsSI`` (representing the requested outputs) can be either a character vector or a cell array thereof, in either syntax. Whenever more than one output is requested, a 3D array will be returned.

The tool also provides a convenience method for calling the low-level ``AbstractState`` API::

    [abState, CoolProp] = AbstractState('HEOS', 'Water');

The returned ``abState`` object can then be used according to relevant `documentation <http://coolprop.sourceforge.net/coolprop/LowLevelAPI.html>`__.

Troubleshooting
===============

If you have problems/questions/bugs/ideas related to the MATLAB wrapper, please open an issue in the repository and tag @Dev-iL.

License
-------

Copyright (C) 2017 Iliya Romm, under the MIT license.