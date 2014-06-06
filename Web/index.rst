
=================
What is CoolProp?
=================

CoolProp is a C++ library that implements:

- Pure and pseudo-pure fluid equations of state and transport properties for 114 components
- Mixture properties using high-accuracy Helmholtz energy formulations (or cubic EOS)
- Correlations of properties of incompressible fluids and brines
- Highest accuracy psychrometric routines

======================
Environments Supported
======================

Programming Languages:

- Fully-featured wrappers: :sfdownloads:`Python (2.x, 3.x) <Python>` , :sfdownloads:`Modelica <Modelica>`, :sfdownloads:`Octave <Octave>`, :sfdownloads:`C# <Csharp>`, :sfdownloads:`MathCAD <MathCAD>`, :sfdownloads:`Java <Java>`
- High-level interface only: :sfdownloads:`Labview <Labview>`, :sfdownloads:`EES <EES>`, :sfdownloads:`Microsoft Excel <Microsoft Excel>`, :sfdownloads:`Javascript <Javascript>`, :sfdownloads:`MATLAB <MATLAB>`
- Rudimentary wrappers exist for FORTRAN, :sfdownloads:`Maple <Maple>`, :sfdownloads:`Mathematica <Mathematica>`, :sfdownloads:`Scilab <Scilab>`, VB.net

Architectures:

- 32-bit/64-bit
- Windows, Linux, OSX, Raspberry PI, VxWorks Compact Rio, etc.

============================
High-Level Interface Example
============================
In most languages, the code to calculate density ``D`` of Nitrogen at a temperature ``T`` of 298 K and a pressure ``P`` of 101325 Pa is something like::

    rho = PropsSI('D', 'T', 298.15, 'P', 101325, 'Nitrogen')
    
See more examples at `examples/examples <Examples>`_
    
===================================
Open-Source Projects Using CoolProp
===================================

- `Thermocycle <http://www.thermocycle.net/>`_
- `PDSim <http://pdsim.sourceforge.net/>`_
- `ACHP <http://achp.sourceforge.net/>`_
- `DWSim <http://sourceforge.net/projects/dwsim/>`_

====
Help
====

- File a `Github issue <https://github.com/ibell/coolprop/issues>`_
- Email the `Google group <https://groups.google.com/d/forum/coolprop-users>`_


===============
Main Developers
===============

The primary developers are:

- `Ian Bell <mailto:ian.h.bell@gmail.com>`_, `Sylvain Quoilin <mailto:squoilin@ulg.ac.be>`_, `Vincent Lemort <mailto:vincent.lemort@ulg.ac.be>`_, University of Liege, Liege, Belgium
- `Jorrit Wronski <mailto:jowr@mek.dtu.dk>`_, Technical University of Denmark, Kgs. Lyngby, Denmark

==========
Supporters
==========

.. image:: _static/labothap.png
   :height: 100px
   :alt: labothap
   :target: http://www.labothap.ulg.ac.be/
   
.. image:: _static/logo_ORCNext.jpg
   :height: 100px
   :alt: ORCNext
   :target: http://www.orcnext.be/

\

.. image:: _static/herrick.png
   :height: 100px
   :alt: Herrick
   :target: https://engineering.purdue.edu/Herrick/index.html
   
.. image:: _static/maplesoft_logo.png
   :height: 100px
   :alt: Maple
   :target: http://maplesoft.com/index.aspx