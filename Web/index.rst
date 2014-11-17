
*******************
Welcome to CoolProp
*******************

These pages help you to get started using CoolProp and provide detailed information for the 
returning user. Please feel free to browse the pages and use the menu on the left to navigate
on this website.


What is CoolProp?
-----------------

CoolProp is a C++ library that implements:

- :ref:`Pure and pseudo-pure fluid equations of state and transport properties for 114 components <list_of_fluids>`
- :ref:`Mixture properties using high-accuracy Helmholtz energy formulations (or cubic EOS) <mixtures>`
- :ref:`Correlations of properties of incompressible fluids and brines <Pure>`
- :ref:`Highest accuracy psychrometric routines <Humid-Air>`


Environments Supported
----------------------

Programming Languages:

- Fully-featured wrappers: :ref:`Python (2.x, 3.x) <Python>` , :sfdownloads:`Modelica <Modelica>`, :ref:`Octave <Octave>`, :ref:`C# <Csharp>`, :ref:`MathCAD <MathCAD>`, :ref:`Java <Java>`, :ref:`MATLAB <MATLAB>`
- High-level interface only: :ref:`Labview <Labview>`, :ref:`EES <EES>`, :ref:`Microsoft Excel <Excel>`, :ref:`LibreOffice <LibreOffice>`, :ref:`Javascript <Javascript>`, :ref:`PHP <PHP>`, :ref:`FORTRAN <FORTRAN>`, :ref:`Maple <Maple>`, :ref:`Mathematica <Mathematica>`, :ref:`Scilab <Scilab>`, VB.net, :ref:`Delphi & Lazarus <Delphi>`

Architectures:

- 32-bit/64-bit
- Windows, Linux, OSX, Raspberry PI, VxWorks Compact Rio, etc. (if you can compile C++ on it, CoolProp will run)


High-Level Interface Example
----------------------------

In most languages, the code to calculate density ``D`` of Nitrogen at a temperature ``T`` of 298 K and a pressure ``P`` of 101325 Pa is something like::

    rho = PropsSI('D', 'T', 298.15, 'P', 101325, 'Nitrogen')
    
See more examples of PropsSI usage at :ref:`High-Level API <high_level_api>` or :ref:`Examples <examples>`
    

Help
----

- File a `Github issue <https://github.com/CoolProp/CoolProp/issues>`_
- Email the `Google group <https://groups.google.com/d/forum/coolprop-users>`_


Open-Source Projects Using CoolProp
-----------------------------------

- `Thermocycle <http://www.thermocycle.net/>`_
- `PDSim <http://pdsim.sourceforge.net/>`_
- `ACHP <http://achp.sourceforge.net/>`_
- `DWSim <http://sourceforge.net/projects/dwsim/>`_


Main Developers
---------------

The primary developers are:

- `Ian Bell <mailto:ian.h.bell@gmail.com>`_, `Sylvain Quoilin <mailto:squoilin@ulg.ac.be>`_, `Vincent Lemort <mailto:vincent.lemort@ulg.ac.be>`_, University of Liege, Liege, Belgium
- `Jorrit Wronski <mailto:jowr@mek.dtu.dk>`_, Technical University of Denmark, Kgs. Lyngby, Denmark


Supporters
----------

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
   :target: http://www.maplesoft.com
