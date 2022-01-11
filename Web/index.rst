
*******************
Welcome to CoolProp
*******************

These pages help you to get started using CoolProp and provide detailed information for the
returning user. Please feel free to browse the pages and use the menu on the left to navigate
on this website.

What is CoolProp?
-----------------

CoolProp is a C++ library that implements:

* :ref:`Pure and pseudo-pure fluid equations of state and transport properties for 122 components <list_of_fluids>`
* :ref:`Mixture properties using high-accuracy Helmholtz energy formulations <mixtures>`
* :ref:`Correlations of properties of incompressible fluids and brines <Pure>`
* :ref:`Computationally efficient tabular interpolation <tabular_interpolation>`
* :ref:`Highest accuracy psychrometric routines <Humid-Air>`
* :ref:`User-friendly interface around the full capabilities of NIST REFPROP <REFPROP>`
* :ref:`Fast IAPWS-IF97 (Industrial Formulation) for Water/Steam <IF97>`
* :ref:`Cubic equations of state (SRK, PR) <cubic_backend>`

Environments Supported
----------------------

Programming Languages:

* Fully-featured wrappers: :ref:`Python (2.x, 3.x) <Python>`, :ref:`C++ (as static library) <static_library>`, :ref:`C++ as shared library <shared_library>`, `Modelica <https://github.com/modelica/ExternalMedia>`_, :ref:`Octave <Octave>`, :ref:`C# <Csharp>`, :ref:`VB.net <VBdotNet>`, :ref:`MathCAD <MathCAD>`, :ref:`Java <Java>`, :ref:`Android <Android>`, :ref:`MATLAB <MATLAB>`
* High-level interface only: :ref:`Labview <Labview>`, :ref:`EES <EES>`, :ref:`Microsoft Excel <Excel>`, :ref:`LibreOffice <LibreOffice>`, :ref:`Javascript <Javascript>`, :ref:`PHP <PHP>`, :ref:`FORTRAN <FORTRAN>`, :ref:`Maple <Maple>`, :ref:`Mathematica <Mathematica>`, :ref:`Scilab <Scilab>`, :ref:`Delphi & Lazarus <Delphi>`, :ref:`Julia <Julia>`

Architectures:

* 32-bit/64-bit
* Windows, Linux, OSX, Raspberry PI, VxWorks Compact Rio, etc. (if you can compile C++ on it, CoolProp will run)


High-Level Interface Example
----------------------------

In most languages, the code to calculate density ``D`` of Nitrogen at a temperature ``T`` of 298 K and a pressure ``P`` of 101325 Pa is something like::

    rho = PropsSI('D', 'T', 298.15, 'P', 101325, 'Nitrogen')

See more examples of PropsSI usage at :ref:`High-Level interface <high_level_api>` or :ref:`Examples <examples>`

.. _help:

Help
----

* (**General Discussion**) Email the `Google group <https://groups.google.com/d/forum/coolprop-users>`_
* (**Bugs, feature requests**) File a `Github issue <https://github.com/CoolProp/CoolProp/issues>`_
* `Docs for v4 of CoolProp <http://www.coolprop.org/v4/>`_
* `Docs for development version of CoolProp <http://www.coolprop.org/dev/>`_

Projects Using CoolProp
-----------------------------------

* `Thermocycle <http://www.thermocycle.net/>`_
* `PDSim <http://pdsim.sourceforge.net/>`_
* `ACHP <http://achp.sourceforge.net/>`_
* `DWSim <http://sourceforge.net/projects/dwsim/>`_
* `StateCalc <https://itunes.apple.com/us/app/statecalc/id891848148?ls=1&mt=8>`_
* `SmoWeb <http://platform.sysmoltd.com>`_
* `T-Props <https://play.google.com/store/apps/details?id=com.innoversetech.tprops>`_
* `PropiedadesDeFluidos <http://jfc.us.es/propiedadesdefluidos/descripcion/>`_
* `CoolPropJavascriptDemo <https://github.com/dvd101x/CoolPropJavascriptDemo>`_

Main Developers
---------------

.. warning:: Please do not email the developers directly, see :ref:`Help` above for assistance (this way the correspondence is google-able)

The primary developers are:

- `Ian Bell <mailto:ian.h.bell@gmail.com>`_, Bell Thermal Consultants
- `Jorrit Wronski <mailto:jowr@ipu.dk>`_, `IPU Refrigeration and Energy Technology <https://www.ipu.dk/expertise/thermodynamics-energy-technology/>`_, Kgs. Lyngby, Denmark
- `Sylvain Quoilin <mailto:squoilin@ulg.ac.be>`_, `Vincent Lemort <mailto:vincent.lemort@ulg.ac.be>`_, Thermodynamics Laboratory, University of Liege, Liege, Belgium

Please be so kind and cite our work in your publication: :ref:`Citation information <citation>`.

Supporters
----------

\ 

.. image:: _static/logo_labothap.png
   :height: 100px
   :alt: labothap
   :target: http://www.labothap.ulg.ac.be/

.. image:: _static/logo_ORCNext.jpg
   :height: 100px
   :alt: ORCNext

\

.. image:: _static/logo_herrick.png
   :height: 100px
   :alt: Herrick
   :target: https://engineering.purdue.edu/Herrick/index.html

.. image:: _static/logo_maplesoft.png
   :height: 100px
   :alt: Maple
   :target: https://www.maplesoft.com

\

.. image:: _static/logo_dtu_mekanik.png
   :height: 50px
   :alt: DTU Mechanical Engineering - Section for Thermal Energy
   :target: https://www.mek.dtu.dk/english/Sections/TES

.. image:: _static/logo_ipu.png
   :height: 50px
   :alt: IPU Refrigeration and Energy Technology
   :target: https://www.ipu.dk
   

License Information
-------------------

CoolProp has flexible licensing terms and you can use it for commercial projects and academic work free of charge. Have a look at the actual `license <https://github.com/CoolProp/CoolProp/blob/master/LICENSE>`_, if you are in doubt. 
