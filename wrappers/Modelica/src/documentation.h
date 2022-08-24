/*!
  \mainpage External Media HowTo

  \section intro_sec Introduction

  The ExternalMedia project was started in 2006 by Francesco Casella and
  Christoph Richter, with the aim of providing a framework for interfacing
  external codes computing fluid properties to Modelica.Media-compatible
  component models. The two main requirements are: maximizing the efficiency
  of the code and minimizing the amount of extra code required to use your
  own external code within the framework.

  The first implementation featured a hidden cache in the C++ layer and used
  integer unique IDs to reference that cache. This architecture worked well
  if the models did not contain implicit algebraic equations involving medium
  properties, but had serious issues when such equations were involved, which
  is often the case when solving steady-state initialization problems.

  The library has been restructured in 2012 by Francesco Casella and Roberto
  Bonifetto. The main idea has been to get rid of the hidden cache and of the
  unique ID references and use the Modelica state records for caching. In this
  way, all optimizations performed by Modelica tools are guaranteed to give
  correct results, which was previously not the case. The current
  version implements the Modelica.Media.Interfaces.PartialTwoPhaseMedium
  interface, so it can handle pure fluids, either one-phase or two-phase,
  and is compatible with Modelica and Modelica Standard Library 3.2. Please note
  that the paths of the medium packages have been changed from the previous
  versions, so you might need some small changes if you want to upgrade your
  models from previous versions of the ExternalMedia library.

  There are two ways to use this library. The easiest way is to use the
  releases available on the Modelica website, which include a pre-compiled
  interface to the FluidProp tool (http://www.fluidprop.com). FluidProp features
  many built-in fluid models, and can optionally be used to access the whole
  NIST RefProp database, thus giving easy access to a wide range of fluid models
  with state-of-the-art accuracy. If you want to use your own fluid property
  computation code instead, then you need to check out the source code and
  add the interface to it, as described in this manual.

  Please contact the main developer, Francesco Casella
  (<a href="mailto://casella@elet.polimi.it">casella@elet.polimi.it</a>)
  if you have questions or suggestions for improvement.

  Licensed by the Modelica Association under the Modelica License 2

  Copyright (c) 2006-2012, Politecnico di Milano, TU Braunschweig, Politecnico
  di Torino.

  \section releases_sec Using the pre-packaged releases with FluidProp

  Download and install the latest version of FluidProp from
  <a href="http://www.fluidprop.com">http://www.fluidprop.com</a>.
  If you want to use the RefProp fluid models,
  you need to get the full version of FluidProp, which has an extra license fee.

  Download and unzip the library corresponding to the version of Microsoft Visual Studio
  that you use to compile your Modelica models, in order to avoid linker errors.
  Make sure that you load the ExternalMedia library in your Modelica tool
  workspace, e.g. by opening the main package.mo file.

  You can now define medium models for the different libraries supported by FluidProp,
  by extending the ExternalMedia.Media.FluidPropMedium package. Please note that
  only single-component fluids are supported. Set libraryName
  to "FluidProp.RefProp", "FluidProp.StanMix", "FluidProp.TPSI", or "FluidProp.IF97",
  depending on the specific library you need to use.
  Set substanceNames to a single-element string array containing the name
  of the specific medium, as specified by the FluidProp documentation. Set
  mediumName to a string that describes the medium (this only used for
  documentation purposes but has no effect in selecting the medium model).
  See ExternalMedia.Examples for examples.

  Please note that the medium model IF97 is already available
  natively in Modelica.Media as Water.StandardWater, which is much faster than
  the FluidProp version. If you need ideal gas models (single-component or
  mixtures), use the medium packages contained in Modelica.Media.IdealGases.

  \section architecture_sec Architecture of the package

  This section gives an overview of the package structure, in order to help
  you understand how to interface your own code to Modelica using it.

  At the top level there is a Modelica package (ExternalMedia), which
  contains all the basic infrastructure needed to use external fluid
  properties computation software through a Modelica.Media compliant
  interface. In particular, the
  ExternalMedia.Media.ExternalTwoPhaseMedium package is a full-fledged
  implementation of a two-phase medium model, compliant with the
  Modelica.Media.Interfaces.PartialTwoPhaseMedium interface. The
  ExternalTwoPhaseMedium package can be used with any external fluid
  property computation software; the specific software to be used is
  specified by changing the libraryName package constant, which is then
  handled by the underlying C code to select the appropriate external
  code to use.

  The Modelica functions within ExternalTwoPhaseMedium communicate to a
  C/C++ interface layer (called externalmedialib.cpp) via external C functions
  calls, which in turn make use of C++ objects. This layer takes care of
  initializing the external fluid computation codes, called solvers from now on.
  Every solver is wrapped by a C++ class, inheriting from the BaseSolver C++
  class. The C/C++ layer maintains a set of active solvers, one for each
  different combination of the libraryName and mediumName strings, by means of
  the SolverMap C++ class. The key to each solver in the map is given by those
  strings. It is then possible to use multiple instances of many solvers in the
  same Modelica model at the same time.

  All the external C functions pass the libraryName, mediumName and substanceNames
  strings to the corresponding functions of the interface layer. These in turn
  use the SolverMap object to look for an active solver in the solver map,
  corresponding to those strings. If one is found, the corresponding function
  of the solver is called, otherwise a new solver object is instantiated and
  added to the map, before calling the corresponding function of the solver.

  The default implementation of an external medium model is implemented by the
  ExternalTwoPhaseMedium Modelica package. The setState_xx() and setSat_x()
  function calls are rerouted to the corresponding functions of the solver
  object. These compute all the required properties and return them in the
  ExternalThermodynamicState and  ExternalSaturationProperties C structs, which
  map onto the corresponding ThermodynamicState and SaturationProperties
  records defined in ExternalTwoPhaseMedium. All the functions returning
  properties as a function of the state records are implemented in Modelica and
  simply return the corresponding element in the state record, which acts as
  a cache. This is an efficient implementation for many complex fluid models,
  where most of the CPU time is spent solving the basic equation of state, while
  the computation of all derived properties adds a minor overhead, so it makes
  sense to compute them once and for all when the setState_XX() or setSat_xx()
  functions are called.

  In case some of the thermodynamic properties require a significant amount of
  CPU time on their own, it is possible to override this default implementation.
  On one hand, it is necessary to extend the ExternalTwoPhaseMedium Modelica
  package and redeclare those functions, so that they call the corresponding
  external C functions defined in externalmedium.cpp, instead of returning the
  value cached in the state record. On the other hand, it
  is also necessary to provide an implementation of the corresponding functions
  in the C++ solver object, by overriding the virtual functions of the
  BaseSolver object. In this case, the setState_xx() and setSat_X() functions
  need not compute all the values of the cache state records; uncomputed
  properties might be set to zero. This is not a problem, since Modelica.Media
  compatible models should never access the elements of the state records
  directly, but only through the appropriate functions, so these values should
  never be actually used by component models using the medium package.

  \section development_sec Developing your own external medium package

  The ExternalMedia package has been designed to ease your task, so that
  you will only have to write the mimum amount of code which is strictly
  specific to your external code - everything else is already provided.
  The following instructions apply if you want to develop an external
  medium model which include a (sub)set of the functions defined in
  Modelica.Media.Interfaces.PartialTwoPhaseMedium.

  The most straightforward implementation is the one in which all fluid
  properties are computed at once by the setState_XX() and setSat_X() functions
  and all the other functions return the values cached in the state records.

  Get the source code from the SVN repository of the Modelica Association:
  <a href="https://svn.modelica.org/projects/ExternalMediaLibrary/trunk">
  https://svn.modelica.org/projects/ExternalMediaLibrary/trunk</a>.

  First of all, you have to write you own solver object code: you can
  look at the code of the TestMedium and FluidPropMedium code as
  examples. Inherit from the BaseSolver object, which provides default
  implementations for most of the required functions, and then just add
  your own implementation for the following functions: object
  constructor, object destructor, setMediumConstants(), setSat_p(),
  setSat_T(), setState_ph(), setState_pT(), setState_ps(), setState_dT().
  Note that the setState and setSat functions need to compute and fill in all
  the fields of the corresponding C structs for the library to work correctly.
  On the other hand, you don't necessarily need to implement all of the four
  setState functions: if you know in advance that your models will only
  use certain combinations of variables as inputs (e.g. p, h), then
  you might omit implementing the setState and setSat functions corresponding
  to the other ones.

  Then you must modify the SolverMap::addSolver() function, so that it
  will instantiate your new solver when it is called with the appropriate
  libraryName string. You are free to invent you own syntax for the
  libraryName string, in case you'd like to be able to set up the
  external medium with some additional configuration data from within
  Modelica - it is up to you to decode that syntax within the addSolver()
  function, and within the constructor of your solver object. Look at how the
  FluidProp solver is implemented for an example.

  Finally, add the .cpp and .h files of the solver object to the C/C++ project,
  set the include.h file according to your needs and recompile it to a
  static library (or to a DLL). The compiled libraries and the externalmedialib.h
  files must then be copied into the Include subdirectory of the Modelica package
  so that the Modelica tool can link them when compiling the models.

  As already mentioned in the previous section, you might provide customized
  implementations where some of the properties are not computed by the setState
  and setSat functions and stored in the cache records, but rather computed
  on demand, based on a smaller set of thermodynamic properties computed by the
  setState and setSat functions and stored in the state C struct.

  Please note that compiling ExternalMedia from source code might require
  the professional version of Microsoft Visual Studio, which includes the
  COM libraries used by the FluidProp interface. However, if you remove
  all the FluidProp files and references from the project, then you should be
  able to compile the source code with the Express edition, or possibly also
  with gcc.
*/