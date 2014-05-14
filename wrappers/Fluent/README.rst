CoolProp wrapper for FLUENT
===========================

Contributors
------------
Primary CoolProp Developer: Ian Bell, University of Liege, Belgium (ian.h.bell@gmail.com)
FLUENT experts: Joris Degroote and Iva Papes, University of Gent, Belgium
First release: October 17, 2013

Requirements
------------
A linux version of FLUENT
g++

To Build
--------

Let us call the directory the directory where the Fluent wrapper is (coolprop/wrappers/Fluent) as CUSTOM_DIRECTORY.

1. Make sure you are in the CUSTOM_DIRECTORY

2. Make sure your case/mesh and UDF files are in the same directory as compile.sh (CUSTOM_DIRECTORY).

   a. If they are not, move your case/mesh and UDF files to CUSTOM_DIRECTORY. Do NOT move compile.sh.
   
   b. By default, an UDF example (UDF.c and properties.c) has been included in CUSTOM_DIRECTORY.
   
3. Run the script compile.sh (sh compile.sh SOLVER FLUENT_BIN_PATH), this should generate the libudf folder.

   a. SOLVER is the type of solver you want to run (2d, 2ddp, 3d, 3ddp)
   
   b. FLUENT_BIN_FOLDER is the path of your Fluent's installation bin folder (i.e. /home/ansys_inc/v145/fluent/bin)
   
   c. An example of how to run the shell file: sh compile.sh 2ddp /home/ansys_inc/v145/fluent/bin
   
   d. The script will compile all .c files in CUSTOM_DIRECTORY as Fluent UDFs and all .h files as UDF headers.
   
   e. Several warnings may show up, those should not be a problem.
   
4. Run Fluent.

   a. Make sure it runs the same solver as specified in the SOLVER variable (2d, 2ddp, 3d, 3ddp)
   
5. Open your case/mesh.

6. Compile the UDFs.

   a. Using the Fluent interface, go to Define > User-Defined > Functions > Compiled
   
   b. Select the .c files in CUSTOM_DIRECTORY as UDFs and .h files, if any, as headers
   
   c. Make sure "libudf" is written in the Library field
   
   d. Hit Load
   
7. The default UDF.c provided by the Coolprop wrapper is an EXECUTE_ON_DEMAND file.

   a. To check if it is working: Define > User-defined > Execute on Demand
   
   b. Select "call_coolprop::libudf" and hit execute
   
8. The default properties.c provided integrates thermal conductivity, density, viscosity and specific heat from Coolprop library with the Fluent solver. Do not load UDF.c and properties.c simultaneously.

   a. After loading properties.c, go to the Materials tab and change the model of each property listed above to user-defined and select the corresponding function (libudf::viscosity, libudf::density, libudf::specificHeat or libudf::thermalConductivity)
   
   b. The fluid is, by default, carbon dioxide. To change the working fluid, you have to change the FLUID variable in the properties.c UDF file BEFORE running compile.sh. A full list of fluids supported by Coolprop may be found here: http://coolprop.sourceforge.net/Fluids/FluidInformation.html
   
   c. If the operating pressure in your Fluent case is different than atmospheric (101325 Pa), you will also have to change gauge variable in properties.c
   
   c. Specific heat is currently only a function of temperature in the Fluent wrapper.

   
Note: If no argument is specified when running the shell file (step 3), then the script will assume Fluent can be run from command line (fluent) and the solver is 2d double precision (2ddp)
  
Warning
-------
Absolutely no guarantee of utility or accuracy can be made, although we have done our best to ensure useful and accurate results.  Caveat emptor!
