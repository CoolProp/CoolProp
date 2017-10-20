CoolProp wrapper for FLUENT
===========================

Contributors
------------
Primary CoolProp Developer: Ian Bell, University of Liege, Belgium (ian.h.bell@gmail.com)
FLUENT experts: Joris Degroote and Iva Papes, University of Gent, Belgium
Others Contributors : Frederic Sonnino, AREVA company (frederic.sonnino1@areva.com)
Second release: October , 2017

Requirements
------------
A linux version of FLUENT
g++
python 2.7
cmake >2.8

To Build
--------

Let us call the main directory where the Fluent case and the fluent wrapper is (coolprop/wrappers/Fluent) as CUSTOM_DIRECTORY.

1. Make sure you are in the CUSTOM_DIRECTORY

2. Make sure your case/mesh and UDF files are in the same directory as compile.sh (CUSTOM_DIRECTORY).

   a. If they are not, move your case/mesh and UDF files to CUSTOM_DIRECTORY. Do NOT move compile.sh.
   
   b. By default, several UDF examples (CoolProp_Properties_of_Water.c OR CoolProp_Properties_of_Brine.c) has been included in CUSTOM_DIRECTORY, one and only must be chosen
   
3. Run the script compile.sh (sh compile.sh SOLVER FLUENT_BIN_PATH), this should generate the libudf folder.

   a. SOLVER is the type of solver you want to run (2d, 2ddp, 3d, 3ddp)
   
   b. FLUENT_BIN_FOLDER is the path of your Fluent's installation bin folder (i.e. /home/ansys_inc/v150/fluent/bin)
   
   c. An example of how to run the shell file: sh compile.sh 2ddp /home/ansys_inc/v145/fluent/bin
   
   d. The script will compile all .c files in CUSTOM_DIRECTORY as Fluent UDFs and all .h files as UDF headers.
   
   e. Several warnings may show up, those should not be a problem.
   
   f. another script is present compile_only_udf.sh, it runs the same way (sh compile_only_udf.sh SOLVER FLUENT_BIN_PATH) but it skips the coolprop compilation (in order to save compilation time) and only makes the compilation of the fluent udf into shared objects and link it with shared library CoolProp and output libudf.so. It could be useful while debugging other udfs for instance
   
4. Add in your .bashrc the following line and replace CUSTOM_DIRECTORY by the complete path coolprop/wrappers/Fluent/CoolProp_Build 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:CoolProp_Build_directory 
   
5. Run Fluent.

   a. Make sure it runs the same solver as specified in the SOLVER variable (2d, 2ddp, 3d, 3ddp)
   
6. Open your case/mesh.

7. By default the udf examples ( for instance CoolProp_Properties_of_water.c ) embed EXECUTE_ON_DEMAND  

   a. Using the Fluent interface to check if it is working: Define > User-defined > Execute on Demand
      
   b. Select "call_coolprop_water::libudf" and hit execute
   
8. The default CoolProp_Properties_of_Water.c provided integrates thermal conductivity, density, viscosity and specific heat from CoolProp library with the Fluent solver.

   a. After loading , go to the Materials tab and change the model of each property listed above to user-defined and select the corresponding function (libudf::water_viscosity, libudf::water_density, libudf::water_specificHeat or libudf::water_thermalConductivity)
   
   b. The fluid is, by default, water. To change the working fluid, you have to change the FLUID[] variable in the CoolProp_Properties_of_Water.c UDF file BEFORE running compile.sh. A full list of fluids supported by Coolprop may be found here: http://www.coolprop.org/FluidInformation.html
   Another example is given in CoolProp_Properties_of_Brine.c, the FLUID[] is a mixture between chloride sodium with a 20% concentration
   
   c. If the operating pressure in your Fluent case is different than atmospheric (101325 Pa), you will also have to change gauge pressure variable in CoolProp_Properties_of_Water.c
   
   c. Specific heat is currently only a function of temperature in the Fluent wrapper.

   
Note: If no argument is specified when running the shell file (step 3), then the script will assume Fluent can be run from command line (fluent) and the solver is 3d double precision (2ddp) ---> actually it do not work


  
Warning
-------
Absolutely no guarantee of utility or accuracy can be made, although we have done our best to ensure useful and accurate results.  Caveat emptor!
