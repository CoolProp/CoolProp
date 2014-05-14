A wrapper of CoolProp for the Octave programming language (an open-source version of MATLAB)

This CoolProp.oct module is built using SWIG, and is built for octave for windows built with Visual Studio.  

Installation
============
Separate versions of the .oct file are available for Octave 3.6.1 and 3.6.2.
Put the .oct file for the version you have in somewhere in the octave path, or use the ``addpath`` function to add the folder that contains CoolProp.oct to the Octave path

On Linux systems you can put .oct file in
"/usr/share/octave/?octave.version.number?/m" folder. You will need superuser
privileges to do this.

If you place .oct file somewhere outside octave path, you have to use
"addpath" function at begining of your code.

Example: adding the folder that contains CoolProp.oct file to the Octave path:
    "addpath('/home/?user_name?/Some_folder/CoolProp')"

Developer Notes:
===============

A. 3.6.1 needs to use VS 2008 to build; 3.6.2 needs to use VS 2010

B. In the win32 distro of 3.6.2, the hard-coded path in OCTAVE/bin/include/math.h around line 73 might need to be changed to 

/* Include VC++ original math.h */

#include <c:/Program Files (x86)/Microsoft Visual Studio 10.0/VC/include/math.h>

depending on where your VS is installed


Building on Linux (Ubuntu in this case)
---------------------------------------
1. You will need to run 
      sudo apt-get update
      sudo apt-get install octave liboctave-dev swig
   to install the necesary dependencies.  The install of octave might not be necessary but it cant hurt
2. Check out the full source for coolprop from github
      git clone https://github.com/ibell/coolprop
3. Change into the coolprop-code/wrappers/Octave folder
      cd coolprop/wrappers/Octave
4. Call
      octave _OctaveBuilder_Linux.m
5. Call
      octave sample_code.m
   to run the sample
   
Building on Raspberry PI
------------------------
1. You will need to run
      sudo aptitude update
      sudo aptitude install octave liboctave-dev swig subversion
2. Download all the sources from subversion using
      git clone https://github.com/ibell/coolprop
3. Change into the folder
      cd coolprop/wrappers/Octave
4. Run the build script
      octave _OctaveBuilder_Linux.m
5. Call 
      octave sample_code.m
    to run the sample
    

