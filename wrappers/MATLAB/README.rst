The *.mexw64 files in this folder are wrappers of CoolProp for 64-bit MATLAB on Windows
The *.mexw32 files in this folder are wrappers of CoolProp for 32-bit MATLAB on Windows
You can tell what type of MATLAB you have when you start up MATLAB.  It will say on the splash screen

Place them somewhere on the MATLAB path, or add their containing folder to the MATLAB path.

For Users
=========
There is an example file called MATLAB_sample.m that when run will demonstrate calling both
Props (for fluid properties) and HAProps(for humid air properties)

Due to the poor integration of MATLAB and C/C++ code, the wrapper for MATLAB is 
a mess compared with other languages, and does not follow the 
semantics of the other languages.

The "normal" calls of Props(Ref,Output) and Props(Output,Name1,Val1,Name2,Val2,Fluid)
work just like the other wrappers.  Here is a summary of the other hacky things
that have been added to the MATLAB wrapper:

For a given fluid, you can use one of the following types of calls to get the given output::

    Props(Fluid,'aliases') : Returns a comma separated string of aliases for a given fluid
    Props(Fluid,'CAS') : Returns CAS code for the given fluid
    Props(Fluid,'ASHRAE34') : Returns ASHRAE 34 safety code
    Props(Fluid,'REFPROPName') : Returns REFPROP name for fluid
    Props(Fluid,'TTSE_mode') : Returns TTSE mode, one of 'TTSE' or 'BICUBIC'
    Props(Fluid,'enable_TTSE') : Turn on tabular interpolation for a given fluid
    Props(Fluid,'disable_TTSE') : Turn on tabular interpolation for a given fluid

Or you can get a global constant using a call like::

    Props('FluidsList') : Returns a comma separated string of fluids in CoolProp
    Props('version') : Returns the version of Coolprop in Use
    Props('gitrevision') : Returns the git revision
    Props('set_UNIT_SYSTEM_SI') : Set the unit system to SI units
    Props('set_UNIT_SYSTEM_KSI') : Set the unit system to kilo-SI units
    Props('get_unit_system') : Get the unit system in use, one of 'UNIT_SYSTEM_SI' or 'UNIT_SYSTEM_KSI'

For Developers
==============
To build the mex files on windows, you should enter these commands at the command prompt::

    git clone https://github.com/ibell/coolprop coolprop-code
    cd coolprop-code/wrappers/MATLAB
    matlab -r MATLABBuilder.m

You will need to have git installed (google it).  You will also need a compiler installed (Visual studio express works). You might need to run mex -setup to select the compiler

On OSX/Linux, the same idea.  Do this::

    git clone https://github.com/ibell/coolprop coolprop-code
    cd coolprop-code/wrappers/MATLAB
    # Change line 6 MATLABBuilder_OSX.m according to instructions in file.
    matlab -r MATLABBuilder_OSX.m

Please send an email to ian.h.bell@gmail.com if it works or if you have problems

(c) Ian Bell, 2012-