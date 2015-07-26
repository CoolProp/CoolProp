
Frequently Asked Questions
==========================

Usage
-----

1. My enthalpy and entropy values are not the same as what are used in REFPROP, or EES, or ...

    Enthalpy and entropy are only equivalent to within an additive constant, so you should always compare enthalpy and entropy *differences* rather than enthalpy and entropy values
    
2. My calls to REFPROP fail with a mysterious `Segementation Fault`: 

    Make sure that you have at least REFPROP version 9.1 installed. If you have a license for 9.0, 
    the upgrade is for free: http://www.nist.gov/srd/nist23.cfm
    
3. Transport properties
  1. Fluid XY does not have viscosity/thermal conductivity/ ... 

    Please have a look at 
    http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids 
    if there is no reference for your property, we do not have the information required to 
    implement the functionality. If you file an issue, please provide information and 
    a link to a publication with the required data. 
    
  1. ..but it used to work in an earlier version of CoolProp
    
    For some fluids, version 5 does not have transport properties even though version 4 had 
    that information. We decided to remove some correlations, which had very large uncertainties 
    and we would like to make the user aware of the fact that the old implementation was 
    experimental and results were not reliable at all.
    
Compilation
-----------

1. I'm on linux and I get compilation errors like

    ```
    gcc: error trying to exec 'cc1plus': execvp: No such file or directory
    error: Setup script exited with error: command 'gcc' failed with exit status 1
    ```
    
    This means that the g++ package is missing.  On OpenSUSE you can just install the following packages in software manager:
    
    ```
    gcc-c++
    gcc48-c++
    libstdc++48-devel
    ```
    
    On ubuntu, it should be sufficient to just run
    
    ```
    sudo apt-get install g++
    ```
    
    
