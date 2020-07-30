
Frequently Asked Questions
==========================

Usage
-----

1. **My enthalpy and entropy values are not the same as what are used in REFPROP, or EES, or ...**

    Values for enthalpy, entropy and internal energy are calculated as *differences* 
    with respect to an arbitrarily chosen reference state. 
    As a consequence, absolute values may differ hugely between different tools, 
    but the difference between two states should be equal.  
    At least the following three states are frequently chosen as reference states:  
    IIR: saturated liquid (Q=0) at T=273.15 K  
    ASHRAE: saturated liquid (Q=0) at T=233.15 K  
    NBP: saturated liquid (Q=0) at p=101325.0 Pa
    
    
2. **My calls to REFPROP fail with a mysterious `Segmentation Fault`:** 

    Make sure that you have at least REFPROP version 9.1 installed. If you have a license for 9.0, 
    the upgrade is for free: http://www.nist.gov/srd/nist23.cfm
    
3. **Transport properties missing**

    a. Fluid XY does not have viscosity/thermal conductivity/ ... 

       Please have a look at 
       http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids 
       if there is no reference for your property, we do not have the information required to 
       implement the functionality. If you file an issue, please provide information and 
       a link to a publication with the required data. 
    
    b. ..but it used to work in an earlier version of CoolProp
    
       For some fluids, version 5 does not have transport properties even though version 4 had 
       that information. We decided to remove some correlations, which had very large uncertainties 
       and we would like to make the user aware of the fact that the old implementation was 
       experimental and results were not reliable at all.

4. **Mixture calculation gives error: Could not match the binary pair [0000-00-0,11-11-1] - for now this is an error.**  

    Mixture calculations require binary interaction parameters for each pair in the mixture.  While many binary interaction parameters are available in the CoolProp library, sadly, many are not.  If you get this error message, then the binary interaction parameters for the CAS fluids listed are not available in CoolProp.  
    
    If you have data for the binary interaction parameters, you can enter them interactively using the [set_mixture_binary_pair_data](http://www.coolprop.org/dev/fluid_properties/Mixtures.html#id826) function in CoolProp.  Otherwise, a more sophisticated mixing model is needed, like the ones in NIST RefProp.
    
5. **Where is the MATLAB wrapper?**

    We retired the MATLAB wrapper in favour of a Python-based approach, read more in the [docs](http://www.coolprop.org/coolprop/wrappers/MATLAB/index.html#matlab-wrapper).

6. **Where is the enthalpy of mixing for incompressible fluids?**

    This is a major flaw in the incompressible backend. The incompressible fluids module was designed to handle heat transfer fluids with fixed compositions and you should not compare enthalpy and entropy values for different concentrations. The way the reference state is calculated makes it impossible to properly handle mixing enthalpy etc.

Compilation
-----------

1. **I'm on Linux and I get compilation errors like**

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
    
2. **Building Python wrapper on Windows fails with an error similar to**

    ```
    error: Microsoft Visual C++ 14.0 is required. Get it with "Microsoft Visual C++ Build Tools": https://visualstudio.microsoft.com/visual-cpp-build-tools/
    ```
    
    Different versions of Python (2.7, 3.6, etc.) are built with and require specific versions of the Microsoft Visual C++ compiler.  Please see the [common wrapper prerequisites](http://www.coolprop.org/dev/coolprop/wrappers/index.html#wrapper-common-prereqs) specifically for Windows build requirements.
    
