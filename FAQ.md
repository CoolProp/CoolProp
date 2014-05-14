
Frequently Asked Questions
==========================

Usage
-----

1. My enthalpy and entropy values are not the same as what are used in REFPROP, or EES, or ...

    Enthalpy and entropy are only equivalent to within an additive constant, so you should always compare enthalpy and entropy *differences* rather than enthalpy and entropy values
    
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
    
    