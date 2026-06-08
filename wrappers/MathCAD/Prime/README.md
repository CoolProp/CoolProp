CoolProp Wrapper for MathCAD Prime 7.0 or later (64-bit)
==========================================================

| Copyright Scott Polak and Ian Bell, 2013  
| Updated by Jeff Henning, 2016

To Install the CoolProp Add-in
==============================

* Download the CoolPropMathcadWrapper.dll file from SourceForge and copy to your "C:\\Program Files\\PTC\\Mathcad Prime ``x.0.0.0``\\Custom Functions" folder, where ``x.0.0.0`` is replaced by your working version of Mathcad Prime.

* Alternatively, build the DLL from the Coolprop source following the instructions in the [Mathcad wrapper](https://github.com/CoolProp/CoolProp/tree/master/wrappers/MathCAD) directory on Github.

* (Optional) Install the Custom Functions Add-in, [CustFunc](https://github.com/henningjp/CustFunc) and copy the associated ``CoolProp_EN.xml`` file to the ``Custom Functions\docs`` folder, creating it if it doesn't exist (this file is optional and for features in development). The [CustFunc](https://github.com/henningjp/CustFunc) add-in provides a pop-up (when pressing ``<F3>``) that gives brief descriptions of each implemented function and its input parameters and facilitates inserting these functions into the worksheet.

* Restart Mathcad Prime.

Using the Example File
======================

The Mathcad Prime file ``CoolPropFluidProperties.mcdx`` demonstrates how to use the CoolProp high-level and some low-level API calls from Mathcad Prime.  The file is saved in Prime 7.0 format so that it can be read in any versions from 7.0 onward.  In the GitHub CoolProp repository there is also a PDF of the file for viewing without Mathcad Prime.

* Open the CoolPropFluidProperties.mcdx file in Mathcad Prime; all CoolProp functions should revaluate properly.  

> _**IMPORTANT!**_  
> Once the example file is loaded, press **`<Ctrl>-<F9>`** to force recalculation of the entire workbook.  

* The examples file should cover most functions included in the Mathcad wrapper, but check the [High-Level Interface documentation at Coolprop.org](https://coolprop.org/coolprop/HighLevelAPI.html#) for more details and a full table of string inputs/outputs to **PropsSI** and **Prop1SI**.  