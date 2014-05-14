CoolProp for Excel
------------------

Ian Bell, 2012, ian.h.bell@gmail.com

Installation
------------
(for Excel 2010, though other versions follow a similar protocol)

Part 1:
^^^^^^^
1. Copy the file CoolProp.dll to c:\\CoolProp folder if you have 32-bit Excel, the CoolProp_x64.dll if you have 64-bit Excel.  It can't hurt to copy both if you don't know which version you have.

Part 2:
^^^^^^^
1. Open Excel
2. Go to the menu File-->Options-->Add-Ins
3. At the bottom, select Manage: Excel Add-ins, then click the Go.. button
4. Click the browse button
5. Browse to the file CoolProp.xlam you downloaded, select it
6. Make sure the CoolProp Add-in is selected.
7. Open the file TestExcel.xlsx and try to re-evaluate one of the cells - they should work now

Part 3 (optional):
^^^^^^^^^^^^^^^^^^
1. If you have REFPROP installed, make sure you have a copy of REFPROP installed 
   in c:\\Program Files\\REFPROP, if not, copy your REFPROP to this location.
   Alternatively you can also modify the paths in the TestExcel.xlam file
   
Notes
-----
The DLL is built using the __stdcall calling convention on 32-bit windows