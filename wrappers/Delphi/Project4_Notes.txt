CoolProp Delphi Wrapper Demo.
Bruce Wernick, info@coolit.co.za
13 December 2015

This demo is based on the win32 shared DLL downloaded from the CoolProp site.

There are two parts to the demo:
  1. The actual DLL interface.
  2. A simple ph-chart builder.

The original code was written in C++.  In the source folder, you can find the text file exports.txt 
that details the complete list of functions exported by the CoolProp DLL.  I have only implemented 
the essential ones but it would be easy to add any ones I have omitted.  The interface mimics the 
C header file but there are a few options.  The returning string is an array of char.  In Delphi, 
you could use an array of char parameter.  Instead, I chose to use the PAnsiChar because it makes 
the dll interface simplere and more universal.  One thing you have to be aware of - the calling 
program has to create a buffer for the result.  If this is not sized sufficiently, then a blank string 
is returned.  This could happen, for example, if a lot of fluids are added to the file.

The graphical part draws the saturation curve on conventional log(p) vs enthalpy axes.  The program main 
form has a Delphi VCL TListView and a TPaintBox.  On startup, the program uses the DLL to populate the 
ListView with all of the available fluids.  The PaintBox is a very basic components that is basically 
just a canvas.  The main program creates a TMolChart class that does all the work on a canvas.  When 
the OnPaint fires in the main form, the chart is re-drawn.

The chart is auto scaled around the critical point and atmospheric pressure (where possible).  This 
makes it quite convenient for comparing the shape of the various saturation curves.  What it needs is 
some isotherms and some isentropic lines.  If anyone has the urge to add these I would like to get an 
update.


Lazarus users
With some minor changes, you could quite easily run this code with Lazarus.


