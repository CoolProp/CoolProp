[v, o] = getversion();

// Idiotically, scilab can only handle functions in shared libraries that take
// all arguments by reference.  Why?  Who knows.  I've filed an issue about this
// and there was no reply.
//
// Conveniently, FORTRAN77 also requires that all functions take all arguments by
// reference, and we added FORTRAN77-compatible versions of HAPropsSI and PropsSI
// (but no other functions).  So these reference-only functions will work with
// scilab's crippled shared library interface.
//
// Below we have made small Scilab wrapper functions around these FORTRAN77-compatible
// functions for convenience.
if (getos() == "Windows") then
    if o(2) == "x64" then
        link('CoolProp_x64.dll', ['propssi_','hapropssi_'], 'c');
    else
        link('CoolProp.dll', ['propssi_','hapropssi_'], 'c');
    end
elseif (getos() == "Darwin") then
    link('libCoolProp.dylib', ['propssi_','hapropssi_'], 'c');
else // Linux
    link('libCoolProp.so', ['propssi_','hapropssi_'], 'c');
end

// Uncomment the following line to see what functions were found in the shared library
link('show');

funcprot(0)

function [out]=PropsSI(Output,Input1,Value1,Input2,Value2,Name);
  out = call("propssi_",Output,1,"c",Input1,2,"c",Value1,3,"d",Input2,4,"c",Value2,5,"d",Name,6,"c","out",[1,1],7,"d");
endfunction;

function [out]=HAPropsSI(Output,Input1,Value1,Input2,Value2,Input3,Value3);
  out = call("hapropssi_",Output,1,"c",Input1,2,"c",Value1,3,"d",Input2,4,"c",Value2,5,"d",Input3,6,"c",Value3,7,"d","out",[1,1],8,"d");
endfunction;

[rho] = PropsSI("D","T",298.15,"P",101325.0,"Air");
disp(rho)

[h] = HAPropsSI("H","T",298.15,"P",101325.0,"W",0.003);
disp(h)

//ulink();
