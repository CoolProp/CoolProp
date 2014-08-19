link('CoolProp.dll', ['F77PropsSI','F77HAPropsSI'], 'c');
//link('show');

funcprot(0)

function [out]=PropsSI(Output,Input1,Value1,Input2,Value2,Name)
  out = call("F77PropsSI",Output,1,"c",Input1,2,"c",Value1,3,"d",Input2,4,"c",Value2,5,"d",Name,6,"c","out",[1,1],7,"d")
endfunction

function [out]=HAPropsSI(Output,Input1,Value1,Input2,Value2,Input3,Value3)
  out = call("F77HAPropsSI",Output,1,"c",Input1,2,"c",Value1,3,"d",Input2,4,"c",Value2,5,"d",Input3,6,"c",Value3,7,"d","out",[1,1],8,"d")
endfunction

[rho] = PropsSI("D","T",298.15,"P",101325.0,"Air")
disp(rho)

[h] = HAPropsSI("H","T",298.15,"P",101325.0,"W",0.003)
disp(h)

ulink()
