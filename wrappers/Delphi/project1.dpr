program project1;
{$apptype console}

uses
  SysUtils, Windows,
  cpIntf;

var
  h: Double;
begin

  {
   Very simple example
   Get the vapor enthalpy at a given pressure.
  }

  h := PropsSI('H', 'P', 101e3, 'Q', 1, 'R22');
  writeln(format('h = %0.3f kJ/kg',[1e-3*h]));
  readln;
end.
