program project6;
{$apptype console}

{ Not using cpIntf
  This is to demo the problem when the exception is not disabled.
  It happens with the Delphi compiler.
}

function PropsSI(spec: PAnsiChar;
                 prop1: PAnsiChar; val1: Double;
                 prop2: PAnsiChar; val2: Double;
                 ref: PAnsiChar): Double; stdcall;
                 external 'CoolProp.dll' name '_PropsSI@32';

const
  MCW_EM = Word($133f); // Disable all fpu exceptions
var
  Saved8087CW: Word;
  p, h, t: double;
begin
  Saved8087CW := Default8087CW;
  Set8087CW(MCW_EM);
  h := 330;
  p := 4.3; // no problem
  t := PropsSI('T', 'P', 1e6*p, 'H', 1e3*h, 'R22');
  writeln(t);
  p := 4.4; // crash if fpu not disabled
  t := PropsSI('T', 'P', 1e6*p, 'H', 1e3*h, 'R22');
  Set8087CW(Saved8087CW);
  writeln(t);
  readln;
end.
