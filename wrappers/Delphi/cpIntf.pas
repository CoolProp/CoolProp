{
  Delphi Interface for CoolProp
  Bruce Wernick, 13 December 2015

  Notes:
  All of the char parameters are PAnsiChar, the floating point numbers are
  Double and the integer types are either LongInt or Integer depending on the
  type in the dll.  Where there is an output string, the calling program must
  create the space for the result.  Be aware here that insufficient space will
  result in a empty string.  So for example, you should ensure sufficient space
  for the refrigerant list.
}

unit cpIntf;

interface

const
  cpdll = 'CoolProp.dll';

function get_global_param_string(
  param: PAnsiChar;
  res: PAnsiChar;
  n: Integer): LongInt; stdcall;
  external cpdll name '_get_global_param_string@12';

function get_fluid_param_string(
  fluid: PAnsiChar;
  param: PAnsiChar;
  res: PAnsiChar;
  n: Integer): Longint; stdcall;
  external cpdll name '_get_fluid_param_string@12';

function PropsSI(
  spec: PAnsiChar;
  prop1: PAnsiChar; val1: Double;
  prop2: PAnsiChar; val2: Double;
  fluid: PAnsiChar): Double; stdcall;
  external cpdll name '_PropsSI@32';

function Props1SI(
  fluid: PAnsiChar;
  res: PAnsiChar): Double; stdcall;
  external cpdll name '_Props1SI@8';

function HAPropsSI(
  spec: PAnsiChar;
  prop1: PAnsiChar; val1: Double;
  prop2: PAnsiChar; val2: Double;
  prop3: PAnsiChar; val3: Double): Double; stdcall;
  external cpdll name '_HAPropsSI@40';

implementation

initialization
  {Disable floating point exception
   This is a problem with the Delphi compiler only}
  Set8087CW($133f);
end.
