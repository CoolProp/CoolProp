{
  Delphi Interface for CoolProp
  Bruce Wernick, 10 December 2015
}

unit cpIntf;

interface

const
  cpdll = 'CoolProp.dll';

type
  CoolChar = AnsiChar;
  PCoolChar = PAnsiChar;
  CoolStr = AnsiString;
  CoolVec = array[0..1023] of CoolChar;

function get_global_param_string(
  param: PCoolChar;
  var res: CoolVec;
  n: Integer): LongInt; stdcall;
  external cpdll name '_get_global_param_string@12';

function get_fluid_param_string(
  fluid: PCoolChar; param: PCoolChar;
  var Output: CoolVec;
  n: Integer): Longint; stdcall;
  external cpdll name '_get_fluid_param_string@12';

function PropsSI(spec: PCoolChar;
  prop1: PCoolChar; val1: Double;
  prop2: PCoolChar; val2: Double;
  ref: CoolStr): Double; stdcall;
  external cpdll name '_PropsSI@32';

function Props1SI(spec: PCoolChar;
  ref: CoolStr): Double; stdcall;
  external cpdll name '_Props1SI@8';

function HAPropsSI(spec: PCoolChar;
  prop1: PCoolChar; val1: Double;
  prop2: PCoolChar; val2: Double;
  prop3: PCoolChar; val3: Double): Double; stdcall;
  external cpdll name '_HAPropsSI@40';

implementation

begin
end.
