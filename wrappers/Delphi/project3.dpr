program project3;
{$apptype console}

uses
  SysUtils, Windows,
  cpIntf;

var
  res: AnsiString;
  T: Double;
begin
  SetLength(res, 1024);
  get_global_param_string('FluidsList', PAnsiChar(res), 1024);
  writeln(res);
  
  T := PropsSI('T', 'P', 4.7e6, 'H', 330e3, 'R22')-273.15;
  writeln(T:6:2);
  readln;
end.