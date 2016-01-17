program Project5;
{$apptype console}

uses
  cpIntf;

var
  dat: AnsiString;
begin
  SetLength(dat, 1024);
  get_global_param_string('FluidsList', PAnsiChar(dat), 1024);
  writeln(dat);
  readln;
end.  