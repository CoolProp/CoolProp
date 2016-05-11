program project2;
{$apptype console}

uses
  SysUtils, Windows,
  cpIntf;

var
  T, h, p, D: Double;
begin

  {
  This example comes from Ian Bell.
  When I ran the original with Delphi, there were a few problems.
  In my cpIntf interface, the fpu exceptions are disabled
  and Ian's example works fine. 
  }

  Writeln('FLUID STATE INDEPENDENT INPUTS');

  Writeln('TWO PHASE INPUTS (Pressure)');
  Writeln(Format('Density of saturated liquid Propane at 101325 Pa: %g kg/m^3', [PropsSI('D', 'P', 101325, 'Q', 0, 'Propane')]));
  Writeln(Format('Density of saturated vapor R290 at 101325 Pa:     %g kg/m^3'+ sLineBreak, [PropsSI('D', 'P', 101325, 'Q', 1, 'R290')]));

  Writeln('TWO PHASE INPUTS (Temperature)');
  Writeln(Format('Density of saturated liquid Propane at 300 K: %g kg/m^3', [PropsSI('D', 'T', 300, 'Q', 0, 'Propane')]));
  Writeln(Format('Density of saturated vapor R290 at 300 K:     %g kg/m^3' + sLineBreak, [PropsSI('D', 'T', 300, 'Q', 1, 'R290')]));

  Writeln('SINGLE PHASE CYCLE (Propane)');
  p := PropsSI('P', 'T', 300, 'D', 1, 'Propane');
  h := PropsSI('H', 'T', 300, 'D', 1, 'Propane');
  Writeln(Format('T,D -> P,H : 300,1 -> %g,%g', [p, h]));
  T := PropsSI('T', 'P', p, 'H', h, 'Propane');
  D := PropsSI('D', 'P', p, 'H', h, 'Propane');
  Writeln(Format('P,H -> T,D : %g, %g -> %g, %g' + sLineBreak, [p, h, T, D]));

  Writeln('************ HUMID AIR PROPERTIES *************' + sLineBreak);
  Writeln(Format('Humidity ratio of 50%% rel. hum. air at 300 K, 101.325 kPa: %g kg_w/kg_da', [HAPropsSI('W', 'T', 300, 'P', 101325, 'R', 0.5)]));
  Writeln(Format('Relative humidity from last calculation: %g (fractional)', [HAPropsSI('R', 'T', 300, 'P', 101325, 'W', HAPropsSI('W', 'T', 300, 'P', 101325, 'R', 0.5))]));
  readln;
end.
