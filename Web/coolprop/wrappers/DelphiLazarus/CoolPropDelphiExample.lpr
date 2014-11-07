program CoolPropDelphiExample;

{$APPTYPE CONSOLE}
{$R *.res}

uses
  SysUtils, Windows;

CONST
  COOL_PROP = 'CoolProp.dll';

type
  TArrayOfAnsiChar = array [0 .. 1023] of AnsiChar;

function get_global_param_string(param: PAnsichar; var Output: TArrayOfAnsiChar): Longint; stdcall; external COOL_PROP name '_get_global_param_string@8';
function get_fluid_param_string(fluid: PAnsichar; param: PAnsichar; var Output: TArrayOfAnsiChar): Longint; stdcall; external COOL_PROP name '_get_fluid_param_string@12';
function Props1(FluidName: PAnsichar; Output: PAnsichar): Double; stdcall; external COOL_PROP name '_Props1@8';
function PropsSI(Output: PAnsichar; Name1: PAnsichar; Prop1: Double; Name2: PAnsichar; Prop2: Double; FluidName: PAnsichar): Double; stdcall; external COOL_PROP name '_PropsSI@32';
function HAProps(Output: PAnsichar; Name1: PAnsichar; Prop1: Double; Name2: PAnsichar; Prop2: Double; Name3: PAnsichar; Prop3: Double): Double; stdcall; external COOL_PROP name '_HAProps@40';
function enable_TTSE_LUT(FluidName: PAnsichar): Boolean; stdcall; external COOL_PROP name '_enable_TTSE_LUT@4';
function disable_TTSE_LUT(FluidName: PAnsichar): Boolean; stdcall; external COOL_PROP name '_disable_TTSE_LUT@4';

var
  OutputArrayOfAnsiChar: TArrayOfAnsiChar;
  OutputString: string;
  Saved8087CW: Word;
//  RPName: PAnsichar;
  T, h, p, D: Double;

begin
    Saved8087CW := Default8087CW;
    get_global_param_string('version', OutputArrayOfAnsiChar);
    OutputString := string(OutputArrayOfAnsiChar);
    Writeln(Format('CoolProp version:       %s', [OutputString]));
    get_global_param_string('gitrevision', OutputArrayOfAnsiChar);
    OutputString := string(OutputArrayOfAnsiChar);
    Writeln(Format('CoolProp gitrevision:   %s', [OutputString]));
    get_global_param_string('FluidsList', OutputArrayOfAnsiChar);
    OutputString := string(OutputArrayOfAnsiChar);
    Writeln(Format('CoolProp fluids:        %s' + sLineBreak, [OutputString]));

    Writeln('************ USING EOS *************' + sLineBreak);
    Writeln('FLUID STATE INDEPENDENT INPUTS');
    Writeln(Format('Critical Density Propane: %g kg/m^3' + sLineBreak, [Props1('Propane', 'rhocrit')]));
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

    Writeln('************ BRINES AND SECONDARY WORKING FLUIDS *************' + sLineBreak);
    Writeln(Format ('Density of 50%% (mass) ethylene glycol/water at 300 K, 101325 Pa: %g kg/m^3', [PropsSI('D', 'T', 300, 'P', 101325, 'EG-50%')]));
    Writeln(Format('Viscosity of Therminol D12 at 350 K, 101325 Pa: %g Pa-s' + sLineBreak, [PropsSI('V', 'T', 350, 'P', 101325, 'TD12')]));

    Writeln('************ HUMID AIR PROPERTIES *************' + sLineBreak);
    Writeln(Format('Humidity ratio of 50%% rel. hum. air at 300 K, 101.325 kPa: %g kg_w/kg_da', [HAProps('W', 'T', 300, 'P', 101.325, 'R', 0.5)]));
    Writeln(Format('Relative humidity from last calculation: %g (fractional)', [HAProps('R', 'T', 300, 'P', 101.325, 'W', HAProps('W', 'T', 300, 'P', 101.325, 'R', 0.5))]));
end.
