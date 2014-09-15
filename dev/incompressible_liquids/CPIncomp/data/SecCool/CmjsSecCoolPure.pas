unit CmjsSecCoolPure;

interface

uses
  SysUtils,CmjsSecCoolFluid;

type
  TDowtherm_J =
    class(TmjsSecCoolFluid)
      private
      protected
        function  GetRhoT(T : Double): Double; override;
        function  GetCpT(T : Double): Double; override;
        function  GetCondT(T : Double): Double; override;
        function  GetMuT(T : Double): Double; override;
      public
        constructor Create; override;
    end;

  TDowtherm_Q =
    class(TmjsSecCoolFluid)
      private
      protected
        function  GetRhoT(T : Double): Double; override;
        function  GetCpT(T : Double): Double; override;
        function  GetCondT(T : Double): Double; override;
        function  GetMuT(T : Double): Double; override;
      public
        constructor Create; override;
    end;

  TSyltherm_XLT =
    class(TmjsSecCoolFluid)
      private
      protected
        function  GetRhoT(T : Double): Double; override;
        function  GetCpT(T : Double): Double; override;
        function  GetCondT(T : Double): Double; override;
        function  GetMuT(T : Double): Double; override;
      public
        constructor Create; override;
    end;

  TThermogen_VP1869 =
    class(TmjsSecCoolFluid)
      private
      protected
        function  GetRhoT(T : Double): Double; override;
        function  GetCpT(T : Double): Double; override;
        function  GetCondT(T : Double): Double; override;
        function  GetMuT(T : Double): Double; override;
      public
        constructor Create; override;
    end;

  TTemper_10 =
    class(TmjsSecCoolFluid)
      private
      protected
        function  GetRhoT(T : Double): Double; override;
        function  GetCpT(T : Double): Double; override;
        function  GetCondT(T : Double): Double; override;
        function  GetMuT(T : Double): Double; override;
      public
        constructor Create; override;
    end;

  TTemper_20 =
    class(TmjsSecCoolFluid)
      private
      protected
        function  GetRhoT(T : Double): Double; override;
        function  GetCpT(T : Double): Double; override;
        function  GetCondT(T : Double): Double; override;
        function  GetMuT(T : Double): Double; override;
      public
        constructor Create; override;
    end;

  TTemper_30 =
    class(TmjsSecCoolFluid)
      private
      protected
        function  GetRhoT(T : Double): Double; override;
        function  GetCpT(T : Double): Double; override;
        function  GetCondT(T : Double): Double; override;
        function  GetMuT(T : Double): Double; override;
      public
        constructor Create; override;
    end;

  TTemper_40 =
    class(TmjsSecCoolFluid)
      private
      protected
        function  GetRhoT(T : Double): Double; override;
        function  GetCpT(T : Double): Double; override;
        function  GetCondT(T : Double): Double; override;
        function  GetMuT(T : Double): Double; override;
      public
        constructor Create; override;
    end;

  TTemper_55 =
    class(TmjsSecCoolFluid)
      private
      protected
        function  GetRhoT(T : Double): Double; override;
        function  GetCpT(T : Double): Double; override;
        function  GetCondT(T : Double): Double; override;
        function  GetMuT(T : Double): Double; override;
      public
        constructor Create; override;
    end;

  THyCool_20 =
    class(TmjsSecCoolFluid)
      private
      protected
        function  GetRhoT(T : Double): Double; override;
        function  GetCpT(T : Double): Double; override;
        function  GetCondT(T : Double): Double; override;
        function  GetMuT(T : Double): Double; override;
      public
        constructor Create; override;
    end;

  THyCool_30 =
    class(TmjsSecCoolFluid)
      private
      protected
        function  GetRhoT(T : Double): Double; override;
        function  GetCpT(T : Double): Double; override;
        function  GetCondT(T : Double): Double; override;
        function  GetMuT(T : Double): Double; override;
      public
        constructor Create; override;
    end;

  THyCool_40 =
    class(TmjsSecCoolFluid)
      private
      protected
        function  GetRhoT(T : Double): Double; override;
        function  GetCpT(T : Double): Double; override;
        function  GetCondT(T : Double): Double; override;
        function  GetMuT(T : Double): Double; override;
      public
        constructor Create; override;
    end;

  THyCool_45 =
    class(TmjsSecCoolFluid)
      private
      protected
        function  GetRhoT(T : Double): Double; override;
        function  GetCpT(T : Double): Double; override;
        function  GetCondT(T : Double): Double; override;
        function  GetMuT(T : Double): Double; override;
      public
        constructor Create; override;
    end;

  THyCool_50 =
    class(TmjsSecCoolFluid)
      private
      protected
        function  GetRhoT(T : Double): Double; override;
        function  GetCpT(T : Double): Double; override;
        function  GetCondT(T : Double): Double; override;
        function  GetMuT(T : Double): Double; override;
      public
        constructor Create; override;
    end;

  THydroFluoroEther =
    class(TmjsSecCoolFluid)
      private
      protected
        function  GetRhoT(T : Double): Double; override;
        function  GetCpT(T : Double): Double; override;
        function  GetCondT(T : Double): Double; override;
        function  GetMuT(T : Double): Double; override;
      public
        constructor Create; override;
    end;

implementation

uses
  CmjsSecCool;

type
  coef = array[0..26] of Double;

function evalcratl2(order:integer; var x: Double; var c: coef) : Double;
var
  i, j, k: integer;
  xs, n, d, tmp, t2, t1 : Double;
begin
  xs := (x-394.1499999999999773)/(194.0000000000000000);
  i  := order;
  j  := order;
  if (order mod 2)=1 then
    i :=i-1
  else
    j :=j-1;
  t2 := 0.0;
  t1 := 0.0;
  k  := i;
  while k>1 do
    begin
      tmp := t1;
      t1  := 2*xs*t1-t2+c[k];
      t2  := tmp;
      k   := k-2;
    end;
  n  := xs*t1-t2+c[k];
  t2 := 0.0;
  t1 := 0.0;
  k  := j;
  while k>=1 do
    begin
      tmp := t1;
      t1  := 2*xs*t1-t2+c[k];
      t2  := tmp;
      k   := k-2;
    end;
  d := xs*t1-t2+1.0;
  if d=0.0 then
    Result := 0.0
  else
    Result := n/d;
end;

{ TDowtherm_J }

constructor TDowtherm_J.Create;
begin
  inherited;
  TFreeze   := -86;
  FluidType := 'Diethylbenzene mixture';
  TradeName := 'Dowtherm J';
  Reference := 'DOW Chemical Company';
  TMin      := -73;
  TMax      := 315;
end;

function TDowtherm_J.GetCondT(T: Double): Double;
begin
  Result := 0.1567891965299059-8.197929309283446E-05*(T+273.15);
end;

function TDowtherm_J.GetCpT(T: Double): Double;
{ Eqn# 6004  y=a+bx+cx^2+dx^3+ex^4+fx^5+gx^6+hx^7 }
const
 a= -5476.30952146975;      b= 146.3997303809014;     c= -1.309135572360928;
 d= 0.006423251883522429;   e= -1.84410593468011E-05; f= 3.111472792972415E-08;
 g= -2.860712306793638E-11; h= 1.106640602427616E-14;
begin
  T := T+273.15;
  Result := a+T*(b+T*(c+T*(d+T*(e+T*(f+T*(g+T*h))))));
end;

function TDowtherm_J.GetMuT(T: Double): Double;
{ Eqn# 7613  Chebyshev Rational Order 6/7 }
var
  c : coef;
begin
  c[0]   := 0.001064371051088026542;   c[1]  := 0.4332761878651155856;
  c[2]   := -0.001428722989359673563;  c[3]  := 0.2955792635389909284;
  c[4]   := 0.0009713993787354970530;  c[5]  := 0.3504194915584169806;
  c[6]   := -0.0003829302988434106282; c[7]  := 0.04935576384982437492;
  c[8]   := 0.0001456080745522824655;  c[9]  := 0.1538488081319953276;
  c[10]  := 8.907959978832755140E-06;  c[11] := 0.10452520257749200515;
  c[12]  := 1.283975792724080212E-05;  c[13] := 0.004023249322418941525;
  T      := T+273.15;
  Result := (evalcratl2(13,T,c))*1e3;
end;

function TDowtherm_J.GetRhoT(T: Double): Double;
{ Eqn# 6006  y=a+bx+cx^2+dx^3+ex^4+fx^5+gx^6+hx^7+ix^8+jx^9 }
const
 a = -1354.474533119642;     b = 62.37124259565723;      c = -0.7189258383359517;
 d = 0.00470925809684546;    e = -1.952403154224894E-05; f = 5.312254521315433E-08;
 g = -9.487853242180682E-11; h = 1.072667547273402E-13;  i = -6.965808790863649E-17;
 j = 1.979206386602891E-20;
begin
  T      := T+273.15;
  Result := a + T*(b+T*(c+T*(d+T*(e+T*(f+T*(g+T*(h+T*(i+j*T))))))));
end;

{ TDowtherm_Q }

constructor TDowtherm_Q.Create;
begin
  inherited;
  TFreeze   := -40;
  FluidType := 'Diphenylethane/alkylated aromatics';
  TradeName := 'Dowtherm Q';
  Reference := 'DOW Chemical Company';
  TMin      := -35;
  TMax      := 20;
end;

function TDowtherm_Q.GetCondT(T: Double): Double;
begin
  Result := 0.12689093+T*(0.000110253+T*(-1.30003E-06+T*(-3.63756E-08+
            T*(-1.7094E-09+T*(-5.26395E-11-6.53595E-13*T)))))
end;

function TDowtherm_Q.GetCpT(T: Double): Double;
begin
  Result := (1.589770396+0.003191608*T)*1000;
end;

function TDowtherm_Q.GetMuT(T: Double): Double;
begin
  Result := (836.2199774-12.94719372*T)/(1+0.021426081*T)*1e-2;
end;

function TDowtherm_Q.GetRhoT(T: Double): Double;
begin
  Result := 986.2587413-0.732167832*T;
end;

{ TSyltherm_XLT }

constructor TSyltherm_XLT.Create;
begin
  inherited;
  TFreeze   := -93;
  FluidType := 'Polydimethylsiloxan';
  TradeName := 'Syltherm XLT';
  Reference := 'DOW Chemical Company';
  TMin      := -70;
  TMax      := 30;
end;

function TSyltherm_XLT.GetCondT(T: Double): Double;
begin
  Result := 0.107155405-7.09459E-05*T;
end;

function TSyltherm_XLT.GetCpT(T: Double): Double;
begin
  Result := (1.598216216+0.002562162*T)*1000;
end;

function TSyltherm_XLT.GetMuT(T: Double): Double;
begin
  Result := (186.2068182-1.665321134*T)/(1+0.010123277*T)*1e-2;
end;

function TSyltherm_XLT.GetRhoT(T: Double): Double;
begin
  Result := 862.0608108-0.866891892*T;
end;

{ TThermogen_VP1869 }

constructor TThermogen_VP1869.Create;
begin
  inherited;
  TFreeze   := -95;
  FluidType := 'Thermogen VP 1869';
  TradeName := 'Thermogen VP 1869';
  Reference := 'Hoechst';
  TMin      := -80;
  TMax      := 20;
end;

function TThermogen_VP1869.GetCondT(T: Double): Double;
begin
  Result := 0.15-0.000154545*T;
end;

function TThermogen_VP1869.GetCpT(T: Double): Double;
begin
  Result := (2.322218182+0.003843636*T)*1000
end;

function TThermogen_VP1869.GetMuT(T: Double): Double;
begin
  Result := (341.3688975+T*(-0.713408301+0.017723992*T))/
            (1+T*(0.034502393+T*(0.000401319+1.57288E-06*T)))*1e-2;
end;

function TThermogen_VP1869.GetRhoT(T: Double): Double;
begin
  Result := 945.5454545-1.054545455*T;
end;

{ TTemper_10 }

constructor TTemper_10.Create;
begin
  inherited;
  TFreeze   := -10;
  FluidType := 'Potassium acetate/formate';
  TradeName := 'Aspen Temper -10';
  Reference := 'Aspen Petroleum AB';
  TMin      := -10;
  TMax      := 30;
end;

function TTemper_10.GetCondT(T: Double): Double;
begin
  Result := 0.001483*(273+T)+0.1094;
end;

function TTemper_10.GetCpT(T: Double): Double;
begin
  Result := 1000*(3.54183+T*(0.00201-0.0000132589*T));
end;

function TTemper_10.GetMuT(T: Double): Double;
begin
  Result := (2.80154921+T*(-0.10631376+T*(0.00235381-0.0000211481*T)));
end;

function TTemper_10.GetRhoT(T: Double): Double;
begin
  Result := 1090+T*(-0.2+T*(1.66533E-16-1.89735E-18*T));
end;

{ TTemper_20 }

constructor TTemper_20.Create;
begin
  inherited;
  TFreeze   := -20;
  FluidType := 'Potassium acetate/formate';
  TradeName := 'Aspen Temper -20';
  Reference := 'Aspen Petroleum AB';
  TMin      := -20;
  TMax      := 30;
end;

function TTemper_20.GetCondT(T: Double): Double;
begin
  Result := 0.001342*(273+T)+0.1144;
end;

function TTemper_20.GetCpT(T: Double): Double;
begin
  Result := 1000*(3.26252+T*(0.00286-0.0000122673*T));
end;

function TTemper_20.GetMuT(T: Double): Double;
begin
  Result := (0.97244+(7.14199*EXP((-20-T)/18.6014)));
end;

function TTemper_20.GetRhoT(T: Double): Double;
begin
  Result := 1147+T*(-0.22142857+T*(-0.00142857+2.98156E-19*T));
end;

{ TTemper_30 }

constructor TTemper_30.Create;
begin
  inherited;
  TFreeze   := -30;
  FluidType := 'Potassium acetate/formate';
  TradeName := 'Aspen Temper -30';
  Reference := 'Aspen Petroleum AB';
  TMin      := -30;
  TMax      := 30;
end;

function TTemper_30.GetCondT(T: Double): Double;
begin
  Result := 0.001256*(273+T)+0.1175;
end;

function TTemper_30.GetCpT(T: Double): Double;
begin
  Result := 1000*(3.07504+T*(0.00299-0.0000270232*T));
end;

function TTemper_30.GetMuT(T: Double): Double;
begin
  Result := (1.30143+(16.01368*EXP((-30-T)/16.69989)));
end;

function TTemper_30.GetRhoT(T: Double): Double;
begin
  Result := 1183.85930736+T*(-0.31085859+T*(-0.00103896+0.0000176768*T));
end;

{ TTemper_40 }

constructor TTemper_40.Create;
begin
  inherited;
  TFreeze   := -40;
  FluidType := 'Potassium acetate/formate';
  TradeName := 'Aspen Temper -40';
  Reference := 'Aspen Petroleum AB';
  TMin      := -40;
  TMax      := 30;
end;

function TTemper_40.GetCondT(T: Double): Double;
begin
  Result := 0.001099*(273+T)+0.1433;
end;

function TTemper_40.GetCpT(T: Double): Double;
begin
  Result := 1000*(2.97788+T*(0.00228-0.0000387227*T));
end;

function TTemper_40.GetMuT(T: Double): Double;
begin
  Result := (39.11536*EXP((-40-T)/9.99495)+(12.41564*EXP((-40-T)/38.45577)));
end;

function TTemper_40.GetRhoT(T: Double): Double;
begin
  Result := 1214.83982684+T*(-0.37819865+T*(-0.00109307+0.000016835*T));
end;

{ TTemper_55 }

constructor TTemper_55.Create;
begin
  inherited;
  TFreeze   := -55;
  FluidType := 'Potassium acetate/formate';
  TradeName := 'Aspen Temper -55';
  Reference := 'Aspen Petroleum AB';
  TMin      := -55;
  TMax      := 30;
end;

function TTemper_55.GetCondT(T: Double): Double;
begin
  Result := (273+T)*(0.000002287*(273+T)-0.0003108)+0.3402;
end;

function TTemper_55.GetCpT(T: Double): Double;
begin
  Result := 1000*(2.83985+T*(0.00229-0.0000248618*T));
end;

function TTemper_55.GetMuT(T: Double): Double;
begin
  Result := (317.40673*EXP((-55-T)/7.24125)+(51.22151*EXP((-55-T)/26.28052)));
end;

function TTemper_55.GetRhoT(T: Double): Double;
begin
  Result := 1249.7534665+T*(-0.47629615+T*(-0.00117189+0.0000198824*T));
end;

{ THyCool_20 }

constructor THyCool_20.Create;
begin
  inherited;
  TFreeze   := -20;
  FluidType := 'Potassium Formate';
  TradeName := 'HYCOOL 20';
  Reference := 'Hydro Chemicals';
  TMin      := -20;
  TMax      := 50;
end;

function THyCool_20.GetCondT(T: Double): Double;
begin
  if T <= 20 then
    Result := 0.001978*T+0.5172
  else
    Result := 0.001005*T+0.5368;
end;

function THyCool_20.GetCpT(T: Double): Double;
begin
  Result := 1000*(0.0023*T+2.955);
end;

function THyCool_20.GetMuT(T: Double): Double;
begin
  if T <= 20 then
    Result := 0.07190*exp(524.75/(T+142.05))
  else
    Result := T*(0.0005524*T - 0.06281)+2.8536;
end;

function THyCool_20.GetRhoT(T: Double): Double;
begin
  Result := -0.429180*T+1202.2;
end;

{ THyCool_30 }

constructor THyCool_30.Create;
begin
  inherited;
  TFreeze   := -30;
  FluidType := 'Potassium Formate';
  TradeName := 'HYCOOL 30';
  Reference := 'Hydro Chemicals';
  TMin      := -30;
  TMax      := 50;
end;

function THyCool_30.GetCondT(T: Double): Double;
begin
  if T <= 20 then
    Result := 0.001840*T+0.4980
  else
    Result := 0.001000*T+0.514;
end;

function THyCool_30.GetCpT(T: Double): Double;
begin
  Result := 1000*(0.0023*T+2.783);
end;

function THyCool_30.GetMuT(T: Double): Double;
begin
  if T <= 20 then
    Result := 0.11100*exp(408.17/(T+123.15))
  else
    Result := T*(0.000295*T - 0.0441)+2.6836;
end;

function THyCool_30.GetRhoT(T: Double): Double;
begin
  Result := -0.475350*T+1257.5;
end;

{ THyCool_40 }

constructor THyCool_40.Create;
begin
  inherited;
  TFreeze   := -40;
  FluidType := 'Potassium Formate';
  TradeName := 'HYCOOL 40';
  Reference := 'Hydro Chemicals';
  TMin      := -40;
  TMax      := 20;
end;

function THyCool_40.GetCondT(T: Double): Double;
begin
  Result := 0.001730*T+0.4826;
end;

function THyCool_40.GetCpT(T: Double): Double;
begin
  Result := 1000*(0.0023*T+2.646);
end;

function THyCool_40.GetMuT(T: Double): Double;
begin
  Result := 0.07830*exp(498.13/(T+130.25));
end;

function THyCool_40.GetRhoT(T: Double): Double;
begin
  Result := -0.512290*T+1304.5;
end;

{ THyCool_45 }

constructor THyCool_45.Create;
begin
  inherited;
  TFreeze   := -45;
  FluidType := 'Potassium Formate';
  TradeName := 'HYCOOL 45';
  Reference := 'Hydro Chemicals';
  TMin      := -40;
  TMax      := 20;
end;

function THyCool_45.GetCondT(T: Double): Double;
begin
  Result := 0.001674*T+0.4750;
end;

function THyCool_45.GetCpT(T: Double): Double;
begin
  Result := 1000*(0.0023*T+2.578);
end;

function THyCool_45.GetMuT(T: Double): Double;
begin
  Result := 0.08990*exp(479.09/(T+126.55));
end;

function THyCool_45.GetRhoT(T: Double): Double;
begin
  Result := -0.530754*T+1328.7;
end;

{ THyCool_50 }

constructor THyCool_50.Create;
begin
  inherited;
  TFreeze   := -50;
  FluidType := 'Potassium Formate';
  TradeName := 'HYCOOL 50';
  Reference := 'Hydro Chemicals';
  TMin      := -50;
  TMax      := 20;
end;

function THyCool_50.GetCondT(T: Double): Double;
begin
  Result := 0.001610*T+0.4660;
end;

function THyCool_50.GetCpT(T: Double): Double;
begin
  Result := 1000*(0.0023*T+2.498);
end;

function THyCool_50.GetMuT(T: Double): Double;
begin
  Result := 0.0491*exp(581.12/(T+129.05));
  if T > -10 then
    Result := Result+0.2;
end;

function THyCool_50.GetRhoT(T: Double): Double;
begin
  Result := -0.552300*T+1359.0;
end;

{ THydroFluoroEther }

constructor THydroFluoroEther.Create;
begin
  inherited;
  TFreeze   := -135;
  FluidType := 'Hydrofluoroether';
  TradeName := 'HFE-7100';
  Reference := '3M Novec';
  TMin      := -80;
  TMax      := 60;
end;

function THydroFluoroEther.GetCondT(T: Double): Double;
begin
  Result := -0.00019548*T+0.073714;
end;

function THydroFluoroEther.GetCpT(T: Double): Double;
begin
 Result := 2.00*T+1133;
end;

function THydroFluoroEther.GetMuT(T: Double): Double;
{ TableCurve Function:D:\Fluid properties\3M Novec Engineered Fluid HFE-7100\NuTcSt.pas Aug 27, 2003 3:30:50 PM }
{ D:\Fluid properties\3M Novec Engineered Fluid HFE-7100\Nu.xls }
{ X= Tc }
{ Y= Mu }
{ Eqn# 7104  lny=(a+cx+ex^2)/(1+bx+dx^2+fx^3) }
{ r2=0.9999682367141715 }
{  r2adj=0.9999523550712573 }
{  StdErr=0.006930163661963493 }
{  Fstat=81852.91123515403 }
const
  a = -0.6799203395116083;     b = 0.01648053800797168;
  c = -0.02549521321469719;    d = 6.358175288949973E-05;
  e = -0.0002027441598398813;  f = 8.339855339776066E-08;
begin
  Result := (a+T*(c+T*e))/(1+T*(b+T*(d+T*f)));
  Result := Exp(Result)*GetRhoT(T)*1e-3;
end;

function THydroFluoroEther.GetRhoT(T: Double): Double;
begin
 Result := -2.2690*T+1538.3;
end;

initialization
  //TmjsSecCool.RegisterFluid(False,'Diethylbenzene mixture','','Dowtherm J','','',TDowtherm_J);
  //TmjsSecCool.RegisterFluid(False,'Diphenylethane/alkylated aromatics','','Dowtherm Q','','',TDowtherm_Q);
  //TmjsSecCool.RegisterFluid(False,'Polydimethylsiloxan','','Syltherm XLT','','',TSyltherm_XLT);
  //TmjsSecCool.RegisterFluid(False,'Thermogen VP 1869','','Thermogen VP 1869','','',TThermogen_VP1869);
  TmjsSecCool.RegisterFluid(False,'Potassium acetate/formate','CAS 127-08-2 590-29-4','Aspen Temper -10','','',TTemper_10);
  TmjsSecCool.RegisterFluid(False,'Potassium acetate/formate','','Aspen Temper -20','','',TTemper_20);
  TmjsSecCool.RegisterFluid(False,'Potassium acetate/formate','','Aspen Temper -30','','',TTemper_30);
  TmjsSecCool.RegisterFluid(False,'Potassium acetate/formate','','Aspen Temper -40','','',TTemper_40);
  TmjsSecCool.RegisterFluid(False,'Potassium acetate/formate','','Aspen Temper -55','','',TTemper_55);
  TmjsSecCool.RegisterFluid(False,'Potassium Formate','','HYCOOL 20','','',THyCool_20);
  TmjsSecCool.RegisterFluid(False,'Potassium Formate','','HYCOOL 30','','',THyCool_30);
  TmjsSecCool.RegisterFluid(False,'Potassium Formate','','HYCOOL 40','','',THyCool_40);
  TmjsSecCool.RegisterFluid(False,'Potassium Formate','','HYCOOL 45','','',THyCool_45);
  TmjsSecCool.RegisterFluid(False,'Potassium Formate','','HYCOOL 50','','',THyCool_50);
  //TmjsSecCool.RegisterFluid(False,'Hydrofluoroether','','HFE-7100','','',THydroFluoroEther);
end.
