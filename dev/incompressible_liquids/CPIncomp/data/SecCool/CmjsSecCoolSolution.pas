unit CmjsSecCoolSolution;

interface

uses
  SysUtils,CmjsSecCoolFluid;

type
  TFreezium =
    class(TmjsSecCoolFluid)
      private
      protected
        function  GetTFreezeX(X : Double) : Double; override;
        function  GetRhoTX(T,X : Double): Double; override;
        function  GetCpTX(T,X : Double): Double; override;
        function  GetCondTX(T,X : Double): Double; override;
        function  GetMuTX(T,X : Double): Double; override;
      public
        constructor Create; override;
    end;

implementation

uses
  Math,CmjsSecCool;

{ TFreezium }

constructor TFreezium.Create;
begin
  inherited;
  FluidKind      := fkSolution;
  ConcBase       := cobMass;
  EqConcBase     := cobMass;
  ConcBaseChange := False;
  FluidType      := 'Potassium Formate';
  TradeName      := 'Freezium';
  Reference      := 'Kemira Chemicals OY';
  TMin           := -40;//-59;
  TMax           := 40;
  XMin           := 19;
  XMax           := 50;
end;

function TFreezium.GetCondTX(T, X: Double): Double;
begin
  T      := T/100;
  X      := X/100;
  Result := 0.55-0.15*X+0.18*T-0.16*X*T;
end;

function TFreezium.GetCpTX(T, X: Double): Double;
begin
  T      := T/100;
  X      := X/100;
  Result := (4.15*exp(-0.9*X)+0.63*T*X)*1000;
end;

function TFreezium.GetMuTX(T, X: Double): Double;
var
  Tr,c : Double;
begin
  Tr     := T/100;
  c      := X/100;
  Result := 0.32+c*(-0.70+c*2.26)+Tr*(-1.26+Tr*(1.12-Tr*0.894));
  Result := GetRhoTX(T,X)*Power(10,Result)*1E-3;
end;

function TFreezium.GetRhoTX(T, X: Double): Double;
begin
  T      := T/100;
  X      := X/100;
  Result := (1.015+X*(0.462+X*0.406)-0.04*T)*1000;
end;

function TFreezium.GetTFreezeX(X: Double): Double;
const
  a = 0.03422039835160944;    b = -0.05425629002714395;
  c = -0.007991085555390726;  d = 0.001036937163846157;
  e = 0.0003268582531827402;  f = -7.721483884155584E-06;
  g = -4.841293026057464E-06; h = 1.216929917247388E-08;
  i = 2.41704396074865E-08;   j = 4.314861246570078E-11;
begin
  Result := (a+X*(c+X*(e+X*(g+i*X))))/
            (1+X*(b+X*(d+X*(f+X*(h+j*X)))));
  Result := Result*100;
  if Result < TMin then
    Result := TMin;
end;

initialization
  {RegisterFluid(IsFit : Boolean;AName,ACAS,ATradeName,AFileName,Formula : string;
          AFluidClass : TmjsSecCoolFluidClass);}
  TmjsSecCool.RegisterFluid(False,'Potassium Formate','','Freezium','','',TFreezium);

end.
