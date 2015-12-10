unit uMolChart;

interface

uses
  SysUtils, Windows, Graphics, Classes,
  superobject,
  cpIntf;

type
  fvec = array of Double;
  TMolChart = class
  private
    FCanvas: TCanvas;
    R: TRect;
    FMargin: Integer;
    x0, x1, y0, y1: Double;
    lny0, lny1: Double;
    xf, yf: Double;
    procedure CalcScaleFactor;
  public
    constructor Create(Canvas: TCanvas);
    function X2P(x: Double): Integer;
    function Y2P(y: Double): Integer;
    function P2X(p: Integer): Double;
    function P2Y(p: Integer): Double;
    procedure MoveTo(x, y: Double);
    procedure Lineto(x, y: Double);
    procedure DrawPoint(x, y: Double; Color: TColor);
    procedure DrawLine(x0, y0, x1, y1: Double; Color: TColor; PenWidth: Integer);
    procedure DrawVec(x, y: fvec; Color: TColor; PenWidth: Integer);
    procedure DrawChart(ref: string; Width, Height: Integer);
  end;

implementation

uses
  Math;

constructor TMolChart.Create(Canvas: TCanvas);
begin
  FMargin := 25;
  FCanvas := Canvas;
end;

procedure TMolChart.CalcScaleFactor;
{calculate the scale factor}
begin
  xf := (x1-x0)/(R.Right-R.Left);
  lny0 := ln(y0);
  lny1 := ln(y1);
  yf := (lny1-lny0)/(R.Top-R.Bottom);
end;

function TMolChart.X2P(x: Double): Integer;
{X scale value to chart pixel value}
begin
  Result := Trunc(R.Left+(x-x0)/xf);
end;

function TMolChart.Y2P(y: Double): Integer;
{Y scale value to chart pixel value}
begin
  Result := Trunc(R.Bottom+(ln(y)-lny0)/yf);
end;

function TMolChart.P2X(p: Integer): Double;
{pixel point to x scale value}
begin
  Result := x0+(p-R.Left)*xf
end;

function TMolChart.P2Y(p: Integer): Double;
{pixel point to y scale value}
begin
  Result := exp(lny0+(p-R.bottom)*yf);
end;

procedure TMolChart.MoveTo(x, y: Double);
{move to (x,y) scale position}
begin
  FCanvas.MoveTo(X2P(x), Y2P(y));
end;

procedure TMolChart.LineTo(x, y: Double);
{line to (x,y) scale position}
begin
  FCanvas.LineTo(X2P(x), Y2P(y));
end;

procedure TMolChart.DrawPoint(x, y: Double; Color: TColor);
{draw a square at (x,y) scale position}
var
  px, py: Integer;
begin
  FCanvas.Pen.Color := Color;
  px := X2P(x);
  py := Y2P(y);
  FCanvas.Rectangle(px-1, py+1, px+1, py-1);
end;

procedure TMolChart.DrawLine(x0, y0, x1, y1: Double; Color: TColor; PenWidth: Integer);
{draw a line from (x0,y0) to (x1,y1) scale co-ordinates}
begin
  FCanvas.Pen.Color := Color;
  FCanvas.Pen.Width := PenWidth;
  MoveTo(x0, y0);
  LineTo(x1, y1);
end;

procedure TMolChart.DrawVec(x, y: fvec; Color: TColor; PenWidth: Integer);
{draw a (xi,yi) vector}
var
  i: Integer;
begin
  FCanvas.Pen.Color := Color;
  FCanvas.Pen.Width := PenWidth;
  MoveTo(x[0], y[0]);
  for i := 1 to length(x)-1 do begin
    LineTo(x[i], y[i]);
    //DrawPoint(x[i], y[i], Color);
  end;
end;

procedure TMolChart.DrawChart(ref: string; Width, Height: Integer);
{draw the mollier chart}
const
  n=17;
var
  i: Integer;
  pc, tc, dc, hc: Double;
  p, hf, hg: fvec;
  dx, dy: Double;
begin
  R := Rect(0, 0, Width, Height);
  InflateRect(R, -FMargin, -FMargin);
  FCanvas.Pen.Color := clGray;
  FCanvas.Pen.Width := 1;
  FCanvas.Rectangle(R);
  pc := 1e-6*Props1SI('PCRIT', CoolStr(ref));
  tc := Props1SI('TCRIT', CoolStr(ref));
  dc := Props1SI('RHOCRIT', CoolStr(ref));
  try
    // CoolProp crashes here! (with some fluids)
    hc := 1e-3*PropsSI('H', 'T', tc, 'D', dc, CoolStr(ref));
  except
    Exit;
  end;
  y0 := 1;
  y1 := pc;
  if y1 <= 1 then begin
    // fluid critical pressure < 1 Mpa
    Exit;
  end;

  // generate n points
  SetLength(p, n);
  SetLength(hf, n);
  SetLength(hg, n);
  for i := 0 to n-1 do begin
    p[i] := y0 + i*(y1-y0)/n;
    hf[i] := 1e-3*PropsSI('H', 'P', 1e6*p[i], 'Q', 0, CoolStr(ref));
    hg[i] := 1e-3*PropsSI('H', 'P', 1e6*p[i], 'Q', 1, CoolStr(ref));
  end;

  // adjust the scale
  x0 := MinValue(hf);
  x1 := MaxValue(hg);
  dx := x1-x0;
  dy := y1-y0;
  x0 := x0-0.1*dx;
  x1 := x1+0.4*dx;
  y1 := y1+0.1*dy;
  CalcScaleFactor;

  // draw the liquid line
  DrawVec(hf, p, clGreen, 2);

  // draw the vapor line
  DrawVec(hg, p, clRed, 2);
  DrawPoint(hc, pc, clBlue);
end;

end.
