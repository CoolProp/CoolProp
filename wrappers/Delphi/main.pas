unit main;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics,
  Controls, Forms, Dialogs, StdCtrls, ExtCtrls,
  superobject, uMolChart;

type
  TMainForm = class(TForm)
    ListBox1: TListBox;
    PaintBox1: TPaintBox;
    StaticText1: TStaticText;
    procedure FormActivate(Sender: TObject);
    procedure ListBox1Click(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure PaintBox1Paint(Sender: TObject);
    procedure PaintBox1MouseMove(Sender: TObject; Shift: TShiftState; X,
      Y: Integer);
  private
    cp: TCanvas;
    ref: string;
    MolChart: TMolChart;
  public
  end;

var
  MainForm: TMainForm;

implementation

{$R *.dfm}

uses
  cpIntf;

procedure TMainForm.FormCreate(Sender: TObject);
var
  res, rev: CoolVec;
begin
  Caption := 'CoolProp Demo';
  res := '';
  rev := '';
  try
    // CoolProp crashes here!
    get_global_param_string('version', res);
    Caption := Caption + ' ' + 'Version ' + string(res);

    get_global_param_string('gitrevision', rev);
    Caption := Caption + ' ' + 'Revision ' + string(rev);
  except
  end;

  cp := PaintBox1.Canvas;
  MolChart := TMolChart.Create(cp);
  ref := 'R22';
end;

procedure TMainForm.FormActivate(Sender: TObject);
var
  dat: CoolVec;
begin
  get_global_param_string('FluidsList', dat);
  try
    // CoolProp crashes here!
    ListBox1.Items.CommaText := string(dat);
  except
  end;
  ListBox1.ItemIndex := ListBox1.Items.IndexOf(ref);
end;

procedure TMainForm.ListBox1Click(Sender: TObject);
begin
  ref := ListBox1.Items[ListBox1.ItemIndex];
  PaintBox1.Refresh;
end;

procedure TMainForm.PaintBox1MouseMove(Sender: TObject; Shift: TShiftState; X,
  Y: Integer);
var
  h, p, t: Double;
  i, j: Integer;
begin
  h := MolChart.P2X(X);
  p := MolChart.P2Y(Y);
  try
    // CoolProp crashes here! (at high pressure liquid)
    t := PropsSI('T', 'P', 1e6*p, 'H', 1e3*h, CoolStr(ref))-273.15;
  except
    t := -99;
    // paint a yellow blob at the position of the error
    for j := -1 to 1 do
      for i := -1 to 1 do
        PaintBox1.Canvas.Pixels[X+i,Y+j] := clYellow;
  end;
  StaticText1.Caption := format('h=%0.3f, p=%0.3f, t=%0.2f', [h,p,t]);
end;

procedure TMainForm.PaintBox1Paint(Sender: TObject);
begin
  MolChart.DrawChart(ref, PaintBox1.Width, PaintBox1.Height);
end;

end.
