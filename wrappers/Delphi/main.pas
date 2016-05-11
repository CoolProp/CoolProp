{
  CoolProp Graphic Demo - Main form
  Bruce Wernick, 13 December 2015
}

unit main;

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics,
  Controls, Forms, Dialogs, StdCtrls, ExtCtrls,
  uMolChart;

type
  TMainForm = class(TForm)
    ListBox1: TListBox;
    PaintBox1: TPaintBox;
    StaticText1: TStaticText;
    Button1: TButton;
    procedure FormActivate(Sender: TObject);
    procedure ListBox1Click(Sender: TObject);
    procedure FormCreate(Sender: TObject);
    procedure PaintBox1Paint(Sender: TObject);
    procedure PaintBox1MouseMove(Sender: TObject; Shift: TShiftState; X, Y: Integer);
    procedure PaintBox1MouseLeave(Sender: TObject);
    procedure Button1Click(Sender: TObject);
    procedure PaintBox1MouseDown(Sender: TObject; Button: TMouseButton; Shift: TShiftState; X, Y: Integer);
  private
    fluid: AnsiString;
    ref: PAnsiChar;
    MolChart: TMolChart;
  public
  end;

var
  MainForm: TMainForm;

implementation

{$R *.dfm}

uses
  clipbrd,
  cpIntf;

procedure TMainForm.FormCreate(Sender: TObject);
const
  n = 80;
var
  res: AnsiString;
begin
  Caption := 'CoolProp Demo';
  SetLength(res, n);
  get_global_param_string('version', PAnsiChar(res), n);
  Caption := Caption + ' ' + 'Version ' + string(res);
  get_global_param_string('gitrevision', PAnsiChar(res), n);
  Caption := Caption + ' ' + 'Revision ' + string(res);

  MolChart := TMolChart.Create(PaintBox1.Canvas);
  fluid := 'R22'; // initial fluid
  ref := PAnsiChar(fluid)
end;

procedure TMainForm.FormActivate(Sender: TObject);
const
  n = 2048;
var
  res: AnsiString;
begin
  SetLength(res, n);
  get_global_param_string('FluidsList', PAnsiChar(res), n);
  ListBox1.Items.CommaText := string(res);
  ListBox1.ItemIndex := ListBox1.Items.IndexOf(string(ref))
end;

procedure TMainForm.ListBox1Click(Sender: TObject);
begin
  fluid := AnsiString(ListBox1.Items[ListBox1.ItemIndex]);
  ref := PAnsiChar(fluid);
  PaintBox1.Refresh
end;

procedure TMainForm.PaintBox1MouseDown(Sender: TObject; Button: TMouseButton; Shift: TShiftState; X, Y: Integer);
var
  h, p: Double;
begin
  h := MolChart.P2X(X); p := MolChart.P2Y(Y);
  Clipboard.AsText := format('h:%0.3f, p:%0.3f', [h, p]);
end;

procedure TMainForm.PaintBox1MouseLeave(Sender: TObject);
begin
  StaticText1.Caption := string(fluid)
end;

procedure TMainForm.PaintBox1MouseMove(Sender: TObject; Shift: TShiftState; X, Y: Integer);
const
  k0 = 273.15;
  fmt = '%s: t=%0.2fC, p=%0.3f MPa, h=%0.3f kJ/kg, s=%0.4f kJ/kg-K';
var
  h, p, t, s: Double;
begin
  h := MolChart.P2X(X); p := MolChart.P2Y(Y);
  t := PropsSI('T', 'P', 1e6*p, 'H', 1e3*h, ref);
  s := 1e-3*PropsSI('S', 'P', 1e6*p, 'T', t+k0, ref);
  StaticText1.Caption := format(fmt, [fluid, t-k0, p, h, s])
end;

procedure TMainForm.PaintBox1Paint(Sender: TObject);
begin
  StaticText1.Caption := string(fluid);
  MolChart.DrawChart(ref, PaintBox1.Width, PaintBox1.Height)
end;

procedure TMainForm.Button1Click(Sender: TObject);
var
  MyFormat: Word;
  AData: THandle;
  APalette: HPALETTE;
  bmp: TBitmap;
begin
  bmp := self.GetFormImage;
  try
    bmp.SaveToClipBoardFormat(MyFormat, AData, APalette);
    ClipBoard.SetAsHandle(MyFormat,AData);
  finally
    bmp.Free;
  end;
end;

end.
