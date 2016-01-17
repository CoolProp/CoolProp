object MainForm: TMainForm
  Left = 0
  Top = 0
  Caption = 'CoolProp Demo'
  ClientHeight = 372
  ClientWidth = 595
  Color = clBtnFace
  DoubleBuffered = True
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  OnActivate = FormActivate
  OnCreate = FormCreate
  DesignSize = (
    595
    372)
  PixelsPerInch = 96
  TextHeight = 13
  object PaintBox1: TPaintBox
    Left = 175
    Top = 40
    Width = 412
    Height = 301
    Anchors = [akLeft, akTop, akRight, akBottom]
    OnMouseDown = PaintBox1MouseDown
    OnMouseLeave = PaintBox1MouseLeave
    OnMouseMove = PaintBox1MouseMove
    OnPaint = PaintBox1Paint
  end
  object ListBox1: TListBox
    Left = 8
    Top = 40
    Width = 161
    Height = 324
    Anchors = [akLeft, akTop, akBottom]
    ItemHeight = 13
    TabOrder = 0
    OnClick = ListBox1Click
  end
  object StaticText1: TStaticText
    Left = 175
    Top = 347
    Width = 412
    Height = 17
    Anchors = [akLeft, akRight, akBottom]
    AutoSize = False
    BorderStyle = sbsSunken
    Caption = 'Chart position...'
    TabOrder = 1
  end
  object Button1: TButton
    Left = 8
    Top = 8
    Width = 49
    Height = 25
    Caption = '&Copy'
    TabOrder = 2
    OnClick = Button1Click
  end
end
