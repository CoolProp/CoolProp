object MainForm: TMainForm
  Left = 0
  Top = 0
  Caption = 'CoolProp Demo'
  ClientHeight = 372
  ClientWidth = 595
  Color = clBtnFace
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
    Top = 8
    Width = 412
    Height = 333
    Anchors = [akLeft, akTop, akRight, akBottom]
    OnMouseMove = PaintBox1MouseMove
    OnPaint = PaintBox1Paint
    ExplicitWidth = 361
    ExplicitHeight = 321
  end
  object ListBox1: TListBox
    Left = 8
    Top = 8
    Width = 161
    Height = 356
    Anchors = [akLeft, akTop, akBottom]
    ItemHeight = 13
    TabOrder = 0
    OnClick = ListBox1Click
    ExplicitHeight = 344
  end
  object StaticText1: TStaticText
    Left = 175
    Top = 347
    Width = 412
    Height = 17
    Anchors = [akLeft, akRight, akBottom]
    Caption = 'Chart position...'
    TabOrder = 1
  end
end
