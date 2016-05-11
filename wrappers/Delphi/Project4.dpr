{
  Project 4: CoolProp Graphic Demo
  Bruce Wernick, 13 December 2015
}

program Project4;

uses
  Forms,
  main in 'main.pas' {MainForm},
  uMolChart in 'uMolChart.pas',
  cpIntf in 'cpIntf.pas';

{$R *.res}

begin
  Application.Initialize;
  Application.MainFormOnTaskbar := True;
  Application.CreateForm(TMainForm, MainForm);
  Application.Run;
end.
