Attribute VB_Name = "Module1"
'Information on calling DLL from Excel:
'http://msdn.microsoft.com/en-us/library/office/bb687915.aspx
'
'Information on compiling DLL:
'http://msdn.microsoft.com/en-us/library/office/bb687850.aspx

Private Declare Function Props_private Lib "CoolProp_stdcall.dll" Alias "_Props@32" (ByVal Output As String, ByVal Name1 As Long, ByVal Value1 As Double, ByVal Name2 As Long, ByVal Value2 As Double, ByVal Ref As String) As Double
Private Declare Function Props1_private Lib "CoolProp_stdcall.dll" Alias "_Props1@8" (ByVal Ref As String, ByVal Output As String) As Double
Private Declare Function HAProps_private Lib "CoolProp_stdcall.dll" Alias "_HAProps@40" (ByVal Output As String, ByVal Input1Name As String, ByVal Value1 As Double, ByVal Input2Name As String, ByVal Value2 As Double, ByVal Input3name As String, ByVal Value3 As Double) As Double

Public Function Props(ByVal Output As String, ByVal Name1 As String, ByVal Value1 As Double, ByVal Name2 As String, ByVal Value2 As Double, ByVal Ref As String) As Double
    'Get the single character forms
    N1 = Asc(Left(Name1, 1))
    N2 = Asc(Left(Name2, 1))
    Props = Props_private(Output, N1, Value1, N2, Value2, Ref)
End Function

Public Function Props1(ByVal Ref As String, ByVal Output As String) As Double
    Props1 = Props1_private(Ref, Output)
End Function

Public Function HAProps(ByVal Output As String, ByVal Input1Name As String, ByVal Value1 As Double, ByVal Input2Name As String, ByVal Value2 As Double, ByVal Input3name As String, ByVal Value3 As Double) As Double
    HAProps = HAProps_private(Output, Input1Name, Value1, Input2Name, Value2, Input3name, Value3)
End Function


