Option Explicit On
Imports System.Math
Module coolprop
    Private Declare Function errstring_private Lib "C:\CoolProp\CoolProp.dll" Alias "_get_errstring_copy@4" (ByVal errstring As String) As Long
    Private Declare Function Props_private Lib "C:\CoolProp\CoolProp.dll" Alias "_Props@32" (ByVal Output As String, ByVal Name1 As Integer, ByVal Value1 As Double, ByVal Name2 As Integer, ByVal Value2 As Double, ByVal Ref As String) As Double
    Private Declare Function Props1_private Lib "C:\CoolProp\CoolProp.dll" Alias "_Props1@8" (ByVal Ref As String, ByVal Output As String) As Double
    Private Declare Function HAProps_private Lib "C:\CoolProp\CoolProp.dll" Alias "_HAProps@40" (ByVal Output As String, ByVal Input1Name As String, ByVal Value1 As Double, ByVal Input2Name As String, ByVal Value2 As Double, ByVal Input3name As String, ByVal Value3 As Double) As Double

    Public Function Props(ByVal Output As String, ByVal Name1 As String, ByVal Value1 As Double, ByVal Name2 As String, ByVal Value2 As Double, ByVal Ref As String) As Double
        Dim N1 As Integer
        Dim N2 As Integer
        Dim Props_temp As Double
        Dim errstring As String

        '1/ Get the single character forms
        N1 = Asc(Left(Name1, 1))
        N2 = Asc(Left(Name2, 1))

        '2/ Call the function coolprop
        Props_temp = Props_private(Output, N1, Value1, N2, Value2, Ref)

        If Abs(Props_temp) > 100000000 Then
            '3/ Make an empty string
            errstring = StrDup(1000, vbNullChar)
            '4/ Get the error
            errstring_private(errstring)
            MsgBox(errstring)
            Props = Props_temp '5/ A été ajouté afin de compiler 
        Else
            Props = Props_temp
        End If
    End Function

    Public Function Props1(ByVal Ref As String, ByVal Output As String) As Double
        Dim Props1_temp As Double
        Dim errstring As String

        '1/ Call the function coolprop
        Props1_temp = Props1_private(Ref, Output)

        If Abs(Props1_temp) > 100000000 Then
            '2/ Make an empty string
            errstring = StrDup(1000, vbNullChar)
            '3/ Get the error
            errstring_private(errstring)
            MsgBox(errstring)
            Props1 = Props1_temp '5/ A été ajouté afin de compiler 
        Else
            Props1 = Props1_temp
        End If
    End Function

    Public Function HAProps(ByVal Output As String, ByVal Input1Name As String, ByVal Value1 As Double, ByVal Input2Name As String, ByVal Value2 As Double, ByVal Input3name As String, ByVal Value3 As Double) As Double
        HAProps = HAProps_private(Output, Input1Name, Value1, Input2Name, Value2, Input3name, Value3)
    End Function

End Module
