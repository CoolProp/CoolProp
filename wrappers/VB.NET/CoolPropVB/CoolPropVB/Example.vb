Module Example

    Sub Main()
        ' High-level interface example
        Dim Tbp As Double
        Tbp = CoolProp.CoolProp.PropsSI("T", "P", 101325, "Q", 0, "Water")
        Console.Write(Tbp)
        Console.Write(vbCrLf)

        ' Low-level interface example
        Dim Water As CoolProp.AbstractState
        Water = CoolProp.AbstractState.factory("HEOS", "Water")
        Water.update(CoolProp.input_pairs.PQ_INPUTS, 101325, 0)
        Console.Write(Water.T())
        Console.Write(vbCrLf)

    End Sub

End Module
