Class MainWindow

    Private Sub Button_Click(sender As Object, e As RoutedEventArgs)
        Dim resultat As Double

        resultat = Props1("R134a", "Tcrit")
        MsgBox(resultat)

        resultat = HAProps("H", "T", 40 + 273.15, "R", 0.4, "P", 101.325)
        MsgBox(resultat)

        ' ne fonctionne pas :
        resultat = Props("P", "T", 273.15 + 2, "Q", 0, "R407C")
        MsgBox(resultat)

    End Sub
End Class
