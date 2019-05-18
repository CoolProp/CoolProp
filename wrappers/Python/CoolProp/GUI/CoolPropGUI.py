import wx
import wx.grid
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as WXCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx as WXToolbar
import matplotlib as mpl
import CoolProp as CP
from CoolProp.Plots.Plots import Ph, Ts
from CoolProp.Plots import PsychChart
import numpy as np

# Munge the system path if necessary to add the lib folder (only really needed
# for packaging using cx_Freeze)
# if os.path.exists('lib') and os.path.abspath(os.path.join(os.curdir,'lib')) not in os.:


class PlotPanel(wx.Panel):
    def __init__(self, parent, **kwargs):
        wx.Panel.__init__(self, parent, **kwargs)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.figure = mpl.figure.Figure(dpi=100)
        self.canvas = WXCanvas(self, -1, self.figure)
        self.ax = self.figure.add_axes((0.15, 0.15, 0.8, 0.8))
        #self.toolbar = WXToolbar(self.canvas)
        # self.toolbar.Realize()
        sizer.Add(self.canvas, 1, wx.EXPAND)
        # sizer.Add(self.toolbar)
        self.SetSizer(sizer)
        sizer.Layout()


class TSPlotFrame(wx.Frame):
    def __init__(self, Fluid):
        wx.Frame.__init__(self, None, title='T-s plot: ' + Fluid)

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.PP = PlotPanel(self, size=(-1, -1))

        sizer.Add(self.PP, 1, wx.EXPAND)
        self.SetSizer(sizer)
        Ts(str(Fluid),
           axis=self.PP.ax,
           Tmin=CP.CoolProp.Props(str(Fluid), 'Ttriple') + 0.01)
        sizer.Layout()

        self.add_menu()

    def add_menu(self):
        # Menu Bar
        self.MenuBar = wx.MenuBar()
        self.File = wx.Menu()

        mnuItem = wx.MenuItem(self.File, -1, "Edit...", "", wx.ITEM_NORMAL)

        self.File.AppendItem(mnuItem)
        self.MenuBar.Append(self.File, "File")

        self.SetMenuBar(self.MenuBar)


class PsychOptions(wx.Dialog):
    def __init__(self, parent):
        wx.Dialog.__init__(self, parent)

        self.build_contents()
        self.layout()

    def build_contents(self):
        self.p_label = wx.StaticText(self, label='Pressure [kPa (absolute)]')
        self.p = wx.TextCtrl(self, value='101.325')
        self.Tmin_label = wx.StaticText(self, label='Minimum dry bulb temperature [\xb0 C]')
        self.Tmin = wx.TextCtrl(self, value='-10')
        self.Tmax_label = wx.StaticText(self, label='Maximum dry bulb temperature [\xb0 C]')
        self.Tmax = wx.TextCtrl(self, value='60')
        self.GoButton = wx.Button(self, label='Accept')
        self.GoButton.Bind(wx.EVT_BUTTON, self.OnAccept)

    def OnAccept(self, event):
        self.EndModal(wx.ID_OK)

    def layout(self):
        sizer = wx.FlexGridSizer(cols=2)
        sizer.AddMany([self.p_label, self.p, self.Tmin_label, self.Tmin, self.Tmax_label, self.Tmax])
        sizer.Add(self.GoButton)
        sizer.Layout()
        self.Fit()


class PsychPlotFrame(wx.Frame):
    def __init__(self, Tmin=263.15, Tmax=333.15, p=101.325, **kwargs):

        wx.Frame.__init__(self, None, title='Psychrometric plot', **kwargs)

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.PP = PlotPanel(self)

        self.PP.figure.delaxes(self.PP.ax)
        self.PP.ax = self.PP.figure.add_axes((0.1, 0.1, 0.85, 0.85))

        sizer.Add(self.PP, 1, wx.EXPAND)
        self.SetSizer(sizer)

        PsychChart.p = p
        PsychChart.Tdb = np.linspace(Tmin, Tmax)

        SL = PsychChart.SaturationLine()
        SL.plot(self.PP.ax)

        RHL = PsychChart.HumidityLines([0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
        RHL.plot(self.PP.ax)

        HL = PsychChart.EnthalpyLines(range(-20, 100, 10))
        HL.plot(self.PP.ax)

        PF = PsychChart.PlotFormatting()
        PF.plot(self.PP.ax)

        sizer.Layout()

        self.add_menu()

        self.PP.toolbar = WXToolbar(self.PP.canvas)
        self.PP.toolbar.Realize()
        self.PP.GetSizer().Add(self.PP.toolbar)

        self.PP.Layout()

    def add_menu(self):
        # Menu Bar
        self.MenuBar = wx.MenuBar()
        self.File = wx.Menu()

        mnuItem = wx.MenuItem(self.File, -1, "Edit...", "", wx.ITEM_NORMAL)

        self.File.AppendItem(mnuItem)
        self.MenuBar.Append(self.File, "File")

        self.SetMenuBar(self.MenuBar)


class PHPlotFrame(wx.Frame):
    def __init__(self, Fluid):
        wx.Frame.__init__(self, None, title='p-h plot: ' + Fluid)

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.PP = PlotPanel(self, size=(-1, -1))

        sizer.Add(self.PP, 1, wx.EXPAND)
        self.SetSizer(sizer)
        Ph(str(Fluid),
           axis=self.PP.ax,
           Tmin=CP.CoolProp.Props(str(Fluid), 'Ttriple') + 0.01)
        sizer.Layout()

        self.add_menu()

    def add_menu(self):
        # Menu Bar
        self.MenuBar = wx.MenuBar()
        self.File = wx.Menu()

        mnuItem = wx.MenuItem(self.File, -1, "Edit...", "", wx.ITEM_NORMAL)

        self.File.AppendItem(mnuItem)
        self.MenuBar.Append(self.File, "File")

        self.SetMenuBar(self.MenuBar)

    def overlay_points(self):
        pass

    def overlay_cycle(self):
        pass


class SimpleGrid(wx.grid.Grid):
    def __init__(self, parent, ncol=20, nrow=8):
        wx.grid.Grid.__init__(self, parent)

        self.CreateGrid(ncol, nrow)
        [self.SetCellValue(i, j, '0.0') for i in range(20) for j in range(8)]


class SaturationTableDialog(wx.Dialog):
    def __init__(self, parent):
        wx.Dialog.__init__(self, parent)

        self.FluidLabel = wx.StaticText(self, label="Fluid")
        self.FluidCombo = wx.ComboBox(self)
        self.FluidCombo.AppendItems(sorted(CP.__fluids__))
        self.FluidCombo.SetEditable(False)
        self.TtripleLabel = wx.StaticText(self, label="Critical Temperature [K]")
        self.TtripleValue = wx.TextCtrl(self)
        self.TtripleValue.Enable(False)
        self.TcritLabel = wx.StaticText(self, label="Critical Temperature [K]")
        self.TcritValue = wx.TextCtrl(self)
        self.TcritValue.Enable(False)
        self.NvalsLabel = wx.StaticText(self, label="Number of values")
        self.NvalsValue = wx.TextCtrl(self)
        self.TminLabel = wx.StaticText(self, label="Minimum Temperature [K]")
        self.TminValue = wx.TextCtrl(self)
        self.TmaxLabel = wx.StaticText(self, label="Maximum Temperature [K]")
        self.TmaxValue = wx.TextCtrl(self)

        self.Accept = wx.Button(self, label="Accept")

        sizer = wx.FlexGridSizer(cols=2)
        sizer.AddMany([self.FluidLabel, self.FluidCombo,
                       self.TtripleLabel, self.TtripleValue,
                       self.TcritLabel, self.TcritValue])
        sizer.AddSpacer(10)
        sizer.AddSpacer(10)
        sizer.AddMany([self.NvalsLabel, self.NvalsValue,
                       self.TminLabel, self.TminValue,
                       self.TmaxLabel, self.TmaxValue])
        sizer.Add(self.Accept)

        self.Bind(wx.EVT_COMBOBOX, self.OnSelectFluid)
        self.Bind(wx.EVT_BUTTON, self.OnAccept)

        self.SetSizer(sizer)
        sizer.Layout()
        self.Fit()

        # Bind a key-press event to all objects to get Esc
        children = self.GetChildren()
        for child in children:
            child.Bind(wx.EVT_KEY_UP, self.OnKeyPress)

    def OnKeyPress(self, event=None):
        """ cancel if Escape key is pressed """
        event.Skip()
        if event.GetKeyCode() == wx.WXK_ESCAPE:
            self.EndModal(wx.ID_CANCEL)

    def get_values(self):
        Fluid = str(self.FluidCombo.GetStringSelection())
        if Fluid:
            N = float(self.NvalsValue.GetValue())
            Tmin = float(self.TminValue.GetValue())
            Tmax = float(self.TmaxValue.GetValue())
            Tvals = np.linspace(Tmin, Tmax, N)
            return Fluid, Tvals
        else:
            return '', []

    def OnCheckTmin(self):
        pass

    def OnCheckTmax(self):
        pass

    def OnAccept(self, event=None):
        self.EndModal(wx.ID_OK)

    def OnSelectFluid(self, event=None):
        Fluid = str(self.FluidCombo.GetStringSelection())
        if Fluid:
            Tcrit = CP.CoolProp.Props(Fluid, 'Tcrit')
            Ttriple = CP.CoolProp.Props(Fluid, 'Ttriple')
            self.TcritValue.SetValue(str(Tcrit))
            self.TtripleValue.SetValue(str(Ttriple))
            self.NvalsValue.SetValue('100')
            self.TminValue.SetValue(str(Ttriple + 0.01))
            self.TmaxValue.SetValue(str(Tcrit - 0.01))


class SaturationTable(wx.Frame):
    def __init__(self, parent):
        wx.Frame.__init__(self, parent)
        self.Fluid, self.Tvals = self.OnSelect()
        if self.Fluid:
            self.tbl = SimpleGrid(self,
                                  ncol=len(self.Tvals)
                                  )
            sizer = wx.BoxSizer(wx.VERTICAL)
            sizer.Add(self.tbl, 1, wx.EXPAND)
            self.SetSizer(sizer)
            sizer.Layout()
            self.build()
            self.add_menu()
        else:
            self.Destroy()

    def OnSelect(self, event=None):
        dlg = SaturationTableDialog(None)
        if dlg.ShowModal() == wx.ID_OK:
            Fluid, Tvals = dlg.get_values()
            cancel = False
        else:
            cancel = True
        dlg.Destroy()
        if not cancel:
            return Fluid, Tvals
        else:
            return None, None

    def build(self):
        self.SetTitle('Saturation Table: ' + self.Fluid)
        self.tbl.SetColLabelValue(0, "Temperature\n[K]")
        self.tbl.SetColLabelValue(1, "Liquid Pressure\n[kPa]")
        self.tbl.SetColLabelValue(2, "Vapor Pressure\n[kPa]")
        self.tbl.SetColLabelValue(3, "Liquid Density\n[kg/m3]")
        self.tbl.SetColLabelValue(4, "Vapor Density\n[kg/m3]")

        for i, T in enumerate(self.Tvals):
            Fluid = self.Fluid
            pL = CP.CoolProp.Props('P', 'T', T, 'Q', 0, Fluid)
            pV = CP.CoolProp.Props('P', 'T', T, 'Q', 1, Fluid)
            rhoL = CP.CoolProp.Props('D', 'T', T, 'Q', 0, Fluid)
            rhoV = CP.CoolProp.Props('D', 'T', T, 'Q', 1, Fluid)

            self.tbl.SetCellValue(i, 0, str(T))
            self.tbl.SetCellValue(i, 1, str(pL))
            self.tbl.SetCellValue(i, 2, str(pV))
            self.tbl.SetCellValue(i, 3, str(rhoL))
            self.tbl.SetCellValue(i, 4, str(rhoV))

    def add_menu(self):
        # Menu Bar
        self.MenuBar = wx.MenuBar()
        self.File = wx.Menu()

        mnuItem0 = wx.MenuItem(self.File, -1, "Select All \tCtrl+A", "", wx.ITEM_NORMAL)
        mnuItem1 = wx.MenuItem(self.File, -1, "Copy selected data \tCtrl+C", "", wx.ITEM_NORMAL)
        mnuItem2 = wx.MenuItem(self.File, -1, "Copy table w/ headers \tCtrl+H", "", wx.ITEM_NORMAL)

        self.File.AppendItem(mnuItem0)
        self.File.AppendItem(mnuItem1)
        self.File.AppendItem(mnuItem2)
        self.MenuBar.Append(self.File, "Edit")
        self.Bind(wx.EVT_MENU, lambda event: self.tbl.SelectAll(), mnuItem0)
        self.Bind(wx.EVT_MENU, self.OnCopy, mnuItem1)
        self.Bind(wx.EVT_MENU, self.OnCopyHeaders, mnuItem2)

        self.SetMenuBar(self.MenuBar)

    def OnCopy(self, event=None):

        # Number of rows and cols
        rows = self.tbl.GetSelectionBlockBottomRight()[0][0] - self.tbl.GetSelectionBlockTopLeft()[0][0] + 1
        cols = self.tbl.GetSelectionBlockBottomRight()[0][1] - self.tbl.GetSelectionBlockTopLeft()[0][1] + 1

        # data variable contain text that must be set in the clipboard
        data = ''

        # For each cell in selected range append the cell value in the data variable
        # Tabs '\t' for cols and '\r' for rows
        for r in range(rows):
            for c in range(cols):
                data = data + str(self.tbl.GetCellValue(self.tbl.GetSelectionBlockTopLeft()[0][0] + r, self.tbl.GetSelectionBlockTopLeft()[0][1] + c))
                if c < cols - 1:
                    data = data + '\t'
            data = data + '\n'
        # Create text data object
        clipboard = wx.TextDataObject()
        # Set data object value
        clipboard.SetText(data)
        # Put the data in the clipboard
        if wx.TheClipboard.Open():
            wx.TheClipboard.SetData(clipboard)
            wx.TheClipboard.Close()
        else:
            wx.MessageBox("Can't open the clipboard", "Error")
        event.Skip()

    def OnCopyHeaders(self, event=None):
        self.tbl.SelectAll()
        # Number of rows and cols
        rows = self.tbl.GetSelectionBlockBottomRight()[0][0] - self.tbl.GetSelectionBlockTopLeft()[0][0] + 1
        cols = self.tbl.GetSelectionBlockBottomRight()[0][1] - self.tbl.GetSelectionBlockTopLeft()[0][1] + 1

        # data variable contain text that must be set in the clipboard
        data = ''

        # Add the headers
        for c in range(cols):
            data += str(self.tbl.GetColLabelValue(c).replace('\n', ' '))
            if c < cols - 1:
                data += '\t'
        data = data + '\n'
        # For each cell in selected range append the cell value in the data variable
        # Tabs '\t' for cols and '\r' for rows
        for r in range(rows):
            for c in range(cols):
                data = data + str(self.tbl.GetCellValue(self.tbl.GetSelectionBlockTopLeft()[0][0] + r, self.tbl.GetSelectionBlockTopLeft()[0][1] + c))
                if c < cols - 1:
                    data = data + '\t'
            data = data + '\n'
        # Create text data object
        clipboard = wx.TextDataObject()
        # Set data object value
        clipboard.SetText(data)
        # Put the data in the clipboard
        if wx.TheClipboard.Open():
            wx.TheClipboard.SetData(clipboard)
            wx.TheClipboard.Close()
        else:
            wx.MessageBox("Can't open the clipboard", "Error")
        event.Skip()


class MainFrame(wx.Frame):

    def __init__(self):
        wx.Frame.__init__(self, None)

        self.build()

    def build(self):
        # Menu Bar
        self.MenuBar = wx.MenuBar()

        self.plots = wx.Menu()
        self.PHPlot = wx.Menu()
        self.TSPlot = wx.Menu()
        self.tables = wx.Menu()
        self.PsychPlot = wx.MenuItem(self.plots, -1, 'Psychrometric Plot')
        self.SatTable = wx.MenuItem(self.tables, -1, ' Saturation Table', "", wx.ITEM_NORMAL)

        for Fluid in sorted(CP.__fluids__):
            mnuItem = wx.MenuItem(self.PHPlot, -1, Fluid, "", wx.ITEM_NORMAL)
            self.PHPlot.AppendItem(mnuItem)
            self.Bind(wx.EVT_MENU, lambda event: self.OnPHPlot(event, mnuItem), mnuItem)

            mnuItem = wx.MenuItem(self.TSPlot, -1, Fluid, "", wx.ITEM_NORMAL)
            self.TSPlot.AppendItem(mnuItem)
            self.Bind(wx.EVT_MENU, lambda event: self.OnTSPlot(event, mnuItem), mnuItem)

        self.MenuBar.Append(self.plots, "Plots")
        self.plots.AppendItem(self.PsychPlot)
        self.plots.AppendMenu(-1, 'p-h plot', self.PHPlot)
        self.plots.AppendMenu(-1, 'T-s plot', self.TSPlot)
        self.MenuBar.Append(self.tables, "Tables")
        self.tables.AppendItem(self.SatTable)
        self.Bind(wx.EVT_MENU, self.OnSatTable, self.SatTable)
        self.Bind(wx.EVT_MENU, self.OnPsychPlot, self.PsychPlot)

        self.SetMenuBar(self.MenuBar)

    def OnPsychPlot(self, event=None):

        # Load the options
        dlg = PsychOptions(None)
        if dlg.ShowModal() == wx.ID_OK:
            Tmin = float(dlg.Tmin.GetValue()) + 273.15
            Tmax = float(dlg.Tmax.GetValue()) + 273.15
            p = float(dlg.p.GetValue())
            PPF = PsychPlotFrame(Tmin=Tmin, Tmax=Tmax, p=p, size=(1000, 700))
            PPF.Show()
        dlg.Destroy()

    def OnSatTable(self, event):
        TBL = SaturationTable(None)
        TBL.Show()

    def OnPHPlot(self, event, mnuItem):
        # Make a p-h plot instance in a new frame
        # Get the label (Fluid name)
        Fluid = self.PHPlot.FindItemById(event.Id).Label
        PH = PHPlotFrame(Fluid)
        PH.Show()

    def OnTSPlot(self, event, mnuItem):
        # Make a p-h plot instance in a new frame
        # Get the label (Fluid name)
        Fluid = self.TSPlot.FindItemById(event.Id).Label
        TS = TSPlotFrame(Fluid)
        TS.Show()


if __name__ == '__main__':
    app = wx.App(False)
    wx.InitAllImageHandlers()

    frame = MainFrame()
    frame.Show(True)
    app.MainLoop()
