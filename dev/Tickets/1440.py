from CoolProp.Plots import PropertyPlot
plot = PropertyPlot("REFPROP::ISOBUTAN[0.8]&PROPANE[0.2]", 'PH', unit_system='EUR', tp_limits='ACHP')
plot.calc_isolines()
plot.show()
