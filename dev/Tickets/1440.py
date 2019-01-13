from CoolProp.Plots import PropertyPlot

#plot = PropertyPlot("HEOS::ISOBUTAN[0.8]&PROPANE[0.2]", 'PH', unit_system='EUR', tp_limits='ACHP')
# plot.calc_isolines()
# plot.show()

#plot = PropertyPlot("HEOS::ISOBUTAN[0.8]&PROPANE[0.2]", 'TS', unit_system='EUR', tp_limits='ORC')
# plot.calc_isolines()
# plot.show()

plot = PropertyPlot("REFPROP::ISOBUTAN[0.8]&PROPANE[0.2]", 'PH', unit_system='EUR', tp_limits='ACHP')
plot.calc_isolines()
plot.show()

#plot = PropertyPlot("REFPROP::ISOBUTAN[0.8]&PROPANE[0.2]", 'TS', unit_system='EUR', tp_limits='ORC')
# plot.calc_isolines()
# plot.show()

plot = PropertyPlot("HEOS::R134a", 'PH', unit_system='EUR', tp_limits='ACHP')
plot.calc_isolines()
plot.show()
