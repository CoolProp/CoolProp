disp 'Critical temperature of R410A'
Tc = Props('R410A','Tcrit')
disp 'should yield 344.494'

disp 'Density of carbon dioxide (R744) at 100 bar and 25C'
rho = Props('D','T',298.15,'P',10000,'R744')
disp 'should yield 817.6342521183478'

disp 'Saturated vapor enthalpy [kJ/kg] of R134a at STP'
rho = Props('H','T',298.15,'Q',1,'R134a')
disp 'should yield 412.333953231868'

disp 'Enthalpy (kJ per kg dry air) as a function of temperature, pressure, and 50% relative humidity at STP'
h = HAProps('H','T',298.15,'P',101.325,'R',0.5)
disp 'should yield 50.4323691981'

disp('Temperature of saturated air at the previous enthalpy')
T = HAProps('T','P',101.325,'H',h,'R',1.0)
disp 'should yield 290.962139387'