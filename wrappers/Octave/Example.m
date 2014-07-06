% Example of CoolProp for Octave
% Ian Bell, 2013

CoolProp
disp(['CoolProp version: ', CoolProp.get_global_param_string('version')])
disp(['CoolProp gitrevision: ', CoolProp.get_global_param_string('gitrevision')])
disp(['CoolProp fluids: ', CoolProp.get_global_param_string('FluidsList')])

disp(' ')
disp('************ USING EOS *************')
disp(' ')
disp('FLUID STATE INDEPENDENT INPUTS')
disp(['Critical Density Propane: ', num2str(CoolProp.Props1('Propane','rhocrit')), ' kg/m^3'])
disp(['TWO PHASE INPUTS (Pressure)'])
disp(['Density of saturated liquid Propane at 101.325 kPa: ', num2str(CoolProp.Props('D','P',101.325,'Q',0,'Propane')), ' kg/m^3'])
disp(['Density of saturated vapor R290 at 101.325 kPa: ', num2str(CoolProp.Props('D','P',101.325,'Q',1,'R290')), ' kg/m^3'])
disp(['TWO PHASE INPUTS (Temperature)'])
disp(['Density of saturated liquid Propane at 300 K: ', num2str(CoolProp.Props('D','T',300,'Q',0,'Propane')), ' kg/m^3'])
disp(['Density of saturated vapor R290 at 300 K: ', num2str(CoolProp.Props('D','T',300,'Q',1,'R290')), ' kg/m^3'])
disp(['SINGLE PHASE CYCLE (propane)'])
p = CoolProp.Props('P','T',300,'D',1,'Propane'); 
h = CoolProp.Props('H','T',300,'D',1,'Propane');
disp(['T,D -> P,H ', num2str(300),',',num2str(1), ' --> ',num2str(p),',',num2str(h)])
T = CoolProp.Props('T','P',p,'H',h,'Propane'); 
D = CoolProp.Props('D','P',p,'H',h,'Propane');
disp(['P,H -> T,D', num2str(p),',',num2str(h),'-->',num2str(T),',',num2str(D)])
 
disp([' '])
disp(['************ USING TTSE ***************'])
disp([' '])
CoolProp.enable_TTSE_LUT('Propane');
disp(['TWO PHASE INPUTS (Pressure)'])
disp(['Density of saturated liquid Propane at 101.325 kPa: ', num2str(CoolProp.Props('D','P',101.325,'Q',0,'Propane')), ' kg/m^3'])
disp(['Density of saturated vapor R290 at 101.325 kPa: ', num2str(CoolProp.Props('D','P',101.325,'Q',1,'R290')), ' kg/m^3'])
disp(['TWO PHASE INPUTS (Temperature)'])
disp(['Density of saturated liquid Propane at 300 K: ', num2str(CoolProp.Props('D','T',300,'Q',0,'Propane')), ' kg/m^3'])
disp(['Density of saturated vapor R290 at 300 K: ', num2str(CoolProp.Props('D','T',300,'Q',1,'R290')), ' kg/m^3'])
disp(['SINGLE PHASE CYCLE (propane)'])
p = CoolProp.Props('P','T',300,'D',1,'Propane'); 
h = CoolProp.Props('H','T',300,'D',1,'Propane');
disp(['T,D -> P,H ', num2str(300),',',num2str(1), ' --> ',num2str(p),',',num2str(h)])
T = CoolProp.Props('T','P',p,'H',h,'Propane'); 
D = CoolProp.Props('D','P',p,'H',h,'Propane');
disp(['P,H -> T,D ', num2str(p),',',num2str(h),' --> ',num2str(T),',',num2str(D)])
CoolProp.disable_TTSE_LUT('Propane');

try
    disp(' ')
    disp('************ USING REFPROP ***************')
    disp(' ')
    disp('FLUID STATE INDEPENDENT INPUTS')
    disp(['Critical Density Propane:', num2str(CoolProp.Props('REFPROP-Propane','rhocrit')), 'kg/m^3'])
    disp(['TWO PHASE INPUTS (Pressure)'])
    disp(['Density of saturated liquid Propane at 101.325 kPa: ', num2str(CoolProp.Props('D','P',101.325,'Q',0,'REFPROP-Propane')), ' kg/m^3'])
    disp(['Density of saturated vapor R290 at 101.325 kPa: ', num2str(CoolProp.Props('D','P',101.325,'Q',1,'REFPROP-Propane')), ' kg/m^3'])
    disp(['TWO PHASE INPUTS (Temperature)'])
    disp(['Density of saturated liquid Propane at 300 K: ', num2str(CoolProp.Props('D','T',300,'Q',0,'REFPROP-Propane')), ' kg/m^3'])
    disp(['Density of saturated vapor R290 at 300 K: ', num2str(CoolProp.Props('D','T',300,'Q',1,'REFPROP-Propane')), ' kg/m^3'])
    disp(['SINGLE PHASE CYCLE (propane)'])
    p = CoolProp.Props('P','T',300,'D',1,'REFPROP-Propane');
    h = CoolProp.Props('H','T',300,'D',1,'REFPROP-Propane');
    disp(['T,D -> P,H ', num2str(300),',',num2str(1), ' --> ',num2str(p),',',num2str(h)])
    T = CoolProp.Props('T','P',p,'H',h,'REFPROP-Propane'); 
    D = CoolProp.Props('D','P',p,'H',h,'REFPROP-Propane');
    disp(['P,H -> T,D ', num2str(p),',',num2str(h),' --> ',num2str(T),',',num2str(D)])
catch
    disp(' ')
    disp('************ CANT USE REFPROP ************')
    disp(' ')
end
 
 disp([' '])
 disp('************ CHANGE UNIT SYSTEM (default is kSI) *************')
 disp(' ')
 CoolProp.set_standard_unit_system(CoolProp.UNIT_SYSTEM_SI)
 disp(['Vapor pressure of water at 373.15 K in SI units (Pa):',num2str(CoolProp.Props('P','T',373.15,'Q',0,'Water'))]);
 CoolProp.set_standard_unit_system(CoolProp.UNIT_SYSTEM_KSI)
 disp(['Vapor pressure of water at 373.15 K in kSI units (kPa):',num2str(CoolProp.Props('P','T',373.15,'Q',0,'Water'))]);
 
disp(' ')
disp('************ BRINES AND SECONDARY WORKING FLUIDS *************')
disp(' ')
disp(['Density of 50% (mass) ethylene glycol/water at 300 K, 101.325 kPa: ', num2str(CoolProp.Props('D','T',300,'P',101.325,'EG-50%')), 'kg/m^3'])
disp(['Viscosity of Therminol D12 at 350 K, 101.325 kPa: ', num2str(CoolProp.Props('V', 'T', 350, 'P', 101.325, 'TD12')), 'Pa-s'])

disp(' ')
disp('************ HUMID AIR PROPERTIES *************')
disp(' ')
disp(['Humidity ratio of 50% rel. hum. air at 300 K, 101.325 kPa: ', num2str(CoolProp.HAProps('W','T',300,'P',101.325,'R',0.5)), ' kg_w/kg_da'])
disp(['Relative humidity from last calculation: ', num2str(CoolProp.HAProps('R','T',300,'P',101.325,'W',HAProps('W','T',300,'P',101.325,'R',0.5))), '(fractional)'])
