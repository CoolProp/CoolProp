from __future__ import division, print_function
import numpy as np
from .DataObjects import PureData, SolutionData, DigitalData,\
    CoefficientData


class PureExample(PureData):
    def __init__(self):
        PureData.__init__(self)
        self.name = "ExamplePure"
        self.description = "Heat transfer fluid TherminolD12 by Solutia"
        self.reference = "Solutia data sheet"
        self.Tmax = 150 + 273.15
        self.Tmin = 50 + 273.15
        self.TminPsat = self.Tmax

        self.density.source = self.density.SOURCE_DATA
        self.specific_heat.source = self.specific_heat.SOURCE_DATA
        self.conductivity.source = self.conductivity.SOURCE_DATA
        self.viscosity.source = self.viscosity.SOURCE_DATA
        self.saturation_pressure.source = self.saturation_pressure.SOURCE_DATA

        self.temperature.data = np.array([50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150]) + 273.15  # Kelvin
        self.density.data = np.array([740, 733, 726, 717, 710, 702, 695, 687, 679, 670, 662])        # kg/m3
        self.specific_heat.data = np.array([2235, 2280, 2326, 2361, 2406, 2445, 2485, 2528, 2571, 2607, 2645])        # J/kg-K
        self.viscosity.data = np.array([0.804, 0.704, 0.623, 0.556, 0.498, 0.451, 0.410, 0.374, 0.346, 0.317, 0.289])        # Pa-s
        self.conductivity.data = np.array([0.105, 0.104, 0.102, 0.100, 0.098, 0.096, 0.095, 0.093, 0.091, 0.089, 0.087])        # W/m-K
        self.saturation_pressure.data = np.array([0.5, 0.9, 1.4, 2.3, 3.9, 6.0, 8.7, 12.4, 17.6, 24.4, 33.2])        # Pa
        self.reshapeAll()


class SolutionExample(SolutionData):
    def __init__(self):
        SolutionData.__init__(self)
        self.name = "ExampleSolution"
        self.description = "Ethanol ice slurry"
        self.reference = "SecCool software,Skovrup2013"

        self.temperature.data = np.array([-45, -40, -35, -30, -25, -20, -15, -10]) + 273.15  # Kelvin
        self.concentration.data = np.array([5, 10, 15, 20, 25, 30, 35]) / 100.0  # mass fraction

        self.density.data = np.array([
          [1064.0, 1054.6, 1045.3, 1036.3, 1027.4, 1018.6, 1010.0],
          [1061.3, 1052.1, 1043.1, 1034.3, 1025.6, 1017.0, 1008.6],
          [1057.6, 1048.8, 1040.1, 1031.5, 1023.1, 1014.8, 1006.7],
          [1053.1, 1044.6, 1036.2, 1028.0, 1019.9, 1012.0, 1004.1],
          [1047.5, 1039.4, 1031.5, 1023.7, 1016.0, 1008.4, 1000.9],
          [1040.7, 1033.2, 1025.7, 1018.4, 1011.2, 1004.0, 997.0],
          [1032.3, 1025.3, 1018.5, 1011.7, 1005.1, 998.5, 992.0],
          [1021.5, 1015.3, 1009.2, 1003.1, 997.1, 991.2, 985.4]])  # kg/m3

        self.specific_heat.data = np.copy(self.density.data)

        self.density.source = self.density.SOURCE_DATA
        self.specific_heat.source = self.specific_heat.SOURCE_DATA

        self.Tmax = np.max(self.temperature.data)
        self.Tmin = np.min(self.temperature.data)
        self.xmax = np.max(self.concentration.data)
        self.xmin = np.min(self.concentration.data)
        self.xid = self.ifrac_mass
        self.TminPsat = self.Tmax


class DigitalExample(DigitalData):
    def __init__(self):
        DigitalData.__init__(self)

        self.name = "ExampleDigital"
        self.description = "some fluid"
        self.reference = "none"

        self.Tmin = 273.00;
        self.Tmax = 500.00;
        self.xmax = 1.0
        self.xmin = 0.0
        self.xid = self.ifrac_mass
        self.TminPsat = self.Tmin;

        self.temperature.data = self.getTrange()
        self.concentration.data = self.getxrange()

        def funcRho(T, x):
            return T + x * 100.0 + T * (x + 0.5)
        self.density.xData, self.density.yData, self.density.data = self.getArray(dataID="D", func=funcRho, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.density.DEBUG)
        self.density.source = self.density.SOURCE_EQUATION

        def funcCp(T, x):
            return T + x * 50.0 + T * (x + 0.6)
        self.specific_heat.xData, self.specific_heat.yData, self.specific_heat.data = self.getArray(dataID="C", func=funcCp, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.specific_heat.DEBUG)
        self.specific_heat.source = self.specific_heat.SOURCE_EQUATION


class DigitalExamplePure(PureData, DigitalData):
    def __init__(self):
        DigitalData.__init__(self)
        PureData.__init__(self)

        self.name = "ExampleDigitalPure"
        self.description = "water at 100 bar"
        self.reference = "none"

        self.Tmin = 280.00;
        self.Tmax = 500.00;

        self.TminPsat = self.Tmin;

        self.temperature.data = self.getTrange()
        self.concentration.data = self.getxrange()

        import CoolProp.CoolProp as CP

        def funcD(T, x):
            return CP.PropsSI('D', 'T', T, 'P', 1e7, 'water')

        def funcC(T, x):
            return CP.PropsSI('C', 'T', T, 'P', 1e7, 'water')

        def funcL(T, x):
            return CP.PropsSI('L', 'T', T, 'P', 1e7, 'water')

        def funcV(T, x):
            return CP.PropsSI('V', 'T', T, 'P', 1e7, 'water')

        def funcP(T, x):
            return CP.PropsSI('P', 'T', T, 'Q', 0.0, 'water')

        self.density.xData, self.density.yData, self.density.data = self.getArray(dataID="D", func=funcD, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.density.DEBUG)
        self.density.source = self.density.SOURCE_EQUATION

        self.specific_heat.xData, self.specific_heat.yData, self.specific_heat.data = self.getArray(dataID="C", func=funcC, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.specific_heat.DEBUG)
        self.specific_heat.source = self.specific_heat.SOURCE_EQUATION

        self.conductivity.xData, self.conductivity.yData, self.conductivity.data = self.getArray(dataID="L", func=funcL, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.conductivity.DEBUG)
        self.conductivity.source = self.conductivity.SOURCE_EQUATION

        self.viscosity.xData, self.viscosity.yData, self.viscosity.data = self.getArray(dataID="V", func=funcV, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.viscosity.DEBUG)
        self.viscosity.source = self.viscosity.SOURCE_EQUATION

        self.saturation_pressure.xData, self.saturation_pressure.yData, self.saturation_pressure.data = self.getArray(dataID="P", func=funcP, x_in=self.temperature.data, y_in=self.concentration.data, DEBUG=self.saturation_pressure.DEBUG)
        self.saturation_pressure.source = self.saturation_pressure.SOURCE_EQUATION


class SecCoolExample(CoefficientData):
    """
    Ethanol-Water mixture according to Melinder book
    Source: SecCool Software
    """

    def __init__(self):
        CoefficientData.__init__(self)
        self.name = "ExampleSecCool"
        self.description = "Methanol solution"
        #self.reference = "SecCool software"
        self.Tmax = 20 + 273.15
        self.Tmin = -50 + 273.15
        self.xmax = 0.5
        self.xmin = 0.0
        self.xid = self.ifrac_mass
        self.TminPsat = 20 + 273.15

        self.Tbase = -4.48 + 273.15
        self.xbase = 31.57 / 100.0

        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = self.convertSecCoolArray(np.array([
           960.24665800,
          -1.2903839100,
          -0.0161042520,
          -0.0001969888,
           1.131559E-05,
           9.181999E-08,
          -0.4020348270,
          -0.0162463989,
           0.0001623301,
           4.367343E-06,
           1.199000E-08,
          -0.0025204776,
           0.0001101514,
          -2.320217E-07,
           7.794999E-08,
           9.937483E-06,
          -1.346886E-06,
           4.141999E-08]))

        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = self.convertSecCoolArray(np.array([
           3822.9712300,
          -23.122409500,
           0.0678775826,
           0.0022413893,
          -0.0003045332,
          -4.758000E-06,
           2.3501449500,
           0.1788839410,
           0.0006828000,
           0.0002101166,
          -9.812000E-06,
          -0.0004724176,
          -0.0003317949,
           0.0001002032,
          -5.306000E-06,
           4.242194E-05,
           2.347190E-05,
          -1.894000E-06]))

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = self.convertSecCoolArray(np.array([
           0.4082066700,
          -0.0039816870,
           1.583368E-05,
          -3.552049E-07,
          -9.884176E-10,
           4.460000E-10,
           0.0006629321,
          -2.686475E-05,
           9.039150E-07,
          -2.128257E-08,
          -5.562000E-10,
           3.685975E-07,
           7.188416E-08,
          -1.041773E-08,
           2.278001E-10,
           4.703395E-08,
           7.612361E-11,
          -2.734000E-10]))

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
        self.viscosity.coeffs = self.convertSecCoolArray(np.array([
           1.4725525500,
           0.0022218998,
          -0.0004406139,
           6.047984E-06,
          -1.954730E-07,
          -2.372000E-09,
          -0.0411841566,
           0.0001784479,
          -3.564413E-06,
           4.064671E-08,
           1.915000E-08,
           0.0002572862,
          -9.226343E-07,
          -2.178577E-08,
          -9.529999E-10,
          -1.699844E-06,
          -1.023552E-07,
           4.482000E-09]))

        self.T_freeze.type = self.T_freeze.INCOMPRESSIBLE_POLYOFFSET
        self.T_freeze.coeffs = np.array([
           27.755555600 / 100.0,
          -22.973221700 + 273.15,
          -1.1040507200 * 100.0,
          -0.0120762281 * 100.0 * 100.0,
          -9.343458E-05 * 100.0 * 100.0 * 100.0])

        self.density.source = self.density.SOURCE_COEFFS
        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.viscosity.source = self.viscosity.SOURCE_COEFFS
        self.T_freeze.source = self.T_freeze.SOURCE_COEFFS


class MelinderExample(CoefficientData):
    """
    Methanol-Water mixture according to Melinder book
    Source: Book
    """

    def __init__(self):
        CoefficientData.__init__(self)
        self.name = "ExampleMelinder"
        self.description = "Methanol solution"
        self.reference = "Melinder2010"
        self.Tmax = 40 + 273.15
        self.Tmin = -50 + 273.15
        self.xmax = 0.6
        self.xmin = 0.0
        self.xid = self.ifrac_mass
        self.TminPsat = self.Tmax

        self.Tbase = 3.5359 + 273.15;
        self.xbase = 30.5128 / 100.0

        coeffs = np.array([
        [-26.29, 958.1, 3887, 0.4175, 1.153],
        [-0.000002575, -0.4151, 7.201, 0.0007271, -0.03866],
        [-0.000006732, -0.002261, -0.08979, 0.0000002823, 0.0002779],
        [0.000000163, 0.0000002998, -0.000439, 0.000000009718, -0.000001543],
        [-1.187, -1.391, -18.5, -0.004421, 0.005448],
        [-0.00001609, -0.0151, 0.2984, -0.00002952, 0.0001008],
        [0.000000342, 0.0001113, -0.001865, 0.00000007336, -0.000002809],
        [0.0000000005687, -0.0000003264, -0.00001718, 0.0000000004328, 0.000000009811],
        [-0.01218, -0.01105, -0.03769, 0.00002044, -0.0005552],
        [0.0000003865, 0.0001828, -0.01196, 0.0000003413, 0.000008384],
        [0.000000008768, -0.000001641, 0.00009801, -0.000000003665, -0.00000003997],
        [-0.0000000002095, 0.0000000151, 0.000000666, -0.00000000002791, -0.0000000003466],
        [-0.00006823, -0.0001208, -0.003776, 0.0000002943, 0.000003038],
        [0.00000002137, 0.000002992, -0.00005611, -0.0000000009646, -0.00000007435],
        [-0.0000000004271, 0.000000001455, -0.0000007811, 0.00000000003174, 0.0000000007442],
        [0.0000001297, 0.000004927, -0.0001504, -0.0000000008666, 0.00000006669],
        [-0.0000000005407, -0.0000001325, 0.000007373, -0.0000000000004573, -0.0000000009105],
        [0.00000002363, -0.00000007727, 0.000006433, -0.0000000002033, -0.0000000008472]
        ])

        self.setMelinderMatrix(coeffs)
