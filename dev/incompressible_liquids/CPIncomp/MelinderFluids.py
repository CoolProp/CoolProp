from __future__ import division, print_function
import numpy as np
from .DataObjects import PureData, CoefficientData
from .BaseObjects import IncompressibleFitter


class DEBLiquidClass(CoefficientData, PureData):
    """
    Pure fluid according to Melinder's book
    """

    def __init__(self):
        CoefficientData.__init__(self)
        PureData.__init__(self)
        self.name = "DEB"
        self.description = "Diethylbenzene mixture - Dowtherm J"
        self.reference = "Melinder2010"

        self.Tmin = -80.0 + 273.15
        self.Tmax = 100.0 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.0

        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.density.coeffs = IncompressibleFitter.shapeArray(np.array([1076.5, -0.731182]))

        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.specific_heat.coeffs = IncompressibleFitter.shapeArray(np.array([999.729, 2.87576]))

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
        _, _, self.viscosity.coeffs = IncompressibleFitter.shapeArray(np.array([3.5503, -0.0566396, 7.03331e-05]))

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.conductivity.coeffs = IncompressibleFitter.shapeArray(np.array([0000.189132, -2.06364e-04]))

        self.density.source = self.density.SOURCE_COEFFS
        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.viscosity.source = self.viscosity.SOURCE_COEFFS


class HCMLiquidClass(CoefficientData, PureData):
    """
    Pure fluid according to Melinder's book
    """

    def __init__(self):
        CoefficientData.__init__(self)
        PureData.__init__(self)
        self.name = "HCM"
        self.description = "Hydrocarbon mixture - Gilotherm D12"
        self.reference = "Melinder2010"

        self.Tmin = -80.0 + 273.15
        self.Tmax = 100.0 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.0

        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.density.coeffs = IncompressibleFitter.shapeArray(np.array([971.725, -0.718788]))

        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.specific_heat.coeffs = IncompressibleFitter.shapeArray(np.array([844.023, 4.31212]))

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
        _, _, self.viscosity.coeffs = IncompressibleFitter.shapeArray(np.array([18.3237, -0.14706, 0.000209096]))

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.conductivity.coeffs = IncompressibleFitter.shapeArray(np.array([0000.153716, -1.51212e-04]))

        self.density.source = self.density.SOURCE_COEFFS
        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.viscosity.source = self.viscosity.SOURCE_COEFFS


class HFELiquidClass(CoefficientData, PureData):
    """
    Pure fluid according to Melinder's book
    """

    def __init__(self):
        CoefficientData.__init__(self)
        PureData.__init__(self)
        self.name = "HFE"
        self.description = "Hydrofluoroether - HFE-7100 3M Novec"
        self.reference = "Melinder2010"

        self.Tmin = -80.0 + 273.15
        self.Tmax = 100.0 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.0

        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.density.coeffs = IncompressibleFitter.shapeArray(np.array([1822.37, -0.918485]))

        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.specific_heat.coeffs = IncompressibleFitter.shapeArray(np.array([871.834, 0.858788]))

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
        _, _, self.viscosity.coeffs = IncompressibleFitter.shapeArray(np.array([-4.22878, -0.0114765, 7.39823e-06]))

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.conductivity.coeffs = IncompressibleFitter.shapeArray(np.array([9.92958e-01, -8.33333e-05]))

        self.density.source = self.density.SOURCE_COEFFS
        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.viscosity.source = self.viscosity.SOURCE_COEFFS


class PMS1LiquidClass(CoefficientData, PureData):
    """
    Pure fluid according to Melinder's book
    """

    def __init__(self):
        CoefficientData.__init__(self)
        PureData.__init__(self)
        self.name = "PMS1"
        self.description = "Polydimethylsiloxan 1 - Baysilone KT3"
        self.reference = "Melinder2010"

        self.Tmin = -80.0 + 273.15
        self.Tmax = 100.0 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.0

        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.density.coeffs = IncompressibleFitter.shapeArray(np.array([1172.35, -0.9025]))

        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.specific_heat.coeffs = IncompressibleFitter.shapeArray(np.array([1223.69, 1.48417]))

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
        _, _, self.viscosity.coeffs = IncompressibleFitter.shapeArray(np.array([6.36183, -0.0636352, 7.51428e-05]))

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.conductivity.coeffs = IncompressibleFitter.shapeArray(np.array([0000.207526, -2.84167e-04]))

        self.density.source = self.density.SOURCE_COEFFS
        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.viscosity.source = self.viscosity.SOURCE_COEFFS


class PMS2LiquidClass(CoefficientData, PureData):
    """
    Pure fluid according to Melinder's book
    """

    def __init__(self):
        CoefficientData.__init__(self)
        PureData.__init__(self)
        self.name = "PMS2"
        self.description = "Polydimethylsiloxan 2 - Syltherm XLT"
        self.reference = "Melinder2010"

        self.Tmin = -80.0 + 273.15
        self.Tmax = 100.0 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.0

        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.density.coeffs = IncompressibleFitter.shapeArray(np.array([1155.94, -1.02576]))

        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.specific_heat.coeffs = IncompressibleFitter.shapeArray(np.array([1153.55, 2.10788]))

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
        _, _, self.viscosity.coeffs = IncompressibleFitter.shapeArray(np.array([5.66926, -0.065582, 8.09988e-05]))

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.conductivity.coeffs = IncompressibleFitter.shapeArray(np.array([0000.172305, -2.11212e-04]))

        self.density.source = self.density.SOURCE_COEFFS
        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.viscosity.source = self.viscosity.SOURCE_COEFFS


class SABLiquidClass(CoefficientData, PureData):
    """
    Pure fluid according to Melinder's book
    """

    def __init__(self):
        CoefficientData.__init__(self)
        PureData.__init__(self)
        self.name = "SAB"
        self.description = "Synthetic alkyl benzene - Marlotherm X"
        self.reference = "Melinder2010"

        self.Tmin = -80.0 + 273.15
        self.Tmax = 100.0 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.0

        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.density.coeffs = IncompressibleFitter.shapeArray(np.array([1102.34, -0.801667]))

        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.specific_heat.coeffs = IncompressibleFitter.shapeArray(np.array([1360.94, 1.51667]))

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
        _, _, self.viscosity.coeffs = IncompressibleFitter.shapeArray(np.array([5.21288, -0.0665792, 8.5066e-05]))

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.conductivity.coeffs = IncompressibleFitter.shapeArray(np.array([0000.208374, -2.61667e-04]))

        self.density.source = self.density.SOURCE_COEFFS
        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.viscosity.source = self.viscosity.SOURCE_COEFFS


class HCBLiquidClass(CoefficientData, PureData):
    """
    Pure fluid according to Melinder's book
    """

    def __init__(self):
        CoefficientData.__init__(self)
        PureData.__init__(self)
        self.name = "HCB"
        self.description = "Hydrocarbon blend - Dynalene MV"
        self.reference = "Melinder2010"

        self.Tmin = -80.0 + 273.15
        self.Tmax = 100.0 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.0

        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.density.coeffs = IncompressibleFitter.shapeArray(np.array([1071.78, -0.772024]))

        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.specific_heat.coeffs = IncompressibleFitter.shapeArray(np.array([761.393, 3.52976]))

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
        _, _, self.viscosity.coeffs = IncompressibleFitter.shapeArray(np.array([7.16819, -0.0863212, 0.000130604]))

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.conductivity.coeffs = IncompressibleFitter.shapeArray(np.array([0000.203186, -2.3869e-04]))

        self.density.source = self.density.SOURCE_COEFFS
        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.viscosity.source = self.viscosity.SOURCE_COEFFS


class TCOLiquidClass(CoefficientData, PureData):
    """
    Pure fluid according to Melinder's book
    """

    def __init__(self):
        CoefficientData.__init__(self)
        PureData.__init__(self)
        self.name = "TCO"
        self.description = "Citrus oil terpene - d-Limonene"
        self.reference = "Melinder2010"

        self.Tmin = -80.0 + 273.15
        self.Tmax = 100.0 + 273.15
        self.TminPsat = self.Tmax
        self.Tbase = 0.0

        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.density.coeffs = IncompressibleFitter.shapeArray(np.array([1071.02, -0.778166]))

        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.specific_heat.coeffs = IncompressibleFitter.shapeArray(np.array([223.775, 5.2159]))

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPPOLYNOMIAL
        _, _, self.viscosity.coeffs = IncompressibleFitter.shapeArray(np.array([-3.47971, -0.0107031, 1.14086e-06]))

        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        _, _, self.conductivity.coeffs = IncompressibleFitter.shapeArray(np.array([0000.174156, -1.85052e-04]))

        self.density.source = self.density.SOURCE_COEFFS
        self.specific_heat.source = self.specific_heat.SOURCE_COEFFS
        self.conductivity.source = self.conductivity.SOURCE_COEFFS
        self.viscosity.source = self.viscosity.SOURCE_COEFFS


class EGSolution(CoefficientData):
    def __init__(self):
        CoefficientData.__init__(self)
        self.name = "MEG"
        self.description = "Ethylene Glycol - aq"
        self.reference = "Melinder2010"

        self.Tmin = -100 + 273.15
        self.Tmax = 100 + 273.15
        self.TminPsat = self.Tmax

        self.xmin = 0.0
        self.xmax = 0.6
        self.xid = self.ifrac_mass

        self.Tbase = 31.728 + 273.15
        self.xbase = 30.8462 / 100.0

        coeffs = np.array([
            [-15.25, 1034, 3737, 0.472, 0.4705],
            [-0.000001566, -0.4781, 2.93, 0.0008903, -0.0255],
            [-0.0000002278, -0.002692, -0.004675, -0.000001058, 0.0001782],
            [0.000000002169, 0.000004725, -0.00001389, -0.000000002789, -0.0000007669],
            [-0.808, 1.311, -17.99, -0.004286, 0.02471],
            [-0.000001339, -0.006876, 0.1046, -0.00001473, -0.0001171],
            [0.00000002047, 0.00004805, -0.0004147, 0.0000001059, 0.000001052],
            [-0.00000000002717, 0.0000000169, 0.0000001847, -0.0000000001142, -0.00000001634],
            [-0.01334, 0.0000749, -0.09933, 0.00001747, 0.000003328],
            [0.00000006322, 0.00007855, 0.0003516, 0.00000006814, 0.000001086],
            [0.0000000002373, -0.0000003995, 0.000005109, -0.000000003612, 0.00000001051],
            [-0.000000000002183, 0.000000004982, -0.00000007138, 0.000000000002365, -0.0000000006475],
            [-0.00007293, -0.0001062, 0.00261, 0.00000003017, 0.000001659],
            [0.000000001764, 0.000001229, -0.000001189, -0.000000002412, 0.000000003157],
            [-0.00000000002442, -0.00000001153, -0.0000001643, 0.00000000004004, 0.0000000004063],
            [0.000001006, -0.0000009623, 0.00001537, -0.000000001322, 0.00000003089],
            [-0.00000000007662, -0.00000007211, -0.0000004272, 0.00000000002555, 0.0000000001831],
            [0.00000000114, 0.00000004891, -0.000001618, 0.00000000002678, -0.000000001865]
        ])

        self.setMelinderMatrix(coeffs)


class PGSolution(CoefficientData):
    def __init__(self):
        CoefficientData.__init__(self)
        self.name = "MPG"
        self.description = "Propylene Glycol - aq"
        self.reference = "Melinder2010"

        self.Tmin = -100 + 273.15
        self.Tmax = 100 + 273.15
        self.TminPsat = self.Tmax

        self.xmin = 0.0
        self.xmax = 0.6
        self.xid = self.ifrac_mass

        self.Tbase = 32.7083 + 273.15
        self.xbase = 30.7031 / 100.0

        coeffs = np.array([
            [-13.25, 1018, 3882, 0.4513, 0.6837],
            [-0.0000382, -0.5406, 2.699, 0.0007955, -0.03045],
            [0.0000007865, -0.002666, -0.001659, 0.00000003482, 0.0002525],
            [-0.000000001733, 0.00001347, -0.00001032, -0.000000005966, -0.000001399],
            [-0.6631, 0.7604, -13.04, -0.004795, 0.03328],
            [0.000006774, -0.00945, 0.0507, -0.00001678, -0.0003984],
            [-0.00000006242, 0.00005541, -0.00004752, 0.00000008941, 0.000004332],
            [-0.0000000007819, -0.0000001343, 0.000001522, 0.0000000001493, -0.0000000186],
            [-0.01094, -0.002498, -0.1598, 0.00002076, 0.00005453],
            [0.00000005332, 0.000027, 0.00009534, 0.0000001563, -0.000000086],
            [-0.000000004169, -0.0000004018, 0.00001167, -0.000000004615, -0.00000001593],
            [0.00000000003288, 0.000000003376, -0.0000000487, 0.000000000009897, -0.00000000004465],
            [-0.0002283, -0.000155, 0.0003539, -0.00000009083, -0.0000039],
            [-0.00000001131, 0.000002829, 0.00003102, -0.000000002518, 0.0000001054],
            [0.0000000001918, -0.000000007175, -0.000000295, 0.00000000006543, -0.000000001589],
            [-0.000003409, -0.000001131, 0.00005, -0.0000000005952, -0.00000001587],
            [0.00000000008035, -0.00000002221, -0.0000007135, -0.00000000003605, 0.0000000004475],
            [0.00000001465, 0.00000002342, -0.0000004959, 0.00000000002104, 0.000000003564]
        ])

        self.setMelinderMatrix(coeffs)


class EASolution(CoefficientData):
    def __init__(self):
        CoefficientData.__init__(self)
        self.name = "MEA"
        self.description = "Ethyl Alcohol (Ethanol) - aq"
        self.reference = "Melinder2010"

        self.Tmin = -100 + 273.15
        self.Tmax = 40 + 273.15
        self.TminPsat = self.Tmax

        self.xmin = 0.0
        self.xmax = 0.6
        self.xid = self.ifrac_mass

        self.Tbase = 8.1578 + 273.15
        self.xbase = 29.2361 / 100.0

        coeffs = np.array([
            [-19.41, 961.9, 4204, 0.4067, 1.474],
            [-0.0003668, -0.5222, 2.319, 0.0006775, -0.04745],
            [-0.00004005, -0.003281, -0.03042, 0.0000003105, 0.0004314],
            [0.000001524, 0.00001569, 0.000686, -0.00000002, -0.000003023],
            [-0.954, -1.433, -21.02, -0.005008, 0.01565],
            [-0.00001209, -0.01989, 0.4927, -0.00002377, -0.00004106],
            [0.000002877, 0.000187, -0.003072, -0.00000003216, -0.000005135],
            [-0.00000004394, -0.0000009154, -0.0000569, 0.00000000008362, 0.00000007004],
            [-0.002648, -0.0226, -0.3714, 0.00002801, -0.0008435],
            [-0.0000003173, 0.0002281, -0.002335, 0.0000002669, 0.0000164],
            [0.000000008652, -0.00000008581, -0.0000196, -0.000000003606, -0.0000001091],
            [-0.0000000003717, 0.000000004056, 0.0000007461, 0.00000000001552, -0.000000001967],
            [0.0003851, -0.000169, 0.01743, -0.00000002009, 0.000007552],
            [0.0000000134, 0.000008594, -0.0002969, -0.000000006813, -0.0000001118],
            [-0.000000002091, -0.00000009607, 0.000001901, 0.0000000001429, 0.000000001899],
            [-0.0000002858, 0.00001291, -0.00006292, -0.000000001506, 0.0000001529],
            [0.0000000009312, -0.000000159, 0.000005353, 0.0000000001167, -0.0000000009481],
            [-0.000000167, -0.00000008318, -0.00000829, -0.00000000001653, -0.00000000413]
        ])

        self.setMelinderMatrix(coeffs)


class MASolution(CoefficientData):
    def __init__(self):
        CoefficientData.__init__(self)
        self.name = "MMA"
        self.description = "Methyl Alcohol (Methanol) - aq"
        self.reference = "Melinder2010"

        self.Tmin = -100 + 273.15
        self.Tmax = 40 + 273.15
        self.TminPsat = self.Tmax

        self.xmin = 0.0
        self.xmax = 0.6
        self.xid = self.ifrac_mass

        self.Tbase = 3.5359 + 273.15
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


class GLSolution(CoefficientData):
    def __init__(self):
        CoefficientData.__init__(self)
        self.name = "MGL"
        self.description = "Glycerol - aq"
        self.reference = "Melinder2010"

        self.Tmin = -100 + 273.15
        self.Tmax = 40 + 273.15
        self.TminPsat = self.Tmax

        self.xmin = 0.0
        self.xmax = 0.6
        self.xid = self.ifrac_mass

        self.Tbase = 8.9110 + 273.15
        self.xbase = 36.1905 / 100.0

        coeffs = np.array([
            [-13, 1093, 3486, 0.4532, 1.52],
            [-0.0008638, -0.3624, 3.766, 0.0009782, -0.03729],
            [0.000006895, -0.002451, -0.0001222, -0.0000001196, 0.0003572],
            [0.0000005229, 0.00001547, 0.00003219, -0.000000008778, -0.000003648],
            [-0.5742, 2.74, -22.5, -0.003223, 0.0445],
            [-0.000007991, -0.008373, 0.1811, -0.00001951, -0.0002688],
            [-0.0000007515, 0.00009596, 0.0003766, 0.0000001178, 0.0000008876],
            [0.0000000171, -0.0000006999, -0.0000173, 0.0000000001048, -0.00000002209],
            [-0.009119, 0.004081, -0.03258, 0.000005539, 0.0003633],
            [0.000002973, 0.00001808, 0.0007249, -0.00000006878, -0.000004088],
            [0.000000002379, -0.0000008516, -0.00002133, -0.000000002587, 0.00000006219],
            [-0.000000001237, 0.00000001441, 0.0000004907, 0.00000000002262, 0.0000000006331],
            [-0.0001641, -0.00004744, 0.002922, -0.0000002073, 0.000001069],
            [-0.000000005313, 0.000001833, -0.00005346, -0.0000000002235, -0.0000001248],
            [0.0000000004546, -0.00000003077, 0.00000023, -0.00000000001421, 0.000000003019],
            [-0.000002408, -0.000001415, 0.00002238, 0.000000003, 0.00000005729],
            [-0.000000001682, 0.00000001863, -0.0000005383, 0.0000000001604, -0.00000000178],
            [-0.000000007734, -0.000000006097, -0.0000005944, 0.0000000002221, 0.000000002116]
        ])

        self.setMelinderMatrix(coeffs)


class AMSolution(CoefficientData):
    def __init__(self):
        CoefficientData.__init__(self)
        self.name = "MAM"
        self.description = "Ammonia (NH3) - aq"
        self.reference = "Melinder2010"

        self.Tmin = -100 + 273.15
        self.Tmax = 30 + 273.15
        self.TminPsat = self.Tmax

        self.xmin = 0.0
        self.xmax = 0.3
        self.xid = self.ifrac_mass

        self.Tbase = -4.6490 + 273.15
        self.xbase = 16.0784 / 100.0

        coeffs = np.array([
            [-25.76, 944.5, 4233, 0.4551, 0.9255],
            [-0.0001817, -0.2743, -1.618, 0.001673, -0.03439],
            [0.00001204, -0.003113, 0.0161, -0.000002214, 0.0003217],
            [0.0000005567, 0.000003349, 0.00001662, 0.0000001228, -0.000004544],
            [-2.385, -2.914, 1.145, -0.005216, 0.01327],
            [0.00002315, -0.02307, 0.02715, 0.000003544, 0.0001856],
            [0.0000001415, 0.0001341, 0.001072, 0.000001057, -0.00001646],
            [-0.00000004244, -0.0000005151, -0.00005266, -0.00000003474, 0.0000003004],
            [-0.07636, 0.02262, -0.001965, 0.00008909, -0.0005979],
            [-0.00000229, 0.00006645, 0.003472, 0.000003526, -0.00002184],
            [-0.000000262, -0.0000087, -0.00009051, -0.0000001782, 0.000001297],
            [0.0000000001786, 0.00000008999, 0.000002106, 0.000000001858, -0.00000001141],
            [-0.00261, -0.0006685, 0.002131, 0.000006807, -0.0001097],
            [-0.000000376, -0.000001002, 0.00004117, -0.0000003394, 0.00000167],
            [0.0000000136, 0.0000003309, 0.0000004446, 0.000000008315, -0.00000003377],
            [-0.000073, 0.000001635, 0.0002136, -0.0000001898, 0.000003065],
            [0.00000003524, 0.0000004465, -0.00001354, 0.000000006304, -0.00000006166],
            [-0.000001075, 0.0000006298, -0.000008551, -0.00000001361, 0.0000003244]
        ])

        self.setMelinderMatrix(coeffs)


class KCSolution(CoefficientData):
    def __init__(self):
        CoefficientData.__init__(self)
        self.name = "MKC"
        self.description = "Potassium Carbonate (K2CO3) - aq"
        self.reference = "Melinder2010"

        self.Tmin = -100 + 273.15
        self.Tmax = 40 + 273.15
        self.TminPsat = self.Tmax

        self.xmin = 0.0
        self.xmax = 0.4
        self.xid = self.ifrac_mass

        self.Tbase = 11.2422 + 273.15
        self.xbase = 22.0833 / 100.0

        coeffs = np.array([
            [-10.3, 1216, 3217, 0.5622, 0.8063],
            [-0.0001575, -0.4114, 1.492, 0.001656, -0.02362],
            [0.000003598, -0.003383, -0.001539, 0.000002108, 0.0001851],
            [0.000000004324, 0.0000299, -0.00003546, -0.00000004482, -0.000002372],
            [-0.7786, 10.64, -37.33, -0.000892, 0.03624],
            [0.00001112, -0.007614, -0.01775, 0.000002031, -0.00001262],
            [-0.0000003479, 0.00005214, 0.0005416, -0.0000002616, -0.0000003022],
            [-0.000000001244, -0.0000008087, 0.00000659, 0.0000000007965, -0.000000007761],
            [-0.02766, 0.04413, 0.2023, -0.0000005233, 0.0006659],
            [0.000001616, 0.00007806, 0.004553, -0.0000002962, -0.00001611],
            [-0.00000001681, -0.000001173, 0.00003587, -0.000000009174, 0.000000153],
            [0.00000000003847, 0.00000005658, -0.0000003707, 0.0000000001027, 0.000000001061],
            [-0.0008226, -0.0001333, 0.01971, -0.0000009283, 0.00001077],
            [-0.00000004913, -0.000002381, 0.0001367, -0.00000001814, -0.00000009796],
            [0.000000001395, 0.0000001696, -0.000003667, 0.0000000008767, 0.00000000307],
            [-0.000002372, -0.00001246, 0.0003405, -0.00000001011, -0.0000001097],
            [-0.000000002886, 0.0000002708, -0.00001676, 0.0000000008471, 0.000000007825],
            [0.0000003251, 0.0000002035, -0.00003488, 0.000000001311, -0.000000008453]
        ])

        self.setMelinderMatrix(coeffs)


class CASolution(CoefficientData):
    def __init__(self):
        CoefficientData.__init__(self)
        self.name = "MCA"
        self.description = "Calcium Chloride (CaCl2) - aq"
        self.reference = "Melinder2010"

        self.Tmin = -100.0 + 273.15
        self.Tmax = 40.0 + 273.15
        self.TminPsat = self.Tmax

        self.xmin = 0.0
        self.xmax = 0.3
        self.xid = self.ifrac_mass

        self.Tbase = 7.52570 + 273.15
        self.xbase = 18.7414 / 100.0

        coeffs = np.array([
            [-16.21, 1171, 3133, 0.558, 0.8939],
            [-0.0001344, -0.1463, 2.81, 0.00146, -0.02647],
            [0.000005073, -0.001302, -0.01563, 0.0000003861, 0.0001718],
            [-0.0000000482, -0.0001871, -0.00001233, 0.00000001307, -0.0000007918],
            [-1.555, 9.847, -44.8, -0.00115, 0.04389],
            [0.00002146, -0.02488, 0.03271, -0.00001008, 0.0002102],
            [-0.0000015, -0.000553, -0.001205, -0.0000000654, -0.0000008688],
            [0.00000002219, 0.00001665, 0.000009346, 0.0000000004728, -0.00000004353],
            [-0.05496, 0.03389, 0.9511, -0.00001784, 0.0009121],
            [0.0000009415, -0.002302, -0.005191, 0.0000003496, 0.000003993],
            [0.00000006185, 0.00003526, 0.0002282, -0.00000000484, 0.0000003785],
            [-0.000000001723, 0.0000002788, -0.000000929, -0.0000000002011, -0.000000009979],
            [-0.002624, 0.001062, 0.01472, -0.0000004415, 0.000001963],
            [-0.0000001082, 0.00006291, 0.0001615, -0.000000003462, -0.0000004892],
            [0.000000003036, 0.000001806, 0.000005073, 0.0000000003922, 0.00000000001526],
            [-0.0001726, 0.00002785, -0.001346, 0.00000004699, 0.0000002997],
            [-0.000000004396, 0.000006423, -0.000009754, 0.0000000006179, -0.00000003631],
            [-0.000004494, -0.000001811, -0.00006674, 0.000000002929, 0.00000003435]
        ])

        self.setMelinderMatrix(coeffs)


class MGSolution(CoefficientData):
    def __init__(self):
        CoefficientData.__init__(self)
        self.name = "MMG"
        self.description = "MgCl2 - aq"
        self.reference = "Melinder2010"

        self.Tmin = -100.0 + 273.15
        self.Tmax = 40.0 + 273.15
        self.TminPsat = self.Tmax

        self.xmin = 0.0
        self.xmax = 0.3
        self.xid = self.ifrac_mass

        self.Tbase = 9.3163 + 273.15
        self.xbase = 14.1327 / 100.0

        coeffs = np.array([
            [-15.12, 1124, 3365, 0.5461, 0.9573],
            [-0.0004843, -0.3072, 2.229, 0.001784, -0.03065],
            [0.00001113, -0.003295, -0.004627, -0.0000008171, 0.0001115],
            [0.0000001858, 0.00001015, 0.00009186, -0.00000006594, -0.000002923],
            [-1.885, 9.071, -52.22, -0.00273, 0.04968],
            [-0.00005461, -0.006513, 0.1162, -0.00001483, 0.0001559],
            [0.000003579, 0.00004664, 0.001249, 0.000000385, -0.00001796],
            [-0.00000003999, 0.000002287, 0.000002421, -0.000000005556, 0.0000003051],
            [-0.05455, 0.02449, 0.6202, 0.000008675, -0.002722],
            [0.00001887, 0.00003574, 0.002337, -0.000001489, -0.00001753],
            [-0.0000003171, 0.000004337, 0.0000724, -0.00000003882, 0.000002021],
            [-0.000000006246, 0.0000006044, -0.000003613, 0.0000000009282, 0.00000002614],
            [-0.0007257, 0.003402, 0.01052, 0.0000008651, 0.00009351],
            [0.0000007588, -0.0001409, -0.000459, 0.0000001992, -0.000008353],
            [-0.00000004102, 0.000001857, -0.00002477, -0.000000001196, 0.0000001901],
            [-0.0003208, 0.0003344, -0.001067, -0.0000004779, 0.00007364],
            [-0.0000000492, -0.00000683, -0.0001048, 0.00000001797, -0.0000004014],
            [-0.00001794, 0.000007239, -0.0000696, -0.00000002503, 0.000003717]
        ])

        self.setMelinderMatrix(coeffs)


class NASolution(CoefficientData):
    def __init__(self):
        CoefficientData.__init__(self)
        self.name = "MNA"
        self.description = "Sodium Chloride (NaCl) - aq"
        self.reference = "Melinder2010"

        self.Tmin = -100.0 + 273.15
        self.Tmax = 40.0 + 273.15
        self.TminPsat = self.Tmax

        self.xmin = 0.0
        self.xmax = 0.23
        self.xid = self.ifrac_mass

        self.Tbase = 12.6179 + 273.15
        self.xbase = 13.3897 / 100.0

        coeffs = np.array([
            [-9.383, 1099, 3593, 0.5736, 0.4369],
            [-0.00002581, -0.3758, 1.669, 0.001595, -0.02666],
            [0.000001423, -0.002023, -0.02019, -0.0000004003, 0.0002035],
            [0, 0, 0, 0, 0],
            [-0.9039, 7.723, -32.48, -0.0009383, 0.02346],
            [0.000002578, -0.01426, -0.03331, -0.00001248, -0.00005368],
            [-0.00000003318, 0.0001535, -0.001164, 0.0000003353, 0.000002871],
            [0, 0, 0, 0, 0],
            [-0.02204, 0.02567, 0.6453, -0.00001057, 0.0004276],
            [0.0000001192, 0.0003994, -0.009314, 0.0000004158, -0.000004526],
            [-0.000000008993, -0.000007281, 0.0002236, -0.00000002032, 0.0000001838],
            [0, 0, 0, 0, 0],
            [-0.0004827, 0.0001108, -0.01629, -0.0000004853, 0.000007386],
            [-0.00000001467, 0.0000003522, 0.0007927, -0.00000001587, 0.0000005437],
            [0, 0, 0, 0, 0],
            [0.000002247, -0.00001686, 0.0002331, -0.000000004654, 0.0000004688],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0]
        ])

        self.setMelinderMatrix(coeffs)


class KASolution(CoefficientData):
    def __init__(self):
        CoefficientData.__init__(self)
        self.name = "MKA"
        self.description = "Potassium Acetate (CH3CO2K) - aq"
        self.reference = "Melinder2010"

        self.Tmin = -100.0 + 273.15
        self.Tmax = 40.0 + 273.15
        self.TminPsat = self.Tmax

        self.xmin = 0.0
        self.xmax = 0.45
        self.xid = self.ifrac_mass

        self.Tbase = 6.7757 + 273.15
        self.xbase = 25.6757 / 100.0

        coeffs = np.array([
            [-17.04, 1138, 3327, 0.4958, 1.042],
            [-0.0001082, -0.3565, 1.806, 0.00134, -0.03071],
            [0.000006892, -0.00202, -0.001766, 0.00000006837, 0.0002819],
            [-0.0000001397, 0.000004205, 0.00004357, 0.000000002637, -0.00000219],
            [-1.228, 5.796, -28.95, -0.002931, 0.03405],
            [0.0000002302, -0.0079, 0.04131, -0.00001477, -0.0001511],
            [-0.0000008711, 0.00002753, 0.0004429, 0.00000004659, 0.000001172],
            [0.00000002016, 0.0000000514, 0.00001125, 0.0000000002886, -0.00000002379],
            [-0.03862, 0.01306, 0.04663, 0.00001032, 0.0005017],
            [0.000001565, 0.00006845, 0.0007775, -0.0000002396, -0.000007779],
            [0.000000007565, -0.00000113, 0.00003463, -0.000000004352, 0.00000009125],
            [-0.0000000003063, 0.00000002433, -0.0000007261, -0.0000000000223, -0.000000001888],
            [-0.0004571, -0.001427, -0.001249, -0.0000002024, 0.000005637],
            [-0.00000003734, 0.0000008304, 0.00005115, -0.000000005541, 0.00000002534],
            [0.000000001268, 0.00000001303, -0.000002987, 0.00000000008984, 0.000000001596],
            [0.00002969, 0.000009353, 0.0001659, -0.000000002371, -0.0000002922],
            [-0.000000002817, 0.00000002322, -0.000005193, 0.0000000005573, 0.000000004601],
            [0.000000869, 0.000002285, -0.0000004612, 0.0000000005515, -0.00000000796]
        ])

        self.setMelinderMatrix(coeffs)


class KFSolution(CoefficientData):
    def __init__(self):
        CoefficientData.__init__(self)
        self.name = "MKF"
        self.description = "Potassium Formate (CHKO2) - aq"
        self.reference = "Melinder2010"

        self.Tmin = -100.0 + 273.15
        self.Tmax = 40.0 + 273.15
        self.TminPsat = self.Tmax

        self.xmin = 0.0
        self.xmax = 0.48
        self.xid = self.ifrac_mass

        self.Tbase = 5.89080 + 273.15
        self.xbase = 29.1447 / 100.0

        coeffs = np.array([
            [-20.19, 1189, 3144, 0.5253, 0.8088],
            [-0.0001703, -0.3515, 1.698, 0.001241, -0.02556],
            [-0.000007478, -0.001918, -0.001303, -0.00000003799, 0.0002195],
            [0.0000003761, 0.00003132, 0.00005177, -0.000000002951, -0.000001667],
            [-1.106, 7.044, -29.94, -0.001972, 0.01758],
            [-0.000004203, -0.00656, 0.0229, -0.000006322, 0.00008603],
            [-0.0000005737, 0.00007018, 0.000003898, 0.00000002654, 0.000002498],
            [0.00000001474, -0.000001459, 0.000007391, -0.0000000007324, -0.00000003569],
            [-0.021, 0.01889, 0.2944, -0.00002596, -0.00008372],
            [-0.0000003802, 0.0001034, -0.002482, -0.00000004693, -0.000001601],
            [0.00000006504, -0.000003889, 0.000033, -0.000000004362, 0.000000004815],
            [-0.000000001877, 0.000000009225, -0.0000004871, 0.00000000002899, -0.000000001861],
            [-0.0002638, -0.0002262, 0.001161, -0.0000005062, -0.000001184],
            [0.00000004027, 0.0000003692, 0.00006758, -0.000000005229, -0.0000001331],
            [-0.0000000004789, 0.00000006609, -0.000001277, 0.0000000001035, -0.000000005489],
            [0.0000008327, 0.00002368, -0.0001429, 0.00000001501, 0.000001088],
            [0.0000000009507, 0.00000008766, 0.0000009949, 0.0000000005376, -0.000000007003],
            [0.00000006345, 0.000001061, -0.000003221, 0.0000000005562, 0.00000003098]
        ])

        self.setMelinderMatrix(coeffs)


class LISolution(CoefficientData):
    def __init__(self):
        CoefficientData.__init__(self)
        self.name = "MLI"
        self.description = "Lithium Chloride (LiCl) - aq"
        self.reference = "Melinder2010"

        self.Tmin = -100.0 + 273.15
        self.Tmax = 40.0 + 273.15
        self.TminPsat = self.Tmax

        self.xmin = 0.0
        self.xmax = 0.24
        self.xid = self.ifrac_mass

        self.Tbase = 1.4895 + 273.15
        self.xbase = 14.8000 / 100.0

        coeffs = np.array([
            [-23.29, 1088, 3383, 0.5362, 1.013],
            [0.0006555, -0.1772, 3.958, 0.001454, -0.03062],
            [-0.0001208, -0.002619, -0.0003035, -0.0000000326, 0.000294],
            [0.000002616, 0.000006209, -0.000003477, -0.0000000142, -0.000002719],
            [-3.051, 6.056, -50.36, -0.001855, 0.0392],
            [-0.0003972, -0.008588, 0.4415, -0.00001405, 0.00006246],
            [0.00003674, 0.0001567, -0.0002609, -0.000000005424, -0.000001752],
            [-0.0000005569, -0.000001847, 0.000003651, 0.0000000009821, 0.00000008346],
            [-0.179, 0.02556, 0.6298, 0.00001017, 0.000332],
            [0.00001391, 0.00007194, -0.004384, 0.0000006821, 0.000000784],
            [-0.000001997, -0.00001053, 0.0001039, -0.000000008674, -0.00000031],
            [0.00000002931, 0.00000009716, -0.000001076, 0.0000000001095, 0.00000001381],
            [-0.002917, 0.0009407, 0.04544, 0.0000007084, 0.000002206],
            [0.00000149, -0.000007253, -0.0008787, -0.00000007434, -0.0000006011],
            [-0.00000006904, 0.0000003144, -0.000008457, 0.0000000006988, 0.000000004023],
            [0.0005715, -0.00008105, 0.002527, -0.0000001273, 0.000001745],
            [0.0000001186, 0.000001072, -0.0000358, -0.000000003058, -0.00000007094],
            [0.00002757, -0.000003974, 0.00004058, -0.000000009124, 0.00000006699]
        ])

        self.setMelinderMatrix(coeffs)
