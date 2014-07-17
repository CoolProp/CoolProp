import numpy as np
from CPIncomp.CoefficientObjects import CoefficientData

class MelinderExample(CoefficientData):
    """ 
    Methanol-Water mixture according to Melinder book
    Source: Book
    """ 
    def __init__(self):
        CoefficientData.__init__(self) 
        self.name = "ExampleMelinder"
        self.description = "Methanol solution"
        self.reference = "Melinder-BOOK-2010"
        self.Tmax =  40 + 273.15
        self.Tmin = -50 + 273.15
        self.xmax = 0.6
        self.xmin = 0.0
        self.TminPsat =  self.Tmax
    
        self.Tbase =   3.5359 + 273.15;
        self.xbase =  30.5128 / 100.0
       
        coeffs = np.array([
        [-26.29           , 958.1           ,3887           ,   0.4175            ,   1.153           ],
        [ -0.000002575    ,  -0.4151        ,   7.201       ,   0.0007271         ,  -0.03866         ],
        [ -0.000006732    ,  -0.002261      ,  -0.08979     ,   0.0000002823      ,   0.0002779       ],
        [  0.000000163    ,   0.0000002998  ,  -0.000439    ,   0.000000009718    ,  -0.000001543     ],
        [ -1.187          ,  -1.391         , -18.5         ,  -0.004421          ,   0.005448        ],
        [ -0.00001609     ,  -0.0151        ,   0.2984      ,  -0.00002952        ,   0.0001008       ],
        [  0.000000342    ,   0.0001113     ,  -0.001865    ,   0.00000007336     ,  -0.000002809     ],
        [  0.0000000005687,  -0.0000003264  ,  -0.00001718  ,   0.0000000004328   ,   0.000000009811  ],
        [ -0.01218        ,  -0.01105       ,  -0.03769     ,   0.00002044        ,  -0.0005552       ],
        [  0.0000003865   ,   0.0001828     ,  -0.01196     ,   0.0000003413      ,   0.000008384     ],
        [  0.000000008768 ,  -0.000001641   ,   0.00009801  ,  -0.000000003665    ,  -0.00000003997   ],
        [ -0.0000000002095,   0.0000000151  ,   0.000000666 ,  -0.00000000002791  ,  -0.0000000003466 ],
        [ -0.00006823     ,  -0.0001208     ,  -0.003776    ,   0.0000002943      ,   0.000003038     ],
        [  0.00000002137  ,   0.000002992   ,  -0.00005611  ,  -0.0000000009646   ,  -0.00000007435   ],
        [ -0.0000000004271,   0.000000001455,  -0.0000007811,   0.00000000003174  ,   0.0000000007442 ],
        [  0.0000001297   ,   0.000004927   ,  -0.0001504   ,  -0.0000000008666   ,   0.00000006669   ],
        [ -0.0000000005407,  -0.0000001325  ,   0.000007373 ,  -0.0000000000004573,  -0.0000000009105 ],
        [  0.00000002363  ,  -0.00000007727 ,   0.000006433 ,  -0.0000000002033   ,  -0.0000000008472 ]
        ])
        
        coeffs = self.convertMelinderMatrix(coeffs).T
        
        self.T_freeze.type = self.T_freeze.INCOMPRESSIBLE_POLYNOMIAL
        self.T_freeze.coeffs = self.convertMelinderArray(coeffs[0])
        
        self.density.type = self.density.INCOMPRESSIBLE_POLYNOMIAL
        self.density.coeffs = self.convertMelinderArray(coeffs[1])
        
        self.specific_heat.type = self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL
        self.specific_heat.coeffs = self.convertMelinderArray(coeffs[2])
        
        self.conductivity.type = self.conductivity.INCOMPRESSIBLE_POLYNOMIAL
        self.conductivity.coeffs = self.convertMelinderArray(coeffs[3])
        
        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_POLYNOMIAL
        self.viscosity.coeffs = self.convertMelinderArray(coeffs[4])   

