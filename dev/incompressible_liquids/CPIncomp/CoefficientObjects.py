import numpy as np
from CPIncomp.DataObjects import SolutionData

class CoefficientData(SolutionData):
    """ 
    A class to convert parameter arrays from different other sources
    """ 
    def __init__(self):
        SolutionData.__init__(self) 
        self.reference = "Some other software"
        
    def convertSecCoolArray(self, array):
        if len(array)!=18:
            raise ValueError("The lenght is not equal to 18!")
        
        self.reference = "SecCool software"
        array = np.array(array)        
        tmp = np.zeros((4,6))
        
        tmp[0][0] = array[0]
        tmp[0][1] = array[1]
        tmp[0][2] = array[2]
        tmp[0][3] = array[3]
        tmp[0][4] = array[4]
        tmp[0][5] = array[5]
        
        tmp[1][0] = array[6]
        tmp[1][1] = array[7]
        tmp[1][2] = array[8]
        tmp[1][3] = array[9]
        tmp[1][4] = array[10]
        #tmp[1][5] = array[11]
        
        tmp[2][0] = array[11]
        tmp[2][1] = array[12]
        tmp[2][2] = array[13]
        tmp[2][3] = array[14]
        #tmp[2][4] = array[4]
        #tmp[2][5] = array[5]
        
        tmp[3][0] = array[15]
        tmp[3][1] = array[16]
        tmp[3][2] = array[17]
        #tmp[3][3] = array[3]
        #tmp[3][4] = array[4]
        #tmp[3][5] = array[5]
        
        # Concentration is no longer handled in per cent!
        for i in range(6):
            tmp.T[i] *= 100.0**i 
                
        return tmp


    def convertMelinderArray(self, array):
        """The same function as the SecCool converter, 
        the original source code is slightly different though.
        That is why the implementation is in a transposed form..."""
        
        if len(array)!=18:
            raise ValueError("The lenght is not equal to 18!")
        
        self.reference = "Melinder Book"
        array = np.array(array)        
        tmp = np.zeros((6,4))
        
        tmp[0][0] = array[0] 
        tmp[0][1] = array[6] 
        tmp[0][2] = array[11]
        tmp[0][3] = array[15]
        
        tmp[1][0] = array[1] 
        tmp[1][1] = array[7] 
        tmp[1][2] = array[12]
        tmp[1][3] = array[16]
        
        tmp[2][0] = array[2] 
        tmp[2][1] = array[8] 
        tmp[2][2] = array[13]
        tmp[2][3] = array[17]
        
        tmp[3][0] = array[3] 
        tmp[3][1] = array[9] 
        tmp[3][2] = array[14]
        
        tmp[4][0] = array[4] 
        tmp[4][1] = array[10]
        
        tmp[5][0] = array[5] 
        
        # Concentration is no longer handled in per cent!
        for i in range(6):
            tmp[i] *= 100.0**i 
                
        return tmp.T
    
    def convertMelinderMatrix(self, array):
        """Function to convert the full coefficient array
        from the very first CoolProp implementation
        based on the book by Melinder"""
        if len(array)!=18:
            raise ValueError("The lenght is not equal to 18!")
        if len(array[0])!=5:
            raise ValueError("The lenght is not equal to 5!")
        array = np.array(array)
        tmp = np.zeros((18,5))
        
        for j in range(5):
            tmp[ 0][j] = array[ 0][j]
            tmp[ 1][j] = array[ 4][j]
            tmp[ 2][j] = array[ 8][j]
            tmp[ 3][j] = array[12][j]
            tmp[ 4][j] = array[15][j]
            tmp[ 5][j] = array[17][j]
            tmp[ 6][j] = array[ 1][j]
            tmp[ 7][j] = array[ 5][j]
            tmp[ 8][j] = array[ 9][j]
            tmp[ 9][j] = array[13][j]
            tmp[10][j] = array[16][j]
            tmp[11][j] = array[ 2][j]
            tmp[12][j] = array[ 6][j]
            tmp[13][j] = array[10][j]
            tmp[14][j] = array[14][j]
            tmp[15][j] = array[ 3][j]
            tmp[16][j] = array[ 7][j]
            tmp[17][j] = array[11][j]
            
        return tmp
    
    
    
    
    

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
        self.Tmax =  20 + 273.15
        self.Tmin = -50 + 273.15
        self.xmax = 0.5
        self.xmin = 0.0
        self.TminPsat =  20 + 273.15
    
        self.Tbase =  -4.48 + 273.15
        self.xbase =  31.57 / 100.0
       
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

        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_POLYNOMIAL
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
        self.T_freeze.coeffs = np.array([[
           27.755555600/100.0,
          -22.973221700, 
          -1.1040507200*100.0, 
          -0.0120762281*100.0*100.0, 
          -9.343458E-05*100.0*100.0*100.0]])  
        
        
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