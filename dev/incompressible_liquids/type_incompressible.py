import numpy as np
import itertools,scipy.interpolate
from quantities.units import viscosity

# Here we define the types. This is done to keep the definitions at one
# central place instead of hiding them somewhere in the data.
 
class IncompressibleData(object):
    """ 
    The work horse of the incompressible classes. 
    Implements both data structures and fitting 
    procedures.
    """ 
    def __init__(self):
        self.INCOMPRESSIBLE_NOT_SET       = 'notdefined'
        self.INCOMPRESSIBLE_POLYNOMIAL    = 'polynomial'
        self.INCOMPRESSIBLE_EXPONENTIAL   = 'exponential'
        self.INCOMPRESSIBLE_EXPPOLYNOMIAL = 'exppolynomial'
        self.INCOMPRESSIBLE_EXPOFFSET     = 'expoffset'
        self.INCOMPRESSIBLE_POLYOFFSET    = 'polyoffset'
        self.INCOMPRESSIBLE_CHEBYSHEV     = 'chebyshev'
        self.type   = self.INCOMPRESSIBLE_NOT_SET
        self.coeffs = None #np.zeros((4,4))
        self.data   = None #np.zeros((10,10))
        
        self.maxLog = numpy.log(numpy.finfo(numpy.float64).max-1)
        self.minLog = -self.maxLog
        
        
        
    ### Base functions that handle the custom data type, just a place holder to show the structure.
    def baseFunction(self, x, y=0.0, xbase=0.0, ybase=0.0, c=None):
        if c==None:
            c = self.coeffs
               
        if self.type==self.INCOMPRESSIBLE_POLYNOMIAL:
            return np.polynomial.polynomial.polyval2d(x-xbase, y-ybase, c)
        
        elif self.type==self.INCOMPRESSIBLE_POLYOFFSET:
            return self.data.basePolyOffset(c, x-xbase, y-ybase)
        
        elif self.type==self.INCOMPRESSIBLE_EXPONENTIAL:
            return self.data.baseExponential(c, x, xbase)
        
        elif self.type==self.INCOMPRESSIBLE_EXPPOLYNOMIAL:
            return np.exp(np.polynomial.polynomial.polyval2d(x-xbase, y-ybase, c))
        
        elif self.type==self.INCOMPRESSIBLE_EXPOFFSET:
            return self.data.baseExponentialOffset(c, x)
        
        else:  raise ValueError("Unknown function.")

    
    def baseExponential(self, data, y, ybase=0.0):
        r=len(data.coeffs)
        c=len(data.coeffs[0]);
        if (r!=3 or c!=1):
            raise ValueError("You have to provide a 3,1 matrix of coefficients, not ({0},{1}).".format(r,c))
        return numpy.exp(numpy.clip( (coefficients[0][0]/ ( y-ybase+coefficients[1][0]) - coefficients[2][0]),self.minLog,self.maxLog))
    
    def baseExponentialOffset(self, data, y):
        r=len(data.coeffs)
        c=len(data.coeffs[0])
        if (r!=4 or c!=1):
            raise ValueError("You have to provide a 4,1 matrix of coefficients, not ({0},{1}).".format(r,c))
        return numpy.exp(numpy.clip( (coefficients[1][0]/ ( y-data.coeffs[0][0]+coefficients[2][0]) - coefficients[3][0]),self.minLog,self.maxLog))

    def basePolyOffset(self, data, y, z=0.0):
        r=len(data.coeffs)
        c=len(data.coeffs[0])
        offset = 0.0
        inV    = 0.0
        coeffs = np.array()
        if (r>0 and c>0):
            offset = data.coeffs[0][0]
            if (r==1 and c>1): # row vector -> function of z
                coeffs = data.coeffs.T[1:c].T
                inV = z
            elif (r>1 and c==1): # column vector -> function of y
                coeffs = data.coeffs[1:r].T
                inV = y
            else:
                raise ValueError("You have to provide a vector (1D matrix) of coefficients, not ({0},{1}).".format(r,c))
            return np.polynomial.polynomial.polyval(inV-offset, data.coeffs) 
        raise ValueError("You have to provide a vector (1D matrix) of coefficients, not ({0},{1}).".format(r,c))


        
    def fit(self, T, x=0.0, Tbase=0.0, xbase=0.0):
        (cr,cc) = self.coeffs.shape
        (dr,dc) = self.data.shape
        (Tr,Tc) = (len(T),1) #T.shape #(len(T),1)
        if (Tc!=1): raise ValueError("Temperature has to be a 2D array with one column.")
        if (Tr!=dr): raise ValueError("Temperature and fitting data have to have the same number of rows.")
        
        if self.type==self.INCOMPRESSIBLE_POLYNOMIAL or self.type==self.INCOMPRESSIBLE_EXPPOLYNOMIAL:
            if (np.max(x)<1e-10): # x not given
                x = np.array([0.0])
                #if (cc==1): # requires 1D coefficients
                if (dc==1):
                    if self.type==self.INCOMPRESSIBLE_POLYNOMIAL:
                        self.coeffs = self.getCoeffs2d(T-Tbase, x-xbase, self.data, cr-1, 0)
                        
                    elif self.type==self.INCOMPRESSIBLE_EXPPOLYNOMIAL:
                        
                        self.coeffs = self.getCoeffs2d(T-Tbase, x-xbase, np.log(self.data), cr-1, 0)
                    else:  raise ValueError("Unknown function.")
                else:
                    raise ValueError("Cannot use 2D data with 1D references")
                #else:
                #    raise ValueError("Cannot use 2D coefficients without concentration input")
            else: # Assume 2D input
                (xr,xc) = (1,len(x))#x.shape#(1,len(x))
                if (xr!=1): raise ValueError("Concentration has to be a 2D array with one row.")
                if (xc!=dc): raise ValueError("Concentration and fitting data have to have the same number of columns.")
                if self.type==self.INCOMPRESSIBLE_POLYNOMIAL:
                    self.coeffs = self.getCoeffs2d(T-Tbase, x-xbase, self.data, cr-1, cc-1)
                    
                elif self.type==self.INCOMPRESSIBLE_EXPPOLYNOMIAL:
                    self.coeffs = self.getCoeffs2d(T-Tbase, x-xbase, np.log(self.data), cr-1, cc-1)
                    
                else:  raise ValueError("Unknown function.")
                #print self.coeffs
        else:
            raise ValueError("Cannot fit that function.")
        

    def getCoeffs1d(self, x, z, order):
        if (len(x)<order+1): 
            raise ValueError("You have only {0} elements and try to fit {1} coefficients, please reduce the order.".format(len(x),order+1))
        A = np.vander(x,order+1)[:,::-1]
        #Anew = np.dot(A.T,A)
        #znew = np.dot(A.T,z)
        #coeffs = np.linalg.solve(Anew, znew)
        coeffs, resids, rank, singulars  = np.linalg.lstsq(A, z)
        return np.reshape(coeffs, (len(x),1))
            

    def getCoeffs2d(self, x_in, y_in, z_in, x_order, y_order):
        x_order += 1
        y_order += 1
        #To avoid overfitting, we only use the upper left triangle of the coefficient matrix
        x_exp = range(x_order)
        y_exp = range(y_order)
        limit = max(x_order,y_order)
        xy_exp = []
        # Construct the upper left triangle of coefficients    
        for i in x_exp:
            for j in y_exp:
                if(i+j<limit): xy_exp.append((i,j))
        
        # Construct input pairs   
        xx, yy = np.meshgrid(x_in,y_in,indexing='ij')
        xx = np.array(xx.flat)
        yy = np.array(yy.flat)
        zz = np.array(z_in.flat)
        
        # TODO: Check for rows with only nan values
        x_num = len(x_in)
        y_num = len(y_in)      
                    
        cols = len(xy_exp)
        eqns = x_num * y_num
        #if (eqns<cols):
        #    raise ValueError("You have only {0} equations and try to fit {1} coefficients, please reduce the order.".format(eqns,cols))   
        if (x_num<x_order):
            raise ValueError("You have only {0} x-entries and try to fit {1} x-coefficients, please reduce the x_order.".format(x_num,x_order))
        if (y_num<y_order):
            raise ValueError("You have only {0} y-entries and try to fit {1} y-coefficients, please reduce the y_order.".format(y_num,y_order))
        
        # Build the functional matrix
        A = np.zeros((eqns,cols))
        for i in range(eqns): # row loop
            for j, (xj,yj) in enumerate(xy_exp): # makes columns
                A[i][j] = xx[i]**xj * yy[i]**yj
        
        # Remove np.nan elements
        mask = np.isfinite(zz)
        A = A[mask]
        zz = zz[mask]
        
        if (len(A) < cols):
            raise ValueError("Your matrix has only {0} valid rows and you try to fit {1} coefficients, please reduce the order.".format(len(A),cols))
        
        coeffs, resids, rank, singulars  = np.linalg.lstsq(A, zz)
        #print resids
        
        #Rearrange coefficients to a matrix shape
        C = np.zeros((x_order,y_order))
        for i, (xi,yi) in enumerate(xy_exp): # makes columns
            C[xi][yi] = coeffs[i]
            
        return C
        
        
    def toJSON(self):
        j = {}
        try:
            j['coeffs'] = self.coeffs.tolist()
        except:
            j['coeffs'] = 'null'
            
        j['type']   = self.type
        return j
    
    
class SolutionData(object):
    """ 
    A base class that defines all the variables needed 
    in order to make a proper fit. You can copy this code
    put in your data and add some documentation for where the
    information came from. 
    """
    def __init__(self):
        self.name        = None # Name of the current fluid
        self.description = None # Description of the current fluid
        self.reference   = None # Reference data for the current fluid
        
        self.Tmax        = None # Maximum temperature in K
        self.Tmin        = None # Minimum temperature in K
        self.xmax        = 1.0 # Maximum concentration
        self.xmin        = 0.0 # Minimum concentration
        self.TminPsat    = None # Minimum saturation temperature in K
        self.Tbase       = 0.0  # Base value for temperature fits
        self.xbase       = 0.0  # Base value for concentration fits
    
        self.temperature   = IncompressibleData() # Temperature for data points in K
        self.concentration = IncompressibleData() # Concentration data points in weight fraction
        self.density       = IncompressibleData() # Density in kg/m3
        self.specific_heat = IncompressibleData() # Heat capacity in J/(kg.K)
        self.viscosity     = IncompressibleData() # Dynamic viscosity in Pa.s
        self.conductivity  = IncompressibleData() # Thermal conductivity in W/(m.K)
        self.saturation_pressure = IncompressibleData() # Saturation pressure in Pa
        self.T_freeze      = IncompressibleData() # Freezing temperature in K
        self.volume2mass   = IncompressibleData() # dd
        self.mass2mole     = IncompressibleData() # dd
        
        self.xref = None
        self.Tref = None
        self.pref = None
        self.href = None
        self.sref = None
        self.uref = None
        self.rhoref = None
        
    
    def rho (self, T, p, x=0.0, c=None):
        if c==None: 
            c=self.density.coeffs
        if self.density.type==self.density.INCOMPRESSIBLE_POLYNOMIAL:
            return np.polynomial.polynomial.polyval2d(T-self.Tbase, x-self.xbase, c)
        else:  raise ValueError("Unknown function.")
    
    def c   (self, T, p, x=0.0, c=None):
        if c==None: 
            c = self.specific_heat.coeffs
        if self.specific_heat.type==self.specific_heat.INCOMPRESSIBLE_POLYNOMIAL:
            return np.polynomial.polynomial.polyval2d(T-self.Tbase, x-self.xbase, c)
        else:  raise ValueError("Unknown function.")
    
    def cp  (self, T, p, x=0.0, c=None):
        return self.c(T,p,x,c)
    
    def cv  (self, T, p, x=0.0, c=None):
        return self.c(T,p,x,c)
    
    #def s   (T, p, x);
    #def u   (T, p, x);
    
    def h   (self, T, p, x=0.0):
        return h_u(T,p,x)
    
    def visc(self, T, p, x=0.0, c=None):
        return self.viscosity.baseFunction(T, x, self.Tbase, self.xbase, c=c)
    
    def cond(self, T, p, x=0.0, c=None):
        return self.conductivity.baseFunction(T, x, self.Tbase, self.xbase, c=c)
        
    def psat(self, T, x=0.0, c=None):
        return self.saturation_pressure.baseFunction(T, x, self.Tbase, self.xbase, c=c)
        
    def Tfreeze(self, p, x=0.0, c=None):
        return self.T_freeze.baseFunction(x, 0.0, self.xbase, 0.0, c=c)
    
    #def V2M (T,           y);
    #def M2M (T,           x);
    
    
    def h_u(self, T, p, x):
        return u(T,p,x)+p/rho(T,p,x)-self.href

    def u_h(self, T, p, x):
        return h(T,p,x)-p/rho(T,p,x)+self.href
    
    def set_reference_state(self, T0, p0, x0=0.0, h0=0.0, s0=0.0):
        self.rhoref = rho(T0,p0,x0)
        self.pref = p0
        self.uref = h0 - p0/rhoref
        self.uref = u(T0,p0,x0)
        self.href = h0
        self.sref = s0
        self.href = h(T0,p0,x0)
        self.sref = s(T0,p0,x0)
        
        
        
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
        
        
        
class PureExample(SolutionData):
    def __init__(self):
        SolutionData.__init__(self) 
        self.name = "ExamplePure"
        self.description = "Heat transfer fluid TherminolD12 by Solutia"
        self.reference = "Solutia data sheet"
        self.Tmax = 150 + 273.15
        self.Tmin =  50 + 273.15
        self.TminPsat =  self.Tmax
        
        self.temperature.data         = np.array([    50 ,     60 ,     70 ,     80 ,     90 ,    100 ,    110 ,    120 ,    130 ,    140 ,    150 ])+273.15 # Kelvin
        self.density.data             = np.array([[  740],[   733],[   726],[   717],[   710],[   702],[   695],[   687],[   679],[   670],[   662]])        # kg/m3
        self.specific_heat.data       = np.array([[ 2235],[  2280],[  2326],[  2361],[  2406],[  2445],[  2485],[  2528],[  2571],[  2607],[  2645]])        # J/kg-K
        self.viscosity.data           = np.array([[0.804],[ 0.704],[ 0.623],[ 0.556],[ 0.498],[ 0.451],[ 0.410],[ 0.374],[ 0.346],[ 0.317],[ 0.289]])        # Pa-s
        self.conductivity.data        = np.array([[0.105],[ 0.104],[ 0.102],[ 0.100],[ 0.098],[ 0.096],[ 0.095],[ 0.093],[ 0.091],[ 0.089],[ 0.087]])        # W/m-K
        self.saturation_pressure.data = np.array([[  0.5],[   0.9],[   1.4],[   2.3],[   3.9],[   6.0],[   8.7],[  12.4],[  17.6],[  24.4],[  33.2]])        # Pa


class SolutionExample(SolutionData):
    def __init__(self):
        SolutionData.__init__(self) 
        self.name = "ExampleSolution"
        self.description = "Ethanol ice slurry"
        self.reference = "SecCool software"
        self.Tmax = -10 + 273.15
        self.Tmin = -45 + 273.15
        self.TminPsat =  self.Tmax
        
        self.temperature.data         = np.array([   -45 ,    -40 ,    -35 ,    -30 ,    -25 ,    -20 ,    -15 ,    -10])+273.15 # Kelvin
        self.concentration.data       = np.array([     5 ,     10 ,     15 ,     20 ,     25 ,     30 ,     35 ])/100.0 # mass fraction
        
        self.density.data             = np.array([
          [1064.0,    1054.6,    1045.3,    1036.3,    1027.4,    1018.6,    1010.0],
          [1061.3,    1052.1,    1043.1,    1034.3,    1025.6,    1017.0,    1008.6],
          [1057.6,    1048.8,    1040.1,    1031.5,    1023.1,    1014.8,    1006.7],
          [1053.1,    1044.6,    1036.2,    1028.0,    1019.9,    1012.0,    1004.1],
          [1047.5,    1039.4,    1031.5,    1023.7,    1016.0,    1008.4,    1000.9],
          [1040.7,    1033.2,    1025.7,    1018.4,    1011.2,    1004.0,     997.0],
          [1032.3,    1025.3,    1018.5,    1011.7,    1005.1,     998.5,     992.0],
          [1021.5,    1015.3,    1009.2,    1003.1,     997.1,     991.2,     985.4]]) # kg/m3
        
#        self.density.data             = np.array([
#          [np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan],
#          [np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan],
#          [np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan],
#          [np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan],
#          [np.nan,    np.nan,    np.nan,    np.nan,    1016.0,    1008.4,    np.nan],
#          [np.nan,    1033.2,    np.nan,    np.nan,    np.nan,    np.nan,    np.nan],
#          [np.nan,    1025.3,    1018.5,    np.nan,    np.nan,     998.5,     992.0],
#          [np.nan,    np.nan,    1009.2,    np.nan,    np.nan,    np.nan,    np.nan]]) # kg/m3

