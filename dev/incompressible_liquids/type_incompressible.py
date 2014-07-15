import numpy as np
import itertools,scipy.interpolate
from scipy.optimize._minimize import minimize
from scipy.optimize.minpack import curve_fit

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
        
        self.maxLog = np.log(np.finfo(np.float64).max-1)
        self.minLog = -self.maxLog
        
        
    ### Base functions that handle the custom data type, just a place holder to show the structure.
    def baseFunction(self, x, y=0.0, xbase=0.0, ybase=0.0, c=None):
        if c==None:
            c = self.coeffs
               
        if self.type==self.INCOMPRESSIBLE_POLYNOMIAL:
            return np.polynomial.polynomial.polyval2d(x-xbase, y-ybase, c)
        
        elif self.type==self.INCOMPRESSIBLE_POLYOFFSET:
            if y!=0.0: raise ValueError("This is 1D only, use x not y.")
            return self.basePolyOffset(c, x) # offset included in coeffs
        
        elif self.type==self.INCOMPRESSIBLE_EXPONENTIAL:
            if y!=0.0: raise ValueError("This is 1D only, use x not y.")
            return self.baseExponential(c, x-xbase)
        
        elif self.type==self.INCOMPRESSIBLE_EXPPOLYNOMIAL:
            return np.exp(np.polynomial.polynomial.polyval2d(x-xbase, y-ybase, c))
        
        elif self.type==self.INCOMPRESSIBLE_EXPOFFSET:
            if y!=0.0: raise ValueError("This is 1D only, use x not y.")
            return self.baseExponentialOffset(c, x-xbase)
        
        else:
            raise ValueError("Unknown function: {0}.".format(self.type))

    
    def baseExponential(self, co, x):
        r,c,coeffs = self.shapeArray(co)
        offset = 0.0
        inV    = 0.0
        if not ( (r==3 and c==1) or (r==1 and c==3) ):
            raise ValueError("You have to provide a 3,1 matrix of coefficients, not ({0},{1}).".format(r,c))
        coeffs_tmp = np.array(coeffs.flat)
        return np.exp(np.clip( (coeffs_tmp[0]/ ( x+coeffs_tmp[1] ) - coeffs_tmp[2]),self.minLog,self.maxLog))
    
    
    def baseExponentialOffset(self, c, x):
        raise ValueError("Function not implemented.")
        r,c,coeffs = self.shapeArray(c)
        offset = 0.0
        inV    = 0.0
        if not ( (r==4 and c==1) or (r==1 and c==4) ):
            raise ValueError("You have to provide a 4,1 matrix of coefficients, not ({0},{1}).".format(r,c))
        coeffs_tmp = np.array(coeffs.flat)
        return np.exp(np.clip( (coeffs_tmp[1]/ ( x-coeffs_tmp[0]+coeffs_tmp[2]) - coeffs_tmp[3]),self.minLog,self.maxLog))


    def basePolyOffset(self, co, x):
        r,c,coeffs = self.shapeArray(co)
        if not ( c==1 or r==1 ):
           raise ValueError("You have to provide a 1D vector of coefficients, not ({0},{1}).".format(r,c))
        offset = coeffs[0][0]
        coeffs = np.array(coeffs.flat)[1:c]
        return np.polynomial.polynomial.polyval(x-offset, coeffs)
    

    def shapeArray(self, array, axs=0):
        """
        A function that promotes a 1D array to 2D and 
        also returns the columns and rows.
        1D vectors are interpreted as column vectors.
        """
        r = 0
        c = 0
        array = np.array(array)
        if array.ndim==0:
            array = np.array([[array]])
            r = 1
            c = 1
        elif array.ndim==1:
            if axs==0:
                r = len(array)
                c = 1
            elif axs==1:
                r = 1
                c = len(array)
            else:
                raise ValueError("You have to provide 0 or 1 to the axs parameter, not {0}.".format(axs))
        elif array.ndim==2:
            (r,c) = array.shape
        else:
            print array
            raise ValueError("You have to provide a 1D-vector or a 2D-matrix.")
        return (r,c,np.reshape(array,(r,c)))
        
        
    def fit(self, x, y=0.0, xbase=0.0, ybase=0.0):
        """ 
        A function that selects the correct equations and
        fits coefficients. Some functions require a start
        guess for the coefficients to work properly.
        """
        cr,cc,_ = self.shapeArray(self.coeffs)
        dr,dc,_ = self.shapeArray(self.data)
        xr,xc,x = self.shapeArray(x)
        yr,yc,y = self.shapeArray(y, axs=1)
                
        if (xc!=1): raise ValueError("The first input has to be a 2D array with one column.")
        if (yr!=1): raise ValueError("The second input has to be a 2D array with one row.")
        if (xr!=dr): raise ValueError("First independent vector and result vector have to have the same number of rows, {0} is not {1}.".format(xr,dr))
        if (yc!=dc): raise ValueError("Second iIndependent vector and result vector have to have the same number of columns, {0} is not {1}.".format(yc,dc))
        
        # Polynomial fitting works for both 1D and 2D functions
        if self.type==self.INCOMPRESSIBLE_POLYNOMIAL or \
           self.type==self.INCOMPRESSIBLE_EXPPOLYNOMIAL:
            print "Polynomial detected, fitting {0}".format(self.type)
            if (xr<cr): raise ValueError("Less data points than coefficients in first dimension ({0} < {1}), aborting.".format(xr,cr))
            if (yc<cc): raise ValueError("Less data points than coefficients in second dimension ({0} < {1}), aborting.".format(yc,cc))
            x_input = x-xbase
            y_input = y-ybase
            z_input = np.copy(self.data)
            if self.type==self.INCOMPRESSIBLE_EXPPOLYNOMIAL:
                z_input = np.log(z_input)
            self.coeffs = self.getCoeffs2d(x_input, y_input, z_input, cr-1, cc-1)

        
        # Select if 1D or 2D fitting
        if yc==1: # 1D fitting, only one input
            print "1D function detected, fitting {0}".format(self.type)
            x_input = x-xbase
            print self.getCoeffsIterative1D(x_input, coeffs_start=None)

        elif yc>1: # 2D fitting
            raise ValueError("There are no other 2D fitting functions than polynomials, cannot use {0}.".format(self.type))
        
        else:  
            raise ValueError("Unknown function.")           
        

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
    
    
    def getCoeffsIterative1D(self, xData, coeffs_start=None):
        #fit = "Powell" # use Powell's algorithm
        #fit = "BFGS" # use Broyden-Fletcher-Goldfarb-Shanno
        #fit = "LMA" # use the Levenberg-Marquardt algorithm from curve_fit 
        fit  = ["LMA","Powell","BFGS"] # First try LMA, use others as fall-back

        # Basic settings to keep track of our efforts
        success = False
        counter = -1
        
        if coeffs_start==None and \
          self.coeffs!=None:
            coeffs_start = self.coeffs
        
        # make sure that we use other routines for polynomials
        if (self.type==self.INCOMPRESSIBLE_POLYNOMIAL) or \
           (self.type==self.INCOMPRESSIBLE_EXPPOLYNOMIAL) :
            raise ValueError("Please use the specific polynomial functions, they are much better.")
        
        expLog = False
        # Preprocess the exponential data
        if (self.type==self.INCOMPRESSIBLE_EXPONENTIAL) or \
           (self.type==self.INCOMPRESSIBLE_EXPOFFSET) :
            expLog = True
        
        xData = np.array(xData.flat)
        if expLog:
            yData = np.log(self.data.flat)
        else:
            yData = np.array(self.data.flat)
            
        # The residual function            
        def fun(coefficients,xArray,yArray):
            """ 
            This function takes the coefficient array and
            x and y data. It evaluates the function and returns
            the sum of the squared residuals if yArray is not
            equal to None
            """
            calculated = self.baseFunction(xArray, 0.0, 0.0, 0.0, coefficients)
            if expLog: calculated = np.log(calculated)
            if yArray==None: return calculated
            data       = yArray
            res        = np.sum(np.square(calculated-data))           
            return res
        
        # Loop through the list of algorithms
        while (not success):
            counter += 1
            
            #fit = "LMA" # use the Levenberg-Marquardt algorithm from curve_fit 
            if fit[counter]=="LMA": 
                
                def func(xVector, *coefficients):
                    return fun(np.array(coefficients), xVector, None)
                    #return self.baseFunction(xVector, 0.0, 0.0, 0.0, np.array(coefficients)) #np.array([self._PropsFit(coefficients,xName,T=Ti) for Ti in T])
                
                try:
                    #print func(xData, coeffs_start)
                    # Do the actual fitting
                    popt, pcov = curve_fit(func, xData, yData, p0=coeffs_start)
                    #print popt 
                    #print pcov
                    if np.any(popt!=coeffs_start):
                        success = True
#                        print "Fit succeeded for "+fit[counter]+": "
#                        print "data: {0}, func: {1}".format(yData[ 2],func(xData[ 2], popt))
#                        print "data: {0}, func: {1}".format(yData[ 6],func(xData[ 6], popt))
#                        print "data: {0}, func: {1}".format(yData[-1],func(xData[-1], popt))
                        return popt
                    else:
                        print "Fit failed for "+fit[counter]+": "
                        success = False 
                    
                except RuntimeError as e:
                    print "Exception using "+fit[counter]+": "+str(e)
                    print "Using "+str(fit[counter+1])+" as a fall-back."
                    success = False
        
            #fit = "MIN" # use a home-made minimisation with Powell and Broyden-Fletcher-Goldfarb-Shanno
            elif fit[counter]=="Powell" or fit[counter]=="BFGS":
                
                arguments  = (xData,yData)
                #options    = {'maxiter': 1e2, 'maxfev': 1e5}
                
                tolStart   = 1e-13
                tol        = tolStart
                res = minimize(fun, coeffs_start, method=fit[counter], args=arguments, tol=tol)
                
                while ((not res.success) and tol<1e-6):
                    tol *= 1e2
                    print "Fit did not succeed, reducing tolerance to "+str(tol)
                    res = minimize(fun, coeffs_start, method=fit[counter], args=arguments, tol=tol)
                
                # Include these lines for an additional fit with new guess values. 
                #if res.success and tol>tolStart:
                #    print "Refitting with new guesses and original tolerance of "+str(tolStart)
                #    res = minimize(fun, res.x, method=method, args=arguments, tol=tolStart)
                
                if res.success:
                    success = True
                    return res.x
                else:
                    print "Fit failed for "+fit[counter]+": "
                    print "Using "+str(fit[counter+1])+" as a fall-back."
                    print res
                    success = False
                    
            # Something went wrong, probably a typo in the algorithm selector
            else:
                raise (ValueError("Error: You used an unknown fit method."))
        
        
        
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
        
        # Some of the functions might need a guess array
        self.viscosity.type = self.viscosity.INCOMPRESSIBLE_EXPONENTIAL
        self.viscosity.coeffs = np.array([+7e+2, -6e+1, +1e+1])
        self.saturation_pressure.type = self.saturation_pressure.INCOMPRESSIBLE_EXPONENTIAL
        self.saturation_pressure.coeffs = np.array([-5e+3, +3e+1, -1e+1])
        
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

if __name__ == '__main__':
    obj = PureExample()
    obj.saturation_pressure.type = obj.saturation_pressure.INCOMPRESSIBLE_EXPONENTIAL
    print obj.saturation_pressure.coeffs
    xData = obj.temperature.data
    coeffs = obj.saturation_pressure.getCoeffsIterative1D(xData)
    obj.saturation_pressure.coeffs = coeffs
    print obj.saturation_pressure.coeffs
    
    obj = PureExample()
    print obj.saturation_pressure.coeffs
    obj.saturation_pressure.fit(obj.temperature.data)
    print obj.saturation_pressure.coeffs
