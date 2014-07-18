from __future__ import division, absolute_import, print_function
import numpy as np
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
        self.data   = None # None #np.zeros((10,10))
        
        self.maxLog = np.log(np.finfo(np.float64).max-1)
        self.minLog = -self.maxLog
        
        self.DEBUG = False
        
        
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
        if not ( (r==3 and c==1) or (r==1 and c==3) ):
            raise ValueError("You have to provide a 3,1 matrix of coefficients, not ({0},{1}).".format(r,c))
        coeffs_tmp = np.array(coeffs.flat)
        return np.exp(np.clip( (coeffs_tmp[0]/ ( x+coeffs_tmp[1] ) - coeffs_tmp[2]),self.minLog,self.maxLog))
    
    
    def baseExponentialOffset(self, c, x):
        raise ValueError("Function not implemented.")
        r,c,coeffs = self.shapeArray(c)
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
    

    @staticmethod
    def shapeArray(array, axs=0):
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
            print(array)
            raise ValueError("You have to provide a 1D-vector or a 2D-matrix.")
        return (r,c,np.reshape(array,(r,c)))
        
        
    def fitCoeffs(self, x=0.0, y=0.0, xbase=0.0, ybase=0.0):
        """ 
        A function that selects the correct equations and
        fits coefficients. Some functions require a start
        guess for the coefficients to work properly.
        """      
        dr,dc,_ = self.shapeArray(self.data)
        xr,xc,x = self.shapeArray(x)
        yr,yc,y = self.shapeArray(y, axs=1)
        
        if self.DEBUG: print("Data        : ({0},{1})".format(dr,dc))
        if self.DEBUG: print("x-axis      : ({0},{1})".format(xr,xc))
        if self.DEBUG: print("y-axis      : ({0},{1})".format(yr,yc))
        
        if dr==1 and dc==1: # 
            if self.DEBUG: print("Data no set, we cannot fit the coefficients")
            self.coeffs = None
            self.type = self.INCOMPRESSIBLE_NOT_SET
            return 
                
        if (xc!=1): raise ValueError("The first input has to be a 2D array with one column.")
        if (yr!=1): raise ValueError("The second input has to be a 2D array with one row.")
        if (xr!=dr): raise ValueError("First independent vector and result vector have to have the same number of rows, {0} is not {1}.".format(xr,dr))
        if (yc!=dc): raise ValueError("Second independent vector and result vector have to have the same number of columns, {0} is not {1}.".format(yc,dc))
        
        cr,cc,_ = self.shapeArray(self.coeffs)
        if self.DEBUG: print("Coefficients: ({0},{1})".format(cr,cc))
        if (xr==1 and xc==1 and cr>1):
            if self.DEBUG: print("Discarding coefficient rows, {0} -> {1}".format(cr,xr))
            self.coeffs = self.coeffs[0]
        if (yr==1 and yc==1 and cc>1):
            if self.DEBUG: print("Discarding coefficient columns, {0} -> {1}".format(cc,yc))
            self.coeffs = self.coeffs.T[0].T
        cr,cc,_ = self.shapeArray(self.coeffs)
                
        if self.DEBUG: print("Coefficients before fitting: \n{0}".format(self.coeffs))
        
        # Polynomial fitting works for both 1D and 2D functions
        if self.type==self.INCOMPRESSIBLE_POLYNOMIAL or self.type==self.INCOMPRESSIBLE_EXPPOLYNOMIAL:
            if self.DEBUG: print("polynomial detected, fitting {0}".format(self.type))
            if cr==1 and cc==1: 
                if self.DEBUG: print("No coefficients left to fit, aborting procedure.")
                self.coeffs = np.array([[0.0]])
                return
            if (xr<cr): raise ValueError("Less data points than coefficients in first dimension ({0} < {1}), aborting.".format(xr,cr))
            if (yc<cc): raise ValueError("Less data points than coefficients in second dimension ({0} < {1}), aborting.".format(yc,cc))
            x_input = np.array(x.flat)-xbase
            y_input = np.array(y.flat)-ybase
            z_input = np.copy(self.data)
            if self.type==self.INCOMPRESSIBLE_EXPPOLYNOMIAL:
                z_input = np.log(z_input)
            self.coeffs = self.getCoeffs2d(x_input, y_input, z_input, cr-1, cc-1)
            if self.DEBUG: print("Coefficients after fitting: \n{0}".format(self.coeffs))
            return
        
        
        # Select if 1D or 2D fitting
        if yc==1: # 1D fitting, only one input
            if self.DEBUG: print("1D function detected, fitting {0}".format(self.type))      
            x_input = x-xbase
            self.coeffs = self.getCoeffsIterative1D(x_input, coeffs_start=None)
            if self.DEBUG: print("Coefficients after fitting: \n{0}".format(self.coeffs))
            return 

        elif yc>1: # 2D fitting
            raise ValueError("There are no other 2D fitting functions than polynomials, cannot use {0}.".format(self.type))
        
        else:  
            raise ValueError("Unknown function.")           
        

#    def getCoeffs1d(self, x, z, order):
#        if (len(x)<order+1): 
#            raise ValueError("You have only {0} elements and try to fit {1} coefficients, please reduce the order.".format(len(x),order+1))
#        A = np.vander(x,order+1)[:,::-1]
#        #Anew = np.dot(A.T,A)
#        #znew = np.dot(A.T,z)
#        #coeffs = np.linalg.solve(Anew, znew)
#        coeffs, resids, rank, singulars  = np.linalg.lstsq(A, z)
#        return np.reshape(coeffs, (len(x),1))
            

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
        if self.DEBUG: print("Linear algebra solver returned:")
        if self.DEBUG: print(coeffs)
        if self.DEBUG: print(resids)
        if self.DEBUG: print(rank)
        if self.DEBUG: print(singulars)
        
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
                    if np.any(popt!=coeffs_start):
                        success = True
#                        print "Fit succeeded for "+fit[counter]+": "
#                        print "data: {0}, func: {1}".format(yData[ 2],func(xData[ 2], popt))
#                        print "data: {0}, func: {1}".format(yData[ 6],func(xData[ 6], popt))
#                        print "data: {0}, func: {1}".format(yData[-1],func(xData[-1], popt))
                        if self.DEBUG: print("Estimated covariance of parameters: ".format(pcov))
                        return popt
                    else:
                        print("Fit failed for "+fit[counter]+": ")
                        success = False 
                    
                except RuntimeError as e:
                    print("Exception using "+fit[counter]+": "+str(e))
                    print("Using "+str(fit[counter+1])+" as a fall-back.")
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
                    print("Fit did not succeed, reducing tolerance to "+str(tol))
                    res = minimize(fun, coeffs_start, method=fit[counter], args=arguments, tol=tol)
                
                # Include these lines for an additional fit with new guess values. 
                #if res.success and tol>tolStart:
                #    print "Refitting with new guesses and original tolerance of "+str(tolStart)
                #    res = minimize(fun, res.x, method=method, args=arguments, tol=tolStart)
                
                if res.success:
                    success = True
                    return res.x
                else:
                    print("Fit failed for "+fit[counter]+": ")
                    print("Using "+str(fit[counter+1])+" as a fall-back.")
                    print(res)
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

            
    