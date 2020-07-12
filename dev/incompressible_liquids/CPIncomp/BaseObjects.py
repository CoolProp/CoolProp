from __future__ import division, print_function
import numpy as np
from scipy.optimize._minimize import minimize
from scipy.optimize.minpack import curve_fit
import sys

# Here we define the types. This is done to keep the definitions at one
# central place instead of hiding them somewhere in the data.


class IncompressibleData(object):
    """
    The work horse of the incompressible classes.
    Implements both data structures and fitting
    procedures.
    """
    INCOMPRESSIBLE_NOT_SET = 'notdefined'
    INCOMPRESSIBLE_POLYNOMIAL = 'polynomial'
    INCOMPRESSIBLE_EXPPOLYNOMIAL = 'exppolynomial'
    INCOMPRESSIBLE_EXPONENTIAL = 'exponential'
    INCOMPRESSIBLE_LOGEXPONENTIAL = 'logexponential'
    INCOMPRESSIBLE_POLYOFFSET = 'polyoffset'
    INCOMPRESSIBLE_CHEBYSHEV = 'chebyshev'

    SOURCE_DATA = 'data'
    SOURCE_EQUATION = 'equation'
    SOURCE_COEFFS = 'coefficients'
    SOURCE_NOT_SET = 'notdefined'

    maxLin = np.finfo(np.float64).max - 1
    minLin = -maxLin

    maxLog = np.log(maxLin)
    minLog = -maxLog

    def __init__(self):
        self.source = self.SOURCE_NOT_SET
        self.type = self.INCOMPRESSIBLE_NOT_SET
        self.coeffs = None  # np.zeros((4,4))
        self.data = None  # None #np.zeros((10,10))
        self.xData = None  # In case you need a customised first data set (temperature?)
        self.yData = None  # In case you need a customised second data set (concentration?)
        self.sErr = None  # Coefficient of determination
        self.NRMS = None  # Normalised RMS: (square root of the mean of the squares of the deviations)/(max-min)
        self.DEBUG = False

    @staticmethod
    def baseFunc(x, y=0.0, xbase=0.0, ybase=0.0, eqnType=None, c=None):

        if eqnType is None: raise ValueError("You did not provide data for eqnType.")
        if c is None: raise ValueError("You did not provide data for the coefficients.")

        if eqnType == IncompressibleData.INCOMPRESSIBLE_POLYNOMIAL:
            return np.polynomial.polynomial.polyval2d(x - xbase, y - ybase, c)

        elif eqnType == IncompressibleData.INCOMPRESSIBLE_POLYOFFSET:
            #if y!=0.0: raise ValueError("This is 1D only, use x not y.")
            return IncompressibleData.basePolyOffset(c, x)  # offset included in coeffs

        elif eqnType == IncompressibleData.INCOMPRESSIBLE_EXPONENTIAL:
            #if y!=0.0: raise ValueError("This is 1D only, use x not y.")
            return IncompressibleData.baseExponential(c, x)

        elif eqnType == IncompressibleData.INCOMPRESSIBLE_LOGEXPONENTIAL:
            #if y!=0.0: raise ValueError("This is 1D only, use x not y.")
            return IncompressibleData.baseLogexponential(c, x)

        elif eqnType == IncompressibleData.INCOMPRESSIBLE_EXPPOLYNOMIAL:
            return np.exp(np.polynomial.polynomial.polyval2d(x - xbase, y - ybase, c))

        else:
            raise ValueError("Unknown function: {0}.".format(eqnType))

    @staticmethod
    def baseExponential(co, x):
        r, c, coeffs = IncompressibleFitter.shapeArray(co)
        if not ((r == 3 and c == 1) or (r == 1 and c == 3)):
            raise ValueError("You have to provide a 3,1 matrix of coefficients, not ({0},{1}).".format(r, c))
        coeffs_tmp = np.array(coeffs.flat)
        return np.exp(np.clip((coeffs_tmp[0] / (x + coeffs_tmp[1]) - coeffs_tmp[2]), IncompressibleData.minLog, IncompressibleData.maxLog))

    @staticmethod
    def baseLogexponential(co, x):
        r, c, coeffs = IncompressibleFitter.shapeArray(co)
        if not ((r == 3 and c == 1) or (r == 1 and c == 3)):
            raise ValueError("You have to provide a 4,1 matrix of coefficients, not ({0},{1}).".format(r, c))
        coeffs_tmp = np.array(coeffs.flat)
        return np.exp(np.clip(np.log(np.clip(1.0 / (x + coeffs_tmp[0]) + 1.0 / (x + coeffs_tmp[0])**2, 1e-10, IncompressibleData.maxLin)) * coeffs_tmp[1] + coeffs_tmp[2], IncompressibleData.minLog, IncompressibleData.maxLog))

    @staticmethod
    def basePolyOffset(co, x):
        r, c, coeffs = IncompressibleFitter.shapeArray(co)
        if not (c == 1 or r == 1):
            raise ValueError("You have to provide a 1D vector of coefficients, not ({0},{1}).".format(r, c))
        offset = coeffs[0][0]
        coeffs = np.array(coeffs.flat)[1:]
        return np.polynomial.polynomial.polyval(x - offset, coeffs)

        # Base functions that handle the custom data type, just a place holder to show the structure.
    def baseFunction(self, x, y=0.0, xbase=0.0, ybase=0.0, c=None):
        if c is None: c = self.coeffs
        return self.baseFunc(x, y, xbase, ybase, self.type, c)

    def fitCoeffs(self, xbase, ybase, x=None, y=None):

        if (not x is None and not self.xData is None and not IncompressibleFitter.allClose(x, self.xData)) \
        or (x is None and self.xData is None): raise ValueError("I do not know which x-value you would like to use. Define either x or self.xData.")
        if (not y is None and not self.yData is None and not IncompressibleFitter.allClose(y, self.yData)) \
        or (y is None and self.yData is None): raise ValueError("I do not know which y-value you would like to use. Define either y or self.yData.")

        if x is None and not self.xData is None: x = self.xData
        if y is None and not self.yData is None: y = self.yData

        #res = None
        #r2 = None

        res, sErr = IncompressibleFitter.fitter(x=x, y=y, z=self.data, \
                  xbase=xbase, ybase=ybase, \
                  eqnType=self.type, \
                  coeffs=self.coeffs, DEBUG=self.DEBUG)

        bestCoeffs = np.copy(res)
        bestType = self.type
        bestsErr = np.copy(sErr)
        bestRMS = np.sqrt(np.square(bestsErr).mean()).sum()

        count = 0
        while bestRMS > 0.03 and count < 2:
            #if self.DEBUG: print("Poor solution found, trying once more with more coefficients.")
            if self.type == IncompressibleData.INCOMPRESSIBLE_EXPONENTIAL:
                if self.DEBUG: print("Poor solution found with exponential, trying once more with log exponential.")
                self.type = IncompressibleData.INCOMPRESSIBLE_LOGEXPONENTIAL
                self.coeffs = np.array([-250, 1.5, 10])
                res, sErr = IncompressibleFitter.fitter(x=x, y=y, z=self.data, \
                      xbase=xbase, ybase=ybase, \
                      eqnType=self.type, \
                      coeffs=self.coeffs, DEBUG=self.DEBUG)

            elif self.type == IncompressibleData.INCOMPRESSIBLE_LOGEXPONENTIAL:
                xLen = np.round([len(x) / 1.5])
                yLen = np.round([len(y) / 1.5])
                xLen = int(np.min([xLen, 4]))
                yLen = int(np.min([yLen, 6]))

                if (xLen + yLen) > 2:
                    if self.DEBUG: print("Poor solution found with log exponential, trying once more with exponential polynomial.")
                    self.type = IncompressibleData.INCOMPRESSIBLE_EXPPOLYNOMIAL

                    self.coeffs = np.zeros((xLen, yLen))
                    res, sErr = IncompressibleFitter.fitter(x=x, y=y, z=self.data, \
                          xbase=xbase, ybase=ybase, \
                          eqnType=self.type, \
                          coeffs=self.coeffs, DEBUG=self.DEBUG)

#            elif self.type==IncompressibleData.INCOMPRESSIBLE_EXPPOLYNOMIAL:
#                if self.DEBUG: print("Poor solution found with exponential polynomial, trying once more with normal polynomial.")
#                self.type=IncompressibleData.INCOMPRESSIBLE_POLYNOMIAL
#                self.coeffs = np.zeros((4,6))
#                res,sErr = IncompressibleFitter.fitter(x=x, y=y, z=self.data, \
#                      xbase=xbase, ybase=ybase, \
#                      eqnType=self.type, \
#                      coeffs=self.coeffs, DEBUG=self.DEBUG)

            RMS = np.sqrt(np.square(sErr).mean()).sum()
            if RMS < bestRMS:  # Better fit
                bestCoeffs = np.copy(res)
                bestType = self.type
                bestsErr = np.copy(sErr)
                bestRMS = RMS

            count += 1

        if bestCoeffs is None:
            if self.DEBUG: print("There was a fitting error, no solution found.")
        elif IncompressibleFitter.allClose(bestCoeffs, self.coeffs):
            if self.DEBUG: print("Coefficients did not change.")
        else:
            if self.DEBUG: print("Best fit for: {0}".format(bestType))
            self.coeffs = bestCoeffs
            self.type = bestType
            self.sErr = bestsErr
            if self.type == IncompressibleData.INCOMPRESSIBLE_POLYNOMIAL:
                self.NRMS = bestRMS / (np.nanmax(self.data) - np.nanmin(self.data))
            elif self.type == IncompressibleData.INCOMPRESSIBLE_EXPONENTIAL or \
              self.type == IncompressibleData.INCOMPRESSIBLE_EXPPOLYNOMIAL or \
              self.type == IncompressibleData.INCOMPRESSIBLE_LOGEXPONENTIAL:
                self.NRMS = bestRMS / (np.log(np.nanmax(self.data)) - np.log(np.nanmin(self.data)))
            else:
                raise ValueError("Unknown function.")

            #if self.DEBUG: print("Fitting statistics:")
            # SSE = np.square(self.sErr).sum() # Sum of squares due to error
            #SST = ((zData-zData.mean())**2).sum()
            #R2  = 1-(ssErr/ssTot )

    def setxData(self, xData):
        if self.xData is None:
            self.xData = xData
        else:
            if self.DEBUG: print("Cannot update xData, value is set already.")

    def setyData(self, yData):
        if self.yData is None:
            self.yData = yData
        else:
            if self.DEBUG: print("Cannot update yData, value is set already.")

    def setxyData(self, xData, yData):
        self.setxData(xData)
        self.setyData(yData)

    def toJSON(self):
        j = {}
        try:
            j['coeffs'] = self.coeffs.tolist()
        except:
            j['coeffs'] = 'null'

        try:
            j['NRMS'] = self.NRMS
        except:
            j['NRMS'] = 'null'

        j['type'] = self.type
        return j

    def fromJSON(self, j):
        try:
            self.coeffs = np.array(j['coeffs'])
            self.type = j['type']
        except:
            self.coeffs = None
            self.type = IncompressibleData.INCOMPRESSIBLE_NOT_SET
        try:
            self.NRMS = j['NRMS']
        except:
            self.NRMS = None
        return


class IncompressibleFitter(object):

    @staticmethod
    def allClose(a, b):
        if np.array(a).shape == np.array(b).shape:
            return np.allclose(a, b)
        else: return False

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
        if array.ndim == 0:
            array = np.array([[array]])
            r = 1
            c = 1
        elif array.ndim == 1:
            if axs == 0:
                r = len(array)
                c = 1
            elif axs == 1:
                r = 1
                c = len(array)
            else:
                raise ValueError("You have to provide 0 or 1 to the axs parameter, not {0}.".format(axs))
        elif array.ndim == 2:
            (r, c) = array.shape
        else:
            print(array)
            raise ValueError("You have to provide a 1D-vector or a 2D-matrix.")
        return (r, c, np.reshape(array, (r, c)))

    @staticmethod
    def fitter(x=None, y=None, z=None, \
                  xbase=0.0, ybase=0.0, \
                  eqnType=None, \
                  coeffs=None, DEBUG=False):
        """ The entry point to the fitting routines

        :param x       : a 1D array in x direction or 2D with one column, most likely temperature
        :param y       : a 1D array in y direction or 2D with one row, most likely cocentration
        :param z       : a 2D array of data, rows = len(x[:,0]) and cols = len(y[0])
        :param xbase   : a value to be subtracted from x, might not be used
        :param ybase   : a value to be subtracted from y, might not be used
        :param eqnType : an instance of IncompressibleData.INCOMPRESSIBLE_ ...
        :param coeffs  : the initial guess and shape (!) for the desired coefficients, can be zeros
        :param DEBUG   : message to display

        :returns       : None if failed or coeffs filled with the right values


        A function that selects the correct equations and
        fits coefficients. Some functions require a start
        guess for the coefficients to work properly.
        """

        if x is None: raise ValueError("You did not provide data for the x-values.")
        if y is None: raise ValueError("You did not provide data for the y-values.")
        if z is None: raise ValueError("You did not provide data for the z-values.")
        if xbase is None: raise ValueError("You did not provide data for xbase.")
        if ybase is None: raise ValueError("You did not provide data for ybase.")
        if eqnType is None: raise ValueError("You did not provide data for eqnType.")
        if coeffs is None: raise ValueError("You did not provide data for the coefficients.")
        if DEBUG is None: raise ValueError("You did not provide data for DEBUG.")

        zr, zc, _ = IncompressibleFitter.shapeArray(z)
        xr, xc, x = IncompressibleFitter.shapeArray(x)
        yr, yc, y = IncompressibleFitter.shapeArray(y, axs=1)

        if DEBUG: print("Data        : ({0},{1})".format(zr, zc))
        if DEBUG: print("x-axis      : ({0},{1})".format(xr, xc))
        if DEBUG: print("y-axis      : ({0},{1})".format(yr, yc))

        if zr == 1 and zc == 1:
            if DEBUG: print("Data no set, we cannot fit the coefficients")
            return None, None

        if (xc != 1): raise ValueError("The first input has to be a 2D array with one column.")
        if (yr != 1): raise ValueError("The second input has to be a 2D array with one row.")
        if (xr != zr): raise ValueError("First independent vector and result vector have to have the same number of rows, {0} is not {1}.".format(xr, zr))
        if (yc != zc): raise ValueError("Second independent vector and result vector have to have the same number of columns, {0} is not {1}.".format(yc, zc))

        if DEBUG: print("Coefficients before fitting: \n{0}".format(coeffs))

        # Polynomial fitting works for both 1D and 2D functions
        if eqnType == IncompressibleData.INCOMPRESSIBLE_POLYNOMIAL or eqnType == IncompressibleData.INCOMPRESSIBLE_EXPPOLYNOMIAL:

            cr, cc, _ = IncompressibleFitter.shapeArray(coeffs)
            if DEBUG: print("Coefficients: ({0},{1})".format(cr, cc))
            if (xr == 1 and xc == 1 and cr > 1):
                if DEBUG: print("Discarding coefficient rows, {0} -> {1}".format(cr, xr))
                coeffs = coeffs[0]
                coeffs = coeffs.reshape((1, cc))
            if (yr == 1 and yc == 1 and cc > 1):
                if DEBUG: print("Discarding coefficient columns, {0} -> {1}".format(cc, yc))
                coeffs = coeffs.T[0]
                coeffs = coeffs.reshape((cr, 1))
            cr, cc, _ = IncompressibleFitter.shapeArray(coeffs)

            if DEBUG: print("polynomial detected, fitting {0}".format(eqnType))
            if cr == 1 and cc == 1:
                if DEBUG: print("No coefficients left to fit, aborting procedure.")
                coeffs = np.array([[0.0]])
                return coeffs, None
            if (xr < cr):
                if DEBUG: print("Less data points than coefficients in first dimension ({0} < {1}), reducing coefficient matrix.".format(xr, cr))
                coeffs = coeffs[:xr, :]
            if (yc < cc):
                if DEBUG: print("Less data points than coefficients in second dimension ({0} < {1}), reducing coefficient matrix.".format(yc, cc))
                coeffs = coeffs[:, :yc]
            cr, cc, _ = IncompressibleFitter.shapeArray(coeffs)
            x_input = np.array(x.flat) - xbase
            y_input = np.array(y.flat) - ybase
            z_input = np.copy(z)
            if eqnType == IncompressibleData.INCOMPRESSIBLE_EXPPOLYNOMIAL:
                z_input = np.log(z_input)
            coeffs, sErr = IncompressibleFitter.getCoeffs2d(x_input, y_input, z_input, cr - 1, cc - 1, DEBUG=DEBUG)
            if DEBUG: print("Coefficients after fitting: \n{0}".format(coeffs))
            if DEBUG: print("Standard deviation: {0}".format(np.nanstd(sErr)))
            if DEBUG: print("Sum of squared errors: {0}".format(np.square(sErr).sum()))
            if DEBUG: print("Root mean squared errors: {0}".format(np.sqrt(np.square(sErr).mean()).sum()))
            return coeffs, sErr

        # Select if 1D or 2D fitting
        if yc == 1 or xr == 1:  # 1D fitting, only one input
            if DEBUG: print("1D function detected.")
            if yc == 1:
                if DEBUG: print("Fitting {0} in x-direction.".format(eqnType))
                coeffs, sErr = IncompressibleFitter.getCoeffsIterative1D(x, z, eqnType=eqnType, coeffs=coeffs, DEBUG=DEBUG)
            elif xr == 1:
                if DEBUG: print("Fitting {0} in y-direction.".format(eqnType))
                coeffs, sErr = IncompressibleFitter.getCoeffsIterative1D(y, z, eqnType=eqnType, coeffs=coeffs, DEBUG=DEBUG)
            else: raise ValueError("Unknown error in matrix shapes.")
            if DEBUG: print("Coefficients after fitting: \n{0}".format(coeffs))
            if DEBUG: print("Standard deviation:    {0}".format(np.nanstd(sErr)))
            if DEBUG: print("Sum of squared errors: {0}".format(np.square(sErr).sum()))
            if DEBUG: print("Root mean squared errors: {0}".format(np.sqrt(np.square(sErr).mean()).sum()))
            return coeffs, sErr

        elif yc > 1:  # 2D fitting
            raise ValueError("There are no other 2D fitting functions than polynomials, cannot use {0}.".format(eqnType))

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

    @staticmethod
    def getCoeffs2d(x_in, y_in, z_in, x_order, y_order, DEBUG=False):

        if x_in is None: raise ValueError("You did not provide data for the x-values.")
        if y_in is None: raise ValueError("You did not provide data for the y-values.")
        if z_in is None: raise ValueError("You did not provide data for the z-values.")
        if x_order is None: raise ValueError("You did not provide data for x_order.")
        if y_order is None: raise ValueError("You did not provide data for y_order.")
        if DEBUG is None: raise ValueError("You did not provide data for DEBUG.")

        x_order += 1
        y_order += 1
        # To avoid overfitting, we only use the upper left triangle of the coefficient matrix
        x_exp = range(x_order)
        y_exp = range(y_order)
        limit = max(x_order, y_order)
        xy_exp = []
        # Construct the upper left triangle of coefficients
        for i in x_exp:
            for j in y_exp:
                if(i + j < limit): xy_exp.append((i, j))

        # Construct input pairs
        xx, yy = np.meshgrid(x_in, y_in, indexing='ij')
        xx = np.array(xx.flat)
        yy = np.array(yy.flat)
        zz = np.array(z_in.flat)

        # TODO: Check for rows with only nan values
        x_num = len(x_in)
        y_num = len(y_in)

        cols = len(xy_exp)
        eqns = x_num * y_num
        # if (eqns<cols):
        #    raise ValueError("You have only {0} equations and try to fit {1} coefficients, please reduce the order.".format(eqns,cols))
        if (x_num < x_order):
            raise ValueError("You have only {0} x-entries and try to fit {1} x-coefficients, please reduce the x_order.".format(x_num, x_order))
        if (y_num < y_order):
            raise ValueError("You have only {0} y-entries and try to fit {1} y-coefficients, please reduce the y_order.".format(y_num, y_order))

        # Build the functional matrix
        A = np.zeros((eqns, cols))
        for i in range(eqns):  # row loop
            for j, (xj, yj) in enumerate(xy_exp):  # makes columns
                A[i][j] = xx[i]**xj * yy[i]**yj

        # Remove np.nan elements
        mask = np.isfinite(zz)
        A = A[mask]
        xx = xx[mask]
        yy = yy[mask]
        zz = zz[mask]

        if (len(A) < cols):
            raise ValueError("Your matrix has only {0} valid rows and you try to fit {1} coefficients, please reduce the order.".format(len(A), cols))

        coeffs, resids, rank, singulars = np.linalg.lstsq(A, zz, rcond=None)
        if DEBUG: print("Linear algebra solver returned:")
        if DEBUG: print(coeffs)
        if DEBUG: print(resids)
        if DEBUG: print(rank)
        if DEBUG: print(singulars)

        # if resids.size>0:
        #    r2 = 1 - resids / (zz.size * zz.var())
        # else:
        #    r2 = 0
        #print("\n r2 2d: ",r2.shape,r2,"\n")

        # Rearrange coefficients to a matrix shape
        C = np.zeros((x_order, y_order))
        for i, (xi, yi) in enumerate(xy_exp):  # makes columns
            C[xi][yi] = coeffs[i]

        sErr = zz - np.polynomial.polynomial.polyval2d(xx, yy, C)
        return C, sErr

    @staticmethod
    def getCoeffsIterative1D(x_in, z_in, eqnType, coeffs, DEBUG=False):

        if x_in is None: raise ValueError("You did not provide data for the x-values.")
        if z_in is None: raise ValueError("You did not provide data for the z-values.")
        if eqnType is None: raise ValueError("You did not provide data for eqnType.")
        if coeffs is None: raise ValueError("You did not provide data for the coefficients.")
        if DEBUG is None: raise ValueError("You did not provide data for DEBUG.")

        sErr = None

        # fit = "Powell" # use Powell's algorithm
        # fit = "BFGS" # use Broyden-Fletcher-Goldfarb-Shanno
        # fit = "LMA" # use the Levenberg-Marquardt algorithm from curve_fit
        fit = ["LMA", "Powell", "BFGS"]  # First try LMA, use others as fall-back

        # make sure that we use other routines for polynomials
        if (eqnType == IncompressibleData.INCOMPRESSIBLE_POLYNOMIAL) or \
           (eqnType == IncompressibleData.INCOMPRESSIBLE_EXPPOLYNOMIAL):
            raise ValueError("Please use the specific polynomial functions, they are much better.")

        expLog = False
        # Fitting the logarithm of z_in?
        if (eqnType == IncompressibleData.INCOMPRESSIBLE_EXPONENTIAL or eqnType == IncompressibleData.INCOMPRESSIBLE_LOGEXPONENTIAL):
            expLog = True

        xData = np.array(x_in.flat)
        if expLog: zData = np.log(np.clip(z_in.flat, 1e-10, IncompressibleData.maxLin))
        else: zData = np.array(z_in.flat)

        # Remove np.nan elements
        mask = np.isfinite(zData)
        xData = xData[mask]
        zData = zData[mask]

        # The residual function
        def fun(coefficients, xArray, yArray):
            """
            This function takes the coefficient array and
            x and y data. It evaluates the function and returns
            the sum of the squared residuals if yArray is not
            equal to None
            """
            # No offset and no Tbase etc for 1D functions!
            calculated = IncompressibleData.baseFunc(xArray, y=0.0, xbase=0.0, ybase=0.0, eqnType=eqnType, c=coefficients)
            if expLog: calculated = np.log(calculated)
            if yArray is None: return calculated
            data = yArray
            res = np.sum(np.square(calculated - data))
            return res

        # Loop through the list of algorithms with basic settings to keep track of our efforts
        success = False
        counter = 0
        tolerance = 1e-16
        while (not success):
            algorithm = fit[counter]
            # fit = "LMA" # use the Levenberg-Marquardt algorithm from curve_fit
            if algorithm == "LMA":

                def func(xVector, *coefficients):
                    return fun(np.array(coefficients), xVector, None)
                    # return self.baseFunction(xVector, 0.0, 0.0, 0.0, np.array(coefficients)) #np.array([self._PropsFit(coefficients,xName,T=Ti) for Ti in T])

                try:
                    # print func(xData, coeffs_start)
                    # Do the actual fitting
                    popt, pcov = curve_fit(func, xData, zData, p0=coeffs, ftol=tolerance)
                    if np.any(popt != coeffs):
                        success = True
                        if DEBUG: print("Fit succeeded with: {0}".format(algorithm))
                        sErr = zData - func(xData, popt)
                        # print "Fit succeeded for "+fit[counter]+": "
                        # print "data: {0}, func: {1}".format(yData[ 2],func(xData[ 2], popt))
                        # print "data: {0}, func: {1}".format(yData[ 6],func(xData[ 6], popt))
                        # print "data: {0}, func: {1}".format(yData[-1],func(xData[-1], popt))
                        #if DEBUG: print("Estimated covariance of parameters: {0}".format(pcov))
                        #ssErr = np.sqrt(np.diag(pcov)).sum()
                        #ssTot = ((zData-zData.mean())**2).sum()
                        #r2 = 1-(ssErr/ssTot )
                        #print("\n r2 FMA: ",r2.shape,r2,"\n")
                        return popt, sErr
                    else:
                        if DEBUG: print("Fit failed for {0}.".format(algorithm))
                        if DEBUG: sys.stdout.flush()
                        success = False

                except RuntimeError as e:
                    if DEBUG: print("Exception using " + algorithm + ": " + str(e))
                    if DEBUG: sys.stdout.flush()
                    success = False

            # fit = "MIN" # use a home-made minimisation with Powell and Broyden-Fletcher-Goldfarb-Shanno
            elif algorithm == "Powell" or algorithm == "BFGS":

                arguments = (xData, zData)
                #options    = {'maxiter': 1e2, 'maxfev': 1e5}

                try:
                    res = minimize(fun, coeffs, method=algorithm, args=arguments, tol=tolerance)
                    if res.success:
                        success = True
                        if DEBUG: print("Fit succeeded with: {0}".format(algorithm))
                        sErr = zData - fun(np.array(res.x), xData, None)
                        # if res.has_key('fvec'):
                            #ssErr = (res['fvec']**2).sum()
                            #ssTot = ((zData-zData.mean())**2).sum()
                            #r2 = 1-(ssErr/ssTot )
                        #print("\n r2 : ",r2.shape,r2,algorithm,"\n")
                        return res.x, sErr
                    else:
                        if DEBUG: print("Fit failed for {0}.".format(algorithm))
                        if DEBUG: sys.stdout.flush()
                        success = False
                except RuntimeError as e:
                    if DEBUG: print("Exception using " + algorithm + ": " + str(e))
                    if DEBUG: sys.stdout.flush()
                    success = False

            # Something went wrong, probably a typo in the algorithm selector
            else:
                raise (ValueError("Error: You used an unknown fit method."))

            if counter < len(fit) - 1:
                #print("Fit did not succeed with {0}, reducing tolerance to {1}.".format(algorithm,tol))
                success = False
                counter += 1
            elif tolerance < 1e-3:
                tolerance *= 1e2
                if DEBUG: print("Fit did not succeed, reducing tolerance to {0}.".format(tolerance))
                success = False
                counter = 0
            else:
                if DEBUG: print("--------------------------------------------------------------")
                if DEBUG: print("Fit failed for {0}. ".format(fit))
                if DEBUG: print("--------------------------------------------------------------")
                return coeffs, 1
