#cython: embedsignature = True

cdef extern from "CoolProp.h":
    double _Props "Props" (char* ,char,double,char,double,char*)
    double _Props1 "Props" (char*,char*)
    void UseSinglePhaseLUT(bool)
    double DerivTerms(char *, double, double, char*)
    char * get_errstringc()
    void get_errstring(char*)
    void _debug "debug" (int)
    double _IProps "IProps" (long,long,double,long,double,long)
    long _get_param_index "get_param_index" (string param)
    long _get_Fluid_index "get_Fluid_index" (string param)
    void _set_phase "set_phase"(string Phase_str)
   
from libc.math cimport pow, sin, cos, exp
from math import pow as pow_
cdef bint _LUT_Enabled
import CoolProp as CP
cimport cython

## Variables for each of the possible variables
cdef long iMM = _get_param_index('M')
cdef long iT = _get_param_index('T')
cdef long iD = _get_param_index('D')
cdef long iH = _get_param_index('H')
cdef long iP = _get_param_index('P')
cdef long iC = _get_param_index('C')
cdef long iC0 = _get_param_index('C0')
cdef long iO = _get_param_index('O')
cdef long iV = _get_param_index('V')
cdef long iL = _get_param_index('L')
cdef long iS = _get_param_index('S')
cdef long iU = _get_param_index('U')
cdef long iDpdT = _get_param_index('dpdT')

#A dictionary mapping parameter index to string for use with non-CoolProp fluids
cdef dict paras = {iMM : 'M',
                   iT : 'T',
                   iH : 'H',
                   iP : 'P',
                   iC : 'C',
                   iC0 : 'C0',
                   iO : 'O',
                   iV : 'V',
                   iL : 'L',
                   iS : 'S',
                   iU : 'U',
                   iDpdT : 'dpdT'}

cpdef int debug(int level):
    """
    Sets the debug level
    
    Parameters
    ----------
    level : int
        Flag indicating how verbose the debugging should be.
            0 : no debugging output
            ...
            ...
            10 : very annoying debugging output - every function call debugged
    """
    _debug(level)

cpdef int LUT(bint LUTval):
    """
    
    LUTval : boolean
        If ``True``, turn on the use of lookup table.  Parameters must have 
        already been set through the use of set_1phase_LUT_params

    """
    if LUTval:
        _LUT_Enabled = True
        print 'Turning on singlephase LUT'
        UseSinglePhaseLUT(True)
    else:
        _LUT_Enabled = False
        UseSinglePhaseLUT(False)
        
cpdef double Props(bytes Parameter, bytes param1, float value1, bytes param2, float value2, bytes Fluid):
    """
    Expose the Props() function.  Uses the same call signature as the Props() function in CoolProp.CoolProp
    """
    cdef char _param1 = param1[0]
    cdef char _param2 = param2[0]  
    return _Props(Parameter, _param1, value1, _param2, value2, Fluid)

@cython.final
cdef class State: 
    """
    A class that contains all the code that represents a thermodynamic state
    """
    
    def __init__(self,bytes Fluid, dict StateDict, double xL=-1.0, bytes Liquid=str(''), phase = str('Gas')):
        self.Fluid = Fluid
        self.iFluid = _get_Fluid_index(Fluid)
        #Try to get the fluid from CoolProp
        if self.iFluid >= 0:
            #It is a CoolProp Fluid so we can use the faster integer passing function
            self.is_CPFluid = True
        else:
            self.is_CPFluid = False
        self.xL = xL
        self.Liquid = Liquid
        self.phase = phase
        #Parse the inputs provided
        self.update(StateDict)
        #Set the phase flag
        if self.phase == str('Gas') or self.phase == str('Liquid') or self.phase == str('Supercritical'):
            _set_phase(self.phase)
            
    def __reduce__(self):
        d={}
        d['xL']=self.xL
        d['Liquid']=self.Liquid
        d['Fluid']=self.Fluid
        d['T']=self.T_
        d['rho']=self.rho_
        return rebuildState,(d,)
          
          
    cpdef update_Trho(self, double T, double rho):
        """
        Just use the temperature and density directly
        """
        self.T_ = T
        self.rho_ = rho
        cdef double p
        
        if self.is_CPFluid:
            p = _IProps(iP,iT,T,iD,rho,self.iFluid)
        else:
            p = _Props('P','T',T,'D',rho,self.Fluid)
        
        if abs(p)<1e90:
            self.p_=p
        else:
            errstr = get_errstringc()
            raise ValueError(errstr)
        
    cpdef update(self,dict params, double xL=-1.0):
        """
        *params* is a list(or tuple) of strings that represent the parameters 
        that have been updated and will be used to fix the rest of the state. 
        ['T','P'] for temperature and pressure for instance
        """
            
        cdef double p
        cdef bytes errstr
        
        # If no value for xL is provided, it will have a value of -1 which is 
        # impossible, so don't update xL
        if xL > 0:
            #There is liquid
            self.xL=xL
            self.hasLiquid=True
        else:
            #There's no liquid
            self.xL=0.0
            self.hasLiquid=False
        
        #You passed in a dictionary, use the values to update the state
        if 'T' not in params:
            raise AttributeError('T must be provided in params dict in State.update')
            
        #Consume the 'T' key since it is required (TODO?)
        self.T_=float(params.pop('T'))
            
        #Given temperature and pressure, determine density of gas 
        # (or gas and oil if xL is provided)
        if abs(self.xL)<=1e-15:
            #Get the density if T,P provided, or pressure if T,rho provided
            if 'P' in params:
                self.p_=params['P']
                
                if self.is_CPFluid:
                    rho = _IProps(iD,iT,self.T_,iP,self.p_,self.iFluid)
                else:
                    rho = _Props('D','T',self.T_,'P',self.p_,self.Fluid)
                
                if abs(rho)<1e90:
                    self.rho_=rho
                else:
                    errstr = get_errstringc()
                    raise ValueError(errstr)
            elif 'D' in params:
                
                self.rho_=params['D']
                
                if self.is_CPFluid:
                    p = _IProps(iP,iT,self.T_,iD,self.rho_,self.iFluid)
                else:
                    p = _Props('P','T',self.T_,'D',self.rho_,self.Fluid)
                
                if abs(p)<1e90:
                    self.p_=p
                else:
                    errstr = get_errstringc()
                    raise ValueError(errstr+str(params))
            else:
                raise KeyError("Dictionary must contain the key 'T' and one of 'P' or 'D'")
            
        elif self.xL>0 and self.xL<=1:
            raise ValueError('xL is out of range - value for xL is [0,1]')
        else:
            raise ValueError('xL must be between 0 and 1')
        
    cpdef double Props(self, long iOutput):
        if iOutput<0:
            raise ValueError('Your output is invalid') 
        
        if self.is_CPFluid:
            return _IProps(iOutput,iT,self.T_,iD,self.rho_,self.iFluid)
        else:
            return _Props(paras[iOutput],'T',self.T_,'D',self.rho_,self.Fluid)
            
    cpdef double get_MM(self):
        return _Props1(self.Fluid,'molemass')
    
    cpdef double get_rho(self): 
        return self.rho_
    property rho:
        def __get__(self):
            return self.rho_
            
    cpdef double get_p(self): 
        return self.p_
    property p:
        def __get__(self):
            return self.p_
    
    cpdef double get_T(self): 
        return self.T_
    property T:
        def __get__(self):
            return self.T_
    
    cpdef double get_h(self): 
        return self.Props(iH)
    property h:
        def __get__(self):
            return self.get_h()
          
    cpdef double get_u(self): 
        return self.Props(iU)
    property u:
        def __get__(self):
            return self.get_u()
            
    cpdef double get_s(self): 
        return self.Props(iS)
    property s:
        def __get__(self):
            return self.get_s()
    
    cpdef double get_cp0(self):
        return self.Props(iC0)
    
    cpdef double get_cp(self): 
        return self.Props(iC)
    property cp:
        def __get__(self):
            return self.get_cp()
            
    cpdef double get_cv(self): 
        return self.Props(iO)
    property cv:
        def __get__(self):
            return self.get_cv()
            
    cpdef double get_visc(self):
        return self.Props(iV)
    property visc:
        def __get__(self):
            return self.get_visc()

    cpdef double get_cond(self):
        return self.Props(iL)
    property k:
        def __get__(self):
            return self.get_cond()
            
    property Prandtl:
        def __get__(self):
            return self.cp * self.visc / self.k
            
    cpdef double get_dpdT(self):
        return self.Props(iDpdT)
    property dpdT:
        def __get__(self):
            return self.get_dpdT()
        
    cpdef speed_test(self, int N):
        from time import clock
        cdef int i
        cdef char * k
        cdef char * Fluid = self.Fluid
        cdef long IT = 'T'
        cdef long ID = 'D'
        
        print 'Direct c++ call to CoolProp without the Python call layer'
        print "'M' involves basically no computational effort and is a good measure of the function call overhead"
        keys = ['H','P','S','U','C','O','V','L','M','C0']
        for key in keys:
            t1=clock()
            for i in range(N):
                _Props(key,'T',self.T_,'D',self.rho_,Fluid)
            t2=clock()
            print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,key,(t2-t1)/N*1e6)
        
        print 'Call to the c++ layer using integers'
        keys = [iH,iP,iS,iU,iC,iO,iV,iL,iMM,iC0]
        for key in keys:
            t1=clock()
            for i in range(N):
                _IProps(key,iT,self.T_,iD,self.rho_,self.iFluid)
            t2=clock()
            print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,paras[key],(t2-t1)/N*1e6)
            
        print 'Call to the Python call layer'
        print "'M' involves basically no computational effort and is a good measure of the function call overhead"
        keys = ['H','P','S','U','C','O','V','L','M','C0']
        for key in keys:
            t1=clock()
            for i in range(N):
                CP.Props(key,'T',self.T_,'D',self.rho_,Fluid)
            t2=clock()
            print 'Elapsed time for {0:d} calls for "{1:s}" at {2:g} us/call'.format(N,key,(t2-t1)/N*1e6)
    
    def __str__(self):
        """
        Return a string representation of the state
        """
        units={'T': 'K', 
               'p': 'kPa', 
               'rho': 'kg/m^3',
               'h':'kJ/kg',
               'u':'kJ/kg',
               's':'kJ/kg/K',
               'visc':'Pa-s',
               'k':'kW/m/K',
               'cp':'kJ/kg/K',
               'cv':'kJ/kg/K',
               'dpdT':'kPa/K'}
        s=''
        for k in ['T','p','rho','h','u','s','visc','k','cp','cv','dpdT','Prandtl']:
            if k in units:
                s+=k+' = '+str(getattr(self,k))+' '+units[k]+'\n'
            else:
                s+=k+' = '+str(getattr(self,k))+' NO UNITS'+'\n'
        return s.rstrip()
        
    cpdef copy(self):
        cdef double T = self.T_*(1.0+1e-20)
        cdef double rho = self.rho_*(1.0+1e-20)
        ST=State(self.Fluid,{'T':T,'D':rho})
        return ST
    
def rebuildState(d):
    S=State(d['Fluid'],{'T':d['T'],'D':d['rho']})
    S.xL = d['xL']
    S.Liquid=d['Liquid']
    return S
