#This file gets directly included in CoolProp.pyx, separate here for cleanness of code

cpdef HAProps(str OutputName, str Input1Name, Input1, str Input2Name, Input2, str Input3Name, Input3):
    """
    Copyright Ian Bell, 2011 email: ian.h.bell@gmail.com

    The function is called like

    HAProps('H','T',298.15,'P',101.325,'R',0.5)

    which will return the enthalpy of the air for a set of inputs of dry bulb temperature of 25C, atmospheric pressure, and a relative humidity of 50%.

    This function implements humid air properties based on the analysis in ASHRAE RP-1845 which is available online: http://rp.ashrae.biz/page/ASHRAE-D-RP-1485-20091216.pdf

    It employs real gas properties for both air and water, as well as the most accurate interaction parameters and enhancement factors.  The IAPWS-95 formulation for the properties of water is used throughout in preference to the industrial formulation.  It is unclear why the industrial formulation is used in the first place.

    Since humid air is nominally a binary mixture, three variables are needed to fix the state.  At least one of the input parameters must be dry-bulb temperature, relative humidity, dew-point temperature, or humidity ratio.  The others will be calculated.  If the output variable is a transport property (conductivity or viscosity), the state must be able to be found directly - i.e. make sure you give temperature and relative humidity or humidity ratio.  The list of possible input variables are

    ========  ========    ========================================
    String    Aliases     Description
    ========  ========    ========================================
    T         Tdb         Dry-Bulb Temperature [K]
    B         Twb         Wet-Bulb Temperature [K]
    D         Tdp         Dew-Point Temperature [K]
    P                     Pressure [kPa]
    V         Vda         Mixture volume [m3/kg dry air]
    R         RH          Relative humidity in (0,1) [-]
    W         Omega       Humidity Ratio [kg water/kg dry air]
    H         Hda         Mixture enthalpy [kJ/kg dry air]
    S         Sda         Mixture entropy [kJ/kg dry air/K]
    C         cp          Mixture specific heat [kJ/kg dry air/K]
    M         Visc        Mixture viscosity [Pa-s]
    K                     Mixture thermal conductivity [kW/m/K]
    ========  ========    ========================================

    There are also strings for the mixture volume and mixture enthalpy that will return the properties on a total humid air flow rate basis, they are given by 'Vha' [units of m^3/kg humid air] and 'Cha' [units of kJ/kg humid air/K] and 'Hha' [units of kJ/kg humid air] respectively.

    For more information, go to http://coolprop.sourceforge.net
    """
    #Convert all strings to byte-strings
    cdef bytes _OutputName = OutputName.encode('ascii')
    cdef bytes _Input1Name = Input1Name.encode('ascii')
    cdef bytes _Input2Name = Input2Name.encode('ascii')
    cdef bytes _Input3Name = Input3Name.encode('ascii')
    
    if isinstance(Input1, (int, long, float, complex)) and isinstance(Input2, (int, long, float, complex)) and isinstance(Input3, (int, long, float, complex)):
        val = _HAProps(_OutputName,_Input1Name,Input1,_Input2Name,Input2,_Input3Name,Input3)
    
        if math.isinf(val) or math.isnan(val):
            err_string = _get_global_param_string('errstring')
            if not len(err_string) == 0:
                raise ValueError("{err:s} :: inputs were:\"{out:s}\",\'{in1n:s}\',{in1:0.16e},\'{in2n:s}\',{in2:0.16e},\'{in3n:s}\',{in3:0.16e} ".format(err=err_string,out=_OutputName,in1n=_Input1Name,in1=Input1,in2n=_Input2Name,in2=Input2,in3n=_Input3Name,in3=Input3))
            else:
                raise ValueError("HAProps failed ungracefully with inputs: \"{out:s}\",\'{in1n:s}\',{in1:0.16e},\'{in2n:s}\',{in2:0.16e},\'{in3n:s}\',{in3:0.16e} ".format(out=_OutputName,in1n=_Input1Name,in1=Input1,in2n=_Input2Name,in2=Input2,in3n=_Input3Name,in3=Input3))
        
        return val #Error raised by HAProps on failure
        
    # At least one is iterable, convert non-iterable to a list of the same length
    elif isinstance(Input1, (int, long, float, complex)) or isinstance(Input2, (int, long, float, complex)):
        
        iterable_lengths = []
        if not isinstance(Input1, (int, long, float, complex)):
            iterable_lengths.append(len(Input1))
        if not isinstance(Input2, (int, long, float, complex)):
            iterable_lengths.append(len(Input2))
        if not isinstance(Input3, (int, long, float, complex)):
            iterable_lengths.append(len(Input3))
        
        if not len(set(iterable_lengths)) == 1:
            raise TypeError("Iterable inputs are not all the same length.  Lengths: "+str(iterable_lengths))
        else:
            L = iterable_lengths[0]
            
        
        if isinstance(Input1, (int, long, float, complex)):
            Input1vec = [Input1]*L
        else:
            Input1vec = Input1
            
        if isinstance(Input2, (int, long, float, complex)):
            Input2vec = [Input2]*L
        else:
            Input2vec = Input2
            
        if isinstance(Input3, (int, long, float, complex)):
            Input3vec = [Input3]*L
        else:
            Input3vec = Input3
        
        vals = []
        
        for _Input1, _Input2, _Input3 in zip(Input1vec, Input2vec, Input3vec):
            val = _HAProps(_OutputName,_Input1Name,_Input1,_Input2Name,_Input2,_Input3Name,_Input3)
        
            if math.isinf(val) or math.isnan(val):
                err_string = _get_global_param_string('errstring')
                if not len(err_string) == 0:
                    raise ValueError("{err:s} :: inputs were:\"{out:s}\",\'{in1n:s}\',{in1:0.16e},\'{in2n:s}\',{in2:0.16e},\'{in3n:s}\',{in3:0.16e}".format(err=err_string,out=_OutputName,in1n=_Input1Name,in1=_Input1,in2n=_Input2Name,in2=_Input2,in3n=_Input3Name,in3=_Input3))
                else:
                    raise ValueError("HAProps failed ungracefully with inputs: \"{out:s}\",\'{in1n:s}\',{in1:0.16e},\'{in2n:s}\',{in2:0.16e},\'{in3n:s}\',{in3:0.16e} ".format(out=_OutputName,in1n=_Input1Name,in1=Input1,in2n=_Input2Name,in2=Input2,in3n=_Input3Name,in3=Input3))
            
            vals.append(val)
            
        if _numpy_supported and isinstance(Input1, np.ndarray):
            return np.array(vals).reshape(Input1.shape)
        elif _numpy_supported and isinstance(Input2, np.ndarray):
            return np.array(vals).reshape(Input2.shape)
        elif _numpy_supported and isinstance(Input3, np.ndarray):
            return np.array(vals).reshape(Input3.shape)
        else:
            return vals
    
    else:
        raise TypeError('Numerical inputs to Props must be ints, floats, lists, or 1D numpy arrays.')

cpdef tuple HAProps_Aux(str OutputName, double T, double p, double w):
    """
    Allows low-level access to some of the routines employed in HumidAirProps

    Returns tuples of the form ``(Value, Units)`` where ``Value`` is the actual value and ``Units`` is a string that describes the units

    The list of possible inputs is

    * Baa [First virial air-air coefficient]
    * Caaa [Second virial air coefficient]
    * Bww [First virial water-water coefficient]
    * Cwww [Second virial water coefficient]
    * Baw [First cross virial coefficient]
    * Caww [Second air-water-water virial coefficient]
    * Caaw [Second air-air-water virial coefficient]
    * beta_H 
    * kT
    * vbar_ws [Molar saturated volume of water vapor]
    * p_ws [Saturated vapor pressure of pure water (>=0.01C) or ice (<0.01 C)]
    * f [Enhancement factor]
    """
    #Convert all strings to byte-strings
    cdef bytes units = (' '*100).encode('ascii')
    cdef bytes _OutputName = OutputName.encode('ascii')
    
    output = _HAProps_Aux(_OutputName,T,p,w,units)
    units = units.strip()
    units = units[0:len(units)-1] #Leave off the null character
    return output, units
    
cpdef double cair_sat(double T):
    """
    The derivative of the saturation enthalpy cair_sat = d(hsat)/dT
    """
    return _cair_sat(T)