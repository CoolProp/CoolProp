from libcpp cimport bool 
from libcpp.string cimport string

cdef class PureFluidClass:
    cdef CoolPropStateClass CPS     # hold a C++ instance which we're wrapping
    cpdef update(self, long iInput1, double Value1, long iInput2, double Value2)
    cpdef double rhoL(self)
    cpdef double rhoV(self)
    cpdef double pL(self)
    cpdef double pV(self)
    cpdef double TL(self)
    cpdef double TV(self)
    cpdef double sL(self)
    cpdef double sV(self)
    cpdef double hL(self)
    cpdef double hV(self)
    
    ## ---------------------------------------- 
    ##        Fluid property accessors
    ## ----------------------------------------
    
    cpdef double T(self)
    cpdef double rho(self)
    cpdef double p(self)
    cpdef double h(self)
    cpdef double s(self)
    cpdef double cp(self)
    cpdef double cv(self)
    cpdef double speed_sound(self)

    ## ---------------------------------------- 
    ##        TTSE LUT things
    ## ----------------------------------------

    
    cpdef enable_TTSE_LUT(self) # Enable the TTSE
    cpdef bool isenabled_TTSE_LUT(self) # Check if TTSE is enabled
    cpdef disable_TTSE_LUT(self) # Disable the TTSE

    cpdef double dTdp_along_sat(self)
    cpdef double d2Tdp2_along_sat(self)

    cpdef double dhdp_along_sat_vapor(self)
    cpdef double dhdp_along_sat_liquid(self)
    cpdef double d2hdp2_along_sat_vapor(self)
    cpdef double d2hdp2_along_sat_liquid(self)

    cpdef double dsdp_along_sat_vapor(self)
    cpdef double dsdp_along_sat_liquid(self)
    cpdef double d2sdp2_along_sat_vapor(self)
    cpdef double d2sdp2_along_sat_liquid(self)

    cpdef double drhodp_along_sat_vapor(self)
    cpdef double drhodp_along_sat_liquid(self)
    cpdef double d2rhodp2_along_sat_vapor(self)
    cpdef double d2rhodp2_along_sat_liquid(self)

    cpdef double drhodT_along_sat_vapor(self)
    cpdef double drhodT_along_sat_liquid(self)
    
    cpdef double drhodT_constp(self)
    cpdef double drhodp_constT(self)
    cpdef double d2rhodp2_constT(self)
    cpdef double d2rhodTdp(self)
    cpdef double d2rhodT2_constp(self)
    cpdef double d2rhodhdQ(self)
    cpdef double d2rhodpdQ(self)
    cpdef double d2rhodhdp(self)
    cpdef double d2rhodh2_constp(self)

    cpdef double dpdrho_constT(self)
    cpdef double dpdrho_consth(self)
    cpdef double dpdT_constrho(self)
    cpdef double dpdT_consth(self)
    cpdef double d2pdrho2_constT(self)
    cpdef double d2pdrhodT(self)
    cpdef double d2pdT2_constrho(self)

    cpdef double dhdrho_constT(self)
    cpdef double dhdrho_constp(self)
    cpdef double dhdT_constrho(self)
    cpdef double dhdT_constp(self)
    cpdef double dhdp_constT(self)
    cpdef double d2hdrho2_constT(self)
    cpdef double d2hdrhodT(self)
    cpdef double d2hdT2_constrho(self)
    cpdef double d2hdT2_constp(self)
    cpdef double d2hdp2_constT(self)
    cpdef double d2hdTdp(self)

    cpdef double dsdrho_constT(self)
    cpdef double dsdT_constrho(self)
    cpdef double dsdrho_constp(self)
    cpdef double dsdT_constp(self)
    cpdef double dsdp_constT(self)
    cpdef double d2sdrho2_constT(self)
    cpdef double d2sdrhodT(self)
    cpdef double d2sdT2_constrho(self)
    cpdef double d2sdT2_constp(self)
    cpdef double d2sdp2_constT(self)
    cpdef double d2sdTdp(self)