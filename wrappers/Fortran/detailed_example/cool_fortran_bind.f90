
program hello

    USE cpinterface
    use iso_c_binding
    implicit none
      
    !Initialize the variables used in the example
    real(C_DOUBLE) T, Q
    character(LEN=32) fluid, out1, n1, n2
    real(C_DOUBLE) dens1, dens2, c_p_h2o

    CHARACTER(100) mesgbfr
    integer(C_long) handle
    integer(C_long) errcode, bfrlen
    integer(C_long) PT_INPUTS, icp

    mesgbfr(1:100) = ' '
    bfrlen = 100

    !----------------------
    !Example calculates density of saturated liquid propane at 300 K:
    !---------------------
      
    T = 300                 ! Temperature [K]
    Q = 0                   ! Quality [-]

    out1 = "D"//CHAR(0)      ! String with of the output Property
    n1 = "T"//CHAR(0)       ! String with of the input Property #1
    n2 = "Q"//CHAR(0)       ! String with of the input Property #2
    fluid = "Propane"//CHAR(0)   ! String with the fluid name

    dens1 = PropsSI(out1, n1, T, n2, Q, fluid)                                !calling props, strings are Variables
    dens2 = PropsSI("D"//CHAR(0), "T"//CHAR(0), T, "Q"//CHAR(0), Q, fluid)    !calling props defining the strings directly in the arguments 
    Print *, dens1, dens2

    !try with the low level interface
    
    PT_INPUTS = get_input_pair_index("PT_INPUTS"//CHAR(0))
    Print *, PT_INPUTS
    icp = get_param_index("C"//CHAR(0))
    Print *, icp

    handle = AbstractState_factory("HEOS"//CHAR(0),"Water"//CHAR(0),errcode,mesgbfr,bfrlen)
    if (errcode .ne. 0) then
        Print *, handle, errcode, mesgbfr
        call exit()
    endif

    Print *, "handle", handle
    call AbstractState_update(handle,PT_INPUTS,101325.D0, 300.D0,errcode,mesgbfr,bfrlen)
    if (errcode .ne. 0) then
        Print *, errcode, mesgbfr
        call exit()
    endif
    
    c_p_h2o = AbstractState_keyed_output(handle,icp,errcode,mesgbfr,bfrlen)
    Print *, c_p_h2o, errcode, bfrlen, mesgbfr

end program hello

