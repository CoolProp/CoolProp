
program hello

    USE cpinterface
    use iso_c_binding
    implicit none
      
    !Initialize the variables used in the example
    real(C_DOUBLE) T, Q
    character(LEN=32) fluid, out1, n1, n2
    real(C_DOUBLE) dens1, dens2, dens3

    CHARACTER(100) mesgbfr
    integer(C_long) handle
    integer(C_long) errcode, bfrlen
    integer(C_long) QT_INPUTS, iDmass

    mesgbfr(1:100) = ' '
    bfrlen = 100

    ! ----------------------
    ! Example calculates density of saturated liquid propane at 300 K:
    ! ---------------------
      
    T = 300                 ! Temperature [K]
    Q = 0                   ! Quality [-]

    out1 = "D"//C_NULL_CHAR      ! String with of the output Property
    n1 = "T"//C_NULL_CHAR      ! String with of the input Property #1
    n2 = "Q"//C_NULL_CHAR       ! String with of the input Property #2
    fluid = "Propane"//C_NULL_CHAR   ! String with the fluid name

    dens1 = PropsSI(out1, n1, T, n2, Q, fluid)                                !calling props, strings are Variables
    dens2 = PropsSI(C_CHAR_"D"//C_NULL_CHAR, C_CHAR_"T"//C_NULL_CHAR, T, C_CHAR_"Q"//C_NULL_CHAR, Q, fluid)    !calling props defining the strings directly in the arguments 
    Print *, dens1, dens2

    ! --------
    ! Do the same calculation with the low-level interface
    
    QT_INPUTS = get_input_pair_index(C_CHAR_"QT_INPUTS"//C_NULL_CHAR)
    iDmass = get_param_index("Dmass"//C_NULL_CHAR)

    handle = AbstractState_factory(C_CHAR_"HEOS"//C_NULL_CHAR,fluid,errcode,mesgbfr,bfrlen)
    if (errcode .ne. 0) then
        Print *, handle, errcode, mesgbfr
        call exit()
    endif

    call AbstractState_update(handle,QT_INPUTS, 0.D0, 300.D0, errcode,mesgbfr,bfrlen)
    if (errcode .ne. 0) then
        Print *, errcode, mesgbfr
        call exit()
    endif
    
    dens3 = AbstractState_keyed_output(handle,iDmass,errcode,mesgbfr,bfrlen)
    if (errcode .ne. 0) then
        Print *, errcode, mesgbfr
        call exit()
    endif
    Print *, dens3

end program hello

