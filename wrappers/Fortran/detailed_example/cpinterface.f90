MODULE CPINTERFACE
    INTERFACE
        FUNCTION PropsSI (output, name1, prop1, name2, prop2, fluidname) BIND(C, NAME='PropsSI')
            use iso_c_binding
            real(C_DOUBLE) :: PropsSI
            character(KIND=c_char) :: output(*)
            character(c_char) :: name1(*)
            real(C_DOUBLE), VALUE :: prop1
            character(c_char) :: name2(*)
            real(C_DOUBLE), VALUE :: prop2
            character(kind=c_char) :: fluidname(*)
                
        END FUNCTION PropsSI
    END INTERFACE
END MODULE CPINTERFACE