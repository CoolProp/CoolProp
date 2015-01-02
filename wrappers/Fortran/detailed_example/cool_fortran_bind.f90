
program hello

    USE cpinterface
    
    implicit none

    !Initialize the variables used in the example
    double precision T, Q
    character(LEN=32) fluid, out1, n1, n2
    double precision dens1, dens2

    !----------------------
    !Example calculates density of saturated liquid propane at 300 K:
    !---------------------
      
    T = 300                 ! Temperature [K]
    Q = 0                   ! Quality [-]

    out1 = "D"//CHAR(0)      ! String with of the output Property
    n1  = "T"//CHAR(0)       ! String with of the input Property #1
    n2  = "Q"//CHAR(0)       ! String with of the input Property #2
    fluid    = "Propane"//CHAR(0)   ! String with the fluid name
      
    dens1 = PropsSI(out1, n1, T, n2, Q, fluid)                                !calling props, strings are Variables
      
    dens2 = PropsSI("D"//CHAR(0), "T"//CHAR(0), T, "Q"//CHAR(0), Q, fluid)    !calling props defining the strings directly in the arguments 

    Print *, dens1, dens2  

end program hello