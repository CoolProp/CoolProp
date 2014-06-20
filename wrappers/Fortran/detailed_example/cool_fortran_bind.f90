program hello

    !Use TestData       ! Include a Module with the interface, Optional (to not re-write the Interface every at Subroutine or Program)
    
    use iso_c_binding
    implicit none 

!The interface to use the PropsSI Function
INTERFACE
     FUNCTION PropsSI  &
        (output, name1, prop1, name2, prop2, fluidname)    &
            BIND(C, NAME='PropsSI')
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


!Here Start the programa


!Initializate the variables used in the example
   !implicit none   
double precision T,Q,D,h,s
character(LEN=32) fluid, PropName
character(LEN=32) out1, n1, n2
double precision dens1, dens2, prop1, prop2


!----------------------
!Example calculates: (Example Code for MATLAB)
!
!Density of saturated liquid Propane at 300 K:
!---------------------

  
T = 300                 ! Defining the temperature
Q = 0                   ! Defining the quality
!D = 1250;


out1 = "D"//CHAR(0)      ! String with of the output Propertie   
n1  = "T"//CHAR(0)       ! String with of the input Propertie #1
n2  = "Q"//CHAR(0)       ! String with of the input Propertie #2
fluid    = "Propane"//CHAR(0)   ! String with the fluid name

 
! From the CoolProp 4.2.1 documentation ->> CoolProp.CoolProp.Props
! Call Type #2:
!   Props(OutputName,InputName1,InputProp1,InputName2,InputProp2,Fluid) --> float
!
  
dens1 = PropsSI(out1, n1, T, n2, Q, fluid)                                !calling props, strings are Variables
  
dens2 = PropsSI("D"//CHAR(0), "T"//CHAR(0), T, "Q"//CHAR(0), Q, fluid)    !calling props deffinig the strings directlly in the arguments 


Print *, dens1
  

end program hello