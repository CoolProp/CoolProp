! Example fortran file 
! 
      double precision T,Q,D,h,s
      character(LEN=32) Ref
      character(LEN=2) Output, Name1, Name2
      double precision outVal,Prop1,Prop2

      integer inte

      T = 285
      Q = 0
      D = 1250;

      Output = "P"//CHAR(0)
      Name1  = "T"//CHAR(0)
      Prop1  = T
      Name2  = "Q"//CHAR(0)
      Prop2  = Q
      Ref    = "R134a"//CHAR(0)

!       Output(LEN(Output):LEN(Output)) = CHAR(0)
!       Name1(LEN(Name1):LEN(Name1)) = CHAR(0)
!       Name2(LEN(Name2):LEN(Name2)) = CHAR(0)
!       Ref(LEN(Ref):LEN(Ref)) = CHAR(0)

      write(*,*) "Saturation pressure for R134a: "
      outVal = props(Output, Name1, Prop1, Name2, Prop2, Ref)
      write(*,*) "Result was: ", outVal/1e2, " bar"
      write(*,*) "-----------------------------------------------"

      Output = "S"//CHAR(0)



      outVal = props1(Ref , "Tcrit"//CHAR(0))
      write(*,*) "Tcrit from props1     : ", outVal 
      outVal = props(Output, Name1, Prop1, Name2, Prop2, Ref)
      write(*,*) "S from standard props : ", outVal
      outVal = 0.0
      CALL propsmod(Output, Name1, Prop1, Name2, Prop2, Ref, outVal)
      write(*,*) "S from modified props : ", outVal
      outVal = derivterms("dpdrho"//CHAR(0), Prop1, D, Ref)
      write(*,*) "dpdrho from derivterms: ", outVal
      inte = setreferencestates(Ref, "IIR"//CHAR(0))
      write(*,*) "reference to IIR      : ", inte
      outVal = props(Output, Name1, Prop1, Name2, Prop2, Ref)
      write(*,*) "S from standard props : ", outVal
      inte = setreferencestates(Ref, "ASHRAE"//CHAR(0))
      write(*,*) "reference to ASHRAE   : ", inte
      outVal = props(Output, Name1, Prop1, Name2,Prop2, Ref)
      write(*,*) "S from standard props : ", outVal
      inte = enablettselut(Ref)
      write(*,*) "enabling TTSE         : ", inte
      inte = isenabledttselut(Ref)
      write(*,*) "state of TTSE         : ", inte
      outVal = props(Output, Name1, Prop1, Name2, Prop2, Ref)
      write(*,*) "S from TTSE props     : ", outVal
      inte = setttsemode(Ref , "bicubic"//CHAR(0))
      write(*,*) "TTSE mode to bicubic  : ", inte
      outVal = props(Output, Name1, Prop1, Name2, Prop2, Ref)
      write(*,*) "S from bicubic props  : ", outVal
      inte = disablettselut(Ref)
      write(*,*) "Disabling TTSE        : ", inte
      outVal = props(Output, Name1, Prop1, Name2, Prop2, Ref)
      write(*,*) "S from standard props : ", outVal

      end program 
