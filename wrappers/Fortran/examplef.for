! Example fortran file 
! 
      double precision T,Q,D,h,s
      character(LEN=32) Ref,Output, Name1, Name2
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
      outval = 9999999

!       Output(LEN(Output):LEN(Output)) = CHAR(0)
!       Name1(LEN(Name1):LEN(Name1)) = CHAR(0)
!       Name2(LEN(Name2):LEN(Name2)) = CHAR(0)
!       Ref(LEN(Ref):LEN(Ref)) = CHAR(0)

      write(*,*) "Saturation pressure for R134a: "
      call propssi(Output, Name1, Prop1, Name2, Prop2, Ref, outVal)
      write(*,*) "Result was: ", outVal/1e2, " bar"
      write(*,*) "-----------------------------------------------"
c
c      Output = "S"//CHAR(0)
c
c
c
c      outVal = props1(Ref , "Tcrit"//CHAR(0))
c      write(*,*) "Tcrit from props1     : ", outVal 
c      outVal = props(Output, Name1, Prop1, Name2, Prop2, Ref)
c      write(*,*) "S from standard props : ", outVal
c      outVal = 0.0
c      CALL propsmod(Output, Name1, Prop1, Name2, Prop2, Ref, outVal)
c      write(*,*) "S from modified props : ", outVal
c      outVal = derivterms("dpdrho"//CHAR(0), Prop1, D, Ref)
c      write(*,*) "dpdrho from derivterms: ", outVal
c      inte = setreferencestates(Ref, "IIR"//CHAR(0))
c      write(*,*) "reference to IIR      : ", inte
c      outVal = props(Output, Name1, Prop1, Name2, Prop2, Ref)
c      write(*,*) "S from standard props : ", outVal
c      inte = setreferencestates(Ref, "ASHRAE"//CHAR(0))
c      write(*,*) "reference to ASHRAE   : ", inte
c      outVal = props(Output, Name1, Prop1, Name2,Prop2, Ref)
c      write(*,*) "S from standard props : ", outVal
c      inte = enablettselut(Ref)
c      write(*,*) "enabling TTSE         : ", inte
c      inte = isenabledttselut(Ref)
c      write(*,*) "state of TTSE         : ", inte
c      outVal = props(Output, Name1, Prop1, Name2, Prop2, Ref)
c      write(*,*) "S from TTSE props     : ", outVal
c      inte = setttsemode(Ref , "bicubic"//CHAR(0))
c      write(*,*) "TTSE mode to bicubic  : ", inte
c      outVal = props(Output, Name1, Prop1, Name2, Prop2, Ref)
c      write(*,*) "S from bicubic props  : ", outVal
c      inte = disablettselut(Ref)
c      write(*,*) "Disabling TTSE        : ", inte
c      outVal = props(Output, Name1, Prop1, Name2, Prop2, Ref)
c      write(*,*) "S from standard props : ", outVal

      end program 
