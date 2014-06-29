subroutine myfunc(a, b, c) bind(C, NAME="myfunc")!
!
use iso_c_binding!
!
integer(c_int), value :: a!
real(c_double), value :: b!
real(c_double), dimension(*) :: c!
!
integer :: i!
!
print *, a!
print *, b!
!
do i=1,20!
print*, c(i)!
end do!
!
end subroutine myfunc