module const_m
  implicit none
  integer, parameter :: SI = selected_int_kind(8)
  integer, parameter :: DI = selected_int_kind(16)
  integer, parameter :: SR = selected_real_kind(5)
  integer, parameter :: DR = selected_real_kind(15)
  real(DR), parameter :: PI = atan(1.0_DR)*4
  real(DR), parameter :: TWOPI = 2*PI

  integer, parameter :: FILE_STANDARD_OUT = 6

contains

  subroutine const__print
    print *, ' SI = ', SI
    print *, ' DI = ', DI
    print *, ' SR = ', SR
    print *, ' DR = ', DR
    print *, ' PI = ', PI
    print *, ' TWOPI = ', TWOPI
  end subroutine const__print

end module const_m
