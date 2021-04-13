module vector_m
  use const_m
  implicit none
  private

  public :: operator(+)
  public :: operator(-)
  public :: operator(.dot.)
  public :: operator(.cross.)
  public :: assignment(=)

  interface operator(+)
    module procedure operator_add
  end interface

  interface operator(-)
    module procedure operator_subtract
  end interface

  interface operator(.dot.)
    module procedure operator_dot_product
  end interface

  interface operator(.cross.)
    module procedure operator_outer_product
  end interface

  interface assignment(=)
    module procedure operator_assign
  end interface

  type, public :: vector_t
    real(DR) :: x, y, z
  contains
    procedure :: print => vector_print
    procedure :: amplitude => vector_amplitude
    procedure :: normalize => vector_normalize
    procedure :: vector_set_components
    procedure :: vector_set_a_real
    generic :: set => vector_set_components, &
                      vector_set_a_real
    procedure :: amplify => vector_amplify
  end type vector_t


contains


  function operator_add( v1, v2 )
    type(vector_t), intent(in) :: v1, v2
    type(vector_t) :: operator_add
    operator_add%x = v1%x + v2%x
    operator_add%y = v1%y + v2%y
    operator_add%z = v1%z + v2%z
  end function operator_add


  function operator_subtract( v1, v2 )
    type(vector_t), intent(in) :: v1, v2
    type(vector_t) :: operator_subtract
    operator_subtract%x = v1%x - v2%x
    operator_subtract%y = v1%y - v2%y
    operator_subtract%z = v1%z - v2%z
  end function operator_subtract


  function operator_dot_product( a, b )
    type(vector_t), intent(in) :: a, b
    real(DR) :: operator_dot_product

    operator_dot_product = a%x*b%x + a%y*b%y + a%z*b%z
  end function operator_dot_product


  function operator_outer_product( a, b )
    type(vector_t), intent(in) :: a, b
    type(vector_t) :: operator_outer_product

    operator_outer_product%x = a%y * b%z - a%z * b%y
    operator_outer_product%y = a%z * b%x - a%x * b%z
    operator_outer_product%z = a%x * b%y - a%y * b%x
  end function operator_outer_product


  subroutine operator_assign( v1, v2 )
    type(vector_t), intent(out) :: v1
    type(vector_t), intent(in) :: v2

    v1%x = v2%x
    v1%y = v2%y
    v1%z = v2%z
  end subroutine operator_assign


  function vector_amplitude( self )
    class(vector_t), intent(in) :: self
    real(DR) :: vector_amplitude

    vector_amplitude = sqrt ( self%x**2 &
                            + self%y**2 &
                            + self%z**2 )
  end function vector_amplitude


  subroutine vector_amplify( self, factor )
    class(vector_t), intent(inout) :: self
    real(DR), intent(in) :: factor

    self%x = self%x * factor
    self%y = self%y * factor
    self%z = self%z * factor
  end subroutine vector_amplify


  subroutine vector_normalize( self )
    class(vector_t), intent(inout) :: self
   
    call self%amplify( 1.0_DR / self%amplitude() )
  end subroutine vector_normalize


  subroutine vector_set_components( self, x, y, z )
    class(vector_t), intent(out) :: self
    real(DR), intent(in) :: x, y, z

    self%x = x
    self%y = y
    self%z = z
  end subroutine vector_set_components


  subroutine vector_print( self )
    class(vector_t), intent(in) :: self
    print *, ' <vector_t> x, y, z = ', self%x, &
                                       self%y, &
                                       self%z
  end subroutine vector_print


  subroutine vector_set_a_real( self, v )
    class(vector_t), intent(out) :: self
    real(DR), intent(in) :: v

    self%x = v
    self%y = v
    self%z = v
  end subroutine vector_set_a_real

!
! Public
!

end module vector_m
