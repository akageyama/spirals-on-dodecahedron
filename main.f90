!
! multiple spiral pattern on regular dodecahedron
!
!   developed by Akira Kageyama,
!                Kobe University,
!             on Feb. 2021,
!            for multiple spiral convection.
!
program main
  use const_m
  use vector_m
  use spiral_m
  implicit none

  integer(SI), parameter :: NTHT = 401
  integer(SI), parameter :: NPHI = 800

  real(DR), parameter :: DTHT =    PI / (NTHT-1)
  real(DR), parameter :: DPHI = TWOPI / NPHI

  type(spiral_t), dimension(6) :: spirals
  type(vector_t), dimension(12) :: work_centers

  integer, parameter :: NUM_ARM_SPIRAL = 4

  real(DR) :: tht, phi, x, y, z
  type(vector_t) :: point
  integer(SI) :: j, k
  real(DR), dimension(6) :: amps ! amplitudes
  real(DR), dimension(6) :: closest_dist_to_arm_in_spirals
  integer :: n
  real(DR), parameter :: EPSILON_FACTOR = 0.28_DR
                        ! A control parameter to enhance the
                        ! separation between neigbhoring spirals. 
                        ! Have found by trials and erros.
  real(DR) :: epsilon = EPSILON_FACTOR / NUM_ARM_SPIRAL

  real(DR), parameter :: GR  = (1+sqrt(5.0_DR))/2  ! golden ratio

  call const__print

  ! central points of the 12 surfaces of dodecahedron
  call work_centers( 1)%set(    +1.0_DR,        +GR,  0+epsilon )
  call work_centers( 2)%set(    -1.0_DR,        +GR,  0-epsilon )
  call work_centers( 3)%set(    +1.0_DR,        -GR,  0-epsilon )
  call work_centers( 4)%set(    -1.0_DR,        -GR,  0+epsilon )
  call work_centers( 5)%set(  0+epsilon,    +1.0_DR,        +GR )
  call work_centers( 6)%set(  0-epsilon,    -1.0_DR,        +GR )
  call work_centers( 7)%set(  0-epsilon,    +1.0_DR,        -GR )
  call work_centers( 8)%set(  0+epsilon,    -1.0_DR,        -GR )
  call work_centers( 9)%set(        +GR,  0+epsilon,    +1.0_DR )
  call work_centers(10)%set(        +GR,  0-epsilon,    -1.0_DR )
  call work_centers(11)%set(        -GR,  0-epsilon,    +1.0_DR )
  call work_centers(12)%set(        -GR,  0+epsilon,    -1.0_DR )

  do n = 1, 12 ! 12 surfaces in dodecahedron
    call work_centers(n)%normalize
  end do

  do n = 1, 6  ! 6 pairs (spirals) 
    call spirals(n)%initialize( NUM_ARM_SPIRAL,  &
                                work_centers(2*n-1), &
                                work_centers(2*n  ) )
  end do

  open( 10, file='output.data' )
  do k = 1, NPHI
    phi = DPHI*(k-1)
    do j = 1, NTHT
      tht = DTHT*(j-1)
      x = sin(tht)*cos(phi)
      y = sin(tht)*sin(phi)
      z = cos(tht)
      call point%set( x, y, z )
      do n = 1, 6
        closest_dist_to_arm_in_spirals(n)  &
                = spirals(n)%dist_to_closest_arm( NUM_ARM_SPIRAL, &
                                                  point )
        amps(n) = delta( spirals(n)%width,  &
                         closest_dist_to_arm_in_spirals(n) )
      end do
      write(10,*) phi, tht, sum( amps(:) )
    end do
    write(10,*) ' '
  end do
  close(10)
       
contains


  function delta( width, x )
    real(DR), intent(in) :: width
    real(DR), intent(in) :: x
    real(DR) :: delta
    
    real(DR) :: x2

    !             |            
    !         1  *|*           
    !             |            
    !           * | *          
    !             |            
    !          *  |  *         
    !             |            
    !   -*****----|----****---------> x
    !             0             

    x2 = (x/width)**2

    if ( x2 > 100.0_DR ) then
      delta = 0.0_DR
    else
      delta = exp( -x2/2 )
    end if
  end function delta

end program main
