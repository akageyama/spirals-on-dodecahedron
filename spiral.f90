module spiral_m
  use const_m
  use vector_m
  implicit none
  private
  public :: spiral__distance
  public :: spiral__angle
  public :: spiral__tau_from_distance
  public :: spiral__tau_from_angle

  type, public :: spiral_axis_t
    type(vector_t) :: position
    type(vector_t) :: base01, base02
  end type spiral_axis_t

  type, public :: spiral_t
    type(spiral_axis_t), dimension(2) :: axis
    type(vector_t) :: center
    real(DR) :: radius    ! of single spiral around an axis
    real(DR) :: width     ! of spiral arm
  contains
    procedure :: initialize => spiral_initialize
    procedure :: dist_to_closest_arm => spiral_dist_to_closest_arm
  end type spiral_t

  ! The spiral equation with parameter tau 
  !
  ! for 0<=tau<=1 
  !   distance_from_axis1(tau) = tau / TETRAHEDRON__DISTANCE_H
  !   angle_around_axis1(tau) = TWOPI * SPIRAL__MODE * tau
  !
  ! for 1<=tau<=2 
  !   distance_from_axis2(tau) = (2-tau) / TETRAHEDRON__DISTANCE_H
  !   angle_around_axis2(tau) = TWOPI * SPIRAL__MODE * (2-tau)
  !
  !            o o               o   o   o --------+
  !         o       o          o            o      |
  !       o    o o   o      o      o  o      o     |radius
  !           o   o   o    o     o      o          |
  !   ---o----o--*----o----o----o----*--o-----o----+--
  !            o |   o     o    o   o|  o     
  !        o     o o      o      o   o o     o 
  !           o  |       o          o|      o
  !              o  o  o             | o  o
  !              |                   |
  !            axis1               axis2
  !              
  !            base02
  !             /|\         base01<--+
  !              |                   |
  !              +-->base01         \|/
  !                                base02
  !
  !                      angle(tau)
  !       distance(tau) .             [tau=1]
  !          o   o   \ .         end of the fist half
  !      o            o                 |
  !    o      o o    .    o       AXES_DISTANCE_H
  !         o     o .       o           |
  !  o     o       o         o          |
  ! o-----o-------*----------o----------o-------> base01
  ! o     o     [tau=0]     o           o
  !        o    axis 01    o            
  !         o            o  \          o
  !   o       o        o     arm1       \   
  !              o  o               o    arm2
  !       o                      o
  !           o             o
  !               o    o


contains


  subroutine spiral_initialize( self, num_arm, axis_pos1, axis_pos2 )
    class(spiral_t), intent(out) :: self
    integer, intent(in) :: num_arm
    type(vector_t), intent(in) :: axis_pos1, axis_pos2

    type(vector_t) :: work_vec12, work_vec21
    real(DR) :: wave_length_along_radius
    integer :: n
    
    self%center = iCalc_average( axis_pos1, axis_pos2 )

    self%axis(1)%position = axis_pos1
    self%axis(2)%position = axis_pos2
    work_vec12 = axis_pos2 - axis_pos1
    work_vec21 = axis_pos1 - axis_pos2

    self%radius = acos( axis_pos1 .dot. axis_pos2 ) / 2

    call work_vec12%normalize
    call work_vec21%normalize

    self%axis(1)%base01 = work_vec12
    self%axis(2)%base01 = work_vec21

    do n = 1, 2
      self%axis(n)%base02 = self%axis(n)%position   &
                                .cross.             &
                            self%axis(n)%base01
    end do

    wave_length_along_radius = self%radius / num_arm
    self%width = wave_length_along_radius * 0.1_DR

  contains

    function iCalc_average( v1, v2 ) result(ans)
      type(vector_t), intent(in) :: v1, v2
      type(vector_t) :: ans

      ans = v1 + v2
      call ans%amplify(0.5_DR)
    end function iCalc_average

  end subroutine spiral_initialize


  function calc_angle_around_axis( point, axis ) result(angle)
    type(vector_t), intent(in) :: point
    type(spiral_axis_t), intent(in) :: axis
    real(DR) :: angle

    type(vector_t) :: axis_to_point
    real(DR) :: work_x, work_y
    real(DR), parameter :: ALMOST_ZERO = 1.e-12

    axis_to_point = point - axis%position

    if ( axis_to_point%amplitude() < ALMOST_ZERO ) then
      angle = 0.0_DR
    else
      call axis_to_point%normalize
      !             axis%base02
      !                 |
      !                 |    axis_to_point
      !                 |   /
      !                 |  /o
      !                 | /   o  angle
      !                 |/     o 
      !       ----------*-------o---------> axis%base01
      !            axis%position
      work_x = axis_to_point .dot. axis%base01
      work_y = axis_to_point .dot. axis%base02
      angle = atan2( work_y, work_x )
    end if
  end function calc_angle_around_axis


  function find_index_min( n, array ) result(index_min)
    integer, intent(in) :: n
    real(DR), dimension(n), intent(in) :: array
    integer :: index_min

    integer :: i
    real(DR) :: minimum_value_so_far
    
    index_min = 1
    minimum_value_so_far = array(1)
    do i = 2, n
      if ( array(i) < minimum_value_so_far ) then
        index_min = i
        minimum_value_so_far = array(i)
      end if
    end do
  end function find_index_min


  function spiral_dist_to_closest_arm( self, num_arm, point ) result(ans)
    class(spiral_t), intent(in) :: self
    type(vector_t), intent(in) :: point
    integer, intent(in) :: num_arm
    real(DR) :: ans

    integer :: m, n, closest_arm_id
    integer :: id_near, id_far   ! 1 or 2
    real(DR), dimension(2) :: dists
    real(DR), dimension(2) :: angles
    real(DR) :: over_run_radius, angle, work_tau, work_dist
    real(DR) :: smoother_factor = 1.5_DR
    real(DR), dimension(num_arm) :: possible_taus
    real(DR), dimension(num_arm) :: dist_to_arm
    
    do n = 1, 2
      dists(n) = acos( point .dot.               &
                       self%axis(n)%position )
      angles(n) = calc_angle_around_axis( point, &
                                          self%axis(n) )
    end do

    id_near = find_index_min( 2, dists )
    id_far  = iToggle( id_near )

    call iSpecial_care

    angle = angles( id_near )
    possible_taus(:) = spiral__tau_from_angle( num_arm, angle )

    do m = 1, num_arm
      work_tau = possible_taus(m)
      work_dist = spiral__distance( self%radius, work_tau )
      dist_to_arm(m) = abs( dists( id_near ) - work_dist )
    end do

    closest_arm_id = find_index_min( num_arm, dist_to_arm )
    ans = dist_to_arm( closest_arm_id )
    
  contains
  
    subroutine iSpecial_care
      !                 |   . 
      !                 |.   *  .
      !                .  *  . 
      !               . |*  .
      !              .  *  .
      !   --F-----------*-----------N--
      !      \       .  *XX <-- Regeion where
      !       \     .  *|XX <-- you need a
      !self%radius .  * |X  <-- special care.
      !         \.  *  .|
      !        . \*  .  | distance between F and X
      !          / .    |           = dist(id_far) 
      !         / /       and angle(id_near) > 0
      !  self%wid/th       
      over_run_radius = dists(id_far) - self%radius
      if ( over_run_radius > 0.0_DR                      &
                         .and.                           &
           over_run_radius < self%width*smoother_factor  &
                         .and.                           &
           angles(id_near) > 0.0_DR ) then
        id_near = iToggle( id_near ) ! switch
        id_far  = iToggle( id_near )
      end if
    end subroutine iSpecial_care

    function iToggle( id )
      integer, intent(in) :: id
      integer :: iToggle
      iToggle = 3 - id    ! {1,2} <==> {2,1}
    end function iToggle

  end function spiral_dist_to_closest_arm


! Private
!-----------
! Public

  function spiral__distance( distance_max, tau )
    real(DR), intent(in) :: distance_max, tau
    real(DR) :: spiral__distance
    spiral__distance = tau * distance_max
  end function spiral__distance


  function spiral__tau_from_distance( distance_max,  &
                                       distance ) result(tau)
    real(DR), intent(in) :: distance_max, distance
    real(DR) :: tau
    tau = distance / distance_max
  end function spiral__tau_from_distance


  function spiral__angle( num_arm, tau )
    integer, intent(in) :: num_arm
    real(DR), intent(in) :: tau
    real(DR) :: spiral__angle
    spiral__angle = TWOPI * num_arm * tau
  end function spiral__angle


  function spiral__tau_from_angle( num_arm, angle ) result(taus)
    integer, intent(in) :: num_arm
    real(DR), intent(in) :: angle  ! supposing [-pi:pi]
    real(DR), dimension(num_arm) :: taus

    real(DR) :: angle_  ! [0:twopi]
    integer :: n

    angle_ = angle  ! to make it [0:twopi]
    if ( angle_ < 0.0_DR ) angle_ = angle_ + TWOPI

    do n = 1, num_arm
      taus(n) = ( angle_ + TWOPI*(n-1) ) / ( TWOPI * num_arm )
    end do
  end function spiral__tau_from_angle

end module spiral_m
