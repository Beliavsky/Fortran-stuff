module diffusion_mod
use kind_mod, only: dp
use random_mod, only: random_normal
implicit none
private
public :: zero_cross
contains
function zero_cross(n) result(k)
! Given n normally distributed increments, return in k the first time for 
! which there was a value both above and below the initial value. Return
! n+1 if this never occurred.
integer, intent(in) :: n
integer :: k
real(kind=dp) :: x
logical :: init_is_min, init_is_max
logical, parameter :: debug=.false.
init_is_min = .true.
init_is_max = .true.
x = 0.0_dp
if (debug) print*
do k=1,n
   x = x + random_normal()
   if (x < 0.0_dp) then
      init_is_min = .false.
   else
      init_is_max = .false.
   end if
   if (debug) print "(f12.4)", x
   if (.not. init_is_min .and. .not. init_is_max) then
      if (debug) print*,k
      return
   end if
end do
if (debug) print*,k
end function zero_cross
end module diffusion_mod
