module random
! Derived from https://jblevins.org/mirror/amiller/random.f90 of Alan Miller
implicit none
public
integer, parameter :: dp = selected_real_kind(12, 60)
real, parameter, private :: half = 0.5, s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,    &
                            r1 = 0.27597, r2 = 0.27846
interface random_normal
   module procedure random_normal_scalar,random_normal_vec
end interface
contains
!
function random_normal_vec(n) result(vec)
integer, intent(in) :: n
real(kind=dp)       :: vec(n)
integer             :: i
do i=1,n
   vec(i) = random_normal_scalar()
end do
end function random_normal_vec
!
function random_normal_scalar() result(fn_val)

! adapted from the following fortran 77 code
!      algorithm 712, collected algorithms from acm.
!      this work published in transactions on mathematical software,
!      vol. 18, no. 4, december, 1992, pp. 434-435.

!  the function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  the algorithm uses the ratio of uniforms method of a.j. kinderman
!  and j.f. monahan augmented with quadratic bounding curves.

real :: fn_val

!     local variables
real :: u, v, x, y, q

!     generate p = (u,v) uniform in rectangle enclosing acceptance region

do
  call random_number(u)
  call random_number(v)
  v = 1.7156 * (v - half)

!     evaluate the quadratic form
  x = u - s
  y = abs(v) - t
  q = x**2 + y*(a*y - b*x)

!     accept p if inside inner ellipse
  if (q < r1) exit
!     reject p if outside outer ellipse
  if (q > r2) cycle
!     reject p if outside acceptance region
  if (v**2 < -4.0*log(u)*u**2) exit
end do

!     return ratio of p's coordinates as the normal deviate
fn_val = v/u
end function random_normal_scalar
!
function random_normal_vec_inline(n) result(fn_val)
integer, intent(in) :: n
real                :: fn_val(n)
!     local variables
real :: u, v, x, y, q
integer :: i
!     generate p = (u,v) uniform in rectangle enclosing acceptance region

do i=1,n
do
  call random_number(u)
  call random_number(v)
  v = 1.7156 * (v - half)
!     evaluate the quadratic form
  x = u - s
  y = abs(v) - t
  q = x**2 + y*(a*y - b*x)
!     accept p if inside inner ellipse
  if (q < r1) exit
!     reject p if outside outer ellipse
  if (q > r2) cycle
!     reject p if outside acceptance region
  if (v**2 < -4.0*log(u)*u**2) exit
end do
!     return ratio of p's coordinates as the normal deviate
fn_val(i) = v/u
end do
end function random_normal_vec_inline
end module random

