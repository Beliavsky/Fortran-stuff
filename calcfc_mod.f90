module calcfc_mod
private
public :: calcfc
contains
subroutine calcfc (n, m, x, f, con)
use common_nprob, only: nprob
implicit none

integer, parameter :: dp = selected_real_kind(14, 60)

integer, intent(in)     :: n  ! size of x(:)
integer, intent(in)     :: m  ! size of con(:) -- # of constraints
real (dp), intent(in)   :: x(:)   ! function argument
real (dp), intent(out)  :: f      ! function value
real (dp), intent(out)  :: con(:) ! (m) constraint values

if (nprob == 1) then

!  test problem 1 (simple quadratic)

  f = 10.0_dp*(x(1) + 1.0_dp)**2 + x(n)**2
else if (nprob == 2) then

!    test problem 2 (2d unit circle calculation)

  f = x(1)*x(2)
  con(m) = 1.0_dp - x(1)**2 - x(n)**2
else if (nprob == 3) then

!  test problem 3 (3d ellipsoid calculation)

  f = x(1)*x(2)*x(3)
  con(1) = 1.0_dp - x(1)**2 - 2.0_dp*x(2)**2 - 3.0_dp*x(n)**2
else if (nprob == 4) then

!  test problem 4 (weak rosenbrock)

  f = (x(1)**2 - x(2))**2 + (1.0_dp + x(1))**2
else if (nprob == 5) then

!  test problem 5 (intermediate rosenbrock)

  f = 10.0_dp*(x(1)**2 - x(n))**2 + (1.0_dp + x(1))**2
else if (nprob == 6) then

!  test problem 6 (equation (9.1.15) in fletcher's book)

  f = - x(1) - x(2)
  con(1) = x(2) - x(1)**2
  con(2) = 1.0_dp - x(1)**2 - x(2)**2
else if (nprob == 7) then

!  test problem 7 (equation (14.4.2) in fletcher's book)

  f = x(3)
  con(1) = 5.0_dp*x(1) - x(2) + x(3)
  con(2) = x(3) - x(1)**2 - x(2)**2 - 4.0_dp*x(2)
  con(m) = x(3) - 5.0_dp*x(1) - x(2)
else if (nprob == 8) then

!  test problem 8 (rosen-suzuki)

  f = x(1)**2 + x(2)**2 + 2.0*x(3)**2 + x(4)**2 - 5.0_dp*x(1) - 5.0_dp*x(2)  &
      - 21.0_dp*x(3) + 7.0_dp*x(4)
  con(1) = 8.0_dp - x(1)**2 - x(2)**2 - x(3)**2 - x(4)**2 - x(1) + x(2)  &
           - x(3) + x(4)
  con(2) = 10._dp - x(1)**2 - 2._dp*x(2)**2 - x(3)**2 - 2._dp*x(4)**2 + x(1) + x(4)
  con(m) = 5.0_dp - 2.0*x(1)**2 - x(2)**2 - x(3)**2 - 2.0_dp*x(1) + x(2) + x(4)
else if (nprob == 9) then

!  test problem 9 (hock and schittkowski 100)

  f = (x(1) - 10._dp)**2 + 5._dp*(x(2) - 12._dp)**2 + x(3)**4 + 3._dp*(x(4)  &
       - 11._dp)**2 + 10._dp*x(5)**6 + 7._dp*x(6)**2 + x(7)**4 - 4._dp*x(6)*x(7) &
       - 10._dp*x(6) - 8._dp*x(7)
  con(1) = 127._dp - 2._dp*x(1)**2 - 3._dp*x(2)**4 - x(3) - 4._dp*x(4)**2   &
           - 5._dp*x(5)
  con(2) = 282._dp - 7._dp*x(1) - 3._dp*x(2) - 10._dp*x(3)**2 - x(4) + x(5)
  con(3) = 196._dp - 23._dp*x(1) - x(2)**2 - 6._dp*x(6)**2 + 8._dp*x(7)
  con(4) = - 4._dp*x(1)**2 - x(2)**2 + 3._dp*x(1)*x(2) - 2._dp*x(3)**2 - 5._dp*x(6)  &
           + 11._dp*x(7)
else if (nprob == 10) then

!  test problem 10 (hexagon area)

  f = - 0.5_dp*(x(1)*x(4) - x(2)*x(3) + x(3)*x(n) - x(5)*x(n) + x(5)*x(8) &
                - x(6)*x(7))
  con(1) = 1.0_dp - x(3)**2 - x(4)**2
  con(2) = 1.0_dp - x(n)**2
  con(3) = 1.0_dp - x(5)**2 - x(6)**2
  con(4) = 1.0_dp - x(1)**2 - (x(2) - x(n))**2
  con(5) = 1.0_dp - (x(1) - x(5))**2 - (x(2) - x(6))**2
  con(6) = 1.0_dp - (x(1) - x(7))**2 - (x(2) - x(8))**2
  con(7) = 1.0_dp - (x(3) - x(5))**2 - (x(4) - x(6))**2
  con(8) = 1.0_dp - (x(3) - x(7))**2 - (x(4) - x(8))**2
  con(9) = 1.0_dp - x(7)**2 - (x(8) - x(n))**2
  con(10) = x(1)*x(4) - x(2)*x(3)
  con(11) = x(3)*x(n)
  con(12) = - x(5)*x(n)
  con(13) = x(5)*x(8) - x(6)*x(7)
  con(m) = x(n)
end if

return
end subroutine calcfc
end module calcfc_mod