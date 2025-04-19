module champ_mod
! functions related to the Champernowne distribution
use kind_mod, only: dp
implicit none
private
public :: champ_integral, champ_integral_cmplx
contains
elemental function champ_integral(x, c) result(y)
! valid for abs(c) < 2
real(kind=dp), intent(in) :: x, c
real(kind=dp)             :: y
real(kind=dp)             :: d
d = sqrt(4.0_dp - c**2)
y = 2 * atan((c + 2*exp(x)) / d) / d
end function champ_integral
!
elemental function champ_integral_cmplx(x, c) result(y)
complex(kind=dp), intent(in) :: x, c
complex(kind=dp)             :: y
complex(kind=dp)             :: d
d = sqrt(4.0_dp - c**2)
y = 2 * atan((c + 2*exp(x)) / d) / d
end function champ_integral_cmplx
end module champ_mod