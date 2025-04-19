module quantile_mod
use kind_mod     , only: dp
use constants_mod, only: two_over_pi,pi_over_2,pi
implicit none
private
public :: quantile_sech,quantile_cauchy,quantile_logistic,quantile_laplace
contains
elemental function quantile_sech(x) result(y)
! quantile of the hyperbolic secant distribution https://en.wikipedia.org/wiki/Hyperbolic_secant_distribution
real(kind=dp), intent(in) :: x
real(kind=dp)             :: y
y = two_over_pi * log(tan(x*pi_over_2))
end function quantile_sech
!
elemental function quantile_logistic(x) result(y)
! quantile of the logistic distribution
real(kind=dp), intent(in) :: x
real(kind=dp)             :: y
y = log(x/(1-x))
end function quantile_logistic
!
elemental function quantile_cauchy(x) result(y)
! quantile of the hyperbolic secant distribution https://en.wikipedia.org/wiki/Hyperbolic_secant_distribution
real(kind=dp), intent(in) :: x
real(kind=dp)             :: y
y = tan(pi*(x-0.5_dp))
end function quantile_cauchy
!
elemental function quantile_laplace(x) result(y)
! quantile of the Laplace distribution
real(kind=dp), intent(in) :: x
real(kind=dp)             :: y
if (x <= 0.5_dp) then
   y = log(2*x)
else
   y = -log(2 - 2*x)
end if
end function quantile_laplace
end module quantile_mod