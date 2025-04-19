module log_lik_mod
use kind_mod, only: dp
implicit none
private
public :: minus_normal_log_lik
contains
elemental function minus_normal_log_lik(x,mu,sigma) result(minus_log_lik)
! -log_likelihood of x given normal distribution with mean mu and standard deviation sigma
real(kind=dp), intent(in) :: x
real(kind=dp), intent(in) :: mu
real(kind=dp), intent(in) :: sigma
real(kind=dp), parameter  :: pi = 4*atan(1.0_dp), half_log_two_pi = log(2*pi)/2
real(kind=dp)             :: minus_log_lik
minus_log_lik = half_log_two_pi + log(sigma) + 0.5_dp*((x-mu)/sigma)**2
end function minus_normal_log_lik
end module log_lik_mod