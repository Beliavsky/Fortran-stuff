module sv_sim_mod
use kind_mod  , only: dp
use random_mod, only: random_normal
private
public :: arsv_sim
contains
subroutine arsv_sim(sd_eq,sd_ar1,sd_sd,xx,xsd)
! simulate from an autoregressive log stochastic volatility model with zero correlation 
! between return and volatility innovations
real(kind=dp), intent(in)  :: sd_eq  ! equilibrium standard deviation
real(kind=dp), intent(in)  :: sd_ar1 ! 1st order autoregressive coefficient of log standard deviation
real(kind=dp), intent(in)  :: sd_sd  ! standard deviation of standard deviation
real(kind=dp), intent(out) :: xsd(:) ! simulated standard deviations
real(kind=dp), intent(out) :: xx(:)  ! simulated time series with stochastic volatility
real(kind=dp)              :: log_sd,log_sd_old,log_sd_eq
integer                    :: i,n
n = size(xx)
if (size(xsd) /= n) then
   write (*,*) "in arsv_sim, size(xsd), size(xx) =",size(xsd),n," must be equal, STOPPING"
   stop
end if
if (n < 1) return
log_sd_eq = log(sd_eq)
do i=1,n
   if (i > 1) then
      log_sd_old = log(xsd(i-1))
      log_sd = (1-sd_ar1)*log_sd_eq + sd_ar1*log_sd_old + sd_sd*random_normal()
      xsd(i) = exp(log_sd)
   else
      xsd(1) = sd_eq
   end if   
   xx(i)  = xsd(i)*random_normal()
end do
end subroutine arsv_sim
end module sv_sim_mod
