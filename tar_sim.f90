module tar_sim_mod
use kind_mod  , only: dp
use random_mod, only: random_normal
implicit none
private
public :: tar1_sim
contains
function tar1_sim(n,ar1_up,ar1_dn,xsd,x1,xconst) result(xar)
! simulate a TAR1 (threshold autoregression of order 1) process with Gaussian innovations
integer      , intent(in)           :: n              ! # of observations
real(kind=dp), intent(in)           :: ar1_up,ar1_dn  ! AR(1) coefficients in up and down regimes
real(kind=dp)                       :: xar(n)         ! simulated time series
real(kind=dp), intent(in), optional :: xsd            ! standard deviation of noise
real(kind=dp), intent(in), optional :: x1             ! initial observation
real(kind=dp), intent(in), optional :: xconst         ! unconditional mean
real(kind=dp)                       :: xsd_,xconst_
integer                             :: i
if (present(xsd)) then
   xsd_ = xsd
else
   xsd_ = 1.0_dp
end if
if (present(xconst)) then
   xconst_ = xconst
else
   xconst_ = 0.0_dp
end if
xar  = random_normal(n)
if (n < 1) return
if (present(x1)) then
   xar(1) = x1
else
   xar(1) = xconst_
end if
do i=2,n
   xar(i) = xsd_*xar(i) + merge(ar1_up,ar1_dn,xar(i-1) > 0) * xar(i-1) + xconst_
end do
end function tar1_sim
end module tar_sim_mod
