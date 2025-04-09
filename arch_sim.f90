module arch_sim_mod
use kind_mod  , only: dp
use random_mod, only: random_normal
use util_mod  , only: assert_equal
implicit none
private
public :: arch_1_sim,arch_sim
contains
function arch_1_sim(n,wgt_1,vol_const,x0) result(xx)
! simulate from an ARCH(1) process with normal innovations
integer      , intent(in)  :: n          ! # of observations to simulate
real(kind=dp), intent(in)  :: wgt_1      ! ARCH(1) coefficient -- weight on previous squared observation
real(kind=dp), intent(in)  :: vol_const  ! constant term in variance prediction
real(kind=dp), intent(in)  :: x0         ! pre-sample observation
real(kind=dp)              :: xx(n)      ! simulated ARCH time series
real(kind=dp)              :: ee(n)
real(kind=dp)              :: vol
integer                    :: i
ee = random_normal(n)
vol = sqrt(max(0.0_dp,vol_const**2 + wgt_1*x0**2))
xx(1) = vol*ee(1)
do i=2,n
   vol   = sqrt(max(0.0_dp,vol_const**2 + wgt_1*xx(i-1)**2))
   xx(i) = vol*ee(i)
end do
end function arch_1_sim
!
function arch_sim(n,wgt,vol_const,vol_init) result(xx)
! simulate from an ARCH(n) process with normal innovations
integer      , intent(in)  :: n          ! # of observations to simulate
real(kind=dp), intent(in)  :: wgt(:)     ! ARCH coefficients -- weight on previous squared observation
real(kind=dp), intent(in)  :: vol_const  ! constant term in variance prediction
real(kind=dp), intent(in)  :: vol_init   ! volatility for first size(wgt) observations
real(kind=dp)              :: xx(n)      ! simulated ARCH time series
real(kind=dp)              :: ee(n)
real(kind=dp)              :: vol,past_sum_squares
integer                    :: i,nwgt
nwgt = size(wgt)
ee = random_normal(n)
do i=1,n
   if (i <= nwgt) then
      vol = vol_init
   else
      past_sum_squares = sum(wgt(nwgt:1:-1)*xx(i-nwgt:i-1)**2)
      vol = sqrt(max(0.0_dp,vol_const**2 + past_sum_squares))
   end if
   xx(i) = vol*ee(i)
end do
end function arch_sim
end module arch_sim_mod
