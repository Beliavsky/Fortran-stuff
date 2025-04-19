module garch_mod
use kind_mod, only: dp
implicit none
contains
subroutine garch_1_1(x,omega,alpha,beta,var1,var)
! given a time series and GARCH(1,1) parameters, compute the conditional variance at each time
real(kind=dp), intent(in)  :: x(:)   ! (n) -- time series of returns
real(kind=dp), intent(in)  :: omega  ! constant term in variance
real(kind=dp), intent(in)  :: alpha  ! coefficient of previous squared return
real(kind=dp), intent(in)  :: beta   ! coefficient of previous variance
real(kind=dp), intent(in)  :: var1   ! initial variance
real(kind=dp), intent(out) :: var(:) ! (n) -- time series of conditional variance
integer                    :: i,n
character (len=*), parameter :: msg="in garch_1_1, "
n = size(x)
if (size(var) /= n) then
   print*,msg,"size(x), size(var) =",n,size(var)," must be equal, STOPPING"
   error stop
else if (n < 1) then
   return
end if
var(1) = var1
do i=2,n
   var(i) = omega + alpha*x(i-1)**2 + beta*var(i-1)
end do
end subroutine garch_1_1
!
subroutine garch_1_1_sim(omega,alpha,beta,noise,var1,x,var)
! simulate from a GARCH(1,1) process
real(kind=dp), intent(in)  :: omega    ! constant term in variance
real(kind=dp), intent(in)  :: alpha    ! coefficient of previous squared return
real(kind=dp), intent(in)  :: beta     ! coefficient of previous variance
real(kind=dp), intent(in)  :: var1     ! initial variance
real(kind=dp), intent(in)  :: noise(:) ! (n) innovations with unit variance
real(kind=dp), intent(out) :: x(:)     ! (n) simulated time series
real(kind=dp), intent(out) :: var(:)   ! (n) time series of conditional variance
integer                    :: i,n
character (len=*), parameter :: msg="in garch_1_1_sim "
n = size(noise)
if (size(x) /= n .or. size(var) /= n) then
   print*,msg,"size(noise), size(x), size(var) =",n,size(x),size(var)," must all be equal, STOPPING"
   error stop
end if
if (n < 1) return
var(1) = var1
x(1) = sqrt(var(1)) * noise(i)
do i=2,n
   var(i) = omega + alpha*x(i-1)**2 + beta*var(i-1)
   x(i)   = sqrt(var(i))*noise(i)
end do
end subroutine garch_1_1_sim
end module garch_mod