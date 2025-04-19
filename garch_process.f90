module garch_process_mod
use kind_mod, only: dp
private
public :: garch_sim,gjr_garch_sim,agarch_sim
contains
function garch_sim(e,w,a,b,var0) result(y)
! simulate a GARCH(1,1) process
real(kind=dp), intent(in) :: e(:)       ! noise with mean 0 and variance 1
real(kind=dp), intent(in) :: w          ! constant term in variance
real(kind=dp), intent(in) :: a          ! weight of previous squared return in variance
real(kind=dp), intent(in) :: b          ! weight of previous variance in variance
real(kind=dp), intent(in) :: var0       ! initial variance
real(kind=dp)             :: y(size(e)) ! simulated observations
real(kind=dp)             :: yvar,yvar_old
integer                   :: i,n
n = size(e)
if (n < 1) return
y(1) = e(1)*sqrt(var0)
yvar_old = var0
do i=2,n
   yvar     = w + a*y(i-1)**2 + b*yvar_old ! updated conditional variance
   y(i)     = sqrt(yvar)*e(i)              ! observation
   yvar_old = yvar
end do
end function garch_sim
!
function gjr_garch_sim(e,w,a,adn,b,var0) result(y)
! simulate a GJR-GARCH(1,1) process
real(kind=dp), intent(in) :: e(:)       ! noise with mean 0 and variance 1
real(kind=dp), intent(in) :: w          ! constant term in variance
real(kind=dp), intent(in) :: a          ! weight of previous squared return in variance
real(kind=dp), intent(in) :: adn        ! extra weight of previous squared return in variance when return is negative
real(kind=dp), intent(in) :: b          ! weight of previous variance in variance
real(kind=dp), intent(in) :: var0       ! initial variance
real(kind=dp)             :: y(size(e)) ! simulated observations
real(kind=dp)             :: yvar,yvar_old
integer                   :: i,n
n = size(e)
if (n < 1) return
y(1) = e(1)*sqrt(var0)
yvar_old = var0
do i=2,n
   yvar     = w + merge(a+adn,a,y(i-1)<0)*y(i-1)**2 + b*yvar_old ! updated conditional variance
   y(i)     = sqrt(yvar)*e(i)              ! observation
   yvar_old = yvar
end do
end function gjr_garch_sim
!
function agarch_sim(e,w,a,yshift,adn,b,var0) result(y)
! simulate an AGARCH(1,1) process with shift and twist effects
real(kind=dp), intent(in) :: e(:)       ! noise with mean 0 and variance 1
real(kind=dp), intent(in) :: w          ! constant term in variance
real(kind=dp), intent(in) :: a          ! weight of previous squared return in variance
real(kind=dp), intent(in) :: adn        ! extra weight of previous squared return in variance when return is negative
real(kind=dp), intent(in) :: b          ! weight of previous variance in variance
real(kind=dp), intent(in) :: var0       ! initial variance
real(kind=dp), intent(in) :: yshift     ! shift of return that is squared in news impact curve
real(kind=dp)             :: y(size(e)) ! simulated observations
real(kind=dp)             :: yvar,yvar_old
integer                   :: i,n
n = size(e)
if (n < 1) return
y(1) = e(1)*sqrt(var0)
yvar_old = var0
do i=2,n
   yvar     = w + merge(a+adn,a,y(i-1)<0)*(y(i-1)-yshift)**2 + b*yvar_old ! updated conditional variance
   y(i)     = sqrt(yvar)*e(i)              ! observation
   yvar_old = yvar
end do
end function agarch_sim
end module garch_process_mod
