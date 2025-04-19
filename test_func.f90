module test_func_mod
use kind_mod     , only: dp
use constants_mod, only: pi
implicit none
private
public :: test_func,test_func_scalar
real (kind=dp), parameter :: bad_real = -999.0_dp
contains
function test_func(xx,fname) result(yy) 
! test functions for nonlinear multiple regression
real(kind=dp), intent(in)     :: xx(:,:)
character (len=*), intent(in) :: fname
real(kind=dp)                 :: yy(size(xx,1))
integer                       :: i,j,nobs,nvar
nobs = size(xx,1)
nvar = size(xx,2)
select case (fname)
   case ("tanh")    ; yy = tanh(sum(xx,dim=2))
   case ("exp")     ; yy = exp(sum(xx,dim=2))
   case ("r")       ; yy = sqrt(sum(xx**2,dim=2))
   case ("r^2")     ; yy = sum(xx**2,dim=2)
   case ("sin")     ; yy = sin(pi*sum(xx,dim=2))
   case ("cos")     ; yy = cos(pi*sum(xx,dim=2))
   case ("sum")     ; yy = sum(xx,dim=2)
   case ("weighted"); forall (i=1:nobs) yy(i) = sum([(j*xx(i,j),j=1,nvar)])
   case ("const")   ; yy = 10.0_dp
   case ("hyperbolic"); yy = sqrt(sum(xx**2,dim=2)+1.0_dp)
   case default     ; write (*,*) "invalid function name " // trim(fname); stop
end select
end function test_func
!
elemental function test_func_scalar(xx,fname) result(yy) 
! test functions for nonlinear multiple regression
real(kind=dp), intent(in)     :: xx
character (len=*), intent(in) :: fname
real(kind=dp)                 :: yy
select case (fname)
   case ("x-2*tanh(x)"); yy = xx - 2*tanh(xx)
   case ("tanh")    ; yy = tanh(xx)
   case ("exp")     ; yy = exp(xx)
   case ("r")       ; yy = abs(xx)
   case ("r^2")     ; yy = xx**2
   case ("sin")     ; yy = sin(pi*xx)
   case ("cos")     ; yy = cos(pi*xx)
   case ("const")   ; yy = 10.0_dp
   case ("hyperbolic"); yy = sqrt(xx**2+1.0_dp)
   case default     ; yy = bad_real
end select
end function test_func_scalar
end module test_func_mod
