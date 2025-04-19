module convolution_mod
use kind_mod, only: dp
private
public :: convolution
contains
pure function convolution(xx,wgt) result(yy)
real(kind=dp), intent(in) :: xx(:)        ! (nx)
real(kind=dp), intent(in) :: wgt(:)       ! (nwgt)
real(kind=dp)             :: yy(size(xx)) ! (nx) 
integer                   :: i,nx,nwgt
nx   = size(xx)
nwgt = size(wgt)
do i=1,nx
   if (i < nwgt) then
      yy(i) = sum(xx(:i)*wgt(:i))
   else
      yy(i) = sum(xx(i-nwgt+1:i)*wgt)
   end if
end do
end function convolution
end module convolution_mod
