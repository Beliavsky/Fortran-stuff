module trig_reg_mod
use kind_mod     , only: dp
use constants_mod, only: pi
implicit none
private
public :: poly_trig_basis
contains
function poly_trig_basis(xx,npow,nk) result(basis)
real(kind=dp), intent(in)  :: xx(:)
integer      , intent(in)  :: npow,nk
real(kind=dp)              :: basis(size(xx),npow+2*nk)
integer                    :: ik,ipow,j
forall (ipow=1:npow) basis(:,ipow) = xx**ipow
do ik=1,nk
   j            = npow + 2*ik - 1
   basis(:,j)   = sin(ik*pi*xx)
   basis(:,j+1) = cos(ik*pi*xx)
end do
end function poly_trig_basis 
end module trig_reg_mod
