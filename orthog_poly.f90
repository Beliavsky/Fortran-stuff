module orthog_poly_mod
use kind_mod, only: dp
implicit none
private
public :: hermite,hermite_basis
contains
function hermite_basis(x,norder) result(y)
! return a Hermite polynomial basis at points x(:) up to order norder
real(kind=dp), intent(in) :: x(:)
integer      , intent(in) :: norder
real(kind=dp)             :: y(size(x),norder)
integer                   :: iorder
do iorder=1,norder
   y(:,iorder) = hermite(x,iorder)
end do
end function hermite_basis
!
elemental function hermite(x,norder) result(y)
real(kind=dp), intent(in) :: x
integer      , intent(in) :: norder ! order of hermite polynomial
real(kind=dp)             :: y
select case (norder)
   case  (0); y = 1.0_dp
   case  (1); y = x
   case  (2); y = x**2  - 1
   case  (3); y = x**3  - 3*x
   case  (4); y = x**4  - 6*x**2  +  3
   case  (5); y = x**5  - 10*x**3 +  15*x
   case  (6); y = x**6  - 15*x**4 +  45*x**2 -   15
   case  (7); y = x**7  - 21*x**5 + 105*x**3 -  105*x
   case  (8); y = x**8  - 28*x**6 + 210*x**4 -  420*x**2 +  105
   case  (9); y = x**9  - 36*x**7 + 378*x**5 - 1260*x**3 +  945*x
   case (10); y = x**10 - 45*x**8 + 630*x**6 - 3150*x**4 + 4725*x**2 - 945
end select
end function hermite
end module orthog_poly_mod
