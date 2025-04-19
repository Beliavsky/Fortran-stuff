module product_mod
use kind_mod, only: dp
implicit none
private
public :: prod, dp
contains
elemental function prod(i, j, bad) result(ij)
! return the product of i and j, unless it overflows, in which case return bad
integer, intent(in) :: i, j, bad
integer             :: ij
real(kind=dp)       :: xij
xij = real(i, kind=dp) * real(j, kind=dp)
if (abs(xij) > huge(i)) then
   ij = bad
else
   ij = i*j
end if
end function prod
end module product_mod
