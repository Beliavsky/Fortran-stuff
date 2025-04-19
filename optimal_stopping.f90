module optimal_stopping_mod
implicit none
private
public :: secretary
contains
function secretary(iperm) result(ival)
integer, intent(in) :: iperm(:) ! successive observations, with values from 1 to size(iperm)
integer             :: ival(size(iperm)) ! value obtained by "leaping" after various numbers of observations
integer             :: i, ibest(size(iperm)), nobs
logical             :: best_so_far
nobs = size(iperm)
if (nobs < 1) return
ival       = 0
ibest(1)   = iperm(1)
ival(1)    = iperm(1)
do i=2,nobs ! loop over observations
   best_so_far = iperm(i) > ibest(i-1)
   ibest(i) = merge(iperm(i), ibest(i-1), best_so_far)
   if (i == nobs .or. best_so_far) where (ival(2:i) == 0) ival(2:i) = iperm(i)
end do
end function secretary
end module optimal_stopping_mod