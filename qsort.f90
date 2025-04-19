module qsort_mod
use kind_mod, only: dp
implicit none
private
public :: indexx
interface indexx
   module procedure indexx_real,indexx_int
end interface indexx
contains
function indexx_int(xx) result(iord)
! return iord(:) such that xx(iord) is in ascending order
integer, intent(in) :: xx(:)
integer             :: iord(size(xx))
iord = indexx_real(1.0_dp*xx)
end function indexx_int
!
function indexx_real(xx) result(iord)
! return iord(:) such that xx(iord) is sorted
real(kind=dp), intent(in) :: xx(:)
integer                   :: iord(size(xx))
real(kind=dp)             :: xcp(size(xx))
xcp = xx
call quick_sort(xcp,iord)
end function indexx_real
!
recursive subroutine quick_sort(list, order)
! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which gives
! the positions of the elements in the original order.
real(kind=dp), intent(in out) :: list(:)
integer      , intent(out)    :: order(:)
! Local variable
integer :: i
do i = 1, size(list)
  order(i) = i
end do
if (any(isnan(list))) return
call quick_sort_1(1, size(list))
contains
recursive subroutine quick_sort_1(left_end, right_end)
integer, intent(in) :: left_end, right_end
!     Local variables
integer             :: i, j, itemp
real(kind=dp)                :: reference, temp
integer, parameter  :: max_simple_sort_size = 6
if (right_end < left_end + max_simple_sort_size) then
  ! Use interchange sort for small lists
  call interchange_sort(left_end, right_end)
else
  ! Use partition ("quick") sort
  reference = list((left_end + right_end)/2)
  i = left_end - 1; j = right_end + 1
  do
    ! Scan list from left end until element >= reference is found
    do
      i = i + 1
!      if (i > size(list)) exit ! added 01/03/2021 10:54 AM by Vivek Rao to avoid out-of-bounds error
      if (list(i) >= reference) exit
    end do
    ! Scan list from right end until element <= reference is found
    do
      j = j - 1
!      if (j < 1) exit ! added 01/03/2021 10:55 AM by Vivek Rao to avoid out-of-bounds error
      if (list(j) <= reference) exit
    end do
    if (i < j) then
      ! Swap two out-of-order elements
      temp = list(i); list(i) = list(j); list(j) = temp
      itemp = order(i); order(i) = order(j); order(j) = itemp
    else if (i == j) then
      i = i + 1
      exit
    else
      exit
    end if
  end do
  if (left_end < j) call quick_sort_1(left_end, j)
  if (i < right_end) call quick_sort_1(i, right_end)
end if
end subroutine quick_sort_1

subroutine interchange_sort(left_end, right_end)
integer, intent(in) :: left_end, right_end
!     Local variables
integer             :: i, j, itemp
real(kind=dp)                :: temp
do i = left_end, right_end - 1
  do j = i+1, right_end
    if (list(i) > list(j)) then
      temp = list(i); list(i) = list(j); list(j) = temp
      itemp = order(i); order(i) = order(j); order(j) = itemp
    end if
  end do
end do
end subroutine interchange_sort
end subroutine quick_sort
end module qsort_mod
