module stats_mod
use kind_mod, only: dp
implicit none
public :: mode, mean, median, sd, bin_counts, acf
interface median
   module procedure median_real, median_int
end interface median
integer, parameter :: bad_int = -999
contains
pure function mean(x) result(y)
real(kind=dp), intent(in) :: x(:)
real(kind=dp)             :: y
y = sum(x)/max(1,size(x))
end function mean
!
pure function sd(x) result(y)
! sample standard deviation
real(kind=dp), intent(in) :: x(:)
real(kind=dp)             :: y
real(kind=dp)             :: xmean
integer                   :: n
n = size(x)
if (n < 2) then
   y = -1.0_dp
   return
end if
xmean = sum(x)/n
y = sqrt(sum((x-xmean)**2)/(n-1))
end function sd
!
pure function mode(ivec) result(imode)
integer, intent(in) :: ivec(:)
integer             :: imode
integer             :: i,freq_max,freq,imin
if (size(ivec) < 1) then
   imode = bad_int
   return
end if
imin = minval(ivec)
imode = imin
freq_max = count(ivec == imode)
do i=minval(ivec)+1,maxval(ivec)
   freq = count(ivec == i)
   if (freq > freq_max) then
      freq_max = freq
      imode    = i
   end if
end do
end function mode
!
function median_real(xx) result(xmed)
! return the median of xx(:)
real(kind=dp), intent(in) :: xx(:)
real(kind=dp)             :: xmed
real(kind=dp)             :: xcopy(size(xx))
xcopy = xx
call median_sub(xcopy,size(xx),xmed)
end function median_real
!
function median_int(xx) result(xmed)
! return the median of xx(:)
integer, intent(in) :: xx(:)
real(kind=dp)       :: xmed
xmed = median_real(real(xx,kind=dp))
end function median_int
!
subroutine median_sub(x,n,xmed)
! Find the median of X(1), ... , X(N), using as much of the quicksort
! algorithm as is needed to isolate it.
! N.B. On exit, the array X is partially ordered.
! By Alan Miller
!     Latest revision - 26 November 1996
implicit none
integer, intent(in)  :: n
real(kind=dp), intent(in out) :: x(:)
real(kind=dp), intent(out)    :: xmed
! Local variables
real(kind=dp)    :: temp, xhi, xlo, xmax, xmin
logical :: odd
integer :: hi, lo, nby2, nby2p1, mid, i, j, k
nby2 = n / 2
nby2p1 = nby2 + 1
odd = .true.
!     HI & LO are position limits encompassing the median.
if (n == 2 * nby2) odd = .false.
lo = 1
hi = n
if (n < 3) then
  if (n < 1) then
    xmed = 0.0_dp
    return
  end if
  xmed = x(1)
  if (n == 1) return
  xmed = 0.5_dp*(xmed + x(2))
  return
end if

!     Find median of 1st, middle & last values.

10 mid = (lo + hi)/2
xmed = x(mid)
xlo = x(lo)
xhi = x(hi)
if (xhi < xlo) then          ! Swap xhi & xlo
  temp = xhi
  xhi = xlo
  xlo = temp
end if
if (xmed > xhi) then
  xmed = xhi
else if (xmed < xlo) then
  xmed = xlo
end if

! The basic quicksort algorithm to move all values <= the sort key (XMED)
! to the left-hand end, and all higher values to the other end.

i = lo
j = hi
50 do
  if (x(i) >= xmed) exit
  i = i + 1
end do
do
  if (x(j) <= xmed) exit
  j = j - 1
end do
if (i < j) then
  temp = x(i)
  x(i) = x(j)
  x(j) = temp
  i = i + 1
  j = j - 1

!     Decide which half the median is in.

  if (i <= j) go to 50
end if

if (.not. odd) then
  if (j == nby2 .and. i == nby2p1) go to 130
  if (j < nby2) lo = i
  if (i > nby2p1) hi = j
  if (i /= j) go to 100
  if (i == nby2) lo = nby2
  if (j == nby2p1) hi = nby2p1
else
  if (j < nby2p1) lo = i
  if (i > nby2p1) hi = j
  if (i /= j) go to 100

! Test whether median has been isolated.

  if (i == nby2p1) return
end if
100 if (lo < hi - 1) go to 10

if (.not. odd) then
  xmed = 0.5_dp*(x(nby2) + x(nby2p1))
  return
end if
temp = x(lo)
if (temp > x(hi)) then
  x(lo) = x(hi)
  x(hi) = temp
end if
xmed = x(nby2p1)
return

! Special case, N even, J = N/2 & I = J + 1, so the median is
! between the two halves of the series.   Find max. of the first
! half & min. of the second half, then average.

130 xmax = x(1)
do k = lo, j
  xmax = max(xmax, x(k))
end do
xmin = x(n)
do k = i, hi
  xmin = Min(xmin, x(k))
end do
xmed = 0.5_dp*(xmin + xmax)
end subroutine median_sub
!
function bin_counts(x, thresh) result(counts)
! return # of elements of x(:) in bins defined by thresholds
real(kind=dp), intent(in) :: x(:)
real(kind=dp), intent(in) :: thresh(:) ! must be ascending
integer                   :: counts(size(thresh) + 1)
integer                   :: i, j, n, nthresh
n = size(x)
nthresh = size(thresh)
counts = 0
do i=1,n
   do_j: do j=1,nthresh
      if (x(i) < thresh(j)) exit do_j
   end do do_j
   counts(j) = counts(j) + 1
end do
end function bin_counts
!
function acf(x, nacf) result(xacf)
        integer, intent(in) :: nacf
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp) :: xacf(nacf)
        integer :: i, n
        real(kind=dp) :: mean_x, var_x
        n = size(x)
        mean_x = sum(x) / n
        var_x = sum((x - mean_x)**2) / n

        do i = 1, nacf
            xacf(i) = sum((x(1:n-i) - mean_x) * (x(i+1:) - mean_x)) / ((n - i) * var_x)
        end do
    end function acf
end module stats_mod
