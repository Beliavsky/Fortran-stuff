module find_subsets

! A module of routines for finding and recording best-fitting subsets of
! regression variables

!   Version 1.10, 17 February 2004
!   Author: Alan Miller
!           formerly of CSIRO Division of Mathematical & Information Sciences
!   Phone: (+61) 3 9592-5085
!   e-mail: amiller @ bigpond.net.au
!   WWW-page: http://users.bigpond.net.au/amiller/

! 17 Feb 2004 Correction to subroutine EFROYM for the case in which all the
!             variables are selected. Thanks to David Jones.
! 12 Nov 1999 Made changes to routines exadd1 & add2 to prevent the calculation
!             of negative residual sums of squares which could occur in cases
!             in which the true RSS is zero.   Routine seq2 changed to avoid
!             cycling.
! 24 May 2000 Changed lsq_kind to dp (cosmetic change only)
! 4 June 2000 Added routine random_pick which picks a random set of variables.
! 29 Aug 2002 Set value of size in subroutine EFROYM when IER /= 0.
use lsq, only: dp, ss, vorder, ncol, nobs, initialized, tol, rss, rhs, d, r, row_ptr, &
               reordr, vmove, startup, tolset, sing, sserr, includ, regcf
use reg_sub_mod, only: ranord, calc_penalty
implicit none
public :: init_subsets, bakwrd, forwrd, random_pick, xhaust, seqrep, seq2, efroym, &
          lopt, ress, nbest, max_size, shell, replace2, cross_validation
integer, save                 :: max_size, nbest, lopt_dim1
real (dp), allocatable, save  :: bound(:), ress(:,:)
integer, allocatable, save    :: lopt(:,:)


contains

subroutine init_subsets(nvar_max, fit_const, nvar)

integer, intent(in)           :: nvar_max
logical, intent(in)           :: fit_const
integer, intent(in), optional :: nvar

!     Local variables

integer    :: i, ier
real (dp)  :: eps = 1.e-14
logical    :: lindep(ncol)

!     The LSQ module has probably already been initialized, but just in case ..

if (.not. initialized) then
  if (present(nvar)) call startup(nvar, fit_const)
end if

if (fit_const) then
  max_size = nvar_max + 1
else
  max_size = nvar_max
end if

lopt_dim1 = max_size * (max_size + 1) / 2
if (allocated(bound)) deallocate(bound, ress, lopt)
allocate (bound(max_size), ress(max_size,nbest), lopt(lopt_dim1,nbest))

bound = huge(eps)
ress  = huge(eps)
lopt  = 0

call tolset(eps)
call sing(lindep, ier)

call ss()
do i = 1, max_size
  call report(i, rss(i))
end do

return
end subroutine init_subsets


subroutine add1(first, last, ss, smax, jmax, ier)

! Calculate the reduction in residual sum of squares when one variable,
! selected from those in positions FIRST .. LAST, is added in position FIRST,
! given that the variables in positions 1 .. FIRST-1 (if any) are already
! included.

integer, intent(in)     :: first, last
integer, intent(out)    :: jmax, ier
real (dp), intent(out)  :: ss(:), smax

!     Local variables

integer    :: j, inc, pos, row, col
real (dp)  :: zero = 0.0_dp, diag, dy, ssqx, sxx(ncol), sxy(ncol)

!     Check call arguments

jmax = 0
smax = zero
ier = 0
if (first > ncol) ier = 1
if (last < first) ier = ier + 2
if (first < 1) ier = ier + 4
if (last > ncol) ier = ier + 8
if (ier /= 0) return

!     Accumulate sums of squares & products from row FIRST

sxx(first:last) = zero
sxy(first:last) = zero
inc = ncol - last
pos = row_ptr(first)
do row = first, last
  diag = d(row)
  dy = diag * rhs(row)
  sxx(row) = sxx(row) + diag
  sxy(row) = sxy(row) + dy
  do col = row+1, last
    sxx(col) = sxx(col) + diag * r(pos)**2
    sxy(col) = sxy(col) + dy * r(pos)
    pos = pos + 1
  end do
  pos = pos + inc
end do

!     Incremental sum of squares for a variable = SXY * SXY / SXX.
!     Calculate whenever sqrt(SXX) > TOL for that variable.

do j = first, last
  ssqx = sxx(j)
  if (sqrt(ssqx) > tol(j)) then
    ss(j) = sxy(j)**2 / sxx(j)
    if (ss(j) > smax) then
      smax = ss(j)
      jmax = j
    end if
  else
    ss(j) = zero
  end if
end do

return
end subroutine add1


subroutine add2(first, last, smax, j1, j2, ier)

!     Calculate the maximum reduction in residual sum of squares when 2
!     variables, selected from those in positions FIRST .. LAST, are
!     added, given that the variables in positions 1 .. FIRST-1 (if
!     any) are already included.    J1, J2 are the positions of the two
!     best variables.   N.B. J2 < J1.

integer, intent(in)     :: first, last
integer, intent(out)    :: j1, j2, ier
real (dp), intent(out)  :: smax

!     Local variables

integer    :: start, i1, i2, row, pos1, pos2, inc
real (dp)  :: zero = 0.0_dp, temp, det, two = 2.0, sxx(ncol), sxy(ncol), sx1x2

!     Check call arguments

smax = zero
j1 = 0
j2 = 0
ier = 0
if (first > ncol) ier = 1
if (last <= first) ier = ier + 2
if (first < 1) ier = ier + 4
if (last > ncol) ier = ier + 8
if (ier /= 0) return

start = row_ptr(first)

!     Cycle through all pairs of variables from those between FIRST & LAST.

do i1 = first, last
  sxx(i1) = d(i1)
  sxy(i1) = d(i1) * rhs(i1)
  pos1 = start + i1 - first - 1
  do row = first, i1-1
    temp = d(row) * r(pos1)
    sxx(i1) = sxx(i1) + temp*r(pos1)
    sxy(i1) = sxy(i1) + temp*rhs(row)
    pos1 = pos1 + ncol - row - 1
  end do

  do i2 = first, i1-1
    pos1 = start + i1 - first - 1
    pos2 = start + i2 - first - 1
    sx1x2 = zero
    do row = first, i2-1
      sx1x2 = sx1x2 + d(row)*r(pos1)*r(pos2)
      inc = ncol - row - 1
      pos1 = pos1 + inc
      pos2 = pos2 + inc
    end do
    sx1x2 = sx1x2 + d(i2)*r(pos1)

!     Calculate reduction in RSS for pair I1, I2.
!     The sum of squares & cross-products are in:
!              ( SXX(I1)  SX1X2   )      ( SXY(I1) )
!              ( SX1X2    SXX(I2) )      ( SXY(I2) )

    det = max( (sxx(i1) * sxx(i2) - sx1x2**2), zero)
    temp = sqrt(det)
    if (temp < tol(i1)*sqrt(sxx(i2)) .or.             &
        temp < tol(i2)*sqrt(sxx(i1))) cycle
    temp = ((sxx(i2)*sxy(i1) - two*sx1x2*sxy(i2))*sxy(i1) + sxx(i1)*sxy(i2)**2) &
           / det
    if (temp > smax) then
      smax = temp
      j1 = i1
      j2 = i2
    end if
  end do ! i2 = first, i1-1
end do   ! i1 = first, last

return
end subroutine add2


subroutine bakwrd(first, last, ier)
!     Backward elimination from variables in positions FIRST .. LAST.
!     If FIRST > 1, variables in positions prior to this are forced in.
!     If LAST < ncol, variables in positions after this are forced out.
!     On exit, the array VORDER contains the numbers of the variables
!     in the order in which they were deleted.
integer, intent(in)  :: first, last
integer, intent(out) :: ier
!     Local variables
integer    :: pos, jmin, i
real (dp)  :: ss(last), smin

!     Check call arguments

ier = 0
if (first >= ncol) ier = 1
if (last <= 1) ier = ier + 2
if (first < 1) ier = ier + 4
if (last > ncol) ier = ier + 8
if (ier /= 0) return

!     For POS = LAST, ..., FIRST+1 call DROP1 to find best variable to
!     find which variable to drop next.

do pos = last, first+1, -1
  call drop1(first, pos, ss, smin, jmin, ier)
  call exdrop1(first, pos, ss, smin, jmin)
  if (jmin > 0 .and. jmin < pos) then
    call vmove(jmin, pos, ier)
    if (nbest > 0) then
      do i = jmin, pos-1
        call report(i, rss(i))
      end do
    end if
  end if
end do

return
end subroutine bakwrd


subroutine drop1(first, last, ss, smin, jmin, ier)

! Calculate the increase in the residual sum of squares when the variable in
! position J is dropped from the model (i.e. moved to position LAST),
! for J = FIRST, ..., LAST-1.

integer, intent(in)     :: first, last
integer, intent(out)    :: jmin, ier
real (dp), intent(out)  :: ss(:), smin

!     Local variables

integer    :: j, pos1, inc, pos, row, col, i
real (dp)  :: large = huge(1.0_dp), zero = 0.0_dp, d1, rhs1, d2, x, wk(last), &
              vsmall = tiny(1.0_dp)

!     Check call arguments

jmin = 0
smin = large
ier = 0
if (first > ncol) ier = 1
if (last < first) ier = ier + 2
if (first < 1) ier = ier + 4
if (last > ncol) ier = ier + 8
if (ier /= 0) return

!     POS1 = position of first element of row FIRST in r.

pos1 = row_ptr(first)
inc = ncol - last

!     Start of outer cycle for the variable to be dropped.

do j = first, last
  d1 = d(j)
  if (sqrt(d1) < tol(j)) then
    ss(j) = zero
    smin = zero
    jmin = j
    go to 50
  end if
  rhs1 = rhs(j)
  if (j == last) go to 40

!     Copy row J of R into WK.

  pos = pos1
  do i = j+1, last
    wk(i) = r(pos)
    pos = pos + 1
  end do
  pos = pos + inc

!     Lower the variable past each row.

  do row = j+1, last
    x = wk(row)
    d2 = d(row)
    if (abs(x) * sqrt(d1) < tol(row) .or. d2 < vsmall) then
      pos = pos + ncol - row
      cycle
    end if
    d1 = d1 * d2 / (d2 + d1 * x**2)
    do col = row+1, last
      wk(col) = wk(col) - x * r(pos)
      pos = pos + 1
    end do
    rhs1 = rhs1 - x * rhs(row)
    pos = pos + inc
  end do
  40 ss(j) = rhs1 * d1 * rhs1
  if (ss(j) < smin) then
    jmin = j
    smin = ss(j)
  end if

!     Update position of first element in row of r.

  50 if (j < last) pos1 = pos1 + ncol - j
end do

return
end subroutine drop1


subroutine efroym(first, last, fin, fout, size, ier, lout)

!     Efroymson's stepwise regression from variables in positions FIRST,
!     ..., LAST.  If FIRST > 1, variables in positions prior to this are
!     forced in.  If LAST < ncol, variables in positions after this are
!     forced out.

!     A report is written to unit LOUT if LOUT >= 0.

integer, intent(in)    :: first, last, lout
integer, intent(out)   :: size, ier
real (dp), intent(in)  :: fin, fout

!     Local variables

integer    :: jmax, jmin, i
real (dp)  :: one = 1.0, eps, zero = 0.0, ss(last), smax, base, var, f, smin

!     Check call arguments

ier = 0
if (first >= ncol) ier = 1
if (last <= 1) ier = ier + 2
if (first < 1) ier = ier + 4
if (last > ncol) ier = ier + 8
if (fin < fout .or. fin <= zero) ier = ier + 256
if (nobs <= ncol) ier = ier + 512
if (ier /= 0) then
  size = 0
  return
end if

!     EPS approximates the smallest quantity such that the calculated value of
!     (1 + EPS) is > 1.   It is used to test for a perfect fit (RSS = 0).

eps = epsilon(one)

!     SIZE = number of variables in the current subset

size = first - 1

!     Find the best variable to add next

20 call add1(size+1, last, ss, smax, jmax, ier)
if (nbest > 0) call exadd1(size+1, smax, jmax, ss, last)

!     Calculate 'F-to-enter' value

if (size > 0) then
  base = rss(size)
else
  base = rss(1) + ss(1)
end if
var = (base - smax) / (nobs - size - 1)
if (var < eps*base) then
  ier = -1
  f = zero
else
  f = smax / var
end if
if (lout >= 0) write(lout, 900) vorder(jmax), f
900 format(' best variable to add:  ', i4, '  f-to-enter = ', f10.2)

!     Exit if F < FIN or IER < 0 (perfect fit)

if (f < fin .or. ier < 0) return

!     Add the variable to the subset (in position FIRST).

if (lout >= 0) write(lout, '(50x, "variable added")')
size = size + 1
if (jmax > first) call vmove(jmax, first, ier)
do i = first, min(jmax-1, max_size)
  call report(i, rss(i))
end do

!     See whether a variable entered earlier can be deleted now.

30 if (size <= first) go to 20
call drop1(first+1, size, ss, smin, jmin, ier)
call exdrop1(first+1, size, ss, smin, jmin)
var = rss(size) / (nobs - size)
f = smin / var
if (lout >= 0) write(lout, 910) vorder(jmin), f
910 format(' best variable to drop: ', i4, '  f-to-drop  = ', f10.2)

if (f < fout) then
  if (lout >= 0) write(lout, '(50x, "variable dropped")')
  call vmove(jmin, size, ier)
  if (nbest > 0) then
    do i = jmin, size-1
      call report(i, rss(i))
    end do
  end if
  size = size - 1
  go to 30
end if

if (size >= last) return
go to 20
end subroutine efroym


subroutine exadd1(ivar, smax, jmax, ss, last)

!     Update the NBEST subsets of IVAR variables found from a call
!     to subroutine ADD1.

integer, intent(in)    :: ivar, jmax, last
real (dp), intent(in)  :: smax, ss(:)

!     Local variables

real (dp)  :: zero = 0.0_dp, ssbase, sm, temp, wk(last)
integer    :: i, j, ltemp, jm

if (jmax == 0) return
if (ivar <= 0) return
if (ivar > max_size) return
ltemp = vorder(ivar)
jm = jmax
sm = smax
if (ivar > 1) ssbase = rss(ivar-1)
if (ivar == 1) ssbase = rss(ivar) + ss(1)
wk(ivar:last) = ss(ivar:last)

do i = 1, nbest
  temp = max(ssbase - sm, zero)
  if (temp >= bound(ivar)) exit
  vorder(ivar) = vorder(jm)
  if (jm == ivar) vorder(ivar) = ltemp
  call report(ivar, temp)
  if (i >= nbest) exit
  wk(jm) = zero
  sm = zero
  jm = 0
  do j = ivar, last
    if (wk(j) <= sm) cycle
    jm = j
    sm = wk(j)
  end do
  if (jm == 0) exit
end do

!     Restore VORDER(IVAR)

vorder(ivar) = ltemp

return
end subroutine exadd1


subroutine exdrop1(first, last, ss, smin, jmin)
! Record any new subsets of (LAST-1) variables found from a call to DROP1

integer, intent(in)    :: first, last, jmin
real (dp), intent(in)  :: ss(:), smin

! Local variables
integer    :: list(1:last), i
real (dp)  :: rss_last, ssq

if (jmin == 0 .or. last < 1 .or. last-1 > max_size) return

rss_last = rss(last)
if (rss_last + smin > bound(last-1)) return

list = vorder(1:last)
do i = first, last-1
  vorder(i:last-1) = list(i+1:last)
  ssq = rss_last + ss(i)
  call report(last-1, ssq)
  vorder(i) = list(i)
end do

return
end subroutine exdrop1


subroutine forwrd(first, last, ier)

!     Forward selection from variables in positions FIRST .. LAST.
!     If FIRST > 1, variables in positions prior to this are forced in.
!     If LAST < ncol, variables in positions after this are forced out.
!     On exit, the array VORDER contains the numbers of the variables
!     in the order in which they were added.

integer, intent(in)  :: first, last
integer, intent(out) :: ier

!     Local variables

integer    :: pos, jmax
real (dp)  :: ss(last), smax

!     Check call arguments

ier = 0
if (first >= ncol) ier = 1
if (last <= 1) ier = ier + 2
if (first < 1) ier = ier + 4
if (last > ncol) ier = ier + 8
if (ier /= 0) return

!     For POS = FIRST .. max_size, call ADD1 to find best variable to put
!     into position POS.

do pos = first, max_size
  call add1(pos, last, ss, smax, jmax, ier)
  if (nbest > 0) call exadd1(pos, smax, jmax, ss, last)

!     Move the best variable to position POS.

  if (jmax > pos) call vmove(jmax, pos, ier)
end do

return
end subroutine forwrd


subroutine report(nv, ssq)

!     Update record of the best NBEST subsets of NV variables, if
!     necessary, using SSQ.

integer, intent(in)    :: nv
real (dp), intent(in)  :: ssq

!     Local variables

integer    :: rank, pos1, j, list(nv)
real (dp)  :: under1 = 0.99999_dp, above1 = 1.00001_dp

!     If residual sum of squares (SSQ) for the new subset > the
!     appropriate bound, return.

if(nv > max_size) return
if(ssq >= bound(nv)) return
pos1 = (nv*(nv-1))/2 + 1

!     Find rank of the new subset

do rank = 1, nbest
  if(ssq < ress(nv,rank)*above1) then
    list = vorder(1:nv)
    call shell(list, nv)

!     Check list of variables if ssq is almost equal to ress(nv,rank) -
!     to avoid including the same subset twice.

    if (ssq > ress(nv,rank)*under1) then
      if (same_vars(list, lopt(pos1:,rank), nv)) return
    end if

!     Record the new subset, and move the others down one place.

    do j = nbest-1, rank, -1
      ress(nv,j+1) = ress(nv,j)
      lopt(pos1:pos1+nv-1, j+1) = lopt(pos1:pos1+nv-1, j)
    end do
    ress(nv,rank) = ssq
    lopt(pos1:pos1+nv-1, rank) = list(1:nv)
    bound(nv) = ress(nv,nbest)
    return
  end if
end do

return
end subroutine report



subroutine shell(l, n)

!      Perform a SHELL-sort on integer array L, sorting into increasing order.

!      Latest revision - 5 July 1995

integer, intent(in)     :: n
integer, intent(in out) :: l(:)

!     Local variables
integer   :: start, finish, temp, new, i1, i2, incr, it

incr = n
do
  incr = incr/3
  if (incr == 2*(incr/2)) incr = incr + 1
  do start = 1, incr
    finish = n

!      TEMP contains the element being compared; IT holds its current
!      location.   It is compared with the elements in locations
!      IT+INCR, IT+2.INCR, ... until a larger element is found.   All
!      smaller elements move INCR locations towards the start.   After
!      each time through the sequence, the FINISH is decreased by INCR
!      until FINISH <= INCR.

    20 i1 = start
    temp = l(i1)
    it = i1

!      I2 = location of element new to be compared with TEMP.
!      Test I2 <= FINISH.

    do
      i2 = i1 + incr
      if (i2 > finish) then
        if (i1 > it) l(i1) = temp
        finish = finish - incr
        exit
      end if
      new = l(i2)

!     If TEMP > NEW, move NEW to lower-numbered position.

      if (temp > new) then
        l(i1) = new
        i1 = i2
        cycle
      end if

!     TEMP <= NEW so do not swap.
!     Use NEW as the next TEMP.

      if (i1 > it) l(i1) = temp
      i1 = i2
      temp = new
      it = i1

!     Repeat until FINISH <= INCR.
    end do

    if (finish > incr) go to 20
  end do

!      Repeat until INCR = 1.

  if (incr <= 1) return
end do

return
end subroutine shell



function same_vars(list1, list2, n) result(same)

logical              :: same
integer, intent(in)  :: n, list1(:), list2(:)

same = all(list1(1:n) == list2(1:n))

return
end function same_vars



subroutine seq2(first, last, ier)

! Sequential replacement algorithm applied to the variables in positions
! FIRST, ..., LAST.   2 variables at a time are added or replaced.
! If FIRST > 1, variables in positions prior to this are forced in.
! If LAST < NP, variables in positions after this are left out.

integer, intent(in)  :: first, last
integer, intent(out) :: ier

!     Local variables

integer  :: nv, nsize

!     Check call arguments

ier = 0
if (first >= ncol) ier = 1
if (last <= 1) ier = ier + 2
if (first < 1) ier = ier + 4
if (last > ncol) ier = ier + 8
if (ier /= 0 .or. nbest <= 0) return

nv = min(max_size, last-1)

!     Outer loop; SIZE = current size of subset being considered.

do nsize = first+1, nv
  call replace2(first, last, nsize)
end do

return
end subroutine seq2



subroutine replace2(first, last, nsize)
! Replace 2 variables at a time from those in positions first, ..., nsize
! with 2 from positions nsize, .., last - if they reduce the RSS.

integer, intent(in)  :: first, last, nsize

! Local variables

integer              :: ier, j1, j2, pos1, pos2, best(2), i, iwk(last)
real (dp)            :: smax, rssnew, rssmin, save_rss
real (dp), parameter :: zero = 0.0_dp

10 best(1) = 0
best(2) = 0
rssmin = rss(nsize)

!     Two loops to place all pairs of variables in positions nsize-1 and nsize.
!     POS1 = destination for variable from position nsize.
!     POS2 = destination for variable from position nsize-1.

do pos1 = first, nsize
  do pos2 = pos1, nsize-1
    call add2(nsize-1, last, smax, j1, j2, ier)

    if (j1+j2 > nsize + nsize - 1) then
      rssnew = max(rss(nsize-2) - smax, zero)
      if (rssnew < rssmin) then
        best(1) = vorder(j1)
        best(2) = vorder(j2)
        iwk(1:nsize-2) = vorder(1:nsize-2)
        rssmin = rssnew
      end if
    end if

    call vmove(nsize-1, pos2, ier)
  end do
  call vmove(nsize, pos1, ier)
  do i = pos1, nsize
    call report(i, rss(i))
  end do
end do

!     If any replacement reduces the RSS, make the best one.

if (best(1) + best(2) > 0) then
  iwk(nsize-1) = best(2)
  iwk(nsize) = best(1)
  save_rss = rss(nsize)
  call reordr(iwk, nsize, 1, ier)
  do i = first, nsize
    call report(i, rss(i))
  end do

!    The calculated value of rssmin above is only a rough approximation to
!    the real residual sum of squares, thiugh usually good enough.
!    The new value of rss(nsize) is more accurate.   It is used below
!    to avoid cycling when several subsets give the same RSS.

  if (rss(nsize) < save_rss) go to 10
end if

return
end subroutine replace2



subroutine seqrep(first, last, ier)

!     Sequential replacement algorithm applied to the variables in
!     positions FIRST, ..., LAST.
!     If FIRST > 1, variables in positions prior to this are forced in.
!     If LAST < ncol, variables in positions after this are forced out.

integer, intent(in)  :: first, last
integer, intent(out) :: ier

!     Local variables

integer    :: nv, size, start, best, from, i, jmax, count, j
real (dp)  :: zero = 0.0_dp, ssred, ss(last), smax

!     Check call arguments

ier = 0
if (first >= ncol) ier = 1
if (last <= 1) ier = ier + 2
if (first < 1) ier = ier + 4
if (last > ncol) ier = ier + 8
if (ier /= 0 .or. nbest <= 0) return

nv = min(max_size, last-1)

!     Outer loop; SIZE = current size of subset being considered.

do size = first, nv
  count = 0
  start = first
  10 ssred = zero
  best = 0
  from = 0

!     Find the best variable from those in positions SIZE+1, ..., LAST
!     to replace the one in position SIZE.   Then rotate variables in
!     positions START, ..., SIZE.

  do i = start, size
    call add1(size, last, ss, smax, jmax, ier)
    if (jmax > size) then
      call exadd1(size, smax, jmax, ss, last)
      if (smax > ssred) then
        ssred = smax
        best = jmax
        if (i < size) then
          from = size + start - i - 1
        else
          from = size
        end if
      end if
    end if
    if (i < size) call vmove(size, start, ier)
    do j = start, size-1
      call report(j, rss(j))
    end do
  end do ! i = start, size

!     If any replacement reduces the RSS, make the best one.
!     Move variable from position FROM to SIZE.
!     Move variable from position BEST to FIRST.

  if (best > size) then
    if (from < size) call vmove(from, size, ier)
    call vmove(best, first, ier)
    do j = first, best-1
      call report(j, rss(j))
    end do
    count = 0
    start = first + 1
  else
    count = count + 1
  end if

!     Repeat until COUNT = SIZE - START + 1

  if (count <= size - start) go to 10
end do

return
end subroutine seqrep


subroutine xhaust(first, last, ier)

!     Exhaustive search algorithm, using leaps and bounds, applied to
!     the variables in positions FIRST, ..., LAST.
!     If FIRST > 1, variables in positions prior to this are forced in.
!     If LAST < ncol, variables in positions after this are forced out.

integer, intent(in)  :: first, last
integer, intent(out) :: ier

!     Local variables

integer    :: row, i, jmax, ipt, newpos, iwk(max_size)
real (dp)  :: ss(last), smax, temp

!     Check call arguments

ier = 0
if (first >= ncol) ier = 1
if (last <= 1) ier = ier + 2
if (first < 1) ier = ier + 4
if (last > ncol) ier = ier + 8
if (ier /= 0 .or. nbest <= 0) return

!     Record subsets contained in the initial ordering, including check
!     for variables which are linearly related to earlier variables.
!     This should be redundant if the user has first called SING and
!     init_subsets.

do row = first, max_size
  if (d(row) <= tol(row)) then
    ier = -999
    return
  end if
  call report(row, rss(row))
end do

!     IWK(I) contains the upper limit for the I-th simulated DO-loop for
!     I = FIRST, ..., max_size-1.
!     IPT points to the current DO loop.

iwk(first:max_size) = last

!     Innermost loop.
!     Find best possible variable for position max_size from those in
!     positions max_size, .., IWK(max_size).

30 call add1(max_size, iwk(max_size), ss, smax, jmax, ier)
call exadd1(max_size, smax, jmax, ss, iwk(max_size))

!     Move to next lower numbered loop which has not been exhausted.

ipt = max_size - 1
40 if (ipt >= iwk(ipt)) then
  ipt = ipt - 1
  if (ipt >= first) go to 40
  return
end if

!     Lower variable from position IPT to position IWK(IPT).
!     Record any good new subsets found by the move.

newpos = iwk(ipt)
call vmove(ipt, newpos, ier)
do i = ipt, min(max_size, newpos-1)
  call report(i, rss(i))
end do

!     Reset all ends of loops for I >= IPT.

iwk(ipt:max_size) = newpos - 1

!     If residual sum of squares for all variables above position NEWPOS
!     is greater than BOUND(I), no better subsets of size I can be found
!     inside the current loop.

temp = rss(newpos-1)
do i = ipt, max_size
  if (temp > bound(i)) go to 80
end do
if (iwk(max_size) > max_size) go to 30
ipt = max_size - 1
go to 40

80 ipt = i - 1
if (ipt < first) return
go to 40

end subroutine xhaust



subroutine random_pick(first, last, npick)
! Pick npick variables at random from those in positions first, ..., last
! and move them to occupy positions starting from first.

integer, intent(in)  :: first, last, npick

! Local variables

integer  :: first2, i, ilist(1:last), j, k, navail
real     :: r

navail = last + 1 - first
if (npick >= navail .or. npick <= 0) return
do i = first, last
  ilist(i) = vorder(i)
end do

first2 = first
do i = 1, npick
  call random_number(r)
  k = first2 + r*navail
  if (k > first2) then
    j = ilist(first2)
    ilist(first2) = ilist(k)
    ilist(k) = j
  end if
  first2 = first2 + 1
  navail = navail - 1
end do

call reordr(ilist(first:), npick, first, i)

return
end subroutine random_pick
!

!


subroutine cross_validation(unit_data, line1, ypos, fit_const, nvar, first,  &
                            last, search_method, criterion, nrepl, seed,   &
                            unit_rpt, msep, ier)
! Cross-validation routine excluding 10% at a time
! Search_method
! 1  Efroymson stepwise (F-to-enter = F-to-delete = 4.0)
! 2  Sequential replacement
! 3  Two-at-a-time sequential replacement
! 4  Best subsets (exhaustive search with branch-and-bound)
! Criterion (set to 5 for Efroymson search)
! 1  AIC
! 2  BIC
! 3  Mallows Cp
! 4  Hannan-Quinn
! 5  F-ratio = 4.0
! Other arguments:
! unit_data Unit number from which to read the data
! line1     Number of the first line of data in the input file (in case file
!           contains variable names & other header information)
! ypos      Position of the dependent variable in each line of data
! fit_const .true. if a constant is to be included in the model
! nvar      Number of predictor variables, excluding any constant
! first     Position of the first variable available for selection.   Any
!           variables in earlier positions are to be forced into all subsets.
! last      Position of the last variable available for selection.   Variables
!           to be forced out of subsets (if any) should be after position last.
! nrepl     Number of complete replications of 10 subsets of 10% of the data
! seed()    Array of random number seeds (dimension depends upon compiler)
! unit_rpt  Unit number for output of report
! msep      Average mean squared error of the predictions
! ier       Error indicator
!           = 0 if no error detected
! This version - 15 August 2002
integer, intent(in)      :: unit_data, line1, ypos, nvar, first, last,  &
                            search_method, nrepl, seed(:), unit_rpt
integer, intent(in out)  :: criterion
logical, intent(in)      :: fit_const
integer, intent(out)     :: ier
real (dp), intent(out)   :: msep
! Local variables
integer              :: replicate, i, nobs_full, percentile, i1, i2, case,   &
                        iostatus, num_seeds, nvar_max, nsize, ipos, lout = 6, j
integer, allocatable :: order(:), seeds(:), list(:), vorder_cpy(:),  &
                        init_vorder(:)
real (dp)            :: sumsq_ls, total_ls, zero = 0.0_dp, minus1 = -1.0_dp, y,   &
                  weight, fin = 4.0_dp, fout = 4.0_dp, variance, one = 1.0_dp, &
                        fit_ls, penalty, crit_val, min_crit_val
real (dp), allocatable :: x(:), beta_ls(:)
if (search_method == 1) criterion = 5
allocate( order(nobs), x(0:nvar), init_vorder(ncol) )
nobs_full = nobs
msep = zero
init_vorder = vorder
allocate( list(1:nvar) )
do i = 1, nvar
  list(i) = i
end do
call random_seed(size=num_seeds)
allocate( seeds(num_seeds) )
seeds = seed(1:num_seeds)
call random_seed(put=seeds)
! Return the QR factorization back to the order of variables in the data set.
do i = 1, nvar
  call reordr(list, i, 2, ier)
end do
do replicate = 1, nrepl
! Choose a random permutation of the integers 1, 2, ..., nobs
  do i = 1, nobs_full
    order(i) = i
  end do
  call ranord(order, nobs_full)
  total_ls = zero
! ------------------- Cycle through subsets of the data ----------------
! Delete 10% of the observations using array `order'
  do percentile = 1, 10
    if (percentile == 1) then
      i1 = 1
      i2 = 0.1 * real(nobs_full) + 0.5
    else
      i1 = i2 + 1
      if (percentile == 10) then
      i2 = nobs_full
      else
        i2 = 0.1 * real(nobs_full) * percentile + 0.5
      end if
    end if
    write (*, "(a, i4, a, i4) ") " replicate: ", replicate,  &
                              "  percentile: ", percentile
    sumsq_ls = zero
    rewind(unit_data)
    if (line1 > 1) then
      do i = 1, line1-1
        read(unit_data,*)
      end do
    end if
    case = 1
    weight = minus1
    x(0) = one
    do
      if (ypos > nvar) then
        read(unit_data,*, iostat=iostatus) x(1:nvar), y
      else if (ypos == 1) then
        read(unit_data,*, iostat=iostatus) y, x(1:nvar)
      else
        read(unit_data,*, iostat=iostatus) x(1:ypos-1), y, x(ypos:nvar)
      end if
      if (iostatus > 0) cycle                   ! Error in data
      if (iostatus < 0) exit                    ! End of file
      if(any(case == order(i1:i2))) then        ! Delete case if in this 10%
        call includ(weight, x, y)
      end if
      case = case + 1
    end do
! N.B. Subroutine INCLUD increases nobs even when weights are negative.
! Calculate correct value for the present nobs.
    nobs = nobs_full - (i2 + 1 - i1)
! Find subsets which fit well
    if (fit_const) then
      nvar_max = max_size - 1
    else
      nvar_max = max_size
    end if
    call init_subsets(nvar_max, fit_const)
                             ! Re-order the QR-factorization if variables are
                             ! to be forced in or out.
    if (first > 1) then
      call reordr(init_vorder, first-1, 1, ier)
    end if
    if (last < ncol) then
      call reordr(init_vorder, last, 1, ier)
    end if
    select case(search_method)
      case(1)
        call efroym(first, last, fin, fout, nsize, ier, lout)
      case(2)
        call seqrep(first, last, ier)
      case(3)
        call seq2(first, last, ier)
      case(4)
        call seq2(first, last, ier)
        call xhaust(first, last, ier)
    end select
! Pick best subset
    min_crit_val = huge(one)
    if (search_method > 1) then
      if (criterion == 3) variance = sserr / (nobs - ncol)
      nsize = 0
      do i = first, max_size
        if (criterion /= 3) variance = ress(i,1) / (nobs - i)
        call calc_penalty(i, first, variance, criterion, penalty)
        crit_val = ress(i,1) + penalty
        if (nsize == 0 .or. crit_val < min_crit_val) then
          min_crit_val = crit_val
          nsize = i
        end if
      end do
    end if
    ipos = ((nsize-1)*nsize)/2 + 1
    call reordr(lopt(ipos:,1), nsize, 1, ier)
! Estimate the regression coefficients using the LS-projections
    allocate( beta_ls(1:nsize), vorder_cpy(1:nsize) )
    vorder_cpy = vorder(1:nsize)
    call shell(vorder_cpy, nsize)                 ! Shell sort from find_sub
    write (unit_rpt, 970) percentile, vorder_cpy(1:nsize)
    970 format("percentile no.", i3, "   selected variables:"/(" ", 15i5))
    call regcf(beta_ls, nsize, ier)
!    WRITE (unit_rpt, 980) beta_LS
!    980 FORMAT('LS regression coefficients:'/(' ', 6g13.5))
    vorder_cpy = vorder(1:nsize)
! Return the order of variables to the original order, as in the data set
    do i = 1, nvar
      call reordr(list, i, 2, ier)
    end do
! Estimate the 10% omitted, and re-instate the deleted cases
    rewind(unit_data)
    if (line1 > 1) then
      do i = 1, line1-1
        read(unit_data,*)
      end do
    end if
    case = 1
    weight = one
    x(0) = one
    do
      if (ypos > nvar) then
        read(unit_data,*, iostat=iostatus) x(1:nvar), y
      else if (ypos == 1) then
        read(unit_data,*, iostat=iostatus) y, x(1:nvar)
      else
        read(unit_data,*, iostat=iostatus) x(1:ypos-1), y, x(ypos:nvar)
      end if
      if (iostatus > 0) cycle                   ! Error in data
      if (iostatus < 0) exit                    ! End of file
      if(any(case == order(i1:i2))) then        ! Restore case if in this 10%
        fit_ls = zero
        do i = 1, nsize
          j = vorder_cpy(i)
          fit_ls = fit_ls + beta_ls(i) * x(j)
        end do
        sumsq_ls = sumsq_ls + (y - fit_ls)**2
        call includ(weight, x, y)                ! INCLUD destroys x
      end if
      case = case + 1
    end do
    write (unit_rpt,1000) sumsq_ls
    1000 format("sums of sq. (ls) = ", g12.4)
    total_ls = total_ls + sumsq_ls
!   CALL print_QR
    deallocate( beta_ls, vorder_cpy )
  end do             ! percentile = 1, 10
  write (unit_rpt, 1030) total_ls
  write (*, 1030) total_ls
  1030 format(/" total sum of squares (ls) = ", g13.5)
  write (unit_rpt, "(/) ")
  msep = msep + total_ls
end do             ! replicate = 1, nrepl
msep = msep / (nrepl * nobs_full)
write (*, 900) msep
write (unit_rpt, 900) msep
900 format(" overall mean squared error of prediction = ", g13.5)
write (*, 910) sqrt(msep)
write (unit_rpt, 910) sqrt(msep)
910 format(" rms (prediction error) = ", g13.5)
deallocate( order, x, list, seeds )
stop
end subroutine cross_validation
!
end module find_subsets
