module reg_sub_mod
! Alan Miller's codes
use kind_mod, only: dp
use lsq     , only: rss,reordr
private
public :: calc_penalty,ranord,crossvalidation,update,f1max
contains
subroutine calc_penalty(size1, size2, variance, penalty_num, penalty_val)
! Calculate value of penalty for size of subset.
! Currently the penalties available are:
! Number Name
!    1   AIC
!    2   BIC
!    3   Mallows' Cp
!    4   Hannan-Quinn
!    5   Efroymson (F-to-delete = F-to-add = 4.0)
integer , intent(in)   :: size1, size2, penalty_num
real(dp), intent(in)   :: variance
real(dp), intent(out)  :: penalty_val
! Local variables
real(dp)  :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp, four = 4.0_dp
if (size1 == size2) then
  penalty_val = zero
  return
end if
if (penalty_num < 1 .or. penalty_num > 5) then
  penalty_val = zero
  return
end if
select case(penalty_num)
  case(1)
    penalty_val = rss(size1) * (one - exp(two * (size2 - size1) / nobs))
  case(2)
    penalty_val = rss(size1) *    &
                  (one - exp((size2 - size1) * log(real(nobs)) / nobs))
  case(3)
    penalty_val = two * variance * (size1 - size2)
  case(4)
    penalty_val = rss(size1) *    &
                  (one - exp(two * (size2 - size1) * log(log(real(nobs))) / nobs))
  case(5)
    penalty_val = four * variance * (size1 - size2)
end select
end subroutine calc_penalty
!
subroutine ranord(order, n)
! Generate a random ordering of the integers 1 ... n.
integer, intent(in)  :: n
integer, intent(out) :: order(:)
! Local variables
real                 :: wk(n)
integer              :: i
do i = 1, n
  order(i) = i
end do
call random_number(wk)
call sqsort(wk, n, order) ! in out: wk,order
end subroutine ranord
!
subroutine sqsort(a, n, t)
!   NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTH'S PASCAL
!   BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.
!   SINGLE PRECISION, ALSO CHANGES THE ORDER OF THE ASSOCIATED ARRAY T.
integer, intent(in)     :: n
integer, intent(in out) :: t(:)
real, intent(in out)    :: a(:)
! Local variables
integer :: i, j, k, l, r, s, stackl(15), stackr(15), ww
real    :: w, x
s = 1
stackl(1) = 1
! KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.
stackr(1) = n
10 l = stackl(s)
r = stackr(s)
! KEEP SPLITTING A(L), ... , A(R) UNTIL L>=R.
s = s - 1
20 i = l
j = r
k = (l + r)/2
! REPEAT UNTIL I > J.
x = a(k)
30 if(a(i) >= x) go to 40
i = i + 1
go to 30
40 if(x >= a(j)) go to 50
j = j - 1
go to 40
50 if(i > j) go to 60
w = a(i)
ww = t(i)
a(i) = a(j)
t(i) = t(j)
a(j) = w
t(j) = ww
i = i + 1
j = j - 1
if(i <= j) go to 30
60 if(j - l < r - i) go to 75
if(l >= j) go to 65
s = s + 1
stackl(s) = l
stackr(s) = j
65 l = i
go to 90
75 if(i >= r) go to 80
s = s + 1
stackl(s) = i
stackr(s) = r
80 r = j
90 if(l < r) go to 20
if(s /= 0) go to 10
return
end subroutine sqsort
!
subroutine update(x, n, xmin, xmax, xmean, sxx)
! Update statistics for x
! n = observation number ( = 1 for the first call)
! (sxx is the sum of squares of deviations from the mean)
real(kind=dp), intent(in)     :: x
integer      , intent(in)     :: n
real(kind=dp), intent(in out) :: xmin, xmax, xmean, sxx
! Local variables
real(kind=dp)  :: dev
if (n == 1) then
  xmin = x
  xmax = x
  xmean = x
  sxx = 0.0d0
  return
end if
if (x < xmin) then
  xmin = x
else if (x > xmax) then
  xmax = x
end if
dev = x - xmean
xmean = xmean + dev/n
sxx = sxx + dev * (x - xmean)
end subroutine update
!

!
elemental subroutine f1max(ndf, nf, f1, f5, f10, ier)
! Calculates approximations to the 1%, 5% & 10% points of the distribution
! of the maximum of NF values of an F-ratio with 1 d.f. for the numerator
! and NDF d.f. for the denominator.
! An approximation is used to the values given in table 2 of:
!     Gilmour, S.G. (1996) `The interpretation of Mallows's Cp-statistic',
!     The Statistician, vol.45, pp.49-56
implicit none
integer      , intent(in)  :: ndf, nf
real(kind=dp), intent(out) :: f1, f5, f10
integer      , intent(out) :: ier
! Local variables
real(kind=dp), parameter ::  a1(6) =  [ 1.67649, 6.94330,  1.22627, 0.25319,  0.06136, -2.41097], &
                             a5(6) =  [ 1.28152, 4.93268, -0.29583, 0.28518, -0.23566, -1.60581], &
                            a10(6) =  [ 1.06642, 3.96276, -0.62483, 0.30228, -0.52843, -1.04499], &
                            one = 1.0_dp
real(kind=dp) :: log_nf
if (ndf < 4 .or. nf < 1) then
  ier = 1
  return
end if
ier = 0
log_nf = log(real(nf,kind=dp))
f1  = one + exp(a1(1) + (a1(3)/ndf + a1(2))/ndf + a1(4)*log_nf + (a1(6)/ndf + a1(5))/nf)
f5  = one + exp(a5(1) + (a5(3)/ndf + a5(2))/ndf + a5(4)*log_nf + (a5(6)/ndf + a5(5))/nf)
f10 = one + exp(a10(1) + (a10(3)/ndf + a10(2))/ndf + a10(4)*log_nf + (a10(6)/ndf + a10(5))/nf)
end subroutine f1max
end module reg_sub_mod