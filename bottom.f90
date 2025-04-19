module bottom_mod
implicit none
integer, parameter :: dp = kind(1.0d0)
private
public :: dp, first_lower, first_lower_look, first_lower_many_thresh
contains
function first_lower(max_steps, thresh, drift, sd) result(ifirst)
! return the first time that the cumulative sum of normal variates
! with given mean (drift) and standard deviation (sd) is below thresh
integer      , intent(in) :: max_steps
real(kind=dp), intent(in) :: thresh
real(kind=dp), intent(in) :: drift
real(kind=dp), intent(in) :: sd
integer                   :: ifirst
integer                   :: i
real(kind=dp)             :: x
x = 0.0_dp
do i=1,max_steps
   x = x + drift + sd*random_normal()
   if (x < thresh) then
      ifirst = i
      return
   end if
end do
ifirst = 0 ! threshold never breached
end function first_lower
!
function first_lower_many_thresh(max_steps, thresh, drift, sd) result(ifirst)
! return the first time that the cumulative sum of normal variates
! with given mean (drift) and standard deviation (sd) is below thresh(:)
integer      , intent(in) :: max_steps
real(kind=dp), intent(in) :: thresh(:)
real(kind=dp), intent(in) :: drift
real(kind=dp), intent(in) :: sd
integer                   :: ifirst(size(thresh))
integer                   :: i,ithresh
real(kind=dp)             :: x, xchange
x = 0.0_dp
ifirst = 0
do i=1,max_steps
   xchange = drift + sd*random_normal()
   x = x + xchange
   if (xchange < 0.0_dp) then
      do ithresh=1,size(thresh)
         if (ifirst(ithresh) == 0 .and. x < thresh(ithresh)) ifirst(ithresh) = i
      end do
   end if
end do
end function first_lower_many_thresh
!
function first_lower_look(max_steps, nlook, thresh, drift, sd) result(ifirst)
! return the first times that the cumulative sum of normal variates
! with given mean (drift) and standard deviation (sd) is below thresh,
! given times between monitoring nlook(:)
integer      , intent(in) :: max_steps
integer      , intent(in) :: nlook(:) ! # of periods between observations
real(kind=dp), intent(in) :: thresh
real(kind=dp), intent(in) :: drift
real(kind=dp), intent(in) :: sd
integer                   :: ifirst(size(nlook))
integer                   :: istep,ilook
real(kind=dp)             :: x
if (any(nlook < 1)) error stop "in first_lower_look need all(nlook > 0)"
x = 0.0_dp
ifirst = 0
do istep=1,max_steps
   x = x + drift + sd*random_normal()
   if (x < thresh) then
      do ilook=1,size(nlook)
         if (ifirst(ilook) == 0 .and. mod(istep,nlook(ilook)) == 0) ifirst(ilook) = istep
      end do
      if (all(ifirst /= 0)) return
   end if
end do
end function first_lower_look
!
function random_normal() result(fn_val)
!   Generate a random normal deviate using the polar method. (Alan Miller)
!   Reference: Marsaglia,G. & Bray,T.A. 'A convenient method for generating
!              normal variables',Siam Rev.,vol.6,260-264,1964.
real(kind=dp)  :: fn_val
! Local variables
real(kind=dp)           :: u,sum
real(kind=dp),save      :: v,sln
logical,save   :: second = .false.
real(kind=dp),parameter :: one = 1.0d0,vsmall = tiny(one)
if (second) then
! If second,use the second random number generated on last call
  second = .false.
  fn_val = v*sln
else
! First call; generate a pair of random normals
  second = .true.
  do
    call random_number(u)
    call random_number(v)
    u = scale(u,1) - one
    v = scale(v,1) - one
    sum = u*u + v*v + vsmall         ! vsmall added to prevent log(zero) / zero
    if(sum < one) exit
  end do
  sln = sqrt(- scale(log(sum),1) / sum)
  fn_val = u*sln
end if
end function random_normal
end module bottom_mod
