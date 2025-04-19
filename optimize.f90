module optimize_mod
implicit none
private
public :: dp,random_search
integer, parameter :: dp=kind(1.0d0)
contains
subroutine foo(xx,yy)
real(kind=dp), intent(in)  :: xx(:)
real(kind=dp), intent(out) :: yy
if (size(xx) >= 2) yy = (xx(1)-3)**2 + (xx(2)+5)**2
end subroutine foo
!
subroutine random_search(xmin,xmax,nran,iopt,xbest,ybest)
real(kind=dp), intent(in)  :: xmin(:),xmax(:)
integer      , intent(in)  :: nran
integer      , intent(in)  :: iopt ! if 1, seach in same neigborhood defined by xmin(:), xmax(:) for all guesses
real(kind=dp), intent(out) :: xbest(:),ybest
real(kind=dp)              :: fval,rng(size(xmin)),xran(size(xmin)),xx(size(xmin))
integer                    :: i,ndim
logical                    :: minim
ndim = size(xmin)
if (size(xmax) /= ndim .or. size(xbest) /= ndim) then
   write (*,*) "in random_search, sizes of xmin, xmax, xbest are ",ndim,size(xmax),size(xbest)," must be equal, STOPPING"
   stop
end if
rng = xmax - xmin
xbest = (xmax + xmin)/2
do i=1,nran
   call random_number(xran)
   if (iopt == 1) then
      xx = xmin + xran*rng
   else
      xx = xbest + (xran-0.5_dp)*rng
   end if
   call foo(xx,fval)
   if (i > 1) then
      minim = (fval < ybest)
   else
      minim = .true.
   end if
   if (minim) then
      ybest = fval
      xbest = xx
   end if
end do
end subroutine random_search
end module optimize_mod
