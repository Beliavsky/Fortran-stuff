module change_point_mod
implicit none
private
public :: dp,change_points_2,change_points_2_slow
integer, parameter :: dp = kind(1.0d0)
contains
subroutine change_points_2_slow(x,xmean,i1,i2)
! find the two positions i1 and i2 that best partition x(:) into segments of minimum variance
real(kind=dp), intent(in)  :: x(:)
real(kind=dp), intent(out) :: xmean(:)
integer      , intent(out) :: i1,i2
integer                    :: j1,j2,n
real(kind=dp)              :: sumsq,sumsq_min,xmean_try(size(x))
n = size(x)
if (n < 3) stop "need size(x) > 2"
if (size(xmean) /= size(x)) stop "need size(x) == size(xmean)"
sumsq_min = huge(x)
do j1=2,n-1
   call set_mean(x,1,j1-1,xmean_try)
   do j2=j1+1,n
      call set_mean(x,j1,j2-1,xmean_try)
      call set_mean(x,j2,n,xmean_try)
      sumsq = sum((x-xmean_try)**2)
      if (sumsq < sumsq_min) then
         sumsq_min = sumsq
         i1 = j1
         i2 = j2
         xmean = xmean_try
      end if  
   end do
end do
end subroutine change_points_2_slow
!
subroutine change_points_2(x,xmean,i1,i2)
! find the two positions i1 and i2 that best partition x(:) into segments of minimum variance
! precomputes the cumulative sum of x(:) for speed
real(kind=dp), intent(in)  :: x(:)
real(kind=dp), intent(out) :: xmean(:)
integer      , intent(out) :: i1,i2
integer                    :: i,j1,j2,n
real(kind=dp)              :: sumsq,sumsq_min,xmean_try(size(x)),cumul_sum(size(x))
n = size(x)
if (n < 3) stop "need size(x) > 2"
if (size(xmean) /= size(x)) stop "need size(x) == size(xmean)"
cumul_sum(1) = x(1)
do i=2,n
   cumul_sum(i) = cumul_sum(i-1) + x(i)
end do
sumsq_min = huge(x)
do j1=2,n-1
   call set_mean(x,1,j1-1,xmean_try,cumul_sum)
   do j2=j1+1,n
      call set_mean(x,j1,j2-1,xmean_try,cumul_sum)
      call set_mean(x,j2,n,xmean_try,cumul_sum)
      sumsq = sum((x-xmean_try)**2)
      if (sumsq < sumsq_min) then
         sumsq_min = sumsq
         i1 = j1
         i2 = j2
         xmean = xmean_try
      end if  
   end do
end do
end subroutine change_points_2
!
subroutine set_mean(x,j1,j2,xmean,cumul_sum)
! set xmean(j1:j2) to mean(x(j1:j2))
real(kind=dp), intent(in)    :: x(:)
integer      , intent(in)    :: j1,j2
real(kind=dp), intent(inout) :: xmean(:) 
real(kind=dp), intent(in), optional :: cumul_sum(:)
integer                      :: n,denom
n = size(x)
if (n > 0 .and. j1 > 0 .and. j2 >= j1 .and. j2 <= n) then
   if (j1 == j2) then
      xmean(j1) = x(j1)
      return
   end if
   denom = j2 - j1 + 1
   if (present(cumul_sum)) then
      if (j1 > 1) then
         xmean(j1:j2) = (cumul_sum(j2)-cumul_sum(j1-1))/denom
      else
         xmean(j1:j2) = cumul_sum(j2)/denom
      end if
   else
      xmean(j1:j2) = sum(x(j1:j2))/denom
   end if
end if
end subroutine set_mean
end module change_point_mod