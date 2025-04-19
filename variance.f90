module variance_mod
use kind_mod, only: dp
implicit none
real, parameter :: bad_real = -999.0_dp
public :: variance,rolling_variance,sd,rolling_sd
interface variance
   module procedure variance_vec,variance_matrix
end interface variance
interface sd
   module procedure sd_vec,sd_matrix
end interface sd
interface rolling_variance
   module procedure rolling_variance_vec,rolling_variance_matrix
end interface rolling_variance
interface rolling_sd
   module procedure rolling_sd_vec,rolling_sd_matrix
end interface rolling_sd
contains
pure function variance_vec(xx,xxmean) result(xvar)
! compute the variance of xx(:)
real(kind=dp), intent(in)           :: xx(:)   ! data for which variance computed
real(kind=dp), intent(in), optional :: xxmean  ! mean to use in oomputing variance
real(kind=dp)                       :: xvar    ! result
integer                             :: n
real(kind=dp)                       :: xxmean_
n   = size(xx)
if (size(xx) < 2) then
   xvar = bad_real
   return
end if
if (present(xxmean)) then
   xxmean_ = xxmean
else
   xxmean_ = sum(xx)/n
end if
xvar = sum((xx-xxmean_)**2)/merge(n,n-1,present(xxmean))
end function variance_vec
!
pure function sd_vec(xx,xxmean) result(xsd)
! compute the standard deviation of xx(:)
real(kind=dp), intent(in)           :: xx(:)   ! data for which standard deviation computed
real(kind=dp), intent(in), optional :: xxmean  ! mean to use in oomputing variance
real(kind=dp)                       :: xsd     ! result
xsd = sqrt(variance_vec(xx,xxmean))
end function sd_vec
!
pure function variance_matrix(xx,xxmean) result(xvar)
! compute the variance of each column of xx(:,:)
real(kind=dp), intent(in)           :: xx(:,:)          ! (n,ncol) -- data for which variance computed
real(kind=dp), intent(in), optional :: xxmean(:)        ! (ncol)   -- means of columns of xx(:,:)
real(kind=dp)                       :: xvar(size(xx,2)) ! (ncol)
integer                             :: icol,ncol
ncol = size(xx,2)
if (present(xxmean)) then
   forall (icol=1:ncol) xvar(icol) = variance(xx(:,icol),xxmean(icol))
else
   forall (icol=1:ncol) xvar(icol) = variance(xx(:,icol))
end if
end function variance_matrix
!
pure function sd_matrix(xx,xxmean) result(xsd)
! compute the standard deviation of each column of xx(:,:)
real(kind=dp), intent(in)           :: xx(:,:)          ! (n,ncol) -- data for which standard deviation computed
real(kind=dp), intent(in), optional :: xxmean(:)        ! (ncol)   -- means of columns of xx(:,:)
real(kind=dp)                       :: xsd(size(xx,2))  ! (ncol)
integer                             :: icol,ncol
ncol = size(xx,2)
if (present(xxmean)) then
   forall (icol=1:ncol) xsd(icol) = sd(xx(:,icol),xxmean(icol))
else
   forall (icol=1:ncol) xsd(icol) = sd(xx(:,icol))
end if
end function sd_matrix
!
pure function rolling_variance_vec(nlook,xx,xxmean) result(xvar)
! compute the rolling covariance of xx(:)
integer      , intent(in)           :: nlook  ! # of trailing observations to use to compute variance
real(kind=dp), intent(in)           :: xx(:)  ! data for which variance computed
real(kind=dp), intent(in), optional :: xxmean ! mean of data
real(kind=dp)                       :: xvar(size(xx))
integer                             :: i,imin
do i=1,size(xx)
   imin    = max(1,i-nlook+1)
   xvar(i) = variance(xx(imin:i),xxmean)
end do
end function rolling_variance_vec
!
pure function rolling_sd_vec(nlook,xx,xxmean) result(xvar)
! compute the rolling standard deviation of xx(:)
Integer      , intent(in)           :: nlook  ! # of trailing observations to use to compute variance
real(kind=dp), intent(in)           :: xx(:)  ! data for which variance computed
real(kind=dp), intent(in), optional :: xxmean ! mean of data
real(kind=dp)                       :: xvar(size(xx))
integer                             :: i,imin
do i=1,size(xx)
   imin    = max(1,i-nlook+1)
   xvar(i) = sd(xx(imin:i),xxmean)
end do
end function rolling_sd_vec
!
pure function rolling_variance_matrix(nlook,xx,xxmean) result(xvar)
! compute the rolling covariance of xx(:)
integer      , intent(in)           :: nlook     ! # of trailing observations to use to compute variance
real(kind=dp), intent(in)           :: xx(:,:)   ! data for which variance computed
real(kind=dp), intent(in), optional :: xxmean(:) ! mean of data
real(kind=dp)                       :: xvar(size(xx,1),size(xx,2))
integer                             :: i,imin
do i=1,size(xx,1)
   imin      = max(1,i-nlook+1)
   xvar(i,:) = variance(xx(imin:i,:),xxmean)
end do
end function rolling_variance_matrix
!
pure function rolling_sd_matrix(nlook,xx,xxmean) result(xvar)
! compute the rolling standard deviation of xx(:,:)
integer      , intent(in)           :: nlook     ! # of trailing observations to use to compute standard deviation
real(kind=dp), intent(in)           :: xx(:,:)   ! data for which standard deviation computed
real(kind=dp), intent(in), optional :: xxmean(:) ! mean of data
real(kind=dp)                       :: xvar(size(xx,1),size(xx,2))
integer                             :: i,imin
do i=1,size(xx,1)
   imin      = max(1,i-nlook+1)
   xvar(i,:) = sd(xx(imin:i,:),xxmean)
end do
end function rolling_sd_matrix
!
end module variance_mod
