module covar_mod
use kind_mod, only: dp
implicit none
private
public :: covar_mat
contains
function covar_mat(xmat) result(cov)
! compute the covariance of the columns of xmat(:,:)
real(kind=dp), intent(in)    :: xmat(:,:) ! data for which covariance matrix is calculated
real(kind=dp)                :: cov(size(xmat,2),size(xmat,2))
integer                      :: n,nvar
integer                      :: i,j
real(kind=dp)                :: xmean(size(xmat,2))
real(kind=dp)                :: ymat(size(xmat,1),size(xmat,2))
n    = size(xmat,1)
nvar = size(xmat,2)
cov  = 0.0_dp
if (n < 2) then
   write (*,*) "in covar_mat, size(xmat,1)=",n," should be > 1"
   return
end if
ymat = xmat
do i=1,nvar
   ymat(:,i) = xmat(:,i) - sum(xmat(:,i))/n
end do
do i=1,nvar
   do j=1,i
      cov(i,j) = sum(ymat(:,i)*ymat(:,j))/(n-1)
   end do
end do
do i=1,nvar
   do j=i+1,nvar
      cov(i,j) = cov(j,i) ! copy above diagonal elements to below diagonal
   end do
end do
end function covar_mat
end module covar_mod
