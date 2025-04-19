module cv_reg_mod
use kind_mod    , only: dp
use linear_solve, only: hat_matrix
implicit none
private
public :: cv_reg_resid
contains
function cv_reg_resid(aa,bb) result(resid)
! return residuals from leave-one-out cross validated (LOOCV) linear regression
real(kind=dp)     , intent(in)             :: aa(:,:)      ! (nrow,ncol) input: LHS matrix
real(kind=dp)     , intent(in)             :: bb(:)        ! (nrow)      input: RHS vector
real(kind=dp)                              :: resid(size(bb))
real(kind=dp)                              :: hat(size(bb),size(bb)),bpred(size(bb)),resid_full(size(bb)),resid_abs,denom
integer                                    :: i
hat = hat_matrix(aa)
bpred = matmul(hat,bb)
resid_full = bb - bpred
do i=1,size(bb)
   denom = 1-hat(i,i)
!   print*,"i, denom=",i,denom
   if (denom > 0) then
      resid_abs = abs(resid_full(i))/denom
      resid(i) = merge(1,-1,resid_full(i) > 0.0_dp) * resid_abs
   else
      resid(i) = resid_full(i)
   end if
end do
end function cv_reg_resid
end module cv_reg_mod