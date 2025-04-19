! extracted from xfind_sub.f90
module reg_mod
use lsq, only: dp
implicit none
contains
subroutine calc_residuals(unit_data, unit_resid, fname_dat, nvar, vname,   &
                          ypos, nreq, beta, list, stdev, autoc1)
integer, intent(in)              :: unit_data, unit_resid, nvar, ypos, nreq,  &
                                    list(:)
character (len = 40), intent(in) :: fname_dat
character (len = 8), intent(in)  :: vname(0:)
real (dp), intent(in)            :: beta(:), stdev
real (dp), intent(out)           :: autoc1
! Local variables
integer    :: i, ier, iostatus
real (dp)  :: fit, lsq_resid, h, std_resid, last_resid = 0.0, x(0:nvar), &
              y, xrow(nreq), one = 1.0, zero = 0.0, sumsq
write(unit_resid) 'data file name: ', fname_dat
write(unit_resid, 900) (list(i), vname(list(i)), i=1,nreq)
900 format('variables in model:' / (5(i3, ' ', a8, ' | ')))
write(unit_resid, 910) beta(1:nreq)
910 format('regression coefficients used:' / (5g15.6))
write(unit_resid) 'actual y  fitted y  ls-residual  std-residual  leverage'
x(0) = one
last_resid = zero
autoc1 = zero
sumsq = zero
do
  read(unit_data, *, iostat=iostatus) x(1:ypos-1), y, x(ypos:nvar)
  if (iostatus > 0) then
    cycle
  else if (iostatus < 0) then
    exit
  end if
  fit = zero
  do i = 1, nreq
    xrow(i) = x(list(i))
    fit = fit + beta(i) * xrow(i)
  end do
  lsq_resid = y - fit
  call hdiag(xrow, nreq, h, ier)
  std_resid = lsq_resid / (stdev * sqrt(one - h))
  write(unit_resid, 920) y, fit, lsq_resid, std_resid, h
  920 format(3g12.4, f9.2, f10.3)
  sumsq = sumsq + lsq_resid**2
  autoc1 = autoc1 + lsq_resid * last_resid
  last_resid = lsq_resid
end do
autoc1 = autoc1 / sumsq
write(*, '(a, f9.3)') 'lag 1 auto-correlation of residuals = ', autoc1
write(unit_resid, '(a, f9.3)') 'lag 1 auto-correlation of residuals = ', autoc1
return
end subroutine calc_residuals
end module reg_mod