module arch_fit_mod
use kind_mod, only: dp
use util_mod, only: istdout,default,write_format
implicit none
private
public :: arch_vol,arch_model,set,display
interface set
   module procedure set_arch_model
end interface set
interface display
   module procedure display_arch_model
end interface display
type :: arch_model
   integer       :: nlags
   real(kind=dp) :: geo_ratio=1.0_dp,casymm=0.0_dp,vol_const=0.0_dp,vol_mult=1.0_dp,log_lik
   real(kind=dp), allocatable :: vols(:)
end type arch_model
contains
subroutine display_arch_model(yy,outu,fmt_header,fmt_trailer)
type(arch_model) , intent(in)           :: yy
integer          , intent(in), optional :: outu
character (len=*), intent(in), optional :: fmt_header,fmt_trailer
integer                                 :: outu_
character (len=*), parameter            :: fmt_ci="(a30,':',100(1x,i10))",fmt_cr="(a30,':',100(1x,f10.4))"
outu_ = default(istdout,outu)
call write_format(fmt_header,outu)
write (outu_,fmt_ci) "#lags",yy%nlags
write (outu_,fmt_cr) "casymm",yy%casymm
write (outu_,fmt_cr) "geo_ratio",yy%geo_ratio
write (outu_,fmt_cr) "vol_const",yy%vol_const
write (outu_,fmt_cr) "vol_mult",yy%vol_mult
write (outu_,fmt_cr) "log_lik",yy%log_lik
call write_format(fmt_trailer,outu)
end subroutine display_arch_model
!
subroutine set_arch_model(yy,nlags,casymm,geo_ratio,vol_const,vol_mult,log_lik,vols)
type(arch_model), intent(in out) :: yy
integer         , intent(in), optional :: nlags
real(kind=dp)   , intent(in), optional :: casymm,geo_ratio,vol_const,vol_mult,log_lik
real(kind=dp)   , intent(in), optional :: vols(:)
if (present(nlags))     yy%nlags = nlags
if (present(casymm))    yy%casymm = casymm
if (present(geo_ratio)) yy%geo_ratio = geo_ratio
if (present(vol_const)) yy%vol_const = vol_const
if (present(vol_mult))  yy%vol_mult  = vol_mult
if (present(log_lik))   yy%log_lik   = log_lik
if (present(vols))      yy%vols = vols
end subroutine set_arch_model
!
function arch_vol(xx,wgt,vol_const,vol_init) result(xvol)
real(kind=dp), intent(in) :: xx(:)
real(kind=dp), intent(in) :: wgt(:)
real(kind=dp), intent(in) :: vol_const
real(kind=dp), intent(in) :: vol_init
real(kind=dp)             :: xvol(size(xx))
real(kind=dp)             :: past_sum_squares
integer                   :: i,n,nwgt
n = size(xx)
nwgt = size(wgt)
do i=1,n
   if (i <= nwgt) then
      xvol(i) = vol_init
   else
      past_sum_squares = sum(wgt(nwgt:1:-1)*xx(i-nwgt:i-1)**2)
      xvol(i) = sqrt(max(0.0_dp,vol_const**2 + past_sum_squares))
   end if
end do
end function arch_vol
end module arch_fit_mod
