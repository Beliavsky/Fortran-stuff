module mix_mod
use kind_mod, only: dp
use random_mod, only: random_normal
use constants_mod, only: one_over_sqrt_two_pi
private
public :: normal_mix_sim_2_components, fit_normal_mix_2_components
contains
elemental function normal_pdf(x, xmean, xsd) result(y)
real(kind=dp), intent(in) :: x, xmean, xsd
real(kind=dp)             :: y
y = (one_over_sqrt_two_pi/xsd) * exp(-0.5*((x-xmean)/xsd)**2)
end function normal_pdf
!
function normal_mix_sim_2_components(n, p1, m1, m2, sd1, sd2) result(x)
integer      , intent(in) :: n
real(kind=dp), intent(in) :: p1       ! probability of 1st component
real(kind=dp), intent(in) :: m1, m2   ! means
real(kind=dp), intent(in) :: sd1, sd2 ! standard deviations
real(kind=dp)             :: x(n)
real(kind=dp)             :: u, z
integer                   :: i
do i=1,n
   call random_number(u)
   z = random_normal()
   if (u < p1) then
      x(i) = m1 + sd1*z
   else
      x(i) = m2 + sd2*z
   end if
end do
end function normal_mix_sim_2_components
!
subroutine fit_normal_mix_2_components(x, p1, m1, m2, sd1, sd2, niter)
real(kind=dp), intent(in)     :: x(:)
real(kind=dp), intent(in out) :: p1
real(kind=dp), intent(in out) :: m1, m2
real(kind=dp), intent(in out) :: sd1, sd2
integer      , intent(in), optional :: niter
real(kind=dp)                 :: p2, r1(size(x)), r2(size(x)), w1(size(x)), w2(size(x)), &
                                 sum_w1
integer                       :: n, niter_
if (present(niter)) then
   niter_ = niter
else
   niter_ = 100
end if
n = size(x)
do i=1,niter_
   p2 = 1.0_dp - p1
   r1 = p1 * normal_pdf(x, m1, sd1)
   r2 = p2 * normal_pdf(x, m2, sd2)
   w1 = r1 / (r1+r2)
   w2 = 1.0_dp - w1
   sum_w1 = sum(w1)
   sum_w2 = n - sum(w1)
   p1 = sum_w1 / n
   m1 = sum(w1*x) / sum_w1
   m2 = sum(w2*x) / sum_w2
   sd1 = sqrt(sum(w1 * (x - m1)**2) / sum_w1)
   sd2 = sqrt(sum(w2 * (x - m2)**2) / sum_w2)
end do
end subroutine fit_normal_mix_2_components
end module mix_mod