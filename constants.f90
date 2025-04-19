module constants_mod
use kind_mod, only: dp
implicit none
private
public :: pi, pi_over_2, log_two_pi, log_two, sqrt_two, &
   one_over_sqrt_two_pi, pi_over_4, two_over_pi, dpmpar
real(kind=dp), parameter :: &
   pi         = 3.141592653589793238462643_dp, &
   pi_over_2  = pi/2.0_dp, &
   pi_over_4  = pi/4.0_dp, &
   two_over_pi = 2.0_dp/pi, &
   one_over_sqrt_two_pi = 0.39894228040143270_dp , &
   log_two_pi = 1.837877066409345483560659_dp, &
   log_two    = 0.69314718055994529_dp, &
   sqrt_two   = 1.4142135623730951_dp, &
   euler      = 0.5772156649015328606065120900824024_dp, &
   polygamma_one_half = 4.93480220054468_dp, & ! polygamma(1, 0.5), &
   psi_half = -1.9635100260214235_dp ! psi(0.5) from SciPy

contains

pure function dpmpar(i) result(fn_val)
!-----------------------------------------------------------------------
!     dpmpar provides the double precision machine constants for
!     the computer being used. it is assumed that the argument
!     i is an integer having one of the values 1, 2, or 3. if the
!     double precision arithmetic being used has m base b digits and
!     its smallest and largest exponents are emin and emax, then
!        dpmpar(1) = b**(1 - m), the machine precision,
!        dpmpar(2) = b**(emin - 1), the smallest magnitude,
!        dpmpar(3) = b**emax*(1 - b**(-m)), the largest magnitude.
integer      , intent(in) :: i
real(kind=dp)             :: fn_val
real(kind=dp), parameter  :: one = 1._dp
select case (i)
  case (1) ; fn_val = epsilon(one)
  case (2) ; fn_val = tiny(one)
  case (3) ; fn_val = huge(one)
end select
end function dpmpar
end module constants_mod
