module kernels_mod
use kind_mod     , only: dp
use constants_mod, only: one_over_sqrt_two_pi,pi_over_2,pi_over_4,sqrt_two,two_over_pi
implicit none
private
public :: weight
real(kind=dp), parameter :: tiny_real = 1.0d0
contains
elemental function weight(x,kernel) result(y)
! kernels for nonparametric regression from https://en.wikipedia.org/wiki/Kernel_(statistics)
real(kind=dp)    , intent(in) :: x
character (len=*), intent(in) :: kernel
real(kind=dp)                 :: y
if (any(["uniform     ","triangular  ","epanechnikov","quartic     ","triweight   ", &
         "tricube     ","cosine      ","1/d         ","1/d^2       "] == kernel)) then
   if (abs(x) >= 1.0_dp) then
      y = 0.0_dp
      return
   end if
end if
select case (kernel)
   case ("uniform")     ; y = 0.5_dp
   case ("triangular")  ; y = 1.0_dp - abs(x)
   case ("epanechnikov"); y = 0.75_dp*(1.0_dp-x**2)
   case ("quartic")     ; y = 0.9375_dp*(1.0_dp-x**2)**2 ! also known as biweight
   case ("triweight")   ; y = 1.09375_dp*(1.0_dp-x**2)**3
   case ("tricube")     ; y = 0.86419753086_dp*(1.0_dp-abs(x)**3)**3 ! 0.86419753086 = 70/81
   case ("gaussian")    ; y = one_over_sqrt_two_pi*(exp(-0.5*x**2))
   case ("cosine")      ; y = pi_over_4*cos(pi_over_2*x)
   case ("logistic")    ; y = 1.0_dp/(exp(x) + 2.0_dp + exp(-x))
   case ("sigmoid")     ; y = two_over_pi/(exp(x) + exp(-x))
   case ("exponential") ; y = exp(-abs(x))
   case ("silverman")   ; y = 0.5_dp*exp(-abs(x)/sqrt_two)*cos(abs(x)/sqrt_two + pi_over_4)
   case ("1/d")         ; y = 1/(abs(x) + tiny_real)
   case ("1/d^2")       ; y = 1/(x**2 + tiny_real)
   case default         ; y = -huge(x) ! should not get here
end select
end function weight
end module kernels_mod
