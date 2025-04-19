module mc_func_mod
use kind_mod, only: dp
implicit none
private
public :: ifunc, f, f_integral, f_name
integer, parameter :: ifunc = 3
contains
elemental function f(x) 
! function to integrate
real(kind=dp), intent(in) :: x
real(kind=dp)             :: f
select case (ifunc)
   case (1) ; f = x**2
   case (2) ; f = exp(x)
   case (3) ; f = log(1+x)
end select
end function f
!
elemental function f_integral(jfunc)
! value of integral from 0 to 1 of f
integer, intent(in) :: jfunc
real(kind=dp)       :: f_integral
select case (jfunc)
   case (1) ; f_integral = 1/3.0_dp
   case (2) ; f_integral = exp(1.0_dp) - 1
   case (3) ; f_integral = 2.0_dp * log(2.0_dp) - 1
end select
end function f_integral
!
elemental function f_name(jfunc)
! name of function
integer, intent(in) :: jfunc
character (len=10)  :: f_name
select case (jfunc)
   case (1) ; f_name = "x^2"
   case (2) ; f_name = "exp(x)"
   case (3) ; f_name = "log(1+x)"
end select
end function f_name
end module mc_func_mod
