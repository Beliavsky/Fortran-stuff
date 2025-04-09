module pdf_mod
use kind_mod     , only: dp
use constants_mod, only: log_two_pi,pi_over_2,sqrt_two,one_over_sqrt_two_pi
use student_t_mod, only: student_t_density
use util_mod     , only: sech,bad_real
implicit none
private
public :: laplace_density,dist_log_lik,base_hyperbolic_secant_density, &
          normal_log_lik,std_density,normal_log_lik_var
character (len=*), parameter, public :: laplace_str = "Laplace", normal_str = "normal", &
   sech_str = "sech", t1_str = "t1", t2_str = "t2", t3_str = "t3", t4_str = "t4", &
   t5_str = "t5", t6_str = "t6", t7_str = "t7", t8_str = "t8", t9_str = "t9", &
   t10_str = "t10"     
contains
elemental function base_hyperbolic_secant_density(xx) result(lik)
! hyperbolic secant density
real(kind=dp), intent(in) :: xx  ! argument of density
real(kind=dp)             :: lik ! likelihood
lik = sech(pi_over_2*xx)/2
end function base_hyperbolic_secant_density
!
elemental function normal_log_lik(xx,mu,xsd) result(log_lik)
! log likelihood of normal density with mean mu and standard deviation xsd
real(kind=dp), intent(in) :: xx      ! argument of density
real(kind=dp), intent(in) :: mu      ! mean of distribution
real(kind=dp), intent(in) :: xsd     ! standard deviation of distribution
real(kind=dp)             :: log_lik ! log likelihood
log_lik = - 0.5_dp*log_two_pi - log(xsd) - 0.5_dp*(((xx-mu)/xsd)**2)
end function normal_log_lik
!
elemental function normal_log_lik_var(xx,xvar) result(log_lik)
! log likelihood of normal density with mean zero and standard deviation xvar
real(kind=dp), intent(in) :: xx      ! argument of density
real(kind=dp), intent(in) :: xvar    ! variance of distribution
real(kind=dp)             :: log_lik ! log likelihood
log_lik = - 0.5_dp*log_two_pi - 0.5_dp*log(xvar) - 0.5_dp*(xx**2/xvar)
end function normal_log_lik_var
!
elemental function laplace_density(xx,mu,xsd) result(lik)
! likelihood of Laplace density with mean mu and standard deviation xsd
real(kind=dp), intent(in) :: xx      ! argument of density
real(kind=dp), intent(in) :: mu      ! mean of distribution
real(kind=dp), intent(in) :: xsd     ! standard deviation of distribution
real(kind=dp)             :: lik     ! log likelihood
lik = 1/(sqrt(2.0_dp)*xsd)*exp(-sqrt(2.0_dp)*abs((xx-mu)/xsd))
end function laplace_density
!
elemental function laplace_log_lik(xx,mu,xsd) result(log_lik)
! log likelihood of normal density with mean mu and standard deviation xsd
real(kind=dp), intent(in) :: xx! argument of density
real(kind=dp), intent(in) :: mu! mean of distribution
real(kind=dp), intent(in) :: xsd     ! standard deviation of distribution
real(kind=dp)             :: log_lik ! log likelihood
log_lik = -0.5_dp*log(2.0_dp) - log(xsd) - sqrt(2.0_dp)*abs((xx-mu)/xsd)
end function laplace_log_lik
!
elemental function hyperbolic_secant_density(xx,mu,xsd) result(lik)
! likelihood of Laplace density with mean mu and standard deviation xsd
real(kind=dp), intent(in) :: xx      ! argument of density
real(kind=dp), intent(in) :: mu      ! mean of distribution
real(kind=dp), intent(in) :: xsd     ! standard deviation of distribution
real(kind=dp)             :: lik     ! log likelihood
lik = base_hyperbolic_secant_density((xx-mu)/xsd)/xsd
end function hyperbolic_secant_density
!
elemental function hyperbolic_secant_log_lik(xx,mu,xsd) result(log_lik)
! log likelihood of hyperbolic secant density with mean mu and standard deviation xsd
real(kind=dp), intent(in) :: xx      ! argument of density
real(kind=dp), intent(in) :: mu      ! mean of distribution
real(kind=dp), intent(in) :: xsd     ! standard deviation of distribution
real(kind=dp)             :: log_lik ! log likelihood
log_lik = log(hyperbolic_secant_density(xx,mu,xsd))
end function hyperbolic_secant_log_lik 
!
elemental function dist_log_lik(xx,dist,mu,xsd) result(log_lik)
! likelihood of dist with mean mu and standard deviation xsd
real(kind=dp)    , intent(in) :: xx      ! argument of density
character (len=*), intent(in) :: dist    ! name of distribution
real(kind=dp)    , intent(in) :: mu      ! mean of distribution
real(kind=dp)    , intent(in) :: xsd     ! standard deviation of distribution
real(kind=dp)                 :: log_lik ! log likelihood
select case (dist)
   case (laplace_str); log_lik = laplace_log_lik(xx,mu,xsd)
   case (normal_str) ; log_lik =  normal_log_lik(xx,mu,xsd)
   case (sech_str)   ; log_lik = hyperbolic_secant_log_lik(xx,mu,xsd)
   case (t1_str)     ; log_lik = log(student_t_density(xx,dof=1.0_dp,mu=mu,xsd=xsd))
   case (t2_str)     ; log_lik = log(student_t_density(xx,dof=2.0_dp,mu=mu,xsd=xsd))
   case (t3_str)     ; log_lik = log(student_t_density(xx,dof=3.0_dp,mu=mu,xsd=xsd))
   case (t4_str)     ; log_lik = log(student_t_density(xx,dof=4.0_dp,mu=mu,xsd=xsd))
   case (t5_str)     ; log_lik = log(student_t_density(xx,dof=5.0_dp,mu=mu,xsd=xsd))
   case (t6_str)     ; log_lik = log(student_t_density(xx,dof=6.0_dp,mu=mu,xsd=xsd))
   case (t7_str)     ; log_lik = log(student_t_density(xx,dof=7.0_dp,mu=mu,xsd=xsd))
   case (t8_str)     ; log_lik = log(student_t_density(xx,dof=8.0_dp,mu=mu,xsd=xsd))
   case (t9_str)     ; log_lik = log(student_t_density(xx,dof=9.0_dp,mu=mu,xsd=xsd))
   case (t10_str)    ; log_lik = log(student_t_density(xx,dof=10.0_dp,mu=mu,xsd=xsd))
   case default      ; log_lik = bad_real
end select
end function dist_log_lik
!
elemental function std_density(xx,dist) result(yy)
real(kind=dp)    , intent(in) :: xx
character (len=*), intent(in) :: dist
real(kind=dp)                 :: yy
select case (dist)
   case ("normal")  ; yy = one_over_sqrt_two_pi * exp(-0.5*xx**2)
   case ("Laplace") ; yy = 1/(sqrt_two)*exp(-sqrt_two*abs(xx))
   case ("sech")    ; yy = sech(pi_over_2*xx)/2
   case ("t5")      ; yy = student_t_density(xx,dof=5.0_dp,mu=0.0_dp,xsd=1.0_dp)
   case ("t10")     ; yy = student_t_density(xx,dof=10.0_dp,mu=0.0_dp,xsd=1.0_dp)
   case default     ; yy = bad_real
end select
end function std_density
end module pdf_mod