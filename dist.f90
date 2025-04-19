module dist_mod
use kind_mod, only: dp
implicit none
private
public :: pdf
real(kind=dp), parameter :: bad_real = -999.0_dp
contains
elemental function pdf(x,dist,c) result(y)
! probability density functions
real(kind=dp)    , intent(in) :: x
real(kind=dp)    , intent(in) :: c
character (len=*), intent(in) :: dist
real(kind=dp)                 :: y
select case (dist)
!  distributions without a shape parameter  
   case ("normal")              ; y = exp(-x**2)               ! kurtosis 0          
   case ("Laplace")             ; y = exp(-abs(x))             ! kurtosis 3
   case ("Cauchy")              ; y = 1/(1+x**2)               
   case ("logistic")            ; y = 1/(exp(x) + 2 + exp(-x)) ! kurtosis 1.2
   case ("sech")                ; y = 1/(exp(x) + exp(-x))     ! kurtosis 2
!  distributions with a shape parameter c
   case ("generalized_error")   ; y = exp(-(abs(x))**c)        ! normal   for c = 2, Laplace for c = 1 
   case ("student_t"   )        ; y = 1/(1+(x**2)/c)**(c+1)/2  ! Cauchy   for c = 1, Normal  for large c -- kurtosis = 6/(c-4)
   case ("symmetric_hyperbolic"); y = exp(-sqrt(x**2 + c**2))  ! Laplace  for c = 0, Normal  for large c
   case ("Champernowne")        ; y = 1/(exp(x) + c + exp(-x)) ! logistic for c = 2, hyperbolic secant for c = 0; need c > -2
   case default     ; y = bad_real
end select
end function pdf
end module dist_mod