module student_t_mod
use constants_mod, only: pi
implicit none
integer, parameter :: dp = kind(1.0d0)
real(kind=dp), parameter :: bad_real = -999.0_dp
public :: dp,student_t_density,student_t_variance,student_t_kurtosis,tstd,tstd_x_pow,t5
contains
elemental function student_t_variance(dof) result(var)
! variance of Student-t distribution
real(kind=dp), intent(in) :: dof
real(kind=dp)             :: var
if (dof > 2) then
   var = dof/(dof-2)
else
   var = -1.0_dp
end if
end function student_t_variance
!
elemental function student_t_kurtosis(dof) result(kurt)
! excess kurtosis of Student-t distribution
real(kind=dp), intent(in) :: dof
real(kind=dp)             :: kurt
if (dof > 4) then
   kurt = 6/(dof-4)
else
   kurt = -999.0_dp
end if
end function student_t_kurtosis
!
elemental function tstd(x,dof) result(y)
! density of standardized Student-t distribution with dof degrees of freedom at x
real(kind=dp), intent(in) :: x   ! value where density computed
real(kind=dp), intent(in) :: dof ! degrees of freedom
real(kind=dp)             :: y   ! density
real(kind=dp)             :: xscale,xscaled
if (dof > 2) then
   xscale  = sqrt(dof/(dof-2))
else
   xscale = 1.0_dp
end if
xscaled = xscale*x
y = xscale*student_t_density(xscaled,dof)
end function tstd
!
elemental function tstd_x_pow(x,dof,ipow) result(y)
! density of standardized Student-t distribution
real(kind=dp), intent(in) :: x
real(kind=dp), intent(in) :: dof ! degrees of freedom
integer      , intent(in) :: ipow
real(kind=dp)             :: y
y = (x**ipow)*tstd(x,dof)
end function tstd_x_pow
!
elemental function student_t_density(x,dof,mu,xsd) result(y)
! likelihood of Student-t density with mean mu and standard deviation xsd
real(kind=dp), intent(in)           :: x   ! argument of density
real(kind=dp), intent(in), optional :: mu  ! mean of distribution
real(kind=dp), intent(in), optional :: xsd ! standard deviation of distribution
real(kind=dp), intent(in)           :: dof ! degrees of freedom
real(kind=dp)                       :: y   ! likelihood
real(kind=dp)                       :: a,b
if (present(mu)) then
   a = mu
else
   a = 0.0_dp
end if
if (present(xsd)) then
   if (xsd <= 0.0_dp) then
      y = bad_real
      return
   end if
   if (dof > 2) then
      b = xsd/sqrt(student_t_variance(dof))
   else
      b = xsd
   end if
else
   b = 1.0_dp
end if
call student_pdf(x,a,b=b,c=dof,pdf=y)
end function student_t_density

pure subroutine student_pdf(x,a,b,c,pdf)
!*****************************************************************************80
!! STUDENT_PDF evaluates the central Student T PDF.
!  Discussion:
!    PDF(A,B,C;X) = Gamma ( (C+1)/2 ) /
!      ( Gamma ( C / 2 ) * Sqrt ( PI * C )
!      * ( 1 + ((X-A)/B)^2/C )^(C + 1/2 ) )
!  Licensing: This code is distributed under the GNU LGPL license.
!  Modified: 02 November 2005
!  Author: John Burkardt
!  Parameters:
!    Input, real ( kind = 8 ) X, the argument of the PDF.
!
!    Input, real ( kind = 8 ) A, B, shape parameters of the PDF,
!    used to transform the argument X to a shifted and scaled
!    value Y = ( X - A ) / B.  It is required that B be nonzero.
!    For the standard distribution, A = 0 and B = 1.
!
!    Input, real ( kind = 8 ) C, is usually called the number of
!    degrees of freedom of the distribution.  C is typically an
!    integer, but that is not essential.  It is required that
!    C be strictly positive.
!
!    Output, real( kind = 8) PDF, the value of the PDF.
real(kind=dp), intent(in)  :: x,a,b,c
real(kind=dp), intent(out) :: pdf
real(kind=dp)              :: y
y = (x-a)/b
pdf = gamma(0.5D+00 *(c + 1.0D+00)) /(b * sqrt(pi * c) * gamma(0.5D+00 * c) &
    * sqrt((1.0D+00 + y * y / c) **(c + 1.0D+00)))
end subroutine student_pdf
!
elemental function t5(x) result(y)
! density of Student's t with 5 degrees of freedom
real(kind=dp), intent(in) :: x
real(kind=dp)             :: y
y = 8/(3*pi*sqrt(5.0_dp)*(1+0.2_dp*x**2)**3)
end function t5
end module student_t_mod
