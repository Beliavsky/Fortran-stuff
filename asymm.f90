module asymm_mod
implicit none
integer, parameter :: dp = kind(1.0d0)
private
public :: asymm
contains
pure function asymm(x,c,abs_x) result(y)
! create asymmetry for symmetric zero-mean, unit-variance data by differentially scaling positive and negative values of x(:)
real(kind=dp), intent(in) :: x(:)
real(kind=dp), intent(in) :: c     ! asymmetry parameter: > 1 to create positively skewed data, 0 < c < 1 for negatively skewed data
real(kind=dp), intent(in) :: abs_x ! known absolute mean
real(kind=dp)             :: y(size(x))
real(kind=dp)             :: ymean,y2,yvar
ymean = (c-(1/c)) * abs_x/2  ! mean of scaled data
y2    = (c**2 + 1/(c**2))/2  ! 2nd moment of scaled data
yvar  = y2 - ymean**2        ! variance of scaled data
y     = (merge(c,1/c,x>0) * x - ymean)/sqrt(yvar) ! shift and scale y to have zero mean and unit variance
end function asymm
end module asymm_mod