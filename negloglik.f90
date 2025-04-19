module negloglik_mod
implicit none
private
public :: negloglik
contains
subroutine negloglik(x, n, mu, sigma, nll, imethod)
  integer         , intent(in)  :: n,imethod
  double precision, intent(in)  :: x(n), mu, sigma
  double precision, intent(out) :: nll
  double precision, parameter   :: pi = 4*atan(1.0d0), sqrt_2_pi = sqrt(2.0d0*pi)
  double precision              :: z(n),ll,zscalar,log_sigma_sqrt_2_pi
  integer                       :: i
  if (imethod < 3) z = (x - mu)/sigma
  if (imethod == 1) then
     ll = 0.0d0
     do i=1,n
        ll = ll - 0.5*z(i)**2 - log(sigma*sqrt_2_pi)
     end do
  else if (imethod == 2) then ! avoid computing log(sigma*sqrt(2*pi)) in loop
     ll = -n * log(sigma*sqrt_2_pi) - 0.5d0*sum(z**2)
  else if (imethod == 3) then ! loop fusion
     ll = 0.0d0
     do i=1,n
       zscalar = (x(i) - mu)/sigma
       ll = ll - 0.5d0*zscalar**2 - log(sigma*sqrt_2_pi)
     end do
  else if (imethod == 4) then ! better loop fusion with log outside loop
     ll = -n * log(sigma*sqrt_2_pi)
     do i=1,n
       zscalar = (x(i) - mu)/sigma
       ll = ll - 0.5d0*zscalar**2
     end do
  else if (imethod == 5) then ! better with 0.5d0/sigma**2 outside loop
     ll = 0.0d0
     do i=1,n
        ll = ll - ((x(i) - mu))**2
     end do
     ll = ll * 0.5d0 / (sigma**2)
     ll = ll - n * log(sigma*sqrt_2_pi)
  end if
  nll = -ll
end subroutine negloglik
end module negloglik_mod