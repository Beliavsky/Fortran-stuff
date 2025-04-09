module activation_mod
! 10/25/2019 09:35 PM Activation functions from https://en.wikipedia.org/wiki/Activation_function
use kind_mod, only: dp
implicit none
private
public :: activation_name,linear,ramp,logistic,elu,relu,leaky_relu,step,gaussian,sinc,soft_sign,soft_clip,inverse_square_root, &
          inverse_square_root_linear,arsinh,sqnl,gelu,softplus,bent_identity,activation_basis,sq_rbf
! activation functions not coded because they are intrinsic functions: atanh, tanh, sin
real(kind=dp), parameter :: sqrt_two = 1.414213562373095_dp, bad_real = -999.0_dp
integer, public, parameter :: ilinear = 1, iramp = 2, ilogistic = 3, ielu = 4, irelu = 5, ileaky_relu = 6, istep = 7,igaussian=8, &
                              isinc = 9, isoft_sign = 10, isoft_clip = 11, iinverse_square_root = 12, &
                              iinverse_square_root_linear = 13, iarsinh = 14, isqnl = 15, &
                              igelu = 16, isoftplus = 17, ibent_identity = 18, isqrbf = 19, ix2 = 20, ix3 = 21, nactivations = 21
character (len=10), parameter :: activation_names(nactivations) = ["linear    ", & !  1 ilinear
                                                         "ramp      ", & !  2 iramp
                                                         "logistic  ", & !  3 ilogistic
                                                         "ELU       ", & !  4 ielu
                                                         "RELU      ", & !  5 irelu
                                                         "leaky_RELU", & !  6 ileaky_relu
                                                         "step      ", & !  7 istep
                                                         "gaussian  ", & !  8 igaussian
                                                         "sinc      ", & !  9 isinc
                                                         "soft_sign ", & ! 10 isoft_sign
                                                         "soft_clip ", & ! 11 isoft_clip
                                                         "invsqrt   ", & ! 12 iinverse_square_root
                                                         "invsqrtlin", & ! 13 iinverse_square_root_linear
                                                         "arsinh    ", & ! 14 iarsinh
                                                         "sqr_nonlin", & ! 15 isqnl
                                                         "GELU      ", & ! 16 igelu
                                                         "softplus  ", & ! 17 isoftplus
                                                         "bent_ident", & ! 18 ibent_identity
                                                         "sqrbf     ", & ! 19 isqrbf
                                                         "x^2+      ", & ! 20 x^2+
                                                         "x^3+      "]   ! 21 x^3+
contains
!
elemental function activation_name(iact) result(aname)
integer, intent(in) :: iact
character (len=10)  :: aname
if (iact < 1 .or. iact > size(activation_names)) then
   aname = "???"
else
   aname = activation_names(iact)
end if
end function activation_name
!
pure function activation_basis(xx,iact,xknots,xscale) result(yy)
real(kind=dp), intent(in) :: xx(:)
integer      , intent(in) :: iact
real(kind=dp), intent(in) :: xknots(:)
real(kind=dp), intent(in) :: xscale
real(kind=dp)             :: yy(size(xx),size(xknots))
integer                   :: iknot
real(kind=dp)             :: xscale_
xscale_ = default(1.0_dp,xscale)
do iknot=1,size(xknots)
   yy(:,iknot) = activation(xscale_*(xx-xknots(iknot)),iact)
end do
end function activation_basis
!
elemental function activation(x,iact) result(y)
real(kind=dp), intent(in) :: x
integer      , intent(in) :: iact
real(kind=dp)             :: y
select case (iact)
   case (ilinear)                    ; y = linear(x)
   case (iramp)                      ; y = ramp(x)
   case (ilogistic)                  ; y = logistic(x)
   case (ielu)                       ; y = elu(x)
   case (irelu)                      ; y = relu(x)
   case (ileaky_relu)                ; y = leaky_relu(x)
   case (istep)                      ; y = step(x)
   case (igaussian)                  ; y = gaussian(x)
   case (isinc)                      ; y = sinc(x)
   case (isoft_sign)                 ; y = soft_sign(x)
   case (isoft_clip)                 ; y = soft_clip(x)
   case (iinverse_square_root)       ; y = inverse_square_root(x)
   case (iinverse_square_root_linear); y = inverse_square_root_linear(x)
   case (iarsinh)                    ; y = arsinh(x)
   case (isqnl)                      ; y = sqnl(x)
   case (igelu)                      ; y = gelu(x)
   case (isoftplus)                  ; y = softplus(x)
   case (ibent_identity)             ; y = bent_identity(x)
   case (isqrbf)                     ; y = sq_rbf(x)
   case (ix2)                        ; y = max(x,0.0_dp)**2
   case (ix3)                        ; y = max(x,0.0_dp)**3
   case default                      ; y = bad_real
end select
end function activation
!
elemental function linear(x) result(y)
real(kind=dp), intent(in) :: x
real(kind=dp)             :: y
y = x
end function linear
!
elemental function step(x) result(y)
real(kind=dp), intent(in) :: x
real(kind=dp)             :: y
y = merge(1.0_dp,0.0_dp,x>=0.0_dp)
end function step
!
elemental function gaussian(x) result(y)
real(kind=dp), intent(in) :: x
real(kind=dp)             :: y
y = exp(-x**2)
end function gaussian
!
elemental function sinc(x) result(y)
real(kind=dp), intent(in) :: x
real(kind=dp)             :: y
if (abs(x) > 1.0d-20) then
   y = sin(x)/x
else
   y = 1.0_dp
end if
end function sinc
!
elemental function soft_sign(x) result(y)
real(kind=dp), intent(in) :: x
real(kind=dp)             :: y
y = 1.0_dp/(1.0_dp  + abs(x))
end function soft_sign
!
elemental function soft_clip(x,alpha) result(y)
real(kind=dp), intent(in) :: x
real(kind=dp)             :: y
real(kind=dp), intent(in), optional :: alpha
real(kind=dp)                       :: alpha_
alpha_ = default(1.0_dp,alpha)
y = (1 + exp(alpha_*x))/(alpha_*(1.0_dp + exp(alpha_*(x-1.0_dp))))
end function soft_clip
!
elemental function inverse_square_root(x,alpha) result(y)
real(kind=dp), intent(in)           :: x
real(kind=dp), intent(in), optional :: alpha
real(kind=dp)                       :: y
y = x/sqrt(1+default(1.0_dp,alpha)*x**2)
end function inverse_square_root
!
elemental function inverse_square_root_linear(x,alpha) result(y)
real(kind=dp), intent(in)           :: x
real(kind=dp), intent(in), optional :: alpha
real(kind=dp)                       :: y
real(kind=dp)                       :: alpha_
if (x >= 0.0_dp) then
   y = x
else
   alpha_ = default(1.0_dp,alpha)
   y      = x/sqrt(1+default(1.0_dp,alpha)*x**2)
end if
end function inverse_square_root_linear
!
elemental function arsinh(x) result(y)
real(kind=dp), intent(in) :: x
real(kind=dp)             :: y
y = log(x + sqrt(1+x**2))
end function arsinh
!
elemental function ramp(x) result(y)
real(kind=dp), intent(in) :: x
real(kind=dp)             :: y
if (x >= 1) then
   y =  1.0_dp
else if (x <= -1.0_dp) then
   y = -1.0_dp
else
   y =  x
end if
end function ramp
!
elemental function sqnl(x) result(y)
! SQuare NonLinearity
real(kind=dp), intent(in) :: x
real(kind=dp)             :: y
if (x >= 2.0_dp) then
   y =  1.0_dp
else if (x <= -2.0_dp) then
   y = -1.0_dp
else if (x >= 0.0_dp) then
   y = x - 0.25_dp*x**2
else
   y = x + 0.25_dp*x**2
end if
end function sqnl
!
elemental function logistic(x) result(y)
real(kind=dp), intent(in) :: x
real(kind=dp)             :: y
y = 1.0_dp/(1.0_dp + exp(-x))
end function logistic
!
elemental function elu(x,alpha) result(y)
real(kind=dp), intent(in)           :: x
real(kind=dp), intent(in), optional :: alpha
real(kind=dp)                       :: y
if (x >= 0) then
   y = x
else
   y = exp(x) - 1.0_dp
   if (present(alpha)) y = alpha*y
end if
end function elu
!
elemental function relu(x) result(y)
! REctified Linear Unit
real(kind=dp), intent(in) :: x
real(kind=dp)             :: y
y = max(x,0.0_dp)
end function relu
!
elemental function gelu(x) result(y)
! Gaussian Error Linear Unit
real(kind=dp), intent(in) :: x
real(kind=dp)             :: y
y = 0.5_dp*x*(1.0_dp+erf(x/sqrt_two))
end function gelu
!
elemental function softplus(x) result(y)
real(kind=dp), intent(in) :: x
real(kind=dp)             :: y
y = log(1.0_dp+exp(x))
end function softplus
!
elemental function bent_identity(x) result(y)
real(kind=dp), intent(in) :: x
real(kind=dp)             :: y
y = 0.5_dp*(sqrt(x**2+1.0_dp) - 1.0_dp) + x
end function bent_identity
!
elemental function leaky_relu(x,alpha) result(y)
! Leaky REctified Linear Unit
real(kind=dp), intent(in)           :: x
real(kind=dp), intent(in), optional :: alpha
real(kind=dp)                       :: y
real(kind=dp), parameter            :: alpha_default = 0.01_dp
if (x >= 0.0_dp) then
   y = x
else
   if (present(alpha)) then
      y = alpha*x
   else
      y = alpha_default*x
   end if
end if
end function leaky_relu
!
elemental function sq_rbf(x) result(y)
! from https://datascience.stackexchange.com/questions/62241/square-law-based-rbf-kernel
real(kind=dp), intent(in) :: x
real(kind=dp)             :: y
real(kind=dp)             :: xabs
xabs = abs(x)
if (xabs <= 1.0_dp) then
   y = 1 - 0.5_dp*xabs**2
else if (xabs <= 2.0_dp) then
   y = 0.5_dp*(xabs-2)**2  
else
   y = 0.0_dp
end if
end function sq_rbf
!
elemental function default(def,opt) result(yy)
! return opt if it is present, otherwise def (real arguments and result)
real(kind=dp), intent(in)           :: def
real(kind=dp), intent(in), optional :: opt
real(kind=dp)                       :: yy
if (present(opt)) then
   yy = opt
else
   yy = def
end if
end function default
end module activation_mod
