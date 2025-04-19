module util_mod
use iso_fortran_env, only: output_unit
use kind_mod, only: dp
implicit none
private
public :: default, assert_equal, write_merge, split_string, display, &
   first_false, sech, bad_real, istdout, write_format, set_optional
interface default
   module procedure default_int, default_real, default_logical, &
      default_character
end interface default
interface display
   module procedure display_matrix
end interface display
interface set_optional
   module procedure set_optional_integer
end interface set_optional
integer, parameter :: istdout = output_unit
real(kind=dp), parameter :: bad_real = -999.0_dp
contains

elemental subroutine set_optional_integer(x, xdef, xopt)
integer, intent(out) :: x
integer, intent(in)  :: xdef
integer, intent(in), optional :: xopt
if (present(xopt)) then
   x = xopt
else
   x = xdef
end if
end subroutine set_optional_integer

elemental function default_int(x, xopt) result(y)
! return xopt if present, otherwise x
integer, intent(in) :: x
integer, intent(in), optional :: xopt
integer             :: y
if (present(xopt)) then
   y = xopt
else
   y = x
end if
end function default_int

elemental function default_real(x, xopt) result(y)
! return xopt if present, otherwise x
real(kind=dp), intent(in) :: x
real(kind=dp), intent(in), optional :: xopt
real(kind=dp)             :: y
if (present(xopt)) then
   y = xopt
else
   y = x
end if
end function default_real

elemental function default_logical(x, xopt) result(y)
! return xopt if present, otherwise x
logical, intent(in) :: x
logical, intent(in), optional :: xopt
logical             :: y
if (present(xopt)) then
   y = xopt
else
   y = x
end if
end function default_logical

elemental function default_character(x, xopt) result(y)
! return xopt if present, otherwise x
character (len=*), intent(in) :: x
character (len=*), intent(in), optional :: xopt
character (len=100) :: y
if (present(xopt)) then
   y = xopt
else
   y = x
end if
end function default_character

subroutine assert_equal(k, kreq, msg)
! check that k == kreq
integer, intent(in) :: k, kreq
character (len=*), intent(in) :: msg
if (k /= kreq) then
   print "(a, i0, a, i0)", msg // " = ", k, ", must equal ", kreq
   stop
end if
end subroutine assert_equal

subroutine write_merge(tf, x, y, outu, fmt)
!> Writes either `x` or `y` to the specified output unit using the given format.
!! If `tf` is true, writes `x`; otherwise, writes `y`.
!! @param tf Logical condition determining whether to write `x` or `y`.
!! @param x The first character string to write if `tf` is true.
!! @param y The second character string to write if `tf` is false.
!! @param outu Optional output unit (defaults to a predefined output unit).
!! @param fmt Optional format specifier (defaults to "(a)").
logical, intent(in) :: tf
character (len=*), intent(in) :: x, y
integer, intent(in), optional :: outu
character (len=*), intent(in), optional :: fmt
integer :: outu_
character (len=100) :: fmt_
outu_ = default(output_unit, outu)
if (present(fmt)) then
   fmt_ = fmt
else
   fmt_ = "(a)"
end if
if (tf) then
   write (outu_, fmt_) x
else
   write (outu_, fmt_) y
end if
end subroutine write_merge

!------------------------------------------------------------------
! Utility: split_string
!
! Splits the input string 'str' at each occurrence of the single-
! character delimiter 'delim' and returns the pieces in the allocatable
! array 'tokens'. To allocate each element (with deferred length)
! properly, we use the length of the input string.
!------------------------------------------------------------------
subroutine split_string(str, delim, tokens)
character(len=*), intent(in)           :: str
character(len=*), intent(in)           :: delim
character(:), allocatable, intent(out) :: tokens(:)
integer :: start, pos, i, count, n

n = len_trim(str)
if (n == 0) then
   allocate(character(len=0) :: tokens(1))
   tokens(1) = ""
   return
end if

! First pass: count tokens.
count = 0
start = 1
do
   pos = index(str(start:), delim)
   if (pos == 0) then
      count = count + 1
      exit
   else
      count = count + 1
      start = start + pos
   end if
end do

! Allocate tokens; each token gets the full length of the input.
allocate(character(len=n) :: tokens(count))

! Second pass: extract tokens.
start = 1
i = 1
do
   pos = index(str(start:), delim)
   if (pos == 0) then
      tokens(i) = adjustl(str(start:))
      exit
   else
      tokens(i) = adjustl(str(start:start+pos-2))
      start = start + pos
      i = i + 1
   end if
end do
end subroutine split_string

subroutine display_matrix(x, outu, fmt_r, fmt_header)
! print a matrix
real(kind=dp)    , intent(in)           :: x(:,:)
integer          , intent(in), optional :: outu
character (len=*), intent(in), optional :: fmt_r, fmt_header
integer                                 :: i, outu_
character (len=100)                     :: fmt_r_
outu_  = default(output_unit, outu)
fmt_r_ = default("(*(1x,f10.4))", fmt_r)
if (present(fmt_header)) write(outu_, fmt_header)
do i=1,size(x,1)
   write(outu_,fmt_r_) x(i,:)
end do
end subroutine display_matrix

pure function first_false(tf) result(i1)
! return the location of the first false element in tf(:), 0 if all .true.
logical, intent(in) :: tf(:)
integer             :: i1
integer :: i
i1 = 0  ! Default to 0, assuming all elements are .true.
do i = 1, size(tf)
    if (.not. tf(i)) then
        i1 = i
        return
    end if
end do
end function first_false

elemental function sech(xx) result(yy)
! hyperbolic secant
real(kind=dp), intent(in) :: xx
real(kind=dp)             :: yy
yy = 1/cosh(xx)
end function sech

subroutine write_format(format_str,iunit,advance)
! write format_str to unit iunit if it is present and not blank
! otherwise, do nothing
character (len=*), intent(in), optional :: format_str
integer          , intent(in), optional :: iunit
character (len=*), intent(in), optional :: advance
integer                                 :: iu
if (.not. present(format_str)) return
if (present(iunit)) then
   iu = iunit
else
   iu = istdout
end if
! iu = default(istdout,iunit)
if (format_str /= "") write (iu,format_str,advance=default("yes",advance))
end subroutine write_format

end module util_mod

