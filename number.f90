module number_mod
use iso_fortran_env
implicit none
private
public :: ikind, sum_of_divisors, perfect, amicable, print_aliquot
integer, parameter :: ikind = int64
contains
elemental function sum_of_divisors(n) result(sum)
  integer(kind=ikind), intent(in) :: n
  integer(kind=ikind) :: sum, i, other_divisor
  if (n < 2) then
     sum = 0
     return
  end if
  sum = 1  ! 1 is always a divisor
  do i = 2, int(sqrt(real(n)))
    if (mod(n, i) == 0) then
      sum = sum + i
      other_divisor = n / i
      if (i /= other_divisor) then  ! Avoid adding the same divisor twice for perfect squares
        sum = sum + other_divisor
      end if
    end if
  end do
end function sum_of_divisors
!
elemental function perfect(n) result(is_perfect)
integer(kind=ikind), intent(in) :: n
logical             :: is_perfect
is_perfect = (n == sum_of_divisors(n))
end function perfect
!
elemental function amicable(i, j) result(are_amicable)
integer(kind=ikind), intent(in) :: i, j
logical :: are_amicable
are_amicable = all([i, j] == sum_of_divisors([j, i]))
end function amicable
!
subroutine print_aliquot(n)
integer(kind=ikind), intent(in) :: n
integer(kind=ikind) :: i, isum, kprint
integer(kind=ikind), parameter :: kprint_max = 30
i = n
write (*,"('aliquot sequence of ',i0,': ',i0)", advance="no") n,n
do kprint=1,kprint_max
   isum = sum_of_divisors(i)
   if (isum == n .or. isum == i) exit
   write (*,"(1x,i0)", advance="no") isum
   if (isum == 0) exit
   i = isum
end do
print*
end subroutine print_aliquot
end module number_mod