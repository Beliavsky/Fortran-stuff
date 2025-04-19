module coupon_mod
use ziggurat_pure_mod, only: dp, random_uni
implicit none
private
public :: dp, avg_tries, avg_tries_sub
contains
impure elemental function avg_tries(n, k) result(avg)
! simulate the average number of tries that it
! takes to collect all n coupons, using k trials
integer, intent(in) :: n, k
real(kind=dp) :: avg
real(kind=dp) :: rn
logical :: collected(n)
integer :: i, j, ncount, total_tries
total_tries = 0
do i = 1, k
   collected = .false.
   ncount = 0
   do while (ncount < n)
      call random_number(rn)
      j = int(rn*n) + 1
      if (.not. collected(j)) then
         ncount = ncount + 1
         collected(j) = .true.
      end if
      total_tries = total_tries + 1
   end do
end do
avg = total_tries / real(k, kind=dp)
end function avg_tries
!
pure elemental subroutine avg_tries_sub(n, k, jsr, avg)
! simulate the average number of tries that it
! takes to collect all n coupons, using k trials
integer, intent(in) :: n, k
integer, intent(in out) :: jsr
real(kind=dp), intent(out) :: avg
real(kind=dp) :: rn
logical :: collected(n)
integer :: i, j, ncount, total_tries
total_tries = 0
do i = 1, k
   collected = .false.
   ncount = 0
   do while (ncount < n)
      call random_uni(jsr, rn)
      j = int(rn*n) + 1
      if (.not. collected(j)) then
         ncount = ncount + 1
         collected(j) = .true.
      end if
      total_tries = total_tries + 1
   end do
end do
avg = total_tries / real(k, kind=dp)
end subroutine avg_tries_sub
end module coupon_mod
