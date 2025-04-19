module primes_mod
implicit none
private
public :: partially_factorize, factorize, factorize_pollard
contains
subroutine partially_factorize(n, primes, powers, factor)
! Partially factors an integer n into powers of the given primes and a remaining factor
integer, intent(in)  :: n  ! integer to factor
integer, intent(in)  :: primes(:)   ! array of primes
integer, intent(out) :: powers(size(primes))  ! powers of the primes in n
integer, intent(out) :: factor  ! Remaining factor after extracting powers of primes
integer :: i, p
factor = n   ! Initialize remaining factor as n
if (size(primes) < 1) return
powers = 0
! Loop over each prime and compute its power in the factorization of n
do i = 1, size(primes)
    p = primes(i)
    ! Divide by the prime as long as the number is divisible
    do while (mod(factor, p) == 0)
        factor = factor / p
        powers(i) = powers(i) + 1
    end do
end do
end subroutine partially_factorize
!
subroutine factorize(n, primes, powers)
! Factors an integer n into powers of primes using trial division
integer, intent(in) :: n  ! integer to factor
integer, intent(out), allocatable :: primes(:)  ! array of prime factors
integer, intent(out), allocatable :: powers(:)  ! powers of the primes in n
integer :: ncount, remaining, p, max_prime, nfac
integer, parameter :: max_fac = 10**5
integer :: prime_factors(max_fac), prime_powers(max_fac)
! Initialize variables
remaining = n
max_prime = 2  ! Start with 2
nfac = 0
! Try to divide n by all numbers starting from 2 upwards
p = 2
loop_p: do while (p * p <= remaining)
   ncount = 0
   ! Divide n by p as long as it's divisible by p
   do while (mod(remaining, p) == 0)
      remaining = remaining / p
      ncount = ncount + 1
   end do   
   ! If p is a factor, store it and its count (power)
   if (ncount > 0) then
      nfac = nfac + 1
      prime_factors(nfac) = p
      prime_powers(nfac) = ncount
      if (nfac == max_fac) exit loop_p
   end if
   ! Move to the next potential factor
   if (p == 2) then
      p = 3
   else
      p = p + 2  ! Skip even numbers
   end if
end do loop_p
! If remaining > 1, it is a prime factor
if (remaining > 1 .and. nfac < max_fac) then
   nfac = nfac + 1
   prime_factors(nfac) = remaining
   prime_powers(nfac) = 1
end if
primes = prime_factors(:nfac)
powers = prime_powers(:nfac)
end subroutine factorize
!
subroutine factorize_pollard(n, primes, powers)
    ! Factors an integer n into powers of primes using Pollard's rho algorithm
    implicit none
    integer, intent(in) :: n  ! integer to factor
    integer, intent(out), allocatable :: primes(:)  ! array of prime factors
    integer, intent(out), allocatable :: powers(:)  ! powers of the primes in n
    integer :: remaining, factor, nfac
    integer, parameter :: max_fac = 10**5
    integer :: prime_factors(max_fac), prime_powers(max_fac)
    ! Initialize variables
!    print*,"entered factorize_pollard, n =", n
    remaining = n
    nfac = 0
    ! Use trial division to remove small primes
    call factorize_small_primes(remaining, prime_factors, prime_powers, nfac, max_fac)
!     print*,"returned from factorize_small_primes, nfac =", nfac
!     print*,"prime_factors =", prime_factors(:nfac)
!     print*,"prime_powers =", prime_powers(:nfac)
!     print*,"remaining =", remaining
    ! Apply Pollard's Rho algorithm for remaining composite number
    do while (remaining > 1 .and. nfac < max_fac)
        factor = pollard_rho(remaining)
!        print*,"remaining, factor:", remaining, factor
        if (factor == remaining) then
            ! If Pollard's Rho fails to find a factor, it's likely prime
            nfac = nfac + 1
            prime_factors(nfac) = remaining
            prime_powers(nfac) = 1
            exit
        else
            ! Apply trial division to fully factorize 'factor'
            call factorize_small_primes(factor, prime_factors, prime_powers, nfac, max_fac)
            remaining = remaining / factor
        end if
    end do

    ! Output the results
    primes = prime_factors(:nfac)
    powers = prime_powers(:nfac)
end subroutine factorize_pollard
!
integer function pollard_rho(n)
    integer, intent(in) :: n
    integer :: x, y, d, c, iteration
    integer, parameter :: iteration_limit = 10000
    ! Return immediately for small values of n
    if (any(n == [1, 2, 3, 5])) then
       pollard_rho = n
       return
    end if
    ! Initialize Pollard's Rho variables
    x = 2
    y = 2
    c = 1  ! constant in f(x) = x^2 + c
    ! Pollard's Rho algorithm loop
    do iteration = 1, iteration_limit
        x = mod(x * x + c, n)  ! f(x) = (x^2 + c) mod n
        y = mod(y * y + c, n)
        y = mod(y * y + c, n)  ! f(f(y)) to make it faster
        d = gcd(abs(x - y), n)

        if (d /= 1 .and. d /= n) then
            pollard_rho = d
            return
        end if
    end do

    ! If the iteration limit is reached, return n as it is likely a prime
    pollard_rho = n
end function pollard_rho
    ! GCD function
    integer function gcd(a, b)
        integer, intent(in) :: a, b
        integer :: temp_a, temp_b
        temp_a = a
        temp_b = b

        do while (temp_b /= 0)
            gcd = temp_b
            temp_b = mod(temp_a, temp_b)
            temp_a = gcd
        end do

        gcd = temp_a
    end function gcd

    ! Small prime factorization (using trial division)
    subroutine factorize_small_primes(n, prime_factors, prime_powers, nfac, max_fac)
        integer, intent(inout) :: n  ! integer to factor
        integer, intent(out) :: prime_factors(:)  ! array of prime factors
        integer, intent(out) :: prime_powers(:)  ! powers of the primes in n
        integer, intent(inout) :: nfac  ! number of factors found
        integer, intent(in) :: max_fac  ! maximum number of factors allowed
        integer :: p, ncount
!        print*,"entered factorize_small_primes, n =",n
        ! Try dividing by small primes starting from 2
        p = 2
        do while (p * p <= n .and. nfac < max_fac)
            ncount = 0
!            print*,"here in loop n, p, nfac =", n, p, nfac
            do while (mod(n, p) == 0)
                n = n / p
                ncount = ncount + 1
            end do
            if (ncount > 0) then
                nfac = nfac + 1
                prime_factors(nfac) = p
                prime_powers(nfac) = ncount
            end if
!            print*,"n, p, ncount:", n, p, ncount
            ! Move to the next potential prime factor
            if (p == 2) then
                p = 3
            else
                p = p + 2
            end if
        end do
!        print*,"exited while loop, n, p, ncount, nfac =", n, p, ncount, nfac
        ! If remaining > 1 and nfac is less than max_fac, the remaining part is prime
        if (n > 1 .and. nfac < max_fac) then
            nfac = nfac + 1
            prime_factors(nfac) = n
            prime_powers(nfac) = 1
        end if
!         print*,"n, nfac =", n, nfac
!         print*,"prime_factors =", prime_factors(:nfac)
!         print*,"prime_powers =", prime_powers(:nfac)
!         print*,"leaving factorize_small_primes"
    end subroutine factorize_small_primes
end module primes_mod
