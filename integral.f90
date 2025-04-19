! NumericalIntegration.f90
module integral_mod
    use kind_mod, only: dp
    implicit none
    public :: integral
    contains

    function integral(f, a, b, eps) result(integral_result)
    ! from ChatGPT
        interface
            real(kind=dp) function f(x)
                import dp
                real(kind=dp), intent(in) :: x
            end function f
        end interface
        real(kind=dp), intent(in) :: a, b, eps
        real(kind=dp) :: integral_result, lastSum, currentSum, trapezoids
        integer :: n, k
        lastSum = 0.0_dp
        n = 1
        do
            currentSum = 0.0_dp
            do k = 0, n
                trapezoids = a + (b - a) * real(k, dp) / real(n, dp)
                currentSum = currentSum + f(trapezoids) * (b - a) / real(n, dp)
            end do
            if (abs(currentSum - lastSum) < eps) exit
            lastSum = currentSum
            n = n * 2
        end do
        integral_result = currentSum
    end function integral

end module integral_mod
