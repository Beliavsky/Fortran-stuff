module date_mdy_mod
    implicit none
    private
    public :: mmddyy_to_mdy
    
contains
    function mmddyy_to_mdy(mmddyy) result(mdy)
        integer, intent(in) :: mmddyy(:)  ! Input array of dates in MMDDYY format
        integer, allocatable :: mdy(:,:)  ! Output array: [month, day, year] for each date
        
        integer :: n, i
        integer :: month, day, year
        
        ! Get size of input array
        n = size(mmddyy)
        
        ! Allocate output array: 3 components (month, day, year) for each input date
        allocate(mdy(3,n))
        
        ! Process each date
        do i = 1, n
            ! Extract components from MMDDYY format
            month = mmddyy(i) / 10000          ! First two digits
            day = mod(mmddyy(i) / 100, 100)    ! Middle two digits
            year = mod(mmddyy(i), 100)         ! Last two digits
            
            ! Basic validation
            if (month < 1 .or. month > 12) month = 0
            if (day < 1 .or. day > 31) day = 0
            if (year < 0) year = 0
            
            ! Assuming years are in 20th/21st century (can be adjusted)
            if (year <= 99) then
                if (year >= 50) then
                    year = year + 1900  ! 1950-1999
                else
                    year = year + 2000  ! 2000-2049
                endif
            endif
            
            ! Store results
            mdy(1,i) = month
            mdy(2,i) = day
            mdy(3,i) = year
        end do
        
    end function mmddyy_to_mdy
    
end module date_mdy_mod
