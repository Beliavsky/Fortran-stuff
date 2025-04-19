module get_unit_mod
    implicit none
    private
    public :: get_unit_open_file, read_lines_file
    
contains
    subroutine get_unit_open_file(filename, unit, action)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: unit
        character(len=*), intent(in) :: action  ! "r" for read, "w" for write
        
        logical :: is_open
        integer :: ios
        
        ! Find an available unit number
        unit = 10  ! Starting point
        do
            inquire(unit=unit, opened=is_open)
            if (.not. is_open) exit
            unit = unit + 1
            if (unit > 999) stop "No available unit numbers"
        end do
        
        ! Open the file with specified action
        if (action == "r") then
            open(unit=unit, file=trim(filename), status='old', action='read', iostat=ios)
        else if (action == "w") then
            open(unit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
        else
            stop "Invalid action in get_unit_open_file: use 'r' or 'w'"
        end if
        
        if (ios /= 0) then
            print *, "Error opening file: ", trim(filename)
            stop
        end if
    end subroutine get_unit_open_file
    
    subroutine read_lines_file(lines, filename, nmax, close_file)
        character(len=*), allocatable, intent(out) :: lines(:)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: nmax
        logical, intent(in) :: close_file
        
        integer :: unit, ios, nlines
        character(len=1000) :: line
        integer :: i
        
        ! Open file for reading
        call get_unit_open_file(filename, unit, "r")
        
        ! Count lines first
        nlines = 0
        do
            read(unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            nlines = nlines + 1
            if (nlines > nmax) stop "File exceeds maximum lines"
        end do
        
        ! Allocate array and reread file
        allocate(lines(nlines))
        rewind(unit)
        
        do i = 1, nlines
            read(unit, '(a)', iostat=ios) lines(i)
            if (ios /= 0) stop "Error reading lines"
        end do
        
        ! Close file if requested
        if (close_file) then
            close(unit)
        end if
    end subroutine read_lines_file
    
end module get_unit_mod
