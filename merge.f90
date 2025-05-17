module m_mod
implicit none
private
public :: mymerge
contains
elemental function mymerge(x, y, tf) result(z)
! merge function that returns a character variable whose LEN is the
! larger of the LENs of the two arguments
character (len=*), intent(in) :: x, y
logical          , intent(in) :: tf
character (len=max(len(x), len(y))) :: z
if (tf) then
   z = x
else
   z = y
end if
end function mymerge
end module m_mod

program main
! test mymerge
use m_mod, only: mymerge
implicit none
character (len=*), parameter :: fmt_c = "(*(1x,a))"
print fmt_c,mymerge("yes", "no", .true.)
print fmt_c,mymerge("yes", "no", [.true., .false.])
print fmt_c,mymerge(["yes ", "sure", "ok  "], "no", .true.)
print fmt_c,mymerge(["yes ", "sure", "ok  "], "no", [.true.,.false.,.true.])
end program main
