module read_xy_mod
use kind_mod, only: dp
use util_mod, only: default
implicit none
private
public :: read_xy
contains
subroutine read_xy(iu,xx,yy,nread,xfirst,comment_char)
integer      , intent(in)  :: iu
real(kind=dp), intent(out) :: xx(:),yy(:)
integer      , intent(out) :: nread
logical      , intent(in), optional :: xfirst
character (len=*), intent(in), optional :: comment_char
integer                    :: iostatus,n
character (len=100)        :: text
logical                    :: xfirst_
xfirst_ = default(.true.,xfirst)
n = size(xx)
nread = 1
do
  read (iu,"(a)",iostat=iostatus) text
  if (iostatus > 0) cycle
  if (iostatus < 0) then
     nread = nread - 1
     exit
  end if
  if (present(comment_char)) then
      if (text(1:1) == comment_char) cycle
  end if
  if (xfirst_) then
    read (text,*,iostat=iostatus) xx(nread),yy(nread)
  else
    read (text,*,iostat=iostatus) yy(nread),xx(nread)
  end if
  if (iostatus /= 0) then
     write (*,*) "in read_xy, could not read two numbers from '" // trim(text) // "', STOPPING"
     stop
  end if
  if (nread == n) exit
  nread = nread + 1
end do
end subroutine read_xy
end module read_xy_mod
