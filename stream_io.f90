module stream_io_mod
use kind_mod, only: dp
implicit none
private
public :: write_vec,read_vec,write_array_real,read_array_real,read_matrix_stream, &
          write_matrix_stream
interface write_vec
   module procedure write_vec_real
end interface write_vec
interface read_vec
   module procedure read_vec_real
end interface read_vec
contains
subroutine write_vec_real(file,vec)
! write 1D array of reals to file as unformatted stream, preceded by size
character (len=*), intent(in)  :: file
real(kind=dp)    , intent(in)  :: vec(:)
integer                        :: iu
open (newunit=iu,file=file,action="write",status="replace",access="stream",form="unformatted")
write (iu) size(vec),vec
close (iu)
end subroutine write_vec_real
!
subroutine read_vec_real(file,vec)
! read 1D array of reals from file as unformatted stream
character (len=*), intent(in)               :: file
real(kind=dp)    , intent(out), allocatable :: vec(:)
integer                                     :: iu,n
open (newunit=iu,file=file,action="read",status="old",access="stream",form="unformatted")
read (iu) n
allocate (vec(n))
read (iu) vec
end subroutine read_vec_real
!
subroutine write_array_real(file,x)
! write array of reals to file as unformatted stream, preceded by shape
character (len=*), intent(in)  :: file
real(kind=dp)    , intent(in)  :: x(..)
integer                        :: iu
open (newunit=iu,file=file,action="write",status="replace",access="stream",form="unformatted")
select rank(x)
   rank(0) ; write (iu) shape(x),x
   rank(1) ; write (iu) shape(x),x
   rank(2) ; write (iu) shape(x),x
   rank(3) ; write (iu) shape(x),x
end select
close (iu)
end subroutine write_array_real
!
subroutine read_array_real(file,x)
! read array of reals from file as unformatted stream
character (len=*), intent(in)               :: file
real(kind=dp)    , intent(out), allocatable :: x(..)
integer                                     :: iu,n1,n2,n3
open (newunit=iu,file=file,action="read",status="replace",access="stream",form="unformatted")
select rank(x)
   rank(0) ; allocate (x)                                ; read (iu) x
   rank(1) ; read (iu) n1       ; allocate (x(n1))       ; read (iu) x
   rank(2) ; read (iu) n1,n2    ; allocate (x(n1,n2))    ; read (iu) x
   rank(3) ; read (iu) n1,n2,n3 ; allocate (x(n1,n2,n3)) ; read (iu) x
end select
close (iu)
end subroutine read_array_real
!
subroutine read_matrix_stream(file,col_labels,row_labels,col_int,row_int,x,n1,n2)
character (len=*), intent(in)                         :: file
character (len=:), intent(out), allocatable, optional :: col_labels(:),row_labels(:)
integer          , intent(out), allocatable, optional :: col_int(:),row_int(:)
real(kind=dp)    , intent(out), allocatable, optional :: x(:,:)
integer          , intent(out),              optional :: n1,n2
integer                                               :: iu,n1_,n2_,nlen
open (newunit=iu,file=file,action="read",status="old",access="stream",form="unformatted")
read (iu) n1_,n2_
if (present(col_labels) .or. present(row_labels)) read (iu) nlen
if (present(col_labels)) allocate (character (len=nlen) :: col_labels(n2_))
if (present(row_labels)) allocate (character (len=nlen) :: row_labels(n1_))
if (present(col_int))    allocate (col_int(n2_))
if (present(row_int))    allocate (row_int(n1_))
if (present(x))          allocate (x(n1_,n2_))
if (present(col_labels)) read (iu) col_labels
if (present(col_int))    read (iu) col_int
if (present(row_labels)) read (iu) row_labels
if (present(row_int))    read (iu) row_int
if (present(x))          read (iu) x
close (iu)
if (present(n1)) n1 = n1_
if (present(n2)) n2 = n2_
end subroutine read_matrix_stream
!
subroutine write_matrix_stream(file,x,col_labels,row_labels,col_int,row_int)
character (len=*), intent(in)           :: file
real(kind=dp)    , intent(in)           :: x(:,:)
character (len=*), intent(in), optional :: col_labels(:),row_labels(:)
integer          , intent(in), optional :: col_int(:),row_int(:)
integer                                 :: iu,nlen
open (newunit=iu,file=file,action="write",status="replace",access="stream",form="unformatted")
write (iu) shape(x)
nlen = -1
if (present(col_labels)) then
   nlen = len(col_labels)
else if (present(row_labels)) then
   nlen = len(row_labels)
end if
if (nlen >= 0) write (iu) nlen
if (present(col_labels)) write (iu) col_labels
if (present(row_labels)) write (iu) row_labels
if (present(col_int))    write (iu) col_int
if (present(row_int))    write (iu) row_int
write (iu) x
close (iu)
end subroutine write_matrix_stream
end module stream_io_mod
