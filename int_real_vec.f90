module int_real_vec_mod
use kind_mod, only: dp,i64
implicit none
private
public :: ikind, int_real_vec
integer, parameter :: ikind = i64
type, public :: int_real_vec
   character (len=100)               :: name
   integer (kind=ikind), allocatable :: ivec(:) ! (n)
   real(kind=dp)       , allocatable :: xx(:)   ! (n)
end type int_real_vec
contains
end module int_real_vec_mod
