module kind_mod
implicit none
private
public :: dp, i64
integer, parameter :: dp = kind(1.0d0), i64 = selected_int_kind(18)
end module kind_mod
