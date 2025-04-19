module loop_mod
  implicit none
  private
  public :: next,partitions
contains
  subroutine partitions(nsum,ni,imin,imax,ipart)
    ! return in ipart(:,:) integer partitions of nsum -- sets of integers within bounds imin(:), imax(:) that sum to nsum
    integer, intent(in)               :: nsum,ni,imin(:),imax(:)
    integer, intent(out), allocatable :: ipart(:,:)
    integer                           :: max_part_,itry(ni),npart
    integer, allocatable              :: jpart(:,:)
    logical,parameter                 :: print_part = .false.
    max_part_ = 1000000
    allocate (jpart(max_part_,ni))
    itry = imin
    npart = 0
    do
       itry = next(itry,imax)
       if (sum(itry) == nsum) then
          npart = npart + 1
          if (print_part) print*,itry
          jpart(npart,:) = itry
       end if
       if (all(itry == imax)) exit
    end do
    allocate (ipart(npart,ni))
    ipart = jpart(:npart,:)
  end subroutine partitions
!  
  function next(ii,imax,ibase) result(inext)
    ! return the next integer vector from an N-dimensional loop with
    ! lower bound ibase and upper bound imax(:) 
    integer, intent(in) :: ii(:),imax(:)
    integer, intent(in), optional :: ibase
    integer             :: inext(size(ii))
    integer             :: idim,ndim,ibase_
    character (len=*), parameter :: msg = "in next, "
    if (present(ibase)) then
       ibase_ = ibase
    else
       ibase_ = 0
    end if
    ndim = size(ii)
    if (size(imax) /= ndim) then
       inext = -1
       print*,msg,"size(ii), size(imax) =",ndim,size(imax)," must be equal, RETURNING"
       return
    else if (any(imax < ibase_)) then
       inext = -2
       print*,msg,"imax = ",imax," need all(imax >=",ibase_,"), RETURNING"
       return
    end if
    if (ii(ndim) < imax(ndim)) then
       inext = ii
       inext(ndim) = inext(ndim) + 1
       return
    end if
    do idim=ndim-1,1,-1
       if (ii(idim) < imax(idim)) then
          inext = ii
          inext(idim) = inext(idim) + 1
          inext(idim+1:) = ibase_
          return
       end if
    end do
    inext = ibase_
  end function next  
end module loop_mod
