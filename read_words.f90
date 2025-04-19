module read_words_mod
use util_mod, only: set_optional
implicit none
private
public :: read_words
contains
subroutine read_words(iu,words,nlines,nwords_write)
integer          , intent(in)           :: iu
character (len=*), intent(in out)       :: words(:)
integer          , intent(in), optional :: nlines,nwords_write
integer                                 :: i,ierr,iline,nlines_,nwords,nwords_write_,nwrite_progress_period_
character (len=10000)                   :: text
nwrite_progress_period_ = 1000
call set_optional(nwords_write_,0,nwords_write)
call set_optional(nlines_,10000000,nlines)
nwords = size(words)
do iline=1,nlines_
   words = ""
   read (iu,"(a)",iostat=ierr) text
   if (ierr /= 0) then
      write (*,*) "could not read line ",iline
      exit
   end if
   read (text,*,iostat=ierr) words
   if (nwords_write_ > 0) then
      do i=1,min(nwords,nwords_write_)
         if (words(i) /= "") write (*,"('line#',i0,1x,i4,1x,a)") iline,i,"'"//trim(words(i))//"'"
      end do
      write (*,*)
   else
      if (modulo(iline,nwrite_progress_period_) == 0) write (*,"('line ',i0)") iline
   end if
end do
end subroutine read_words
end module read_words_mod
