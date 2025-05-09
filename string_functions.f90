! -----------------------------------------------
module string_functions  ! by david frank  dave_frank@hotmail.com
implicit none            ! http://home.earthlink.net/~dave_gemini/strings.f90

! copy (generic) char array to string or string to char array
! clen           returns same as len      unless last non-blank char = null
! clen_trim      returns same as len_trim    "              "
! ctrim          returns same as trim        "              "
! count_items    in string that are blank or comma separated
! reduce_blanks  in string to 1 blank between items, last char not blank
! replace_text   in all occurances in string with replacement string
! spack          pack string's chars == extract string's chars
! tally          occurances in string of text arg
! translate      text arg via indexed code table
! upper/lower    case the text arg

interface copy    ! generic
   module procedure copy_a2s, copy_s2a
end interface copy
public :: replace_text,replace_text_vec
contains
! ------------------------
pure function copy_a2s(a)  result (s)    ! copy char array to string
character,intent(in) :: a(:)
character(size(a)) :: s
integer :: i
do i = 1,size(a)
   s(i:i) = a(i)
end do
end function copy_a2s

! ------------------------
pure function copy_s2a(s)  result (a)   ! copy s(1:clen(s)) to char array
character(*),intent(in) :: s
character :: a(len(s))
integer :: i
do i = 1,len(s)
   a(i) = s(i:i)
end do
end function copy_s2a

! ------------------------
pure integer function clen(s)      ! returns same result as len unless:
character(*),intent(in) :: s       ! last non-blank char is null
integer :: i
clen = len(s)
i = len_trim(s)
if (s(i:i) == char(0)) clen = i-1  ! len of c string
end function clen

! ------------------------
pure integer function clen_trim(s) ! returns same result as len_trim unless:
character(*),intent(in) :: s       ! last char non-blank is null, if true:
integer :: i                       ! then len of c string is returned, note:
                                   ! ctrim is only user of this function
i = len_trim(s) ; clen_trim = i
if (s(i:i) == char(0)) clen_trim = clen(s)   ! len of c string
end function clen_trim

! ----------------
function ctrim(s1)  result(s2)     ! returns same result as trim unless:
character(*),intent(in)  :: s1     ! last non-blank char is null in which
character(clen_trim(s1)) :: s2     ! case trailing blanks prior to null
s2 = s1                            ! are output
end function ctrim

! --------------------
integer function count_items(s1)  ! in string or c string that are blank or comma separated
character(*), intent(in) :: s1
character(clen(s1)) :: s
integer :: i, k

s = s1                            ! remove possible last char null
k = 0  ; if (s /= ' ') k = 1      ! string has at least 1 item
do i = 1,len_trim(s)-1
   if (s(i:i) /= ' '.and.s(i:i) /= ',' &
                    .and.s(i+1:i+1) == ' '.or.s(i+1:i+1) == ',') k = k+1
end do
count_items = k
end function count_items

! --------------------
function reduce_blanks(s)  result (outs)
character(*), intent(in)      :: s
character(len_trim(s)) :: outs
integer           :: i, k, n

n = 0  ; k = len_trim(s)          ! k=index last non-blank (may be null)
do i = 1,k-1                      ! dont process last char yet
   n = n+1 ; outs(n:n) = s(i:i)
   if (s(i:i+1) == '  ') n = n-1  ! backup/discard consecutive output blank
end do
n = n+1  ; outs(n:n)  = s(k:k)    ! last non-blank char output (may be null)
if (n < k) outs(n+1:) = ' '       ! pad trailing blanks
end function reduce_blanks

! ------------------
function replace_text_vec(str,text,rep) result(out_str)
character(len=*), intent(in) :: str(:),text(:),rep(:)
character(len(str)+100)      :: out_str(size(str))     ! provide outs with extra 100 char len
integer                      :: istr,irep,nrep
nrep = size(rep)
if (size(text) /= nrep) then
   write (*,*) "in replace_text_vec, size(text), size(rep) =",size(text),nrep," must be equal, STOPPING"
   stop
end if
do istr=1,size(str)
   out_str(istr) = str(istr)
   do irep=1,nrep
      out_str(istr) = replace_text(out_str(istr),text(irep),rep(irep))
   end do
end do
end function replace_text_vec

function replace_text(str,text,rep)  result(out_str)
character(*), intent(in) :: str,text,rep
character(len(str)+100)  :: out_str     ! provide outs with extra 100 char len
integer                  :: i, nt, nr
out_str = str
nt = len_trim(text)
nr = len_trim(rep)
do
   i = index(out_str,text(:nt)) ; if (i == 0) exit
   out_str = out_str(:i-1) // rep(:nr) // out_str(i+nt:)
end do
end function replace_text

! ---------------------------------
function spack (s,ex)  result (outs)
character(*), intent(in) :: s,ex
character(len(s)) :: outs
character :: aex(len(ex))   ! array of ex chars to extract
integer   :: i, n

n = 0  ;  aex = copy(ex)
do i = 1,len(s)
   if (.not.any(s(i:i) == aex)) cycle   ! dont pack char
   n = n+1 ; outs(n:n) = s(i:i)
end do
outs(n+1:) = ' '     ! pad with trailing blanks
end function spack

! --------------------
integer function tally (s,text)
character(*), intent(in) :: s, text
integer :: i, nt

tally = 0 ; nt = len_trim(text)
do i = 1,len(s)-nt+1
   if (s(i:i+nt-1) == text(:nt)) tally = tally+1
end do
end function tally

! ---------------------------------
function translate(s1,codes)  result (s2)
character(*), intent(in) :: s1, codes(2)
character(len(s1)) :: s2
character          :: ch
integer            :: i, j

do i = 1,len(s1)
   ch = s1(i:i)
   j = index(codes(1),ch) ; if (j > 0) ch = codes(2)(j:j)
   s2(i:i) = ch
end do
end function translate

! ---------------------------------
function upper(s1)  result (s2)
character(*), intent(in) :: s1
character(len(s1)) :: s2
character          :: ch
integer,parameter  :: duc = ichar('a') - ichar('a')
integer            :: i

do i = 1,len(s1)
   ch = s1(i:i)
   if (ch >= 'a'.and.ch <= 'z') ch = char(ichar(ch)+duc)
   s2(i:i) = ch
end do
end function upper

! ---------------------------------
function lower(s1)  result (s2)
character(*), intent(in)       :: s1
character(len(s1)) :: s2
character          :: ch
integer,parameter  :: duc = ichar('a') - ichar('a')
integer            :: i

do i = 1,len(s1)
   ch = s1(i:i)
   if (ch >= 'a'.and.ch <= 'z') ch = char(ichar(ch)-duc)
   s2(i:i) = ch
end do
end function lower

end module string_functions