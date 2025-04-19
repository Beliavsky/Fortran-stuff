module factors_mod
use kind_mod , only: dp
use qsort_mod, only: indexx
use util_mod , only: first_false
implicit none
private
public :: compress,print_sums_by_factor
contains
pure subroutine compress(words,factors,values)
character (len=*)         , intent(in)               :: words(:)   ! (n) 
integer                   , intent(out)              :: factors(:) ! (n) integer values corresponding to words(:)
character (len=len(words)), allocatable, intent(out) :: values(:)  ! unique values of words(:)
integer :: i,n,nfac,imatch
n = size(words)
if (size(factors) /= n) error stop "in compress, need size(words) == size(factors)"
if (n < 1) then
   allocate (values(0))
   return
end if
allocate (values(n))
values(1) = words(1)
factors(1) = 1
nfac = 1
do i=2,n
   imatch = findloc(values(:nfac),words(i),dim=1)
   if (imatch == 0) then
      nfac = nfac + 1
      factors(i) = nfac
      values(nfac) = words(i)
   else
      factors(i) = imatch
   end if
end do
values = values(:nfac)
end subroutine compress
!
subroutine print_sums_by_factor(dollar_positions,sym_unique,ipos_sym,dollars_cat,symbol_fields,symbol_categories,dollar_unit)
! print sums of dollar positions by various categories
real(kind=dp)    , intent(in)  :: dollar_positions(:)    ! (nunique)
character (len=*), intent(in)  :: sym_unique(:)          ! (nunique)
integer          , intent(in)  :: ipos_sym(:)            ! (nunique)
real(kind=dp)    , intent(in)  :: dollars_cat(:)         ! (nsym_cat)
character (len=*), intent(in)  :: symbol_fields(:)       ! (nsym_fields)
character (len=*), intent(in)  :: symbol_categories(:,:) ! (nsym_cat,nsym_fields+1)
real(kind=dp)    , intent(in)  :: dollar_unit
character (len=len(symbol_categories)), allocatable :: fields_unique(:)
integer                        :: i,ifield,isym,ifac,nfac,j,k,nsym_fields,nunique,nsym_cat,ierr
integer          , allocatable :: factors(:),indx_dollar(:),indx_fac(:)
real(kind=dp)                  :: sum_dollars
logical          , allocatable :: mask_fac(:)
logical          , parameter   :: debug = .false.
real(kind=dp)    , allocatable :: dollars_fac(:)
logical                        :: sort_by_dollars
sort_by_dollars = .true.
if (debug) then
   print*,"entered print_sums_by_factor, size(dollar_positions)=",size(dollar_positions)," size(sym_unique)=",size(sym_unique) !! debug
   print*,"size(ipos_sym)=",size(ipos_sym)," size(dollars_cat)=",size(dollars_cat)," size(symbol_fields)=",size(symbol_fields) !! debug
   print*,"shape(symbol_categories)=",shape(symbol_categories) !! debug
end if
ierr = first_false([size(dollar_positions)==size(sym_unique),size(ipos_sym)==size(sym_unique), &
                    size(symbol_categories,1)==size(dollars_cat),size(symbol_categories,2)==size(symbol_fields)+1])
if (ierr /= 0) then
   print*,"in factors_mod::print_sums_by_factor, ierr =",ierr," STOPPING"
   error stop
end if
nunique     = size(sym_unique)
nsym_fields = size(symbol_fields)
nsym_cat    = size(symbol_categories,dim=1)
sum_dollars = sum(dollar_positions)
indx_dollar = indexx(-dollar_positions) ! allocation
allocate (factors(nsym_cat))
do ifield=1,nsym_fields
   call compress(symbol_categories(:,ifield+1),factors,fields_unique) ! out: factors(:), field_unique(:)
   nfac = size(fields_unique)
   allocate (dollars_fac(nfac))
   write (*,"(/,i6,*(1x,a15))") nfac,trim(symbol_fields(ifield)) // "(s)"
   write (*,"(a10,a6,a10,a8)") "category","#","$","%"
   do ifac=1,nfac
      mask_fac = factors==ifac .and. dollars_cat /= 0.0_dp
      dollars_fac(ifac) = sum(dollars_cat,mask_fac)
      if (.not. sort_by_dollars) &
         write (*,"(a10,i6,f10.1,f8.2)") trim(fields_unique(ifac)),count(mask_fac), &
                                 dollars_fac(ifac)/dollar_unit,100*dollars_fac(ifac)/sum_dollars
   end do
   if (sort_by_dollars) then
      indx_fac = indexx(-dollars_fac)
      do j=1,nfac
         ifac = indx_fac(j)
         mask_fac = factors==ifac .and. dollars_cat /= 0.0_dp
         write (*,"(a10,i6,f10.1,f8.2)") trim(fields_unique(ifac)),count(mask_fac), &
                                         dollars_fac(ifac)/dollar_unit,100*dollars_fac(ifac)/sum_dollars
      end do
   end if
   write (*,"(/,a)") "by " // trim(symbol_fields(ifield))
   do ifac=1,nfac
      mask_fac = factors==ifac .and. dollars_cat /= 0.0_dp
      do i=1,nunique
         j = indx_dollar(i)
         isym = ipos_sym(j)
         if (mask_fac(isym)) &
            write (*,"(a15,f15.1,*(a15))") trim(sym_unique(j)),dollar_positions(j)/dollar_unit, &
                                          (trim(symbol_categories(isym,k)),k=2,nsym_fields+1)
      end do
   end do
   deallocate (dollars_fac)
end do
end subroutine print_sums_by_factor
!
end module factors_mod
