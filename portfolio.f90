module portfolio_mod
use kind_mod,  only: dp
use util_mod,  only: default, write_format
use qsort_mod, only: indexx
implicit none
private
public :: print_portfolio
contains
subroutine print_portfolio(sym,shares,prices,dollar_unit,min_abs_dollar_pos,max_print, &
                           sort_criterion,fmt_header,fmt_trailer)
! print a portfolio, optionally in order from largest to smallest positions
character (len=*), intent(in)           :: sym(:)    ! (npos) stock symbols
real(kind=dp)    , intent(in)           :: shares(:) ! (npos) share quantities
real(kind=dp)    , intent(in)           :: prices(:) ! (npos) share prices
real(kind=dp)    , intent(in), optional :: dollar_unit ! dollar unit in which dollar positions are printed -- 1000 to print $K
real(kind=dp)    , intent(in), optional :: min_abs_dollar_pos ! minimum dollar position to be printed
character (len=*), intent(in), optional :: fmt_header,fmt_trailer
character (len=*), intent(in), optional :: sort_criterion ! "none" or "abs_dollars"
integer          , intent(in), optional :: max_print ! max # of positions to print
integer                                 :: ipos,j,npos,indx(size(sym))
real(kind=dp)                           :: pct,cumul_pct, dollars(size(sym)),sum_dollars
character (len=100)                     :: sort_criterion_
call write_format(fmt_header)
write (*,"(*(a12))") "symbol","shares","prices","$","%","%cumul"
npos = size(sym)
if (size(shares) /= npos .or. size(prices) /= npos) &
   error stop "in print_portfolio, need size(sym) == size(shares) == size(prices)"
cumul_pct = 0.0_dp
dollars = shares*prices
sum_dollars = sum(dollars)
if (present(sort_criterion)) then
   sort_criterion_ = sort_criterion
else
   sort_criterion_ = "abs_dollars"
end if
if (sort_criterion_ == "abs_dollars") then
   indx = indexx(-abs(dollars))
else
   indx = [(j,j=1,npos)]
end if
do j=1,npos
   if (present(max_print)) then
      if (j > max_print) exit
   end if
   ipos = indx(j)
   pct = 100*dollars(ipos)/sum_dollars   
   cumul_pct = cumul_pct + pct
   if (abs(dollars(ipos)) >= default(0.0_dp,min_abs_dollar_pos)) write (*,"(a12,*(f12.2))") trim(sym(ipos)), &
       shares(ipos),prices(ipos),dollars(ipos)/default(1.0_dp,dollar_unit),pct,cumul_pct
end do
write (*,"(3a12,*(f12.2))") "<TOTAL>","---","---",sum(dollars)/default(1.0_dp,dollar_unit),100.0_dp
call write_format(fmt_trailer)
end subroutine print_portfolio
end module portfolio_mod