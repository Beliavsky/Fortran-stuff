module spend_mod
use kind_mod  , only: dp
use random_mod, only: random_normal
use util_mod  , only: istdout,default,write_format
implicit none
contains
subroutine spending_sim(spend_frac_init,ann_mean,ann_vol,nyears,npaths,nsim,outu)
integer, intent(in) :: nyears,npaths,nsim
! integer      , parameter :: nyears = 30, npaths = 100000, nsim = 3
real(kind=dp)            :: xx(0:nyears),xret(nyears),xlast(npaths),starting_wealth_year
real(kind=dp), intent(in) :: spend_frac_init,ann_mean,ann_vol
integer      , intent(in), optional :: outu
real(kind=dp), parameter :: x0 = 1000.0_dp
real(kind=dp)            :: ann_spend
integer                  :: i,ipath,nyears_lasted(npaths),num_lasted(0:nyears),isim,outu_
logical      , parameter :: print_terminal_wealth = .false.
character (len=*), parameter :: fmt_cr = "(a25,':',1x,f10.4)",fmt_ci = "(a25,':',1x,i10)"
outu_ = default(istdout,outu)
ann_spend = spend_frac_init*x0
call write_format("()",outu_)
write (outu_,fmt_cr) "initial spending fraction",spend_frac_init,"annual return",ann_mean,"annual vol",ann_vol
write (outu_,fmt_ci) "#paths",npaths,"#sim",nsim
do_sim: do isim=1,nsim
   do ipath=1,npaths
      xx(0) = x0
      xret = ann_mean + ann_vol*random_normal(nyears)
      nyears_lasted(ipath) = nyears
      do i=1,nyears
         starting_wealth_year = xx(i-1) - ann_spend
         if (starting_wealth_year <= 0.0_dp) then
            nyears_lasted(ipath) = i-1
            xx(i:) = 0.0_dp
            exit
         end if
         xx(i) = starting_wealth_year * (1+xret(i))
      end do
      xlast(ipath) = xx(nyears)
   end do
   if (print_terminal_wealth) write (outu_,"(100(f10.4,1x))") xlast
   forall (i=0:nyears) num_lasted(i) = count(nyears_lasted == i)
   write (outu_,"(/,100a8)") "year","BK","survive"
   do i=0,nyears
      write (outu_,"(i8,2f8.4)") i,[num_lasted(i),sum(num_lasted(i:))]/dble(npaths)
   end do
end do do_sim
end subroutine spending_sim
end module spend_mod
