module cobyla2_mod

! the fortran 77 version of this code was by michael powell
! (m.j.d.powell @ damtp.cam.ac.uk)

! this fortran 90 version by alan miller
! this version is in a subset of f90 and is compatible with version 3.0a
! of lahey's free elf90 compiler.
! alan.miller @ vic.cmis.csiro.au
! latest revision - 28 june 2000
use calcfc_mod, only: calcfc
implicit none
private
public :: cobyla,cobyla_dt,set,display
integer, parameter, public :: dp = selected_real_kind(14, 60)
type, public :: cobyla_param
   real(kind=dp) :: rhobeg, rhoend
   integer       :: iprint
   integer       :: maxfun
end type cobyla_param
integer, parameter :: istdout = 6
contains
subroutine display(par,outu,fmt_header,fmt_trailer)
type(cobyla_param), intent(in)           :: par
integer           , intent(in), optional :: outu
character (len=*) , intent(in), optional :: fmt_header,fmt_trailer
integer                                  :: outu_
character  (len=*), parameter            :: fmt_cr="(a20,':',1x,f12.6)",fmt_ci="(a20,':',1x,i0)"
if (present(outu)) then
   outu_ = outu
else
   outu_ = istdout
end if
if (present(fmt_header)) write (outu_,fmt_header)
write (outu_,fmt_cr) "rhobeg",par%rhobeg,"rhoend",par%rhoend
write (outu_,fmt_ci) "maxfun",par%maxfun
call write_format(fmt_trailer)
end subroutine display
!
subroutine write_format(format_str,iunit)
! write format_str to unit iunit if it is present and not blank
! otherwise, do nothing
character (len=*), intent(in), optional :: format_str
integer          , intent(in), optional :: iunit
integer                                 :: iu
if (.not. present(format_str)) return
if (present(iunit)) then
   iu = iunit
else
   iu = istdout
end if
! iu = default(istdout,iunit)
if (format_str /= "") write (iu,format_str)
end subroutine write_format
!
elemental subroutine set(par,rhobeg,rhoend,iprint,maxfun)
type(cobyla_param) , intent(in out)        :: par
real(kind=dp)      , intent(in), optional  :: rhobeg ! size of reasonable initial changes
real(kind=dp)      , intent(in), optional  :: rhoend ! required accuracy in variables
integer            , intent(in), optional  :: iprint
integer            , intent(in), optional  :: maxfun
if (present(rhobeg)) par%rhobeg = rhobeg
if (present(rhoend)) par%rhoend = rhoend
if (present(iprint)) par%iprint = iprint
if (present(maxfun)) par%maxfun = maxfun
end subroutine set
!
subroutine cobyla_dt(ncon, x, par, fmin, nfunc_eval)
integer           , intent(in)            :: ncon   ! # of constraints
real (dp)         , intent(in out)        :: x(:)   ! input: intitial guess; output: location of minimum
type(cobyla_param), intent(in)            :: par
integer           , intent(out), optional :: nfunc_eval
integer                                   :: maxfun
real(kind=dp)     , intent(out), optional :: fmin
maxfun = par%maxfun
call cobyla(size(x), ncon, x, par%rhobeg, par%rhoend, par%iprint, maxfun, fmin=fmin)
if (present(nfunc_eval)) nfunc_eval = maxfun
end subroutine cobyla_dt

!------------------------------------------------------------------------

subroutine cobyla(n, m, x, rhobeg, rhoend, iprint, maxfun, fmin, con)
integer, intent(in)        :: n      ! size of x(:)
integer, intent(in)        :: m      ! # of constraints
real (dp), intent(in out)  :: x(:)   ! input: intitial guess; output: location of minimum
real (dp), intent(in)      :: rhobeg ! size of reasonable initial changes
real (dp), intent(in)      :: rhoend ! required accuracy in variables
integer, intent(in)        :: iprint
integer, intent(in out)    :: maxfun
real(kind=dp), intent(out), optional :: fmin
real(kind=dp), intent(out), optional :: con(:)
!  this subroutine minimizes an objective function f(x) subject to m
!  inequality constraints on x, where x is a vector of variables that has
!  n components.  the algorithm employs linear approximations to the
!  objective and constraint functions, the approximations being formed by
!  linear interpolation at n+1 points in the space of the variables.
!  we regard these interpolation points as vertices of a simplex.  the
!  parameter rho controls the size of the simplex and it is reduced
!  automatically from rhobeg to rhoend.  for each rho the subroutine tries
!  to achieve a good vector of variables for the current size, and then
!  rho is reduced until the value rhoend is reached.  therefore rhobeg and
!  rhoend should be set to reasonable initial changes to and the required
!  accuracy in the variables respectively, but this accuracy should be
!  viewed as a subject for experimentation because it is not guaranteed.
!  the subroutine has an advantage over many of its competitors, however,
!  which is that it treats each constraint individually when calculating
!  a change to the variables, instead of lumping the constraints together
!  into a single penalty function.  the name of the subroutine is derived
!  from the phrase constrained optimization by linear approximations.

!  the user must set the values of n, m, rhobeg and rhoend, and must
!  provide an initial vector of variables in x.  further, the value of
!  iprint should be set to 0, 1, 2 or 3, which controls the amount of
!  printing during the calculation. specifically, there is no output if
!  iprint=0 and there is output only at the end of the calculation if
!  iprint=1.  otherwise each new value of rho and sigma is printed.
!  further, the vector of variables and some function information are
!  given either when rho is reduced or when each new value of f(x) is
!  computed in the cases iprint=2 or iprint=3 respectively. here sigma
!  is a penalty parameter, it being assumed that a change to x is an
!  improvement if it reduces the merit function
!             f(x)+sigma*max(0.0, - c1(x), - c2(x),..., - cm(x)),
!  where c1,c2,...,cm denote the constraint functions that should become
!  nonnegative eventually, at least to the precision of rhoend. in the
!  printed output the displayed term that is multiplied by sigma is
!  called maxcv, which stands for 'maximum constraint violation'.  the
!  argument maxfun is an integer variable that must be set by the user to a
!  limit on the number of calls of calcfc, the purpose of this routine being
!  given below.  the value of maxfun will be altered to the number of calls
!  of calcfc that are made.

!  in order to define the objective and constraint functions, we require
!  a subroutine that has the name and arguments
!             subroutine calcfc (n,m,x,f,con)
!             dimension x(:),con(:)  .
!  the values of n and m are fixed and have been defined already, while
!  x is now the current vector of variables. the subroutine should return
!  the objective and constraint functions at x in f and con(1),con(2),
!  ...,con(m).  note that we are trying to adjust x so that f(x) is as
!  small as possible subject to the constraint functions being nonnegative.

!  n.b. if the starting value for any x(i) is set to zero, that value will
!       not be changed.   this can be a useful feature in comparing
!       nested models.   if all the x(i)'s are set to zero, an error message
!       will result.

! local variable
real(kind=dp) :: fmin_,con_(m)
integer :: mpp

mpp = m + 2
! print*,"n, m, size(x)=",n,m,size(x)," x=",x !! debug
call cobylb (n, m, mpp, x, rhobeg, rhoend, iprint, maxfun)
if (present(fmin) .or. present(con)) call calcfc(n,m,x,fmin_,con_)
if (present(fmin)) fmin = fmin_
if (present(con))   con = con_
end subroutine cobyla

!------------------------------------------------------------------------------

subroutine cobylb (n, m, mpp, x, rhobeg, rhoend, iprint, maxfun)

! n.b. arguments con, sim, simi, datmat, a, vsig, veta, sigbar, dx, w & iact
!   have been removed.

integer, intent(in)        :: n
integer, intent(in)        :: m
integer, intent(in)        :: mpp
real (dp), intent(in out)  :: x(:)
real (dp), intent(in)      :: rhobeg
real (dp), intent(in)      :: rhoend
integer, intent(in)        :: iprint
integer, intent(out)       :: maxfun ! this variable is used in a comparison in line 171 without being previously set
                                     ! i think maxfund should have intent(in out)

!  set the initial values of some parameters. the last column of sim holds
!  the optimal vertex of the current simplex, and the preceding n columns
!  hold the displacements from the optimal vertex to the other vertices.
!  further, simi holds the inverse of the matrix that is contained in the
!  first n columns of sim.

! local variables

real (dp) :: con(mpp), sim(n,n+1), simi(n,n), datmat(mpp,n+1), a(n,m+1),      &
             vsig(n), veta(n), sigbar(n), dx(n), w(n)
real (dp) :: alpha, barmu, beta, cmin, cmax, cvmaxm, cvmaxp, delta, denom,    &
             dxsign, edgmax, error, f, gamma, pareta, parmu, parsig, phi,     &
             phimin, prerec, prerem, ratio, resmax, resnew, rho, temp, tempa, &
             total, trured, vmnew, vmold, weta, wsig
integer   :: i, ibrnch, iflag, ifull, iptem, iptemp, j, jdrop, k, l, mp,  &
             nbest, nfvals, np

iptem = min(n,5)
iptemp = iptem + 1
np = n + 1
mp = m + 1
alpha = 0.25_dp
beta = 2.1_dp
gamma = 0.5_dp
delta = 1.1_dp
rho = rhobeg
parmu = 0.0_dp
if (iprint >= 2) write(*, 10) rho
10 format (/'   the initial value of rho is', g13.6,   &
           '  and parmu is set to zero.')
nfvals = 0
temp = 1.0_dp/rho
do i=1,n
  sim(i,np) = x(i)
  do j=1,n
    sim(i,j) = 0.0_dp
    simi(i,j) = 0.0_dp
  end do
  sim(i,i) = rho
  simi(i,i) = temp
end do
jdrop = np
ibrnch = 0

!  make the next call of the user-supplied subroutine calcfc. these
!  instructions are also used for calling calcfc during the iterations of
!  the algorithm.

40 if (nfvals >= maxfun .and. nfvals > 0) then
  if (iprint >= 1) write(*, 50)
  50 format (/'   return from subroutine cobyla because the ',  &
             'maxfun limit has been reached.')
  go to 600
end if
nfvals = nfvals + 1
call calcfc (n, m, x, f, con)
resmax = 0.0_dp
if (m > 0) then
  do k=1,m
    resmax = max(resmax, - con(k))
  end do
end if
if (nfvals == iprint-1 .or. iprint == 3) then
  write(*, 70) nfvals, f, resmax, x(1:iptem)
  70 format (/'   nfvals = ', i5, '   f = ', g13.6, '    maxcv = ',  &
             g13.6/ ('   x = ', 5g14.6))
  if (iptem < n) write(*, 80) x(iptemp:n)
  80 format (g19.6, g15.6)
end if
con(mp) = f
con(mpp) = resmax
if (ibrnch == 1) go to 440

!  set the recently calculated function values in a column of datmat. this
!  array has a column for each vertex of the current simplex, the entries of
!  each column being the values of the constraint functions (if any)
!  followed by the objective function and the greatest constraint violation
!  at the vertex.

do k=1,mpp
  datmat(k,jdrop) = con(k)
end do
if (nfvals > np) go to 130

!  exchange the new vertex of the initial simplex with the optimal vertex if
!  necessary. then, if the initial simplex is not complete, pick its next
!  vertex and calculate the function values there.

if (jdrop <= n) then
  if (datmat(mp,np) <= f) then
    x(jdrop) = sim(jdrop,np)
  else
    sim(jdrop,np) = x(jdrop)
    do k=1,mpp
      datmat(k,jdrop) = datmat(k,np)
      datmat(k,np) = con(k)
    end do
    do k=1,jdrop
      sim(jdrop,k) = -rho
      temp = -sum( simi(k:jdrop, k) )
      simi(jdrop,k) = temp
    end do
  end if
end if
if (nfvals <= n) then
  jdrop = nfvals
  x(jdrop) = x(jdrop) + rho
  go to 40
end if
130 ibrnch = 1

!  identify the optimal vertex of the current simplex.

140 phimin = datmat(mp,np) + parmu*datmat(mpp,np)
nbest = np
do j=1,n
  temp = datmat(mp,j) + parmu*datmat(mpp,j)
  if (temp < phimin) then
    nbest = j
    phimin = temp
  else if (temp == phimin .and. parmu == 0.0_dp) then
    if (datmat(mpp,j) < datmat(mpp,nbest)) nbest = j
  end if
end do

!  switch the best vertex into pole position if it is not there already,
!  and also update sim, simi and datmat.

if (nbest <= n) then
  do i=1,mpp
    temp = datmat(i,np)
    datmat(i,np) = datmat(i,nbest)
    datmat(i,nbest) = temp
  end do
  do i=1,n
    temp = sim(i,nbest)
    sim(i,nbest) = 0.0_dp
    sim(i,np) = sim(i,np) + temp
    tempa = 0.0_dp
    do k=1,n
      sim(i,k) = sim(i,k) - temp
      tempa = tempa - simi(k,i)
    end do
    simi(nbest,i) = tempa
  end do
end if

!  make an error return if sigi is a poor approximation to the inverse of
!  the leading n by n submatrix of sig.

error = 0.0_dp
do i=1,n
  do j=1,n
    temp = 0.0_dp
    if (i == j) temp = temp - 1.0_dp
    temp = temp + dot_product( simi(i,1:n), sim(1:n,j) )
    error = max(error, abs(temp))
  end do
end do
if (error > 0.1_dp) then
  if (iprint >= 1) write(*, 210)
  210 format (/'   return from subroutine cobyla because ',  &
              'rounding errors are becoming damaging.')
  go to 600
end if

!  calculate the coefficients of the linear approximations to the objective
!  and constraint functions, placing minus the objective function gradient
!  after the constraint gradients in the array a. the vector w is used for
!  working space.

do k=1,mp
  con(k) = -datmat(k,np)
  do j=1,n
    w(j) = datmat(k,j) + con(k)
  end do
  do i=1,n
    temp = dot_product( w(1:n), simi(1:n,i) )
    if (k == mp) temp = -temp
    a(i,k) = temp
  end do
end do

!  calculate the values of sigma and eta, and set iflag = 0 if the current
!  simplex is not acceptable.

iflag = 1
parsig = alpha*rho
pareta = beta*rho
do j=1,n
  wsig = sum( simi(j,1:n)**2 )
  weta = sum( sim(1:n,j)**2 )
  vsig(j) = 1.0_dp/sqrt(wsig)
  veta(j) = sqrt(weta)
  if (vsig(j) < parsig .or. veta(j) > pareta) iflag = 0
end do

!  if a new vertex is needed to improve acceptability, then decide which
!  vertex to drop from the simplex.

if (ibrnch == 1 .or. iflag == 1) go to 370
jdrop = 0
temp = pareta
do j=1,n
  if (veta(j) > temp) then
    jdrop = j
    temp = veta(j)
  end if
end do
if (jdrop == 0) then
  do j=1,n
    if (vsig(j) < temp) then
      jdrop = j
      temp = vsig(j)
    end if
  end do
end if

!  calculate the step to the new vertex and its sign.

temp = gamma*rho*vsig(jdrop)
dx(1:n) = temp*simi(jdrop,1:n)
cvmaxp = 0.0_dp
cvmaxm = 0.0_dp
do k=1,mp
  total = dot_product( a(1:n,k), dx(1:n) )
  if (k < mp) then
    temp = datmat(k,np)
    cvmaxp = max(cvmaxp, -total - temp)
    cvmaxm = max(cvmaxm, total - temp)
  end if
end do
dxsign = 1.0_dp
if (parmu*(cvmaxp - cvmaxm) > total + total) dxsign = -1.0_dp

!  update the elements of sim and simi, and set the next x.

temp = 0.0_dp
do i=1,n
  dx(i) = dxsign*dx(i)
  sim(i,jdrop) = dx(i)
  temp = temp + simi(jdrop,i)*dx(i)
end do
simi(jdrop,1:n) = simi(jdrop,1:n) / temp
do j=1,n
  if (j /= jdrop) then
    temp = dot_product( simi(j,1:n), dx(1:n) )
    simi(j,1:n) = simi(j,1:n) - temp*simi(jdrop,1:n)
  end if
  x(j) = sim(j,np) + dx(j)
end do
go to 40

!  calculate dx = x(*)-x(0).
!  branch if the length of dx is less than 0.5*rho.

370 call trstlp (n, m, a, con, rho, dx, ifull)
if (ifull == 0) then
  temp = sum( dx(1:n)**2 )
  if (temp < 0.25_dp*rho*rho) then
    ibrnch = 1
    go to 550
  end if
end if

!  predict the change to f and the new maximum constraint violation if the
!  variables are altered from x(0) to x(0) + dx.

resnew = 0.0_dp
con(mp) = 0.0_dp
do k=1,mp
  total = con(k) - dot_product( a(1:n,k), dx(1:n) )
  if (k < mp) resnew = max(resnew, total)
end do

!  increase parmu if necessary and branch back if this change alters the
!  optimal vertex. otherwise prerem and prerec will be set to the predicted
!  reductions in the merit function and the maximum constraint violation
!  respectively.

barmu = 0.0_dp
prerec = datmat(mpp,np) - resnew
if (prerec > 0.0_dp) barmu = total/prerec
if (parmu < 1.5_dp*barmu) then
  parmu = 2.0_dp*barmu
  if (iprint >= 2) write(*, 410) parmu
  410 format (/'   increase in parmu to', g13.6)
  phi = datmat(mp,np) + parmu*datmat(mpp,np)
  do j=1,n
    temp = datmat(mp,j) + parmu*datmat(mpp,j)
    if (temp < phi) go to 140
    if (temp == phi .and. parmu == 0.0) then
      if (datmat(mpp,j) < datmat(mpp,np)) go to 140
    end if
  end do
end if
prerem = parmu*prerec - total

!  calculate the constraint and objective functions at x(*).
!  then find the actual reduction in the merit function.

x(1:n) = sim(1:n,np) + dx(1:n)
ibrnch = 1
go to 40

440 vmold = datmat(mp,np) + parmu*datmat(mpp,np)
vmnew = f + parmu*resmax
trured = vmold - vmnew
if (parmu == 0.0_dp .and. f == datmat(mp,np)) then
  prerem = prerec
  trured = datmat(mpp,np) - resmax
end if

!  begin the operations that decide whether x(*) should replace one of the
!  vertices of the current simplex, the change being mandatory if trured is
!  positive. firstly, jdrop is set to the index of the vertex that is to be
!  replaced.

ratio = 0.0_dp
if (trured <= 0.0) ratio = 1.0
jdrop = 0
do j=1,n
  temp = dot_product( simi(j,1:n), dx(1:n) )
  temp = abs(temp)
  if (temp > ratio) then
    jdrop = j
    ratio = temp
  end if
  sigbar(j) = temp*vsig(j)
end do

!  calculate the value of ell.

edgmax = delta*rho
l = 0
do j=1,n
  if (sigbar(j) >= parsig .or. sigbar(j) >= vsig(j)) then
    temp = veta(j)
    if (trured > 0.0_dp) then
      temp = sum( (dx(1:n) - sim(1:n,j))**2 )
      temp = sqrt(temp)
    end if
    if (temp > edgmax) then
      l = j
      edgmax = temp
    end if
  end if
end do
if (l > 0) jdrop = l
if (jdrop == 0) go to 550

!  revise the simplex by updating the elements of sim, simi and datmat.

temp = 0.0_dp
do i=1,n
  sim(i,jdrop) = dx(i)
  temp = temp + simi(jdrop,i)*dx(i)
end do
simi(jdrop,1:n) = simi(jdrop,1:n) / temp
do j=1,n
  if (j /= jdrop) then
    temp = dot_product( simi(j,1:n), dx(1:n) )
    simi(j,1:n) = simi(j,1:n) - temp*simi(jdrop,1:n)
  end if
end do
datmat(1:mpp,jdrop) = con(1:mpp)

!  branch back for further iterations with the current rho.

if (trured > 0.0_dp .and. trured >= 0.1_dp*prerem) go to 140
550 if (iflag == 0) then
  ibrnch = 0
  go to 140
end if

!  otherwise reduce rho if it is not at its least value and reset parmu.

if (rho > rhoend) then
  rho = 0.5_dp*rho
  if (rho <= 1.5_dp*rhoend) rho = rhoend
  if (parmu > 0.0_dp) then
    denom = 0.0_dp
    do k=1,mp
      cmin = datmat(k,np)
      cmax = cmin
      do i=1,n
        cmin = min(cmin, datmat(k,i))
        cmax = max(cmax, datmat(k,i))
      end do
      if (k <= m .and. cmin < 0.5_dp*cmax) then
        temp = max(cmax,0.0_dp) - cmin
        if (denom <= 0.0_dp) then
          denom = temp
        else
          denom = min(denom,temp)
        end if
      end if
    end do
    if (denom == 0.0_dp) then
      parmu = 0.0_dp
    else if (cmax - cmin < parmu*denom) then
      parmu = (cmax - cmin)/denom
    end if
  end if
  if (iprint >= 2) write(*, 580) rho,parmu
  580 format (/'   reduction in rho to ', g13.6, '  and parmu = ', g13.6)
  if (iprint == 2) then
    write(*, 70) nfvals, datmat(mp,np), datmat(mpp,np), sim(1:iptem,np)
    if (iptem < n) write(*, 80) x(iptemp:n)
  end if
  go to 140
end if

!  return the best calculated values of the variables.

if (iprint >= 1) write(*, 590)
590 format (/'   normal return from subroutine cobyla')
if (ifull == 1) go to 620

600 x(1:n) = sim(1:n,np)
f = datmat(mp,np)
resmax = datmat(mpp,np)
620 if (iprint >= 1) then
  write(*, 70) nfvals, f, resmax, x(1:iptem)
  if (iptem < n) write(*, 80) x(iptemp:n)
end if
maxfun = nfvals

return
end subroutine cobylb
!------------------------------------------------------------------------------

subroutine trstlp (n, m, a, b, rho, dx, ifull)

! n.b. arguments z, zdota, vmultc, sdirn, dxnew, vmultd & iact have been removed.

integer, intent(in)     :: n
integer, intent(in)     :: m
real (dp), intent(in)   :: a(:,:)
real (dp), intent(in)   :: b(:)
real (dp), intent(in)   :: rho
real (dp), intent(out)  :: dx(:)
integer, intent(out)    :: ifull

!  this subroutine calculates an n-component vector dx by applying the
!  following two stages. in the first stage, dx is set to the shortest
!  vector that minimizes the greatest violation of the constraints
!    a(1,k)*dx(1)+a(2,k)*dx(2)+...+a(n,k)*dx(n) .ge. b(k), k = 2,3,...,m,
!  subject to the euclidean length of dx being at most rho. if its length is
!  strictly less than rho, then we use the resultant freedom in dx to
!  minimize the objective function
!           -a(1,m+1)*dx(1) - a(2,m+1)*dx(2) - ... - a(n,m+1)*dx(n)
!  subject to no increase in any greatest constraint violation. this
!  notation allows the gradient of the objective function to be regarded as
!  the gradient of a constraint. therefore the two stages are distinguished
!  by mcon .eq. m and mcon .gt. m respectively. it is possible that a
!  degeneracy may prevent dx from attaining the target length rho. then the
!  value ifull = 0 would be set, but usually ifull = 1 on return.

!  in general nact is the number of constraints in the active set and
!  iact(1),...,iact(nact) are their indices, while the remainder of iact
!  contains a permutation of the remaining constraint indices.  further, z
!  is an orthogonal matrix whose first nact columns can be regarded as the
!  result of gram-schmidt applied to the active constraint gradients.  for
!  j = 1,2,...,nact, the number zdota(j) is the scalar product of the j-th
!  column of z with the gradient of the j-th active constraint.  dx is the
!  current vector of variables and here the residuals of the active
!  constraints should be zero. further, the active constraints have
!  nonnegative lagrange multipliers that are held at the beginning of
!  vmultc. the remainder of this vector holds the residuals of the inactive
!  constraints at dx, the ordering of the components of vmultc being in
!  agreement with the permutation of the indices of the constraints that is
!  in iact. all these residuals are nonnegative, which is achieved by the
!  shift resmax that makes the least residual zero.

!  initialize z and some other variables. the value of resmax will be
!  appropriate to dx = 0, while icon will be the index of a most violated
!  constraint if resmax is positive. usually during the first stage the
!  vector sdirn gives a search direction that reduces all the active
!  constraint violations by one simultaneously.

! local variables

real (dp) :: z(n,n), zdota(m+1), vmultc(m+1), sdirn(n), dxnew(n), vmultd(m+1)
real (dp) :: acca, accb, alpha, beta, dd, optnew, optold, ratio, resmax,   &
             resold, sd, sp, spabs, ss, step, stpful, sumabs, temp, tempa, &
             tot, total, vsave, zdotv, zdotw, zdvabs, zdwabs
integer   :: i, iact(m+1), icon, icount, isave, k, kk, kl, kp, kw, mcon,   &
             nact, nactx

ifull = 1
mcon = m
nact = 0
resmax = 0.0_dp
do i=1,n
  z(i,1:n) = 0.0_dp
  z(i,i) = 1.0_dp
  dx(i) = 0.0_dp
end do
if (m >= 1) then
  do k=1,m
    if (b(k) > resmax) then
      resmax = b(k)
      icon = k
    end if
  end do
  do k=1,m
    iact(k) = k
    vmultc(k) = resmax - b(k)
  end do
end if
if (resmax == 0.0_dp) go to 480
sdirn(1:n) = 0.0_dp

!  end the current stage of the calculation if 3 consecutive iterations
!  have either failed to reduce the best calculated value of the objective
!  function or to increase the number of active constraints since the best
!  value was calculated. this strategy prevents cycling, but there is a
!  remote possibility that it will cause premature termination.

60 optold = 0.0_dp
icount = 0
70 if (mcon == m) then
  optnew = resmax
else
  optnew = - dot_product( dx(1:n), a(1:n,mcon) )
end if
if (icount == 0 .or. optnew < optold) then
  optold = optnew
  nactx = nact
  icount = 3
else if (nact > nactx) then
  nactx = nact
  icount = 3
else
  icount = icount - 1
  if (icount == 0) go to 490
end if

!  if icon exceeds nact, then we add the constraint with index iact(icon) to
!  the active set. apply givens rotations so that the last n-nact-1 columns
!  of z are orthogonal to the gradient of the new constraint, a scalar
!  product being set to zero if its nonzero value could be due to computer
!  rounding errors. the array dxnew is used for working space.

if (icon <= nact) go to 260
kk = iact(icon)
dxnew(1:n) = a(1:n,kk)
tot = 0.0_dp
k = n
100 if (k > nact) then
  sp = 0.0_dp
  spabs = 0.0_dp
  do i=1,n
    temp = z(i,k)*dxnew(i)
    sp = sp + temp
    spabs = spabs + abs(temp)
  end do
  acca = spabs + 0.1_dp*abs(sp)
  accb = spabs + 0.2_dp*abs(sp)
  if (spabs >= acca .or. acca >= accb) sp = 0.0_dp
  if (tot == 0.0_dp) then
    tot = sp
  else
    kp = k + 1
    temp = sqrt(sp*sp + tot*tot)
    alpha = sp/temp
    beta = tot/temp
    tot = temp
    do i=1,n
      temp = alpha*z(i,k) + beta*z(i,kp)
      z(i,kp) = alpha*z(i,kp) - beta*z(i,k)
      z(i,k) = temp
    end do
  end if
  k = k - 1
  go to 100
end if

!  add the new constraint if this can be done without a deletion from the
!  active set.

if (tot /= 0.0_dp) then
  nact = nact + 1
  zdota(nact) = tot
  vmultc(icon) = vmultc(nact)
  vmultc(nact) = 0.0_dp
  go to 210
end if

!  the next instruction is reached if a deletion has to be made from the
!  active set in order to make room for the new active constraint, because
!  the new constraint gradient is a linear combination of the gradients of
!  the old active constraints.  set the elements of vmultd to the multipliers
!  of the linear combination.  further, set iout to the index of the
!  constraint to be deleted, but branch if no suitable index can be found.

ratio = -1.0_dp
k = nact
130 zdotv = 0.0_dp
zdvabs = 0.0_dp
if (k < 1) return ! added by Vivek
! print*,"in trstlp, k, size(z,2)=",k,size(z,2)," lbound(z)=",lbound(z) !! debug
do i=1,n
  temp = z(i,k)*dxnew(i)
  zdotv = zdotv + temp
  zdvabs = zdvabs + abs(temp)
end do
acca = zdvabs + 0.1_dp*abs(zdotv)
accb = zdvabs + 0.2_dp*abs(zdotv)
if (zdvabs < acca .and. acca < accb) then
  temp = zdotv/zdota(k)
  if (temp > 0.0_dp .and. iact(k) <= m) then
    tempa = vmultc(k)/temp
    if (ratio < 0.0_dp .or. tempa < ratio) then
      ratio = tempa
    end if
  end if
  if (k >= 2) then
    kw = iact(k)
    dxnew(1:n) = dxnew(1:n) - temp*a(1:n,kw)
  end if
  vmultd(k) = temp
else
  vmultd(k) = 0.0_dp
end if
k = k - 1
if (k > 0) go to 130
if (ratio < 0.0_dp) go to 490

!  revise the lagrange multipliers and reorder the active constraints so
!  that the one to be replaced is at the end of the list. also calculate the
!  new value of zdota(nact) and branch if it is not acceptable.

do k=1,nact
  vmultc(k) = max(0.0_dp,vmultc(k) - ratio*vmultd(k))
end do
if (icon < nact) then
  isave = iact(icon)
  vsave = vmultc(icon)
  k = icon
  170 kp = k + 1
  kw = iact(kp)
  sp = dot_product( z(1:n,k), a(1:n,kw) )
  temp = sqrt(sp*sp + zdota(kp)**2)
  alpha = zdota(kp)/temp
  beta = sp/temp
  zdota(kp) = alpha*zdota(k)
  zdota(k) = temp
  do i=1,n
    temp = alpha*z(i,kp) + beta*z(i,k)
    z(i,kp) = alpha*z(i,k) - beta*z(i,kp)
    z(i,k) = temp
  end do
  iact(k) = kw
  vmultc(k) = vmultc(kp)
  k = kp
  if (k < nact) go to 170
  iact(k) = isave
  vmultc(k) = vsave
end if
temp = dot_product( z(1:n,nact), a(1:n,kk) )
if (temp == 0.0_dp) go to 490
zdota(nact) = temp
vmultc(icon) = 0.0_dp
vmultc(nact) = ratio

!  update iact and ensure that the objective function continues to be
!  treated as the last active constraint when mcon>m.

210 iact(icon) = iact(nact)
iact(nact) = kk
if (mcon > m .and. kk /= mcon) then
  k = nact - 1
  sp = dot_product( z(1:n,k), a(1:n,kk) )
  temp = sqrt(sp*sp + zdota(nact)**2)
  alpha = zdota(nact)/temp
  beta = sp/temp
  zdota(nact) = alpha*zdota(k)
  zdota(k) = temp
  do i=1,n
    temp = alpha*z(i,nact) + beta*z(i,k)
    z(i,nact) = alpha*z(i,k) - beta*z(i,nact)
    z(i,k) = temp
  end do
  iact(nact) = iact(k)
  iact(k) = kk
  temp = vmultc(k)
  vmultc(k) = vmultc(nact)
  vmultc(nact) = temp
end if

!  if stage one is in progress, then set sdirn to the direction of the next
!  change to the current vector of variables.

if (mcon > m) go to 320
kk = iact(nact)
temp = dot_product( sdirn(1:n), a(1:n,kk) )
temp = temp - 1.0_dp
temp = temp/zdota(nact)
sdirn(1:n) = sdirn(1:n) - temp*z(1:n,nact)
go to 340

!  delete the constraint that has the index iact(icon) from the active set.

260 if (icon < nact) then
  isave = iact(icon)
  vsave = vmultc(icon)
  k = icon
  do
    kp = k + 1
    kk = iact(kp)
    sp = dot_product( z(1:n,k), a(1:n,kk) )
    temp = sqrt(sp*sp + zdota(kp)**2)
    alpha = zdota(kp)/temp
    beta = sp/temp
    zdota(kp) = alpha*zdota(k)
    zdota(k) = temp
    do i=1,n
      temp = alpha*z(i,kp) + beta*z(i,k)
      z(i,kp) = alpha*z(i,k) - beta*z(i,kp)
      z(i,k) = temp
    end do
    iact(k) = kk
    vmultc(k) = vmultc(kp)
    k = kp
    if (k >= nact) exit
  end do
  iact(k) = isave
  vmultc(k) = vsave
end if
nact = nact - 1

!  if stage one is in progress, then set sdirn to the direction of the next
!  change to the current vector of variables.

if (mcon > m) go to 320
temp = dot_product( sdirn(1:n), z(1:n,nact+1) )
sdirn(1:n) = sdirn(1:n) - temp*z(1:n,nact+1)
go to 340

!  pick the next search direction of stage two.
! print*,"nact, lbound(zdota), ubound(zdota)=",nact, lbound(zdota), ubound(zdota) !! debug
320 temp = 1.0_dp/zdota(nact)
sdirn(1:n) = temp*z(1:n,nact)

!  calculate the step to the boundary of the trust region or take the step
!  that reduces resmax to zero. the two statements below that include the
!  factor 1.0e-6 prevent some harmless underflows that occurred in a test
!  calculation. further, we skip the step if it could be zero within a
!  reasonable tolerance for computer rounding errors.

340 dd = rho*rho
sd = 0.0_dp
ss = 0.0_dp
do i=1,n
  if (abs(dx(i)) >= 1.0e-6*rho) dd = dd - dx(i)**2
  sd = sd + dx(i)*sdirn(i)
  ss = ss + sdirn(i)**2
end do
if (dd <= 0.0_dp) go to 490
temp = sqrt(ss*dd)
if (abs(sd) >= 1.0e-6*temp) temp = sqrt(ss*dd + sd*sd)
stpful = dd/(temp + sd)
step = stpful
if (mcon == m) then
  acca = step + 0.1_dp*resmax
  accb = step + 0.2_dp*resmax
  if (step >= acca .or. acca >= accb) go to 480
  step = min(step,resmax)
end if

!  set dxnew to the new variables if step is the steplength, and reduce
!  resmax to the corresponding maximum residual if stage one is being done.
!  because dxnew will be changed during the calculation of some lagrange
!  multipliers, it will be restored to the following value later.

dxnew(1:n) = dx(1:n) + step*sdirn(1:n)
if (mcon == m) then
  resold = resmax
  resmax = 0.0_dp
  do k=1,nact
    kk = iact(k)
    temp = b(kk) - dot_product( a(1:n,kk), dxnew(1:n) )
    resmax = max(resmax,temp)
  end do
end if

!  set vmultd to the vmultc vector that would occur if dx became dxnew. a
!  device is included to force vmultd(k) = 0.0 if deviations from this value
!  can be attributed to computer rounding errors. first calculate the new
!  lagrange multipliers.

k = nact
390 zdotw = 0.0_dp
zdwabs = 0.0_dp
do i=1,n
  temp = z(i,k)*dxnew(i)
  zdotw = zdotw + temp
  zdwabs = zdwabs + abs(temp)
end do
acca = zdwabs + 0.1_dp*abs(zdotw)
accb = zdwabs + 0.2_dp*abs(zdotw)
if (zdwabs >= acca .or. acca >= accb) zdotw = 0.0_dp
vmultd(k) = zdotw / zdota(k)
if (k >= 2) then
  kk = iact(k)
  dxnew(1:n) = dxnew(1:n) - vmultd(k)*a(1:n,kk)
  k = k - 1
  go to 390
end if
if (mcon > m) vmultd(nact) = max(0.0_dp,vmultd(nact))

!  complete vmultc by finding the new constraint residuals.

dxnew(1:n) = dx(1:n) + step*sdirn(1:n)
if (mcon > nact) then
  kl = nact + 1
  do k=kl,mcon
    kk = iact(k)
    total = resmax - b(kk)
    sumabs = resmax + abs(b(kk))
    do i=1,n
      temp = a(i,kk)*dxnew(i)
      total = total + temp
      sumabs = sumabs + abs(temp)
    end do
    acca = sumabs + 0.1*abs(total)
    accb = sumabs + 0.2*abs(total)
    if (sumabs >= acca .or. acca >= accb) total = 0.0
    vmultd(k) = total
  end do
end if

!  calculate the fraction of the step from dx to dxnew that will be taken.

ratio = 1.0_dp
icon = 0
do k=1,mcon
  if (vmultd(k) < 0.0_dp) then
    temp = vmultc(k)/(vmultc(k) - vmultd(k))
    if (temp < ratio) then
      ratio = temp
      icon = k
    end if
  end if
end do

!  update dx, vmultc and resmax.

temp = 1.0_dp - ratio
dx(1:n) = temp*dx(1:n) + ratio*dxnew(1:n)
do k=1,mcon
  vmultc(k) = max(0.0_dp,temp*vmultc(k) + ratio*vmultd(k))
end do
if (mcon == m) resmax = resold + ratio*(resmax - resold)

!  if the full step is not acceptable then begin another iteration.
!  otherwise switch to stage two or end the calculation.

if (icon > 0) go to 70
if (step == stpful) go to 500
480 mcon = m + 1
icon = mcon
iact(mcon) = mcon
vmultc(mcon) = 0.0_dp
go to 60

!  we employ any freedom that may be available to reduce the objective
!  function before returning a dx whose length is less than rho.

490 if (mcon == m) go to 480
ifull = 0

500 return
end subroutine trstlp

end module cobyla2_mod