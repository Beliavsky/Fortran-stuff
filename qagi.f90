module adapt_quad_infinite
use constants_mod, only: dpmpar
use kind_mod     , only: dp
implicit none
private 
public :: qagi
contains
!
subroutine qagi(f,bound,inf,epsabs,epsrel,result,abserr,neval,ier,limit,last)
!                   Integration over infinite intervals
!  PURPOSE
!     THE RoutinE CALCULATES AN APPROXIMATION  result  to A GIVEN
!     inTEGRAL    I = inTEGRAL OF  F  OVER (BOUND,+inFinITY)
!              OR I = inTEGRAL OF  F  OVER (-inFinITY,BOUND)
!              OR I = inTEGRAL OF  F  OVER (-inFinITY,+inFinITY),
!     HOPEFULLY SATISFYinG FOLLOWinG CLAIM FOR ACCURACY
!     abs(I-result).LE.max(EPSabs,EPSREL*abs(I)).

!  parameterS
!   ON ENTRY
!      F      - real
!               function SUBPROGRAM DEFininG THE inTEGRAND
!               function F(X). THE ACTUAL NAME FOR F NEEDS to BE
!               DECLARED E X T E R N A L in THE DRIVER PROGRAM.

!      BOUND  - real
!               FinITE BOUND OF inTEGRATION RANGE
!               (HAS NO MEANinG if inTERVAL IS doUBLY-inFinITE)

!      inF    - integer
!               inDICATinG THE kind OF inTEGRATION RANGE inVOLVED
!               inF = 1 CORRESPONDS to  (BOUND,+inFinITY),
!               inF = -1            to  (-inFinITY,BOUND),
!               inF = 2             to (-inFinITY,+inFinITY).

!      EPSabs - real
!               absOLUTE ACCURACY REQUESTED

!      EPSREL - real
!               RELATIVE ACCURACY REQUESTED

!   ON return
!      result - real
!               APPROXIMATION to THE inTEGRAL

!      absERR - real
!               ESTIMATE OF THE MODULUS OF THE absOLUTE ERROR,
!               WHICH SHOULD EQUAL OR EXCEED abs(I-result)

!      NEVAL  - integer
!               NUMBER OF inTEGRAND EVALUATIONS

!      IER    - integer
!               IER = 0 NORMAL AND RELIABLE TERMinATION OF THE RoutinE.
!                       IT IS ASSUMED THAT THE REQUESTED
!                       ACCURACY HAS BEEN ACHIEVED.
!             - IER.GT.0 ABNORMAL TERMinATION OF THE RoutinE.  THE ESTIMATES
!                       FOR result AND ERROR ARE LESS RELIABLE. IT IS ASSUMED
!                       THAT THE REQUESTED ACCURACY HAS NOT BEEN ACHIEVED.

!               IER = 1 maxIMUM NUMBER OF SUBDIVISIONS ALLOWED HAS BEEN
!                       ACHIEVED.  ONE CAN ALLOW MORE SUBDIVISIONS BY
!                       inCREASinG THE VALUE OF LIMIT (AND TAKinG THE
!                       ACCORDinG dimension ADJUSTMENTS into ACCOUNT).
!                       HOWEVER, if THIS YIELDS NO IMPROVEMENT IT IS ADVISED
!                       to ANALYZE THE inTEGRAND in ORDER to DETERMinE THE
!                       inTEGRATION DifFICULTIES.
!                       if THE POSITION OF A LOCAL DifFICULTY CAN BE DETERMinED
!                       (E.G. SinGULARITY, DISCONTinUITY WITHin THE inTERVAL)
!                       ONE WILL PROBABLY GAin FROM SPLITTinG UP THE inTERVAL
!                       AT THIS POinT AND CALLinG THE inTEGRAtoR ON THE
!                       SUBRANGES.
!                       if POSSIBLE, AN APPROPRIATE SPECIAL-PURPOSE inTEGRAtoR
!                       SHOULD BE useD, WHICH IS DESIGNED FOR HANDLinG THE
!                       TYPE OF DifFICULTY inVOLVED.
!                   = 2 THE OCCURRENCE OF ROUNdoFF ERROR IS DETECTED, WHICH
!                       PREVENTS THE REQUESTED toLERANCE FROM BEinG ACHIEVED.
!                       THE ERROR MAY BE UNDER-ESTIMATED.
!                   = 3 EXTREMELY BAD inTEGRAND BEHAVIOUR OCCURS AT SOME POinTS
!                       OF THE inTEGRATION inTERVAL.
!                   = 4 THE ALgoRITHM doES NOT CONVERGE.
!                       ROUNdoFF ERROR IS DETECTED in THE EXTRAPOLATION TABLE.
!                       IT IS ASSUMED THAT THE REQUESTED toLERANCE CANNOT BE
!                       ACHIEVED, AND THAT THE returnED result IS THE BEST
!                       WHICH CAN BE OBTAinED.
!                   = 5 THE inTEGRAL IS PROBABLY DIVERGENT, OR SLOWLY
!                       CONVERGENT. IT MUST BE NOTED THAT DIVERGENCE CAN OCCUR
!                       WITH ANY OTHER VALUE OF IER.
!                   = 6 THE inPUT IS inVALID BECAuse EPSabs OR EPSREL IS
!                       NEGATIVE, LIMIT .LT. 1, OR LENW .LT. 4 * LIMIT.
!                       result, absERR, NEVAL, LAST ARE SET to ZERO.

!   dimensioninG parameterS
!      LIMIT - integer
!              dimensioninG parameter FOR IWORK
!              LIMIT DETERMinES THE maxIMUM NUMBER OF SUBinTERVALS in THE
!              PARTITION OF THE GIVEN inTEGRATION inTERVAL (A,B), LIMIT.GE.1.
!              if LIMIT.LT.1, THE RoutinE WILL end WITH IER = 6.

!      LENW  - integer
!              dimensioninG parameter FOR WORK
!              LENW MUST BE AT LEAST LIMIT*4.
!              if LENW.LT.LIMIT*4, THE RoutinE WILL end WITH IER = 6.

!      LAST  - integer
!              ON return, LAST EQUALS THE NUMBER OF SUBinTERVALS PRODUCED in
!              THE SUBDIVISION PROCESS, WHICH DETERMinES THE NUMBER OF
!              SIGNifICANT ELEMENTS ACTUALLY in THE WORK ARRAYS.

!  subroutineS OR functionS NEEDED
!        - QAGIE
!        - QK15I
!        - QPSRT
!        - QELG
!        - F (useR PROVIDED function)
!        - SPMPAR

!-----------------------------------------------------------------------

real (dp), intent(in)   :: bound  ! finite bound of integration range -- has no meaning if integral is double infinite
real (dp), intent(in)   :: epsabs ! absolute accuracy requested
real (dp), intent(in)   :: epsrel ! relative accuracy requested
integer  , intent(in)   :: inf    ! indicates integration range: (1,(bound,Inf)) (-1,(-Inf,bound)) (2,-Inf,Inf)
integer  , intent(in)   :: limit  ! determines maximum number of subintervals in the partition of the integration region
real (dp), intent(out)  :: result ! approximation to the integral
real (dp), intent(out)  :: abserr ! estimate of the modulus of the absolute error
integer  , intent(out)  :: neval  ! # of integrand evaluations
integer  , intent(out)  :: ier    ! error flag -- returned zero if no error
integer  , intent(out)  :: last   ! # of subintervals produced in the subdivision process

interface
  function f(x) result(fx)
    use constants_mod
    use kind_mod, only: dp
    implicit none
    real (dp), intent(in)  :: x
    real (dp)              :: fx
  end function f
end interface

!         CHECK VALIDITY OF LIMIT.

ier = 6
neval = 0
last = 0
result = 0.0D0
abserr = 0.0D0
if (limit < 1) return

!         PREPARE CALL FOR QAGIE.

CALL qagie (f, bound, inf, epsabs, epsrel, limit, result, abserr,  &
            neval, ier, last)
return
end subroutine qagi



subroutine qagie (f, bound, inf, epsabs, epsrel, limit, result, abserr,  &
                  neval, ier, last)
!-----------------------------------------------------------------------

!                   inTEGRATION OVER inFinITE inTERVALS

!-----------------------------------------------------------------------

!  PURPOSE
!     THE RoutinE CALCULATES AN APPROXIMATION result to A GIVEN inTEGRAL
!                 I = inTEGRAL OF  F  OVER (BOUND,+inFinITY)
!              OR I = inTEGRAL OF  F  OVER (-inFinITY,BOUND)
!              OR I = inTEGRAL OF  F  OVER (-inFinITY,+inFinITY),
!     HOPEFULLY SATISFYinG FOLLOWinG CLAIM FOR ACCURACY
!     abs(I-result).LE.max(EPSabs,EPSREL*abs(I)).

!  parameterS
!   ON ENTRY
!      F      - real
!               function SUBPROGRAM DEFininG THE inTEGRAND
!               function F(X). THE ACTUAL NAME FOR F NEEDS to BE
!               DECLARED E X T E R N A L in THE DRIVER PROGRAM.

!      BOUND  - real
!               FinITE BOUND OF inTEGRATION RANGE
!               (HAS NO MEANinG if inTERVAL IS doUBLY-inFinITE)

!      inF    - integer
!               inDICATinG THE kind OF inTEGRATION RANGE inVOLVED
!               inF = 1 CORRESPONDS to  (BOUND,+inFinITY),
!               inF = -1            to  (-inFinITY,BOUND),
!               inF = 2             to (-inFinITY,+inFinITY).

!      EPSabs - real
!               absOLUTE ACCURACY REQUESTED

!      EPSREL - real
!               RELATIVE ACCURACY REQUESTED

!      LIMIT  - integer
!               GIVES AN UPPER BOUND ON THE NUMBER OF SUBinTERVALS in THE
!               PARTITION OF (A,B),
!               LIMIT.GE.1

!   ON return
!      result - real
!               APPROXIMATION to THE inTEGRAL

!      absERR - real
!               ESTIMATE OF THE MODULUS OF THE absOLUTE ERROR,
!               WHICH SHOULD EQUAL OR EXCEED abs(I-result)

!      NEVAL  - integer
!               NUMBER OF inTEGRAND EVALUATIONS

!      IER    - integer
!               IER = 0 NORMAL AND RELIABLE TERMinATION OF THE RoutinE.
!                       IT IS ASSUMED THAT THE REQUESTED ACCURACY HAS BEEN
!                       ACHIEVED.
!             - IER.GT.0 ABNORMAL TERMinATION OF THE RoutinE.  THE ESTIMATES
!                       FOR result AND ERROR ARE LESS RELIABLE.  IT IS ASSUMED
!                       THAT THE REQUESTED ACCURACY HAS NOT BEEN ACHIEVED.

!               IER = 1 maxIMUM NUMBER OF SUBDIVISIONS ALLOWED HAS BEEN
!                       ACHIEVED.  ONE CAN ALLOW MORE SUBDIVISIONS BY
!                       inCREASinG THE VALUE OF LIMIT (AND TAKinG THE
!                       ACCORDinG dimension ADJUSTMENTS into ACCOUNT).
!                       HOWEVER, if THIS YIELDS NO IMPROVEMENT IT IS ADVISED
!                       to ANALYZE THE inTEGRAND in ORDER to DETERMinE THE
!                       inTEGRATION DifFICULTIES.
!                       if THE POSITION OF A LOCAL DifFICULTY CAN BE
!                       DETERMinED (E.G. SinGULARITY, DISCONTinUITY WITHin THE
!                       inTERVAL) ONE WILL PROBABLY GAin FROM SPLITTinG UP THE
!                       inTERVAL AT THIS POinT AND CALLinG THE inTEGRAtoR ON
!                       THE SUBRANGES.  if POSSIBLE, AN APPROPRIATE SPECIAL-
!                       PURPOSE inTEGRAtoR SHOULD BE useD, WHICH IS DESIGNED
!                       FOR HANDLinG THE TYPE OF DifFICULTY inVOLVED.
!                   = 2 THE OCCURRENCE OF ROUNdoFF ERROR IS DETECTED, WHICH
!                       PREVENTS THE REQUESTED toLERANCE FROM BEinG ACHIEVED.
!                       THE ERROR MAY BE UNDER-ESTIMATED.
!                   = 3 EXTREMELY BAD inTEGRAND BEHAVIOUR OCCURS AT SOME
!                       POinTS OF THE inTEGRATION inTERVAL.
!                   = 4 THE ALgoRITHM doES NOT CONVERGE.
!                       ROUNdoFF ERROR IS DETECTED in THE EXTRAPOLATION TABLE.
!                       IT IS ASSUMED THAT THE REQUESTED toLERANCE CANNOT BE
!                       ACHIEVED, AND THAT THE returnED result IS THE BEST
!                       WHICH CAN BE OBTAinED.
!                   = 5 THE inTEGRAL IS PROBABLY DIVERGENT, OR SLOWLY
!                       CONVERGENT.  IT MUST BE NOTED THAT DIVERGENCE CAN
!                       OCCUR WITH ANY OTHER VALUE OF IER.
!                   = 6 THE inPUT IS inVALID BECAuse EPSabs OR EPSREL IS
!                       NEGATIVE.
!                       result, absERR, NEVAL, LAST, RLIST(1), ELIST(1) AND
!                       IORD(1) ARE SET to ZERO.
!                       ALIST(1) AND BLIST(1) ARE SET to 0 AND 1 RESPECTIVELY.

!      ALIST  - real
!               VECtoR OF dimension AT LEAST LIMIT, THE FIRST
!                LAST  ELEMENTS OF WHICH ARE THE LEFT
!               end POinTS OF THE SUBinTERVALS in THE PARTITION
!               OF THE TRANSFORMED inTEGRATION RANGE (0,1).

!      BLIST  - real
!               VECtoR OF dimension AT LEAST LIMIT, THE FIRST
!                LAST  ELEMENTS OF WHICH ARE THE RIGHT
!               end POinTS OF THE SUBinTERVALS in THE PARTITION
!               OF THE TRANSFORMED inTEGRATION RANGE (0,1).

!      RLIST  - real
!               VECtoR OF dimension AT LEAST LIMIT, THE FIRST
!                LAST  ELEMENTS OF WHICH ARE THE inTEGRAL
!               APPROXIMATIONS ON THE SUBinTERVALS

!      ELIST  - real
!               VECtoR OF dimension AT LEAST LIMIT,  THE FIRST
!                LAST  ELEMENTS OF WHICH ARE THE MODULI
!               OF THE absOLUTE ERROR ESTIMATES ON THE
!               SUBinTERVALS

!      IORD   - integer
!               VECtoR OF dimension LIMIT, THE FIRST K ELEMENTS OF WHICH ARE
!               POinTERS to THE ERROR ESTIMATES OVER THE SUBinTERVALS,
!               SUCH THAT ELIST(IORD(1)), ..., ELIST(IORD(K)) FORM A
!               DECREASinG SEQUENCE, WITH K = LAST if LAST.LE.(LIMIT/2+2),
!               AND K = LIMIT+1-LAST OTHERWISE

!      LAST   - integer
!               NUMBER OF SUBinTERVALS ACTUALLY PRODUCED in THE SUBDIVISION
!               PROCESS

!  subroutineS OR functionS NEEDED
!        - QK15I
!        - QPSRT
!        - QELG
!        - F (useR-PROVIDED function)
!        - SPMPAR

!-----------------------------------------------------------------------

real (dp), intent(in)   :: bound, epsabs, epsrel
real (dp), intent(out)  :: result, abserr
integer, intent(in)     :: inf, limit
integer, intent(out)    :: neval, ier, last

real (dp) :: abseps, area, area1, area12, area2, a1, a2, boun, b1, b2,     &
             correc, defabs, defab1, defab2, dres, epmach, erlarg, erlast, &
             errbnd, errmax, error1, error2, erro12, errsum, ertest,       &
             oflow, rerr, resabs, reseps, res3la(3), rlist2(52), small, t, &
             uflow
integer :: id, ierro, iroff1, iroff2, iroff3, jupbnd, k, ksgn, ktmin,  &
           maxerr, nres, nrmax, numrl2
logical :: extrap, noext

real (dp), dimension(limit)  :: alist, blist, elist, rlist
integer, dimension(limit)    :: iord

inTERFACE
  function f(x) result(fx)
    use constants_mod
    use kind_mod, only: dp
    implicit none
    real (dp), intent(in) :: x
    real (dp)             :: fx
  end function f
end inTERFACE

!      THE dimension OF RLIST2 IS DETERMinED BY THE VALUE OF
!      LIMEXP in subroutine QELG.


!      LIST OF MAJOR VARIABLES
!      -----------------------

!     ALIST     - LIST OF LEFT end POinTS OF ALL SUBinTERVALS CONSIDERED UP
!                 to NOW
!     BLIST     - LIST OF RIGHT end POinTS OF ALL SUBinTERVALS CONSIDERED UP
!                 to NOW
!     RLIST(I)  - APPROXIMATION to THE inTEGRAL OVER (ALIST(I),BLIST(I))
!     RLIST2    - ARRAY OF dimension AT LEAST (LIMEXP+2), CONTAininG THE PART
!                 OF THE epsilon TABLE WHICH IS STILL NEEDED FOR FURTHER
!                 COMPUTATIONS
!     ELIST(I)  - ERROR ESTIMATE APPLYinG to RLIST(I)
!     maxERR    - POinTER to THE inTERVAL WITH LARGEST ERROR ESTIMATE
!     ERRmax    - ELIST(maxERR)
!     ERLAST    - ERROR ON THE inTERVAL CURRENTLY SUBDIVIDED
!                 (BEFORE THAT SUBDIVISION HAS TAKEN PLACE)
!     AREA      - SUM OF THE inTEGRALS OVER THE SUBinTERVALS
!     ERRSUM    - SUM OF THE ERRORS OVER THE SUBinTERVALS
!     ERRBND    - REQUESTED ACCURACY max(EPSabs, EPSREL*abs(result))
!     *****1    - VARIABLE FOR THE LEFT SUBinTERVAL
!     *****2    - VARIABLE FOR THE RIGHT SUBinTERVAL
!     LAST      - inDEX FOR SUBDIVISION
!     NRES      - NUMBER OF CALLS to THE EXTRAPOLATION RoutinE
!     NUMRL2    - NUMBER OF ELEMENTS CURRENTLY in RLIST2.  if AN APPROPRIATE
!                 APPROXIMATION to THE COMPOUNDED inTEGRAL HAS BEEN OBTAinED,
!                 IT IS PUT in RLIST2(NUMRL2) AFTER NUMRL2 HAS BEEN inCREASED
!                 BY ONE.
!     SMALL     - LENGTH OF THE SMALLEST inTERVAL CONSIDERED UP to NOW,
!                 MULTIPLIED BY 1.5
!     ERLARG    - SUM OF THE ERRORS OVER THE inTERVALS LARGER THAN THE
!                 SMALLEST inTERVAL CONSIDERED UP to NOW
!     EXTRAP    - logical VARIABLE DENOTinG THAT THE RoutinE IS ATTEMPTinG to
!                 PERFORM EXTRAPOLATION. I.E. BEFORE SUBDIVIDinG THE SMALLEST
!                 inTERVAL WE TRY to DECREASE THE VALUE OF ERLARG.
!     NOEXT     - logical VARIABLE DENOTinG THAT EXTRAPOLATION IS NO LONGER
!                 ALLOWED (TRUE-VALUE)

!      MACHinE DEPendENT CONSTANTS
!      ---------------------------

!     EPMACH IS THE LARGEST RELATIVE SPACinG.
!     UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
!     OFLOW IS THE LARGEST POSITIVE MAGNITUDE.

epmach = dpmpar(1)
uflow = dpmpar(2)
oflow = dpmpar(3)

!           CHECK EPSabs AND EPSREL
!           -----------------------

neval = 0
last = 0
result = 0.0D0
abserr = 0.0D0
alist(1) = 0.0D0
blist(1) = 1.0D0
rlist(1) = 0.0D0
elist(1) = 0.0D0
iord(1) = 0
ier = 6
if (epsabs < 0.0D0 .OR. epsrel < 0.0D0) go to 999
ier = 0
rerr = max(epsrel, 50.0D0*epmach, 0.5D-14)

!     FIRST APPROXIMATION to THE inTEGRAL
!     -----------------------------------

!     DETERMinE THE inTERVAL to BE MAPPED ONto (0,1).
!     if inF = 2 THE inTEGRAL IS COMPUTED AS I = I1+I2, WHERE
!     I1 = inTEGRAL OF F OVER (-inFinITY,0),
!     I2 = inTEGRAL OF F OVER (0,+inFinITY).

boun = bound
if (inf == 2) boun = 0.0D0
CALL qk15i (f, boun, inf, 0.0D0, 1.0D0, result, abserr, defabs, resabs, &
            epmach, uflow)

!           TEST ON ACCURACY

last = 1
rlist(1) = result
elist(1) = abserr
iord(1) = 1
dres = abs(result)
errbnd = max(epsabs,rerr*dres)
if (abserr <= 100.0D0*epmach*defabs .AND. abserr > errbnd)  &
ier = 2
if (limit == 1) ier = 1
if (ier /= 0 .OR. (abserr <= errbnd .AND. abserr /= resabs)  &
.OR. abserr == 0.0D0) go to 130

!           inITIALIZATION
!           --------------

rlist2(1) = result
errmax = abserr
maxerr = 1
area = result
errsum = abserr
abserr = oflow
correc = 0.0D0
nrmax = 1
nres = 0
ktmin = 0
numrl2 = 2
extrap = .false.
noext = .false.
ierro = 0
iroff1 = 0
iroff2 = 0
iroff3 = 0
ksgn = -1
if (dres >= (1.0D0 - 50.0D0*epmach)*defabs) ksgn = 1
t = 1.0D0 + 100.0D0*epmach

!           MAin do-LOOP
!           ------------

do last = 2,limit
  
!           BISECT THE SUBinTERVAL WITH NRmax-TH LARGEST ERROR ESTIMATE.
  
  a1 = alist(maxerr)
  b1 = 0.5D0*(alist(maxerr) + blist(maxerr))
  a2 = b1
  b2 = blist(maxerr)
  erlast = errmax
  CALL qk15i (f, boun, inf, a1, b1, area1, error1,  &
              resabs, defab1, epmach, uflow)
  CALL qk15i (f, boun, inf, a2, b2, area2, error2,  &
              resabs, defab2, epmach, uflow)
  
!           IMPROVE PREVIOUS APPROXIMATIONS to inTEGRAL AND ERROR
!           AND TEST FOR ACCURACY.
  
  area12 = area1 + area2
  erro12 = error1 + error2
  errsum = errsum + erro12 - errmax
  area = area + area12 - rlist(maxerr)
  if (defab1 == error1 .OR. defab2 == error2) go to 15
  if (abs(rlist(maxerr) - area12) > 0.1D-04*abs(area12)  &
      .OR. erro12 < 0.99D0*errmax) go to 10
  if (extrap) iroff2 = iroff2 + 1
  if (.NOT.extrap) iroff1 = iroff1 + 1
  10 if (last > 10 .AND. erro12 > errmax) iroff3 = iroff3 + 1
  15 rlist(maxerr) = area1
  rlist(last) = area2
  errbnd = max(epsabs,rerr*abs(area))
  
!           TEST FOR ROUNdoFF ERROR AND EVENTUALLY SET ERROR FLAG.
  
  if (iroff1 + iroff2 >= 10 .OR. iroff3 >= 20) ier = 2
  if (iroff2 >= 5) ierro = 3
  
!           SET ERROR FLAG in THE case THAT THE NUMBER OF
!           SUBinTERVALS EQUALS LIMIT.
  
  if (last == limit) ier = 1
  
!           SET ERROR FLAG in THE case OF BAD inTEGRAND BEHAVIOUR
!           AT SOME POinTS OF THE inTEGRATION RANGE.
  
  
  if (max(abs(a1),abs(b2)) <= t*(abs(a2) + 0.1D+04*uflow)) ier = 4
  
!           APPend THE NEWLY-CREATED inTERVALS to THE LIST.
  
  if (error2 > error1) go to 20
  alist(last) = a2
  blist(maxerr) = b1
  blist(last) = b2
  elist(maxerr) = error1
  elist(last) = error2
  go to 30
  20 alist(maxerr) = a2
  alist(last) = a1
  blist(last) = b1
  rlist(maxerr) = area2
  rlist(last) = area1
  elist(maxerr) = error2
  elist(last) = error1
  
!           CALL subroutine QPSRT to MAinTAin THE DESCendinG ORDERinG
!           in THE LIST OF ERROR ESTIMATES AND select THE SUBinTERVAL WITH
!           NRmax-TH LARGEST ERROR ESTIMATE (to BE BISECTED NEXT).
  
  30 CALL qpsrt (limit,last,maxerr,errmax,elist,iord,nrmax)
  if (errsum <= errbnd) go to 115
  if (ier /= 0) go to 100
  if (last == 2) go to 80
  if (noext) CYCLE
  erlarg = erlarg - erlast
  if (abs(b1 - a1) > small) erlarg = erlarg + erro12
  if (extrap) go to 40
  
!           TEST WHETHER THE inTERVAL to BE BISECTED NEXT IS THE
!           SMALLEST inTERVAL.
  
  if (abs(blist(maxerr) - alist(maxerr)) > small) CYCLE
  extrap = .true.
  nrmax = 2
  40 if (ierro == 3 .OR. erlarg <= ertest) go to 60
  
!           THE SMALLEST inTERVAL HAS THE LARGEST ERROR.
!           BEFORE BISECTinG DECREASE THE SUM OF THE ERRORS
!           OVER THE LARGER inTERVALS (ERLARG) AND PERFORM EXTRAPOLATION.
  
  id = nrmax
  jupbnd = last
  if (last > (2 + limit/2)) jupbnd = limit + 3 - last
  do k = id,jupbnd
    maxerr = iord(nrmax)
    errmax = elist(maxerr)
    if (abs(blist(maxerr) - alist(maxerr)) > small) CYCLE
    nrmax = nrmax + 1
  end do
  
!           PERFORM EXTRAPOLATION.
  
  60 numrl2 = numrl2 + 1
  rlist2(numrl2) = area
  CALL qelg (numrl2, rlist2, reseps, abseps, res3la, nres,epmach, oflow)
  ktmin = ktmin + 1
  if (ktmin > 5 .AND. abserr < 0.1D-02*errsum) ier = 5
  if (abseps >= abserr) go to 70
  ktmin = 0
  abserr = abseps
  result = reseps
  correc = erlarg
  ertest = max(epsabs,rerr*abs(reseps))
  if (abserr <= ertest) go to 100
  
!            PREPARE BISECTION OF THE SMALLEST inTERVAL.
  
  70 if (numrl2 == 1) noext = .true.
  if (ier == 5) go to 100
  maxerr = iord(1)
  errmax = elist(maxerr)
  nrmax = 1
  extrap = .false.
  small = small*0.5D0
  erlarg = errsum
  CYCLE
  80 small = 0.375D0
  erlarg = errsum
  ertest = errbnd
  rlist2(2) = area
end do

!           SET FinAL result AND ERROR ESTIMATE.
!           ------------------------------------

100 if (abserr == oflow) go to 115
if (ier + ierro == 0) go to 110
if (ierro == 3) abserr = abserr + correc
if (ier == 0) ier = 3
if (result /= 0.0D0 .AND. area /= 0.0D0) go to 105
if (abserr > errsum) go to 115
if (area == 0.0D0) go to 130
go to 110
105 if (abserr/abs(result) > errsum/abs(area)) go to 115

!           TEST ON DIVERGENCE

110 if (ksgn == -1 .AND. max(abs(result),abs(area)) <=  &
defabs*0.1D-01) go to 130
if (0.1D-01 > (result/area) .OR. (result/area) > 0.1D+03  &
.OR. errsum > abs(area)) ier = 6
go to 130

!           COMPUTE GLOBAL inTEGRAL SUM.

115 result = 0.0D0
do k = 1,last
  result = result + rlist(k)
end do
abserr = errsum
130 neval = 30*last - 15
if (inf == 2) neval = 2*neval
if (ier > 2) ier = ier - 1
999 return
end subroutine qagie



subroutine qk15i (f, boun, inf, a, b, result, abserr, resabs,  &
                  resasc, epmach, uflow)
!-----------------------------------------------------------------------

! 1.  PURPOSE
!        THE ORIGinAL (inFinITE) inTEGRATION RANGE IS MAPPED
!        ONto THE inTERVAL (0,1) AND (A,B) IS A PART OF (0,1).
!        IT IS THE PURPOSE to COMPUTE
!        I = inTEGRAL OF TRANSFORMED inTEGRAND OVER (A,B),
!        J = inTEGRAL OF abs(TRANSFORMED inTEGRAND) OVER (A,B).

! 2.  parameterS
!      ON ENTRY
!        F      - real
!                 function SUBPROGRAM DEFininG THE inTEGRAND function F(X).
!                 THE ACTUAL NAME FOR F NEEDS to BE DECLARED E X T E R N A L
!                 in THE CALLinG PROGRAM.

!        BOUN   - real
!                 FinITE BOUND OF ORIGinAL inTEGRATION
!                 RANGE (SET to ZERO if inF = +2)

!        inF    - integer
!                 if inF = -1, THE ORIGinAL inTERVAL IS
!                             (-inFinITY,BOUND),
!                 if inF = +1, THE ORIGinAL inTERVAL IS
!                             (BOUND,+inFinITY),
!                 if inF = +2, THE ORIGinAL inTERVAL IS
!                             (-inFinITY,+inFinITY) AND
!                 THE inTEGRAL IS COMPUTED AS THE SUM OF TWO inTEGRALS,
!                 ONE OVER (-inFinITY,0) AND ONE OVER (0,+inFinITY).

!        A      - real
!                 LOWER LIMIT FOR inTEGRATION OVER SUBRANGE OF (0,1)

!        B      - real
!                 UPPER LIMIT FOR inTEGRATION OVER SUBRANGE OF (0,1)

!        EPMACH - real
!                 THE RELATIVE PRECISION OF THE FLOATinG ARITHMETIC BEinG useD.

!        UFLOW  - real
!                 THE SMALLEST POSITIVE MAGNITUDE.

!      ON return
!        result - real
!                 APPROXIMATION to THE inTEGRAL I
!                 result IS COMPUTED BY APPLYinG THE 15-POinT KRONROD RULE
!                 (RESK) OBTAinED BY OPTIMAL ADDITION OF absCISSAE to THE
!                 7-POinT GAUSS RULE(RESG).

!        absERR - real
!                 ESTIMATE OF THE MODULUS OF THE absOLUTE ERROR,
!                 WHICH SHOULD EQUAL OR EXCEED abs(I-result)

!        RESabs - real
!                 APPROXIMATION to THE inTEGRAL J

!        RESASC - real
!                 APPROXIMATION to THE inTEGRAL OF
!                 abs((TRANSFORMED inTEGRAND)-I/(B-A)) OVER (A,B)

! 3.  subroutineS OR functionS NEEDED
!           - F (useR-PROVIDED function)

!-----------------------------------------------------------------------

real (dp), intent(in)   :: boun, a, b, epmach, uflow
integer, intent(in)     :: inf
real (dp), intent(out)  :: result, abserr, resabs, resasc

real (dp) :: fv1(7), fv2(7)

inTERFACE
  function f(x) result(fx)
    use constants_mod
    use kind_mod, only: dp
    implicit none
    real (dp), intent(in) :: x
    real (dp)             :: fx
  end function f
end inTERFACE

!     THE absCISSAE AND WEIGHTS ARE SUPPLIED FOR THE inTERVAL
!     (-1,1).  BECAuse OF SYMMETRY ONLY THE POSITIVE absCISSAE AND
!     THEIR CORRESPONDinG WEIGHTS ARE GIVEN.

!     XGK    - absCISSAE OF THE 15-POinT KRONROD RULE
!              XGK(2), XGK(4), ... absCISSAE OF THE 7-POinT GAUSS RULE
!              XGK(1), XGK(3), ...  absCISSAE WHICH ARE OPTIMALLY
!              ADDED to THE 7-POinT GAUSS RULE

!     WGK    - WEIGHTS OF THE 15-POinT KRONROD RULE

!     WG     - WEIGHTS OF THE 7-POinT GAUSS RULE, CORRESPONDinG to THE
!              absCISSAE XGK(2), XGK(4), ... WG(1), WG(3), ...
!              ARE SET to ZERO.

real (dp) :: absc, absc1, absc2,  centr, dinf, fc, fsum, fval1, fval2,  &
             hlgth, resg, resk, reskh, tabsc1, tabsc2, tol
integer   :: j
real (dp), dimension(8) :: xgk = (/  &
                           0.9914553711208126D+00, 0.9491079123427585D+00,  &
                           0.8648644233597691D+00, 0.7415311855993944D+00,  &
                           0.5860872354676911D+00, 0.4058451513773972D+00,  &
                           0.2077849550078985D+00, 0.0000000000000000D+00 /), &
                           wgk = (/  &
                           0.2293532201052922D-01, 0.6309209262997855D-01,  &
                           0.1047900103222502D+00, 0.1406532597155259D+00,  &
                           0.1690047266392679D+00, 0.1903505780647854D+00,  &
                           0.2044329400752989D+00, 0.2094821410847278D+00 /), &
                           wg = (/   &
                           0.0000000000000000D+00, 0.1294849661688697D+00,  &
                           0.0000000000000000D+00, 0.2797053914892767D+00,  &
                           0.0000000000000000D+00, 0.3818300505051189D+00,  &
                           0.0000000000000000D+00, 0.4179591836734694D+00 /)

!     LIST OF MAJOR VARIABLES
!     -----------------------

!     CENTR  - MID POinT OF THE inTERVAL
!     HLGTH  - HALF-LENGTH OF THE inTERVAL
!     absC*  - absCISSA
!     TabsC* - TRANSFORMED absCISSA
!     FVAL*  - function VALUE
!     RESG   - result OF THE 7-POinT GAUSS FORMULA
!     RESK   - result OF THE 15-POinT KRONROD FORMULA
!     RESKH  - APPROXIMATION to THE MEAN VALUE OF THE TRANSFORMED
!              inTEGRAND OVER (A,B), I.E. to I/(B-A)

dinf = Min(1,inf)

centr = 0.5D0*(a + b)
hlgth = 0.5D0*(b - a)
tabsc1 = boun + dinf*(1.0D0 - centr)/centr
fval1 = f(tabsc1)
if (inf == 2) fval1 = fval1 + f(-tabsc1)
fc = (fval1/centr)/centr

!           COMPUTE THE 15-POinT KRONROD APPROXIMATION to THE inTEGRAL,
!           AND ESTIMATE THE ERROR.

resg = wg(8)*fc
resk = wgk(8)*fc
resabs = abs(resk)
do j = 1,7
  absc = hlgth*xgk(j)
  absc1 = centr - absc
  absc2 = centr + absc
  tabsc1 = boun + dinf*(1.0D0 - absc1)/absc1
  tabsc2 = boun + dinf*(1.0D0 - absc2)/absc2
  fval1 = f(tabsc1)
  fval2 = f(tabsc2)
  if (inf == 2) fval1 = fval1 + f(-tabsc1)
  if (inf == 2) fval2 = fval2 + f(-tabsc2)
  fval1 = (fval1/absc1)/absc1
  fval2 = (fval2/absc2)/absc2
  fv1(j) = fval1
  fv2(j) = fval2
  fsum = fval1 + fval2
  resg = resg + wg(j)*fsum
  resk = resk + wgk(j)*fsum
  resabs = resabs + wgk(j)*(abs(fval1) + abs(fval2))
end do
reskh = resk / 2
resasc = wgk(8)*abs(fc - reskh)
do j = 1,7
  resasc = resasc + wgk(j)*(abs(fv1(j)-reskh) + abs(fv2(j)-reskh))
end do
result = resk*hlgth
resasc = resasc*hlgth
resabs = resabs*hlgth
abserr = abs((resk - resg)*hlgth)
if (resasc /= 0.0D0 .AND. abserr /= 0.0D0) abserr = resasc*  &
                           Min(1.0D0, (0.2D+03*abserr/resasc)**1.5D0)
tol = 50.0D0*epmach
if (resabs > uflow/tol) abserr = max(abserr, tol*resabs)

return
end subroutine qk15i



subroutine qpsrt(limit, last, maxerr, ermax, elist, iord, nrmax)
!     ..................................................................

! 1.  QPSRT
!     ORDERinG RoutinE
!        STANDARD FORTRAN subroutine
!        real VERSION

! 2.  PURPOSE
!        THIS RoutinE MAinTAinS THE DESCendinG ORDERinG in THE LIST OF THE
!        LOCAL ERROR ESTIMATES resultinG FROM THE inTERVAL SUBDIVISION
!        PROCESS.  AT EACH CALL TWO ERROR ESTIMATES ARE inSERTED USinG THE
!        SEQUENTIAL SEARCH METHOD, toP-doWN FOR THE LARGEST ERROR ESTIMATE
!        AND BOTtoM-UP FOR THE SMALLEST ERROR ESTIMATE.

! 3.  CALLinG SEQUENCE
!        CALL QPSRT(LIMIT, LAST, maxERR, ERmax, ELIST, IORD, NRmax)

!     parameterS (MEANinG AT outPUT)
!        LIMIT  - integer
!                 maxIMUM NUMBER OF ERROR ESTIMATES THE LIST CAN CONTAin

!        LAST   - integer
!                 NUMBER OF ERROR ESTIMATES CURRENTLY in THE LIST

!        maxERR - integer
!                 maxERR POinTS to THE NRmax-TH LARGEST ERROR ESTIMATE
!                 CURRENTLY in THE LIST

!        ERmax  - real
!                 NRmax-TH LARGEST ERROR ESTIMATE
!                 ERmax = ELIST(maxERR)

!        ELIST  - real
!                 VECtoR OF dimension LAST CONTAininG THE ERROR ESTIMATES

!        IORD   - integer
!                 VECtoR OF dimension LAST, THE FIRST K ELEMENTS OF
!                 WHICH CONTAin POinTERS to THE ERROR ESTIMATES,
!                 SUCH THAT ELIST(IORD(1)), ... , ELIST(IORD(K))
!                 FORM A DECREASinG SEQUENCE, WITH K = LAST if
!                 LAST <= (LIMIT/2+2), AND K = LIMIT+1-LAST OTHERWISE

!        NRmax  - integer
!                 maxERR = IORD(NRmax)

! 4.  NO subroutineS OR functionS NEEDED

!     ..................................................................


integer, intent(in)                  :: limit, last
real (dp), dimension(:), intent(in)  :: elist
integer, intent(in out)              :: nrmax
integer, dimension(:), intent(out)   :: iord
integer, intent(out)                 :: maxerr
real (dp), intent(out)               :: ermax

real (dp) :: errmax, errmin
integer   :: i, ibeg, ido, isucc, j, jbnd, jupbn, k

!           CHECK WHETHER THE LIST contains MORE THAN TWO ERROR ESTIMATES.

!***FIRST EXECUTABLE STATEMENT  QPSRT
if(last > 2) go to 10
iord(1) = 1
iord(2) = 2
go to 90

!           THIS PART OF THE RoutinE IS ONLY EXECUTED if,
!           DUE to A DifFICULT inTEGRAND, SUBDIVISION inCREASED
!           THE ERROR ESTIMATE.   in THE NORMAL case THE inSERT PROCEDURE
!           SHOULD START AFTER THE NRmax-TH LARGEST ERROR ESTIMATE.

10 errmax = elist(maxerr)
if(nrmax == 1) go to 30
ido = nrmax-1
do i = 1, ido
  isucc = iord(nrmax-1)
! ***JUMP out OF do-LOOP
  if(errmax <= elist(isucc)) EXIT
  iord(nrmax) = isucc
  nrmax = nrmax-1
end do

!           COMPUTE THE NUMBER OF ELEMENTS in THE LIST to
!           BE MAinTAinED in DESCendinG ORDER. THIS NUMBER
!           DEPendS ON THE NUMBER OF SUBDIVISIONS STILL ALLOWED.

30 jupbn = last
if(last > (limit/2+2)) jupbn = limit+3-last
errmin = elist(last)

!           inSERT ERRmax BY TRAVERSinG THE LIST toP-doWN,
!           STARTinG COMPARISON FROM THE ELEMENT ELIST(IORD(NRmax+1)).

jbnd = jupbn-1
ibeg = nrmax+1
do i=ibeg, jbnd
  isucc = iord(i)
! ***JUMP out OF do-LOOP
  if(errmax >= elist(isucc)) go to 60
  iord(i-1) = isucc
end do
iord(jbnd) = maxerr
iord(jupbn) = last
go to 90

!           inSERT ERRMin BY TRAVERSinG THE LIST BOTtoM-UP.

60 iord(i-1) = maxerr
k = jbnd
do j=i, jbnd
  isucc = iord(k)
! ***JUMP out OF do-LOOP
  if(errmin < elist(isucc)) go to 80
  iord(k+1) = isucc
  k = k-1
end do
iord(i) = last
go to 90
80 iord(k+1) = last

!           SET maxERR AND ERmax.

90 maxerr = iord(nrmax)
ermax = elist(maxerr)
return
end subroutine qpsrt



subroutine qelg (n, epstab, result, abserr, res3la, nres, epmach, oflow)
!-----------------------------------------------------------------------

! 1.  PURPOSE
!        THE RoutinE DETERMinES THE LIMIT OF A GIVEN SEQUENCE OF
!        APPROXIMATIONS, BY MEANS OF THE epsilon ALgoRITHM OF P. WYNN.
!        AN ESTIMATE OF THE absOLUTE ERROR IS ALSO GIVEN.
!        THE CONDENSED epsilon TABLE IS COMPUTED. ONLY THOSE ELEMENTS NEEDED
!        FOR THE COMPUTATION OF THE NEXT DIAgoNAL ARE PRESERVED.

! 2.  parameterS
!        N      - integer
!                 EPSTAB(N) contains THE NEW ELEMENT in THE
!                 FIRST COLUMN OF THE epsilon TABLE.

!        EPSTAB - real
!                 VECtoR OF dimension 52 CONTAininG THE ELEMENTS OF THE TWO
!                 LOWER DIAgoNALS OF THE TRIANGULAR epsilon TABLE.
!                 THE ELEMENTS ARE NUMBERED STARTinG AT THE RIGHT-HAND
!                 CORNER OF THE TRIANGLE.

!        result - real
!                 resultinG APPROXIMATION to THE inTEGRAL

!        absERR - real
!                 ESTIMATE OF THE absOLUTE ERROR COMPUTED FROM
!                 result AND THE 3 PREVIOUS resultS

!        RES3LA - real
!                 VECtoR OF dimension 3 CONTAininG THE LAST 3 resultS

!        NRES   - integer
!                 NUMBER OF CALLS to THE RoutinE
!                 (SHOULD BE ZERO AT FIRST CALL)

!        EPMACH - real
!                 THE RELATIVE PRECISION OF THE FLOATinG ARITHMETIC BEinG useD.

!        OFLOW  - real
!                 THE LARGEST POSITIVE MAGNITUDE.

! 3.  NO subroutineS OR functionS useD

!-----------------------------------------------------------------------

integer, intent(in out)                  :: n, nres
real (dp), intent(in)                    :: epmach, oflow
real (dp), intent(out)                   :: abserr, result
real (dp), dimension(:), intent(in out)  :: epstab, res3la
!---------------------

!     LIST OF MAJOR VARIABLES
!     -----------------------

!     E0     - THE 4 ELEMENTS ON WHICH THE
!     E1       COMPUTATION OF A NEW ELEMENT in
!     E2       THE epsilon TABLE IS BASED
!     E3                 E0
!                  E3    E1    NEW
!                        E2
!     NEWELM - NUMBER OF ELEMENTS to BE COMPUTED in THE NEW DIAgoNAL
!     ERROR  - ERROR = abs(E1-E0)+abs(E2-E1)+abs(NEW-E2)
!     result - THE ELEMENT in THE NEW DIAgoNAL WITH LEAST VALUE OF ERROR

!     LIMEXP IS THE maxIMUM NUMBER OF ELEMENTS THE epsilon TABLE CAN CONTAin.
!     if THIS NUMBER IS REACHED, THE UPPER DIAgoNAL OF THE epsilon TABLE IS
!     DELETED.

real (dp) :: delta1, delta2, delta3, epsinf, error, err1, err2, err3, e0, &
             e1, e1abs, e2, e3, res, ss, tol1, tol2, tol3
integer   :: i, ib, ib2, ie, indx, k1, k2, k3, limexp, newelm, num

nres = nres + 1
abserr = oflow
result = epstab(n)
if (n < 3) go to 100
limexp = 50
epstab(n + 2) = epstab(n)
newelm = (n - 1)/2
epstab(n) = oflow
num = n
k1 = n
do i = 1, newelm
  k2 = k1 - 1
  k3 = k1 - 2
  res = epstab(k1 + 2)
  e0 = epstab(k3)
  e1 = epstab(k2)
  e2 = res
  e1abs = abs(e1)
  delta2 = e2 - e1
  err2 = abs(delta2)
  tol2 = max(abs(e2),e1abs)*epmach
  delta3 = e1 - e0
  err3 = abs(delta3)
  tol3 = max(e1abs,abs(e0))*epmach
  if (err2 > tol2 .OR. err3 > tol3) go to 10

!           if E0, E1 AND E2 ARE EQUAL to WITHin MACHinE ACCURACY,
!           CONVERGENCE IS ASSUMED.
!           result = E2
!           absERR = abs(E1-E0) + abs(E2-E1)

  result = res
  abserr = err2 + err3
! ***JUMP out OF do-LOOP
  go to 100
  10 e3 = epstab(k1)
  epstab(k1) = e1
  delta1 = e1 - e3
  err1 = abs(delta1)
  tol1 = max(e1abs,abs(e3))*epmach

!           if TWO ELEMENTS ARE VERY CLOSE to EACH OTHER, OMIT
!           A PART OF THE TABLE BY ADJUSTinG THE VALUE OF N

  if (err1 <= tol1 .OR. err2 <= tol2 .OR. err3 <= tol3) go to 20
  ss = 1.0D0/delta1 + 1.0D0/delta2 - 1.0D0/delta3
  epsinf = abs(ss*e1)

!           TEST to DETECT IRREGULAR BEHAVIOUR in THE TABLE, AND EVENTUALLY
!           OMIT A PART OF THE TABLE ADJUSTinG THE VALUE OF N.

  if (epsinf > 0.1D-03) go to 30
  20 n = i + i - 1
! ***JUMP out OF do-LOOP
  go to 50

!           COMPUTE A NEW ELEMENT AND EVENTUALLY ADJUST THE VALUE OF result.

  30 res = e1 + 1.0D0/ss
  epstab(k1) = res
  k1 = k1 - 2
  error = err2 + abs(res - e2) + err3
  if (error > abserr) CYCLE
  abserr = error
  result = res
end do

!           SHifT THE TABLE.

50 if (n == limexp) n = 2*(limexp/2) - 1
ib = 1
if ((num/2)*2 == num) ib = 2
ie = newelm + 1
do i = 1, ie
  ib2 = ib + 2
  epstab(ib) = epstab(ib2)
  ib = ib2
end do
if (num == n) go to 80
indx = num - n + 1
do i = 1, n
  epstab(i) = epstab(indx)
  indx = indx + 1
end do
80 if (nres >= 4) go to 90
res3la(nres) = result
abserr = oflow
go to 100

!           COMPUTE ERROR ESTIMATE

90 abserr = abs(result - res3la(3)) + abs(result - res3la(2)) +  &
            abs(result - res3la(1))
res3la(1) = res3la(2)
res3la(2) = res3la(3)
res3la(3) = result
100 abserr = max(abserr,5.0D0*epmach*abs(result))
return
end subroutine qelg

end module adapt_quad_infinite
