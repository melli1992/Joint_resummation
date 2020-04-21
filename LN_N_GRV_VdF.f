CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C	I must remove
C	Some thousands of these logs and pile them up,
C
C	Ferdinand,The Tempest, Act 3 Scene 1.
C
C
C...DIRECT PHOTON : "DOUBLE RESUMMATION" FOR dSIGMA/dP_T
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCC
      PROGRAM MAIN
CCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE PRECISION (A-I,L-Z)

      include "int_set.inc"
      include "CUBA_init.inc"

      DIMENSION Q2THR(3), LAMBDA (3:6), QTA(7), C0A(39)
      
      INTEGER IINIP, IRES, IINIF, MP, IREC, ILNNON, ICHI, IFLAG,IPROFILE
      INTEGER intflag, bstint_cuba, reg_div, PDFINPUT,EVOLMODE, WHICHONE
      INTEGER ICONF, FLNR, FLNR2, IEV, IEVFS, CHECKPDF, KK
      INTEGER FLAV, ikar
      INTEGER QQS, QGS, H2Q, H2G
      INTEGER GAMQS, GAMGS
      INTEGER ISQS, IAR, IIAR, component
      INTEGER strlength,lfconf,lmbar,ltsmf,lrecoil,lfname,lfname2, 
     x     lsmf, linteg, lfname3,lfname6, lsmf1,
     x     ldirfrag, npoint, ii, jj, lcmp, lwhich

      double precision xx0(5), evsw, avgv,sigv,chiv,auxv, res, cmp,
     x     x_pdf, xjes(37), VARFLAG,mf1,smf1

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   NB! the following is specific for gfortran/f90 and may not work 
ccc   with other compilers
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PTA
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer offset1, offset2,pdflength
      CHARACTER*70 FNAME, FNAME2, outname
      CHARACTER*16 FCONF
      CHARACTER*6 FNAME3, SWHICHONE
      CHARACTER*8 FNAME4
      CHARACTER*5 FNAME5
      CHARACTER*10 FNAME6
      character*5 order
      character*5 evex
      character*6 pdfcfg
c      character*7 sqS
      character*9 mbar, ccmp
      character*9 smftmp, SMFtmp1
      character*7 recoil
      character*6 integ
      character*8 dirfrag, pdf
      LOGICAL EVOLIS, EVOLFS

      COMMON / RESULT / ERG1,ERG2,ERG3,ERG4
      common /resultv/ avgv,sigv,chiv,auxv
c      COMMON / SCALES / Q2, Q2START, Q2S, Q20, Q21, Q2MUR, Q2MUF, 
c     1     Q2THR, Q20F, Q21F, LAMBDA
      COMMON / SCALES / Q2, Q2START, Q20, Q21, Q2MUR, 
     1                  Q2MUF, Q2THR, LAMBDA, Q2S, Q20F, Q21F
      common / scales2 / mu2, muf2, mur2, mud2
      COMMON / FOURIER / BST, QT
      COMMON / KINE / FAC, PT, SQS
      COMMON / COUPL  / ALPS, ALP0, ALP1, ALPC, ALPB, ALPT, ALPQ, 
     1     ALPQR, LMQ, LMQR
      COMMON / COUPL1 / ALPQQ,ALPQF
      COMMON / CONT / C, CO, SI, XT2, PI, BQLIM, VP
      COMMON / CUT / MUBAR, XNP1, XNP2
      COMMON / PL / RXN0, MP
      COMMON / RESUM / IRES
      COMMON / QTINT / IPROFILE
      COMMON / RECOIL / IREC, ICHI
      COMMON / SVCAPP / ILNNON
      COMMON / FRAG / IFLAG
      COMMON / INTINIP / IINIP
      COMMON / INTINIF / IINIF      
      COMMON / FLAVORS / FLAV
      COMMON / CONFIG / EVOLSW,ICONF,FLNR,FLNR2,IEV, IEVFS
      COMMON / CONFIG2 /ISLL,ISNLL,ISLNON,ISLNON2,FSLL,FSNLL,FSLNON,FLAG
       COMMON / CONFIG3 / EVOLIS, EVOLFS,EVOLMODE 
      COMMON / CONFIG4 / IDIRECT, IHADR, REG_DIV, PDFINPUT
        COMMON / CONFIG5 / x_pdf,  WHICHONE
c mbeekveld changed QTA

      DATA QTA / 1.d0,2.d0,5.d0,7.5d0,10.0d0,20.0d0,50.0d0 /
      DATA XJES /1.d-5,2.d-5,3.5d-5,5d-5,6.d-5,7.5d-5,9.d-5,
     1           1.d-4,2.d-4,3.5d-4,5d-4,6.d-4,7.5d-4,9.d-4,
     1           1.d-3,2.d-3,3.5d-3,5d-3,6.d-3,7.5d-3,9.d-3,
     1           1.d-2,1.5d-2,2.d-2,3.5d-2,5d-2,6.d-2,7.5d-2,9.d-2,
     1           1.d-1,1.5d-1,2.d-1,3.5d-1,5d-1,6.d-1,7.5d-1,9.d-1/

c$$$      DATA C0A / 1.8d0 ,2.3D0, 3.3D0, 4.3D0, 5.3D0, 6.3D0, 7.3D0,8.3D0,
c$$$     x     9.3d0 ,10.3D0, 11.3D0, 12.3D0, 13.3D0, 14.3D0, 15.3D0,
c$$$     x     16.3D0 /
c$$$      DATA C0A / 1.5d0 ,1.75D0, 2.D0, 2.25D0, 2.5D0, 2.75D0, 3.D0,
c$$$     x     3.25D0, 3.5d0, 3.75d0, 4.d0, 4.25d0, 4.5d0, 4.75d0, 5.d0,
c$$$     x     5.25d0, 5.5d0, 5.75d0, 6.d0, 6.25d0, 6.5d0, 6.75d0, 7.d0,
c$$$     x     7.25d0, 7.5d0, 7.75d0, 8.d0, 8.25d0, 8.5d0, 8.75d0, 9.d0,
c$$$     x     9.25d0, 9.5d0, 9.75d0, 10.d0, 10.25d0, 10.5d0, 10.75d0, 
c$$$     x     11.d0 /

      EXTERNAL BSTINT, BSTINT0, bstint_cuba, PDF_CALC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   C0de begins here...
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
       AEM=1./137.
       PI=DACOS(-1.D0)
	ndim = 5
     	write(*,*) "Write 2 to check PDFs only"
	read(*,*) CHECKPDF
	iF (CHECKPDF.eq.2) THEN
	write(*,*) "whichone?"
	read(*,*)  SWHICHONE
    	call stringmod(SWHICHONE,LWHICH)

      	 read(SWHICHONE,*)WHICHONE
	endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 10    write(*,*) "======================================="
       write(*,*) "||   BLMM code for LN/N resummation  ||"
       write(*,*) "||   in prompt photon production     ||"
       write(*,*) "======================================="       
       write(*,*) "PRESS (IS radiation)"
       write(*,*) "0  for no resummation" 
       write(*,*) "1  for LL" 
       write(*,*) "2  for NLL (to muF, then resummation)"
       write(*,*) "3  for NLL (to Q/X, then resummation)" 
       write(*,*) "4  for NLL  + ln(N)/N (resummation h')" 
       write(*,*) "5  for NLL  + ln(N)/N (evolution, pure ln(N)/N)" 
       write(*,*) "6  for NLL  + O(1/N) (evolution)" 
       write(*,*) "7  for NLL  + O(1/N, 1/N^2 ...) (evolution)" 
       write(*,*) "======================================="

       READ(*,*) ICONF       

       if ((iconf.lt.0).or.(iconf.gt.7)) then
          goto 10
	else if (iconf.gt.3) then
	   IEV = ICONF
       endif
	

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ... EVOL or EXP
 61    write(*,*) "======================================="       
       write(*,*) "PRESS"
       write(*,*) "1 For no lnN/N terms in the final state "
       write(*,*) "2 for inclusion of Fragmenting lnN/N by exp.    "
       write(*,*) "3 For inclusion of Fragmenting lnN/N by diagonal ev."
       write(*,*) "4 For inclusion of Fragmenting lnN/N by full ev."
       write(*,*) "======================================="
C     Determine whether pp-, ppbar- initial state
       READ(*,*) IEVFS
       if (IEVFS.eq.1) then
	FNAME6="FS_NOLNN"
	EVOLFS=.FALSE.
	else if (IEVFS.eq.2) then
	FNAME6="FS_EXPON"
	EVOLFS=.FALSE.
	else if (IEVFS.eq.3) then
	FNAME6="FS_PARTEVOL"
	EVOLFS=.TRUE.
	else if (IEVFS.eq.4) then
	FNAME6="FS_FULLEVOL"
	EVOLFS=.TRUE.
	endif

 62    write(*,*) "======================================="       
       write(*,*) "PRESS"
       write(*,*) "1 GRV input   "
       write(*,*) "2 MMHT2014 input    "
       write(*,*) "3 MSTW2008 input    "
       write(*,*) "======================================="
C     Determine whether pp-, ppbar- initial state
       READ(*,*) PDFINPUT
       if ((PDFINPUT.lt.0) .or. (PDFINPUT.gt.3)) then
          goto 62
       endif
	if (PDFINPUT .eq. 1) THEN
		PDF = "GRV"
	elseif (PDFINPUT .eq. 2) THEN
		PDF = "MMHT"
	elseif (PDFINPUT .eq. 3) THEN
		PDF = "MSTW"
	endif
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 20    write(*,*) "======================================="       
       write(*,*) "PRESS"
       WRITE(*,*) "1  for P-P collisions    " 
       WRITE(*,*) "2  for P-Pbar collisions "
       WRITE(*,*) "3  for P-Be(4,8) collisions "
       write(*,*) "======================================="
C     Determine whether pp-, ppbar- initial state
       READ(*,*) IFLAG
       if ((iflag.lt.0).or.(iflag.gt.3)) then
          goto 20
       endif

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 30    write(*,*) "======================================="       
       write(*,*) ""
       IF(IFLAG.EQ.1)THEN
          WRITE(*,*) "1  for sqrt(S) =    63.0 GeV" 
          WRITE(*,*) "2  for sqrt(S) =  7000.0 GeV" 
          WRITE(*,*) "3  for sqrt(S) = 13000.0 GeV" 
       ELSEIF(IFLAG.EQ.2)THEN
          WRITE(*,*) "1  for sqrt(S) =    24.3 GeV"
          WRITE(*,*) "2  for sqrt(S) =  1960.0 GeV"
       ELSEIF(IFLAG.EQ.3)THEN
          WRITE(*,*) "1  for E_Beam  =   530.0 GeV"
       ENDIF
       write(*,*) "======================================="
C     Determine whether pp-, ppbar- initial state
       READ(*,*) ISQS
       if ((isqs.lt.0).or.(isqs.gt.5)) then
          goto 30
       endif

       IAR = 0

       if(isqs.eq.1.and.iflag.eq.1) then 
          FCONF = 'pp_63GeV'
          SQS =  63.0d0
          open (80, file='pT_sqrtS_63.in')
          do while (.true.)
             read (80, *, end=999) DUMMY
             IAR = IAR + 1 
          enddo 
       elseif(isqs.eq.2.and.iflag.eq.1) then 
          FCONF = 'pp_7TeV'
          SQS =  7000.0d0
          open (80, file='pT_sqrtS_7000.in')
          do while (.true.)
             read (80, *, end=999) DUMMY
             IAR = IAR + 1 
          enddo 
       elseif(isqs.eq.3.and.iflag.eq.1) then
          FCONF = 'pp_14TeV'
          SQS = 13000.0d0
          open (80, file='pT_sqrtS_13002.in')
          do while (.true.)
             read (80, *, end=999) DUMMY
             IAR = IAR + 1 
          enddo
       elseif(isqs.eq.1.and.iflag.eq.2) then 
          FCONF = 'ppbar_24.3GeV'
          SQS =    24.3d0
          open (80, file='pT_sqrtS_24-3.in')
          do while (.true.)
             read (80, *, end=999) DUMMY
             IAR = IAR + 1 
          enddo
       elseif(isqs.eq.2.and.iflag.eq.2) then
          FCONF = 'ppbar_1.96TeV'
          SQS =  1960.0d0
          open (80, file='pT_sqrtS_1960.in')
          do while (.true.)
             read (80, *, end=999) DUMMY
             IAR = IAR + 1 
          enddo
       elseif(isqs.eq.1.and.iflag.eq.3) then
          FCONF = 'pBe_Ebeam_530GeV'
          EBEAM=530.D0 
          SQS = DSQRT(2.D0*0.938D0*EBEAM)
          open (80, file='pT_sqrtS_pBe_530.in')
          do while (.true.)
             read (80, *, end=999) DUMMY
             IAR = IAR + 1 
          enddo
       endif

 999   continue

       if(sqs.ge.1000.0d0) then
          FLAV = 5
       else if(sqs.lt. 1000.0d0) then
          FLAV = 3
       endif

       REWIND(80)
       ALLOCATE(PTA(IAR))
       WRITE(*,*)"The number of pT values is: ", IAR
       DO IIAR = 1, IAR
          READ (80, *, end=888) PTA(IIAR)
          WRITE(*,*) "pt =",  pta(iiar)
       ENDDO         
       
 888      continue

       call stringmod(fconf,lfconf)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

cmbeekveld
 38    write(*,*) "======================================="       
       write(*,*) "Compute dsigma/dpTdQt [0] or"
       write(*,*) "Compute dsigma/dpT [1]"
       write(*,*) "======================================="

       read(*,*) IPROFILE
	IQTINT = IPROFILE
       if (IPROFILE.eq.0) then 
	FNAME3 = "PROFILE2"
	ndim = 3
	goto 90 
	endif
	if (IPROFILE.eq.1) then 
	FNAME3 = ""
	goto 39
	endif

cccccccccccccccccccccccccccccccccccccccccccccc

cmbeekveld
 39    write(*,*) "======================================="       
       write(*,*) "Choose way to regulate kin. div.:"
       write(*,*) "1. for mubar cutoff"
       write(*,*) "2. for contour deformation"
       write(*,*) "======================================="

       read(*,*) reg_div
        if (reg_div.eq.1) then
        goto 40 
        endif
        if (reg_div.eq.2) then 
        mbar="CONTDIF"
        call stringmod(mbar,lmbar)
        ndim = 1
        goto 90 
        endif
cccccccccccccccccccccccccccccccccc

 40    write(*,*) "======================================="       
       write(*,*) "CHOOSE THE VALUE OF MUBAR:"
       write(*,*) "NB! Especially for lower values of"
       write(*,*) "    ECM the value of MUBAR should "
       write(*,*) "    be chosen with care."
       write(*,*) "IF mubar = pT choose -1"
       write(*,*) "======================================="
       read(*,*) mbar

       call stringmod(mbar,lmbar)

       read(mbar,*) mubar
       mubar = dble(mubar)
	varflag = mubar
c       if (mubar.eq.-1) then
c       		varflag = 1
c       else
c       		varflag = 0
c       endif
       if ((mubar.lt.-1)) then
          goto 40
       endif
c       write(*,*)"!!!!!!!!!!!!!!!!!!!!!!!!!!!"
c       write(*,*) mubar
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 90    write(*,*) "======================================="       
       write(*,*) "CHOOSE THE VALUE OF SCALES RATIOS:"
       write(*,*) "(r = muR/pT)"
       write(*,*) "======================================="
C     Determine whether pp-, ppbar- initial state
c       READ(*,*) IMUBAR
       read(*,*) SMFtmp

       call stringmod(smftmp,lsmf)

       read(SMFtmp,*)smf
       smf = dble(smf)

       if ((SMF.lt.0)) then
          goto 90
       endif

c       call stringmod(TSMF,LTSMF)
c
c       read(TSMF,*) SMF
       SMF = dble(SMF)
       
 91    write(*,*) "======================================="       
       write(*,*) "CHOOSE THE VALUE OF SCALES RATIOS muF:"
       write(*,*) "(r = muF/pT)"
       write(*,*) "======================================="
C     Determine whether pp-, ppbar- initial state
c       READ(*,*) IMUBAR
       read(*,*) SMFtmp1

       call stringmod(smftmp1,lsmf1)

       read(SMFtmp1,*)smf1
       mf1 = dble(smf1)
        write(*,*) mf1
       if ((SMF.lt.0)) then
          goto 91
       endif
       
c       write(*,*)"!!!!!!!!!!!!!!!!!!!!!!!!!!!"
c       write(*,*) mubar
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...  IFRAG = 0 fragmentation OFF
C...  IFRAG = 1 fragmentation ON
       IFRAG = 1

CRB    THRESHOLD ONLY
c       IRES=2
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ... RECOIL OFF (IREC=0) OR ON (IREC=1) ?
C
C ... IRES=0: UNRESUMMED: PLAIN dSIG/dp_T at LO !
C ... IRES=1: RESUMMED
C ... IRES=2: THRESHOLD RESUMMED WITH PURE GAUSSIAN
C ... IRES=3: RESUMMED WITH (FINITE) C TERMS 
C

C mbeekveld modified for cmp
 49    write(*,*) "======================================="       
       write(*,*) "CHOOSE CMP variation:"
       write(*,*) "======================================="
       read(*,*) ccmp
       call stringmod(ccmp,lcmp)

       read(ccmp,*) cmp
       cmp = dble(cmp)

cccccccccccccccc


 50    write(*,*) "======================================="       
       write(*,*) "PRESS"
       write(*,*) "1 for no recoil effects (TR)"
       write(*,*) "2 for recoil effects    (JR)"
       write(*,*) "======================================="
C     Determine whether pp-, ppbar- initial state
       READ(*,*) IREC
       if ((irec.lt.0) .or. (irec.gt.2)) then
          goto 50
       endif

       if(irec.eq.1)then
          RECOIL = 'TR'
c          irec   = 0
          ires   = 2
       endif
       if(irec.eq.2) then
          RECOIL = 'JR'
c          irec   = 1
          ires   = 3
       endif

       irec = irec - 1

       call stringmod(recoil,lrecoil)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ... intflag=1 : VEGAS      (single core)
C ... intflag=2 : CUBA-VEGAS (multi core) 
 70    write(*,*) "======================================="       
       write(*,*) "PRESS"
       write(*,*) "1 for single core integration (VEGAS)"
       write(*,*) "2 for single core integration (VEGAN)"
       write(*,*) "3 for multi core integration (CUBA lib)"
       write(*,*) "======================================="
C     Determine whether pp-, ppbar- initial state
       READ(*,*) intflag
       if ((intflag.lt.0) .or. (intflag.gt.3)) then
          goto 70
       endif

       if (intflag.eq.1) INTEG = 'VEGAS'
       if (intflag.eq.2) INTEG = 'VEGAN'
       if (intflag.eq.3) INTEG = 'CUBA'

       call stringmod(integ,linteg)

c$$$       if (intflag.eq.1 .or. intflag.eq.2) then
c$$$ 72       write(*,*) "======================================="       
c$$$          write(*,*) "PRESS"
c$$$          write(*,*) "1 for single pT per core "
c$$$          write(*,*) "2 for all pT values read in from file "
c$$$          write(*,*) "======================================="
c$$$       endif
c$$$       
c$$$       if (intflag2.eq.1) 


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 80    write(*,*) "======================================="       
       write(*,*) "PRESS"
       write(*,*) "1 for direct photon production         "
       write(*,*) "2 for photon fragmentation only        "
       write(*,*) "3 for both                             "
       write(*,*) "======================================="
C     Determine whether pp-, ppbar- initial state
       READ(*,*) component
       if ((component.lt.0) .or. (component.gt.3)) then
          goto 80
       endif
       if(component.eq.1)then
          IDIRECT = 1.d0
          IHADR   = 0.d0
          DIRFRAG = 'DIR_ONLY'
       elseif(component.eq.2)then
          IDIRECT = 0.d0
          IHADR   = 1.d0
          DIRFRAG = 'FRA_ONLY'
       elseif(component.eq.3)then
          IDIRECT = 1.d0
          IHADR   = 1.d0
          DIRFRAG = 'DIR+FRAG'
       endif

       call stringmod(dirfrag,ldirfrag)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...  ILNNON = 0 lnN/N terms not included
C...  ILNNON = 1 include lnN/N terms via evolution
C...
c       ILNNON = 1

C...  ICHI = 1 include lnN/N terms via chi=N+b
C...  ICHI = 2 include lnN/N terms via chi=b+Nb/(1+b/4Nb)       
C...
       ICHI = 2
C
C...  EVOLUTION MATRIX BETWEEN TWO COMPLEX-VALUED SCALES, ONLY
C...  AT A*LNNBAR + B LEVEL, NO 1/N

ccc   call to subroutine configuring setup
      CALL  DEFSWITCH(
     X     ICONF,FNAME,LFNAME,
     X     IEV,IEVFS, EVOLIS,EVOLFS,EVOLSW,EVOLMODE,
     X     ISLL,ISNLL,ISLNON,ISLNON2,
     X     FSLL,FSNLL,FSLNON,ILNNON,
     X     FLAG)

	call stringmod(FNAME3, lfname3)
	call stringmod(PDF, pdflength)
	call stringmod(FNAME6, lfname6)
c mbeekveld changed the stored files
	IF(CHECKPDF.eq.2) THEN
	OPEN(FLNR, FILE = "PDFCHECK"//'_'//
     x    PDF(1:pdflength)//'_'//SWHICHONE(1:lWHICH),
     x      ACCESS = 'APPEND',STATUS = 'REPLACE')
	ELSE
       OPEN (FLNR, FILE = 'Results_2019_13TeV_lnNoN/'//
     x      FNAME3(1:lfname3)//'_'//
     x      FCONF(1:lfconf)//'_'//
     x      PDF(1:pdflength)//'_'//
     x      'mubar_'//mbar(1:lmbar)//'_'//
     x      'SMF_'//smftmp(1:lsmf)//'_'//
     x      'muFoQ_'//smftmp1(1:lsmf1)//'_'//
     x      RECOIL(1:lrecoil)//'_'//
     x      INTEG(1:linteg)//'_'//
     x      DIRFRAG(1:ldirfrag)//'_'//
     x      'CMP='//CCMP(1:lcmp)//'_'//
     x      FNAME(1:lfname)//'_'//
     x      FNAME6(1:lfname6)//'.out', 
     x      ACCESS = 'APPEND',STATUS = 'REPLACE')
	ENDIF

C
C...PT, X AND Q**2 VALUES :
C
C
       II = 1
       JJ = 1
       DO 1 WHILE (JJ.le.iar)
       PT    = PTA(JJ)  
       XT2   = ( 2.D0*PT/SQS )**2   
       Q2    =   4.D0*PT**2               
ccc   set scales
c$$$       MU2   = Q2
c$$$       MUF2  = Q2
c$$$       MUR2  = Q2
c$$$       MUD2  = Q2
       MF    = SMF**2
       MU2   = MF * PT**2
       MUF2  = MF1**2* PT**2
       MUR2  = MF * PT**2
       MUD2  = MF * PT**2
cccccccccccccccccccccccccccccc

       write(*,*)"iar = ", iar 
       IF (varflag < 0) THEN
            mubar = -1*PT*varflag
            
        ENDIF
	write(*,*) "mubar", mubar
c mbeekveld
c		if (varflag.eq.1) then
c        mubar = PT
c        endif


C
C ... IF QT INTEGRATION ACTIVE, CHOOSE CUT ON QT:
C       MUBAR = dble(IMUBAR)

C 
C ... NON-PERTURBATIVE PARAMETER :
C ... ITS MEANING IS SQRT(AVERAGE VALUE OF KT**2)
       XNP1  = 1.0D0
       XNP2  = 0.5D0
C
C ... PARAMETERS FOR LOGARITHM :
       EMC   = 0.5772156649D0
       MP    = 1
       DMP   = DBLE(MP) 
       RXN0  = 1.D0

C
C ... CHOICE OF SCALES :
C

       Q2MUF = MF1**2*PT**2
       Q2MUR = MF*PT**2
C
        write(*,*)"MF       =", mf
        write(*,*)"MF1       =", mf1
       write(*,*)"q2       = ", q2
       write(*,*)"q2muf    = ", q2muf
       write(*,*)"q2mur    = ", q2mur

       LMQ   = LOG(Q2/Q2MUF)
       LMQR  = LOG(Q2/Q2MUR)
	IF (PDFINPUT == 1) THEN
       INPUT = ALPHAS (1.d0)*4*PI
       CALL INITALPHAS(1, 1.d0, 1.d0, INPUT, 1.4d0, 4.75d0, 1.73D10)
	ELSE IF (PDFINPUT == 2) THEN
       CALL INITALPHAS(1, 1.d0, 1.d0, 0.49128d0, 1.4d0, 4.75d0, 1.73D10)
	ELSE IF (PDFINPUT == 3) THEN
       CALL INITALPHAS(1, 1.d0, 1.d0, 0.49128d0, 1.4d0, 4.75d0, 1.73D10)
	ENDIF

       ALPS  = ALPHAS_MSTW(Q2START**0.5)/(4*PI)
       ALP0  = ALPHAS_MSTW(Q20**0.5)/(4*PI)
       ALP1  = ALPHAS_MSTW(Q21**0.5)/(4*PI)
       ALPC  = ALPHAS_MSTW(Q2THR(1)**0.5)/(4*PI)
       ALPB  = ALPHAS_MSTW(Q2THR(2)**0.5)/(4*PI)
       ALPT  = ALPHAS_MSTW(Q2THR(3)**0.5)/(4*PI)
       ALPQ  = ALPHAS_MSTW(Q2MUF**0.5)/(4*PI)       
       ALPQR = ALPHAS_MSTW(Q2MUR**0.5)/(4*PI)
       ALPQF = ALPQ
       ALPQQ = ALPHAS_MSTW(Q2**0.5)/(4*PI)
       write(*,*)"------ alpha(Q^2) ------"
       write(*,*)"lmq      = ", lmq
       write(*,*)"lmqr     = ", lmqr

       write(*,*)"alps     = ", alps
       write(*,*)"alp0     = ", alp0
       write(*,*)"alp1     = ", alp1
       write(*,*)"alpc     = ", alpc
       write(*,*)"alpb     = ", alpb
       write(*,*)"alpt     = ", alpt
       write(*,*)"alpq     = ", alpq
       write(*,*)"alpqr    = ", alpqr
       write(*,*)"alpqf    = ", alpqf
       write(*,*)"alpqq    = ", alpqq

       write(*,*)"q2start  = ", q2start
       write(*,*)"q20      = ", q20
       write(*,*)"q21      = ", q21
       write(*,*)"q2thr(1) = ", q2thr(1)
       write(*,*)"q2thr(2) = ", q2thr(2)
       write(*,*)"q2thr(3) = ", q2thr(3)
       write(*,*)"Q2       = ", q2
       write(*,*)"PT       = ", pt
       write(*,*)"----------------------------"

C
c conversion factor between GEV-2 and pb, end result is in pb/GeV
c mbeekveld modified this since is was wrong!
       FAC1  = AEM*3.893793656043442025902826102928969744d8
       FAC   = FAC1*4.D0*PI*ALPQR/PT**3 * XT2**2
c 4*PI*ALPQR = alpas in notes
c XT2**2 = 16*pT**4/S**2
c so indeed, dsigma/dpT is calculated (pT**3 is divided out!)
c and the partonic functions are given by 1/16*Integrate[(4*v*(1-v)w)^N* s/2*dsigma/(dv dw),{v,0,1}]
c = 1/16*1/(8*Pi)*Integrate[xT2^N*<|M|^2>/sqrt(1-xT2),{xT2,0,1}]


C...PARAMETERS OF THE INTEGRATION CONTOURS IN THE COMPLEX PLANE :
C
c       PHI = PI * 3./4.
c       PHI = PI * 5./6.
c mbeekveld
	IF(REG_DIV.eq.2) THEN
       PHI   = PI * 3.d0/4.d0
	
	ELSE 
       PHI   = PI * 25.d0/32.d0
	
	ENDIF
       CO    = COS (PHI)
       SI    = SIN (PHI)       
c     C0 = C0A(JJ)
	    IF(REG_DIV.eq.2) THEN
       CADD  = 1.D0
       ELSE
       CADD  = -3.D0/LOG(XT2)
	ENDIF
C
       npoint = 50000

       if(sqs.gt.1000.d0) c0 = 1.8d0
       if(sqs.le.1000.d0) c0 = 2.3d0

c       do 2 ii=1,19
c          c0 = c0a(ii)
       IF(IRES.EQ.0) THEN
         C     = C0-1.D0 + CADD+DBLE(cmp)
       ELSE IF(IRES.GE.1) THEN
         B0    = B0CALC(Q2MUR,Q2THR)       
C .... LANDAU SINGULARITY: 
         XNLAN = EXP(1.D0/B0/ALPQR/4.D0/PI/2.D0 - EMC)
         XLANL = 0.7D0 * (XNLAN-C0) + C0
         C     = C0-1.0D0 + CADD/(1.D0 + CADD/(XLANL-C0))+DBLE(cmp)
       ENDIF

       print*, "main: xnlan = ", xnlan
       print*, "main: xlanl = ", xlanl
       print*, "main: c     = ", c
       print*, "main: c0    = ", c0
       print*, "main: cadd  = ", cadd
C
C ... BRANCHING POINT FOR B CONTOUR : 
C ... NOTE, BQLIM HAS TO BE SMALLER THAN, BUT CLOSE TO, XNLAN-C   
C ... BQLIM IS BLIM*Q/2 
C
       BQLIM = 1.D0

C
C ... CONTOUR PARAMETER V (IN UNITS OF PI) :
       VP = 2.D0
C 
       WRITE(6,*) XNLAN, C, BQLIM
C 
       
       
       if(IPROFILE.eq.0) then
        mubar = QTA(II)
        II = II + 1
        JJ = JJ - 1
        if (II.eq.7) then
            II = 1
            JJ = JJ + 1
        endif
        endif
        write(*,*) "MUBAR", mubar
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       write(*,*) "======================================="       
       write(*,*) "||   Configuration output            ||"
       write(*,*) "======================================="       
       write(*,*) "     Flags for physical config.        "
       write(*,*) "IRES    = ", IRES
       write(*,*) "ICONF   = ", ICONF
       write(*,*) "ILNNON  = ", ILNNON
       write(*,*) "IREC    = ", IREC
       write(*,*) "ICHI    = ", ICHI
       write(*,*) "IFLAG   = ", IFLAG
       write(*,*) "IQTINT  = ", IQTINT
       write(*,*) " "
       write(*,*) "     Physical variables                "
       write(*,*) "Q2      = ", Q2
       write(*,*) "MU2/Q2  = ", SMF
       write(*,*) "MUBAR   = ", MUBAR
       write(*,*) "Q2MUF   = ", Q2MUF
       write(*,*) "Q2MUR   = ", Q2MUR
       write(*,*) "======================================="       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C ... CROSS SECTION, NORMALIZATION IS dSIG/dp_T :

        write(*,*) "ndim", ndim
        if (CHECKPDF.eq.2) then
            goto 102
        else
            continue
        endif
       if(intflag.eq.1)then
ccc   "plain" VEGAS 
          CALL VEGAS(BSTINT,1.D-6,ndim,npoint,8,1,0)
          SIG = ERG1
          STD = ERG2   
       elseif(intflag.eq.2) then          
ccc   VEGAS with w grid adaptions and 4 integration iterations
          res = intveg(bstint,4+iqtint,npoint,4)
          SIG = avgv
          STD = sigv
          CHI = chiv
       elseif(intflag.eq.3) then
ccc   CUBA-VEGAS multi-core integration
          call cuba(1, ndim, ncomp, bstint_cuba, integral, error, prob)
c          write(*,*)"Currently disabled."
c          stop
          sig = integral(1)
          std = error(1)
       endif

       JJ = JJ + 1
ccc   write output to file
        IF (IPROFILE.EQ.1) THEN
       WRITE(FLNR,13) PT, SIG, STD, C0 
C
  13   FORMAT (4(1PE13.5))
        ENDIF
       IF(IPROFILE.EQ.0) THEN
       WRITE(FLNR,13) MUBAR, PT, SIG, STD
C
  14   FORMAT (5(1PE13.5))
        ENDIF

c  2    CONTINUE
c  3    CONTINUE
	
102	IF (CHECKPDF.eq.2) THEN 
		KK = 1
		ikar = 37
       		DO 101 WHILE (KK.le.ikar)
		write(*,*) xjes(KK)
		x_pdf = xjes(KK)
		CALL VEGAS(PDF_CALC, 1.D-10, 1, npoint, 8, 1, 0)
          	SIG = ERG1
          	STD = ERG2   
    		WRITE(FLNR,13) x_pdf, SIG, STD, Q2
		KK= KK + 1
101		CONTINUE
	 	JJ = JJ + 1
		
	ENDIF
  1    CONTINUE
C
       END


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	FUNCTION PDF_CALC(XX, WGT)
	IMPLICIT DOUBLE COMPLEX (A-Z)
	DOUBLE PRECISION PI, CO, SI, PHI,Z,Y, C, x_pdf, WGT,PDF_CALC
	DIMENSION XX(1)
	INTEGER whichone
	LOGICAL misnan
        external misnan
        COMMON / CONFIG5 / x_pdf, whichone
        PI   = DACOS(-1.D0)
        PHI   = PI * 3.d0/4.d0 
	Z       = XX(1)
        Y      = Z/(1.D0-Z)
	C =1.2d0
	EXPO = DCMPLX(COS(PHI), SIN(PHI))
	XN = DCMPLX(C + COS(PHI)*Y+1.D0, SIN(PHI)*Y)
c	write(*,*) Y
c	x_pdf = 0.5D0
C       
        CALL RENO (VAN, NS3M, NS8M, NS3N,NS8N,NS15N, NS15M,
     1       NS35N, NS35M, NS24N, NS24M, SINTMP, GLN, XN) 
        THRD = DCMPLX(1.D0/3.D0,0.D0)
        TVN = (VAN - NS35M) * 0.5D0 * THRD
        BVN = TVN + (NS35M - NS24M) * 0.2D0  
        CVN = BVN + (NS24M - NS15M) * 0.25D0 
        SVN = CVN + (NS15M - NS8M) * THRD 
        DVN = SVN + (NS8M - NS3M) * 0.5D0 
        UVN = SVN + (NS8M + NS3M) * 0.5D0 
        TP = (SINTMP - NS35N) * 0.5D0 * THRD
        BP = TP + (NS35N - NS24N) * 0.2D0
        CP = BP + (NS24N - NS15N) * 0.25D0
        SP = CP + (NS15N - NS8N) * THRD
        DP = SP + (NS8N - NS3N) * 0.5D0
        UP = SP + (NS8N + NS3N) * 0.5D0
        TSN = 0.5D0*(TP-TVN)
        BSN = 0.5D0*(BP-BVN)
        CSN = 0.5D0*(CP-CVN)
        SSN = 0.5D0*(SP-SVN)
        DSN = 0.5D0*(DP-DVN)
        USN = 0.5D0*(UP-UVN)
	IF (WHICHONE.eq.1) THEN
	RES = UVN
	ELSEIF (WHICHONE.eq.2) THEN
	RES = DVN
	ELSEIF (WHICHONE.eq.3) THEN
	RES = SVN
	ELSEIF (WHICHONE.eq.4) THEN
	RES = USN
	ELSEIF (WHICHONE.eq.5) THEN
	RES = DSN
	ELSEIF (WHICHONE.eq.6) THEN
	RES = SSN
	ELSEIF (WHICHONE.eq.7) THEN
	RES = CSN
	ELSEIF (WHICHONE.eq.8) THEN
	RES = BSN
	ELSEIF (WHICHONE.eq.9) THEN
	RES = GLN
	ENDIF
	PDF_CALC = DIMAG(EXPO*x_pdf**(-XN+1.D0)*RES)
     1       *1.0D0/(PI*(1.0D0-Z)**2.0D0)
c	PDFCALC = REAL(1.D0-XX(1))
 	if(misnan(PDF_CALC) ) then
		PDF_CALC = 0.d0
		WGT = 0.d0
	endif
c	write(*,*) PDF_CALC
	RETURN
	END



       SUBROUTINE CDGRES2(ALPSMU,F,XN,LMQ,LMQR,BBQ,CCQQ,CCQG,CCQQ1
     1     ,CCQQ2,CCQQ3,CCQQ4,CCQQ5,CCQQ6,CCQG1,CCQG2,CCGG1,CCGG2
     1     ,ENS_GRES)
       IMPLICIT DOUBLE COMPLEX (A-Z)
       DOUBLE PRECISION PI,ZET2,CF,CA,TF,AQ1,AQ2,BQ1,B0,B1,EMC,LMQ,
     1      LMQR,ALPSMU,ALPQQ,ALPQF,ALRAT,DMP, AKQ,
     2      Q2, Q2START, Q20, Q21, Q2MUR, Q2MUF, 
     3      Q2THR(3), LAMBDA(3:6), MUBAR, AK, RXN0,
     4      AG1,AG2,BG1, XNP1, XNP2, CFHQ, CFHG, 
     5      CFHQ1,CFHQ2,CFHQ3,CFHQ4,CFHQ5,CFHQ6, CFHG1,
     6      CFHG2,CFHG3,CFHG4, q2s, q20f, q21f,
     7      ISLL, ISNLL, ISLNON, ISLNON2, FSLL, FSNLL, FSLNON, 
     8      FLAG, isthreshold
ccc   external functions
       DOUBLE PRECISION c1_qqp_qqp, c1_qqbp_qqbp, c1_qqb_qpqbp, 
     x      c1_qq_qq, c1_qqb_qqb, c1_qg_qg, c1_qg_gq, c1_gg_qqb
     x      ,q2tst, mu2, muf2, mur2, mud2,IDIRECT, EVOLSW,IHADR
       EXTERNAL c1_qqp_qqp, c1_qqbp_qqbp, c1_qqb_qpqbp, c1_qq_qq, 
     x      c1_qqb_qqb, c1_qg_qg, c1_qg_gq, c1_gg_qqb

       INTEGER ICONF,FLNR,FLNR2,IEV, IEVFS, REG_DIV
       INTEGER MP, F, IRES, ILNNON, IREC, ICHI
       INTEGER HS,FS,FS2, PDFINPUT
       INTEGER QQS, QGS, H2Q, H2G
       INTEGER GAMQS,GAMGS,EVOLMODE 
       INTEGER FLAV, NNFF
       LOGICAL EVOLIS, EVOLFS

       COMMON / RESUM / IRES
       COMMON / COUPL1 / ALPQQ,ALPQF
       COMMON / CUT / MUBAR, XNP1, XNP2
       COMMON / PL / RXN0, MP
c       COMMON / SCALES / Q2, Q2START, Q2S, Q20, Q21, Q2MUR, Q2MUF,
c     1                   Q2THR, Q20F, Q21F, LAMBDA 
       COMMON / SCALES / Q2, Q2START, Q20, Q21, Q2MUR, 
     1                  Q2MUF, Q2THR, LAMBDA, Q2S, Q20F, Q21F
       common / scales2 / mu2, muf2, mur2, mud2
       COMMON / RECOIL / IREC, ICHI
       COMMON / SVCAPP / ILNNON
      COMMON / CONFIG / EVOLSW,ICONF,FLNR,FLNR2,IEV, IEVFS
       COMMON / CONFIG2 /ISLL,ISNLL,ISLNON,ISLNON2,
     x      FSLL,FSNLL,FSLNON,FLAG
       COMMON / CONFIG3 / EVOLIS, EVOLFS,EVOLMODE 
      COMMON / CONFIG4 / IDIRECT, IHADR, REG_DIV, PDFINPUT
c       COMMON / FLAVORS / FLAV

CCC   Number of flavors hard wired for now...
       NNFF = 5

       PI   = DACOS(-1.D0)
       ZET2 = PI**2/6.D0
       CF   = 4.D0/3.D0
       CA   = 3.D0
c       TF   = DBLE(F)/2.D0
       TF   = DBLE(NNFF)/2.D0

       B0   = (11.D0*CA-4.D0*TF)/12.D0/PI
       B1   = (17.D0*CA**2-10.D0*CA*TF-6.D0*CF*TF)/24.D0/PI**2
       EMC  = 0.5772156649D0

       AQ1  = CF
       AG1  = CA
       AK   = CA*(67.D0/18.D0-ZET2) - 10.D0/9.D0*TF
       AKQ  = CF * (7.D0/2.D0-ZET2)
       AQ2  = 0.5D0*CF*AK
       AG2  = 0.5D0*CA*AK
       BQ1  = -3.D0/2.D0 * CF
       BG1  = -2.D0*PI * B0


CAM ... color weights as in eqs (A.1)-(A.10) in PRD71,114004
       G11  =  1.D0/ 3.D0
       G12  =  2.D0/ 3.D0
       G21  =  1.D0/ 9.D0
       G22  =  8.D0/ 9.D0
       G31  =  1.D0
       G32  =  0.D0
       G41  =  9.D0/11.D0
       G42  =  2.D0/11.D0
       G51  =  5.D0/21.D0
       G52  = 16.D0/21.D0
       G61  =  5.D0/ 7.D0
       G62  =  2.D0/ 7.D0
       G71  = 45.D0/88.D0
       G72  = 25.D0/88.D0
       G73  = 18.D0/88.D0
       G81  = 45.D0/88.D0
       G82  = 25.D0/88.D0
       G83  = 18.D0/88.D0
       G91  =  1.D0/ 3.D0
       G92  =  1.D0/ 2.D0
       G93  =  1.D0/ 6.D0
       G101 =  5.D0/ 7.D0
       G102 =  2.D0/ 7.D0
       G103 =  0.D0
       
       XNBAR  = XN*EXP(EMC)
       L      = B0*ALPSMU*LOG(XNBAR) * 4.D0*PI
c mbeekveld if we do deformation, we don't need beta
        IF (REG_DIV.eq.2) THEN
        BETA = L
        ELSE 
       CBQ    = BBQ * EXP(EMC)
       IF(ICHI.EQ.1) THEN
          CHI = XN + BBQ
       ELSEIF(ICHI.EQ.2) THEN
          CHI = CBQ + XNBAR/(1.D0 + CBQ/4.D0/XNBAR)
       ENDIF
       BETA   = B0*ALPSMU*LOG(CHI) * 4.D0*PI
          ENDIF
	ISTHRESHOLD = 1.D0
       IF (IRES.EQ.2) THEN 
C ... (Threshold resummation)     
          BETA = L
	  ISTHRESHOLD = 0.D0
       ENDIF
ccc   Kulesza-Sterman-Vogelsang exponents
c mbeekveld if reg_div is used, we need to modify this 
       IF (REG_DIV.eq.2) THEN
       QIN       =  ISLL   * H0L(L,L,B0,AQ1)/ALPSMU/4.D0/PI 
     1            + ISNLL  * H1L(L,L,B0,B1,AQ1,AQ2,BQ1,LMQ,LMQR) 
     2            + ISLNON * HPRIME(L,B0,AQ1,XN)
     3	+ 1.644934*ALPSMU*4.D0/2.D0*AQ1/(1.D0-2.D0*L)*ISTHRESHOLD
     4	+ ALPSMU*4.D0/2.D0*AQ1*LOG(XNBAR)/XN*ISTHRESHOLD
c mbeekveld QIN had an error in ISTHRESHOLD!
     
	
ccc
       GIN       =  ISLL   * H0L(L,L,B0,AG1)/ALPSMU/4.D0/PI 
     1            + ISNLL  * H1L(L,L,B0,B1,AG1,AG2,BG1,LMQ,LMQR)
c mbeekveld wrong     2            + ISLNON * HPRIMEIS(L,L,B0,3.d0*AG1,ALPSMU*4.D0*PI)    
     2            + ISLNON * HPRIME(L,B0,AG1,XN)    
     3 + 1.644934*ALPSMU*4.D0/2.D0*AG1/(1.D0-2.D0*L)*ISTHRESHOLD 
     4	+ ALPSMU*4.D0/2.D0*AG1*LOG(XNBAR)/XN*ISTHRESHOLD

	    ElSE
       QIN       =  ISLL   * H0L(L,BETA,B0,AQ1)/ALPSMU/4.D0/PI 
     1            + ISNLL  * H1L(L,BETA,B0,B1,AQ1,AQ2,BQ1,LMQ,LMQR) 
     2            + ISLNON * HPRIME(BETA,B0,AQ1,XN)
	
ccc
       GIN       =  ISLL   * H0L(L,BETA,B0,AG1)/ALPSMU/4.D0/PI 
     1            + ISNLL  * H1L(L,BETA,B0,B1,AG1,AG2,BG1,LMQ,LMQR)
cc mbeekveld this is wrong     2            + ISLNON * HPRIMEIS(L,BETA,B0,3.d0*AG1,ALPSMU*4.D0*PI) 
     2            + ISLNON * HPRIME(BETA,B0,AG1,XN)     
        ENDIF
ccc
ccc   
ccc   QINFRAG, GINFRAG are the exponents for the fragmenting partons
ccc   We are not including recoil effects for the fragmenting parton
       QINFRAG   =  ISLL   * H0L_NOEV(L,L,B0,AQ1)/ALPSMU/4.D0/PI 
     1            + ISNLL  * H1L_NOEV(L,L,B0,B1,AQ1,AQ2,BQ1,LMQ,LMQR) 
     2            + ISLNON2* HPRIME(L,B0,AQ1,XN)
ccc
       GINFRAG   =  ISLL   * H0L_NOEV(L,L,B0,AG1)/ALPSMU/4.D0/PI 
     1            + ISNLL  * H1L_NOEV(L,L,B0,B1,AG1,AG2,BG1,LMQ,LMQR)
cmbeekveld this is wrong (factor 3):     2            + ISLNON2* HPRIME(L,B0,3.d0*AG1,ALPSMU*4.D0*PI)     
     2            + ISLNON2* HPRIME(L,B0,AG1,XN)      
CCC   --- ORIGINALS ---
c$$$       QINFRAG   =  ISLL   * H0L_NOEV(L,BETA,B0,AQ1)/ALPSMU/4.D0/PI 
c$$$     1            + ISNLL  * H1L_NOEV(L,BETA,B0,B1,AQ1,AQ2,BQ1,LMQ,LMQR) 
c$$$     2            + ISLNON2* HPRIME(L,B0,AQ1,ALPSMU*4.D0*PI)
c$$$ccc
c$$$       GINFRAG   =  ISLL   * H0L_NOEV(L,BETA,B0,AG1)/ALPSMU/4.D0/PI 
c$$$     1            + ISNLL  * H1L_NOEV(L,BETA,B0,B1,AG1,AG2,BG1,LMQ,LMQR)
c$$$     2            + ISLNON2* HPRIME(L,B0,3.d0*AG1,ALPSMU*4.D0*PI)      
ccc
c       GIN  = GIN * FLAG 

CEL Drop the F1L as this is LL for both initial and final
CEL Note that H1L has an initial argument L again, for possible scale 
CEL logs
    
       QOUT =  FSLL   * F0L(L,B0,AQ1)/ALPSMU/4.D0/PI 
     1       + FSNLL  * F1L(L,B0,B1,AQ1,AQ2,BQ1,LMQR)
     2       + FSLNON * F2L(L,B0,AQ1,XN)
c       QOUT = QOUT * FLAG 
ccc 
       GOUT =  FSLL   * F0L(L,B0,AG1)/ALPSMU/4.D0/PI 
     1       + FSNLL  * F1L(L,B0,B1,AG1,AG2,BG1,LMQR)
ccc     2       + FSLNON * F2L(L,B0,AG1,ALPSMU) * 3.d0 mbeekveld factor 3 is wrong
     2       + FSLNON * F2L(L,B0,AG1,XN)
c       GOUT = GOUT * FLAG 



c$$$       write(*,*) lmq, lmqr, 
c$$$     x      qin, gin, qinfrag, ginfrag, qout, gout
c$$$     x      HPRIME(L,B0,AQ1,ALPSMU*4.D0*PI),
c$$$     x      HPRIME(L,B0,3.d0*AG1,ALPSMU*4.D0*PI),   
c$$$     x      F2L(L,B0,AQ1,ALPSMU),
c$$$     x      F2L(L,B0,AG1,ALPSMU) * 3.d0,
c$$$     x      H1L(L,BETA,B0,B1,AQ1,AQ2,BQ1,LMQ,LMQR),
c$$$     x      H1L(L,BETA,B0,B1,AG1,AG2,BG1,LMQ,LMQR),
c$$$     x      H1L_NOEV(L,L,B0,B1,AQ1,AQ2,BQ1,LMQ,LMQR),
c$$$     x      H1L_NOEV(L,L,B0,B1,AG1,AG2,BG1,LMQ,LMQR),
c$$$     x      F1L(L,B0,B1,AQ1,AQ2,BQ1,LMQR),
c$$$     x      F1L(L,B0,B1,AG1,AG2,BG1,LMQR)

       GHELP   =   1.D0/PI/B0 *LOG(1.-2.*L) *              LOG(2.D0) 
       GAMQQ   = - AG1 * GHELP
       GAMQG   = - AQ1 * GHELP
       GAMQ11  = - 1./PI/B0/2.*LOG(1.-2.*L) *       4.D0 * LOG(2.D0)
       GAMQ12  =   1./PI/B0/2.*LOG(1.-2.*L) *       0.D0
       GAMQ21  = - 1./PI/B0/2.*LOG(1.-2.*L) * 10.D0/3.D0 * LOG(2.D0)
       GAMQ22  =   1./PI/B0/2.*LOG(1.-2.*L) *  8.D0/3.D0 * LOG(2.D0)
       GAMQ31  = - 1./PI/B0/2.*LOG(1.-2.*L) *      10.D0 * LOG(2.D0)
       GAMQ32  =   0.D0
       GAMQ41  = - 1./PI/B0/2.*LOG(1.-2.*L) *       4.D0 * LOG(2.D0)
       GAMQ42  = - 1./PI/B0/2.*LOG(1.-2.*L) *       0.D0 * LOG(2.D0)
       GAMQ51  = - 1./PI/B0/2.*LOG(1.-2.*L) * 10.D0/3.D0 * LOG(2.D0)
       GAMQ52  =   1./PI/B0/2.*LOG(1.-2.*L) *  8.D0/3.D0 * LOG(2.D0)
       GAMQ61  = - 1./PI/B0/2.*LOG(1.-2.*L) * 10.D0/3.D0 * LOG(2.D0)
       GAMQ62  =   1./PI/B0/2.*LOG(1.-2.*L) *  8.D0/3.D0 * LOG(2.D0)
       GAMQ71  = - 1./PI/B0/2.*LOG(1.-2.*L) * 14.D0/3.D0 * LOG(2.D0)
       GAMQ72  =   1./PI/B0/2.*LOG(1.-2.*L) * 10.D0/3.D0 * LOG(2.D0)
       GAMQ73  = - 1./PI/B0/2.*LOG(1.-2.*L) *  2.D0/3.D0 * LOG(2.D0)
       GAMQ81  = - 1./PI/B0/2.*LOG(1.-2.*L) *       8.D0 * LOG(2.D0)
       GAMQ82  = - 1./PI/B0/2.*LOG(1.-2.*L) *       0.D0 * LOG(2.D0)
       GAMQ83  = - 1./PI/B0/2.*LOG(1.-2.*L) *       4.D0 * LOG(2.D0)
       GAMQ91  = - 1./PI/B0/2.*LOG(1.-2.*L) *       0.D0 * LOG(2.D0)
       GAMQ92  = - 1./PI/B0/2.*LOG(1.-2.*L) *      10.D0 * LOG(2.D0)
       GAMQ93  =   1./PI/B0/2.*LOG(1.-2.*L) *       6.D0 * LOG(2.D0)
       GAMQ101 = - 1./PI/B0/2.*LOG(1.-2.*L) *       0.D0 * LOG(2.D0)
       GAMQ102 =   1./PI/B0/2.*LOG(1.-2.*L) *       6.D0 * LOG(2.D0)

       GHELP   =   GHELP   * FSNLL
       GAMQQ   =   GAMQQ   * FSNLL
       GAMQG   =   GAMQG   * FSNLL
       GAMQ11  =   GAMQ11  * FSNLL
       GAMQ12  =   GAMQ12  * FSNLL
       GAMQ21  =   GAMQ21  * FSNLL
       GAMQ22  =   GAMQ22  * FSNLL
       GAMQ31  =   GAMQ31  * FSNLL
       GAMQ32  =   GAMQ32  * FSNLL
       GAMQ41  =   GAMQ41  * FSNLL 
       GAMQ42  =   GAMQ42  * FSNLL
       GAMQ51  =   GAMQ51  * FSNLL
       GAMQ52  =   GAMQ52  * FSNLL
       GAMQ61  =   GAMQ61  * FSNLL 
       GAMQ62  =   GAMQ62  * FSNLL
       GAMQ71  =   GAMQ71  * FSNLL
       GAMQ72  =   GAMQ72  * FSNLL 
       GAMQ73  =   GAMQ73  * FSNLL
       GAMQ81  =   GAMQ81  * FSNLL
       GAMQ82  =   GAMQ82  * FSNLL
       GAMQ83  =   GAMQ83  * FSNLL
       GAMQ91  =   GAMQ91  * FSNLL
       GAMQ92  =   GAMQ92  * FSNLL
       GAMQ93  =   GAMQ93  * FSNLL
       GAMQ101 =   GAMQ101 * FSNLL
       GAMQ102 =   GAMQ102 * FSNLL

c$$$       TMPC    = EXP(AQ1/2.D0/PI/B0**2 * 
c$$$     1      ( -2.D0*L*LOG(1.D0-2.D0*BETA) /ALPSMU/4.D0/PI ))

       IF ((ISNLL.gt.0) .and. (FSNLL.gt.0)) THEN 
CCC   IF NLL, NLL + lnN/N

CCC       IF ((IRES.EQ.2) .OR. (IRES.EQ.3)) THEN 
C
C ... FINITE DELTA-FCT. PART :
C
CEL Set LMQ to zero for the case of using evolmat. These terms
CEL are including the anomalous dimension already.

ccc   testing
c          evolsw =1.d0
ccccccccccccc
CPM   CFHQ and CFHG are taken from JHEP 9903:025,1999
cmbeekveld not exactly the same as in there, but we use Nbar instead of 
c N and we use a different scale, which gives the ln(2) probably
c https://arxiv.org/pdf/hep-ph/9903436.pdf
c if evolsw = 0 then we evolve the inital state to Q^2!
c https://arxiv.org/pdf/hep-ph/9806484.pdf there is explaination sec 4.3
c here are the matching coefficients
         CFHQ   = - 0.5D0*(2.d0*CF-CA)*DLOG(2.D0) 
     1            + AK/2.D0 - AKQ + 2.d0*ZET2*(2.d0*CF-CA/2.d0) 
     2            + 5.d0/4.d0*(2.d0*CF-CA)*DLOG(2.D0)**2 
     3            - PI*B0*(LMQR-DLOG(2.D0)) 
     4            + 3.d0/2.d0*CF*(-DLOG(2.D0))  
     5            + 3.d0/2.d0*CF*LMQ*EVOLSW
c$$$         write(*,*)"CFHQ    = ", cfhq
CPM          CFHQ = 1.D0
C
         CFHG   = - 0.1D0*(CF-2.d0*CA)*DLOG(2.D0) 
     1            - AKQ/2.d0 + ZET2/10.d0*(2.d0*CF+19.d0*CA) 
     2            + 1.d0/2.d0*CF*DLOG(2.D0)**2 
     3            - PI*B0*(LMQR-DLOG(2.D0)) 
     4            + (3.d0/4.d0*CF+PI*B0)*(-DLOG(2.D0))
     5            + (3.d0/4.d0*CF+PI*B0)*LMQ*EVOLSW
c$$$       write(*,*)"CFHG    = ", cfhg
c$$$       write(*,*)"EVSW    = ", evolsw
c$$$       WRITE(*,*)"LMQ     = ", LMQ
c$$$       WRITE(*,*)"LMQR    = ", LMQR
c$$$       WRITE(*,*)"ALPSMU  = ", alpsmu
CPM          CFHG = 1.D0
CCC
         CFHQ1  = 1.D0 + 4.D0 * ALPSMU * CFHQ 
         CFHG1  = 1.D0 + 4.D0 * ALPSMU * CFHG

CAM .... CFHQ11 etc as in (A.1)-(A.10)in PRD71,114004,2005.
ccc   q qp  --> q qp
         CFHQ11 = 1.D0 + 4.d0*ALPSMU * 20.2389d0
c         write(*,*)"(1)  CFHQ11   = ", CFHQ11
ccc   q qbp --> q qbp
         CFHQ12 = 1.D0 + 4.d0*ALPSMU * 21.3772d0
ccc   q qb  --> qp qbp
         CFHQ13 = 1.D0 + 4.d0*ALPSMU * 7.91881d0
ccc   q q   --> q q
         CFHQ14 = 1.D0 + 4.d0*ALPSMU * 19.5350d0
c         write(*,*)"(1)  CFHQ14   = ", CFHQ14
ccc   q qb  --> q qb
         CFHQ15 = 1.D0 + 4.d0*ALPSMU * 19.9643d0
ccc   q qb  --> g g (omit)
         CFHQ16 = 1.D0 + 4.d0*ALPSMU * 12.4329d0
ccc   q g   --> q g: q fragments, g --> jet
         CFHG11 = 1.D0 + 4.d0*ALPSMU * 15.4167d0
cccccccccccccccccccccccccc
ccc   q g   --> g q: g fragments, q --> jet, hence EXCLUDED
         CFHG12 = 1.D0 + 4.d0*ALPSMU * 22.4474d0
cccccccccccccccccccccccccc
ccc   g g   --> g g (omit)
         CFHG13 = 1.D0 + 4.d0*ALPSMU * 21.1977d0
cccccccccccccccccccccccccc
ccc   g g   --> q qb
         CFHG14 = 1.D0 + 4.d0*ALPSMU * 16.7862d0
c$$$
c$$$         write(*,*)"(1)CFHQ13   = ", CFHQ13
c$$$         write(*,*)"(1)raw walue= ", 7.91881d0
c$$$c         write(*,*)"1 CFHQ11  =   ", CFHQ11
ccc   q qp  --> q qp
         CFHQ11 = 1.D0 
     x        + 4.d0*ALPSMU * c1_qqp_qqp(q2,mu2,muf2,mur2,mud2,f)
c         write(*,*)"(2)  CFHQ11   = ", CFHQ11
ccc   q qbp --> q qbp
         CFHQ12 = 1.D0 
     x        + 4.d0*ALPSMU * c1_qqbp_qqbp(q2,mu2,muf2,mur2,mud2,f)
ccc   q qb  --> qp qbp
         CFHQ13 = 1.D0 
     x        + 4.d0*ALPSMU * c1_qqb_qpqbp(q2,mu2,muf2,mur2,mud2,f)
ccc   q q   --> q q
         CFHQ14 = 1.D0 
     x        + 4.d0*ALPSMU * c1_qq_qq(q2,mu2,muf2,mur2,mud2,f)
c         write(*,*)"(2)  CFHQ14   = ", CFHQ14
ccc   q qb  --> q qb
         CFHQ15 = 1.D0 
     x        + 4.d0*ALPSMU * c1_qqb_qqb(q2,mu2,muf2,mur2,mud2,f)
ccc   q g   --> q g
         CFHG11 = 1.D0 
     x        + 4.d0*ALPSMU * c1_qg_qg(q2,mu2,muf2,mur2,mud2,f)
ccccccccccccccccccc
ccc   q g   --> g q (OMITTED in practice, hence no FLAV,nf argument)
         CFHG12 = 1.D0 
     x        + 4.d0*ALPSMU * c1_qg_gq(q2,mu2,muf2,mur2,mud2,f)
ccccccccccccccccccc
ccc   g g   --> q qb   
         CFHG14 = 1.D0 
     x        + 4.d0*ALPSMU *c1_gg_qqb(q2,mu2,muf2,mur2,mud2,f)
c$$$         write(*,*)"(2)CFHQ13   = ", CFHQ13
c$$$         write(*,*)"(2)raw walue= ", c1_qqb_qpqbp(q2,mu2,muf2,mur2,mud2,
c$$$     x        flav)
C
       ELSE 
CCC   ELSE LL
          CFHQ1  = 1.D0
          CFHG1  = 1.D0
          CFHQ11 = 1.D0
          CFHQ12 = 1.D0
          CFHQ13 = 1.D0
          CFHQ14 = 1.D0
          CFHQ15 = 1.D0
          CFHQ16 = 1.D0
          CFHG11 = 1.D0
          CFHG12 = 1.D0
          CFHG13 = 1.D0
          CFHG14 = 1.D0         
       ENDIF

c       write(*,*)"CFHQ1 = ", CFHQ1

C... CONSTRUCT THE SUDAKOV EXPONENTS
C AM...Exponents for new diagrams added
       CCQQ     =  CFHQ1 * EXP(2.D0* QIN   + GOUT +FLAG*GAMQQ )
       CCQG     =  CFHG1 * EXP(QIN + GIN   + QOUT +FLAG*GAMQG )
       
c$$$       write(*,*)"GAMQS   = ", flag
c$$$       write(*,*)"GAMGS   = ", flag
c$$$       write(*,*)"CCQQ    = ", ccqq
c$$$       write(*,*)"CCQG    = ", ccqg 
c$$$       write(*,*)"QIN     = ", qin
c$$$       write(*,*)"GIN     = ", gin 
c$$$       write(*,*)"QOUT    = ", qout
c$$$       write(*,*)"GOUT    = ", gout

ccc   partonic proc.'s for fragmentation 
       CCQQ1    = CFHQ11 * EXP(2.D0*QIN + QINFRAG + 1.D0*QOUT)*
     1            ( G11*EXP(GAMQ11) + G12 * EXP(GAMQ12) )    
       CCQQ2    = CFHQ12 * EXP(2.D0*QIN + QINFRAG + 1.D0*QOUT)*
     1            ( G21*EXP(GAMQ21) + G22 * EXP(GAMQ22) )    
       CCQQ3    = CFHQ13 * EXP(2.D0*QIN + QINFRAG + 1.D0*QOUT)*
     1            ( G31*EXP(GAMQ31) + G32 * EXP(GAMQ32) )    
       CCQQ4    = CFHQ14 * EXP(2.D0*QIN + QINFRAG + 1.D0*QOUT)*
     1            ( G41*EXP(GAMQ41) + G42 * EXP(GAMQ42) )    
       CCQQ5    = CFHQ15 * EXP(2.D0*QIN + QINFRAG + 1.D0*QOUT)*
     1            ( G51*EXP(GAMQ51) + G52 * EXP(GAMQ52) )    

cpm  not needed: gluon fragmentation
cpm       CCQQ6= CFHQ16 * EXP(2.D0*QIN+2.D0*QOUT)*
cpm     1            (G61*EXP(GAMQ61)+ G62 * EXP(GAMQ62) )    
       CCQG1    = CFHG11 * EXP(QIN + GIN + QINFRAG + GOUT) *
     1            (G71 * EXP(GAMQ71)+ G72 * EXP(GAMQ72) 
     2            + G73 * EXP(GAMQ73) )    

ccc   CCQG2 not needed (fragmenting gluon, quark --> jet)
       CCQG2    = CFHG12 * EXP(2.D0 * QIN + GIN + GOUT) *
     1            (G81 * EXP(GAMQ81)+ G82 * EXP(GAMQ82) 
     2            + G83 * EXP(GAMQ83))    
cpm  appears to be erronymous: G83-term was missing
c       CCQG2    = CFHG12 * EXP(2.D0 * QIN + GIN + GOUT) *
c     1            (G81 * EXP(GAMQ81)+ G82 * EXP(GAMQ82) )    
cpm  not needed: gluon fragmentation
cpm       CCGG1= CFHG13 * EXP(2.D0*GIN+2.D0*GOUT)*
cpm     1  (G91*EXP(GAMQ91)+ G92 * EXP(GAMQ92)+G93*EXP(GAMQ93))    
       CCGG2    = CFHG14 * EXP(2.D0*GIN + QINFRAG + QOUT) *
     1            (G101 * EXP(GAMQ101)+ G102 * EXP(GAMQ102) )

C... AND THE NP FACTORS
  
       BRE      = BBQ * 2.D0/DSQRT(Q2)
       EXFAC1   = XNP1**2 * BRE**2 / 2.D0
       CCQQ     = CCQQ  * EXP(-EXFAC1)
       CCQG     = CCQG  * EXP(-EXFAC1)
       CCQQ1    = CCQQ1 * EXP(-EXFAC1) 
       CCQQ2    = CCQQ2 * EXP(-EXFAC1)
       CCQQ3    = CCQQ3 * EXP(-EXFAC1)
       CCQQ4    = CCQQ4 * EXP(-EXFAC1)
       CCQQ5    = CCQQ5 * EXP(-EXFAC1)
cpm   Omitting gluon fragmentation parts
cpm       CCQQ6 = CCQQ6*EXP(-EXFAC1)
       CCQG1    = CCQG1 * EXP(-EXFAC1)
cpm   Omitted: FS gluon fragmenting 
       CCQG2    = CCQG2 * EXP(-EXFAC1)
cpm       CCGG1 = CCGG1*EXP(-EXFAC1)
       CCGG2    = CCGG2 * EXP(-EXFAC1)
       RETURN
       END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   various functions and subroutines
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE COMPLEX FUNCTION H0L(L,B,B0,A1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE COMPLEX (A-Z)
      DOUBLE PRECISION PI,A1,B0,IDIRECT, IHADR
C           
      LOGICAL EVOLIS, EVOLFS, PDFINPUT
	INTEGER REG_DIV,EVOLMODE        
      
       COMMON / CONFIG3 / EVOLIS, EVOLFS,EVOLMODE 
      COMMON / CONFIG4 / IDIRECT, IHADR, REG_DIV, PDFINPUT
C     
      PI=DACOS(-1.D0)
C     
      IF(EVOLIS.EQV..TRUE.) THEN
         H0L = A1/2.D0/PI/B0**2 
     X       * ( 2.d0*B + LOG(1.d0-2.d0*B) )
      ELSE
         H0L = A1/2.D0/PI/B0**2 
     X       * ( 2.d0*B + (1.d0-2.d0*L)*LOG(1.d0-2.d0*B) )
      ENDIF
      
      RETURN
      END
C

c mbeekveld modified for fragmentation function evolution
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE COMPLEX FUNCTION H0L_NOEV(L,B,B0,A1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE COMPLEX (A-Z)
      DOUBLE PRECISION PI,A1,B0
      INTEGER EVOLMODE           
      LOGICAL EVOLIS, EVOLFS       
      
       COMMON / CONFIG3 / EVOLIS, EVOLFS,EVOLMODE 
C     
      PI=DACOS(-1.D0)
C     
      IF((EVOLFS.EQV..TRUE.)) THEN
         H0L_NOEV = A1/2.D0/PI/B0**2 
     X       * ( 2.d0*B + LOG(1.d0-2.d0*B) )
      ELSE
         H0L_NOEV = A1/2.D0/PI/B0**2 
     X       * ( 2.d0*B + (1.d0-2.d0*L)*LOG(1.d0-2.d0*B) )
      ENDIF


      RETURN
      END
C
      DOUBLE COMPLEX FUNCTION H0LOLD(L,B,B0,A1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE COMPLEX (A-Z)
      DOUBLE PRECISION PI,A1,B0
      INTEGER EVOLMODE           
      LOGICAL EVOLIS, EVOLFS       
      
       COMMON / CONFIG3 / EVOLIS, EVOLFS,EVOLMODE 
C     
      PI=DACOS(-1.D0)
C     
         H0LOLD = A1/2.D0/PI/B0**2 
     X       * ( 2.d0*B + (1.d0-2.d0*L)*LOG(1.d0-2.d0*B) )
     


      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE COMPLEX FUNCTION H1L(L,B,B0,B1,A1,A2,BB1,
     1     LOMUF,LOMU)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE COMPLEX (A-Z)
      DOUBLE PRECISION PI,A1,A2,B0,B1,BB1,LOMUF,LOMU, 
     X     EVSW,ISLL,ISNLL,ISLNON,ISLNON2,FSLL,FSNLL,FSLNON,
     X     FLAG,IDIRECT, EVOLSW,IHADR
      
      INTEGER ICONF,FLNR,FLNR2,IEV, IEVFS
      INTEGER HS, FS, FS2, REG_DIV, PDFINPUT
      INTEGER QQS, QGS, H2Q, H2G
      INTEGER GAMQS, GAMGS,EVOLMODE 
      
      LOGICAL EVOLIS, EVOLFS     
      
      COMMON / CONFIG / EVOLSW,ICONF,FLNR,FLNR2,IEV, IEVFS
      COMMON / CONFIG2 /ISLL,ISNLL,ISLNON,ISLNON2,FSLL,FSNLL,FSLNON,FLAG
       COMMON / CONFIG3 / EVOLIS, EVOLFS,EVOLMODE 

      COMMON / CONFIG4 / IDIRECT, IHADR, REG_DIV, PDFINPUT
C     
      PI=DACOS(-1.D0)
C     
CCC       LOLAT = LOG(1.-2.*BT) 
C     Effectively setting BT = 0
      LOLAT     = 0d0
      LOLA      = LOG(1.d0-2.d0*B) 
      LOLA_NOEV = LOG(1.d0-2.d0*L)
      FLH1      = 2.d0*B + LOLA
      FRAT      = 1.d0/(1.d0-2.d0*B)
      FRAT_LSV  = (1.d0-2.d0*L)/(1.d0-2.d0*B)
      FLH2      = 2.d0*B*FRAT + LOLA
      FLH2_LSV  = 2.d0*B*FRAT_LSV + LOLA
C     

      IF(EVOLIS.EQV..TRUE.) THEN
         H1L =  A1*B1/2.D0/PI/B0**3 * (1.D0/2.D0*LOLA**2 + FRAT*FLH1) 
     1        + 1.D0/2.D0/PI/B0 * (-A2/PI/B0+A1*LOMU) * FLH2 
     2        + BB1/2.D0/PI/B0 * LOLA
c     3        - A1/2.D0/PI/B0 * LOG(2.D0) * LOLA_NOEV
ccc     2        + BB1/PI/B0 * LOLA
CCC   Should it be BB1/PI/B0 * LOLA ??? Appears not to be the case
      ELSE
CCC   --- ORIGINAL ---
         H1L =  A1*B1/2.D0/PI/B0**3 * (1.D0/2.D0*LOLA**2 
     1        + FRAT_LSV*FLH1) 
     2        + 1.D0/2.D0/PI/B0 * (-A2/PI/B0+A1*LOMU) * FLH2_LSV
     3        - A1/PI/B0 * L * LOMUF 
     4        + BB1/2.D0/PI/B0 * LOLA * 0.d0
c     5        - A1/2.D0/PI/B0 * LOG(2.D0) * LOLA_NOEV
CCC   BB1 term should be omitted (hence the factor of 0.d0)
      ENDIF

c      write(*,*)"evolnll = ", evolnll

CEL Note that the BB1 term is here extra wrt. LSV, and also hep-ph/0409234
CEL because of we added and subtracted a B-term for the initial state 
CEL to make the anomalous dimension complete.

CEL Note also that the renormalization scale is done ok here.

CEL Note that for the WHEPP method, one does not evolve so that
CEL the following term must be added. L(ambda) is not equal to B(eta)!
cpm       IF(EV .EQV. .FALSE.)  THEN
cpm       H1l = H1L - H2S * A1/PI/B0 * L * LOMUF
cpm       ENDIF


c$$$       H1L = A1*B1/2.D0/PI/B0**3 * (1.D0/2.D0*LOLA**2 + FRAT*FLH1) +
c$$$     1       1.D0/2.D0/PI/B0 * (-A2/PI/B0+A1*LOMU) * FLH2 -
c$$$     2       A1/PI/B0 * L * LOMUF + BB1/2.D0/PI/B0 * LOLAT
       RETURN
       END
C


c mbeekveld modified for fragmentation function evolution
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE COMPLEX FUNCTION H1L_NOEV(L,B,B0,B1,A1,A2,BB1,
     1     LOMUF,LOMU)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE COMPLEX (A-Z)
      DOUBLE PRECISION PI,A1,A2,B0,B1,BB1,LOMUF,EVOLSW,LOMU
      LOGICAL EVOLIS, EVOLFS   
      INTEGER ICONF,FLNR,FLNR2,IEV, IEVFS,EVOLMODE   
      
      COMMON / CONFIG / EVOLSW,ICONF,FLNR,FLNR2,IEV, IEVFS
      COMMON / CONFIG2 /ISLL,ISNLL,ISLNON,ISLNON2,FSLL,FSNLL,FSLNON,FLAG
       COMMON / CONFIG3 / EVOLIS, EVOLFS,EVOLMODE 
     
      PI=DACOS(-1.D0)

      LOLAT     = 0d0
      LOLA      = LOG(1.d0-2.d0*B) 
      LOLA_NOEV = LOG(1.d0-2.d0*L)
      FLH1      = 2.d0*B + LOLA
      FRAT      = 1.d0/(1.d0-2.d0*B)
      FRAT_LSV  = (1.d0-2.d0*L)/(1.d0-2.d0*B)
      FLH2      = 2.d0*B*FRAT + LOLA
      FLH2_LSV  = 2.d0*B*FRAT_LSV + LOLA
    
      IF((EVOLFS.EQV..TRUE.)) THEN
         H1L_NOEV =  A1*B1/2.D0/PI/B0**3 * 
     1        (1.D0/2.D0*LOLA**2 + FRAT*FLH1) 
     2        + 1.D0/2.D0/PI/B0 * (-A2/PI/B0+A1*LOMU) * FLH2 
     3        + BB1/2.D0/PI/B0 * LOLA
c mbeekveld was error
c     4        - A1/2.D0/PI/B0 * LOG(2.D0) * LOLA_NOEV
      ELSE
         H1L_NOEV =  A1*B1/2.D0/PI/B0**3 * (1.D0/2.D0*LOLA**2 
     1        + FRAT_LSV*FLH1) 
     2        + 1.D0/2.D0/PI/B0 * (-A2/PI/B0+A1*LOMU) * FLH2_LSV
     3        - A1/PI/B0 * L * LOMUF 
     4        + BB1/2.D0/PI/B0 * LOLA * 0.d0
c mbeekveld was error
c     5        - A1/2.D0/PI/B0 * LOG(2.D0) * LOLA_NOEV
	ENDIF


       RETURN
       END
C
C     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE COMPLEX FUNCTION HPRIME(L,B0,A1,XN)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE COMPLEX (A-Z)
      DOUBLE PRECISION PI, A1, B0, ALPHAS
      
      PI=DACOS(-1.D0)
c mbeekveld fout in L, dat moet namelijk gewoon 1/N zijn!!!
      HPRIME = -A1/2.D0/PI/B0/XN * 
     1     ( LOG(1.D0 - 2.D0*L) ) 
      RETURN
      END

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       DOUBLE COMPLEX FUNCTION H1LOLD(L,B,BT,B0,B1,A1,A2,BB1,
     1                             LOMUF,LOMU)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       IMPLICIT DOUBLE COMPLEX (A-Z)
       DOUBLE PRECISION PI,A1,A2,B0,B1,BB1,LOMUF,LOMU
C
       PI=DACOS(-1.D0)
C
       LOLA  = LOG(1.D0-2.D0*B) 
       LOLAT = LOG(1.D0-2.D0*BT) 
       FLH1  = 2.D0*B + LOLA
       FRAT  = (1.D0-2.D0*L)/(1.D0-2.D0*B)
       FLH2  = 2.D0*B*FRAT + LOLA
C

C
   

       H1LOLD = A1*B1/2.D0/PI/B0**3 * (1.D0/2.D0*LOLA**2 + FRAT*FLH1) +
     1       1.D0/2.D0/PI/B0 * (-A2/PI/B0+A1*LOMU) * FLH2 -
     2       A1/PI/B0 * L * LOMUF + BB1/2.D0/PI/B0 * LOLAT
CP     3        + BB1/2.D0/PI/B0 * LOLAT
CCC   Should it be BB1/PI/B0 * LOLAT ???
       RETURN
       END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE COMPLEX FUNCTION F0L(L,B0,A1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE COMPLEX (A-Z)
      DOUBLE PRECISION PI,A1,B0
C     
      PI   = DACOS(-1.D0)
      LOG1 = LOG(1.D0 - L)
      LOG2 = LOG(1.D0 - 2.D0*L)
C     
c      F0L = 2.d0*H0L(L/2.d0,L/2.d0,B0,A1) 
c     X    - H0L(L,L,B0,A1)
      F0L = (2.d0*H0LOLD(L/2.D0,L/2.D0,B0,A1) 
     X    -       H0LOLD(L,L,B0,A1))
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       DOUBLE COMPLEX FUNCTION F1L(L,B0,B1,A1,A2,BB1,LOMU)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       IMPLICIT DOUBLE COMPLEX (A-Z)
       DOUBLE PRECISION PI,A1,A2,B0,B1,BB1,LOMU,EMC
C
       PI   = DACOS(-1.D0)
       EMC  = 0.5772156649D0
       LOG1 = LOG(1.D0 - L)
       LOG2 = LOG(1.D0 - 2.D0*L)
       F1L = 2.D0*H1LOLD(L/2.D0,L/2.D0,L/2.D0,B0,B1,A1,A2,BB1,0.D0,LOMU)
     1     -      H1LOLD(L,L,L/2.D0,B0,B1,A1,A2,BB1,0.D0,LOMU) 
     2     + A1 * LOG(2.D0)/PI/B0 * ( LOG2 -  LOG1 )      
c     3     + BB1/2.D0/PI/B0 * LOG1
CPM OLD     2     + A1 * LOG(2.D0)/PI/B0 * ( LOG(1.D0-2.D0*L) - LOG(1.D0-L) )      
CEL Note extra B term at the end

c       F1L=2.D0*H1L_NOEV(L/2.D0,L/2.D0,L/2.D0,B0,B1,A1,A2,BB1,0.D0,LOMU)
c     1    -     H1L_NOEV(L,L,L/2.D0,B0,B1,A1,A2,BB1,0.D0,LOMU) 
c     2    + A1 * LOG(2.D0)/PI/B0 * ( LOG2 - LOG1 )      
c     3    - BB1/PI/B0 * LOG1
CEL Note extra B term at the end
       RETURN
       END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       DOUBLE COMPLEX FUNCTION F2L(L,B0,A1,XN)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       IMPLICIT DOUBLE COMPLEX (A-Z)
       DOUBLE PRECISION PI,A1,B0, ALPSMU
C
       PI   = DACOS(-1.D0)
       LOG1 = LOG(1.D0 - L)
       LOG2 = LOG(1.D0 - 2.D0*L)
C
       F2L = A1/2.D0/PI/B0/XN*(LOG2 - LOG1)

       RETURN
       END

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       DOUBLE PRECISION FUNCTION B0CALC(Q2,Q2THR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       IMPLICIT DOUBLE PRECISION (A-Z)
       INTEGER K
       DIMENSION Q2THR(3)
C
       NF = 3
       DO 10 K = 1, 3
       IF (Q2 .GT. Q2THR (K)) THEN
          NF = NF + 1
       ELSE
          GO TO 20
       END IF
  10   CONTINUE
*
  20   B0CALC = (33.- 2.* NF)/12./DACOS(-1.D0)
       RETURN
       END        

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C
C...Q ** 2 - EVOLUTION OF PARTON DENSITIES AND STRUCTURE FUNCTIONS IN
C...NEXT TO LEADING ORDER :
C
       SUBROUTINE RENO (VAN, NS3M, NS8M, NS3N, NS8N, NS15N, NS15M,
     1       NS35N, NS35M, NS24N, NS24M, SIN, GLN, XN)
c       SUBROUTINE RENO (UVN, DVN, SVN, USN, DSN, SSN, CVN, CSN, BVN, 
c     x     BSN,GLN, XN)
       IMPLICIT DOUBLE COMPLEX (A - Z)
       DOUBLE PRECISION XL, XL1, S, EQS (3:6), Q2, Q2START, Q20, Q21,
     1                  Q2MUR, Q2MUF, ALP, ALPS, ALP0, ALP1, ALPC, 
     2                  ALPB, ALPT, ALPQ, ALPQR, LMQ, LMQR,
     3                  Q2THR(3), LAMBDA(3:6), Q2S, Q20F, Q21F,
     4                 IDIRECT, IHADR
       INTEGER F, FI, IRES, ISWI, REG_DIV, PDFINPUT
       COMMON / COUPL / ALPS, ALP0, ALP1, ALPC, ALPB, ALPT, ALPQ, 
     1                  ALPQR, LMQ, LMQR
c       COMMON / SCALES / Q2, Q2START, Q2S, Q20, Q21, Q2MUR, Q2MUF, 
c     1                   Q2THR, Q20F, Q21F, LAMBDA
      COMMON / SCALES / Q2, Q2START, Q20, Q21, Q2MUR, 
     1                  Q2MUF, Q2THR, LAMBDA, Q2S, Q20F, Q21F
       COMMON / RESUM / IRES       

      COMMON / CONFIG4 / IDIRECT, IHADR, REG_DIV, PDFINPUT
       DATA EQS / 6., 10., 11., 15./
       CALL ANCALC (QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     1              QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, C2QI, C2GF, 
     2              CDYQI, CDYGI, XN)
       
C    
c$$$       write(*,*)"Q2      = ", Q2
c$$$       write(*,*)"Q2start = ", Q2start
c$$$       write(*,*)"Q20     = ", Q20
c$$$       write(*,*)"Q21     = ", Q21
c$$$       write(*,*)"Q2mur   = ", Q2mur
c$$$       write(*,*)"Q2muf   = ", Q2muf
c$$$       write(*,*)"Q2thr   = ", Q2thr
c$$$       write(*,*)"lambda  = ", lambda


       ISWI = 0
       GO TO 20
C...PARAMETERS FOR NLO EVOLUTION :
  10   CONTINUE
       CALL ANOM (ANS, AM, AP, AL, BE, AB, RMIN, RPLUS, RQQ, RQG,
     1            RGQ, RGG, C2Q, C2G, CDYQ, CDYG, XN, F,
     2            QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     3            QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, C2QI, C2GF,
     4            CDYQI, CDYGI)
  11   CONTINUE
       S   = LOG (XL)
       XL1 = 1.- XL
C
       ENS = EXP (-ANS*S)
C
       EM  = EXP (-AM*S)
       EP  = EXP (-AP*S)
       AC  = 1.- AL
       EMP = EXP (S * (-AM+AP))
       EPM = EXP (S * (+AM-AP))
       NMP = 1.- AM + AP
       NPM = 1.- AP + AM
       DMQQ =  AL * RQQ + BE * RGQ
       DMQG =  AL * RQG + BE * RGG
       DMGQ =  AB * RQQ + AC * RGQ
       DMGG =  AB * RQG + AC * RGG
       DPQQ =  AC * RQQ - BE * RGQ
       DPQG =  AC * RQG - BE * RGG
       DPGQ = -AB * RQQ + AL * RGQ
       DPGG = -AB * RQG + AL * RGG
       RMMQQ =   AL * DMQQ + AB * DMQG
       RMMQG =   BE * DMQQ + AC * DMQG
       RMMGQ =   AL * DMGQ + AB * DMGG
       RMMGG =   BE * DMGQ + AC * DMGG
       RMPQQ =  (AC * DMQQ - AB * DMQG) / NMP
       RMPQG = (-BE * DMQQ + AL * DMQG) / NMP
       RMPGQ =  (AC * DMGQ - AB * DMGG) / NMP
       RMPGG = (-BE * DMGQ + AL * DMGG) / NMP
       RPMQQ =  (AL * DPQQ + AB * DPQG) / NPM
       RPMQG =  (BE * DPQQ + AC * DPQG) / NPM
       RPMGQ =  (AL * DPGQ + AB * DPGG) / NPM
       RPMGG =  (BE * DPGQ + AC * DPGG) / NPM
       RPPQQ =   AC * DPQQ - AB * DPQG
       RPPQG =  -BE * DPQQ + AL * DPQG
       RPPGQ =   AC * DPGQ - AB * DPGG
       RPPGG =  -BE * DPGQ + AL * DPGG
C...EVOLUTION OF LIGHT PARTON DESITIES BETWEEN THRESHOLDS :
*       IF(ISWI.EQ.0) THEN
         UVN  = UVN  * ENS * (1.+ ALP * XL1 * RMIN)
         DVN  = DVN  * ENS * (1.+ ALP * XL1 * RMIN)
         VAN = VAN * ENS * (1.+ ALP * XL1 * RMIN)
         NS3N = NS3N * ENS * (1.+ ALP * XL1 * RPLUS)
         NS8N = NS8N * ENS * (1.+ ALP * XL1 * RPLUS)
         NS3M = NS3M * ENS * (1.+ ALP * XL1 * RMIN)
         NS8M = NS8M * ENS * (1.+ ALP * XL1 * RMIN)
         NS15N = NS15N * ENS * (1.+ ALP * XL1 * RPLUS)
         NS24N = NS24N * ENS * (1.+ ALP * XL1 * RPLUS)
         NS35N = NS35N * ENS * (1.+ ALP * XL1 * RPLUS)
         NS15M = NS15M * ENS * (1.+ ALP * XL1 * RMIN)
         NS24M = NS24M * ENS * (1.+ ALP * XL1 * RMIN)
         NS35M = NS35M * ENS * (1.+ ALP * XL1 * RMIN)
c	write(*,*) "ENS(1)", ENS*(1.+ ALP * XL1 * RMIN)
c	write(*,*) "ENS(2)", ENS*(1.+ ALP * XL1 * RPLUS)
*        ENDIF
C...EVOLUTION OF LIGHT PARTON DESITIES BETWEEN THRESHOLDS :
C ... TEST: NEGLECT NLO PART OF EVOLUTION (NS ONLY):
*       IF(ISWI.EQ.1) THEN
c         UVN  = UVN  * ENS 
c         DVN  = DVN  * ENS 
c         NS3N = NS3N * ENS 
c         NS8N = NS8N * ENS 
*       ENDIF
C
       SG = SIN
       GL = GLN
c	write(*,*) "EFF", EM * 
c     1   (AL + ALP * (RMMQQ * XL1 + RMPQQ * (EPM-XL))) + 
c     2   EP *(AC + ALP * (RPPQQ * XL1 + RPMQQ * (EMP-XL))) 
c	write(*,*) "EFG", EM * 
c     1    (BE + ALP * (RMMQG * XL1 + RMPQG * (EPM-XL)))+
c     2    EP *(-BE + ALP * (RPPQG * XL1 + RPMQG * (EMP-XL)))
c	write(*,*) "EGF", EM * 
c     1  (AB + ALP * (RMMGQ * XL1 + RMPGQ * (EPM-XL))) +
c     2   EP *(-AB + ALP * (RPPGQ * XL1 + RPMGQ * (EMP-XL)))
c	write(*,*) "EGG", EM *
c     1   (AC + ALP * (RMMGG * XL1 + RMPGG * (EPM-XL))) +
c     2   EP *(AL + ALP * (RPPGG * XL1 + RPMGG * (EMP-XL)))
       SIN = EM * ((AL + ALP * (RMMQQ * XL1 + RMPQQ * (EPM-XL)))* SG
     1           + (BE + ALP * (RMMQG * XL1 + RMPQG * (EPM-XL))) * GL)
     2     + EP * ((AC + ALP * (RPPQQ * XL1 + RPMQQ * (EMP-XL))) * SG
     3           +(-BE + ALP * (RPPQG * XL1 + RPMQG * (EMP-XL))) * GL)
       GLN = EM * ((AB + ALP * (RMMGQ * XL1 + RMPGQ * (EPM-XL))) * SG
     1           + (AC + ALP * (RMMGG * XL1 + RMPGG * (EPM-XL))) * GL)
     2     + EP *((-AB + ALP * (RPPGQ * XL1 + RPMGQ * (EMP-XL))) * SG
     3           + (AL + ALP * (RPPGG * XL1 + RPMGG * (EMP-XL))) * GL)
       IF ( F .EQ. 3 ) THEN
          GO TO 30
       ELSE IF ( (F .EQ. 4).AND.(FI .EQ. 0) ) THEN
          GO TO 40
       ELSE IF ( (F .EQ. 4).AND.(FI .EQ. 1) ) THEN
          GO TO 45
       ELSE IF ( (F .EQ. 5).AND.(FI .EQ. 0) ) THEN
          GO TO 50
       ELSE IF ( (F .EQ. 5).AND.(FI .EQ. 1) ) THEN
          GO TO 55  
       ELSE
          GO TO 60
       END IF
C...INPUT MOMENTS OF THE PARTON DENSITIES (FOR N(F) = 3) :
  20   CONTINUE
       IF (PDFINPUT == 2) THEN
		CALL MMHT2014INP (UV,DV,LS,DEL,SB,SM,GL,XN)
       ELSE IF (PDFINPUT == 1) THEN
		CALL GRV98INP (UV,DV,LS,DEL,SB,SM,GL,XN)
c       ELSE IF (PDFINPUT == 3) THEN
c		CALL MSTW2008INP (UV,DV,LS,DEL,SB,SM,GL,XN)
	ENDIF
c	write(*,*) "HERE"
c	WRITE(*,*) "xn", XN
c	write(*,*) "UV", UV
c	write(*,*) "DV", DV
c	write(*,*) "LS", LS
c	write(*,*) "DEL", DEL
c	write(*,*) "SB", SB
c	write(*,*) "SM", SM
c	write(*,*) "GL", GL
c	write(*,*) "HERE"
        UVN = UV
        DVN = DV
        NS3N = UVN-DVN-2.*DEL
        NS8N = UVN+DVN+LS-2.*SB
        NS3M = UVN-DVN
        NS8M = UVN+DVN-2.*SM
        VAN = UVN + DVN + SM
        SIN = UVN+DVN+LS+SB	
        GLN = GL
	NS15N = 0.D0
	NS15M = 0.D0
	NS24N = 0.D0
	NS24M = 0.D0
	NS35N = 0.D0
	NS35M = 0.D0

*...EVOLUTION BELOW THE CHARM THRESHOLD :
       F=3
       IF ( ALPQ .GE. ALPC ) THEN
          ALP = ALPQ
          ISWI = 1
       ELSE
          ALP = ALPC
       END IF
       XL = ALPS / ALP
       GO TO 10
  30   CONTINUE
       NS15N = SIN
       NS15M = VAN
       IF ( ALPQ .GE. ALPC ) GO TO 70
C...EVOLUTION BETWEEN CHARM THRESHOLD AND INTERMEDIATE SCALE Q0:
       F = 4
       FI = 0
       IF ( ALPQ .GE. ALP0 ) THEN
          ALP = ALPQ
          ISWI = 1
       ELSE
          ALP = ALP0
       END IF
       XL = ALPC / ALP
       GO TO 10
  40   CONTINUE
C
C...EVOLUTION BETWEEN INTERMEDIATE SCALE Q0 AND BOTTOM THRESHOLD :
       FI = 1
       IF ( ALPQ .GE. ALPB ) THEN
          ALP = ALPQ
          ISWI = 1
       ELSE
          ALP = ALPB
       END IF
       XL = ALP0 / ALP
       GO TO 10
  45   CONTINUE
       NS24N = SIN
       NS24M = VAN
       IF ( ALPQ .GE. ALPB ) GO TO 70
C...EVOLUTION BETWEEN BOTTOM THRESHOLD AND INTERMEDIATE SCALE Q1:
       F = 5
       FI = 0
       IF ( ALPQ .GE. ALP1 ) THEN
          ALP = ALPQ
          ISWI = 1
       ELSE
          ALP = ALP1
       END IF
       XL = ALPB / ALP
       GO TO 10
  50   CONTINUE
C
C...EVOLUTION BETWEEN INTERMEDIATE SCALE Q1 AND TOP THRESHOLD :
       FI = 1
       IF ( ALPQ .GE. ALPT ) THEN
          ALP = ALPQ
          ISWI = 1
       ELSE
          ALP = ALPT
       END IF
       XL = ALP1 / ALP
       GO TO 10
  55   CONTINUE
       NS35M = VAN
       NS35N = SIN
       IF ( ALPQ .GE. ALPT ) GO TO 70
C...EVOLUTION ABOVE THE TOP THRESHOLD :
       F = 6
       ALP = ALPQ
       ISWI = 1       
       XL  = ALPT / ALP
       GO TO 10
  60   CONTINUE
C...FLAVOR DECOMPOSITION OF THE QUARK SEA :
  70   CONTINUE
C

c        write(*,*) "NS3M", NS3M
c        write(*,*) "NS8M", NS8M
c        write(*,*) "VAN", VAN
c        write(*,*) "NS3N", NS3N
c        write(*,*) "NS8N", NS8N
c        write(*,*) "SIN", SIN
c        write(*,*) "GLN", GLN
c        write(*,*) "NS15N", NS15M
c        write(*,*) "P15N", NS15N
c        write(*,*) "M24N", NS24M
c        write(*,*) "P24N", NS24N
c        write(*,*) "M35N", NS35M
c        write(*,*) "P35N", NS35N
c        THRD = DCMPLX(1.D0/3.D0,0.D0)
c        TVN = (VAN - NS35M) * 0.5D0 * THRD
c        BVN = TVN + (NS35M - NS24M) * 0.2D0  
c        CVN = BVN + (NS24M - NS15M) * 0.25D0 
c        SVN = CVN + (NS15M - NS8M) * THRD 
c        DVN = SVN + (NS8M - NS3M) * 0.5D0 
c        UVN = SVN + (NS8M + NS3M) * 0.5D0 
c        TP = (SIN - NS35N) * 0.5D0 * THRD
c        BP = TP + (NS35N - NS24N) * 0.2D0
c        CP = BP + (NS24N - NS15N) * 0.25D0
c        SP = CP + (NS15N - NS8N) * THRD
c        DP = SP + (NS8N - NS3N) * 0.5D0
c        UP = SP + (NS8N + NS3N) * 0.5D0
c        TSN = 0.5D0*(TP-TVN)
c        BSN = 0.5D0*(BP-BVN)
c        CSN = 0.5D0*(CP-CVN)
c        SSN = 0.5D0*(SP-SVN)
c        DSN = 0.5D0*(DP-DVN)
c        USN = 0.5D0*(UP-UVN)

c	read(*,*)
	 RETURN
       END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C...PARAMETERS FOR EVOLUTION OF PARTON DISTRIBUTIONS :
C
       BLOCK DATA
       IMPLICIT DOUBLE PRECISION (A-Z)
       DIMENSION LAMBDA (3:6), Q2THR (3), PVA(13), PSG(11)
       INTEGER N
       COMMON / COUPL  / ALPS, ALP0, ALP1, ALPC, ALPB, ALPT, ALPQ, 
     1                   ALPQR, LMQ, LMQR
c       COMMON / SCALES / Q2,Q2START, Q2S, Q20, Q21, Q2MUR, Q2MUF, 
c     1                   Q2THR, Q20F, Q21F, LAMBDA
      COMMON / SCALES / Q2, Q2START, Q20, Q21, Q2MUR, 
     1                  Q2MUF, Q2THR, LAMBDA, Q2S, Q20F, Q21F 
       COMMON / SCHEME / N
       COMMON / INPAR  / PVA, PSG
C...CHOICE OF THE FACTORIZATION SCHEME : N = 1: B, N = 2: A, N = 3: DIS
       DATA N / 1 /
*...HEAVY QUARK THRESHOLDS AND LAMBDA VALUES :
       DATA Q2START / 0.4D0 /
cpm       DATA Q2S /0.3D0/ 
*...Q0 HAS TO BE BETWEEN CHARM AND BOTTOM THRESHOLDS !
*...Q1 HAS TO BE BETWEEN BOTTOM AND TOP THRESHOLDS !
       DATA Q20 / 1.96D0 /
       DATA Q21 / 20.25D0 /
cpm       DATA Q21 / 30.D0 /
cpm       DATA Q20F / 4.D0/                                              
       DATA Q20F / 1.96D0/                                              
       DATA Q21F / 20.25D0 /
cpm       DATA LAMBDA / 0.248, 0.200, 0.131, 0.053 / 
cpm       DATA Q2THR / 2.25, 20.25, 1.D8 / 

cpm       DATA Q2THR   /  1.960,  20.25,  30625. /
cpm       DATA LAMBDA / 0.2994, 0.2460, 0.1677, 0.0678 /

       DATA Q2THR   /  1.960D0,  20.25D0,  30625.D0 /
       DATA LAMBDA / 0.2994D0, 0.2460D0, 0.1677D0, 0.0678D0 /


       DATA Q2S /0.4D0/ 


*
* ..PARTON DISTRIBUTION INPUT PARAMETERS :
       DATA PVA/          0.43D0,  3.09D0,  0.0D0, 18.2D0, 
     1                    0.43D0,  4.09D0,  0.0D0, 18.2D0,
     1           0.20D0,  0.43D0,  12.4D0, -13.3D0, 60.0D0 /
       DATA PSG/ 2.48D0,  0.2D0,   8.5D0,  -2.3D0,  5.7D0,   
     1           0.0D0, 
     2           1.0D0,   1.6D0,   4.1D0,   0.0D0,  0.0D0/
*
       END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C
C...MOMENTS OF HO-INPUT DISTRIBUTIONS (MRS A) AT  Q**2 = 4 GEV**2
C
       SUBROUTINE GRV98INP (UV, DV, LS, DEL, SS, SM, GL, XN)
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER IU, ID, IUPD, ISMOM
       DOUBLE PRECISION PVA(13), PSG(11)
       COMMON / INPAR  / PVA, PSG
C
       IU=2
       ID=1
       IUPD=0
       ISMOM=0
C
* ..NON-SINGLET INPUT PARAMETERS AT Q2S :
       ETAU1 = PVA(1)
       ETAU2 = PVA(2)
       GAMU1 = PVA(3)
       GAMU2 = PVA(4)
       ETAD1 = PVA(5)
       ETAD2 = PVA(6)
       GAMD1 = PVA(7)
       GAMD2 = PVA(8)
       ADEL  = PVA(9)
       ETAE1 = PVA(10)
       ETAE2 = PVA(11)
       GAME1 = PVA(12)
       GAME2 = PVA(13)
*
* ---------------------------------------------------------------------
*
* ..VALENCE QUARK NORMALIZATION :
*    (IU AND ID ARE THE NUMBER OF U AND D VALENCE QUARKS, RESPECTIVELY)
       IF (IUPD .EQ. 0) THEN
* ..(SEPARATE INPUT FOR UV AND DV) 
         AUD = IU/ (CBETA (ETAU1, ETAU2+1.)
     1          * (1.+ GAMU2 * (ETAU1) / (ETAU1+ETAU2+1.))
     2          + GAMU1 * CBETA (ETAU1+0.5, ETAU2+1.))
       ELSE 
* ..(INPUT FOR UV+DV AND DV)   
         AUD = (IU+ID)/ (CBETA (ETAU1, ETAU2+1.)
     1          * (1.+ GAMU2 * (ETAU1) / (ETAU1+ETAU2+1.))
     2          + GAMU1 * CBETA (ETAU1+0.5, ETAU2+1.))
       END IF
       AD  = ID/ (CBETA (ETAD1, ETAD2+1.)
     1        * (1.+ GAMD2 * (ETAD1) / (ETAD1+ETAD2+1.))
     2        + GAMD1 * CBETA (ETAD1+0.5, ETAD2+1.))
*
  21   CONTINUE
*
* ---------------------------------------------------------------------
*
* ..SECOND MOMENT OF THE VALENCE DISTRIBUTION (FOR MOMENTUM SUM) :
       XXN  = DCMPLX (2.D0, 0.D0)
       AV = AUD * (CBETA (XXN+ETAU1-1., ETAU2+1.)
     1        * (1.+ GAMU2 * (ETAU1+XXN-1.) / (ETAU1+ETAU2+XXN))
     2        + GAMU1 * CBETA (XXN+ETAU1-0.5, ETAU2+1.))
       IF (IUPD .EQ. 0) THEN
       AV = AV + AD * (CBETA (XXN+ETAD1-1., ETAD2+1.)
     1        * (1.+ GAMD2 * (ETAD1+XXN-1.) / (ETAD1+ETAD2+XXN))
     2        + GAMD1 * CBETA (XXN+ETAD1-0.5, ETAD2+1.))
       END IF
*
* --------------------------------------------------------------------- 
*
* ..SEA QUARK AND GLUON SHAPE PARAMETERS AT Q2S :
       ALS = DCMPLX(PSG(2), 0.0D0)    
       DS  = DCMPLX(PSG(3), 0.0D0)  
       ETS1= DCMPLX(PSG(4), 0.0D0)
       ETS2= DCMPLX(PSG(5), 0.0D0)
       AST = DCMPLX(PSG(6), 0.0D0)
*
       ATOT = DCMPLX(PSG(7), 0.0D0)
       ALG  = DCMPLX(PSG(8), 0.0D0)
       DG   = DCMPLX(PSG(9), 0.0D0)    
       ETG1 = DCMPLX(PSG(10), 0.0D0)  
       ETG2 = DCMPLX(PSG(11), 0.0D0)

	ASP = DCMPLX(1.D0,0.D0)
	DSP = DCMPLX(1.D0,0.D0)
	ETSP = DCMPLX(1.D0,0.D0)
*
* ---------------------------------------------------------------------
*
* ..SEA QUARK NORMALIZATION :
       IF (ISMOM .EQ. 0) THEN
* ..(SEA NORMALIZATION FACTOR AS INPUT, AS DESCRIBED ABOVE)
         NS  = DCMPLX(PSG(1), 0.0D0)
       ELSE
* ..(SEA MOMENTUM FRACTION AS INPUT)
         AS  = PSG(1)
       END IF
       X0S = ALS + 1.   
       X1S = DS + 1.  
       IF (ISMOM .EQ. 0) THEN
         AS = NS * (CBETA (X0S, X1S) * (1.+ ETS2 * X0S / (X0S+X1S)) 
     1              + ETS1 *  CBETA (X0S+0.5, X1S)) 
       ELSE
         NS = AS / (CBETA (X0S, X1S) * (1.+ ETS2 * X0S / (X0S+X1S)) 
     1              + ETS1 *  CBETA (X0S+0.5, X1S)) 
       END IF
*
* ..GLUON NORMALIZATION :      
       AG  = ATOT - AV - AS 
       X0G = 1.+ ALG     
       X1G = DG + 1.   
       NG = AG / (CBETA (X0G, X1G) * (1.+ ETG2 * X0G / (X0G+X1G)) 
     1           + ETG1 *  CBETA (X0G+0.5, X1G))
* 
* ---------------------------------------------------------------------
*
* ..VALENCE QUARK DISTRIBUTIONS :
       UDV    = AUD * (CBETA (XN+ETAU1-1., ETAU2+1.)
     1          * (1.+ GAMU2 * (ETAU1+XN-1.) / (ETAU1+ETAU2+XN))
     2          + GAMU1 * CBETA (XN+ETAU1-0.5, ETAU2+1.))
       DV     = AD * (CBETA (XN+ETAD1-1., ETAD2+1.)
     1          * (1.+ GAMD2 * (ETAD1+XN-1.) / (ETAD1+ETAD2+XN))
     2          + GAMD1 * CBETA (XN+ETAD1-0.5, ETAD2+1.))
       IF (IUPD .EQ. 0) THEN
         UV = UDV 
       ELSE 
         UV = UDV - DV
       END IF
       DEL      = ADEL * (CBETA (XN+ETAE1-1., ETAE2+1.)
     1          * (1.+ GAME2 * (ETAE1+XN-1.) / (ETAE1+ETAE2+XN))
     2          + GAME1 * CBETA (XN+ETAE1-0.5, ETAE2+1.))
*
* ..SEA AND GLUON DISTRIBUTIONS : 
       XNS = XN - 1.+ ALS 
       LS = NS * (CBETA (XNS, X1S) * (1.+ ETS2 * XNS / (XNS+X1S)) 
     1               + ETS1 * CBETA (XNS+0.5, X1S)) / (1.+ AST * 0.5)
       SS = LS*1.D0
       SM = 0.D0         
       XNG = XN - 1.+ ALG
       GL = NG * (CBETA (XNG, X1G) * (1.+ ETG2 * XNG / (XNG+X1G)) 
     1               + ETG1 * CBETA (XNG+0.5, X1G))

  1    CONTINUE    
*
       RETURN
       END

ccccccccccccc MMHT input

       SUBROUTINE MMHT2014INP (UV, DV, LS, DEL, SS, SM, GL, XN)
       IMPLICIT DOUBLE COMPLEX (A - Z)
C
	AU = 4.2723268205900551d0
	DELTAU = 0.74687d0
	ETAU = 2.7421d0
	AU1 = 0.26349d0
	AU2 = -0.00256d0
	AU3 = 0.25858d0
	AU4 = 0.05d0
	
	AD = 3.3001500335329315d0
	DELTAD = 0.90012d0
	ETAD = -0.58802d0+ETAU
	AD1 = 1.2898d0
	AD2 = 0.60385d0
	AD3 = 0.33950d0
	AD4 = 0.26150d0
	
	AS = 31.329d0
	DELTAS = -0.13358d0
	ETAS = 11.945d0
	AS1 = -1.602d0
	AS2 = 0.86538d0
	AS3 = -0.29923d0
	AS4 = 0.06022d0

	DELTAINT = 0.09531d0
	ADELTA = 7.1043d0
	ETADELTA = ETAS+2.d0
	DELTADELTA=1.7116d0
	GAMMADELTA=10.659d0
	EPSDELTA = -33.341d0

	AG = 0.88744590534913992d0
	DELTAG = -0.45853d0
	ETAG = 2.8636d0
	AG1 = -0.36317d0
	AG2 = 0.20961d0
	AGP = -1.0187d0
	DELTAGP = -0.42510d0
	ETAGP = 32.614d0

	APL = 4.6779d0
	DELTAPL = DELTAS
	ETAPL = 11.588d0
	APL1 = -1.5910
	APL2 = 0.86501

C mstw	
c	AM = -0.011629d0
c	ETAM = 11.261d0
c	DELTAM = 0.2d0
c	X0 = 0.016050d0	
	AM = -0.01614d0
	ETAM = 7.1599
c mbeekveld was an error in the print of 1412.3989!
 	DELTAM = 0.22208d0
c	DELTAM = -0.26403d0
	X0 = 0.026495d0	


       UV = AU*((AU1+AU2+AU3+AU4+1)*CBETA(DELTAU+XN-1., ETAU+1.)+
     1   (-2.*AU1-8.*AU2-18.*AU3-32.*AU4)*CBETA(DELTAU+XN-0.5,ETAU+1.)+
     2	  (8.*AU2+48.*AU3+160.*AU4)*CBETA(DELTAU+XN,ETAU+1.)+
     3          (-32.*AU3-256.*AU4)*CBETA(DELTAU+XN+0.5,ETAU+1.)+
     4           128.*AU4*CBETA(DELTAU+XN+1.,ETAU+1.))
       DV = AD*((AD1+AD2+AD3+AD4+1)*CBETA(DELTAD+XN-1., ETAD+1.)+
     1   (-2.*AD1-8.*AD2-18.*AD3-32.*AD4)*CBETA(DELTAD+XN-0.5,ETAD+1.)+
     2    (8.*AD2+48.*AD3+160.*AD4)*CBETA(DELTAD+XN,ETAD+1.)+
     3            (-32.*AD3-256.*AD4)*CBETA(DELTAD+XN+0.5,ETAD+1.)+
     4            128.*AD4*CBETA(DELTAD+XN+1.,ETAD+1.))

       SS = AS*((AS1+AS2+AS3+AS4+1)*CBETA(DELTAS+XN-1., ETAS+1.)+
     1   (-2.*AS1-8.*AS2-18.*AS3-32.*AS4)*CBETA(DELTAS+XN-0.5,ETAS+1.)+
     2		 (8.*AS2+48.*AS3+160.*AS4)*CBETA(DELTAS+XN,ETAS+1.)+
     3           (-32.*AS3-256.*AS4)*CBETA(DELTAS+XN+0.5,ETAS+1.)+
     4           128.*AS4*CBETA(DELTAS+XN+1.,ETAS+1.))

       GL = AG*((AG1+AG2+1.)*CBETA(DELTAG+XN-1., ETAG+1.)+
     1     (-2.*AG1-8.*AG2)*CBETA(DELTAG+XN-0.5,ETAG+1.)+
     2		   (8.*AG2)*CBETA(DELTAG+XN,ETAG+1.))+
     3	   AGP*CBETA(DELTAGP+XN-1.,ETAGP+1.)

       SP = APL*((APL1+APL2+AS3+AS4+1)*
     1      CBETA(DELTAPL+XN-1., ETAPL+1.)+
     1   (-2.*APL1-8.*APL2-18.*AS3-32.*AS4)*
     2             CBETA(DELTAPL+XN-0.5,ETAPL+1.)+
     2    (8.*APL2+48.*AS3+160.*AS4)*CBETA(DELTAPL+XN,ETAPL+1.)+
     3      (-32.*AS3-256.*AS4)*CBETA(DELTAPL+XN+0.5,ETAPL+1.)+
     4      128.*AS4*CBETA(DELTAPL+XN+1.,ETAPL+1.))

       DEL      = ADELTA * CBETA (XN+DELTADELTA-1., ETADELTA+1.)
     1   * (1.+ GAMDELTA * (DELTADELTA+XN-1.) 
     1    / (DELTADELTA+ETADELTA+XN)
     2    + EPSDELTA*(DELTADELTA+XN)*(DELTADELTA+XN-1.)/
     3     ((DELTADELTA+ETADELTA+XN)*(DELTADELTA+ETADELTA+XN+1.)))

       SM = AM*(CBETA(DELTAM-1.+XN, ETAM+1.)
     1         -CBETA(DELTAM+XN,ETAM+1.)*1./X0)

       LS = SS - SP
       SS = SP

       RETURN
       END

       SUBROUTINE MSTW2008INP (UV, DV, LS, DEL, SS, SM, GL, XN)
       IMPLICIT DOUBLE COMPLEX (A - Z)
C
	AU = 0.25871d0
	ETAU1 = 0.29065d0
	ETAU2 = 3.2432d0
	EPSU = 4.0603d0
	GAMU = 30.687d0

	AD = 12.288d0
	ETAD1 = 0.96809d0
	ETAD2 = 2.7003d0 + ETAU2
	EPSD = -3.8911d0
	GAMD = 6.0542d0

	AS = 0.31620d0
	DELTAS = -0.21515d0
	ETAS = 9.2726d0
	EPSS = -2.6022d0
	GAMS = 30.785d0


	DELTAINT = 0.087673d0
	ADELTA = 8.1084d0
	ETADELTA = 1.8691d0
	GAMDELTA = 13.609d0
	DELTADELTA = -59.289d0

	AG = 1.0805d0
	DELTAG = -0.42848d0
	ETAG = 3.0225d0
	EPSG = -2.2922d0
	GAMG = 3.4894d0
	AGP = -1.1168d0
	DELTAGP = -0.42776d0
	ETAGP = 32.869d0
	
	AP = 0.047915d0
	ETAP = 9.7466d0
	AM = -0.11629d0
	ETAM = 11.261d0
	DELTAM = 0.2d0
	X0 = 0.016050d0

	UV = AU * (CBETA(ETAU1+XN-1.,ETAU2+1.)+
     1      EPSU*CBETA(ETAU1+XN-0.5,ETAU2+1.)+
     2      GAMU*CBETA(ETAU1+XN, ETAU2+1.))
	DV = AD * (CBETA(ETAD1+XN-1.,ETAD2+1.)+
     1      EPSD*CBETA(ETAD1+XN-0.5,ETAD2+1.)+
     2      GAMD*CBETA(ETAD1+XN, ETAD2+1.))
	SS = AS * (CBETA(DELTAS+XN-1.,ETAS+1.)+
     1      EPSS*CBETA(DELTAS+XN-0.5,ETAS+1.)+
     2      GAMS*CBETA(DELTAS+XN, ETAS+1.))
	DEL = ADELTA*(CBETA(ETADELTA+XN-1., ETAS+3.)+
     1    GAMDELTA*CBETA(ETADELTA+XN, ETAS+3.)+
     2    DELTADELTA*CBETA(ETADELTA+XN+1.,ETAS+3.))
	SP = AP * (CBETA(DELTAS+XN-1.,ETAP+1.)+
     1      EPSS*CBETA(DELTAS+XN-0.5,ETAP+1.)+
     2      GAMS*CBETA(DELTAS+XN, ETAP+1.))
	SM = AM * (CBETA(DELTAM+XN-1.,ETAM+1.)-
     1      1./X0*CBETA(DELTAM+XN,ETAM+1.))
	GL = AG * (CBETA(DELTAG+XN-1.,ETAG+1.)+
     1      EPSG*CBETA(DELTAG+XN-0.5,ETAG+1.)+
     2      GAMG*CBETA(DELTAG+XN, ETAG+1.))+
     3     AGP * (CBETA(DELTAGP+XN-1.,ETAGP+1.))

       LS = SS - SP
       SS = SP
       RETURN
       END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C
C...ANOMALOUS DIMENSIONS FOR LEADING AND NEXT TO LEADING ORDER
C...EVOLUTION OF PARTON DENSITIES AND WILSON COEFFICIENTS FOR
C...NLO STRUCTURE FUNCTIONS :
       SUBROUTINE ANOM (ANS, AM, AP, AL, BE, AB, RMIN, RPLUS, RQQ, RQG,
     1                  RGQ, RGG, C2Q, C2G, CDYQ, CDYG, XN, FR,
     2                  QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     3                  QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, C2QI, C2GF,
     4                  CDYQI, CDYGI)
       IMPLICIT DOUBLE COMPLEX (A - Z)
       DOUBLE PRECISION B0, B1
       INTEGER N, F, FR
       COMMON /SCHEME/ N
C...ANOMALOUS DIMENSIONS AND RELATED QUANTITIES IN LEADING ORDER :
c       F = 3
c mbeekveld F=FR should be variable!!!! NOT FIXED
        F = FR
CPM       F=5
       B0 = 11.- 2./3.* FR
CPM       write(51,*) "B0 = ", B0
       B02 = 2.* B0
       QQ = QQI
       QG = F * QGF
       GQ = GQI
       GG = GGI + F * GGF
       SQ = SQRT ((GG - QQ) * (GG - QQ) + 4.* QG * GQ)
       GP = 0.5 * (QQ + GG + SQ)
       GM = 0.5 * (QQ + GG - SQ)
       ANS = QQ / B02
       AM = GM / B02
       AP = GP / B02
       AL = (QQ - GP) / (GM - GP)
       BE = QG / (GM - GP)
       AB = GQ / (GM - GP)

C...NEXT TO LEADING ORDER : ANOMALOUS DIMENSIONS AND WILSON COEFFICIENTS
C...IN THE MS-BAR FACTORIZATION SCHEME OF BARDEEN ET AL. (1981) :
       NS1M = NS1MI + F * NS1F
       NS1P = NS1PI + F * NS1F
       QQ1 = NS1P + F * QQ1F
       QG1 = F * QG1F
       GQ1 = GQ1I + F * GQ1F
       GG1 = GG1I + F * GG1F
       C2Q = C2QI
       C2G = F * C2GF
       CDYQ = CDYQI
       CDYG = CDYGI
C...CHANGE TO THE SCHEME OF ALTARELLI ET AL. (1979) FOR N = 2 OR TO THE
C...DIS SCHEME FOR N = 3 :
       IF (N .EQ. 2) THEN
          DEL = -0.5 * QG
          C2G = C2G + DEL
          QQ1 = QQ1 - DEL * GQ
          QG1 = QG1 + DEL * (QQ - GG - B02 - QG)
          GQ1 = GQ1 + DEL * GQ
          GG1 = GG1 + DEL * (GQ + B02)
       ELSE IF (N .EQ. 3) THEN
          NS1P = NS1P + B02 * C2Q
          NS1M = NS1M + B02 * C2Q
          QQ1 = QQ1 + C2Q * (QG + B02) + C2G * GQ
          QG1 = QG1 + C2Q * QG + C2G * (GG - QQ + QG + B02)
          GQ1 = GQ1 - C2Q * (QQ - GG + GQ + B02) - C2G * GQ
          GG1 = GG1 - C2Q * QG - C2G * (GQ + B02)
          C2Q = 0.
          C2G = 0.
       END IF
       XN1 = XN + 1.
C       C3Q = C2Q - 8./3.* (1./ XN + 1./ XN1)
C       CLQ = 16./ (3.* XN1)
C       CLG = 8.* F / (XN1 * (XN + 2.))
C...COMBINATIONS OF ANOMALOUS DIMENSIONS FOR NLO - SINGLET EVOLUTION :
       B1 = 102 - 38./3.* FR
       B10 = B1 / B0
       RMIN = (NS1M - QQ * B10) / B02
       RPLUS = (NS1P - QQ * B10) / B02
       RQQ = (QQ1 - QQ * B10) / B02
       RQG = (QG1 - QG * B10) / B02
       RGQ = (GQ1 - GQ * B10) / B02
       RGG = (GG1 - GG * B10) / B02
       RETURN
       END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C
C...CALCULATION OF ANOMALOUS DIMENSIONS AND WILSON COEFFICIENTS
C...UP TO THEIR DEPENDENCE OF THE NUMBER OF ACTIVE FLAVORS F :
       SUBROUTINE ANCALC (QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     1                    QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, 
     2                    C2QI, C2GF, CDYQI, CDYGI, XN)
       IMPLICIT DOUBLE COMPLEX (A - Z)
       DOUBLE PRECISION ZETA2, ZETA3, ALPS, ALP0, ALP1, ALPC, ALPB, 
     1                  ALPT, ALPQ, ALPQR, LMQ, LMQR
       COMMON / COUPL  / ALPS, ALP0, ALP1, ALPC, ALPB, ALPT, ALPQ, 
     1                   ALPQR, LMQ, LMQR
C
       XNS = XN * XN
       XN1 = XN + 1.
       XN2 = XN + 2.
       XNM = XN - 1.
C...LEADING ORDER :
       CPSI = PSIFN (XN1) + 0.577216
       QQI = (8./3.) * (-3.- 2./(XN * XN1) + 4.* CPSI)
       QGF = -4.* (XNS + XN +2.) / (XN * XN1 * XN2)
       GQI = -(16./3.) * (XNS + XN + 2.) / (XN * XN1 * XNM)
       GGI = -22.- 24./(XN * XNM) - 24./(XN1 * XN2) + 24.* CPSI
       GGF = 4./3.
C...NEXT OT LEADING ORDER :
       XNT = XNS * XN
       XNFO = XNT * XN
       XN1S = XN1 * XN1
       XN1T = XN1S * XN1
C...ANALYTIC CONTINUATIONS OF N-SUMS AS GIVEN IN GLUECK ET AL. (1990) :
       ZETA2 = 1.644934
       ZETA3 = 1.202057
       CPSI1 = ZETA2 - PSIFN1 (XN1)
       SPMOM = 1.004D0 / XN1 - 0.846D0 / XN2 + 1.342D0 / (XN+3.) -
     1         1.532D0 / (XN+4.) + 0.839D0 / (XN+5.)
       SLC = -5./8.* ZETA3
       SLV = - ZETA2/2.* (PSIFN (XN1/2.) - PSIFN (XN/2.))
     1       + CPSI/XNS + SPMOM
       SSCHLM = SLC - SLV
       SSTR2M = ZETA2 - PSIFN1 (XN1/2.)
       SSTR3M = 0.5 * PSIFN2 (XN1/2.) + ZETA3
       SSCHLP = SLC + SLV
       SSTR2P = ZETA2 - PSIFN1 (XN2/2.)
       SSTR3P = 0.5 * PSIFN2 (XN2/2.) + ZETA3
C...NON-SINGLET PIECES AS GIVEN IN CURCI ET AL. (1980) :
       NS1MA = 16.* CPSI * (2.* XN + 1.) / (XNS * XN1S) +
     1         16.* (2.* CPSI - 1./(XN * XN1)) * ( CPSI1 - SSTR2M ) +
     2         64.* SSCHLM + 24.* CPSI1 - 3. - 8.* SSTR3M -
     3         8.* (3.* XNT + XNS -1.) / (XNT * XN1T) +
     4         16.* (2.* XNS + 2.* XN +1.) / (XNT * XN1T)
       NS1PA = 16.* CPSI * (2.* XN + 1.) / (XNS * XN1S) +
     1         16.* (2.* CPSI - 1./(XN * XN1)) * ( CPSI1 - SSTR2P ) +
     2         64.* SSCHLP + 24.* CPSI1 - 3. - 8.* SSTR3P -
     3         8.* (3.* XNT + XNS -1.) / (XNT * XN1T) -
     4         16.* (2.* XNS + 2.* XN +1.) / (XNT * XN1T)
       NS1B = CPSI * (536./9. + 8.* (2.* XN + 1.) / (XNS * XN1S)) -
     1        (16.* CPSI + 52./3.- 8./(XN * XN1)) * CPSI1 - 43./6. -
     2        (151.* XNFO + 263.* XNT + 97.* XNS + 3.* XN + 9.) *
     3        4./ (9.* XNT * XN1T)
       NS1C = -160./9.* CPSI + 32./3.* CPSI1 + 4./3. +
     1        16.* (11.* XNS + 5.* XN - 3.) / (9.* XNS * XN1S)
       NS1MI = -2./9.* NS1MA + 4.* NS1B
       NS1PI = -2./9.* NS1PA + 4.* NS1B
       NS1F = 2./3. * NS1C
C...SINGLET PIECES AS GIVEN IN FLORATOS ET AL. (1981) :
       XNFI = XNFO * XN
       XNSI = XNFI * XN
       XNSE = XNSI * XN
       XNE = XNSE * XN
       XNN = XNE * XN
       XNMS = XNM * XNM
       XN2S = XN2 * XN2
       XN2T = XN2S * XN2
       QQ1F = (5.* XNFI + 32.* XNFO + 49.* XNT + 38.* XNS + 28.* XN
     1          + 8.) / (XNM * XNT * XN1T * XN2S) * (-32./3.)
       QG1A = (-2.* CPSI * CPSI + 2.* CPSI1 - 2.* SSTR2P)
     1          * (XNS + XN + 2.) / (XN * XN1 * XN2)
     2        + (8.* CPSI * (2.* XN + 3.)) / (XN1S * XN2S)
     3        + 2.* (XNN + 6.* XNE + 15. * XNSE + 25.* XNSI + 36.* XNFI
     4          + 85.* XNFO + 128.* XNT + 104.* XNS + 64.* XN + 16.)
     5          / (XNM * XNT * XN1T * XN2T)
       QG1B = (2.* CPSI * CPSI - 2.* CPSI1 + 5.) * (XNS + XN + 2.)
     1          / (XN * XN1 * XN2)   -   4.* CPSI / XNS
     2        + (11.* XNFO + 26.* XNT + 15.* XNS + 8.* XN + 4.)
     3          / (XNT * XN1T * XN2)
       QG1F = - 12.* QG1A - 16./3.* QG1B
       GQ1A = (-2.* CPSI * CPSI + 10.* CPSI - 2.* CPSI1)
     1          * (XNS + XN + 2.) / (XNM * XN * XN1)  -  4.* CPSI / XN1S
     2        - (12.* XNSI + 30.* XNFI + 43.* XNFO + 28.* XNT - XNS
     3          - 12.* XN - 4.) / (XNM * XNT * XN1T)
       GQ1B = (CPSI * CPSI + CPSI1 - SSTR2P) * (XNS + XN + 2.)
     1          / (XNM * XN * XN1)
     2        - CPSI * (17.* XNFO + 41.* XNS - 22.* XN - 12.)
     3          / (3.* XNMS * XNS * XN1)
     4        + (109.* XNN + 621.* XNE + 1400.* XNSE + 1678.* XNSI
     5          + 695.* XNFI - 1031.* XNFO - 1304.* XNT - 152.* XNS
     6          + 432.* XN + 144.) / (9.* XNMS * XNT * XN1T * XN2S)
       GQ1C = (CPSI - 8./3.) * (XNS + XN + 2.) / (XNM * XN * XN1)
     1        + 1./ XN1S
       GQ1I = - 64./9.* GQ1A - 32.* GQ1B
       GQ1F = - 64./9.* GQ1C
       GG1A = 16./9.* (38.* XNFO + 76.* XNT + 94.* XNS + 56.* XN + 12.)
     1          / (XNM * XNS * XN1S * XN2)   -   160./9.* CPSI + 32./3.
       GG1B = (2.* XNSI + 4.* XNFI + XNFO - 10.* XNT - 5.* XNS - 4.* XN
     1          - 4.) * 16. / (XNM * XNT * XN1T * XN2)   +   8.
       GG1C = (2.* XNFI + 5.* XNFO + 8.* XNT + 7.* XNS - 2.* XN - 2.)
     1          * 64.* CPSI / (XNMS * XNS * XN1S * XN2S)
     2        + 536./9.* CPSI - 64./3.
     3        + 32.* SSTR2P * (XNS + XN + 1.) / (XNM * XN * XN1 * XN2)
     4        - 16.* CPSI * SSTR2P + 32.* SSCHLP - 4.* SSTR3P
     5        - 4.* (457.* XNN + 2742.* XNE + 6040.* XNSE + 6098.* XNSI
     6          + 1567.* XNFI - 2344.* XNFO - 1632.* XNT + 560.* XNS
     7          + 1488.* XN + 576.) / (9.* XNMS * XNT * XN1T * XN2T)
       GG1I = 9.* GG1C
       GG1F = 3./2.* GG1A + 2./3.* GG1B
C...WILSON COEFFICIENTS :
       C2QI = 4./3.* (2.* CPSI * CPSI - 2.* CPSI1 + 3.* CPSI - 9.
     1       - 2.* CPSI / (XN * XN1) + 3./ XN + 4./ XN1 + 2./ XNS)
       C2QI = C2QI - LMQ * QQI/2.      
       C2GF = - 2.* (CPSI * (XNS + XN + 2.) / (XN * XN1 * XN2)
     1       + 1./ XN - 1./ XNS - 6./ XN1 + 6./ XN2)
       C2GF = C2GF - LMQ * QGF/2.    
C... DRELL-YAN COEFFICIENTS :
       CDYQI = 4./3. * (-8. + 8.*ZETA2 + 2./XNS + 2./XN1S - 
     1                   4.*CPSI/XN/XN1 + 4.*CPSI*CPSI +
     2                   (3. + 2./XN/XN1 - 4.*CPSI)*LMQ )
       CDYGI = 1./2. * ( (4. + 14.*XN + 22.*XNS + 11.*XNT + XNFO)/
     1                   (XNS*XN1S*XN2S) - 
     2                    2.*(2. + XN + XNS)*CPSI/(XN*XN1*XN2) +
     3                   (1./XN - 2./XN1 + 2./XN2)*LMQ )       
       RETURN
       END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C
C...PSI - FUNCTION FOR COMPLEX ARGUMENT
       DOUBLE COMPLEX FUNCTION PSIFN (Z)
       DOUBLE COMPLEX Z, ZZ, RZ, DZ, SUB
       SUB = DCMPLX (0.D0,0.D0)
       ZZ = Z
  1    CONTINUE
       IF (DREAL (ZZ) .LT. 10.) THEN
         SUB = SUB - 1./ ZZ
         ZZ = ZZ + 1.
         GOTO 1
       END IF
       RZ = 1./ ZZ
       DZ = RZ * RZ
       PSIFN = SUB + LOG(ZZ) - RZ/2.- DZ/2520. * ( 210.+ DZ * (-21.+
     1         10.*DZ ))
       RETURN
       END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C...FIRST DERIVATIVE OF THE PSI - FUNCTION FOR COMPLEX ARGUMENT :
       DOUBLE COMPLEX FUNCTION PSIFN1 (Z)
       DOUBLE COMPLEX Z, ZZ, RZ, DZ, SUB
       SUB = DCMPLX (0.D0,0.D0)
       ZZ = Z
  1    CONTINUE
       IF (DREAL (ZZ) .LT. 10.) THEN
         SUB = SUB + 1./ (ZZ * ZZ)
         ZZ = ZZ + 1.
         GOTO 1
       END IF
       RZ = 1./ ZZ
       DZ = RZ * RZ
       PSIFN1 = SUB + RZ + DZ/2. * ( 1 + RZ/630. * ( 210.- DZ * ( 42.-
     1         DZ * ( 30.- 42.*DZ ))))
       RETURN
       END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C...SECOND DERIVATIVE OF THE PSI - FUNCTION FOR COMPLEX ARGUMENT :
       DOUBLE COMPLEX FUNCTION PSIFN2 (Z)
       DOUBLE COMPLEX Z, ZZ, RZ, DZ, SUB
       SUB = DCMPLX (0.D0,0.D0)
       ZZ = Z
  1    CONTINUE
       IF (DREAL (ZZ) .LT. 10.) THEN
         SUB = SUB - 2./ (ZZ * ZZ * ZZ)
         ZZ = ZZ + 1.
         GOTO 1
       END IF
       RZ = 1./ ZZ
       DZ = RZ * RZ
       PSIFN2 = SUB - DZ/60. * ( 60.+ RZ * ( 60.+ RZ * ( 30.- DZ *
     1         ( 10.- DZ * ( 10.- DZ * ( 18.- 50.* DZ ))))))
       RETURN
       END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C...BETA FUNCTION FOR COMPLEX ARGUMENT :
       DOUBLE COMPLEX FUNCTION CBETA (Z1, Z2)
       IMPLICIT DOUBLE COMPLEX (A - Z)
       LNGAM (X) = (X - 0.5) * LOG (X) - X + 0.91893853 + 1./(12.* X)
     1              * (1.- 1./(30.* X*X) * (1.- 1./(3.5 * X*X)
     2              * (1.- 4./(3.* X*X))))
       SUB = DCMPLX (0.D0, 0.D0)
       ZZ1 = Z1
  1    CONTINUE
       IF ( DREAL (ZZ1) .LT. 15.) THEN
          SUB = SUB + LOG ((ZZ1+Z2) / ZZ1)
          ZZ1 = ZZ1 + 1.
          GOTO 1
       END IF
       ZZ2 = Z2
  2    CONTINUE
       IF ( DREAL (ZZ2) .LT. 15.) THEN
          SUB = SUB + LOG ((ZZ1+ZZ2) / ZZ2)
          ZZ2 = ZZ2 + 1.
          GOTO 2
       END IF
       LG1 = LNGAM (ZZ1)
       LG2 = LNGAM (ZZ2)
       LG12 = LNGAM (ZZ1 + ZZ2)
       CBETA = EXP (LG1 + LG2 - LG12 + SUB)
       RETURN
       END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C...FUNCTION GAMMA(Z1+Z2)/GAMMA(Z1) :
       DOUBLE COMPLEX FUNCTION CGRAT(Z1, Z2)
       IMPLICIT DOUBLE COMPLEX (A - Z)
       LNGAM (X) = (X - 0.5) * LOG (X) - X + 0.91893853 + 1./(12.* X)
     1              * (1.- 1./(30.* X*X) * (1.- 1./(3.5 * X*X)
     2              * (1.- 4./(3.* X*X))))
       SUB1 = DCMPLX (0.D0, 0.D0)
       SUB2 = SUB1
       ZZ1 = Z1
  1    CONTINUE
       IF ( DREAL (ZZ1) .LT. 15.) THEN
          SUB1 = SUB1 - LOG(ZZ1)
          ZZ1 = ZZ1 + 1.
          GOTO 1
       END IF
       ZZ2 = Z1 + Z2
  2    CONTINUE
       IF ( DREAL (ZZ2) .LT. 15.) THEN
          SUB2 = SUB2 - LOG(ZZ2)
          ZZ2 = ZZ2 + 1.
          GOTO 2
       END IF
       LG1 = LNGAM (ZZ1)
       LG2 = LNGAM (ZZ2)
       CGRAT = EXP ( - LG1 + LG2 - SUB1 + SUB2 )
       RETURN
       END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C...DOUBLE PRECISE GAUSS INTEGRATION :
       FUNCTION DINTEG (F, ALFA, BETA, EPS)
       IMPLICIT DOUBLE PRECISION (A-H, O-Z)
       DIMENSION W(12), X(12)
       DATA CONST / 1.0 D-12 /
       DATA W
     1  /0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,
     2   0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38647,
     3   0.09515 85116 82492, 0.12462 89712 55533, 0.14959 59888 16576,
     4   0.16915 65193 95002, 0.18260 34150 44923, 0.18945 06104 55068/
       DATA X
     1  /0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,
     2   0.18343 46424 95650, 0.98940 09349 91649, 0.94457 50230 73232,
     3   0.86563 12023 87831, 0.75540 44083 55003, 0.61787 62444 02643,
     4   0.45801 67776 57227, 0.28160 35507 79258, 0.09501 25098 37637/
       DINTEG = 0.0 D0
       IF ( ALFA . EQ. BETA ) RETURN
       A = ALFA
       B = BETA
       DELTA = CONST * (DABS(A-B))
       AA = A
    1  Y = B - AA
       IF( DABS(Y) .LE. DELTA ) RETURN
    2  BB = AA + Y
       C1 = 0.5 D0 * (AA + BB)
       C2 = C1 - AA
       S8 = 0.0 D0
       S16 = 0.0 D0
       DO 15 I = 1, 4
          C3 = X(I) * C2
          S8 = S8 + W(I) * (F(C1+C3) + F(C1-C3))
   15  CONTINUE
       DO 16 I = 5, 12
          C3 = X(I) * C2
          S16 = S16 + W(I) * (F(C1+C3) + F(C1-C3))
  16   CONTINUE
       S8 = S8 * C2
       S16= S16 * C2
       IF( DABS(S16-S8) .GT. EPS * DABS(S8)) THEN
          Y = 0.5 * Y
          IF ( DABS(Y) .LE. DELTA ) THEN
             DINTEG = 0.0
             WRITE (*,10)
  10         FORMAT (1X,' DINTEG : TOO HIGH ACCURACY ')
          ELSE
             GOTO 2
          END IF
       ELSE
          DINTEG = DINTEG + S16
          AA = BB
          GOTO 1
       END IF
       RETURN
       END



