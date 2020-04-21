ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       FUNCTION BSTINT(XX,WGT)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
CAM ..NEW VARIABLES PDGNQQ1,CDGNQQ1, DQGAM,SIGDQQ1,CCQQ1, CFHQ1 ETC ADDED
CAM ..DQGAM AND DGGAM ARE ff's not yet defined
       IMPLICIT DOUBLE PRECISION (A-I,L-Z)
       DOUBLE COMPLEX XN, XNM, XNBAR, CC, ZPII, CEX, PDGNQ, PDGNG,
     1                CDGQ, CDGG, SIGDGQ, SIGDGG, BBC(3), BBQ, BRE, 
     2                XJBC(3), ANGFAC, ANGF1, ANGF2, IC,
     3                PDGNGG,PDGNQQ,F1COMBN,C1QQ,C1GQ,C1QG,C1GG,
     4                E1QQ,E1GQ,E1QG,E1GG,BBQ2, COMMU1,
ccc
     x                PDGNQ1,PDGNQ2,
     x                PDGNQ3,PDGNQ4,PDGNQ5,PDGNQ6,
ccc
     5                PDGNQQ1,PDGNQQ2,PDGNQQ3,PDGNQQ4,PDGNQQ5,
     6                PDGNQQ6,PDGNGMOD,PDGNGG1,PDGNGG2,DQGAM,DGGAM,
     7                SIGDQQ1,SIGDQQ2,SIGDQQ3,SIGDQQ4,SIGDQQ5,
     8                SIGDQQ6,SIGDQG,SIGDGG1,SIGDGG2, DFRAG, ZZ,
     9                FRAGCOMP, Q2SF, Q2THRF, LAMBDAF
       double complex CCQQ1, CCQQ2, CCQQ3, CCQQ4,
     x      CCQQ5, CCQQ6, CCQG1, CCQG2, CCGG1, CCGG2
       
       DOUBLE COMPLEX CMU0,CMU1,ENSP,ENSM,EFF,EFG,EGF,EGG,CHI,
     1      UVN,USN,DVN,DSN,SSN,CSN,GLN,T3N,T8N,SINN,SINTMP,UDSN,UDSBN,
     2      CALCQQ,CALCQG,CALCGQ,CALCGG,CCQQ,CCQG,SVN,VAN, 
     2      NS3M, NS8M, NS3N, NS15N, NS15M,TVN, TSN, BSN,
     1      NS35N, NS35M, NS24N, NS24M,TP, BP, SP, CP, UP, DP,
     3      CCQQSUD,CCQGSUD, ENS_LL, ENS_GRES,NS8N,
     1      CCQQ1SUD,CCQQ2SUD,CCQQ3SUD,CCQQ4SUD,CCQQ5SUD,
     2      CCQQ6SUD,CCQG1SUD,CCQG2SUD,CCGG1SUD,CCGG2SUD,
     4      UN, DN, SN, CN, GLNF, BVN, CVN,
     5      sig_direct, sig_frag

       double precision evsw, FQT

       INTEGER IRES, NB, IREC, ILNNON, ICHI,IFLAG,IPROFILE 
	INTEGER ICONF,FLNR,FLNR2,IEV, IEVFS, REG_DIV
       INTEGER HS,FS,FS2, NBINT
       INTEGER QQS, QGS, H2Q, H2G
       INTEGER GAMQS,GAMGS
       INTEGER N, FLAV,EVOLMODE 

       LOGICAL EVOLIS, misnan, crap,EVOLFS, writeit

CPM    Added one extra dimension in XX: FF integration

       DIMENSION Q2THR(3), LAMBDA (3:6), XX(5)

       external misnan

C  
       common /deb/ crap

       COMMON / COUPL  / ALPS, ALP0, ALP1, ALPC, ALPB, ALPT, ALPQ, 
     1      ALPQR, LMQ, LMQR 

c       COMMON / SCALES / Q2, Q2START, Q2S, Q20, Q21, Q2MUR, Q2MUF, 
c     1      Q2THR, Q20F, Q21F, LAMBDA
      COMMON / SCALES / Q2, Q2START, Q20, Q21, Q2MUR, 
     1                  Q2MUF, Q2THR, LAMBDA, Q2S, Q20F, Q21F
       COMMON / CUT / MUBAR, XNP1, XNP2
       COMMON / CONT / C, CO, SI, XT2, PI, BQLIM, VP
       COMMON / FOURIER / BST, QQT
       COMMON / KINE / FAC, PT, SQS
       COMMON / RESUM / IRES
       COMMON / QTINT / IPROFILE
       COMMON / RECOIL / IREC, ICHI
       COMMON / SVCAPP / ILNNON
C       COMMON / TEST / ICONF,FLNR,HS,FS,FS2,EV,EV2,QQS,QGS,GAMQS,GAMGS,
C     1      H2Q, H2G
       COMMON / CONFIG / EVOLSW,ICONF,FLNR,FLNR2,IEV, IEVFS
       COMMON / CONFIG2 / ISLL,ISNLL,ISLNON,FSLL,FSNLL,FSLNON,FLAG
       COMMON / CONFIG3 / EVOLIS, EVOLFS,EVOLMODE 
       COMMON / CONFIG4 / IDIRECT, IHADR, REG_DIV
       COMMON / FRAG / IFLAG
       COMMON / SCHEME / N 
       COMMON / FLAVORS /  FLAV
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   C0de starts here...
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     For testing: 
c$$$       xx(1) = 0.1
c$$$       xx(2) = 0.1
c$$$       xx(3) = 0.1
c$$$       xx(4) = 0.1
c$$$       xx(5) = 0.1
C     For testing: 
c       xx(1) = 0.1d0
c       xx(2) = 0.2d0
c       xx(3) = 0.3d0
c       xx(4) = 0.4d0
c       xx(5) = 0.5d0
C
C ... INTEGRATION OVER QT = K:
       crap = .false.
       writeit = .false.
c mbeekveld 
	IF(REG_DIV.eq.2) THEN
c mbeekveld if deformation then d2 QT integral is gone:
	QT = 0.D0
	XJ5 = 1.D0
	PHI = 0.D0
	XJ4 = 1.D0
    
	ELSE 
       IF (IPROFILE.EQ.0) THEN
          QT  = MUBAR
          XJ5  = 1.D0   
        PHI = 0.D0
        XJ4 = 1.D0
       ELSE IF (IPROFILE.EQ.1) THEN
          QTDO = 0.D0
          QTUP = MUBAR
          QT   = QTDO + (QTUP-QTDO) * XX(5)
          XJ5  = QTUP-QTDO        
       PHIDO = 0.D0
       PHIUP = 2.D0*PI
       PHI = PHIDO + (PHIUP-PHIDO) * XX(4)
       XJ4 = PHIUP-PHIDO
c     FQT = QT / (2.D0*PI)**2
       ENDIF
C
C ... INTEGRATION OVER phi_p:
	ENDIF
C
C ... CALCULATE   4|Pt - Qt/2|**2/S :
C
       
	IF(REG_DIV.eq.2) THEN
c mbeekveld it is multiplied by xt^2, not x~T^2
         XTP2  = 4.D0 * PT**2 / SQS**2
       ELSE IF (IREC.EQ.0) THEN 
         XTP2  = 4.D0 * PT**2 / SQS**2
       ELSE IF (IREC.EQ.1) THEN 
         XTP2  = 4.D0 * ( PT**2 + QT**2/4.D0 
     x        - PT*QT*COS(PHI) ) / SQS**2
       ENDIF 
C
       AX      = LOG (XTP2)
C
C ... INTEGRATION OVER N:
C
c mbeekveld change X(3) into X(1)!
       Z       = XX(1)
       DZ      = Z/(1.D0-Z)
C       
       XN      = DCMPLX(C+CO*DZ, SI*DZ)
       XNM     = (- XN - 2.d0) * AX
       CC      = DCMPLX(CO, SI)
       ZPII    = 2.d0*PI*DCMPLX(0.d0,1.d0)
       CEX     = EXP (XNM) / ZPII * CC
       EMC     = 0.5772156649D0
C
       XN      = XN + 1.d0
       XNBAR   = XN*EXP(EMC)
C
C ... CALL HERE EVERYTHING RELATED TO PARTON DENSITIES, BORN CROSS
C ... SECTIONS, FF'S (PHFNLO), ETC: 
C ... THIS IS ALL INDEPENDENT OF B AND THETA !
C
CCC   MOMENTS OF BORN CROSS SECTIONS
c       write(*,*)"calling DGDIGLO..."


c mbeekveld catch nans:
c false behaviour in sampling
        IF(ISNAN(XX(3))) THEN
        GOTO 77
        ELSE IF(ISNAN(XX(2))) THEN
        GOTO 77
        ELSE IF(XX(2).eq.1) THEN
        GOTO 77
        ELSE IF(ISNAN(XX(1))) THEN
        GOTO 77
        ELSE IF(XX(1).eq.1) THEN
        GOTO 77
        ELSE IF(ISNAN(XX(4))) THEN
        GOTO 77
        ELSE IF(ISNAN(XX(5))) THEN
        GOTO 77
        ENDIF

       CALL DGDIGLO(XN-1.d0,SIGDGQ,SIGDGG,SIGDQQ1,SIGDQQ2,SIGDQQ3,
     1      SIGDQQ4,SIGDQQ5,SIGDQQ6,SIGDQG,SIGDGG1,SIGDGG2)
C
C ... INTEGRATION OVER B :
C ... NEED THREE DIFFERENT B's 
C
	IF (reg_div.eq.2) THEN
c mbeekveld changed: if contour is deformed, b integrals are gone
c and b is equal to -i*(N+1)/PT
c note that this b definition is slightly different * sqrt(Q2)/2!
	XJBC(1) = DCMPLX(1.D0,0.D0)
	BBC(1) = DCMPLX(0.D0,-1.D0)*(XN+1.D0)/(PT * 2.D0/DSQRT(Q2))
	NBINT = 1
	ELSE
       ZB1     = BQLIM * XX(2)
       BBC(1)  = DCMPLX(ZB1,0.D0)
       XJBC(1) = DCMPLX(BQLIM,0.D0)
C
       ZB0     = XX(2)
       ZB2     = ZB0/(1.D0-ZB0)
       BBC(2)  = DCMPLX(BQLIM,0.D0) - ZB2 / CC 
       XJBC(2) = - 1.D0/CC/(1.D0-ZB0)**2 
C
       ZB3     = ZB2
       BBC(3)  = DCMPLX(BQLIM,0.D0) - ZB3 * CC 
       XJBC(3) = - CC/(1.D0-ZB0)**2 
	NBINT = 3
	ENDIF

C
C ... LOOP OVER B VALUES :
C
       YSUM    = 0.D0
C
       DO 2 NB = 1, NBINT
       BBQ     = BBC(NB)       
       BRE     = BBQ * 2.D0/DSQRT(Q2)    
   
C
c mbeekveld vp is contour parameter
c mbeekveld not needed if NBINT = 1
	
	IF (NBINT.eq.3) THEN
C ... INTEGRATION OVER X_THETA :C
       XTH = XX(3)
       IC = DCMPLX(0.D0,1.D0)
       ANGF1 = DCMPLX(1.D0,-2.D0*VP) * 
     1   EXP(-IC*BRE*QT*SIN(PI*DCMPLX(-XTH,VP*(-1.D0+2.D0*XTH)))) 
       ANGF2 = DCMPLX(1.D0,+2.D0*VP) * 
     1   EXP(-IC*BRE*QT*SIN(PI*DCMPLX( XTH,VP*(-1.D0+2.D0*XTH)))) 
C
C ... ASSIGN ANGULAR FACTORS :
       IF(NB.EQ.1) THEN 
         ANGFAC = ANGF1 + ANGF2
       ELSE IF(NB.EQ.2) THEN 
         ANGFAC = ANGF1
       ELSE IF(NB.EQ.3) THEN 
         ANGFAC = ANGF2
       ENDIF 
	ELSE IF(NBINT.eq.1) THEN
	ANGFAC = 1.D0
	ENDIF

C 
C...PROMPT PHOTON RESUMMATION COEFFICIENT
C
       BBQ2 = BBQ * EXP(EMC) 
       IF(ILNNON.EQ.1) THEN
          IF(ICHI.EQ.1) THEN
             CHI = XN + BBQ
          ELSEIF(ICHI.EQ.2) THEN
             CHI = BBQ2 + XNBAR/(1.D0+BBQ2/4.D0/XNBAR)
          ENDIF
       ELSE
          CHI = (1.0,0.0)
       ENDIF
	IF (NBINT.eq.1) THEN
	CHI = XNBAR
	ENDIF
ccc   complex scales
       CMU0 = DCMPLX(SQRT(Q2MUF))
       CMU1 = DCMPLX(SQRT(Q2)/CHI)
cmbeekveld complex scale added for fragmentation functions
       IF(EVOLFS.eqv. .TRUE.) then
       COMMU1 = DCMPLX(SQRT(Q2), 0)/XNBAR
       ELSE	
       COMMU1 = DCMPLX(SQRT(Q2MUF), 0)
       ENDIF
       if(ihadr.eq. 1.d0) then
          CALL PHFNLO(2.D0*(XN-1.D0)+3.D0, Q2MUF, 
     x        COMMU1*COMMU1, UN, DN, SN, CN, GLNF)
       endif

C... HERE WE IMPLEMENT THE KSV METHOD
C... FIRST CALL INVIDUAL PDF'S AND MAKE (NON)SINGLET COMBINATIONS

         CALL RENO (VAN, NS3M, NS8M, NS3N,NS8N,NS15N, NS15M,
     1       NS35N, NS35M, NS24N, NS24M, SINTMP, GLN, XN) 
 
cCALL RENO (UVN, DVN, SVN, USN, DSN, SSN, CVN, CSN, BVN, 
c     x     BSN,GLN, XN)
c       T3N    = UVN+2.D0*USN-(DVN+2.D0*DSN)
c mbeekveld changed t8n definition to hold svn
c       T8N    = UVN+2.D0*USN+DVN+2.D0*DSN-2.D0*(SVN+2.D0*SSN)
c mbeekveld changed sintmp def to hold svn
c       SINTMP = UVN+2.D0*USN+DVN+2.D0*DSN+SVN+2.D0*SSN

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   Evolving the PDF's
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C... NOW EVOLVE DOWN TO SCALE Q/CHI
       IF    (EVOLIS  .EQV. .TRUE.) THEN
	IF (EVOLMODE .eq. 1) THEN
          CALL EVOLMATNLL(XN,CMU0,CMU1,ENSP,
     x      ENSM,EFF,EFG,EGF,EGG,ENS_LL)
	ELSEIF (EVOLMODE.eq. 2) THEN
          CALL EVOLMATNLLON(XN,CMU0,CMU1,ENSP,
     x       ENSM,EFF,EFG,EGF,EGG,ENS_LL)
	ELSEIF (EVOLMODE.eq. 3) THEN
          CALL EVOLMATLNON(XN,CMU0,CMU1,ENSP,ENSM,EFF,
     x    EFG,EGF,EGG,ENS_LL)
	ELSEIF (EVOLMODE.eq. 4) THEN
          CALL EVOLMAT(XN,CMU0,CMU1,ENSP,ENSM,EFF,EFG,EGF,EGG,ENS_LL)
	ENDIF
       ENDIF
CEL NOTE: when using the WHEPP method and commenting out the
CEL evolution below, don't forget to include the LOMU term 
CEL in H1L if NLL is used, as well as switch-on line 5 in 
CEL CFHQ1, CFHG1
       IF    (EVOLIS.EQV..TRUE.) THEN         
         VAN = VAN * ENSM
         NS3N = NS3N * ENSP
         NS8N = NS8N * ENSP
         NS3M = NS3M * ENSM
         NS8M = NS8M * ENSM
         NS15N = NS15N * ENSP
         NS24N = NS24N * ENSP
         NS35N = NS35N * ENSP
         NS15M = NS15M * ENSM
         NS24M = NS24M * ENSM
         NS35M = NS35M * ENSM 
          SINN = EFF*SINTMP+EFG*GLN
          GLN  = EGF*SINTMP+EGG*GLN          
C..  COMBINE BACK INTO DENSITIES
	ELSE 
	   SINN = SINTMP
       ENDIF
	
        THRD = DCMPLX(1.D0/3.D0,0.D0)
        TVN = (VAN - NS35M) * 0.5D0 * THRD
        BVN = TVN + (NS35M - NS24M) * 0.2D0  
        CVN = BVN + (NS24M - NS15M) * 0.25D0 
        SVN = CVN + (NS15M - NS8M) * THRD 
        DVN = SVN + (NS8M - NS3M) * 0.5D0 
        UVN = SVN + (NS8M + NS3M) * 0.5D0 
        TP = (SINN - NS35N) * 0.5D0 * THRD
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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c... CONSTRUCT FLUXES, KEEP CHARGE^2 FACTORS IN MIND!       
ccc   modify fluxes
       CALL PDGMOD(UVN,DVN,SVN,CVN,BVN,USN,DSN,SSN,CSN,BSN,GLN,
     x     UN,DN,SN,GLNF,
     x     PDGNQ,PDGNG,
     x     PDGNQQ1,PDGNQQ2,PDGNQQ3,PDGNQQ4,PDGNQQ5,PDGNQQ6,
     x     PDGNGG1,PDGNGG2,PDGNGMOD,
     x     IFLAG)
     

C
C... NOW COMBINE WITH COEFFICIENT FUNCTIONS
c..   CCQQ1,...CCGG2 defined below
       CCQQ  = 1.D0
       CCQG  = 1.D0
       CCQQ1 = 1.D0
       CCQQ2 = 1.D0
       CCQQ3 = 1.D0
       CCQQ4 = 1.D0
       CCQQ5 = 1.D0
       CCQQ6 = 1.D0
       CCQG1 = 1.D0
       CCQG2 = 1.D0
       CCGG1 = 1.D0
       CCGG2 = 1.D0
C... GET THE SUDAKOV EXPONENTS (INCL. NP PARTS)

       CALL CDGRES2(ALPQR,5,XN,LMQ,LMQR,BBQ,CCQQSUD,CCQGSUD,
     1  CCQQ1SUD,CCQQ2SUD,CCQQ3SUD,CCQQ4SUD,CCQQ5SUD,
     2  CCQQ6SUD,CCQG1SUD,CCQG2SUD,CCGG1SUD,CCGG2SUD,ENS_GRES)

       CCQQ  = CCQQ *CCQQSUD
       CCQG  = CCQG *CCQGSUD
       CCQQ1 = CCQQ1*CCQQ1SUD
       CCQQ2 = CCQQ2*CCQQ2SUD
       CCQQ3 = CCQQ3*CCQQ3SUD
       CCQQ4 = CCQQ4*CCQQ4SUD
       CCQQ5 = CCQQ5*CCQQ5SUD
       CCQG1 = CCQG1*CCQG1SUD
       CCQG2 = CCQG2*CCQG2SUD
       CCGG2 = CCGG2*CCGG2SUD

CCC   Direct photon cross section
       if(idirect.eq.1.d0) then
          sig_direct = PDGNQ    * SIGDGQ   * CCQQ
     x               + PDGNG    * SIGDGG   * CCQG
       else
          sig_direct = 0.d0
       endif
        
CCC   Fragmentation photon cross section
       if(ihadr.eq.1.d0) then
          sig_frag   = PDGNQQ1  * SIGDQQ1  * CCQQ1 
     x               + PDGNQQ2  * SIGDQQ2  * CCQQ2 
     x               + PDGNQQ3  * SIGDQQ3  * CCQQ3 
     x               + PDGNQQ4  * SIGDQQ4  * CCQQ4 
     x               + PDGNQQ5  * SIGDQQ5  * CCQQ5 
     x               + PDGNGMOD * SIGDQG   * CCQG1 
ccc
ccc   omit: x               + PDGNGMOD * SIGDQG   * CCQG2 
ccc
     x               + PDGNGG2  * SIGDGG2  * CCGG2 
       else
          sig_frag   = 0.d0
       endif
c mbeekveld if contour is deformed, BRE is not needed
	IF (reg_div.eq.2) THEN
C mbeekveld check this, I believe 2pi^2 is needed
	BRE = 1.D0
	ENDIF
CCC   Integrand summed over three different domains (b-integral)
       DGN =  2.D0 * DREAL( 
     x        CEX * ANGFAC   * XJBC(NB) * BRE   * ( 
     x        sig_direct 
     x      + sig_frag 
     x      ) 
     x      )
       
c$$$       write(*,*)"CCQQSUD    = ", ccqqsud
c$$$       write(*,*)"CCQGSUD    = ", ccqgsud
c$$$       write(*,*)"PDGNQ      = ", pdgnq
c$$$       write(*,*)"PDGNG      = ", pdgng
c$$$       write(*,*)"SIGDGQ     = ", sigdgq
c$$$       write(*,*)"SIGDGG     = ", sigdgg
c$$$       write(*,*)"CCQQ       = ", ccqq
c$$$       write(*,*)"CCQG       = ", ccqg
c$$$       write(*,*)"--------------------"
c$$$       write(*,*)"DGN        = ", DGN
c$$$       write(*,*)"CEX        = ", CEX
c$$$       write(*,*)"XTP2       = ", XTP2
c$$$       write(*,*)"ANGFAC     = ", ANGFAC
c$$$       write(*,*)"XJBC       = ", XJBC
c$$$       write(*,*)"BRE        = ", BRE
c$$$       write(*,*)"--------------------"

c$$$       DGN = 2.D0 * DREAL( CEX * ANGFAC   * XJBC(NB) * BRE   * ( 
c$$$     1      IDIRECT * (
c$$$     1                           PDGNQ    * SIGDGQ   * CCQQ 
c$$$     2                        +  PDGNG    * SIGDGG   * CCQG  
c$$$     2      )
c$$$     3      + IHADR * (
c$$$     3                           PDGNQQ1  * SIGDQQ1  * CCQQ1 
c$$$     4                         + PDGNQQ2  * SIGDQQ2  * CCQQ2 
c$$$     5                         + PDGNQQ3  * SIGDQQ3  * CCQQ3 
c$$$     6                         + PDGNQQ4  * SIGDQQ4  * CCQQ4 
c$$$     7                         + PDGNQQ5  * SIGDQQ5  * CCQQ5 
c$$$     9                         + PDGNGMOD * SIGDQG   * CCQG1 
c$$$     4                         + PDGNGG2  * SIGDGG2  * CCGG2 
c$$$cpm   we do not include the following:
c$$$cpm   the gg final states.  
c$$$cpm     8                           PDGNQQ6  * SIGDQQ6  * CCQQ6 +
c$$$cpm   Term proportional to CCQG2 not included as this is the case where
c$$$cpm   the gluon fragments and the quark --> jet
c$$$cpm     1                         + PDGNGMOD * SIGDQG   * CCQG2
c$$$cpm   we do not include the gg final states.  
c$$$cpm   contributions   -- here for clarification
c$$$cpm     3                           PDGNGG1*SIGDGG1*CCGG1+
c$$$     4      )
c$$$     4      ) 
c$$$     5      )

c$$$          write(*,*)"CCQQSUD    = ", ccqqsud
c$$$          write(*,*)"CCQGSUD    = ", ccqgsud
c$$$          write(*,*)"PDGNQ      = ", pdgnq
c$$$          write(*,*)"PDGNG      = ", pdgng
c$$$          write(*,*)"SIGDGQ     = ", sigdgq
c$$$          write(*,*)"SIGDGG     = ", sigdgg
c$$$          write(*,*)"CCQQ       = ", ccqq
c$$$          write(*,*)"CCQG       = ", ccqg
c$$$
c$$$          write(*,*)"SIGDGQ     = ", sigdgq
c$$$          write(*,*)"SIGDGG     = ", sigdgg
c$$$          write(*,*)"CCQQ       = ", ccqq
c$$$          write(*,*)"CCQG       = ", ccqg
c$$$
c$$$          write(*,*)"CCQQ1      = ", ccqq1
c$$$          write(*,*)"CCQQ2      = ", ccqq2
c$$$          write(*,*)"CCQQ3      = ", ccqq3
c$$$          write(*,*)"CCQQ4      = ", ccqq4
c$$$          write(*,*)"CCQQ5      = ", ccqq5
c$$$
c$$$          write(*,*)"PDGNQQ1    = ", pdgnqq1
c$$$          write(*,*)"PDGNQQ2    = ", pdgnqq2
c$$$          write(*,*)"PDGNQQ3    = ", pdgnqq3
c$$$          write(*,*)"PDGNQQ4    = ", pdgnqq4
c$$$          write(*,*)"PDGNQQ5    = ", pdgnqq5
c$$$
c$$$          write(*,*)"SIGDQQ1    = ", sigdqq1
c$$$          write(*,*)"SIGDQQ2    = ", sigdqq2
c$$$          write(*,*)"SIGDQQ3    = ", sigdqq3
c$$$          write(*,*)"SIGDQQ4    = ", sigdqq4
c$$$          write(*,*)"SIGDQQ5    = ", sigdqq5
c$$$
c$$$          write(*,*)"DGN        = ", dgn

c$$$          write(*,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
c$$$          stop
       if(misnan(DGN) ) then
          if(writeit.eqv..TRUE.) then
          write(*,*)"WGT        = ", wgt
          write(*,*)"QT         = ", qt
          write(*,*)"PHI        = ", phi
          write(*,*)"CEX        = ", cex
          write(*,*)"AX         = ", AX
          write(*,*)"XTP2       = ", XTP2
          write(*,*)"XNM        = ", XNM
          write(*,*)"XN         = ", XN
          write(*,*)"CO         = ", CO 
          write(*,*)"SI         = ", SI
          write(*,*)"DZ         = ", DZ
          write(*,*)"Z          = ", Z
          write(*,*)"ZPII       = ", ZPII
          write(*,*)"C          = ", C
          write(*,*)"CC         = ", CC
          

          write(*,*)"XX(1)      = ", xx(1)
          write(*,*)"XX(2)      = ", xx(2)
          write(*,*)"XX(3)      = ", xx(3)
          write(*,*)"XX(4)      = ", xx(4)
          write(*,*)"XX(5)      = ", xx(5)

          write(*,*)"ANGFAC     = ", angfac
          write(*,*)"XJBC(NB)   = ", xjbc(nb)
          write(*,*)"BRE        = ", bre
          write(*,*)"NB         = ", nb
          
          write(*,*)"IDIRECT    = ", idirect
          write(*,*)"IHADR      = ", ihadr

          write(*,*)"CCQQ       = ", ccqq
          write(*,*)"CCQG       = ", ccqg
          write(*,*)"CCQQSUD    = ", ccqqsud
          write(*,*)"CCQGSUD    = ", ccqgsud
          write(*,*)"CCQQ1      = ", ccqq1
          write(*,*)"CCQQ2      = ", ccqq2
          write(*,*)"CCQQ3      = ", ccqq3
          write(*,*)"CCQQ4      = ", ccqq4
          write(*,*)"CCQQ5      = ", ccqq5

          write(*,*)"PDGNQ      = ", pdgnq
          write(*,*)"PDGNG      = ", pdgng
          write(*,*)"PDGNQQ1    = ", pdgnqq1
          write(*,*)"PDGNQQ2    = ", pdgnqq2
          write(*,*)"PDGNQQ3    = ", pdgnqq3
          write(*,*)"PDGNQQ4    = ", pdgnqq4
          write(*,*)"PDGNQQ5    = ", pdgnqq5

          write(*,*)"SIGDGQ     = ", sigdgq
          write(*,*)"SIGDGG     = ", sigdgg
          write(*,*)"SIGDQQ1    = ", sigdqq1
          write(*,*)"SIGDQQ2    = ", sigdqq2
          write(*,*)"SIGDQQ3    = ", sigdqq3
          write(*,*)"SIGDQQ4    = ", sigdqq4
          write(*,*)"SIGDQQ5    = ", sigdqq5
c$$$
c$$$          write(*,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

c          write(*,*)"prod1      = ", PDGNQ    * SIGDGQ   * CCQQ
c          write(*,*)"prod2      = ", PDGNG    * SIGDGG   * CCQG
          write(*,*)"DGN        = ", dgn
          write(*,*)"------------------------------"
          endif
c          crap = .true.
          dgn =  0.d0
          wgt = -1.d0
       endif

C
       YSUM = YSUM + DGN 
C
 2     CONTINUE
C

c mbeekveld if contour is deformed we don't need prefactor
	IF (reg_div.eq.2) THEN
c and we do not need 2/sqrt(q^2)*pi from b integral
	 FQT = 1.D0/(2.D0/DSQRT(Q2) * PI) 
	ELSE 
       FQT = QT / (2.D0*PI)**2
	ENDIF

ccc   the total integrand
ccc   built-in security in case of "pathetic" events where 
ccc   the VEGAS integration weight becomes negative
c       write(*,*)"(2) FQT        = ", FQT
c       write(*,*)"--------"
       if(wgt.gt.0.d0) then
          BSTINT = YSUM / (1.D0-Z)**2 * FQT * FAC * XJ4 * XJ5 *
     1         2.D0/DSQRT(Q2) * PI 
     
       else
 77       BSTINT = 0.d0
       endif

c$$$       write(*,*)"DGN        = ", dgn
c$$$       write(*,*)"FQT        = ", fqt
c$$$       write(*,*)"BSTINT     = ", bstint
c$$$       write(*,*)"--------------------"

c$$$          write(*,*)"DGN        = ", dgn
c$$$          write(*,*)"FQT        = ", fqt
c$$$          write(*,*)"FAC        = ", fac
c$$$          write(*,*)"CCQQSUD    = ", ccqqsud
c$$$          write(*,*)"CCQGSUD    = ", ccqgsud
c$$$          write(*,*)"PDGNQ      = ", pdgnq
c$$$          write(*,*)"PDGNG      = ", pdgng
c$$$          write(*,*)"SIGDGQ     = ", sigdgq
c$$$          write(*,*)"SIGDGG     = ", sigdgg
c$$$          write(*,*)"CCQQ       = ", ccqq
c$$$          write(*,*)"CCQG       = ", ccqg
c$$$          write(*,*)"XJ4        = ", xj4
c$$$          write(*,*)"XJ5        = ", xj5
c$$$          write(*,*)"YSUM       = ", ysum
c$$$          write(*,*)"prod1      = ", PDGNQ    * SIGDGQ   * CCQQ
c$$$          write(*,*)"prod2      = ", PDGNG    * SIGDGG   * CCQG
c$$$          write(*,*)"BSTINT     = ", bstint
c$$$          write(*,*)"------------------------------"

       if(misnan(BSTINT))then
          write(*,*)"BSTINT     = NaN !"          
          write(*,*)"BSTINT     = ", bstint
          write(*,*)"DGN        = ", dgn 
          write(*,*)"FQT        = ", fqt 
          write(*,*)"SIGDGQ     = ", sigdgq
          write(*,*)"SIGDGG     = ", sigdgg
          write(*,*)"PDGNQ      = ", pdgnq
          write(*,*)"PDGNG      = ", pdgng
          write(*,*)"CCQQ       = ", ccqq
          write(*,*)"CCQG       = ", ccqg
          write(*,*)"QT         = ", qt
          write(*,*)"Q2         = ", q2
          write(*,*)"prod1      = ", PDGNQ    * SIGDGQ   * CCQQ
          write(*,*)"prod2      = ", PDGNG    * SIGDGG   * CCQG
          write(*,*)"------------------------------"
c          stop
      endif

       RETURN
       END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
