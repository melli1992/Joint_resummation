      SUBROUTINE ANCALCLNON(QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     1                    QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, 
     2                    C2QI, C2GF, CDYQI, CDYGI, XN)
       IMPLICIT DOUBLE COMPLEX (A - Z)
       DOUBLE PRECISION ZETA2, ZETA3, ALPS, ALP0, ALP1, ALPC, ALPB, 
     1                  ALPT, ALPQ, ALPQR, LMQ, LMQR, CF, CA, TR
       COMMON / COUPL  / ALPS, ALP0, ALP1, ALPC, ALPB, ALPT, ALPQ, 
     1                   ALPQR, LMQ, LMQR
C


       CF = 4.d0/3.d0
       TR = 1.d0/2.d0
       CA = 3.d0
C... Order(1/N), for P^(0):
       LNBAR = LOG(XN) + 0.577216
       QQI = 2.d0*CF*(-3.d0+ 4.d0*LNBAR+2.d0/XN)
       QGF = -8.d0*TR/XN
       GQI = -4.d0*CF/XN
       GGI = 8.d0*CA*LNBAR - 22.d0/3.d0*CA + 4.d0*CA/XN
       GGF = 8.d0/3.d0*TR

C...Order (1/N), for P^(1) :
       ZETA2 = 1.644934
       ZETA3 = 1.202057
C...NON-SINGLET PIECES AS GIVEN IN CURCI ET AL. (1980) :
c check EXP VAL
       NS1MA = (32.D0*LOG(XN)-6.518986723150949)/XN-21.220317620919317
c check EXP VAL
       NS1PA = NS1MA
c check EXP VAL
	NS1B = (16.D0*LNBAR-8.D0*ZETA2+424.D0/9.D0)/XN + 
     1       1.D0/18.D0*(1072.D0*LNBAR- 48.D0*6.D0*ZETA2*LNBAR
     1       -52.D0*6.D0*ZETA2- 129.D0)
c check EXP VAL
	NS1C = -176.D0/(9.D0*XN)-4.D0/9.D0*(40.D0*LNBAR-24.D0*ZETA2-3.D0)
c check EXP
       NS1MI = -2.D0/9.D0* NS1MA + 4.D0* NS1B
       NS1PI = -2.D0/9.D0* NS1PA + 4.D0* NS1B
       NS1F = 2.D0/3.D0 * NS1C

c check EXP VAL
       QQ1F = DCMPLX(0.D0,0.D0)
c check EXP VAL
       QG1A = -2.D0/XN*(LNBAR*LNBAR-1.D0)
c check EXP VAL
       QG1B = (2.D0*LNBAR*LNBAR + 5.D0-ZETA2*2.D0)/XN
c check EXP VAL
       QG1F = - 12.D0* QG1A - 16.D0/3.D0* QG1B
c check EXP VAL
       GQ1A = (-2.D0*LNBAR*LNBAR+10.D0*LNBAR-2.D0*ZETA2-12.D0)/XN
c check EXP VAL
       GQ1B = (LNBAR*LNBAR-17.D0/3.D0*LNBAR+109.D0/9.D0)/XN
c check EXP VAL
       GQ1C = (LNBAR-8.D0/3.D0)/XN
c check EXP
       GQ1I = - 64.D0/9.D0* GQ1A - 32.D0* GQ1B
c check EXP
       GQ1F = - 64.D0/9.D0* GQ1C
c check EXP VAL
       GG1A = -80.D0/(9.D0*XN)-32.D0/9.D0*(5.D0*LNBAR-3.D0)
c check EXP VAL
       GG1B = DCMPLX(8.D0,0.D0)
c check EXP VAL
       GG1C = (34.594262519841024 + 32.*LOG(XN))/XN + 
     1    2.6666666666666665*(12.463728932243974*LOG(XN) -
     2    11.624253271961711)
c check EXP
       GG1I = 9.D0* GG1C
c check EXP
       GG1F = 3.D0/2.D0* GG1A + 2.D0/3.D0* GG1B
C...WILSON COEFFICIENTS :
c check EXP VAL
       C2QI = 4.D0/9.D0*(6.D0*LNBAR*LNBAR+9.D0*LNBAR
     1        -6.D0*ZETA2-27.D0) +
     1  (8.D0/3.D0*LNBAR + 14.D0)/XN
c no idea
       C2QI = C2QI - LMQ * QQI/2.D0   
c check EXP VAL
       C2GF = (-2.D0*LNBAR-2.D0)/XN
c no idea
       C2GF = C2GF - LMQ * QGF/2.D0  

C... DRELL-YAN COEFFICIENTS :
c CHECK VAL
       CDYQI = 16.D0/9.D0*(3.D0*LNBAR*LNBAR+6.D0*ZETA2-6.D0)+ 
     1   16.D0*LNBAR/(3.D0*XN) +
     1    LMQ*(-8.D0/(3.D0*XN)-16.D0/3.D0*LNBAR+4.D0)
c CHECK VAL
       CDYGI = -LNBAR/XN + LMQ/(2.D0*XN)
       RETURN
       END


c for NLL only
       SUBROUTINE ANCALCNLL(QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     1                    QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, 
     2                    C2QI, C2GF, CDYQI, CDYGI, XN)
       IMPLICIT DOUBLE COMPLEX (A - Z)
       DOUBLE PRECISION ZETA2, ZETA3, ALPS, ALP0, ALP1, ALPC, ALPB, 
     1                  ALPT, ALPQ, ALPQR, LMQ, LMQR
       COMMON / COUPL  / ALPS, ALP0, ALP1, ALPC, ALPB, ALPT, ALPQ, 
     1                   ALPQR, LMQ, LMQR
C
       CF = 4.d0/3.d0
       TR = 1.d0/2.d0
       CA = 3.d0
C... Order(1), for P^(0):
       LNBAR = LOG(XN) + 0.577216
       QQI = 2.d0*CF*(-3.d0+ 4.d0*LNBAR)
       QGF = DCMPLX(0.D0,0.D0)
       GQI = DCMPLX(0.D0,0.D0)
       GGI = 8.d0*CA*LNBAR - 22.d0/3.d0*CA
       GGF = 8.d0/3.d0*TR

C...Order (1), for P^(1) :
       ZETA2 = 1.644934
       ZETA3 = 1.202057
C...NON-SINGLET PIECES AS GIVEN IN CURCI ET AL. (1980) :
c check EXP VAL
c mbeekveld       NS1MA = -21.220317620919317
	NS1MA = DCMPLX(0.D0,0.D0)
c check EXP VAL
       NS1PA = NS1MA
c check EXP VAL
	NS1B = 1.D0/18.D0*(1072.D0*LNBAR- 48.D0*6.D0*ZETA2*LNBAR)
c mbeekveld     2     -52.D0*6.D0*ZETA2- 129.D0)
c check EXP VAL
	NS1C = -4.D0/9.D0*(40.D0*LNBAR)
c mbeekveld -24.D0*ZETA2-3.D0)
c check EXP
       NS1MI = -2.D0/9.D0* NS1MA + 4.D0* NS1B
       NS1PI = -2.D0/9.D0* NS1PA + 4.D0* NS1B
       NS1F = 2.D0/3.D0 * NS1C

c check EXP VAL
       QQ1F = DCMPLX(0.D0,0.D0)
c check EXP VAL
       QG1A = DCMPLX(0.D0,0.D0)
c check EXP VAL
       QG1B = DCMPLX(0.D0,0.D0)
c check EXP VAL
       QG1F = - 12.D0* QG1A - 16.D0/3.D0* QG1B
c check EXP VAL
       GQ1A = DCMPLX(0.D0,0.D0)
c check EXP VAL
       GQ1B = DCMPLX(0.D0,0.D0)
c check EXP VAL
       GQ1C = DCMPLX(0.D0,0.D0)
c check EXP
       GQ1I = - 64.D0/9.D0* GQ1A - 32.D0* GQ1B
c check EXP
       GQ1F = - 64.D0/9.D0* GQ1C
c check EXP VAL
       GG1A = -32.D0/9.D0*(5.D0*LNBAR-3.D0)
c check EXP VAL
       GG1B = DCMPLX(8.D0,0.D0)
c check EXP VAL
       GG1C = 2.6666666666666665D0*(12.463728932243974D0*LOG(XN))
c mbeekveld - 11.624253271961711D0)
c check EXP
       GG1I = 9.D0* GG1C
c check EXP
       GG1F = 3.D0/2.D0* GG1A + 2.D0/3.D0* GG1B
C...WILSON COEFFICIENTS :
c check EXP VAL
       C2QI = 4.D0/9.D0*(6.D0*LNBAR*LNBAR+9.D0*LNBAR-6.D0*ZETA2-27.D0)
c no idea
       C2QI = C2QI - LMQ * QQI/2.D0   
c check EXP VAL
       C2GF = DCMPLX(0.D0,0.D0)
c no idea
       C2GF = C2GF - LMQ * QGF/2.D0  

C... DRELL-YAN COEFFICIENTS :
c CHECK VAL
       CDYQI = 16.D0/9.D0*(3.D0*LNBAR*LNBAR+6.D0*ZETA2-6.D0) +
     1    LMQ*(-16.D0/3.D0*LNBAR+4.D0)
c CHECK VAL
       CDYGI = DCMPLX(0.D0,0.D0)
      
       RETURN
       END

c DIAGONAL LN(N)/N 
       SUBROUTINE ANCALCNLLON(QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, 
     1                    NS1F, QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, 
     2                    C2QI, C2GF, CDYQI, CDYGI, XN)
       IMPLICIT DOUBLE COMPLEX (A - Z)
       DOUBLE PRECISION ZETA2, ZETA3, ALPS, ALP0, ALP1, ALPC, ALPB, 
     1                  ALPT, ALPQ, ALPQR, LMQ, LMQR
       COMMON / COUPL  / ALPS, ALP0, ALP1, ALPC, ALPB, ALPT, ALPQ, 
     1                   ALPQR, LMQ, LMQR
C
       CF = 4.d0/3.d0
       TR = 1.d0/2.d0
       CA = 3.d0
C... Order(1/N), for P^(0):
       LNBAR = LOG(XN) + 0.577216
       QQI = 2.d0*CF*(-3.d0+ 4.d0*LNBAR+2.d0/XN)
       QGF = DCMPLX(0.D0,0.D0)
       GQI = DCMPLX(0.D0,0.D0)
       GGI = 8.d0*CA*LNBAR - 22.d0/3.d0*CA + 4.d0*CA/XN
       GGF = 8.d0/3.d0*TR

C...Order (1/N), for P^(1) :
       ZETA2 = 1.644934
       ZETA3 = 1.202057
C...NON-SINGLET PIECES AS GIVEN IN CURCI ET AL. (1980) :
c check EXP VAL
c mbeekveld changed
c       NS1MA = -21.220317620919317
	NS1MA = DCMPLX(0.D0,0.D0)
c check EXP VAL
       NS1PA = NS1MA
c check EXP VAL
	NS1B = 1.D0/18.D0*(1072.D0*LNBAR- 48.D0*6.D0*ZETA2*LNBAR)
c mbeekveld changed
c     2     -52.D0*6.D0*ZETA2- 129.D0)
c check EXP VAL
	NS1C = -4.D0/9.D0*(40.D0*LNBAR)
c mbeekveld add this
c-24.D0*ZETA2-3.D0)
c check EXP
       NS1MI = -2.D0/9.D0* NS1MA + 4.D0* NS1B
       NS1PI = -2.D0/9.D0* NS1PA + 4.D0* NS1B
       NS1F = 2.D0/3.D0 * NS1C

c check EXP VAL
       QQ1F = DCMPLX(0.D0,0.D0)
c check EXP VAL
       QG1A = DCMPLX(0.D0,0.D0)
c check EXP VAL
       QG1B = DCMPLX(0.D0,0.D0)
c check EXP VAL
       QG1F = - 12.D0* QG1A - 16.D0/3.D0* QG1B
c check EXP VAL
       GQ1A = DCMPLX(0.D0,0.D0)
c check EXP VAL
       GQ1B = DCMPLX(0.D0,0.D0)
c check EXP VAL
       GQ1C = DCMPLX(0.D0,0.D0)
c check EXP
       GQ1I = - 64.D0/9.D0* GQ1A - 32.D0* GQ1B
c check EXP
       GQ1F = - 64.D0/9.D0* GQ1C
c check EXP VAL
       GG1A = -32.D0/9.D0*(5.D0*LNBAR-3.D0)
c check EXP VAL
       GG1B = DCMPLX(8.D0,0.D0)
c check EXP VAL
       GG1C = 2.6666666666666665D0*(12.463728932243974D0*LOG(XN))
cmbeekveld     2    -11.624253271961711D0)
c check EXP
       GG1I = 9.D0* GG1C
c check EXP
       GG1F = 3.D0/2.D0* GG1A + 2.D0/3.D0* GG1B
C...WILSON COEFFICIENTS :
c check EXP VAL
       C2QI = 4.D0/9.D0*(6.D0*LNBAR*LNBAR+9.D0*LNBAR-6.D0*ZETA2-27.D0)
c no idea
       C2QI = C2QI - LMQ * QQI/2.D0   
c check EXP VAL
       C2GF = DCMPLX(0.D0,0.D0)
c no idea
       C2GF = C2GF - LMQ * QGF/2.D0  

C... DRELL-YAN COEFFICIENTS :
c CHECK VAL
       CDYQI = 16.D0/9.D0*(3.D0*LNBAR*LNBAR+6.D0*ZETA2-6.D0) +
     1    LMQ*(-16.D0/3.D0*LNBAR+4.D0)
c CHECK VAL
       CDYGI = DCMPLX(0.D0,0.D0)
      
       RETURN
       END


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C...EVOLUTION MATRIX BETWEEN TWO COMPLEX-VALUED SCALES. USED TO 
C...IMPLEMENT LEADING LN(N)/N CORRECTIONS IN JOINT RESUMMATION
C...A LA KULESZA-STERMAN-VOGELSANG

      SUBROUTINE EVOLMAT(XN,CMU0,CMU1,ENSP,ENSM,EFF,EFG,EGF,EGG)

      IMPLICIT DOUBLE COMPLEX(A-Z)
      INTEGER I,J,ETA,F, FLAV

      COMMON / FLAVORS / FLAV


      CMU02 = CMU0*CMU0
      CMU12 = CMU1*CMU1

      F    = FLAV

      ALP0 = ALPHASC_MSTW(CMU0)
      ALP  = ALPHASC_MSTW(CMU1)
       IF( isnan(REAL(ALP))) THEN
       ALP = DCMPLX(1.0,0.0)
        ENDIF
      CALL ANCALC (QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     1              QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, C2QI, C2GF, 
     2              CDYQI, CDYGI, XN)

      
      XL = ALP0/ALP
      XL1 = 1.D0-XL


      CALL ANOM (ANS, AM, AP, AL, BE, AB, RMIN, RPLUS, RQQ, RQG,
     1            RGQ, RGG, C2Q, C2G, CDYQ, CDYG, XN, F,
     2            QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     3            QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, C2QI, C2GF,
     4            CDYQI, CDYGI)


C
       S     = LOG(XL)
       ENS   = EXP (-ANS*S)
C
       EM    = EXP (-AM*S)
       EP    = EXP (-AP*S)
       AC    = 1.D0- AL
       EMP   = EXP (S * (-AM+AP))
       EPM   = EXP (S * (+AM-AP))
       NMP   = 1.D0- AM + AP
       NPM   = 1.D0- AP + AM
       DMQQ  =  AL * RQQ + BE * RGQ
       DMQG  =  AL * RQG + BE * RGG
       DMGQ  =  AB * RQQ + AC * RGQ
       DMGG  =  AB * RQG + AC * RGG
       DPQQ  =  AC * RQQ - BE * RGQ
       DPQG  =  AC * RQG - BE * RGG
       DPGQ  = -AB * RQQ + AL * RGQ
       DPGG  = -AB * RQG + AL * RGG
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
       ENSM  =  ENS * (1.+ ALP * XL1 * RMIN)
       ENSP  =  ENS * (1.+ ALP * XL1 * RPLUS)
       EFF   =  EM * (AL + ALP * (RMMQQ * XL1 + RMPQQ * (EPM-XL))) +
     x          EP * (AC + ALP * (RPPQQ * XL1 + RPMQQ * (EMP-XL))) 
       EFG   =  EM * (BE + ALP * (RMMQG * XL1 + RMPQG * (EPM-XL))) +
     x          EP * (-BE + ALP * (RPPQG * XL1 + RPMQG * (EMP-XL)))
       EGF   =  EM * (AB + ALP * (RMMGQ * XL1 + RMPGQ * (EPM-XL))) +
     x          EP * (-AB + ALP * (RPPGQ * XL1 + RPMGQ * (EMP-XL)))
       EGG   =  EM * (AC + ALP * (RMMGG * XL1 + RMPGG * (EPM-XL))) +
     x          EP * (AL + ALP * (RPPGG * XL1 + RPMGG * (EMP-XL))) 


      RETURN
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      SUBROUTINE EVOLMATLNON(XN,CMU0,CMU1,ENSP,ENSM,EFF,EFG,EGF,EGG,
     1  ENS_LL)

      IMPLICIT DOUBLE COMPLEX(A-Z)
      INTEGER I,J,ETA,F, FLAV

      COMMON /FLAVORS/ FLAV


      CMU02 = CMU0*CMU0
      CMU12 = CMU1*CMU1
      F    = FLAV
      ALP0 = ALPHASC_MSTW(CMU0)
      ALP  = ALPHASC_MSTW(CMU1)
       IF( isnan(REAL(ALP))) THEN
       ALP = DCMPLX(1.0,0.0)
        ENDIF

      CALL ANCALCLNON (QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     1              QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, C2QI, C2GF, 
     2              CDYQI, CDYGI, XN)

      
      XL = ALP0/ALP
      XL1 = 1.-XL

      CALL ANOM (ANS, AM, AP, AL, BE, AB, RMIN, RPLUS, RQQ, RQG,
     1            RGQ, RGG, C2Q, C2G, CDYQ, CDYG, XN, F,
     2            QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     3            QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, C2QI, C2GF,
     4            CDYQI, CDYGI)
      PI=DACOS(-1.D0)

      B0 = (11.- 2./3.* F)/12./PI
C
       S = LOG(XL)
       ENS = EXP (-ANS*S)
       ENS_LL = ENS
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
       ENSM  =  ENS * (1.+ ALP * XL1 * RMIN)
       ENSP  =  ENS * (1.+ ALP * XL1 * RPLUS)
       EFF  = EM * (AL + ALP * (RMMQQ * XL1 + RMPQQ * (EPM-XL))) +
     #        EP * (AC + ALP * (RPPQQ * XL1 + RPMQQ * (EMP-XL))) 
       EFG  = EM * (BE + ALP * (RMMQG * XL1 + RMPQG * (EPM-XL))) +
     #        EP * (-BE + ALP * (RPPQG * XL1 + RPMQG * (EMP-XL)))
       EGF  = EM * (AB + ALP * (RMMGQ * XL1 + RMPGQ * (EPM-XL))) +
     #        EP * (-AB + ALP * (RPPGQ * XL1 + RPMGQ * (EMP-XL)))
       EGG  = EM * (AC + ALP * (RMMGG * XL1 + RMPGG * (EPM-XL))) +
     #        EP * (AL + ALP * (RPPGG * XL1 + RPMGG * (EMP-XL))) 

      RETURN
      END

      SUBROUTINE EVOLMATNLL(XN,CMU0,CMU1,ENSP,ENSM,EFF,EFG,EGF,EGG,
     1  ENS_LL)

      IMPLICIT DOUBLE COMPLEX(A-Z)
      INTEGER I,J,ETA,F, FLAV

      COMMON /FLAVORS/ FLAV


      CMU02 = CMU0*CMU0
      CMU12 = CMU1*CMU1
      F    = FLAV
      ALP0 = ALPHASC_MSTW(CMU0)
      ALP  = ALPHASC_MSTW(CMU1)
      ALP0_LO = ALPHASC_LO(CMU0)
      ALP_LO  = ALPHASC_LO(CMU1)
       IF( isnan(REAL(ALP))) THEN
       ALP = DCMPLX(1.0,0.0)
        ENDIF
       IF( isnan(REAL(ALP_LO))) THEN
       ALP_LO = DCMPLX(1.0,0.0)
        ENDIF

      CALL ANCALCNLL (QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     1              QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, C2QI, C2GF, 
     2              CDYQI, CDYGI, XN)

      
      XL = ALP0/ALP
      XL1 = 1.-XL

      CALL ANOM (ANS, AM, AP, AL, BE, AB, RMIN, RPLUS, RQQ, RQG,
     1            RGQ, RGG, C2Q, C2G, CDYQ, CDYG, XN, F,
     2            QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     3            QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, C2QI, C2GF,
     4            CDYQI, CDYGI)
      PI=DACOS(-1.D0)

      B0 = (11.- 2./3.* F)/12./PI
C
       S = LOG(XL)
       ENS = EXP (-ANS*S)
       EM  = EXP (-AM*S)
       EP  = EXP (-AP*S)
C...EVOLUTION OF LIGHT PARTON DESITIES BETWEEN THRESHOLDS :
       ENSM  =  ENS * EXP((ALP_LO-ALP0_LO) * RMIN)
       ENSP  =  ENS * EXP((ALP_LO-ALP0_LO)* RPLUS)
       EFF  = EM * EXP((ALP_LO-ALP0_LO) * RQQ)
       EFG  = DCMPLX(0.D0,0.D0)
       EGF  = DCMPLX(0.D0,0.D0)
       EGG  = EP * EXP((ALP_LO-ALP0_LO) * RGG)



      RETURN
      END

      SUBROUTINE EVOLMATNLLON(XN,CMU0,CMU1,ENSP,ENSM,EFF,EFG,EGF,EGG,
     1  ENS_LL)

      IMPLICIT DOUBLE COMPLEX(A-Z)
      INTEGER I,J,ETA,F, FLAV

      COMMON /FLAVORS/ FLAV


      CMU02 = CMU0*CMU0
      CMU12 = CMU1*CMU1
      F    = FLAV
      ALP0 = ALPHASC_MSTW(CMU0)
      ALP  = ALPHASC_MSTW(CMU1)
      ALP0_LO = ALPHASC_LO(CMU0)
      ALP_LO  = ALPHASC_LO(CMU1)
       IF( isnan(REAL(ALP))) THEN
       ALP = DCMPLX(1.0,0.0)
        ENDIF
       IF( isnan(REAL(ALP_LO))) THEN
       ALP_LO = DCMPLX(1.0,0.0)
        ENDIF

      CALL ANCALCNLLON (QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     1              QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, C2QI, C2GF, 
     2              CDYQI, CDYGI, XN)

      
      XL = ALP0/ALP
      XL1 = 1.-XL

      CALL ANOM (ANS, AM, AP, AL, BE, AB, RMIN, RPLUS, RQQ, RQG,
     1            RGQ, RGG, C2Q, C2G, CDYQ, CDYG, XN, F,
     2            QQI, QGF, GQI, GGI, GGF, NS1MI, NS1PI, NS1F,
     3            QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, C2QI, C2GF,
     4            CDYQI, CDYGI)
      PI=DACOS(-1.D0)

      B0 = (11.- 2./3.* F)/12./PI
C
       S = LOG(XL)
       ENS = EXP (-ANS*S)
       EM  = EXP (-AM*S)
       EP  = EXP (-AP*S)
C...EVOLUTION OF LIGHT PARTON DESITIES BETWEEN THRESHOLDS :
       ENSM  =  ENS * EXP((ALP_LO-ALP0_LO) * RMIN)
       ENSP  =  ENS * EXP((ALP_LO-ALP0_LO)* RPLUS)
       EFF  = EM * EXP((ALP_LO-ALP0_LO) * RQQ)
       EFG  = DCMPLX(0.D0,0.D0)
       EGF  = DCMPLX(0.D0,0.D0)
       EGG  = EP * EXP((ALP_LO-ALP0_LO) * RGG)



      RETURN
      END


