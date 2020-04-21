C                                                                               
C...PHOTON FRAGMENTATION IN NEXT TO LEADING ORDER                               
C...(POINTLIKE OR WITH VMD INPUT) FOR 3-6 ACTIVE FLAVOURS                       
C                                                                               
       SUBROUTINE PHFNLO(XN, Q2, UN, DN, SN, CN, GLN)
       IMPLICIT DOUBLE PRECISION (A - Z) 
       DOUBLE COMPLEX UN, DN, SN, CN, GLN, XN
       DIMENSION LAMBDA(3:6), Q2THR(3)                          
       INTEGER N
       COMMON / SCALESP / Q2S, Q20, Q2THR, LAMBDA 
       COMMON / SCHEMEP / N 
       COMMON / COUPLP / ALPS, ALP0, ALPC, ALPB, ALPT, ALPQ 
       DATA LAMBDA / 0.248D0, 0.200D0, 0.131D0, 0.053D0 / 
       DATA Q2S, Q20 / 0.3D0, 4.D0/
       DATA Q2THR / 2.25d0, 20.25d0, 30625.d0 / 
cpm       DATA Q2S, Q20 / 0.4D0, 4.D0/       
cpm     1      Q2THR / 2.25, 20.25, 1.D8 / 
C...CHOICE OF THE FACTORIZATION SCHEME: N = 2: MS-BAR, 1: DIS-GAMMA(F1)         
       DATA N / 2 / 
C !!!  NOTE: N=2 MEANS EVOLUTION IS DONE IN DIS-GAMMA. NEVERTHELESS
C      TRANSFORMATION TO MSBAR IS DONE LATER !!!  
C                 
C...COUPLING CONSTANTS AT STATIC POINT, INPUT SCALE AND THRESHOLDS :            
       ALPS = ALPHASP (Q2S)
       ALP0 = ALPHASP (Q20)
       ALPC = ALPHASP (Q2THR(1))
       ALPB = ALPHASP (Q2THR(2))
       ALPT = ALPHASP (Q2THR(3))
       ALPQ = ALPHASP (Q2)
C
       CALL RENOP (UN, DN, SN, CN, GLN, XN) 
C
       RETURN
       END
C
C                                                                               
C...CALCULATION OF ALPHA STRONG (Q**2) DIVIDED BY 4 * PI :                      
       FUNCTION ALPHASP (Q2)

       IMPLICIT DOUBLE PRECISION (A - Z)

       INTEGER NF, K
       DIMENSION LAMBDA (3:6), Q2THR (3)
       COMMON / SCALESP / Q2S, Q20, Q2THR, LAMBDA
 
       NF = 3
       DO 10 K = 1, 3
       IF (Q2 .GT. Q2THR (K)) THEN
          NF = NF + 1
       ELSE
          GO TO 20
       END IF
  10   CONTINUE
  20   B0  = 11.- 2./3.* NF
       B0S = B0 * B0
       B1  = 102.- 38./3.* NF
       LAM2 = LAMBDA (NF) * LAMBDA (NF)
       LQ2  = DLOG (Q2 / LAM2)
       ALPHASP = 1./ (B0 * LQ2) * (1.- B1 / B0S * DLOG (LQ2) / LQ2)
 
       RETURN
       END
C                                                                               
C...MELLIN MOMENTS OF THE FRAGMENTATION DISTRIBUTIONS                           
C...AND STRUCTURE FUNCTIONS :
       SUBROUTINE RENOP (UN, DN, SN, CN, GLN, XN)
       IMPLICIT DOUBLE COMPLEX (A - Z)
       INTEGER N, F
       DOUBLE PRECISION ESQ, EFO, EFS, ALP, ALPS, ALP0, ALPC, ALPB,
     1                  ALPT, ALPQ, XL, S, XL1, ES, PI
       COMMON / SCHEMEP / N
       COMMON / COUPLP / ALPS, ALP0, ALPC, ALPB, ALPT, ALPQ
 
       CALL FRCALCP (QQI, QGF, GQI, GGI, GGF,        NS1PI, NS1F,
     1            QQ1F, QG1F, GQ1F, GQ1I, GG1I, GG1F, QPF, QP1F,
     2            GP1F, C2QI, C1QI, C2GF, C1GF, C2PI, C1PI, XN)
       PI = DACOS(-1.D0)

       GO TO 20
C...PARAMETERS FOR NLO EVOLUTION :
  10   CONTINUE
       CALL ANOMP (ANS, AM, AP, AL, BE, AB,       RPLUS, RQQ, RQG,
     1            RGQ, RGG, KNS, KSI, RKNS, RKSI, RKGL, C2Q, C1Q, C2G,
     2            C1G, C2P, C1P, ESQ, EFO, EFS,
     3            QQI, QGF, GQI, GGI, GGF,        NS1PI, NS1F, QQ1F,
     4            QG1F, GQ1F, GQ1I, GG1I, GG1F, QPF, QP1F, GP1F, C2QI,
     5            C1QI, C2GF, C1GF, C2PI, C1PI, XN, F)
       AC  = 1.- AL
       S   = DLOG (XL)
       XL1 = 1.- XL
       ES  = EXP (-S)
       ENS = EXP (-ANS*S)
       EM  = EXP (-AM*S)
       EP  = EXP (-AP*S)
       EMPL = EXP (S * (-AM+AP)) - XL
       EPML = EXP (S * (+AM-AP)) - XL
       ENS1 = ENS * ES
       EM1  = EM * ES
       EP1  = EP * ES
       EENS = 1.- ENS
       EEN1 = 1.- ENS1
       EEM  = 1.- EM
       EEM1 = 1.- EM1
       EEP  = 1.- EP
       EEP1 = 1.- EP1
       ANS1 = 1.+ ANS
       AM1  = 1.+ AM
       AP1  = 1.+ AP
       NMP  = 1.- AM + AP
       NPM  = 1.- AP + AM
       DMQQ = AL * RQQ + BE * RGQ
       DMQG = AL * RQG + BE * RGG
       DMGQ = AB * RQQ + AC * RGQ
       DMGG = AB * RQG + AC * RGG
       DPQQ = AC * RQQ - BE * RGQ
       DPQG = AC * RQG - BE * RGG
       DPGQ = - AB * RQQ + AL * RGQ
       DPGG = - AB * RQG + AL * RGG
       RMMQQ = AL * DMQQ + AB * DMQG
       RMMQG = BE * DMQQ + AC * DMQG
       RMMGQ = AL * DMGQ + AB * DMGG
       RMMGG = BE * DMGQ + AC * DMGG
       RMPQQ = (AC * DMQQ - AB * DMQG) / NMP
       RMPQG = (-BE * DMQQ + AL * DMQG) / NMP
       RMPGQ = (AC * DMGQ - AB * DMGG) / NMP
       RMPGG = (-BE * DMGQ + AL * DMGG) / NMP
       RPMQQ = (AL * DPQQ + AB * DPQG) / NPM
       RPMQG = (BE * DPQQ + AC * DPQG) / NPM
       RPMGQ = (AL * DPGQ + AB * DPGG) / NPM
       RPMGG = (BE * DPGQ + AC * DPGG) / NPM
       RPPQQ = AC * DPQQ - AB * DPQG
       RPPQG = -BE * DPQQ + AL * DPQG
       RPPGQ = AC * DPGQ - AB * DPGG
       RPPGG = -BE * DPGQ + AL * DPGG
C...INHOMOGENEOUS SOLUTION :
       INSL = EEN1 / ANS1 * KNS / ALP
       IN3L = INSL / (3.* F * EFS)
       IN8L = IN3L
       ISIL = (EEM1 * AL / AM1  +  EEP1 * AC / AP1) * KSI / ALP
       IGLL = (EEM1 / AM1  -  EEP1 / AP1) * AB * KSI / ALP
       INSH =  EEN1 / ANS1 * RPLUS * KNS
     1       + EENS / ANS * (-RPLUS * KNS + RKNS)
       IN3H = INSH / (3.* F * EFS)
       IN8H = IN3H
       ISIH =  (EEM1 / AM1 * RMMQQ  +  EEP1 / AP1 * RPPQQ) * KSI
     2        + EEM / AM * (-RMMQQ * KSI + AL * RKSI + BE * RKGL)
     3        + EEP / AP * (-RPPQQ * KSI + AC * RKSI - BE * RKGL)
     4        + (EEP1 / AP1 - EEM / AM) * RMPQQ * KSI
     5        + (EEM1 / AM1 - EEP / AP) * RPMQQ * KSI
       IGLH =  (EEM1 / AM1 * RMMGQ  +  EEP1 / AP1 * RPPGQ) * KSI
     2        + EEM / AM * (-RMMGQ * KSI + AB * RKSI + AC * RKGL)
     3        + EEP / AP * (-RPPGQ * KSI - AB * RKSI + AL * RKGL)
     4        + (EEP1 / AP1 - EEM / AM) * RMPGQ * KSI
     5        + (EEM1 / AM1 - EEP / AP) * RPMGQ * KSI
C...HOMOGENEOUS SOLUTION :
       HN3L = NS3L * ENS
       HN8L = NS8L * ENS
       HSIL = (EM * AL + EP * AC) * SIL  +  (EM - EP) * BE * GLL
       HGLL = (EM - EP) * AB * SIL  +  (EM * AC + EP * AL) * GLL
       HN3H = NS3L * ENS * XL1 * RPLUS * ALP + NS3H * ENS
       HN8H = NS8L * ENS * XL1 * RPLUS * ALP + NS8H * ENS
       HSIH = (EM * ((RMMQQ * XL1 + RMPQQ * EPML) * SIL
     1              + (RMMQG * XL1 + RMPQG * EPML) * GLL)
     2        + EP * ((RPPQQ * XL1 + RPMQQ * EMPL) * SIL
     3              + (RPPQG * XL1 + RPMQG * EMPL) * GLL)) * ALP
     4        + (EM * AL + EP * AC) * SIH  + (EM - EP) * BE * GLH
       HGLH  = (EM * ((RMMGQ * XL1 + RMPGQ * EPML) * SIL
     1              + (RMMGG * XL1 + RMPGG * EPML) * GLL)
     2        + EP * ((RPPGQ * XL1 + RPMGQ * EMPL) * SIL
     3              + (RPPGG * XL1 + RPMGG * EMPL) * GLL)) * ALP
     4        + (EM - EP) * AB * SIH + (EM * AC + EP * AL) * GLH
C...COMPLETE SOLUTION :
       NS3L = IN3L + HN3L
       NS8L = IN8L + HN8L
       SIL  = ISIL + HSIL
       GLL  = IGLL + HGLL
       NS3H = IN3H + HN3H
       NS8H = IN8H + HN8H
       SIH  = ISIH + HSIH
       GLH  = IGLH + HGLH
C
       IF ( ALP .GE. ALPC ) THEN
          GO TO 30
       ELSE IF ( ALP .GE. ALP0 ) THEN
          GO TO 40
       ELSE IF ( ALP .GE. ALPB ) THEN
          GO TO 45
       ELSE IF ( ALP .GE. ALPT ) THEN
          GO TO 50
       ELSE
          GO TO 60
       END IF
C...EVOLUTION BELOW THE CHARM THRESHOLD (POINTLIKE) :
  20   CONTINUE
       F = 3
       NS3L = 0.D0
       NS8L = 0.D0
       SIL  = 0.D0
       GLL  = 0.D0
       NS3H = 0.D0
       NS8H = 0.D0
       SIH  = 0.D0
       GLH  = 0.D0
       IF ( ALPQ .GE. ALPC ) THEN
          ALP = ALPQ
       ELSE
          ALP = ALPC
       END IF
       XL = ALPS / ALP
       GO TO 10
  30   CONTINUE
       NS15L = SIL
       NS24L = SIL
       NS35L = SIL
       NS15H = SIH
       NS24H = SIH
       NS35H = SIH
       IF ( ALPQ .GE. ALPC ) GO TO 70
C...EVOLUTION BETWEEN CHARM AND INPUT SCALE :
       F = 4
       IF ( ALPQ .GE. ALP0 ) THEN
          ALP = ALPQ
       ELSE
          ALP = ALP0
       END IF
       XL = ALPC / ALP
       GO TO 10
  40   CONTINUE
       HN15L = NS15L * ENS
       NS15L = HN15L - 2.* IN3L
       NS24L = SIL 
       NS35L = SIL
       HN15H = NS15L * ENS * XL1 * RPLUS * ALP + NS15H * ENS
       NS15H = HN15H - 2.* IN3H 
       NS24H = SIH
       NS35H = SIH
       IF ( ALPQ .GT. ALP0 ) GO TO 70
C...INPUT MOMENTS FOR THE HADRONIC PART (MODIFIED OWENS-INPUT) :
       CALL FRHOIN (VAR, XIR, GLR, XN)
       KAP = 1.D0 
C       KAP = 0.  
       SIR  = (2.* VAR + 6.* XIR) * KAP / 2.2D0
       GLR  = GLR * KAP / 2.2D0
       NS3R = 0.0D0
       NS8R = (2.D0* VAR) * KAP / 2.2D0
       NS8L = NS8L + NS8R
       SIL  = SIL + SIR
       GLL  = GLL + GLR
       NS15L = NS15L + SIR
       NS24L = SIL
       NS35L = SIL
C...EVOLUTION BETWEEN INPUT SCALE AND BOTTOM THRESHOLD :
       F = 4
       IF ( ALPQ .GE. ALPB ) THEN
          ALP = ALPQ
       ELSE
          ALP = ALPB
       END IF
       XL = ALP0 / ALP
       GO TO 10
  45   CONTINUE
       HN15L = NS15L * ENS
       NS15L = HN15L - 2.* IN3L
       NS24L = SIL
       NS35L = SIL
       HN15H = NS15L * ENS * XL1 * RPLUS * ALP + NS15H * ENS
       NS15H = HN15H - 2.* IN3H
       NS24H = SIH
       NS35H = SIH
       IF ( ALPQ .GE. ALPB ) GO TO 70
C...EVOLUTION BETWEEN BOTTOM AND TOP THRESHOLD :
       F = 5
       IF ( ALPQ .GE. ALPT ) THEN
          ALP = ALPQ
       ELSE
          ALP = ALPT
       END IF
       XL = ALPB / ALP
       GO TO 10
  50   CONTINUE
       HN15L = NS15L * ENS
       NS15L = HN15L - 2.* IN3L
       HN24L = NS24L * ENS
       NS24L = HN24L + 2.* IN3L
       NS35L = SIL
       HN15H = NS15L * ENS * XL1 * RPLUS * ALP + NS15H * ENS
       NS15H = HN15H - 2.* IN3H
       HN24H = NS24L * ENS * XL1 * RPLUS * ALP + NS24H * ENS
       NS24H = HN24H + 2.* IN3H
       NS35H = SIH
       IF ( ALPQ .GE. ALPT ) GO TO 70
C...EVOLUTION ABOVE THE TOP THRESHOLD :
       F = 6
       ALP = ALPQ
       XL = ALPT / ALP
       GO TO 10
  60   CONTINUE
       HN15L = NS15L * ENS
       NS15L = HN15L - 2.* IN3L
       HN24L = NS24L * ENS
       NS24L = HN24L + 2.* IN3L
       HN35L = NS35L * ENS
       NS35L = HN35L - 3.* IN3L
       HN15H = NS15L * ENS * XL1 * RPLUS * ALP + NS15H * ENS
       NS15H = HN15H - 2.* IN3H
       HN24H = NS24L * ENS * XL1 * RPLUS * ALP + NS24H * ENS
       NS24H = HN24H + 2.* IN3H
       HN35H = NS35L * ENS * XL1 * RPLUS * ALP + NS35H * ENS
       NS35H = HN24H - 3.* IN3H
C...FLAVOUR DECOMPOSITION :
  70   CONTINUE
       TNL = (SIL - NS35L) / 12.D0
       BNL = (5.* SIL + NS35L - 6.* NS24L) / 60.D0
       CNL = (10.* SIL + 2.* NS35L + 3.* NS24L - 15.* NS15L ) / 120.D0
       SNL = (10.* SIL + 2.* NS35L + 3.* NS24L + 5.* NS15L - 20.* NS8L) 
     1      / 120.D0
       DNL = (10.* SIL + 2.* NS35L + 3.* NS24L + 5.* NS15L + 10.* NS8L
     1      - 30.* NS3L) / 120.D0
       UNL = (10.* SIL + 2.* NS35L + 3.* NS24L + 5.* NS15L + 10.* NS8L
     1      + 30.* NS3L) / 120.D0
       TNH = (SIH - NS35H) / 12.D0       
       BNH = (5.* SIH + NS35H - 6.* NS24H) / 60.D0
       CNH = (10.* SIH + 2.* NS35H + 3.* NS24H - 15.* NS15H ) / 120.D0
       SNH = (10.* SIH + 2.* NS35H + 3.* NS24H + 5.* NS15H - 20.* NS8H)
     1       / 120.D0
       DNH = (10.* SIH + 2.* NS35H + 3.* NS24H + 5.* NS15H + 10.* NS8H
     1        - 30.* NS3H) / 120.D0
       UNH = (10.* SIH + 2.* NS35H + 3.* NS24H + 5.* NS15H + 10.* NS8H
     1        + 30.* NS3H) / 120.D0
       TN  = TNL + TNH
       BN  = BNL + BNH
       CN  = CNL + CNH
       SN  = SNL + SNH
       DN  = DNL + DNH
       UN  = UNL + UNH
       SIN = SIL + SIH
       GLN = GLL + GLH
       QNSN = 0.0
       TNL = 0.
       TNH = 0.
C       BNL = 0.
C       BNH = 0.
C...TRANSFORMATION TO THE MS(BAR)-SCHEME:                                       
       IF(N.EQ.2) THEN
        UN  = UN - 1./2./12.566 * 4./9.* C1P
        DN  = DN - 1./2./12.566 * 1./9.* C1P
        SN  = SN - 1./2./12.566 * 1./9.* C1P
        CN  = CN - 1./2./12.566 * 4./9.* C1P
        BN  = BN - 1./2./12.566 * 1./9.* C1P
       ENDIF
C                   
       RETURN
       END
C                                                                               
       SUBROUTINE FRHOIN (VA, XI, GL, XN)
       IMPLICIT DOUBLE COMPLEX (A - Z)
       XNV = XN + 0.5 - 1.0 
       C   = DCMPLX(1.703, 0.D0)
       VA  = (0.146 * (C/XNV - 1./(XNV+1.)) * 0.5)
       XNS = (XN - 0.3 - 1.0) 
       X1S = 3.0
       XI  = (0.188 * CBETA (XNS, X1S) * 0.5)
       XNG = XNS
       X1G = DCMPLX(2.5D0, 0.D0) 
       GL  = 0.225 * CBETA (XNG, X1G) * 0.5
       RETURN 
       END                                                       
C
C
       SUBROUTINE ANOMP (ANS, AM, AP, AL, BE, AB,       RPLUS, RQQ, RQG,
     1            RGQ, RGG, KNS, KSI, RKNS, RKSI, RKGL, C2Q, C1Q, C2G,
     2            C1G, C2P, C1P, ESQ, EFO, EFS,
     3            QQI, QGF, GQI, GGI, GGF,        NS1PI, NS1F, QQ1F,
     4            QG1F, GQ1I, GQ1F, GG1I, GG1F, QPF, QP1F, GP1F, C2QI,
     5            C1QI, C2GF, C1GF, C2PI, C1PI, XN, F)
       IMPLICIT DOUBLE COMPLEX (A - Z)
       DOUBLE PRECISION PI, B0, B02, B1, B10, Q, EQ(6), EQ2K, ESQ, 
     1                  EFO, EFS
       INTEGER F, K1, N
       COMMON / SCHEMEP / N
       DATA EQ / -1., 2., -1., 2., -1., 2./
       PI  = DACOS(-1.D0)
C...PARAMETERS RELATED TO THE NUMBER OF ACTIVE FLAVOURS F:
C....BETA-FUNCTION:
       B0  = 11. - 2./3. * F
       B02 = 2.* B0
       B1  = 102.- 38./3.* F
       B10 = B1 / B0
       Q   = 1./ (PI * B02)
C....AVERAGES OF QUARK CHARGES :
       ESQ = 0.0
       EFO = 0.0
       DO 1 K1 = 1, F
          EQ2K = EQ(K1) * EQ(K1) / 9.
          ESQ  = ESQ + EQ2K
          EFO  = EFO + EQ2K * EQ2K
  1    CONTINUE
       ESQ = ESQ / F
       EFO = EFO / F
       EFS = (EFO - ESQ * ESQ)
C...ANOMALOUS DIMENSIONS AND RELATED QUANTITIES IN LEADING ORDER,
C....HADRONIC PART :
       QQ = QQI
       QG = F * QGF
       GQ = GQI
       GG = GGI + F * GGF
       SQ = SQRT ((GG-QQ) * (GG-QQ) + 4.* QG * GQ)
       GP = 0.5 * (QQ + GG + SQ)
       GM = 0.5 * (QQ + GG - SQ)
       ANS = QQ / B02
       AM = GM / B02
       AP = GP / B02
       AL = (QQ-GP) / (GM-GP)
       BE = QG / (GM-GP)
       AB = GQ / (GM-GP)
C....ADDITIONAL TERMS FOR THE PHOTON STRUCTURE :
       KNS = EFS * QPF * F / 4. * Q
       KSI = ESQ * QPF * F / 4. * Q
C...NEXT TO LEADING ORDER : ANOMALOUS DIMENSIONS AND WILSON COEFFICIENTS
C...IN THE MS-BAR FACTORIZATION SCHEME
C....HADRONIC PART :
C      NS1M = NS1MI + F * NS1F
       NS1P = NS1PI + F * NS1F
       QQ1 = NS1P + F * QQ1F
       QG1 = QG1F * F
       GQ1 = GQ1I + F * GQ1F
       GG1 = GG1I + F * GG1F
       C1Q = C1QI
       C2Q = C2QI
       C1G = F * C1GF
       C2G = F * C2GF
C....COMBINATIONS OF ANOMALOUS DIMENSIONS FOR NLO EVOLUTION :
C      RMIN = (NS1M - QQ * B10) / B02
       RPLUS = (NS1P - QQ * B10) / B02
       RQQ = (QQ1 - QQ * B10) / B02
       RQG = (QG1 - QG * B10) / B02
       RGQ = (GQ1 - GQ * B10) / B02
       RGG = (GG1 - GG * B10) / B02
C...ADDITIONAL TERMS FOR THE PHOTON STRUCTURE :
       KNS1 = F * EFS * QP1F / 8.* Q
       KSI1 = F * ESQ * QP1F / 8.* Q
       KGL1 = F * ESQ * GP1F / 8.* Q
       C1P  = C1PI
       C2P  = C2PI
C....CHANGE OF THE FACTORIZATION SCHEME :
          DEL = 0.0
       IF ( N .EQ. 1 ) THEN
          DEL =  -C2PI
       ELSE IF ( N .EQ. 2 ) THEN
          DEL  = -C1PI
       ELSE IF ( N .EQ. 3 ) THEN
          DEL  = -C2PI - 24./ (XN - 1.) + 24./ XN
       END IF
C      C1P  = C1PI + DEL
C........(SEE RENO FOR AN EXPLANATION)
       C2P  = C2PI + DEL
       KNS1 = KNS1 - (QQ * DEL * F * EFS) / 8.* Q
       KSI1 = KSI1 - (QQ * DEL * F * ESQ) / 8.* Q
       KGL1 = KGL1 - (GQ * DEL * F * ESQ) / 8.* Q
C....COMBINATIONS FOR NLO EVOLUTION :
       RKNS = (KNS1 - KNS * B1 / B02) * 2.
       RKSI = (KSI1 - KSI * B1 / B02) * 2.
       RKGL = KGL1 * 2.
       RETURN
       END
C
C...CALCULATION OF ANOMALOUS DIMENSIONS (AND WILSON COEFFICIENTS)
C...UP TO THEIR DEPENDENCE OF THE NUMBER OF ACTIVE FLAVOURS F :
       SUBROUTINE FRCALCP (QQI, QGF, GQI, GGI, GGF,        NS1PI, NS1F,
     1               QQ1F, QG1F, GQ1I, GQ1F, GG1I, GG1F, QPF, QP1F,
     2               GP1F, C2QI, C1QI, C2GF, C1GF, C2PI, C1PI, N)
       IMPLICIT DOUBLE COMPLEX (A - Z)
       DOUBLE PRECISION ZETA2, ZETA3
       NS = N * N
       N1 = N + 1.
       N2 = N + 2.
       NM = N - 1.
C...LEADING ORDER :
       S1  = PSIFN (N1) + 0.577216
       QQI = (8./3.) * (-3.- 2./(N * N1) + 4.* S1)
       GQI = -2.* (NS + N +2.) / (N * N1 * N2)
       QGF = -(32./3.) * (NS + N + 2.) / (N * N1 * NM)
       GGI = -22.- 24./(N * NM) - 24./(N1 * N2) + 24.* S1
       GGF = 4./3.
C...NEXT TO LEADING ORDER :
       NT  = NS * N
       NFO = NT * N
       N1S = N1 * N1
       N1T = N1S * N1
C...ANALYTIC CONTINUATIONS OF N-SUMS AS GIVEN IN GLUECK ET AL. (1990) :
       ZETA2 = 1.644934
       ZETA3 = 1.202057
       S11 = S1 + 1./N1
       S2  = ZETA2 - PSIFN1 (N1)
       S21 = S2 + 1./N1S
       G1  = 0.5 * (PSIFN1 (N1/2.) - PSIFN1 (N/2.))
       G11 = 0.5 * (PSIFN1 (N2/2.) - PSIFN1 (N1/2.))
       SPMOM = 1.01/N1 - 0.846/N2 + 1.155/(N+3.) - 1.074/(N+4.) +
     1         0.55/(N+5.)
       SLC   = -5./8.* ZETA3
       SLV   = - ZETA2/2.* (PSIFN (N1/2.) - PSIFN (N/2.))
     1         + S1/NS + SPMOM
       SSCHLP = SLC + SLV
       SSTR2P = ZETA2 - PSIFN1 (N2/2.)
       SSTR3P = 0.5 * PSIFN2 (N2/2.) + ZETA3
C
C.......MINUS COMBINATIONS ARE NOT NEEDED HERE (CF. PHOTON STRUCTURE)
       NS1PA = 16.* S1 * (2.* N + 1.) / (NS * N1S) +
     1         16.* (2.* S1 - 1./(N * N1)) * ( S2 - SSTR2P ) +
     2         64.* SSCHLP + 24.* S2 - 3. - 8.* SSTR3P -
     3         8.* (3.* NT + NS -1.) / (NT * N1T) -
     4         16.* (2.* NS + 2.* N +1.) / (NT * N1T)
       NS1B  = S1 * (536./9. + 8.* (2.* N + 1.) / (NS * N1S)) -
     1         (16.* S1 + 52./3.- 8./(N * N1)) * S2 - 43./6. -
     2         (151.* NFO + 263.* NT + 97.* NS + 3.* N + 9.) *
     3         4./ (9.* NT * N1T)
       NS1C  = -160./9.* S1 + 32./3.* S2 + 4./3. +
     1         16.* (11.* NS + 5.* N - 3.) / (9.* NS * N1S)
       DNS1I = - 8.* S1 * S2 + 6.* S2 + 4.* (S2 - ZETA2) / (N * N1)
     1         + 8.* ZETA2 * S1 + 4.* S1 * (2.* N + 1.) / (NS * N1S)
     2         - (6.* NT + 9.* NS + 7.* N + 2.) / (NT * N1T)
     3         - 6.* ZETA2
*      NS1MI = -2./9.* NS1MA + 4.* NS1B
       NS1PI = -2./9.* NS1PA + 4.* NS1B - 128./9.* DNS1I
       NS1F  = 2./3. * NS1C
C...SPACELIKE SINGLET PIECES AS GIVEN IN FLORATOS ET AL. (1981),
C...DIFFERENCE OF THE TIMELIKE PIECES AS INVERTED FROM FURMANSKI AND
C...PETRONZIO (1980) BY GLUECK AND REYA (1992) :
       NFI = NFO * N
       NSI = NFI * N
       NSE = NSI * N
       NE  = NSE * N
       NN  = NE * N
       NMS = NM * NM
       NMT = NMS * NM
       N2S = N2 * N2
       N2T = N2S * N2
       QQ1FS = (5.* NFI + 32.* NFO + 49.* NT + 38.* NS + 28.* N
     1           + 8.) / (NM * NT * N1T * N2S) * (-32./3.)
       DQQ1F = - 80./9./ NM - 12./N + 12./NS + 8./NT - 4./N1
C........ERROR IN THE REPRINT (NOT IN THE PREPRINT) OF FP :
C...........-8./NS (WRONG) INSTEAD OF 12./NS (CORRECT).
     1         + 28./N1S + 8./N1T + 224./9./ N2 + 32./3./ N2S
       QQ1F  = QQ1FS - 16./3.* DQQ1F
       GQ1A = 8./3.* S11 * (NS + N + 2.) / (N * N1 * N2)
     1        + 8./3. * (- 5./(3.* N) + 1./NS - 1./(N * N1)
     2          + 4./(3.* N1) - 2./N1S - 4./(3.* N2) + 4./N2S)
       GQ1B = - (2.* S11*S11 - 2.* S11 - 10.* S21) * (NS + N + 2.)
     1          / (N * N1 * N2)
     2        + 4.* S11 * (1./N - 1./NS + 1./(N * N1) + 2./N1S
     3          - 4./N2S)
     4        - 12./N + 5./NS - 6./(N * N1) - 12./(N * N1S)
     5        + 4./(NS * N1) - 2./NT + 23./N1 - 4./N1S + 4./N1T
     6        - 20./N2
       GQ1C  = (2.* S11 * S11 - 10./3.* S11 - 6.* S21 + 2.* G11
     1           -6.* ZETA2) * (NS + N + 2.) / (N * N1 * N2)
     2         - 4.* S11 * (1./N - 2./NS + 1./(N * N1) + 4./N1S
     3           - 6./N2S)
     4         - 40./9./ NM + 26./9./ N + 8./3./ NS - 190./9./ N1
     5         + 22./3./(N * N1) + 68./3./ N1S + 16./N1T + 4./NT
     6         + 8./(N1S * N2) + 356/9./ N2 - 4./N2S
     7         - 8./(NS * N1S)
C...GQ1 OF FP HAS TO BE DIVIDED BY  2 * F
       GQ1F = - 2.* GQ1A / 2.
       GQ1I = (- 16./3. * GQ1B - 12.* GQ1C) / 2.
       QG1A  = (S1 * S1 - 3.* S2 - 4.* ZETA2) * (2./NM - 2./N + 1./N1)
     1         + 2.* S1 * (-2./(NM * N) -1./N1 + 3./N1S + 4./NMS
     2           - 4./NS)
     3         + 8./(N * NM) - 1./2./ N + 9./2./ N1 - 5./2./ N1S
     4         + 1./N1T - 8./(NMS * N) + 2./NT
       QG1B  = - (S1 * S1 - 5.* S2 + G1 - ZETA2) * (NS + N + 2.)
     1           / (NM * N * N1)
     2         + 2.* S1 * (NFI - NFO + NT - 9.* NS - 2.* N + 2.)
     3           / (NMS * NS * N1S)
     4         - 8./NMT + 6./NMS + 4./(NMS * N) + 17./(9.* NM)
     5         - 12./(NM * NS) + 5./N - 8./NS - 2./(NS * N1)
     6         - 1./N1 - 7./N1S - 2./N1T - 44./9./ N2 - 8./3./ N2S
       QG1F  = (- 128./9.* QG1A - 32.* QG1B) * 2.
C...QG1 OF FP HAS TO BE MULTIPLIED BY 2 * F
       GG1AS = 16./9.* (38.* NFO + 76.* NT + 94.* NS + 56.* N + 12.)
     1           / (NM * NS * N1S * N2)   -   160./9.* S1 + 32./3.
       DGG1A = S2 - 1./NMS + 1./NS - 1./N1S + 1./N2S - ZETA2
       GG1A  = GG1AS + 64./3.* DGG1A
       GG1BS = (2.* NSI + 4.* NFI + NFO - 10.* NT - 5.* NS - 4.* N
     1          - 4.) * 16. / (NM * NT * N1T * N2)   +   8.
       DGG1B = 80./9./ NM - 16./3./ NMS + 12./N - 16./NS + 8./NT
     1         + 4./N1 - 24./N1S + 8./N1T - 224./9./ N2
     2         - 16./3./ N2S
       GG1B  = GG1BS - 8.* DGG1B
       GG1CS = (2.* NFI + 5.* NFO + 8.* NT + 7.* NS - 2.* N - 2.)
     1           * 64.* S1 / (NMS * NS * N1S * N2S)
     2         + 536./9.* S1 - 64./3.
     3         + 32.* SSTR2P * (NS + N + 1.) / (NM * N * N1 * N2)
     4         - 16.* S1 * SSTR2P + 32.* SSCHLP - 4.* SSTR3P
     5         - 4.* (457.* NN + 2742.* NE + 6040.* NSE + 6098.* NSI
     6          + 1567.* NFI - 2344.* NFO - 1632.* NT + 560.* NS
     7          + 1488.* N + 576.) / (9.* NMS * NT * N1T * N2T)
       DGG1C = - S1 * S2 + S1 * (ZETA2 + 1./NMS - 1./NS + 1./N1S
     1         - 1./N2S)
     2         + (S2 - ZETA2) * (11/12.+ 1./NM - 1./N + 1./N1 - 1./N2)
     3         - 1./NMT + 11./12./ NMS - 1./(NMS * N)
     4         - 1./(NM * NS) - 7./12./ NS - 1./NT + 7./12./ N1S
     5         - 1./ N1T - 1./(N1S * N2) - 1./(N1 * N2S)
     6         - 11./12./ N2S - 1./N2T
       GG1C  = GG1CS - 64.* DGG1C
       GG1I  = 9.* GG1C
       GG1F  = 3./2.* GG1A  + 2./3.* GG1B
C...PHOTON PIECES :
       QPF  = -3./4.* QGF
       QP1F = 64./3.* QG1A
       GP1F = 4.* (-4./N + 12./N1 - 164./9./ N2 + 92./9./ NM - 10./NS
     1          -14./N1S - 16./3./ N2S - 16./3./ NMS + 4./NT + 4./N1T)
C...TIMELIKE WILSON COEFFICIENTS :
       C2QI = 4./3.* (2.* S1 * S1 - 2.* S2 + 3.* S1 - 9.
     1       - 2.* S1 / (N * N1) + 3./ N + 4./ N1 + 2./ NS)
       DC2Q = 8./3.* (6.* S2 + 3./N1S - 3./NS - 7./2./ N - 7./2./ N1)
       C2QI = C2QI + DC2Q
       C1QI = C2QI + 16./3./ N
       C2GF = (-S1 * (NS + N + 2.) / (NM * N * N1) - 4./(NM * N)
     1         + 4.* (1.- 2.* N) / (NMS * NS) - 3./ N1S) * 16./3.
       C1GF = C2GF + 16./3.* (4./NM - 4./N)
       C2PI = (- S1 * (NS + N + 2.) / (NM * N * N1) - 4./(NM * N)
     1         + 4.* (1.- 2.* N) / (NMS * NS) - 3./ N1S) * 4.
       C1PI = C2PI + 4.* (4./NM - 4./N)
       RETURN
       END
c$$$C
c$$$C...PSI - FUNCTION FOR COMPLEX ARGUMENT
c$$$       DOUBLE COMPLEX FUNCTION PSIFN (Z)
c$$$       DOUBLE COMPLEX Z, ZZ, RZ, DZ, SUB
c$$$       SUB = DCMPLX (0.D0,0.D0)
c$$$       ZZ = Z
c$$$  1    CONTINUE
c$$$       IF (DREAL (ZZ) .LT. 10.) THEN
c$$$         SUB = SUB - 1./ ZZ
c$$$         ZZ = ZZ + 1.
c$$$         GOTO 1
c$$$       END IF
c$$$       RZ = 1./ ZZ
c$$$       DZ = RZ * RZ
c$$$       PSIFN = SUB + LOG(ZZ) - RZ/2.- DZ/2520. * ( 210.+ DZ * (-21.+
c$$$     1         10.*DZ ))
c$$$       RETURN
c$$$       END
c$$$C
c$$$C...FIRST DERIVATIVE OF THE PSI - FUNCTION FOR COMPLEX ARGUMENT :
c$$$       DOUBLE COMPLEX FUNCTION PSIFN1 (Z)
c$$$       DOUBLE COMPLEX Z, ZZ, RZ, DZ, SUB
c$$$       SUB = DCMPLX (0.D0,0.D0)
c$$$       ZZ = Z
c$$$  1    CONTINUE
c$$$       IF (DREAL (ZZ) .LT. 10.) THEN
c$$$         SUB = SUB + 1./ (ZZ * ZZ)
c$$$         ZZ = ZZ + 1.
c$$$         GOTO 1
c$$$       END IF
c$$$       RZ = 1./ ZZ
c$$$       DZ = RZ * RZ
c$$$       PSIFN1 = SUB + RZ + DZ/2. * ( 1 + RZ/630. * ( 210.- DZ * ( 42.-
c$$$     1         DZ * ( 30.- 42.*DZ ))))
c$$$       RETURN
c$$$       END
c$$$C
c$$$C...SECOND DERIVATIVE OF THE PSI - FUNCTION FOR COMPLEX ARGUMENT :
c$$$       DOUBLE COMPLEX FUNCTION PSIFN2 (Z)
c$$$       DOUBLE COMPLEX Z, ZZ, RZ, DZ, SUB
c$$$       SUB = DCMPLX (0.D0,0.D0)
c$$$       ZZ = Z
c$$$  1    CONTINUE
c$$$       IF (DREAL (ZZ) .LT. 10.) THEN
c$$$         SUB = SUB - 2./ (ZZ * ZZ * ZZ)
c$$$         ZZ = ZZ + 1.
c$$$         GOTO 1
c$$$       END IF
c$$$       RZ = 1./ ZZ
c$$$       DZ = RZ * RZ
c$$$       PSIFN2 = SUB - DZ/60. * ( 60.+ RZ * ( 60.+ RZ * ( 30.- DZ *
c$$$     1         ( 10.- DZ * ( 10.- DZ * ( 18.- 50.* DZ ))))))
c$$$       RETURN
c$$$       END
c$$$C
c$$$C...BETA FUNCTION FOR COMPLEX ARGUMENT :
c$$$       DOUBLE COMPLEX FUNCTION CBETA (Z1, Z2)
c$$$       IMPLICIT DOUBLE COMPLEX (A - Z)
c$$$       LNGAM (X) = (X - 0.5) * LOG (X) - X + 0.91893853 + 1./(12.* X)
c$$$     1              * (1.- 1./(30.* X*X) * (1.- 1./(3.5 * X*X)
c$$$     2              * (1.- 4./(3.* X*X))))
c$$$       SUB = DCMPLX (0.D0, 0.D0)
c$$$       ZZ1 = Z1
c$$$  1    CONTINUE
c$$$       IF ( DREAL (ZZ1) .LT. 15.) THEN
c$$$          SUB = SUB + LOG ((ZZ1+Z2) / ZZ1)
c$$$          ZZ1 = ZZ1 + 1.
c$$$          GOTO 1
c$$$       END IF
c$$$       ZZ2 = Z2
c$$$  2    CONTINUE
c$$$       IF ( DREAL (ZZ2) .LT. 15.) THEN
c$$$          SUB = SUB + LOG ((ZZ1+ZZ2) / ZZ2)
c$$$          ZZ2 = ZZ2 + 1.
c$$$          GOTO 2
c$$$       END IF
c$$$       LG1 = LNGAM (ZZ1)
c$$$       LG2 = LNGAM (ZZ2)
c$$$       LG12 = LNGAM (ZZ1 + ZZ2)
c$$$       CBETA = EXP (LG1 + LG2 - LG12 + SUB)
c$$$       RETURN
c$$$       END
