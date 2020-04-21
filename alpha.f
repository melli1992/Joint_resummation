      FUNCTION ALPHAS (Q2)
*********************************************************************
*                                                                   *
*   THE ALPHA_S ROUTINE.                                            *
*                                                                   *
*   INPUT :  Q2    =  scale in GeV**2  (not too low, of course);    *
*            NAORD =  1 (LO),  2 (NLO).                             *
*                                                                   *
*   OUTPUT:  alphas_s/(4 pi) for use with the GRV(98) partons.      *  
*                                                                   *
*******************************************************i*************
*
      IMPLICIT DOUBLE PRECISION (A - Z)
      INTEGER NF, K, I, NAORD
      DIMENSION LAMBDAN (3:6), Q2THR (3)
c      COMMON / SCALES / QQ2, Q2START, Q2S, Q20, Q21, Q2MUR, 
c     1                  Q2MUF, Q2THR, Q20F, Q21F, LAMBDAN
      COMMON / SCALES / QQ2, Q2START, Q20, Q21, Q2MUR, 
     1                  Q2MUF, Q2THR, LAMBDAN, Q2S, Q20F, Q21F
*
      NAORD=2
*
*...DETERMINATION OF THE APPROPRIATE NUMBER OF FLAVOURS :
c$$$       write(*,*)"(alp) Q2      = ", q2
c$$$       write(*,*)"(alp) Q2START = ", q2start
c$$$       write(*,*)"(alp) Q2THR   = ", q2thr
      NF = 3
      DO 10 K = 1, 3
      IF (Q2 .GT. Q2THR (K)) THEN
         NF = NF + 1
      ELSE
          GO TO 20
       END IF
  10   CONTINUE
*
*...LO ALPHA_S AND BETA FUNCTION FOR NLO CALCULATION :
  20   B0 = 11.- 2./3.* NF
       B1 = 102.- 38./3.* NF
       B10 = B1 / (B0*B0)
       IF (NAORD .EQ. 1) THEN
*         LAM2 = LAMBDAL (NF) * LAMBDAL (NF)
*         ALP  = 1./(B0 * DLOG (Q2/LAM2))
*         GO TO 1
       ELSE IF (NAORD .EQ. 2) then
         LAM2 = LAMBDAN (NF) * LAMBDAN (NF)
         B1 = 102.- 38./3.* NF
         B10 = B1 / (B0*B0)
       ELSE
         WRITE (6,91)
  91     FORMAT ('INVALID CHOICE FOR ORDER IN ALPHA_S')
         STOP
       END IF
c$$$       write(*,*)"(alp) B0      = ", b0
c$$$       write(*,*)"(alp) B1      = ", b1
c$$$       write(*,*)"(alp) B10     = ", b10
c$$$       write(*,*)"(alp) NF      = ", nf
c$$$       write(*,*)"------------------------"
*
*...START VALUE FOR NLO ITERATION :
       LQ2 = DLOG (Q2 / LAM2)
       ALP = 1./(B0*LQ2) * (1.- B10*DLOG(LQ2)/LQ2)
*
*...EXACT NLO VALUE, FOUND VIA NEWTON PROCEDURE :
       DO 2 I = 1, 6
       XL  = DLOG (1./(B0*ALP) + B10)
       XLP = DLOG (1./(B0*ALP*1.01) + B10)
       XLM = DLOG (1./(B0*ALP*0.99) + B10)
       Y  = LQ2 - 1./ (B0*ALP) + B10 * XL
       Y1 = (- 1./ (B0*ALP*1.01) + B10 * XLP
     1       + 1./ (B0*ALP*0.99) - B10 * XLP) / (0.02D0*ALP)
       ALP = ALP - Y/Y1
  2    CONTINUE
*
*...OUTPUT :
  1    ALPHAS = ALP
       RETURN
       END

CCCCCCCCCCCCCCCCCCCCCC
CCC
CCC   Complex alpha
CCC
CCCCCCCCCCCCCCCCCCCCCC

      FUNCTION CALPHAS (CQ2,NAORD,NF)
*********************************************************************
*                                                                   *
*   THE ALPHA_S ROUTINE FOR COMPLEX VARIABLES                               
*                                                                   *
*   INPUT :  CQ2    = complex scale in GeV**2                        *
*            NAORD =  1 (LO)[inactive],  2 (NLO).                             *
*            NF = number of active flavors
*                                                                   *
*   OUTPUT:  alphas_s/(4 pi) 
*                                                                   *
*******************************************************i*************
*
      IMPLICIT DOUBLE COMPLEX (A - Z)
      INTEGER NF, K, I, NAORD
      DOUBLE PRECISION QQ2, Q2START, Q20, Q21, Q2MUR, 
     1     Q2MUF, Q2THR(3), LAMBDAN(3:6), q2s, q20f, q21f
c      COMMON / SCALES / QQ2, Q2START, Q2S, Q20, Q21, Q2MUR, 
c     1                  Q2MUF, Q2THR, Q20F,Q21F, LAMBDAN
      COMMON / SCALES / QQ2, Q2START, Q20, Q21, Q2MUR, 
     1                  Q2MUF, Q2THR, LAMBDAN, Q2S, Q20F, Q21F
*
*
*...LO ALPHA_S AND BETA FUNCTION FOR NLO CALCULATION :
  20   B0 = 11.- 2./3.* NF
       B1 = 102.- 38./3.* NF
       B10 = B1 / (B0*B0)
       IF (NAORD .EQ. 1) THEN
*         LAM2 = LAMBDAL (NF) * LAMBDAL (NF)
*         ALP  = 1./(B0 * DLOG (Q2/LAM2))
*         GO TO 1
          print*,"CALPHAS: NAORD=1 IS NOT IMPLEMENTED; STOP"
          STOP
       ELSE IF (NAORD .EQ. 2) then
         LAM2 = LAMBDAN (NF) * LAMBDAN (NF)
         B1 = 102.- 38./3.* NF
         B10 = B1 / (B0*B0)
       ELSE
         WRITE (6,91)
  91     FORMAT ('INVALID CHOICE FOR ORDER IN CALPHA_S')
         STOP
       END IF

       CLQ2 = LOG (CQ2 / LAM2)
       CALP = 1./(B0*CLQ2) * (1.- B10*LOG(CLQ2)/CLQ2)
*
*...OUTPUT :
  1    CALPHAS = CALP
       RETURN
       END

