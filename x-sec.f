CAM NEW SIGMAS INCLUDED EQS A.1-10 IN FLORIAN 7 WV PRD71,114004
c https://journals.aps.org/prd/pdf/10.1103/PhysRevD.71.114004
       SUBROUTINE DGDIGLO(XN,SIGDGQ,SIGDGG,SIGDQQ1,SIGDQQ2,SIGDQQ3,
     1              SIGDQQ4,SIGDQQ5,SIGDQQ6,SIGDQG,SIGDGG1,SIGDGG2)
       IMPLICIT DOUBLE COMPLEX (A-Z)
       DOUBLE PRECISION PI,CF,CA

       external cbeta

C
       PI=DACOS(-1.D0)
       CF=4.D0/3.D0
       CA=3.D0
       XNH=DCMPLX(0.5D0,0.D0)
C mbeekveld these functions are the functions from catani 9903436 and 0501258 with the phase space integration measure.
       SIGDGQ = PI*CF/CA * CBETA(XN+1.D0,XNH) 
     x      * (XN+2.D0)/(XN+3.D0/2.D0)
       SIGDGG = PI/8.D0/CA * CBETA(XN+1.D0,XNH) 
     x      * (5.D0*XN+7.D0)/(XN+3.D0/2.D0)
       SIGDQQ1= PI*CF/3.D0/CA * CBETA(XN,5.D0*XNH) 
     x      * (5.D0*XN**2 +15.D0*XN+12.D0)
       SIGDQQ2= PI*CF/3.D0/CA * CBETA(XN,5.D0*XNH) 
     x      * (5.D0*XN**2 +15.D0*XN+12.D0)
       SIGDQQ3= PI*CF/6.D0/CA * CBETA(XN+1.D0,5.D0*XNH) 
     x      * (XN+1.D0)*(XN+3.D0)
       SIGDQQ4= 2.D0*PI*CF/3.D0/CA**2 * CBETA(XN,5.D0*XNH) 
     1      * (CA*(5.D0*XN**2 +15.D0*XN+12.D0)
     2      -2.D0*XN*(3.D0+2.D0*XN))
       SIGDQQ5= PI*CF/15.D0/CA**2 * CBETA(XN,7.D0*XNH) 
     1      * (CA*(11.D0*XN**3 +59.D0*XN**2 +102.D0*XN+60.D0)
     2      +XN*(3.D0+XN)*(5.D0+2.D0*XN) )

CCC   Omitted: g-g final state
       SIGDQQ6= PI*CF/3.D0/CA * CBETA(XN+1.D0,5.D0*XNH) *
     1      (2.D0*CF*(XN+2.D0)*(5.D0+2.D0*XN)-CA*(XN+1.D0)*(XN+3.D0))
CPM    erronymous     1      (2.*CF*(XN+2.)*(5.+2.*XN)-CA*(N+1)*(N+3))

       SIGDQG= PI/6.D0/CA * CBETA(XN,5.D0*XNH) *
     1      (CF*XN*(7.D0+5.D0*XN)+2.D0*CA*(5.D0*XN**2 +15.D0*XN+12.D0))

CCC   Omitted: g-g final state
       SIGDGG1= PI*CA/5.D0/CF * CBETA(XN,7.D0*XNH) *
     1      (9.D0*XN**3 +45.D0*XN**2 +72.D0*XN+40.D0)

       SIGDGG2= PI/12.D0/CA/CF * CBETA(XN+1.D0,5.D0*XNH) *
     1      (2.D0*CF*(XN+2.D0)*(5.D0+2.D0*XN)-CA*(XN+1.D0)*(XN+3.D0))
       
C

c$$$       write(*,*)"(dgdiglo)  XN       = ", XN
c$$$       write(*,*)"(dgdiglo)  SIGDGQ   = ", SIGDGQ
c$$$       write(*,*)"(dgdiglo)  SIGDGG   = ", SIGDGG
c$$$       write(*,*)"(dgdiglo)  SIGDQQ1  = ", SIGDQQ1
c$$$       write(*,*)"(dgdiglo)  SIGDQQ2  = ", SIGDQQ2
c$$$       write(*,*)"(dgdiglo)  SIGDQQ3  = ", SIGDQQ3
c$$$       write(*,*)"(dgdiglo)  SIGDQQ4  = ", SIGDQQ4
c$$$       write(*,*)"(dgdiglo)  SIGDQQ5  = ", SIGDQQ5
c$$$       write(*,*)"(dgdiglo)  SIGDQG   = ", SIGDQG
       RETURN
       END
