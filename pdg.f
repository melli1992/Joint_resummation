CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC
CCC   Re-combining the fluxes taking fragmentation 
CCC   components into account (based on the previously
CCC   used subroutine PDG (see below) )
CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CAM.....PDGNQ1,...PDGNQ6  ADDED BELOW IN PDG
       SUBROUTINE PDGMOD(UVN,DVN,SVN,CVN,BVN,USN,DSN,SSN,CSN,
     x     BSN,GLN,
     x     UN,DN,SN,GLNF,
     x     PDGNQ,PDGNG,
     x     PDGNQQ1,PDGNQQ2,PDGNQQ3,PDGNQQ4,PDGNQQ5,PDGNQQ6,
     x     PDGNGG1,PDGNGG2,PDGNGMOD,
     x     FLAG)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     IN : UVN,DVN,USN,DSN,SSN,GLN,
C          UN,DN,SN,CN,GLNF
C     OUT: PDGNQ,PDGNG,
C          PDGNQQ1,PDGNQQ2,PDGNQQ3,PDGNQQ4,PDGNQQ5,PDGNQQ6,
C          PDGNGG1,PGNGG2,PDGNGMOD
C     
C     FLAG = 1 : PP - collisions
C          = 2 : PPAR - collisions
C          = 3 : PN - collision NOT IMPLEMENTED YET!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT DOUBLE COMPLEX (A - Z)
       
       INTEGER FLAG

       DOUBLE PRECISION EQUP, EQDO, ALPS, ALP0, ALP1, ALPC, 
     1                  ALPB, ALPT, ALPQ, ALPQR, LMQ, LMQR

       COMMON / COUPL / ALPS, ALP0, ALP1, ALPC, ALPB, ALPT, ALPQ, 
     1                  ALPQR, LMQ, LMQR
C
       CSN = DCMPLX(0.D0,0.D0)
       CVN  = DCMPLX(0.D0,0.D0)
	   BSN = DCMPLX(0.D0,0.D0)
       BVN  = DCMPLX(0.D0,0.D0)
C
       EQUP = 4.D0/9.D0
       EQDO = 1.D0/9.D0

c mbeekveld SVN added
	
       IF(FLAG .eq. 1) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ***  P-P :
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     fluxes for the direct (non-fragmenting) contribution
       F1COMBN =       EQUP*(UVN+2.D0*USN+CVN+2.D0*CSN) + 
     1                 EQDO*(DVN+2.D0*DSN+SVN+2.D0*SSN
     2                 +BVN+2.D0*BSN)
       PDGNQ = 2.D0 * (EQUP * (UVN+USN)*USN + 
     1                 EQUP * (CVN+CSN)*CSN + 
     1                 EQDO * (DVN+DSN)*DSN +
     2                 EQDO * (SVN+SSN)*SSN +
     2                 EQDO * (BVN+BSN)*BSN )
       PDGNG = 2.D0 *  F1COMBN * GLN 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CPM   FLUX RE-DEFINITION IN ORDER TO INCLUDE FF'S
CPM   construct fluxes incl. fragmentation
CPM   PP-collisions
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
ccc   flux for SIGDQQ1: q qp --> q qp
ccc   treating also \bar{q} \bar{qp} --> \bar{q} \bar{qp}
ccc   taking into account different multiplicity of anti-quarks
       PDGNQQ1 = 2.D0  * (
     x        (UVN+USN)*(DVN+DSN)
     x      + (UVN+USN)*(SVN+SSN) 
     x      + (UVN+USN)*(BVN+BSN)
     x      + (UVN+USN)*(CVN+CSN)   
     x      +      USN * DSN
     x      +      USN * SSN   
     x      +      USN * CSN
     x      +      USN * BSN ) * UN * EQUP 
     x      +    2.D0  * (
     x        (DVN+DSN)*(SVN+SSN) 
     x      + (DVN+DSN)*(UVN+USN)
     x      + (DVN+DSN)*(BVN+BSN) 
     x      + (DVN+DSN)*(CVN+CSN)
     x      +      DSN * SSN  
     x      +      DSN * USN
     x      +      DSN * BSN  
     x      +      DSN * CSN ) * DN * EQDO
     x      +    2.D0  * (
     x        (SVN+SSN) * (UVN+USN)
     x      + (SVN+SSN) * (DVN+DSN)
     x      + (SVN+SSN) * (CVN+CSN)
     x      + (SVN+SSN) * (BVN+BSN)
     x      +      SSN * USN 
     x      +      SSN * DSN 
     x      +      SSN * BSN 
     x      +      SSN * CSN ) * SN * EQDO


ccc   flux for SIGDQQ2: q \bar{qp} --> q \bar{qp}
       PDGNQQ2 = 2.D0 *( (
     x             (UVN+USN)* DSN 
     x      +      (UVN+USN)* SSN 
     x      +      (UVN+USN)* CSN 
     x      +      (UVN+USN)* BSN 
     x      +      (DVN+DSN)* USN 
     x      +      (SVN+SSN)* USN
     x      +      (BVN+BSN)* USN 
     x      +      (CVN+CSN)* USN ) * UN * EQUP
     x      +            (
     x            (DVN+DSN)* SSN 
     x      +      (DVN+DSN)* USN 
     x      +      (DVN+DSN)* BSN 
     x      +      (DVN+DSN)* CSN 
     x      +      (UVN+USN)* DSN 
     x      +      (SVN+SSN) * DSN
     x      +      (CVN+CSN)* DSN 
     x      +      (BVN+BSN) * DSN ) * DN * EQDO
     x      +            (
     x            (SVN+SSN) * USN
     x      +      (SVN+SSN) * DSN
     x      +      (SVN+SSN) * CSN
     x      +      (SVN+SSN) * BSN
     x      +      (UVN+USN)* SSN
     x      +      (DVN+DSN)* SSN 
     x      +      (CVN+CSN)* SSN
     x      +      (BVN+BSN)* SSN  )* SN * EQDO 
     x      )

ccc   flux for SIGDQQ3: q \bar{q} --> qp \bar{qp}
       PDGNQQ3 =   2.D0 *(
     x        (UVN+USN)*USN * 2.D0 * ((DN + SN) * EQDO) 
     x      )
     x      +      2.D0 *(
     x        (DVN+DSN)*DSN * 2.D0 *  UN * EQUP 
     x      +   (DVN+DSN)*DSN * 2.D0 *  SN * EQDO  
     x      )
     x      +      2.D0 *(
     x      (SVN+SSN)*SSN  * 2.D0 *  DN * EQDO
     x      +  (SVN+SSN)*SSN * 2.D0 *  UN * EQUP
     x      )
     x      +      2.D0 *(
     x      (CVN+CSN)*CSN  * 2.D0 *  DN * EQDO
     x      +  (CVN+CSN)*CSN  * 2.D0 *  SN * EQDO
     x      +  (CVN+CSN)*CSN * 2.D0 *  UN * EQUP
     x      )
     x      +      2.D0 *(
     x      (BVN+BSN)*BSN  * 2.D0 *  DN * EQDO
     x      +  (CVN+BSN)*BSN  * 2.D0 *  SN * EQDO
     x      +  (CVN+BSN)*BSN * 2.D0 *  UN * EQUP
     x      )


ccc   flux for SIGDQQ4: q q --> q q
ccc   also taking \bar{q} \bar{q} --> \bar{q} \bar{q} 
ccc   into account
       PDGNQQ4 =  (
     x        (UVN+USN)*(UVN+USN) 
     x      +      USN * USN ) * 2.D0 * UN * EQUP
     x      +     ( 
     x        (DVN+DSN)*(DVN+DSN) 
     x      +      DSN * DSN ) * 2.D0 * DN * EQDO  
     x      +    (
     x         (SVN+SSN)*(SVN+SSN)
     x       +      SSN * SSN ) * 2.D0 * SN * EQDO
 

ccc   flux for SIGDQQ5: q \bar{q} --> q \bar{q}
       PDGNQQ5 = 2.D0 *( 
     x        (UVN+USN)*USN  * 2.D0 * UN * EQUP
     x      ) 
     x      +    2.D0 *( 
     x        (DVN+DSN)*DSN  * 2.D0 * DN * EQDO +
     x        (SVN+SSN)*SSN       * 2.D0 * SN * EQDO 
     x      )

ccc   flux for SIGDQQ6: q \bar{q} --> g g
ccc   NB! not complete (and not needed)

ccc   gluons - not needed for now ccccccccccccccccc
ccc   flux for SIGDGG1: g g --> g g
c$$$       PDGNGG1 = GLN * GLN * 2. * GLNF
ccccccccccccccccccccccccccccccccccccccccccccccccccc

ccc   flux for SIGDGG2: g g --> q \bar{q}
       PDGNGG2 =   GLN * GLN * 2.D0 
     x         * ((UN       * EQUP) 
     x         +  (DN + SN) * EQDO)

ccccccccccccccccccccccccccccccccccccccccccccccccccc


c mbeekveld was error here!
ccc   flux for SIGDQG: q g --> q g / g q
       PDGNGMOD = 2.D0 * ( 
     x         (
     x         (UVN + 2.D0*USN) * UN * EQUP  
     x      )
     x      +  (
     x         (DVN + 2.D0*DSN) * DN * EQDO  
     x      +  (SVN+ 2.D0*SSN)  * SN * EQDO
     x      ) 
     x      ) * GLN 


      ELSEIF(FLAG .EQ. 2) THEN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ***  P-PBAR :
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     fluxes for the direct (non-fragmenting) contribution
        F1COMBN = EQUP * (UVN+2.D0*USN) + 
     1            EQDO * (DVN+2.D0*DSN+SVN+2.D0*SSN)
        PDGNQ   = EQUP * ( (UVN+USN)*(UVN+USN) + USN*USN ) +
     1            EQDO * ( (DVN+DSN)*(DVN+DSN) + DSN*DSN   +
     2                     (SVN+SSN)* SSN +SSN*SSN)
        PDGNG   = 2.D0 * F1COMBN * GLN 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CPM   FLUX RE-DEFINITION IN ORDER TO INCLUDE FF'S
CPM   construct fluxes incl. fragmentation
CPM   PPbar-collisions
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
ccc   flux for SIGDQQ1: q qp --> q qp
ccc   also taking \bar{q} \bar{qp} --> \bar{q} \bar{qp}
ccc   into account
ccc   processes
       PDGNQQ1 = 
ccc   u/\bar{u}-quark coming from proton
     x            ((UVN+USN)* DSN 
     x        +    (UVN+USN)* SSN
     x        +    (UVN+USN)* CSN 
     x        +    (UVN+USN)* BSN
     x        +         USN *(DVN+DSN) 
     x        +         USN *(SSN+SVN) 
     x        +         USN *(BVN+BSN) 
     x        +         USN *(CSN+CVN) 
ccc   u/\bar{u}-quark coming from anti-proton
     x        +    (DVN+DSN)* USN
     x        +    (SVN+SSN)* USN
     x        +    (BVN+BSN)* USN
     x        +    (CVN+CSN)* USN
     x        +         DSN *(UVN+USN)
     x        +         SSN *(UVN+USN)
     x        +         CSN *(UVN+USN)
     x        +         BSN *(UVN+USN) )* UN * EQUP 
ccc   d/\bar{d}-quark coming from proton
     x        +   ((DVN+DSN)* USN
     x        +    (DVN+DSN)* SSN 
     x        +    (DVN+DSN)* BSN
     x        +    (DVN+DSN)* CSN 
     x        +         DSN *(UVN+USN) 
     x        +         DSN *(SSN+SVN) 
     x        +         DSN *(CVN+CSN) 
     x        +         DSN *(BSN+BVN) 
ccc   d/\bar{d}-quark coming from anti-proton
     x        +    (UVN+USN)* DSN
     x        +    (SVN+SSN)* DSN
     x        +    (BVN+BSN)* DSN
     x        +    (CVN+CSN)* DSN
     x        +         USN *(DVN+DSN)
     x        +         SSN *(DVN+DSN)
     x        +         BSN *(DVN+DSN)
     x        +         CSN *(DVN+DSN) )* DN * EQDO 
ccc   s/\bar{s}-quark coming from proton
     x        +  ( (SVN+SSN) * USN 
     x        +    (SVN+SSN) * DSN
     x        +    (SVN+SSN) * BSN 
     x        +    (SVN+SSN) * CSN
     x        +         SSN *(UVN+USN)
     x        +         SSN *(DVN+DSN)
     x        +         SSN *(BVN+BSN)
     x        +         SSN *(CVN+CSN) 
ccc   s/\bar{s}-quark coming from anti-proton
     x        +    (UVN+USN)* SSN
     x        +    (DVN+DSN)* SSN
     x        +    (BVN+BSN)* SSN
     x        +    (CVN+CSN)* SSN
     x        +         USN * (SVN+SSN)
     x        +         DSN * (SVN+SSN)
     x        +         CSN * (SVN+SSN)
     x        +         BSN * (SVN+SSN)  )* SN * EQDO


ccc   flux for SIGDQQ2: q \bar{qp} --> q \bar{qp}
ccc   Then factor of 2.D0 here is for the 
ccc   crossing in the IS
       PDGNQQ2 =      ( 
ccc   u/\bar{u}-quark coming from proton
     x       (UVN+USN)*(DVN+DSN)
     x      +(UVN+USN)*(SVN+SSN)
     x      +(UVN+USN)*(BVN+BSN)
     x      +(UVN+USN)*(CVN+CSN)
     x      +     USN * DSN    
     x      +     USN * SSN 
     x      +     USN * CSN    
     x      +     USN * BSN 
ccc   u/\bar{u}-quark coming from anti-proton
     x      +     DSN * USN
     x      +     SSN * USN
     x      +     BSN * USN
     x      +     CSN * USN
     x      +(DVN+DSN)*(UVN+USN)
     x      +(SVN+SSN)*(UVN+USN)
     x      +(BVN+BSN)*(UVN+USN)
     x      +(CVN+CSN)*(UVN+USN) ) * UN * EQUP
     x      +         (
ccc   d/\bar{d}-quark coming from proton       
     x      +(DVN+DSN)*(UVN+USN) 
     x      +(DVN+DSN)*(SVN+SSN)           
     x      +(DVN+DSN)*(BVN+BSN) 
     x      +(DVN+DSN)*(CVN+CSN)     
     x      +     DSN * USN   
     x      +     DSN * SSN    
     x      +     DSN * CSN   
     x      +     DSN * BSN 
ccc   d/\bar{d}-quark coming from anti-proton       
     x      +     USN * DSN
     x      +     SSN * DSN     
     x      +     BSN * DSN
     x      +     CSN * DSN
     x      +(UVN+USN)*(DVN+DSN)
     x      +(SVN+SSN)*(DVN+DSN)
     x      +(BVN+BSN)*(DVN+DSN)
     x      +(CVN+CSN)*(DVN+DSN) ) * DN * EQDO
     x      +         (
ccc   s/\bar{s}-quark coming from proton       
     x      +(SVN+SSN)*(UVN+USN)
     x      +(SVN+SSN)*(DVN+DSN)    
     x      +(SVN+SSN)*(BVN+BSN)
     x      +(SVN+SSN)*(CVN+CSN)
     x      +     SSN * USN
     x      +     SSN * DSN
     x      +     SSN * CSN
     x      +     SSN * BSN
ccc   s/\bar{s}-quark coming from anti-proton         
     x      +     USN * SSN      
     x      +     DSN * SSN           
     x      +     BSN * SSN      
     x      +     CSN * SSN   
     x      +(UVN+USN)*(SVN+SSN)
     x      +(DVN+DSN)*(SVN+SSN)
     x      +(CVN+CSN)*(SVN+SSN)
     x      +(BVN+BSN)*(SVN+SSN) ) * SN * EQDO  

ccc   flux for SIGDQQ3: q \bar{q} --> qp \bar{qp}
       PDGNQQ3 =  (
ccc   u from proton/\bar{u} from anti-proton
     x      (UVN+USN)*(UVN+USN) 
ccc   u from anti-proton/\bar{u} from proton
     x      +    USN * USN) * 2.D0 * (DN + SN) * EQDO 
     x      +     (
ccc   d from proton/\bar{d} from anti-proton
     x      (DVN+DSN)*(DVN+DSN) 
ccc   d from anti-proton/\bar{d} from proton
     x      +    DSN * DSN) * 2.D0 * (UN * EQUP + SN * EQDO) 
     x      +     (
ccc   s from proton/\bar{s} from anti-proton
     x        (SVN+SSN)*(SVN+SSN)
ccc   s from anti-proton/\bar{s} from proton
     x      + SSN* SSN) * 2.D0 * (UN * EQUP + DN * EQDO) 
     x      +     (
ccc   c from proton/\bar{c} from anti-proton
     x        (CVN+CSN)*(CVN+CSN)
ccc   c from anti-proton/\bar{C} from proton
     x      + CSN* CSN) * 2.D0 * (UN * EQUP + (DN+SN)*EQDO) 
     x      +     (
ccc   b from proton/\bar{b} from anti-proton
     x        (BVN+BSN)*(BVN+BSN)
ccc   b from anti-proton/\bar{b} from proton
     x      + BSN* BSN) * 2.D0 * (UN * EQUP + (DN+SN)*EQDO) 

ccc   flux for SIGDQQ4: q q --> q q
       PDGNQQ4 =      (
     x       (UVN+USN)* USN   
     x      +     USN *(UVN+USN) ) * 2.D0 * UN * EQUP 
     x      +         ( 
     x       (DVN+DSN)* DSN 
     x      +     DSN *(DVN+DSN) ) * 2.D0 * DN * EQDO  
     x      +   (
     x           (SVN+SSN)*(SVN+SSN) 
     x      +     SSN * SSN )      * 2.D0 * SN * EQDO 

ccc   flux for SIGDQQ5: q \bar{q} --> q \bar{q}
       PDGNQQ5 =       ( 
     x        (UVN+USN)*(UVN+USN) 
     x      +      USN * USN )     * 2.D0 * UN * EQUP
     x      +          ( 
     x        (DVN+DSN)*(DVN+DSN) 
     x      +      DSN * DSN )     * 2.D0 * DN * EQDO 
     x      +   (
     x           (SVN+SSN)*(SVN+SSN) 
     x      +       SSN * SSN )     * 2.D0 * SN * EQDO 

ccc   flux for SIGDQQ6: q \bar{q} --> g g
ccc   not complete (and not needed for now)

ccc   gluons - not needed for now ccccccccccccccccc
ccc   flux for SIGDGG1: g g --> g g
c$$$       PDGNGG1 = GLN * GLN * 2. * GLNF
ccccccccccccccccccccccccccccccccccccccccccccccccccc

ccc   flux for SIGDGG2: g g --> q \bar{q}
       PDGNGG2 = GLN * GLN * 2.D0 
     x         * (UN * EQUP + (DN + SN) * EQDO)

ccccccccccccccccccccccccccccccccccccccccccccccccccc

ccc   flux for SIGDQG: q g --> q g / g q
ccc   NB! not taking gluon fragmentation into account
       PDGNGMOD =  GLN* ( ( 
ccc   gluon coming from proton
     x        UVN + 2.D0*USN 
ccc   gluon coming from anti-proton
     x      + UVN + 2.D0*USN   ) * UN * EQUP 
     x      +           (
ccc   gluon coming from proton       
     x        DVN + 2.D0*DSN
ccc   gluon coming from anti-proton
     x      + DVN + 2.D0*DSN   )* DN * EQDO  
     x      +           (
ccc   gluon coming from proton       
     x       + SVN+ 2.D0*SSN  
ccc   gluon coming from anti-proton
     x      + SVN+ 2.D0*SSN        )* SN * EQDO)


      ELSEIF(FLAG .EQ. 3) THEN
	write(*,*) "NOT CORRECT" 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ***  P-N : p-Be(4,9) (E706 Experiment)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
ccc   original F1COMBN
       F1COMBN  = EQUP*(3.D0/2.D0*(UVN+2.D0*USN)+(DVN+2.D0*DSN)/2.D0
     x          + 4.D0*CSN)  
     x          + EQDO*(3.D0/2.D0*(DVN+2.D0*DSN)+(UVN+2.D0*USN)/2.D0
     x          + 4.D0*SSN)

ccc   Beryllium parton densities (auxillariy) 
       UVNBER   = (UVN + DVN)/2.D0
       USNBER   = (USN + DSN)/2.D0
       UBER     = UVNBER + USNBER
       UBARBER  = USNBER
       
       DVNBER   = (DVN + UVN)/2.D0
       DSNBER   = (DSN + USN)/2.D0
       DBER     = DVNBER + DSNBER
       DBARBER  = DSNBER

ccc   QQ initial states
       PDGNQ    = EQUP * ( 
ccc   u (proton) x \bar{u} (Be(4,8))
     x        (UVN+USN)* UBARBER  
ccc   \bar{u} (proton) x u (Be(4,8))
     x      +  USN * UBER )  
     x      +  EQDO * (
ccc   d (proton) x \bar{d} (Be(4,8))
     x        (DVN+DSN)* DBARBER  
ccc   \bar{d} (proton) x d (Be(4,8))
     x      +  DSN * DBER )
     x      +  EQDO * (
ccc   s (proton) x \bar{s} (Be(4,8))
     x         SSN * SSN  
ccc   \bar{s} (proton) x s (Be(4,8))
     x      +  SSN * SSN  ) 

ccc   QG initial states
       PDGNG    = F1COMBN * GLN 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CPM   FLUX RE-DEFINITION IN ORDER TO INCLUDE FF'S
CPM   construct fluxes incl. fragmentation
CPM   P-N-collisions
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
ccc   flux for SIGDQQ1: q qp --> q qp
ccc   also taking \bar{q} \bar{qp} --> \bar{q} \bar{qp}
ccc   into account
ccc   processes
       PDGNQQ1 = 
ccc   u/\bar{u}-quark coming from proton
     x            ((UVN+USN)* DBER 
     x        +    (UVN+USN)* SSN
     x        +         USN * DBARBER
     x        +         USN * SSN
ccc   u/\bar{u}-quark coming from Be(4,8)
     x        +    (DVN+DSN)* UBER
     x        +         SSN * UBER
     x        +         DSN * UBARBER
     x        +         SSN * UBARBER )* UN * EQUP 
ccc   d/\bar{d}-quark coming from proton
     x        +   ((DVN+DSN)* UBER
     x        +    (DVN+DSN)* SSN 
     x        +         DSN * UBARBER
     x        +         DSN * SSN 
ccc   d/\bar{d}-quark coming from Be(4,8)
     x        +    (UVN+USN)* DBER
     x        +         SSN * DBER
     x        +         USN * DBARBER
     x        +         SSN * DBARBER )* DN * EQDO 
ccc   s/\bar{s}-quark coming from proton
     x        +  (      SSN * UBER
     x        +         SSN * DBER
     x        +         SSN * UBARBER
     x        +         SSN * DBARBER
ccc   s/\bar{s}-quark coming from Be(4,8)
     x        +    (UVN+USN)* SSN
     x        +    (DVN+DSN)* SSN
     x        +         USN * SSN
     x        +         DSN * SSN )* SN * EQDO


ccc   flux for SIGDQQ2: q \bar{qp} --> q \bar{qp}
ccc   Then factor of 2.D0 here is for the 
ccc   crossing in the IS
       PDGNQQ2 =      ( 
ccc   u/\bar{u}-quark coming from proton
     x       (UVN+USN)* DBARBER 
     x      +(UVN+USN)* SSN 
     x      +     USN * DBER
     x      +     USN * SSN
ccc   u/\bar{u}-quark coming from Be(4,8)
     x      +     DSN * UBER
     x      +     SSN * UBER
     x      +(DVN+DSN)* UBARBER
     x      +     SSN * UBARBER ) * UN * EQUP
     x      +         (
ccc   d/\bar{d}-quark coming from proton       
     x      +(DVN+DSN)* UBARBER
     x      +(DVN+DSN)* SSN
     x      +     DSN * UBER
     x      +     DSN * SSN
ccc   d/\bar{d}-quark coming from Be(4,8)
     x      +     USN * DBER
     x      +     SSN * DBER
     x      +(UVN+USN)* DBARBER
     x      +     SSN * DBARBER ) * DN * EQDO
     x      +         (
ccc   s/\bar{s}-quark coming from proton       
     x      +     SSN * UBARBER
     x      +     SSN * DBARBER
     x      +     SSN * UBER
     x      +     SSN * DBER
ccc   s/\bar{s}-quark coming from Be(4,8)
     x      +     USN * SSN 
     x      +     DSN * SSN 
     x      +(UVN+USN)* SSN 
     x      +(DVN+DSN)* SSN ) * SN * EQDO  

ccc   flux for SIGDQQ3: q \bar{q} --> qp \bar{qp}
       PDGNQQ3 =  (
ccc   u from proton / \bar{u} from Be(4,8)
     x      (UVN+USN)* UBARBER
ccc   u from Be(4,8) / \bar{u} from proton
     x      + UBER * USN ) 
     x      * 2.D0 * (DN + SN) * EQDO
     x      +     (
ccc   d from proton / \bar{d} from Be(4,8)
     x      (DVN+DSN)* DBARBER
ccc   d from Be(4,8)/\bar{d} from proton
     x      + DBER * DSN ) 
     x      * 2.D0 * (UN * EQUP + SN * EQDO) 
     x      +     (
ccc   s from proton/\bar{s} from Be(4,8)
ccc   s from Be(4,8)/\bar{s} from proton
     x      + 2.D0 * SSN * SSN ) 
     x      * 2.D0 * (UN * EQUP + DN * EQDO) 

ccc   flux for SIGDQQ4: q q --> q q
       PDGNQQ4 =      (
ccc   u-u
     x       (UVN+USN)* UBER
ccc   \bar{u}-\bar{u}
     x      +     USN * UBARBER ) * 2.D0 * UN * EQUP 
     x      +         ( 
ccc   d-d
     x       (DVN+DSN)* DBER
ccc  \bar{d}-\bar{d} 
     x      +     DSN * DBARBER ) * 2.D0 * DN * EQDO  
ccc   s-s/\bar{s}-\bar{s}
     x      +    2.D0*( 
     x            SSN * SSN )  * 2.D0 * SN * EQDO 

ccc   flux for SIGDQQ5: q \bar{q} --> q \bar{q}
       PDGNQQ5 =       ( 
ccc   u (proton) x \bar{u} (Be(4,8))
     x        (UVN+USN) * UBARBER 
ccc   u (Be(4,8)) x \bar{u} (proton)
     x      +      USN  * UBER  )     
     x      * 2.D0 * UN * EQUP
     x      +          ( 
ccc   d (proton) x \bar{d} (Be(4,8))
     x        (DVN+DSN) * DBARBER
ccc   d (Be(4,8)) x \bar{d} (proton) 
     x      +      DSN  * DBER )     
     x      * 2.D0 * DN * EQDO 
ccc   s (proton) x \bar{s} (Be(4,8)) 
ccc   +
ccc   s (Be(4,8)) x \bar{s} (proton)
     x      +    2.D0*(
     x             SSN * SSN )  * 2.D0 * SN * EQDO 


ccc   gluons - not needed for now ccccccccccccccccc
ccc   flux for SIGDGG1: g g --> g g
c$$$       PDGNGG1 = GLN * GLN * 2. * GLNF
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c       gln = 1.d0
ccc   flux for SIGDGG2: g g --> q \bar{q}
       PDGNGG2 = GLN * GLN * 2.D0 
     x         * (UN * EQUP + (DN + SN) * EQDO)

ccccccccccccccccccccccccccccccccccccccccccccccccccc

ccc   flux for SIGDQG: q g --> q g / g q
ccc   NB! not taking gluon fragmentation into account
       PDGNGMOD =  GLN* ( ( 
ccc   gluon coming from proton
     x        (UBER + UBARBER)
ccc   gluon coming from Be(4,8)
     x      + (UVN + 2.D0*USN)   ) * UN * EQUP 
     x      +           (
ccc   gluon coming from proton       
     x        (DBER + DBARBER)
ccc   gluon coming from Be(4,8)
     x      + (DVN + 2.D0*DSN)   )* DN * EQDO  
     x      +           (
ccc   gluon coming from proton       
     x        (2.D0 * SSN)  
ccc   gluon coming from Be(4,8)
     x      + (2.D0 * SSN)       )* SN * EQDO)

      ENDIF

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CAM.....PDGNQ1,...PDGNQ6  ADDED BELOW IN PDG

       SUBROUTINE PDG(PDGNQ, PDGNG, PDGNGG, PDGNQQ,PDGNQ1,PDGNQ2,
     1       PDGNQ3,PDGNQ4,PDGNQ5,PDGNQ6, XN)
       IMPLICIT DOUBLE COMPLEX (A - Z)
       DOUBLE PRECISION EQUP, EQDO, ALPS, ALP0, ALP1, ALPC, 
     1                  ALPB, ALPT, ALPQ, ALPQR, LMQ, LMQR
       COMMON / COUPL / ALPS, ALP0, ALP1, ALPC, ALPB, ALPT, ALPQ, 
     1                  ALPQR, LMQ, LMQR
C
       CALL RENO(UVN, DVN, USN, DSN, SSN, GLN, XN)
       CSN = DCMPLX(0.D0,0.D0)
C
       EQUP = 4./9.
       EQDO = 1./9.
C
C ***  P-P :
       F1COMBN = EQUP*(UVN+2.*USN) + 
     1      EQDO*(DVN+2.*DSN+2.*SSN)
       PDGNQ = 2. * (EQUP * (UVN+USN)*USN + 
     1      EQDO * (DVN+DSN)*DSN +
     2      EQDO * SSN*SSN )
       PDGNG = 2. * F1COMBN * GLN 
       PDGNGQ = 2. * F1COMBN * GLN
       PDGNGG = GLN*GLN
       PDGNQ1 = EQUP*EQDO*((UVN+USN)*DSN+(DVN+DSN)*SSN)
     x      +EQDO*EQDO*(DVN+DSN)*SSN
       PDGNQ2 = EQUP*EQDO*((UVN+USN)*(DVN+DSN)+(UVN+USN)*2*SSN)
     x      +EQDO*EQDO*(DVN+DSN)*2*SSN
       PDGNQ3 = PDGNQ
       PDGNQ4 = PDGQQ
       PDGNQ5 = EQUP*EQUP*(UVN+USN)*USN
     x      +EQDO*EQDO*(DVN+DSN)*DSN+EQDO*EQDO*SSN*SSN
       PDGNQ6 = PDGNQ5


c$$$CCC   Combine with fragmentation components
c$$$       PDGNQ = 2. * (EQUP * (UVN+USN)*USN + 
c$$$     1      EQDO * (DVN+DSN)*DSN +
c$$$     2      EQDO * SSN*SSN )
c$$$       PDGNG = 2. * F1COMBN * GLN 
c$$$       PDGNGQ = 2. * F1COMBN * GLN
c$$$
c$$$ccc   (gg --> qq)
c$$$       PDGNGG1 = GLN*GLN * 2.*(UN + DN + SN + CN)
c$$$
c$$$ccc   (gg --> gg)
c$$$       PDGNGG2 = GLN*GLN * 2.*GLNF
c$$$
c$$$ccc   (
c$$$       PDGNQ1 = EQUP*EQDO*((UVN+USN)*DSN+(DVN+DSN)*SSN)
c$$$     x      +EQDO*EQDO*(DVN+DSN)*SSN
c$$$       PDGNQ2 = EQUP*EQDO*((UVN+USN)*(DVN+DSN)+(UVN+USN)*2*SSN)
c$$$     x      +EQDO*EQDO*(DVN+DSN)*2*SSN
c$$$       PDGNQ3 = PDGNQ
c$$$       PDGNQ4 = PDGQQ
c$$$       PDGNQ5 = EQUP*EQUP*(UVN+USN)*USN
c$$$     x      +EQDO*EQDO*(DVN+DSN)*DSN+EQDO*EQDO*SSN*SSN
c$$$       PDGNQ6 = PDGNQ5       
C
C ***  P-PBAR :
c$$$        F1COMBN = EQUP*(UVN+2.*USN) + 
c$$$     1            EQDO*(DVN+2.*DSN+2.*SSN)
c$$$        PDGNQ = EQUP * ( (UVN+USN)*(UVN+USN) + USN*USN ) +
c$$$     1          EQDO * ( (DVN+DSN)*(DVN+DSN) + DSN*DSN +
c$$$     2       2.*SSN*SSN )
c$$$        PDGNG = 2. * F1COMBN * GLN 
c$$$C	    PDGNGQ = (EQUP*USN+EQDO*(DSN+SSN))* GLN
c$$$        PDGNQQ = 2.*(EQUP*(UVN+USN)*USN+EQDO*((DVN+DSN)*DSN +2.*SSN*SSN)
c$$$     1             + EQUP*USN*(UVN+USN)+EQDO*(DSN*(DVN+DSN)+2.*SSN*SSN))
c$$$     	PDGNGG = 2.D0*GLN*GLN
c$$$        PDGNQ1 = EQUP*EQDO*((UVN+USN)*DVN+(UVN+USN)*2*SSN)
c$$$C                 +EQDO*EQDO*(DVN+DSN)*2*SSN
c$$$        PDGNQ2 = EQUP*EQDO*((UVN+USN)*(DVN+DSN)+(UVN+USN)*SSN)
c$$$C                 +EQDO*EQDO*(DVN+DSN)*SSN
c$$$        PDGNQ3 = PDGNQ
c$$$        PDGNQ4 = PDGQQ
c$$$        PDGNQ5 =EQUP*EQUP*(UVN+USN)*UVN+EQDO*EQDO*(DVN+DSN)*DVN
c$$$     1       +EQDO*EQDO*(SSN*SSN) 
c$$$        PDGNQ6 = PDGNQ5
C
C ***  P-N :
*       F1COMBN = EQUP*(3./2.*(UVN+2.*USN)+(DVN+2.*DSN)/2.+4.*CSN) + 
*     1           EQDO*(3./2.*(DVN+2.*DSN)+(UVN+2.*USN)/2.+4.*SSN)
*       PDGNQ = ( EQUP * ( (UVN+USN)*(USN+DSN)/2. + 
*     1                     USN*(UVN+USN+DVN+DSN)/2. + 2.*CSN*CSN ) +
*     2           EQDO * ( (DVN+DSN)*(USN+DSN)/2. + 
*     3                     DSN*(UVN+USN+DVN+DSN)/2. + 2.*SSN*SSN ) )
*       PDGNG = F1COMBN * GLN 
c        PDGNQ1 =
c        PDGNQ2 = 
c        PDGNQ3 = 
c        PDGNQ4 = 
c        PDGNQ5 = 
c        PDGNQ6 = 
         
C
       RETURN
       END



c$$$        write(*,*)"UVN   = ", uvn
c$$$        write(*,*)"USN   = ", usn
c$$$        write(*,*)"DVN   = ", dvn
c$$$        write(*,*)"DSN   = ", dsn
c$$$        write(*,*)"SSN   = ", ssn
c$$$        write(*,*)"GLN   = ", gln
c$$$
c$$$        write(*,*)"UN    = ", un
c$$$        write(*,*)"DN    = ", dn
c$$$        write(*,*)"SN    = ", sn



