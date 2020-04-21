***************************************************************
	FUNCTION INTVEG(FU,INTDIM,NCALL,ITMX)
***************************************************************
cpm	implicit real*8 (a-z)
	implicit double precision (a-z)
	integer nhist,nbin
	parameter(nhist=7,nbin=2000)
cpm	real*8  bin(nhist,0:20,0:nbin),binerr(nhist,0:20,0:nbin)
cpm	real*8  x(nhist),xmin(nhist),xmax(nhist)
	double precision bin(nhist,0:20,0:nbin),binerr(nhist,0:20,0:nbin)
	double precision x(nhist),xmin(nhist),xmax(nhist)
        integer act,comp,binpop(nhist,0:20,0:nbin),nx(nhist),ln(nhist)
        character*20 fname(nhist)

	integer i,NDIM,NCALL,ITMX,NPRN,INTDIM

	logical indy, prnt

        COMMON/BVEG1/XL(20),XU(20),ACC,NDIM,NPRN
        common/resultv/s1,s2,s3,s4
        common/hist/bin,binerr,xmin,xmax,binpop,nx,ncomp,fname
	common/crapo/ indy, prnt

  
	EXTERNAL FU

C      OPEN(11,FILE='daten.out',STATUS='NEW')

c	random number initialization
c	call rmarin(1802,9373)
	call rmarin(1800,9370)
c	call rmarin(1801,9371)

C--> MONTE CARLO INTEGRATION PARAMETERS:                               
        NPRN = 10
C--> NUMBER OF DIMENSIONS FOR MONTE CARLO INTEGRATION                   
        NDIM = INTDIM                                                   
C--> EVALUATIONS PER ITERATION                                          
C       NCALL = 100000                                                  
C--> MAX. NUMBER OF ITERATIONS                                          
C       ITMX = 8                                                        
C--> DEMANDED ACCURACY                                                  
        ACC = 1d-10
C--> INTEGRATION BOUNDARY
cpm	do 10 i=1,10
	do 10 i=1,20
        XL(i) = 0D0                                                     
        XU(i) = 1D0                                                     
10	continue

c       CALL VEGAN(FU,AVGI,SD,CHIA2,NCALL,ITMX)
	indy = .true.
	CALL VEGAN(FU,AVGI,SD,CHIA2,NCALL,itmx) 

cpm	call histogram(-1,x,comp,add,count)
	indy = .false.
        CALL VEGAN1(FU,AVGI,SD,CHIA2,2*NCALL,ITMX) 

	INTVEG = S1
c	write(*,*) "sigma = ", s2

	end

***************************************************************
	subroutine vegan(fxn,avgi,sd,chi2a,ncall,itmx)
***************************************************************
*      vegas, corrected version from jos vermaseren
*
***************************************************************
*    25.03.91 sa
***************************************************************
        implicit double precision ( a-h,o-z )
	external fxn

	integer ncall, itmx, i, nxi

cpm       dimension xin(204),r(204),dx(20),ia(20),kg(20),dt(20)
       dimension xin(102),r(102),dx(20),ia(20),kg(20),dt(20)
cpm       dimension xin(51),r(51),dx(20),ia(20),kg(20),dt(20)
cpm       dimension xin(51),r(51),dx(10),ia(10),kg(10),dt(10)
cpm       dimension qran(10),x(10)
	dimension qran(20),x(20)
	common/bveg1/xl(20),xu(20),acc,ndim,nprn
csd      common/bveg2/xi(51,10),si,si2,swgt,schi,ndo,it

cpm	common/bveg2/xi(51,20),si,si2,swgt,schi,ndo,it
cpm	common/bveg5/scalls
cpm     + ,d(51,20),di(51,20),nxi(51,20)

	common/bveg2/xi(102,20),si,si2,swgt,schi,ndo,it
	common/bveg5/scalls
     + ,d(102,20),di(102,20),nxi(102,20)

cpm      common/bveg2/xi(204,20),si,si2,swgt,schi,ndo,it
cpm       common/bveg5/scalls
cpm     + ,d(204,20),di(204,20),nxi(204,20)
c$$$       common/bveg5/scalls
c$$$     + ,d(51,10),di(51,10),nxi(51,10)
	common/resultv/s1,s2,s3,s4
cpm	data ndmx/50/,alph/1.5/,one/1./,mds/1/

	data ndmx/100/,alph/1.5/,one/1./,mds/1/

cpm      data ndmx/200/,alph/1.5/,one/1./,mds/1/

	logical dump
	logical indy, prnt

	common /scrdump/ dump
	common/crapo/ indy, prnt

	
	dump = .false.
	prnt = .false.
      ipr=1
      if(nprn.gt.0)ipr=0
      ndo=1
      do 1 j=1,ndim
1     xi(1,j)=one
      entry vegan1(fxn,avgi,sd,chi2a,ncall,itmx)
      it=0
      si=0.
      si2=si
      swgt=si
      schi=si
      scalls=si
      entry vegan2(fxn,avgi,sd,chi2a,ncall,itmx)
      nd=ndmx
      ng=1
      if(mds.eq.0) go to 2
      ng=(ncall*0.5d+00)**(1.d+00/ndim)
      mds=1
      if((2*ng-ndmx).lt.0) go to 2
      mds=-1
      npg=ng/ndmx+1
      nd=ng/npg
      ng=npg*nd
2     k=ng**ndim
      npg=ncall/k
      if(npg.lt.2)npg=2
      calls=npg*k
      dxg=one/ng
      dv2g=dxg**(2*ndim)/npg/npg/(npg-one)
      xnd=nd
      ndm=nd-1
      dxg=dxg*xnd
      xjac=one
      do 3 j=1,ndim
      dx(j)=xu(j)-xl(j)
3     xjac=xjac*dx(j)
      if(nd.eq.ndo) go to 8
      rc=ndo/xnd
      do 7 j=1,ndim
      k=0
      xn=0.
      dr=xn
      i=k
4     k=k+1
      dr=dr+one
      xo=xn
      xn=xi(k,j)
5     if(rc.gt.dr) go to 4
      i=i+1
      dr=dr-rc
      xin(i)=xn-(xn-xo)*dr
      if(i.lt.ndm) go to 5
      do 6  i=1,ndm
6     xi(i,j)=xin(i)
7     xi(nd,j)=one
      ndo=nd
8     if(nprn.ne.0.and.nprn.ne.10)print 200,ndim,calls,it,itmx
     1,acc,mds,nd
      if(nprn.eq.10)print 290,ndim,calls,itmx,acc,mds,nd
      entry vegan3(fxn,avgi,sd,chi2a,ncall,itmx)
9     it=it+1
cpm   printing after 15th iteration
c$$$	write(33,*)"VEGAN iter. no. ", it
	if((it.eq.9 ) .and. 
     1  (indy.eqv. .true. ) ) prnt = .true.
c	if(it.gt.4 ) stop
      ti=0.
      tsi=ti
      do 10 j=1,ndim
      kg(j)=1
      do 10 i=1,nd
      nxi(i,j)=0
      d(i,j)=ti
10    di(i,j)=ti
11    fb=0.
      f2b=fb
      k=0
12    k=k+1
      call aran9(qran,ndim)
      wgt=xjac
      do 15 j=1,ndim
      xn=(kg(j)-qran(j))*dxg+one
      ia(j)=xn
      iaj=ia(j)
      iaj1=iaj-1
      if(iaj.gt.1) go to 13
      xo=xi(iaj,j)
      rc=(xn-iaj)*xo
      go to 14
13    xo=xi(iaj,j)-xi(iaj1,j)
      rc=xi(iaj1,j)+(xn-iaj)*xo
14    x(j)=xl(j)+rc*dx(j)
15    wgt=wgt*xo*xnd
csd start
      do 900 j=1,ndim
	if (x(j).ge.1d0) x(j) = 1d0-1d-15
	if (x(j).le.0d0) x(j) = 1d-15
900   continue
csd end
	f=fxn(x,wgt)*wgt
      f2=f*f
      fb=fb+f
      f2b=f2b+f2
      do 16 j=1,ndim
      iaj=ia(j)
      nxi(iaj,j)=nxi(iaj,j)+1
      di(iaj,j)=di(iaj,j)+f/calls
16    if(mds.ge.0)  d(iaj,j)=d(iaj,j)+f2
      if(k.lt.npg) go to 12
      f2b=f2b*npg
      f2b=dsqrt(f2b)
      f2b=(f2b-fb)*(f2b+fb)
      ti=ti+fb
      tsi=tsi+f2b
      if(mds.ge.0) go to 18
      do 17 j=1,ndim
      iaj=ia(j)
17    d(iaj,j)=d(iaj,j)+f2b
18    k=ndim
19    kg(k)=mod(kg(k),ng)+1
      if(kg(k).ne.1) go to 11
      k=k-1
      if(k.gt.0) go to 19
      ti=ti/calls
      tsi=tsi*dv2g
      ti2=ti*ti
      wgt=ti2/tsi
      si=si+ti*wgt
      si2=si2+ti2
      swgt=swgt+wgt
      schi=schi+ti2*wgt
      scalls=scalls+calls
      avgi=si/swgt
      sd=swgt*it/si2
      chi2a=0.
      if(it.gt.1)chi2a=sd*(schi/swgt-avgi*avgi)/(it-1)
      sd=one/sd
      sd=dsqrt(sd)
      if(nprn.eq.0) go to 21
      tsi=dsqrt(tsi)
      if(nprn.ne.10)print 201,ipr,it,ti,tsi,avgi,sd,chi2a
      if(nprn.eq.10)print 203,it,ti,tsi,avgi,sd,chi2a
      if(nprn.ge.0) go to 21
      do 20 j=1,ndim
      print 202,j
20    print 204,(xi(i,j),di(i,j),d(i,j),i=1,nd)
21    continue
      s1=avgi
      s2=sd
      s3=ti
      s4=tsi
c      do 23 j=1,ndim
c      xo=d(1,j)
c      xn=d(2,j)
c      d(1,j)=(xo+xn)*0.5
c      dt(j)=d(1,j)
c      do 22 i=2,ndm
c      d(i,j)=xo+xn
c      xo=xn
c      xn=d(i+1,j)
c      d(i,j)=(d(i,j)+xn)/3.
c22    dt(j)=dt(j)+d(i,j)
c      d(nd,j)=(xn+xo)*0.5
c23    dt(j)=dt(j)+d(nd,j)
c-----this part of the vegas-algorithm is unstable
c-----it should be replaced by
      do 23 j=1,ndim
      dt(j)=0.
      do 23 i=1,nd
      if(nxi(i,j).gt.0)d(i,j)=d(i,j)/nxi(i,j)
23    dt(j)=dt(j)+d(i,j)
      do 28 j=1,ndim
      rc=0.
      do 24 i=1,nd
      r(i)=0.
      if(d(i,j).le.0.)go to 24
      xo=dt(j)/d(i,j)
      r(i)=((xo-one)/xo/dlog(xo))**alph
24    rc=rc+r(i)
      rc=rc/xnd
      k=0
      xn=0.
      dr=xn
      i=k
25    k=k+1
      dr=dr+r(k)
      xo=xn
      xn=xi(k,j)
26    if(rc.gt.dr) go to 25
      i=i+1
      dr=dr-rc
      xin(i)=xn-(xn-xo)*dr/r(k)
      if(i.lt.ndm) go to 26
      do 27 i=1,ndm
27    xi(i,j)=xin(i)
28    xi(nd,j)=one
      if(it.lt.itmx.and.dabs(acc).lt.dabs(sd/avgi))go to 9
200   format(35h0input parameters for vegas   ndim=,i3
     1,8h  ncall=,f8.0/28x,5h  it=,i5,8h  itmx =,i5/28x
     2,6h  acc=,g9.3/28x,6h  mds=,i3,6h   nd=,i4//)
290   format(13h0vegas  ndim=,i3,8h  ncall=,f8.0,8h  itmx =,i5
     1,6h  acc=,g9.3,6h  mds=,i3,6h   nd=,i4)
201   format(/i1,20hintegration by vegas/13h0iteration no,i3,
     114h.   integral =,g14.8/20x,10hstd dev  =,g10.4/
     234h accumulated results.   integral =,g14.8/
     324x,10hstd dev  =,g10.4 / 24x,18hchi**2 per itn   =,g10.4)
202   format(14h0data for axis,i2 / 7x,1hx,7x,10h  delt i  ,
     12x,11h convce    ,11x,1hx,7x,10h  delt i  ,2x,11h convce
     2,11x,1hx,7x,10h  delt i  ,2x,11h convce     /)
204   format(1x,3g12.4,5x,3g12.4,5x,3g12.4)
203   format(1h ,i3,g20.8,g12.4,g20.8,g12.4,g12.4)
      s1=avgi
      s2=sd
      s3=chi2a
      return
      end
*
      SUBROUTINE ARAN9(QRAN,NDIM)
c      REAL*8 QRAN(10)
c      REAL*8 QRAN(20)
      double precision QRAN(20)
      DO 1 I=1,NDIM
    1 CALL RANMAR(QRAN(I))
c    1 CALL R2455(QRAN(I))
      RETURN
      END
*
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE RANMAR(RVEC)
*     -----------------
* Universal random number generator proposed by Marsaglia and Zaman
* in report FSU-SCRI-87-50
* In this version RVEC is a double precision variable.
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/ RASET1 / RANU(97),RANC,RANCD,RANCM
      COMMON/ RASET2 / IRANMR,JRANMR
      UNI = RANU(IRANMR) - RANU(JRANMR)
      IF(UNI .LT. 0D0) UNI = UNI + 1D0
      RANU(IRANMR) = UNI
      IRANMR = IRANMR - 1
      JRANMR = JRANMR - 1
      IF(IRANMR .EQ. 0) IRANMR = 97
      IF(JRANMR .EQ. 0) JRANMR = 97
      RANC = RANC - RANCD
      IF(RANC .LT. 0D0) RANC = RANC + RANCM
      UNI = UNI - RANC
      IF(UNI .LT. 0D0) UNI = UNI + 1D0
      RVEC = UNI
      END

      SUBROUTINE RMARIN(IJ,KL)
*     -----------------
* Initializing routine for RANMAR, must be called before generating
* any pseudorandom numbers with RANMAR. The input values should be in
* the ranges 0<=ij<=31328 ; 0<=kl<=30081
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/ RASET1 / RANU(97),RANC,RANCD,RANCM
      COMMON/ RASET2 / IRANMR,JRANMR
* This shows correspondence between the simplified input seeds IJ, KL
* and the original Marsaglia-Zaman seeds I,J,K,L.
* To get the standard values in the Marsaglia-Zaman paper (i=12,j=34
* k=56,l=78) put ij=1802, kl=9373
      I = MOD( IJ/177 , 177 ) + 2
      J = MOD( IJ     , 177 ) + 2
      K = MOD( KL/169 , 178 ) + 1
      L = MOD( KL     , 169 )
      WRITE(*,100) IJ,KL,I,J,K,L
  100 FORMAT(' Initialization for RANMAR: ',2I8,4I5)
      DO 300 II = 1 , 97
	S =  0D0
	T = .5D0
	DO 200 JJ = 1 , 24
	  M = MOD( MOD(I*J,179)*K , 179 )
	  I = J
	  J = K
	  K = M
	  L = MOD( 53*L+1 , 169 )
	  IF(MOD(L*M,64) .GE. 32) S = S + T
	  T = .5D0*T
  200   CONTINUE
	RANU(II) = S
  300 CONTINUE
      RANC  =   362436D0 / 16777216D0
      RANCD =  7654321D0 / 16777216D0
      RANCM = 16777213D0 / 16777216D0
      IRANMR = 97
      JRANMR = 33
      END

