cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function c1_qq_qq(q2,mu2,muf2,mur2,mud2,nf)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none

      integer nf

      double precision CA, CF, zeta2, game, pi,
     x     mu, muf, mur, mud, mufp2, q2,
     x     mu2, muf2, mur2, mud2

      parameter(CA = 3.d0)
      parameter(CF = 4.d0/3.d0)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   C0de starts here...
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   scales and constants 
      pi    = dacos(-1.d0)
c      mu2   = mu**2
c      muf2  = muf**2
c      mufp2 = mud**2
      mufp2 = mud2
c      mur2  = mur**2
      zeta2 = pi**2/6.d0
      game  = 0.5772156649d0
      game  = 0.d0 * game
cccccccccccccccccccccccccccc

      c1_qq_qq = 
     x     -(900 - 1125*CA + 460*CA**2 - 575*CA**3 + 216*CA*CF*game 
     x     - 270*CA**2*CF*game + 720*CA*CF*game**2 
     x     - 900*CA**2*CF*game**2 - 160*CA*nf + 200*CA**2*nf + 78*Pi**2 
     x     - 138*CA*Pi**2 + 120*CA**2*Pi**2 - 96*CA**3*Pi**2 
     x     + 720*CA*CF*zeta2 - 900*CA**2*CF*zeta2 + 324*Log(2.d0) 
     x     - 333*CA*Log(2.d0) + 204*CA**2*Log(2.d0) 
     x     - 255*CA**3*Log(2.d0) 
     x     - 1296*game*Log(2.d0) + 2340*CA*game*Log(2.d0) 
     x     + 144*CA**2*game*Log(2.d0) - 900*CA**3*game*Log(2.d0) 
     x     - 96*CA*nf*Log(2.d0) + 120*CA**2*nf*Log(2.d0) 
     x     - 216*Log(2.d0)**2 + 162*CA*Log(2.d0)**2 
     x     + 216*CA**2*Log(2.d0)**2 - 270*CA**3*Log(2.d0)**2 
     x     + 36*CA*(-4 + 5*CA)*CF*(-3 + 4*game)*Log(Q2/muF2) 
     x     + 18*CA*(-4 + 5*CA)*CF*(-3 + 4*game 
     x     + Log(16.d0))*Log(Q2/muFp2) - 528*CA**2*Log(Q2/muR2) 
     x     + 660*CA**3*Log(Q2/muR2) + 96*CA*nf*Log(Q2/muR2) 
     x     - 120*CA**2*nf*Log(Q2/muR2))/(792.*CA)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function c1_qqp_qqp(q2,mu2,muf2,mur2,mud2,nf)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none

      integer nf

      double precision CA, CF, zeta2, game, pi,
     x     mu, muf, mur, mud, mufp2, q2,
     x     mu2, muf2, mur2, mud2

      parameter(CA = 3.d0)
      parameter(CF = 4.d0/3.d0)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   C0de starts here...
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   scales and constants 
      pi    = dacos(-1.d0)
c      mu2   = mu**2
c      muf2  = muf**2
c      mufp2 = mud**2
      mufp2 = mud2
c      mur2  = mur**2
      zeta2 = pi**2/6.d0
      game  = 0.5772156649d0
      game  = 0.d0 * game
cccccccccccccccccccccccccccc
      
      c1_qqp_qqp = 
     x     (1125 + 575*CA**2 + 270*CA*CF*game + 900*CA*CF*game**2 
     x     - 200*CA*nf + 138*Pi**2 + 96*CA**2*Pi**2 + 900*CA*CF*zeta2 
     x     + 333*Log(2.d0) + 255*CA**2*Log(2.d0) - 2340*game*Log(2.d0) 
     x     + 900*CA**2*game*Log(2.d0) - 120*CA*nf*Log(2.d0) 
     x     - 162*Log(2.d0)**2 
     x     + 270*CA**2*Log(2.d0)**2 + 540*CA*CF*Log(Q2/muF2) 
     x     - 720*CA*CF*game*Log(Q2/muF2) 
     x     + 270*CA*CF*Log(Q2/muFp2) 
     x     - 360*CA*CF*game*Log(Q2/muFp2) 
     x     - 360*CA*CF*Log(2.d0)*Log(Q2/muFp2) 
     x     - 660*CA**2*Log(Q2/muR2) 
     x     + 120*CA*nf*Log(Q2/muR2))/(360.d0*CA)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function c1_qqb_qpqbp(q2,mu2,muf2,mur2,mud2,nf)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none

      integer nf

      double precision CA, CF, zeta2, game, pi,
     x     mu, muf, mur, mud, mufp2, s, q2,
     x     mu2, muf2, mur2, mud2

      parameter(CA = 3.d0)
      parameter(CF = 4.d0/3.d0)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   C0de starts here...
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   scales and constants 
      pi    = dacos(-1.d0)
c      mu2   = mu**2
c      muf2  = muf**2
c      mufp2 = mud**2
      mufp2 = mud2
c      mur2  = mur**2
      zeta2 = pi**2/6.d0
      game  = 0.5772156649d0
      game  = 0.d0 * game
cccccccccccccccccccccccccccc

      c1_qqb_qpqbp = 
     x     (225 + 115*CA**2 + 54*CA*CF*game + 180*CA*CF*game**2 
     x     - 40*CA*nf - 30*Pi**2 - 6*CA**2*Pi**2 + 180*CA*CF*zeta2 
     x     - 27*Log(2.d0) - 9*CA**2*Log(2.d0) - 36*game*Log(2.d0) 
     x     + 180*CA**2*game*Log(2.d0) - 18*Log(2.d0)**2 
     x     + 90*CA**2*Log(2.d0)**2 + 108*CA*CF*Log(Q2/muF2) 
     x     - 144*CA*CF*game*Log(Q2/muF2) 
     x     + 54*CA*CF*Log(Q2/muFp2) 
     x     - 72*CA*CF*game*Log(Q2/muFp2) 
     x     - 72*CA*CF*Log(2.d0)*Log(Q2/muFp2) 
     x     - 132*CA**2*Log(Q2/muR2) 
     x     + 24*CA*nf*Log(Q2/muR2))/(72.d0*CA)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function c1_qqb_qqb(q2,mu2,muf2,mur2,mud2,nf)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none 

      integer nf

      double precision CA, CF, zeta2, game, pi,
     x     mu, muf, mur, mud, mufp2, q2,
     x     mu2, muf2, mur2, mud2

      parameter(CA = 3.d0)
      parameter(CF = 4.d0/3.d0)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   C0de starts here...
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   scales and constants 
      pi    = dacos(-1.d0)
c      mu2   = mu**2
c      muf2  = muf**2
c      mufp2 = mud**2
      mufp2 = mud2
c      mur2  = mur**2
      zeta2 = pi**2/6.d0
      game  = 0.5772156649d0
      game  = 0.d0 * game
cccccccccccccccccccccccccccc

      c1_qqb_qqb = 
     x     (450*CF + 2475*CA*CF + 230*CA**2*CF + 1265*CA**3*CF 
     x     + 108*CA*CF**2*game + 594*CA**2*CF**2*game 
     x     + 360*CA*CF**2*game**2 + 1980*CA**2*CF**2*game**2 
     x     - 80*CA*CF*nf - 440*CA**2*CF*nf - 24*CF*Pi**2 
     x     - 186*CA*CF*Pi**2 + 24*CA**2*CF*Pi**2 + 294*CA**3*CF*Pi**2 
     x     + 360*CA*CF**2*zeta2 + 1980*CA**2*CF**2*zeta2 
     x     + 162*CF*Log(2.d0) 
     x     + 927*CA*CF*Log(2.d0) + 78*CA**2*CF*Log(2.d0) 
     x     + 429*CA**3*CF*Log(2.d0) 
     x     - 72*CF*game*Log(2.d0) + 1044*CA*CF*game*Log(2.d0) 
     x     + 360*CA**2*CF*game*Log(2.d0) + 540*CA**3*CF*game*Log(2.d0) 
     x     - 24*CA*CF*nf*Log(2.d0) - 240*CA**2*CF*nf*Log(2.d0) 
     x     - 126*CF*Log(2.d0)**2 - 774*CA*CF*Log(2.d0)**2 
     x     + 144*CA**2*CF*Log(2.d0)**2 + 738*CA**3*CF*Log(2.d0)**2 
     x     + 216*CA*CF**2*Log(Q2/muF2) 
     x     + 1188*CA**2*CF**2*Log(Q2/muF2) 
     x     - 288*CA*CF**2*game*Log(Q2/muF2) 
     x     - 1584*CA**2*CF**2*game*Log(Q2/muF2) 
     x     + 108*CA*CF**2*Log(Q2/muFp2) 
     x     + 594*CA**2*CF**2*Log(Q2/muFp2) 
     x     - 144*CA*CF**2*game*Log(Q2/muFp2) 
     x     - 792*CA**2*CF**2*game*Log(Q2/muFp2) 
     x     - 144*CA*CF**2*Log(2.d0)*Log(Q2/muFp2) 
     x     - 792*CA**2*CF**2*Log(2.d0)*Log(Q2/muFp2) 
     x     - 264*CA**2*CF*Log(Q2/muR2) 
     x     - 1452*CA**3*CF*Log(Q2/muR2) 
     x     + 48*CA*CF*nf*Log(Q2/muR2) 
     x     + 264*CA**2*CF*nf*Log(Q2/muR2))/(3360.d0*CA)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function c1_qqbp_qqbp(q2,mu2,muf2,mur2,mud2,nf)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none

      integer nf

      double precision CA, CF, zeta2, game, pi,
     x     mu, muf, mur, mud, mufp2, q2,
     x     mu2, muf2, mur2, mud2

      parameter(CA = 3.d0)
      parameter(CF = 4.d0/3.d0)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   C0de starts here...
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   scales and constants 
      pi    = dacos(-1.d0)
c      mu2   = mu**2
c      muf2  = muf**2
c      mufp2 = mud**2
      mufp2 = mud2
c      mur2  = mur**2
      zeta2 = pi**2/6.d0
      game  = 0.5772156649d0
      game  = 0.d0 * game
cccccccccccccccccccccccccccc

      c1_qqbp_qqbp = 
     x     (1125 + 575*CA**2 + 270*CA*CF*game + 900*CA*CF*game**2 
     x     - 200*CA*nf - 78*Pi**2 + 150*CA**2*Pi**2 + 900*CA*CF*zeta2 
     x     + 477*Log(2.d0) + 219*CA**2*Log(2.d0) + 540*game*Log(2.d0) 
     x     + 180*CA**2*game*Log(2.d0) - 120*CA*nf*Log(2.d0) 
     x     - 378*Log(2.d0)**2 
     x     + 324*CA**2*Log(2.d0)**2 + 540*CA*CF*Log(Q2/muF2) 
     x     - 720*CA*CF*game*Log(Q2/muF2) 
     x     + 270*CA*CF*Log(Q2/muFp2) 
     x     - 360*CA*CF*game*Log(Q2/muFp2) 
     x     - 360*CA*CF*Log(2.d0)*Log(Q2/muFp2) 
     x     - 660*CA**2*Log(Q2/muR2) 
     x     + 120*CA*nf*Log(Q2/muR2))/(360.d0*CA)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function c1_qg_qg(q2,mu2,muf2,mur2,mud2,nf)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none

      integer nf
      
      double precision CA, CF, zeta2, game, pi,
     x     mu, muf, mur, mud, mufp2, q2,
     x     mu2, muf2, mur2, mud2

      parameter(CA = 3.d0)
      parameter(CF = 4.d0/3.d0)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   C0de starts here...
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   scales and constants 
      pi    = dacos(-1.d0)
c      mu2   = mu**2
c      muf2  = muf**2
c      mufp2 = mud**2
      mufp2 = mud2
c      mur2  = mur**2
      zeta2 = pi**2/6.d0
      game  = 0.5772156649d0
      game  = 0.d0 * game
cccccccccccccccccccccccccccc
      c1_qg_qg =
     x     (-630*CF + 3434*CA**2*CF + 524*CA**4*CF + 990*CA**2*CF*game 
     x     + 330*CA**4*CF*game + 2640*CA**3*CF**2*game - 1080*CF*game**2 
     x     + 1260*CA**2*CF*game**2 + 540*CA**4*CF*game**2 
     x     - 2880*CA*CF**2*game**2 + 4320*CA**3*CF**2*game**2 
     x     - 300*CA*CF*nf - 100*CA**3*CF*nf - 800*CA**2*CF**2*nf 
     x     - 180*CA*CF*game*nf - 60*CA**3*CF*game*nf 
     x     - 480*CA**2*CF**2*game*nf + 36*CF*Pi**2 - 156*CA**2*CF*Pi**2 
     x     + 312*CA**4*CF*Pi**2 - 1080*CF*zeta2 + 1260*CA**2*CF*zeta2 
     x     + 540*CA**4*CF*zeta2 - 2880*CA*CF**2*zeta2 
     x     + 4320*CA**3*CF**2*zeta2 - 288*CF*Log(2.d0) 
     x     + 1146*CA**2*CF*Log(2.d0) + 174*CA**4*CF*Log(2.d0) 
     x     + 720*CF*game*Log(2.d0) - 3960*CA**2*CF*game*Log(2.d0) 
     x     + 4680*CA**4*CF*game*Log(2.d0) - 180*CA*CF*nf*Log(2.d0) 
     x     - 60*CA**3*CF*nf*Log(2.d0) - 480*CA**2*CF**2*nf*Log(2.d0) 
     x     + 360*CF*Log(2.d0)**2 - 1440*CA**2*CF*Log(2.d0)**2 
     x     + 1548*CA**4*CF*Log(2.d0)**2 - 405*CF*Log(Q2/muF2) 
     x     + 1260*CA**2*CF*Log(Q2/muF2) 
     x     + 465*CA**4*CF*Log(Q2/muF2) 
     x     - 1080*CA*CF**2*Log(Q2/muF2) 
     x     + 3720*CA**3*CF**2*Log(Q2/muF2) 
     x     + 540*CF*game*Log(Q2/muF2) 
     x     - 1440*CA**2*CF*game*Log(Q2/muF2) 
     x     - 540*CA**4*CF*game*Log(Q2/muF2) 
     x     + 1440*CA*CF**2*game*Log(Q2/muF2) 
     x     - 4320*CA**3*CF**2*game*Log(Q2/muF2) 
     x     - 180*CA*CF*nf*Log(Q2/muF2) 
     x     - 60*CA**3*CF*nf*Log(Q2/muF2) 
     x     - 480*CA**2*CF**2*nf*Log(Q2/muF2) 
     x     + 810*CA*CF**2*Log(Q2/muFp2) 
     x     + 270*CA**3*CF**2*Log(Q2/muFp2) 
     x     + 2160*CA**2*CF**3*Log(Q2/muFp2) 
     x     - 1080*CA*CF**2*game*Log(Q2/muFp2) 
     x     - 360*CA**3*CF**2*game*Log(Q2/muFp2) 
     x     - 2880*CA**2*CF**3*game*Log(Q2/muFp2) 
     x     - 1080*CA*CF**2*Log(2.d0)*Log(Q2/muFp2) 
     x     - 360*CA**3*CF**2*Log(2.d0)*Log(Q2/muFp2) 
     x     - 2880*CA**2*CF**3*Log(2.d0)*Log(Q2/muFp2) 
     x     - 1980*CA**2*CF*Log(Q2/muR2) 
     x     - 660*CA**4*CF*Log(Q2/muR2) 
     x     - 5280*CA**3*CF**2*Log(Q2/muR2) 
     x     + 360*CA*CF*nf*Log(Q2/muR2) 
     x     + 120*CA**3*CF*nf*Log(Q2/muR2) 
     x     + 960*CA**2*CF**2*nf*Log(Q2/muR2))/(21120.d0*CA)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc
ccc   Omit this one in code: gluon fragmenting, q --> jet
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function c1_qg_gq(q2,mu2,muf2,mur2,mud2,nf)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none

      integer nf

      double precision CA, CF, zeta2, game, pi,
     x     mu, muf, mur, mud, mufp2, q2,
     x     mu2, muf2, mur2, mud2

      parameter(CA = 3.d0)
      parameter(CF = 4.d0/3.d0)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   C0de starts here...
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   scales and constants 
      pi    = dacos(-1.d0)
c      mu2   = mu**2
c      muf2  = muf**2
c      mufp2 = mud**2
      mufp2 = mud2
c      mur2  = mur**2
      zeta2 = pi**2/6.d0
      game  = 0.5772156649d0
      game  = 0.d0 * game
cccccccccccccccccccccccccccc
      c1_qg_gq =
     x     ( -105*CF + 738*CA**2*CF - 417*CA**4*CF 
     x       - 90*CA*CF**2*game + 450*CA**3*CF**2*game 
     x       + 30*CF*game**2 - 420*CA**2*CF*game**2 
     x       + 1350*CA**4*CF*game**2 + 2*CF*Pi**2 
     x       - 12*CA**2*CF*Pi**2 + 154*CA**4*CF*Pi**2 
     x       + 30*CF*zeta2 - 420*CA**2*CF*zeta2 
     x       + 1350*CA**4*CF*zeta2 - 51*CF*Log(2.d0) 
     x       + 222*CA**2*CF*Log(2.d0) - 267*CA**4*CF*Log(2.d0) 
     x       + 60*CF*game*Log(2.d0) - 600*CA**2*CF*game*Log(2.d0) 
     x       + 2460*CA**4*CF*game*Log(2.d0) + 30*CF*Log(2.d0)**2 
     x       - 120*CA**2*CF*Log(2.d0)**2 
     x       + 966*CA**4*CF*Log(2.d0)**2 
     x       + 45*CF*Log(Q2/muF2) 
     x       + 100*CA*CF*Log(Q2/muF2) 
     x       - 380*CA**2*CF*Log(Q2/muF2) 
     x       - 500*CA**3*CF*Log(Q2/muF2) 
     x       + 775*CA**4*CF*Log(Q2/muF2) 
     x       - 60*CF*game*Log(Q2/muF2) 
     x       + 480*CA**2*CF*game*Log(Q2/muF2) 
     x       - 900*CA**4*CF*game*Log(Q2/muF2) 
     x       + 100*CA*CF*Log(Q2/muFp2) 
     x       - 110*CA**2*CF*Log(Q2/muFp2) 
     x       - 500*CA**3*CF*Log(Q2/muFp2) 
     x       + 550*CA**4*CF*Log(Q2/muFp2) 
     x       + 120*CA**2*CF*game*Log(Q2/muFp2) 
     x       - 600*CA**4*CF*game*Log(Q2/muFp2) 
     x       + 120*CA**2*CF*Log(2.d0)*Log(Q2/muFp2) 
     x       - 600*CA**4*CF*Log(2.d0)*Log(Q2/muFp2) 
     x       - 200*CA*CF*Log(Q2/muR2) 
     x       + 220*CA**2*CF*Log(Q2/muR2) 
     x       + 1000*CA**3*CF*Log(Q2/muR2) 
     x       - 1100*CA**4*CF*Log(Q2/muR2))/(7040.*CA)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function c1_gg_qqb(q2,mu2,muf2,mur2,mud2,nf)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none

      integer nf

      double precision CA, CF, zeta2, game, pi,
     x     mu, muf, mur, mud, mufp2, q2,
     x     mu2, muf2, mur2, mud2

      parameter(CA = 3.d0)
      parameter(CF = 4.d0/3.d0)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   C0de starts here...
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   scales and constants 
      pi    = dacos(-1.d0)
c      mu2   = mu**2
c      muf2  = muf**2
c      mufp2 = mud**2
      mufp2 = mud2
c      mur2  = mur**2
      zeta2 = pi**2/6.d0
      game  = 0.5772156649d0
      game  = 0.d0 * game
cccccccccccccccccccccccccccc
      
      c1_gg_qqb = 
     x     (-42*CF + 63*CA**2*CF - 21*CA**4*CF - 18*CA**3*CF**2*game 
     x     + 72*CA**2*CF**3*game + 6*CA**2*CF*game**2 
     x     - 54*CA**4*CF*game**2 - 24*CA*CF**2*game**2 
     x     + 216*CA**3*CF**2*game**2 + 8*CF*Pi**2 - 12*CA**2*CF*Pi**2 
     x     + 4*CA**4*CF*Pi**2 + 6*CA**2*CF*zeta2 - 54*CA**4*CF*zeta2 
     x     - 24*CA*CF**2*zeta2 + 216*CA**3*CF**2*zeta2 - 42*CF*Log(2.d0) 
     x     + 15*CA**2*CF*Log(2.d0) + 3*CA**4*CF*Log(2.d0) 
     x     + 24*CF*game*Log(2.d0) 
     x     - 228*CA**2*CF*game*Log(2.d0) + 60*CA**4*CF*game*Log(2.d0) 
     x     + 72*CF*Log(2.d0)**2 - 36*CA**2*CF*Log(2.d0)**2 
     x     + 18*CA**4*CF*Log(2.d0)**2 - 44*CA**4*CF*Log(Q2/muF2) 
     x     + 176*CA**3*CF**2*Log(Q2/muF2) 
     x     - 288*CA**2*CF*game*Log(Q2/muF2) 
     x     + 336*CA**4*CF*game*Log(Q2/muF2) 
     x     - 768*CA**3*CF**2*game*Log(Q2/muF2) 
     x     + 8*CA**3*CF*nf*Log(Q2/muF2) 
     x     - 32*CA**2*CF**2*nf*Log(Q2/muF2) 
     x     - 18*CA**3*CF**2*Log(Q2/muFp2) 
     x     + 72*CA**2*CF**3*Log(Q2/muFp2) 
     x     + 24*CA**3*CF**2*game*Log(Q2/muFp2) 
     x     - 96*CA**2*CF**3*game*Log(Q2/muFp2) 
     x     + 24*CA**3*CF**2*Log(2.d0)*Log(Q2/muFp2) 
     x     - 96*CA**2*CF**3*Log(2.d0)*Log(Q2/muFp2) 
     x     + 44*CA**4*CF*Log(Q2/muR2) 
     x     - 176*CA**3*CF**2*Log(Q2/muR2) 
     x     - 8*CA**3*CF*nf*Log(Q2/muR2) 
     x     + 32*CA**2*CF**2*nf*Log(Q2/muR2))/(224.d0*CA)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
