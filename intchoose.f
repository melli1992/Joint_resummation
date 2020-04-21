cccccccccccccccccccccccccccccccc
      subroutine intchoose(j)
cccccccccccccccccccccccccccccccc
ccc   choose LO or NLO and 
ccc   integration routine
      implicit none
      integer j, k
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   C0de starting here
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      k = -1
      j = -1

c$$$      do while( (k.gt.2) .or. (k.lt.1) )
c$$$         write(*,*)"##########################################"
c$$$         write(*,*)"||   Select order:                      ||"
c$$$         write(*,*)"||--------------------------------------||"
c$$$         write(*,*)"||   press 1 for LO                     ||"
c$$$         write(*,*)"||   press 2 for NLO                    ||"
c$$$         write(*,*)"##########################################"
c$$$         read(*,*) k
c$$$      enddo
c      do while((k.eq.1) .or. (k.eq.2) ) 
      do while( (j.gt.7) .or. (j.lt.1) ) 
c$$$         if( k.eq.1 ) then
c$$$            write(*,*)"##########################################"
c$$$            write(*,*)"||   Select integration routine:        ||"
c$$$            write(*,*)"||--------------------------------------||"
c$$$            write(*,*)"||   press 1 for VEGAS0                 ||"
c$$$            write(*,*)"||   press 2 for VEGAS7                 ||"
c$$$            write(*,*)"||   press 3 for VEGAN                  ||"            
c$$$            write(*,*)"||   press 4 for CUBA                   ||"
c$$$            write(*,*)"##########################################"
c$$$            read(*,*) j
c$$$         elseif ( k.eq.2 ) then
c$$$c            write(*,*)"not implemented at the moment..."
c$$$c            stop
c$$$            write(*,*)"##########################################"
c$$$            write(*,*)"||   Select integration routine:        ||"
c$$$            write(*,*)"||--------------------------------------||"
c$$$            write(*,*)"||   press 5 for VEGAS0                 ||"
c$$$            write(*,*)"||   press 6 for VEGAS7                 ||"
c$$$            write(*,*)"||   press 7 for VEGAN                  ||"            
c$$$            write(*,*)"##########################################"
c$$$            read(*,*) j
c$$$         endif
c         if(k.eq.1.or.k.eq.2) then
            write(*,*)"##########################################"
            write(*,*)"||   Select integration:                ||"
            write(*,*)"||--------------------------------------||"
            write(*,*)"||   press 1 for full                   ||"
            write(*,*)"||   press 2 for s-chan.                ||"
            write(*,*)"||   press 3 for t-chan.                ||"            
            write(*,*)"||   press 4 for t-chan. (w. same lept.)||"            
            write(*,*)"||   press 5 for full - s-chan          ||"            
            write(*,*)"||   press 6 for full - t-chan          ||"            
            write(*,*)"||   press 7 for full - s- and t-chan.  ||"            
            write(*,*)"##########################################" 
            read(*,*) j
c         endif
c      enddo
      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine integrate(func,j,k,res,sig,chi)
cccccccccccccccccccccccccccccccccccccccccccccccc
ccc   choose integration routine
      implicit none

c      include "myint.inc"

      integer j,k,proc, it, ndo, ndim, ncalls, niter

      double precision 
     x     res, sig, chi,
     x     avg,sigma,chi2,aux,
     x     avgi,sd,ti,tsi,
     x     avgv,sigv,chiv,auxv,
     x     xi(51,20),si,si2,swgt,schi,
     x     intveg

      parameter(ndim   = 5)
      parameter(ncalls = 100000)
      parameter(niter  = 10)

      external func, intveg
      common /result/ avg,sigma,chi2,aux  
      common /myvegasresult/ avgi,sd,ti,tsi    
      common /resultv/ avgv,sigv,chiv,auxv
      common /bveg2/xi,si,si2,swgt,schi,ndo,it
c      common /bveg2/xi(51,20),si,si2,swgt,schi,ndo,it
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   C0de starting here...
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      res = intveg(func,ndim,ncalls,niter)
      res = avgv
      sig = sigv
      chi = chiv

      return 
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
