cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function bstint_cuba(ndim_d, xx, ncomp_d, ff)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc
ccc   Wrapper around psintnlo_t function (used for direct
ccc   VEGAS calls) to be used with the CUBA library
ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none      

      include "int_set.inc"

      integer ndim_d, ncomp_d
      double precision xx(ndim_tmp+iqtint), ff(ncomp), wgt, 
     x     bstint
     
      external bstint
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   C0de begins here...
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   wgt set as dummy
      wgt     = 1.d0

      

c$$$      write(*,*)"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
c$$$      write(*,*)"ndim   = ", ndim
c$$$      write(*,*)"ncomp  = ", ncomp
c$$$      write(*,*)"xx     = ", xx
c$$$      write(*,*)"wgt    = ", wgt
c$$$      write(*,*)"iqtint = ", iqtint
c$$$      write(*,*)"ndimtmp= ", ndim_tmp

      ff(ncomp) = bstint(xx,wgt)

c$$$      write(*,*)"bstint    = ", bstint(xx,wgt)
c$$$      write(*,*)"ff(ncomp) = ", ff(ncomp)

      bstint_cuba = 0
      end
