* cuba.F
* Fortran chooser for the Cuba routines
* last modified 3 Feb 05 th

c$$$#define VEGAS 1
c$$$#define SUAVE 2
c$$$#define DIVONNE 3
c$$$#define CUHRE 4


	subroutine cuba(method, ndim, ncomp, integrand,
     &    integral, error, prob)

	implicit none

	integer ndim, ncomp, ldxgiven, method
	include "CUBA_init.inc"

	external integrand

	ldxgiven = ndim

cpm	if( method .eq. VEGAS ) then
	if( method .eq. 1 ) then
	   write(*,*)"ndim    = ", ndim
	   write(*,*)"ncomp   = ", ncomp
 	call vegas(ndim, ncomp, integrand, userdata,
     &    epsrel, epsabs, verbose, seed,
     &    mineval, maxeval, nstart, nincrease, nbatch,
     &    gridno, statefile,
     &    neval, fail, integral, error, prob)
c	  nregions = 1


cpm	else if( method .eq. SUAVE ) then
	else if( method .eq. 2 ) then

	call suave(ndim, ncomp, integrand, userdata,
     &    epsrel, epsabs, verbose + last, seed,
     &    mineval, maxeval, nnew, flatness,
     &    nregions, neval, fail, integral, error, prob)


cpm	else if( method .eq. DIVONNE ) then
	else if( method .eq. 3 ) then

	call divonne(ndim, ncomp, integrand, userdata,
     &    epsrel, epsabs, verbose, seed,
     &    mineval, maxeval, key1, key2, key3, maxpass,
     &    border, maxchisq, mindeviation,
     &    ngiven, ldxgiven, 0, nextra, 0,
     &    nregions, neval, fail, integral, error, prob)

cpm	else if( method .eq. CUHRE ) then
	else if( method .eq. 4 ) then
	call cuhre(ndim, ncomp, integrand, userdata,
     &    epsrel, epsabs, verbose + last,
     &    mineval, maxeval, key,
     &    nregions, neval, fail, integral, error, prob)

	else

	  print *, "invalid method ", method
	  return

	endif

	print *, "method   =", name(method)
	print *, "nregions =", nregions
	print *, "neval    =", neval
	print *, "fail     =", fail
	print '(G20.12," +- ",G20.12,"   p = ",F8.3)',
     &    (integral(ccuba), error(ccuba), prob(ccuba), ccuba = 1, ncomp)
	end

