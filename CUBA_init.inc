ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc	COMMON variables and parameters
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer verbose, last, seed, mineval, maxeval
	double precision epsrel, epsabs, userdata
	parameter (userdata = 0)
	parameter (epsrel = 1D-2)
	parameter (epsabs = 1D-12)
	parameter (verbose = 6)
	parameter (last = 4)
	parameter (seed = 0)
	parameter (mineval = 1000000)
	parameter (maxeval = 2000000)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc	VEGAS specific variables and parameters 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer nstart, nincrease, nbatch, gridno
	character*7 name(4)
	character*(*) statefile	
	parameter (nstart = 100000)
	parameter (nincrease = 25000)
	parameter (nbatch = 100000)
	parameter (gridno = 0)
	parameter (statefile = "")

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc	SUAVE specific variables and parameters
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer nnew
	double precision flatness
	parameter (nnew = 2000000)
	parameter (flatness = 0.25D0)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc	DIVONNE specific variables and parameters
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer key1, key2, key3, maxpass
	double precision border, maxchisq, mindeviation
	integer ngiven, nextra
	parameter (key1 = 47)
	parameter (key2 = 1)
	parameter (key3 = 1)
	parameter (maxpass = 5)
	parameter (border = 0D0)
	parameter (maxchisq = 10D0)
	parameter (mindeviation = .25D0)
	parameter (ngiven = 0)
c	parameter (ldxgiven = ndim)
	parameter (nextra = 0)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc	CUHRE specific variables and parameters
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer key
	parameter (key = 0)


	double precision integral(ncomp), error(ncomp), prob(ncomp)
c	double precision integral(1), error(1), prob(1)
	integer nregions, neval, fail, ccuba



