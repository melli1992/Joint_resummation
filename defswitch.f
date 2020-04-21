c mbeekveld changed this to have evolfs as well as output

      SUBROUTINE DEFSWITCH(
     X     CONF,FNAME,LFNAME,
     X     IEVOLSW,IEVFS, EVOLIS,EVOLFS,EVOLSW,EVOLMODE ,
     X     ISLL,ISNLL,ISLNON,ISLNON2,
     X     FSLL,FSNLL,FSLNON,ILNNON,
     X     FLAG 
     X     )


      LOGICAL EVOLIS, EVOLFS
      INTEGER CONF,FILENO,FILENO2,IEVFS,
     x     IEVOLSW, ILNNON, LFNAME, LFNAME2, EVOLMODE 
      DOUBLE PRECISION ISLL,ISNLL,ISLNON,ISLNON2,
     x     FSLL,FSNLL,FSLNON,FLAG,EVOLSW
      CHARACTER*30 FNAME, FNAME2
CCC      COMMON / TEST / CONF,HFUNC,FFUNC,EVOLSW       
CCC
CCC      CONF: 
CCC      1      = LL 
CCC      2      = NLL 
CCC      3      = LL + lnNoN 
CCC      4      = NLL + lnNoN 
CCC
CCC      ISLL   : turns on IS LL (on by default)
CCC      ISNLL  : turns on IS NLL 
CCC      ISLNON : turns on IS lnN/N 
CCC      ISLNON2: turns on lnN/N for fragmenting parton  
CCC      FSLL   : turns on FS LL (on by default)
CCC      FSNLL  : turns on FS NLL 
CCC      FSLNON : turns on FS lnN/N 
CCC      FLAG   : (obsolote for now - perhaps used later)
CCC      EVOLSW : (de)activates LMQ-terms in consts
CCC      
ccc   needs modification 
CCC      GAMQS  : turns GAMQQ term on/off
CCC      GAMGS  : turns GAMQG term on/off

      EVOLIS     = .FALSE. 
      EVOLFS     = .FALSE.       
      EVOLSW     = 1.d0
      ILNNON     = 1
      EVOLMODE   = 0
      IF(CONF.eq.0) THEN
         FNAME     = 'noresum'
         ISLL    = 0.d0
         ISNLL   = 0.d0
         ISLNON  = 0.d0
         ISLNON2 = 0.d0
         FSLL    = 0.d0
         FSNLL   = 0.d0
         FSLNON  = 0.d0
      ELSEIF(CONF.eq.1) THEN
         FNAME     = 'LL'
         ISLL    = 1.d0
         ISNLL   = 0.d0
         ISLNON  = 0.d0
         ISLNON2 = 0.d0
         FSLL    = 1.d0
         FSNLL   = 0.d0
         FSLNON  = 0.d0
      ELSEIF(CONF.eq.2) THEN
         FNAME       = 'NLL_EXP'
         ISLL    = 1.d0
         ISNLL   = 1.d0
         ISLNON  = 0.d0
         ISLNON2 = 0.d0
         FSLL    = 1.d0
         FSNLL   = 1.d0
         FSLNON  = 0.d0
      ELSEIF(CONF.eq.3) THEN
         FNAME       = 'NLL_EVOL'
	 EVOLIS = .TRUE.
	 EVOLMODE = 1
         ISLL    = 1.d0
         ISNLL   = 1.d0
         ISLNON  = 0.d0
         ISLNON2 = 0.d0
         FSLL    = 1.d0
         FSNLL   = 1.d0
         FSLNON  = 0.d0     
         EVOLSW     = 0.d0
      ELSEIF(CONF.eq.4) THEN
         FNAME       = 'NLL_EXP_ISLNNON'
	 EVOLIS = .FALSE.
	 EVOLMODE = 0
         ISLL    = 1.d0
         ISNLL   = 1.d0
         ISLNON  = 1.d0
         ISLNON2 = 0.d0
         FSLL    = 1.d0
         FSNLL   = 1.d0
         FSLNON  = 1.d0
      ELSEIF(CONF.eq.5) THEN
         FNAME       = 'NLL_EVOL_ISLNNON'
	 EVOLIS = .TRUE.
	 EVOLMODE = 2
         ISLL    = 1.d0
         ISNLL   = 1.d0
         ISLNON  = 0.d0
         ISLNON2 = 0.d0
         FSLL    = 1.d0
         FSNLL   = 1.d0
         FSLNON  = 1.d0
         EVOLSW     = 0.d0
      ELSEIF(CONF.eq.6) THEN
         FNAME       = 'NLL_EVOL_ISON'
	 EVOLIS = .TRUE.
	 EVOLMODE = 3
         ISLL    = 1.d0
         ISNLL   = 1.d0
         ISLNON  = 0.d0
         ISLNON2 = 0.d0
         FSLL    = 1.d0
         FSNLL   = 1.d0
         FSLNON  = 1.d0
         EVOLSW     = 0.d0
      ELSEIF(CONF.eq.7) THEN
         FNAME       = 'NLL_EVOL_ISFULL'
	 EVOLIS = .TRUE.
	 EVOLMODE = 4
         ISLL    = 1.d0
         ISNLL   = 1.d0
         ISLNON  = 0.d0
         ISLNON2 = 0.d0
         FSLL    = 1.d0
         FSNLL   = 1.d0
         FSLNON  = 1.d0
         EVOLSW  = 0.d0
      ENDIF
	IF (IEVFS.eq.1) THEN
        FSLNON = 0.d0
        ISLNON2 = 0.d0
	   	EVOLFS   = .FALSE.
	ELSEIF (IEVFS.eq.2) THEN
        FSLNON = 1.d0
        ISLNON2 = 1.d0
	   	EVOLFS   = .FALSE.
	ELSEIF (IEVFS.eq.3) THEN
        FSLNON = 1.d0
        ISLNON2 = 0.d0
	   	EVOLFS   = .TRUE.
	ELSEIF (IEVFS.eq.4) THEN
        FSLNON = 1.d0
        ISLNON2 = 0.d0
	   	EVOLFS   = .TRUE.
	ENDIF
        FLAG       = ISNLL * FSNLL  
	FNAME = TRIM(FNAME)

      call stringmod(FNAME, LFNAME)

       write(*,*) "======================================="       
       write(*,*) "||   DEFSWITCH     output            ||"
       write(*,*) "======================================="       
       write(*,*) "EVOLMODE     = ", EVOLMODE 
       write(*,*) "EVOLSW     = ", EVOLSW
       write(*,*) "ILNNON     = ", ILNNON
       write(*,*) "ISLL       = ", ISLL
       write(*,*) "ISNLL      = ", ISNLL
       write(*,*) "ISLNoN     = ", ISLNoN
       write(*,*) "ISLNoN2    = ", ISLNoN2
       write(*,*) "FSLL       = ", FSLL
       write(*,*) "FSNLL      = ", FSNLL
       write(*,*) "FSLNoN     = ", FSLNoN
       write(*,*) "EVOLIS     = ", EVOLIS
       write(*,*) "EVOLFS     = ", EVOLFS
       write(*,*) "FLAG       = ", FLAG
       write(*,*) "======================================="       


      RETURN
      END
