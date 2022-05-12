c----------------------------------------------------------------------
c
c This file contains the common blocks and the data declaration needed 
c for the routines.  It is not necessary to modify the common blocks 
c other than in this file.
c
c----------------------------------------------------------------------
c
c	implicit double precision (a-h,o-z)
	implicit real (a-h,o-z)
c
	common /const / zero,   pt25,    pt33,    pt5,   pt66,   one,    
     &			onept5, two,     three,   four,  sixten, eps7,
     &                  pi,     onept25, onept33, eps8,
     &          	eps9,   eps10
c	
	common /mater / cv, pr, gamma, gamma1, gamma2, gamma3
c
        common /parm  / Rgas, cp, rmu, dmu, rKappa, dKappa, 
     &                  rMach, beta, Te
c
	common /run   / datmat(3), mattyp
c 
c---------------------------------------------------------------------- 
c 
c Scott Collis, Fall 1994
c 
c----------------------------------------------------------------------
