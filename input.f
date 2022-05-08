c----------------------------------------------------------------------
        subroutine input()
c----------------------------------------------------------------------
c
        include "common.h"
c
c.... flow parameters (Assume constant Cp, Cv, Rgas, gamma, Pr)
c
        Rgas = gamma1 * cv
        cp = Rgas + cv
c
c.... Constant mu or Sutherland's law
c
 1      write (*,5)
 5      format(/,'(0) for constant Mu, (1) for Sutherland ==> ',$)
        read (*,*) mattyp
        if (mattyp .ne. 0 .and. mattyp .ne. 1) goto 1
c
c.... Constant mu
c
        if (mattyp .eq. 0) then
          write (*,10) 
 10       format(/,'Enter the viscosity ==> ',$)
          read (*,*) datmat(1)
        else
          write (*,15)
 15       format(/,'Freestream stagnation temperature (T)) ==> ',$)
          read (*,*) T0
c
c.... changed to match N. Adams 1-24-95
c
          datmat(1) = 1.715336725523065e-05     ! 1.716e-5
          datmat(2) = 273.0
          datmat(3) = 110.4                     ! 111.0 
        end if
c
c.... Freestream Mach number
c
        write (*,20)
 20     format(/,'Enter the Mach number ==> ',$)
        read (*,*) rMach 
        Te = T0 / (one + gamma1 * pt5 * rMach**2)
c
c.... Falkner-Skan type parameter (beta = zero for a flate plate)
c
        write (*,"(/,'Enter beta ==> ',$)")
        read(*,*) beta
c
c       beta = zero
c
c.... Echo input
c
        write (*,30) Pr, Cv, Rgas
        write (*,40) gamma, rMach, beta
        write (*,60) Te, T0

 30     format (1p,/,'Pr =    ',e10.3,'   Cv =   ',e10.3,'   Rgas = ',e10.3)
 40     format (1p,  'Gamma = ',e10.3,'   Mach = ',e10.3,'   Beta = ',e10.3)
 60     format (1p,  'Te   =  ',e10.3,'   T0 =   ',e10.3)

        return
        end

