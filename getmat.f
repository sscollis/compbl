c----------------------------------------------------------------------
        subroutine getmat( T, visc, dvisc)
c----------------------------------------------------------------------
c
        include "common.h"
c
c.... Constant viscosity
c
        if ( mattyp .eq. 0) then
          visc = datmat(1)
          dvisc = zero
        end if
c
c.... Sutherland's law
c
        if ( mattyp .eq. 1) then
          visc  = datmat(1) * (t)/datmat(2) * sqrt((t)/datmat(2)) *
     &            (datmat(2) + datmat(3)) / (t + datmat(3))
          dvisc  = (datmat(1) * (3.0d0*datmat(3)+t)*(datmat(3)+datmat(2)) *
     &             sqrt(t/datmat(2)))/(2.0d0 * (datmat(3)+t)**2 * datmat(2))
c
c.... Malik's "Sutherland's Law"
c
c         visc  = 1.458d-6 * sqrt(t) / (one + 110.4d0/t)
c         dvisc = 1.458d-6 * sqrt(t) * ( 3.0d0 * 110.4d0 + t ) / 
c     &           ( two * (110.4d0 + t)**2 )
        end if
c
        return
        end
