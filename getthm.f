c----------------------------------------------------------------------
        subroutine getthm( T, visc, dvisc)
c----------------------------------------------------------------------
c
        include "common.h"
c
c.... Constant viscosity
c
        if ( mattyp .eq. 0) then
          visc = rmu
          dvisc = zero
        end if
c
c.... Sutherland's law
c
        if ( mattyp .eq. 1) then
          visc = 1.46d-6 * T**1.5d0 / (T + 110.3d0)
          dvisc = 1.46d-6 * ( 1.5d0 * T**0.5 / (T + 110.3d0) - 
     &                        T**1.5d0 / (T + 110.3d0)**2 ) 
        end if
c
        return
        end
