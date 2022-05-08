c------------------------------------------------------------------------------
        subroutine derivs(x,y,dy)
c------------------------------------------------------------------------------
c
        include "common.h"
c
        dimension y(5), dy(5)
c
c.... Following Malik, compute temperature from Total enthalpy
c
        T  = y(5) + pt5*rMach**2*gamma1*(y(5) - y(2)**2)
        dT = y(4) + pt5*rMach**2*gamma1*(y(4) - two*y(2)*y(1))

        call getmat( Te, visce, dvisce)

        call getmat( T * Te, visc, dvisc)

        C    = one / T * visc/visce
        
        dC   = ( -one/T**2 * visc/visce + one/T * dvisc/visce * Te ) * dT

        Cinv = one / C
        a1  = C / Pr
        da1 = dc / Pr
        a2  = gamma1 * rMach**2 / (one + pt5*gamma1*rMach**2)*(one-one/Pr) * C
        da2 = a2 / C * dC
        
        dy(1) = Cinv * ( beta * ( y(2)**2 - T ) - y(3) * y(1) - 
     &                   dC * y(1) )
 
        dy(2) = y(1)

        dy(3) = y(2)

        dy(4) = -one/a1 * ( da1*y(4) + da2*y(2)*y(1) + a2*y(1)**2 +
     &                      a2*y(2)*dy(1) + y(3)*y(4) ) 

        dy(5) = y(4)

        return
        end

c------------------------------------------------------------------------------
        subroutine derivs_old(x,y,dy)
c------------------------------------------------------------------------------
c
        include "common.h"
c
        dimension y(5), dy(5)
c
c.... Following Reed
c
        call getmat( Te, visce, dvisce)
c
        call getmat( y(5) * Te, visc, dvisc)
c
c.... the righteous way
c
        C    =  one / y(5) * visc/visce
        
        dC   = -one / y(5)**2 * y(4) * visc/visce + 
     &         (one / y(5)) * (dvisc / visce) * (Te * y(4))
c
c.... the approximate way (see White p. 507)
c
        C    = y(5) ** (-pt33)
        dC   = -pt33 * y(5)**onept33 * y(4)

        Cinv = one / C

        dy(1) = Cinv * ( beta * ( y(2)**2 - y(5) ) - y(3) * y(1) - 
     &                   dC * y(1) )
 
        dy(2) = y(1)

        dy(3) = y(2)

        dy(4) =  -Cinv * (dc * y(4) + pr * y(3) * y(4) + 
     &                    pr * C * gamma1 * rMach**2 * y(1)**2 )

        dy(5) = y(4)

        return
        end

