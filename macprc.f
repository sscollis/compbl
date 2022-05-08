c----------------------------------------------------------------------
        subroutine macprc (epsM)
c
c----------------------------------------------------------------------
c
c This subroutine computes the machine precision.
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        var1 = one
        var2 = one
c
        do 100 i = 1, 1000
          var1 = var1 / two
          if (var1 + var2 .eq. var2) then
            epsM = two * var1
c           write (*,10) epsM
10          format('Machine Precision:  epsM = ',1pe20.13)
            return
          endif
100     continue
c
c.... end
c
        return
        end
