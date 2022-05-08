        subroutine blkdta
c
c----------------------------------------------------------------------
c 
c  Almost all data statements are stored in this block data.
c
c note: hopefully this routine and 'common.h' are the only files 
c       needed to be modified for single or double precision; see 
c       common.h.
c
c----------------------------------------------------------------------
c
        include "common.h"
c
c----------------------------------------------------------------------
c
c.... useful constants
c.... common /const / zero,   pt25,   pt33,   pt5,    pt66,   one,    
c....                 onept5, two,    three,  four,   sixten, eps7,
c....                 pi,     onept25, onept33, eps8, eps9, eps10
c
        data    zero,                  pt25,                 
     &          pt33,                  pt5,                  
     &          pt66,                  one,                  
     &          onept5,                two,                  
     &          three,                 four,                 
     &          sixten,                eps7,
     &          pi,                    onept25,
     &          onept33,               eps8,
     &          eps9,                  eps10
     &  /       0.0000000000000000000d+0,   2.500000000000000000d-1, 
     &          3.3333333333333333333d-1,   5.000000000000000000d-1,   
     &          6.6666666666666666666d-1,   1.000000000000000000d+0,   
     &          1.5000000000000000000d+0,   2.000000000000000000d+0,   
     &          3.0000000000000000000d+0,   4.000000000000000000d+0,   
     &          1.6000000000000000000d+1,   1.000000000000000000d-7,
     &          3.1415926535897932385d+0,   1.250000000000000000d+0,
     &          1.3333333333333333333d+0,   1.000000000000000000d-8,
     &          1.0000000000000000000d-9,   1.000000000000000000d-10  /
c
c----------------------------------------------------------------------
c
c.... common /mater / cv,     pr,     gamma,  gamma1, gamma2, gamma3
c
        data    cv,                   pr,     
     &          gamma,                gamma1, 
     &          gamma2,               gamma3
     &  /       7.16500000000000d+2,  0.70000000000000d+0,
     &          1.40000000000000d+0,  4.00000000000000d-1,
     &          1.01192885125388d-1, -3.50000000000000d+0     /
c
c gamma1 = gamma - 1
c gamma2 = (gamma1) ** (1/gamma1)
c gamma3 = -gamma / gamma1
c
c----------------------------------------------------------------------
        return
        end
