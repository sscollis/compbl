        subroutine calc( etmp, gtmp, nstep, x, z)
c----------------------------------------------------------------------
c
c.... compute the integral of 1/rho and store in z
c
c----------------------------------------------------------------------
c       implicit double precision (a-h,o-z)
        dimension etmp(nstep), gtmp(nstep)
        dimension x(nstep), z(nstep)
c
        dimension eta(10000), g(10000), gs(10000)
        common /yes/ eta, g, gs, num
        external der
c----------------------------------------------------------------------
c
c.... spline g
c       
        num = nstep
        do i = 1, nstep
          g(i) = gtmp(i)
          eta(i) = etmp(i)
        end do
        call SPLINE (nstep, eta, g, gs)
c
c.... setup for the integration
c
        xmin = 0.0d0
        xmax = eta(nstep)

        z(1) = 0.0d0

        x(1) = 0.0d0
        h    = (xmax-xmin)/(nstep-1)    
c
c.... integrate
c
        do i = 1, nstep-1
          call der(x(i),z(i),dz)
          call rk4(z(i),dz,1,x(i),h,z(i+1),der)
c         write(50,"(8(1pe13.6,1x))") x(i), z(i), etmp(i), gtmp(i)
          x(i+1) = x(i) + h
        end do
        
        return
        end
        
c----------------------------------------------------------------------
        subroutine der(x,z,dz)
c----------------------------------------------------------------------
c
c.... given eta in x, compute 1/rho and put in dz using the spline
c
c----------------------------------------------------------------------
c       implicit double precision (a-h,o-z)
        dimension eta(10000), g(10000), gs(10000)
        common /yes/ eta, g, gs, num
c---------------------------------------------------------------------- 
        call SPEVAL (num,eta,g,gs,x,dz)
        dz = 1.0d0 / dz
        
        return
        end
