        program compbl
c------------------------------------------------------------------------------
c
c       Solves the compressible self-similar boundary layer equations
c
c       Using the fixed step-size Runge routine from Numerical Recipes to
c       integrate the boundary-layer equations.
c
c       I have coded the variable step routine, but using it messes up
c       the trapezoid integration that I do to compute the layer thickness
c
c       So, for now make sure you compile with USE_RUNGE defined
c
c       Author: S. Scott Collis
c
c       Revised: 4-10-96
c
c------------------------------------------------------------------------------
c
        include "common.h"
c
c.... These parameters are also in subroutine runge
c
        parameter (NMAX=5,NSTPMX=5000)
c
        dimension q1(NSTPMX), q2(NSTPMX),  q3(NSTPMX),  q4(NSTPMX)
        dimension d1(NSTPMX), d2(NSTPMX),  d3(NSTPMX),  d4(NSTPMX)
        dimension s1(NSTPMX), s2(NSTPMX),  s3(NSTPMX),  s4(NSTPMX)
        dimension x2(NSTPMX), zz(NSTPMX),  zs(NSTPMX)
c
#ifdef USE_RUNGE
        common /NR_RUNGE_PATH/ xx(NSTPMX), y(NMAX,NSTPMX)
#else
        common /NR_ODEINT_PATH/ kmax, kount, dxsav, xx(NSTPMX),
     &                          y(NMAX,NSTPMX)
#endif
        dimension ys(NMAX), key(2)
c
        dimension eta(NSTPMX)
c
        external derivs, NR_RKQC
c
c.... spline integration
c
        dimension ww(10000), g(10000), gs(10000)
        common /yes/ ww, g, gs, num
        external der

        character(1) ans
c------------------------------------------------------------------------------
        kmax  = NSTPMX
        dxsav = 0.0
c
c.... setup block data
c
        call blkdta
c
c.... Note: the Pr is set to 0.7 as default in blkdta
c.... get input
c
        call input
c
c.... compute machine presision
c
        call macprc(epsM)
c
c.... setup domain
c
        n     = 5
        xmin  = zero
        xmax  = 20.0d0
        nstep = 2000
        dx    = (xmax-xmin)/float(nstep)
c
c.... Adiabatic and non-slip boundary conditions
c
        if (.true.) then
          y(1,1) = 4.6659494418796d-01                  ! f'' (u')
          y(2,1) = zero                                 ! f'  (u)
          y(3,1) = zero                                 ! f
          y(4,1) = zero                                 ! g'
          y(5,1) = one                                  ! g   (T guess)
          key(1) = 1
          key(2) = 5
        endif
c
c.... Fixed wall temperature and non-slip boundary conditions
c
        if (.false.) then
          y(1,1) = pt5                                  ! f''
          y(2,1) = zero                                 ! f'
          y(3,1) = zero                                 ! f
          y(4,1) = one                                  ! g' (guess)
          y(5,1) = one / ( one + pt5*rMach**2*gamma1 )  ! g
          key(1) = 1
          key(2) = 4
        end if
c
c.... output some info
c
        write (*,5) nstep, xmin, xmax
 5      format (1p,/,'nstep = ',i5,' xmin = ',g10.3,' xmax =  ',g10.3,/)
c
c.... Begin Newton loop
c
        iter = 0
 100    continue
        iter = iter + 1
c
c.... save the starting values
c
        do i = 1, n
          ys(i) = y(i,1)
        end do
c
c.... integrate
c
#ifdef USE_RUNGE
        call NR_RUNGE(y(1,1),n,xmin,xmax,nstep,derivs)
#else
        call NR_ODEINT(y(1,1),n,xmin,xmax,eps8,1.0d0,0.0d0,nok,nbad,
     &                 derivs,NR_RKQC)
        nstep = kount-1
#endif
        write (*,20) iter, y(key(1),1), y(key(2),1), y(2,nstep+1),
     &               y(5,nstep+1)
 20     format(1p,i5,1x,4(e20.13,1x))
c
c.... Newton solve
c
        b1 = y(2,nstep+1)
        b2 = y(5,nstep+1)
c
c.... Perturb y(1,1)
c
        eps = sqrt(epsM)

        y(key(1),1) = (one + eps) * ys(key(1))
#ifdef USE_RUNGE
        call NR_RUNGE(y(1,1),n,xmin,xmax,nstep,derivs)
#else
        call NR_ODEINT(y(1,1),n,xmin,xmax,eps8,1.0d0,0.0d0,nok,nbad,
     &                derivs,NR_RKQC)
        nstep = kount-1
#endif
        a11 = (y(2,nstep+1)-b1)/(eps*ys(key(1)))
        a21 = (y(5,nstep+1)-b2)/(eps*ys(key(1)))
c
c.... Perturb y(5,1)
c
        y(key(1),1) = ys(key(1))
        y(key(2),1) = (one + eps) * ys(key(2))
#ifdef USE_RUNGE
        call NR_RUNGE(y(1,1),n,xmin,xmax,nstep,derivs)
#else
        call NR_ODEINT(y(1,1),n,xmin,xmax,eps8,1.0d0,0.0d0,nok,nbad,
     &                 derivs,NR_RKQC)
        nstep = kount-1
#endif
        a12 = (y(2,nstep+1)-b1)/(eps*ys(key(2)))
        a22 = (y(5,nstep+1)-b2)/(eps*ys(key(2)))
c
c.... solve for the correction
c
        det = a11*a22-a21*a12

        b1 = one - b1
        b2 = one - b2
        bnorm = sqrt(b1**2 + b2**2)

        dy1 = (b1*a22 - b2*a12) / det
        dy5 = (a11*b2 - a21*b1) / det

        y(key(1),1) = y(key(1),1) + dy1
        y(key(2),1) = y(key(2),1) + dy5

        if (bnorm .gt. eps8 .and. iter .le. 10) goto 100
c
c.... integrate one last time using the latest wall values
c
        iter = iter + 1
#ifdef USE_RUNGE
        call NR_RUNGE(y(1,1),n,xmin,xmax,nstep,derivs)
#else
        call NR_ODEINT(y(1,1),n,xmin,xmax,eps8,1.0d0,0.0d0,nok,nbad,
     &                 derivs,NR_RKQC)
        nstep = kount-1
#endif
c
        write (*,20) iter, y(key(1),1), y(key(2),1), y(2,nstep+1),
     &               y(5,nstep+1)
c
c------------------------------------------------------------------------------
c
c.... Converged
c
        write (*,30) bnorm
 30     format (/,'Converged to ',1pe10.3)
c
c.... compute the boundary layer thickness using Trapezoid integration
c
c     I have converted the integrals from dy to deta so that the results
c     are given in terms of y
c
        i = 1
        if (.false.) then
          T = y(5,i)
        else
          T = y(5,i) + pt5*rMach**2*gamma1*(y(5,i) - y(2,i)**2)
        end if
        rho = one / T
        u = y(2,i)
        delta = pt5 * dx * (T - u)
        theta = pt5 * dx * ( u * ( one - u ) )
c
        i = nstep+1
        if (.false.) then
          T = y(5,i)
        else
          T = y(5,i) + pt5*rMach**2*gamma1*(y(5,i) - y(2,i)**2)
        end if
        rho = one / T
        u = y(2,i)
        delta = delta + pt5 * dx * ( T - u)
        theta = theta + pt5 * dx * ( u * ( one - u ) )
c
        do i = 2, nstep
          if (.false.) then
            T = y(5,i)
          else
            T = y(5,i) + pt5*rMach**2*gamma1*(y(5,i) - y(2,i)**2)
          end if
          rho = one / T
          u = y(2,i)
          delta = delta + (T - u) * dx
          theta = theta + ( u * ( one - u ) ) * dx
        end do

        write (*,40) delta * sqrt(two)
 40     format (/,'delta = ', 1pe20.13)
        write (*,50) theta * sqrt(two)
 50     format ('theta = ', 1pe20.13)
c
c... compute derivatives and save the results
c
        open(unit=13,file='output.dat')
        write (13,*) 1, nstep+1, 5, xx(nstep+1)*sqrt(2.0)

        do i = 1, nstep+1
          T  = y(5,i) + pt5*rMach**2*gamma1*(y(5,i) - y(2,i)**2)
          dT = y(4,i) + pt5*rMach**2*gamma1*(y(4,i) - two*y(2,i)*y(1,i))
c
          rho  = one / T
          drho = -dT/(T**2)
c
          u  = y(2,i)
          du = y(1,i)
c
          v  = T * (xx(i) * y(2,i) - y(3,i))
          dv = dT * (xx(i) * y(2,i) - y(3,i)) + T * (xx(i) * y(1,i))
c
          q1(i) = rho
          q2(i) = u
          q3(i) = v
          q4(i) = T
c
          d1(i) = drho
          d2(i) = du
          d3(i) = dv
          d4(i) = dT
c
          write (13,10) xx(i)*sqrt(2.0), q1(i), q2(i), q3(i), zero,
     &                  q4(i)
c
        end do
        close(13)
c
c.... compute the second derivatives using finite differences
c
        i = 1
        s1(i) = (d1(i+1) - d1(i)) / dx
        s2(i) = (d2(i+1) - d2(i)) / dx
        s3(i) = (d3(i+1) - d3(i)) / dx
        s4(i) = (d4(i+1) - d4(i)) / dx

        do i = 2, nstep
          s1(i) = pt5 * (d1(i+1) - d1(i-1)) / dx
          s2(i) = pt5 * (d2(i+1) - d2(i-1)) / dx
          s3(i) = pt5 * (d3(i+1) - d3(i-1)) / dx
          s4(i) = pt5 * (d4(i+1) - d4(i-1)) / dx
        end do

        i = nstep + 1
        s1(i) = (d1(i) - d1(i-1)) / dx
        s2(i) = (d2(i) - d2(i-1)) / dx
        s3(i) = (d3(i) - d3(i-1)) / dx
        s4(i) = (d4(i) - d4(i-1)) / dx
c
c------------------------------------------------------------------------------
c.... need to compute y as a function of eta
c
c       xx = original eta
c       x2 = eta at the integration points
c       zz = the integral of 1/rho at the integration points
c       y  = sqrt(2*ny*x/U) * zz
c
        call calc ( xx, q1, nstep+1, x2, zz )
c
c.... spline the integral
c
        call SPLINE (nstep+1,x2,zz,zs)
c
c.... compute delta_999, (with linear interpolation)
c
        do i = 1, nstep+1
          if (y(2,i) .gt. 0.999d0) then
            d999 = xx(i-1) + (0.999d0 - y(2,i-1))*(xx(i)-xx(i-1)) /
     &                       (y(2,i)-y(2,i-1))
            goto 200
          end if
        end do
 200    continue
c
c.... convert eta_999 to delta_999
c
c       SSC 200212:  changed z2 to zs
c
c       call SPEVAL (nstep+1,x2,zz,z2,d999,tmp)
        call SPEVAL (nstep+1,x2,zz,zs,d999,tmp)
        d999 = tmp
c
        write (*,60) d999 * sqrt(two)
 60     format (/,'d_999 = ', 1pe20.13)
        write (*,70) d999 / delta
 70     format ('d_999/delta = ',1pe20.13)
c
c.... now write out files in physical coordinates
c
        write(*,"(/,'Enter Re_delta ==> ',$)")
        read(*,*) Red
        rnu = datmat(1)
c
c.... Re_delta = C_delta * sqrt(2) * sqrt(Re_x)
c
        xv = (Red/(delta*sqrt(two)))**2 * rnu

        write (*,*) 'Re_delta = ',Red,' Nu = ',rnu
        write (*,*) 'Virtual origin = ',xv

        call SPEVAL (nstep+1,x2,zz,zs,xx(nstep+1),zi)
        yi = sqrt( two * rnu * xv) * zi

        open(55,file='header.dat')
        write(55,"('{')")
        write(55,"('  ',a,'Ny',a,' : ',i5)") '"', '"', nstep+1
        write(55,"('  ',a,'Ndof',a,' : ',i5)") '"', '"', n
        write(55,"('  ',a,'Ymax',a,' : ',1pe20.13)") '"', '"', yi
        write(55,"('}')")
        close(55)

        open(10,file='profile.dat')
        open(11,file='first.dat')
        open(12,file='second.dat')

        do i = 1, nstep + 1
          call SPEVAL (nstep+1,x2,zz,zs,xx(i),zi)
          yi = sqrt( two * rnu * xv) * zi
c
c.... NOTE that xx(:) has the z integral in it
c
          xx(i) = zi

          dedy   = sqrt( one / (two * rnu * xv) ) * q1(i)
          d2edy2 = sqrt( one / (two * rnu * xv) ) * d1(i) * dedy

          fact = sqrt( pt5 * rnu / xv )

          write (10,10) yi, q1(i), q2(i), fact * q3(i), zero, q4(i)
          write (11,10) yi, dedy*d1(i), dedy*d2(i), fact*dedy*d3(i),
     &                      zero, dedy*d4(i)
          write (12,10) yi, dedy**2 * s1(i) + d2edy2 * d1(i),
     &                      dedy**2 * s2(i) + d2edy2 * d2(i),
     &                      fact * ( dedy**2 * s3(i) + d2edy2 * d3(i) ),
     &                      zero,
     &                      dedy**2 * s4(i) + d2edy2 * d4(i)
        end do
        close(10)
        close(11)
        close(12)
c
c.... compute the displacement thickness (this is a check)
c
        do j = 1, nstep + 1
          u = y(2,j)
          T = y(5,j) + pt5*rMach**2*gamma1*(y(5,j) - y(2,j)**2)
          rho = one / T
          if (j .eq. 1) then
            dy2 = xx(j+1) - xx(j)
            deltap = pt5 * dy2 * (one - rho * u)
          else if ( j .eq. nstep+1) then
            dy1 = xx(j) - xx(j-1)
            deltap = deltap + pt5 * dy1 * (one - rho * u)
          else
            dy1 = xx(j) - xx(j-1)
            dy2 = xx(j+1) - xx(j)
            deltap = deltap + pt5 * (dy1 + dy2) * (one - rho * u)
          end if
        end do
        write (*,*) 'Delta = ', sqrt( two * rnu * xv) * deltap
c
c.... write out the field on a grid
c
        write(*,"(/,'Write on a grid? ==> ',$)")
        read(*,"(a1)") ans
        if (ans.eq.'y' .or. ans.eq.'Y') then
          write (*,*) 'In Grid...'
          call grid ( nstep+1, xx, q1, q2, q3, q4, rnu, xv, delta )
        end if
c
c.... format statements
c
 10     format(1p,7(e20.13,1x))
c
        stop
        end
