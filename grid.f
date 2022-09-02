        subroutine grid( nstep, zz, rho, u, v, T, rnu, xv, delta)
c----------------------------------------------------------------------
c
c.... Put the boundary layer solution on a grid
c
c----------------------------------------------------------------------
        include "common.h"
c
        dimension zz(nstep), rho(nstep), u(nstep), v(nstep), t(nstep)
        
c       parameter (NSTPMX=5000)
        parameter (NSTPMX=1000)
        dimension rhos(NSTPMX), us(NSTPMX), vs(NSTPMX), 
     &            ts(NSTPMX),   zs(NSTPMX)

c       parameter (Nx_max=512,Ny_max=512)
        parameter (Nx_max=256,Ny_max=256)
        real x(Nx_max,Ny_max), y(Nx_max,Ny_max), tmp
        real q(Nx_max,Ny_max,5)
        real Re, Ma, Pr, rgamma, cv

        character(1) ans
        logical temp
c----------------------------------------------------------------------
        if (nstep.gt.NSTPMX) stop 'Increase NSTPMX in grid'
c
c.... Determine the instability type
c
        write(*,"(/,'Enter [T]emporal or [S]patial ==> ',$)")
        read(*,*) ans

        if (ans.eq.'t' .or. ans.eq.'T') then
           temp = .true.
        else
           temp = .false.
        end if
c
c.... set the parameters
c
        Re     = one / rnu
        Ma     = rmach
        Pr     = 0.7
        rgamma = 1.4
        cv     = 716.5
c
c.... read the grid file
c
        open(unit=10,file='grid.dat',form='unformatted',
     &       action='read',status='old',err=1000)
        read(10,err=1010) nx, ny, nz
        write(*,*) nx, ny, nz
        read(10,err=1010) (((x(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &                    (((y(i,j), i=1,nx), j=1,ny), k = 1, nz),
     &                    (((   tmp, i=1,nx), j=1,ny), k = 1, nz)
        close(10)
c
c.... spline the data
c
        call SPLINE (nstep, zz, rho, rhos)
        call SPLINE (nstep, zz,   u,   us)
        call SPLINE (nstep, zz,   v,   vs)
        call SPLINE (nstep, zz,   t,   ts)
c
c.... put boundary layer solution on the grid
c
        xmin = x(1,1)
c
        do i = 1, nx
           if (x(i,1).ge.xmin) then
              if (temp) then
                 xx = xv
              else
                 xx = xv + x(i,1) - xmin
              end if
              do j = 1, ny
                 zi = y(i,j) / sqrt( two * rnu * xx )
                 call SPEVAL (nstep, zz, rho, rhos, zi, q1 )
                 call SPEVAL (nstep, zz,   u,   us, zi, q2 )
                 call SPEVAL (nstep, zz,   v,   vs, zi, q3 )
                 call SPEVAL (nstep, zz,   t,   ts, zi, q5 )
                 q(i,j,1) = q1
                 q(i,j,2) = q2
                 q(i,j,3) = sqrt(rnu / (two * xx)) * q3
                 q(i,j,4) = zero
                 q(i,j,5) = q5
              end do
           else
              q(i,j,1) = zero
              q(i,j,2) = zero
              q(i,j,3) = zero
              q(i,j,4) = zero
              q(i,j,5) = zero
           end if
        end do
c
c.... If temporal, then zero the v-velocity
c       
        if (temp) then
          do i = 1, nx
            do j = 1, ny
              q(i,j,3) = zero
            end do
          end do
        end if
c       
        open(unit=10,file='bl.dat',form='unformatted')
        write(10) nx, ny, nz
        write(10) zero, zero, zero, zero
c       write(10) sngl(zero), sngl(zero), sngl(zero), sngl(zero)
        write(10) (((q(i,j,1), i=1,nx), j=1,ny), k = 1, nz),
     &            (((q(i,j,2), i=1,nx), j=1,ny), k = 1, nz),
     &            (((q(i,j,3), i=1,nx), j=1,ny), k = 1, nz),
     &            (((q(i,j,4), i=1,nx), j=1,ny), k = 1, nz),
     &            (((q(i,j,5), i=1,nx), j=1,ny), k = 1, nz)
        close(10)
c
c.... write out an LNS mean.dat file
c
        open(unit=10,file='mean.dat',form='unformatted')
        write(*,*) 0, sngl(zero), nx, ny, nz, 5, Re, Ma, Pr, rgamma, cv
        write(10) 0, zero, nx, ny, nz, 5, Re, Ma, Pr, rgamma, cv
c       write(10) 0, sngl(zero), nx, ny, nz, 5, Re, Ma, Pr, rgamma, cv
        write(10) (((q(i,j,k), j=1,ny), i=1,nx), k=1,5)
        close(10)       
c
c.... compute the boundary layer displacement thickness
c
        open(10,file='delta.dat',status='unknown')
        do i = 1, nx
           if (x(i,1) .ge. xmin) then
              xx = xv + x(i,1) - xmin
              do j = 1, ny
                 rhol = q(i,j,1)
                 ul   = q(i,j,2)
                 if (j .eq. 1) then
                    dy2 = y(i,j+1) - y(i,j)
                    deltap = pt5 * dy2 * (one - rhol * ul)
                 else if ( j .eq. ny) then
                    dy1 = y(i,j) - y(i,j-1)
                    deltap = deltap + pt5 * dy1 * (one - rhol * ul)
                 else
                    dy1 = y(i,j) - y(i,j-1)
                    dy2 = y(i,j+1) - y(i,j)
                    deltap = deltap + pt5 * (dy1 + dy2) * (one - 
     &                       rhol * ul)
                 end if
              end do
              write(10,"(8(1pe13.6,1x))") x(i,1), deltap/rnu,
     &                                    delta*sqrt(two*xx/rnu)
           end if
        end do
        close(10)

        return
 1000   write(*,"('Error opening file ',a)") "grid.dat" 
        call exit(1)
 1010   write(*,"('Error reading from file ',a)") "grid.dat" 
        call exit(1)
        end
