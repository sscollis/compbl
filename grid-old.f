	subroutine grid_temp( yy, rho, u, v, T, nstep)
c----------------------------------------------------------------------
c
c.... Put the boundary layer solution on a grid
c
c----------------------------------------------------------------------
	include "common.h"
c
	dimension yy(nstep), rho(nstep), u(nstep), v(nstep), t(nstep)
	
        parameter (NSTPMX=5000)
	dimension rhos(NSTPMX), us(NSTPMX), vs(NSTPMX), 
     &            ts(NSTPMX), zs(NSTPMX)

	parameter (Nymax)
	dimension y(Nymax)
c----------------------------------------------------------------------
c
c.... spline the data
c
	call SPLINE (nstep, yy, rho, rhos)
	call SPLINE (nstep, yy,   u,   us)
	call SPLINE (nstep, yy,   v,   vs)
	call SPLINE (nstep, yy,   t,   ts)
c
	xmach = rMach
	write (*,"('Enter length: ',$)") 
	read(*,*) tx
	write (*,"('Enter Number of elements in x: ',$)") 
	read(*,*) nx
c
c.... Note that these are regular T and P, not stagnation
c
        p0 = one / gamma / xmach**2
        T0 = p0 / Rgas
c
	nxp1 = nx + 1
c
c.... get uniform spacing in x
c
        dx   = tx / float(nx)
        xmin = 0.0d0
c
c.... read the y grid from a file
c
	open(33,file='coord.dat')
c
	i = 1
 11	continue
	  read(33,*,end=22) idum, idum, x, y(i)
	  i = i + 1
	  if (i.ge.Nymax) stop 'Error in grid, increase Nymax'
	  goto 11
 22     continue
        ny   = i - 2
        nyp1 = i - 1
        write(*,*) 'Ny = ',ny,'  Ny+1 = ',nyp1
c
c.... put boundary layer on the grid
c
	write (30,*) 'Created by COMBL'
	write (30,*) 0, 0, gamma, cv
	n = 0
	do i = 1, nxp1
	  x  = xmin + float(i-1)*dx
	  do j = 1, nyp1
	    n = n + 1
	    call SPEVAL (nstep, yy, rho, rhos, y(j), rhoi)
	    call SPEVAL (nstep, yy,   u,   us, y(j),   ui)
	    call SPEVAL (nstep, yy,   v,   vs, y(j),   vi)
	    call SPEVAL (nstep, yy,   t,   ts, y(j),   ti)

c	    write (30,10) n, rhoi, ui, sqrt(rnu/(two*xv)) * vi, T0*ti

	    write (30,10) n, rhoi, ui, 0.0d0, T0*ti

 10         format (1p,i5,2x,4(e20.13,1x))
	  end do
	end do
	
	return
	end
