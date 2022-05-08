      program pro
      
      parameter (ny=100)
      parameter (pi=3.1415926535897932385)

      dy = 1.0 / real(ny-1)
      
      do j = 1, ny
        y = real(j-1) * dy
        write(10,*) y, y, 2.0*y - y**2, 1.5*y - 0.5*y**3, sin(pi*0.5*y)
        write(20,*) y**(1.0/7.0), y
      end do
      
      stop
      end
