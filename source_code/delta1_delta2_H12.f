      program integ 20140115

      ! this program reads the output of the simulation (spatial.dat)
      ! and calculates the values of displacement thickness (delta1),
      ! momentum thickness (delta2) and shape factor H12 = delta1 / delta2
      implicit none
      include 'par.for'
      character*15 nome
      integer i, j, k, t, my_rank, shift
      real*8 x, y(jmax), z, a, b, c, det, fc1, fc2, umax,
     &        delta(imax,kphys), delta1(imax,kphys),
     &       delta2(imax,kphys),    H12(imax,kphys),
     &          var1(imax,jmax,kphys), var2(imax,jmax,kphys),
     &           uxp(imax,jmax,kphys),  wxp(imax,jmax,kphys),
     &           uyp(imax,jmax,kphys),  wyp(imax,jmax,kphys),
     &           uzp(imax,jmax,kphys),  wzp(imax,jmax,kphys)
      complex*16 uxt(imax,jmax,kfour),  wxt(imax,jmax,kfour),
     &           uyt(imax,jmax,kfour),  wyt(imax,jmax,kfour),
     &           uzt(imax,jmax,kfour),  wzt(imax,jmax,kfour),
     &            ux(ptsx,jmax,kfour),   wx(ptsx,jmax,kfour),
     &            uy(ptsx,jmax,kfour),   wy(ptsx,jmax,kfour),
     &            uz(ptsx,jmax,kfour),   wz(ptsx,jmax,kfour)
      
      ! Disturbances variables data
      do my_rank = 0, np - 1
        write(nome,'(a,i0.2,a)')'data_',my_rank,'.bin'
        open(2, file = nome, form = 'unformatted')
        read(2) t
        read(2) ux, uy, uz, wx, wy, wz
        close (unit = 2)

        shift = my_rank * (ptsx - inter - 1)
        do k = 1, kfour
          do j = 1, jmax
            do i = 1, ptsx
              uxt(i+shift,j,k) = ux(i,j,k)
              uyt(i+shift,j,k) = uy(i,j,k)
              uzt(i+shift,j,k) = uz(i,j,k)
              wxt(i+shift,j,k) = wx(i,j,k)
              wyt(i+shift,j,k) = wy(i,j,k)
              wzt(i+shift,j,k) = wz(i,j,k)
            end do
          end do
        end do
      end do
      
      call f_to_p(uxp, uxt, imax)
      call f_to_p(uyp, uyt, imax)
      call f_to_p(uzp, uzt, imax)
      call f_to_p(wxp, wxt, imax)
      call f_to_p(wyp, wyt, imax)
      call f_to_p(wzp, wzt, imax)

      do j = 1, jmax
        if (stf .eq. 1.d0) then 
          ! without stretching
          y(j) = dy * dble(j-1)
         else
          ! with stretching
          y(j) = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)
        end if        
      end do

      do k = 1, kphys
        do i = 1, imax
          umax = uxp(i,jmax,k)
          delta(i,k) = 0.d0
          do j = 1, jmax
            var1(i,j,k) = 1.d0 - uxp(i,j,k) / umax
            var2(i,j,k) = uxp(i,j,k) / umax * (1.d0 - uxp(i,j,k) / umax)
            if (uxp(i,j,k)/umax .lt. 0.99d0) then
              delta(i,k) = y(j)
            end if
          end do
        end do

        do i = 1, imax
          fc1 = 0.d0
          fc2 = 0.d0
          do j = 2, jmax - 1, 2
            det = y(j-1)**2 * ( y(j) - y(j+1) )
     &          + y(j-1) * ( y(j+1)**2 - y(j)**2 )
     &          + y(j) * y(j+1) * ( y(j) - y(j+1) )
            a   = (1.d0 / det)
     &          * ( var1(i,j-1,k) * (y(j)   - y(j+1))
     &            + var1(i,j,k)   * (y(j+1) - y(j-1))
     &            + var1(i,j+1,k) * (y(j-1) - y(j)) )
            b   = (1.d0 / det)
     &          * ( var1(i,j-1,k) * (y(j+1)**2 - y(j)**2)
     &            + var1(i,j,k)   * (y(j-1)**2 - y(j+1)**2)
     &            - var1(i,j+1,k) * (y(j-1)**2 - y(j)**2) )
            c   = (1.d0 / det)
     &          * ( var1(i,j-1,k) * y(j)   * y(j+1) * (y(j)   - y(j+1))
     &            + var1(i,j,k)   * y(j-1) * y(j+1) * (y(j+1) - y(j-1))
     &            + var1(i,j+1,k) * y(j)   * y(j-1) * (y(j-1) - y(j)))
            fc1 = fc1 + (a / 3.d0) * (y(j+1)**3 - y(j-1)**3)
     &                + (b / 2.d0) * (y(j+1)**2 - y(j-1)**2)
     &                +  c         * (y(j+1)    - y(j-1))

            a   = (1.d0 / det) * ( var2(i,j-1,k) * (y(j)   - y(j+1))
     &                           + var2(i,j,k)   * (y(j+1) - y(j-1))
     &                           + var2(i,j+1,k) * (y(j-1) - y(j)))
            b   = (1.d0 / det)
     &          * ( var2(i,j-1,k) * (y(j+1)**2 - y(j)**2)
     &            + var2(i,j,k)   * (y(j-1)**2 - y(j+1)**2)
     &            - var2(i,j+1,k) * (y(j-1)**2 - y(j)**2) )
            c   = (1.d0 / det)
     &          * ( var2(i,j-1,k) * y(j)   * y(j+1) * (y(j)   - y(j+1))
     &            + var2(i,j,k)   * y(j-1) * y(j+1) * (y(j+1) - y(j-1))
     &            + var2(i,j+1,k) * y(j)   * y(j-1) * (y(j-1) - y(j) ) )
            fc2 = fc2 + (a / 3.d0) * (y(j+1)**3 - y(j-1)**3)
     &                + (b / 2.d0) * (y(j+1)**2 - y(j-1)**2)
     &                +  c         * (y(j+1)    - y(j-1))
          end do
          if (my_form .eq. 0) then
            delta1(i,k) = fc1 * L_1 * 1.d3 / dsqrt(fac_y)  ! milimeters
            delta2(i,k) = fc2 * L_1 * 1.d3 / dsqrt(fac_y)  ! milimeters
           else
            delta1(i,k) = fc1 / dsqrt(fac_y)
            delta2(i,k) = fc2 / dsqrt(fac_y)
          end if
          H12(i,k) = delta1(i,k) / delta2(i,k)
        end do
      end do

      open (1, file = 'H12.dat',status = 'unknown')
      write(1,*) 'VARIABLES="x","z","delta","delta1","delta2","H12"'
      write(1,*) 'ZONE I=',imax,', K=',kphys,', F=POINT'
      write(*,*) ' The results are stored in the file H12.dat'
      do k = 1, kphys
        z = dble(k-1) * dz
        do i = 1, imax
          x = x0 + dble(i-1) * dx
          write(1,5) x, z, delta(i,k), delta1(i,k), delta2(i,k),
     &               H12(i,k)
        end do
      end do
      close (unit = 1)
 
    5 format(1x, 6d25.17)

      return
      end
