      program energy 05072012

      implicit none
      include 'par.for'
      character*15 nome, nome2
      integer i, j, k
      real*8 fc1, x, var(imax,jmax,kfour), en(imax), y(jmax), a, b, c,
     &       det, dendx(imax)

      call initval(var)
      do j = 1, jmax
        if (stf .ne. 1.d0) then
         y(j) = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)
        else
         y(j) = dy * dble(j-1)
        end if
      end do

      do k = 1, kfour
       do i = 1, imax
        fc1 = 0.d0
        do j = 2, jmax - 1, 2
        det = y(j-1)**2 * ( y(j) - y(j+1) )
     &      + y(j-1) * ( y(j+1)**2 - y(j)**2 )
     &      + y(j) * y(j+1) * ( y(j) - y(j+1) )
        a = (1.d0 / det) * ( var(i,j-1,k) * ( y(j)   - y(j+1) )
     &                     + var(i,j,k)   * ( y(j+1) - y(j-1) )
     &                     + var(i,j+1,k) * ( y(j-1) - y(j)   ) )
        b = (1.d0 / det) * ( var(i,j-1,k) * ( y(j+1)**2 - y(j)**2   )
     &                     + var(i,j,k)   * ( y(j-1)**2 - y(j+1)**2 )
     &                     - var(i,j+1,k) * ( y(j-1)**2 - y(j)**2   ) )
        c = (1.d0 / det) * ( var(i,j-1,k)
     &                     * y(j)   * y(j+1) * ( y(j)   - y(j+1) )
     &                     + var(i,j,k)
     &                     * y(j-1) * y(j+1) * ( y(j+1) - y(j-1) )
     &                     + var(i,j+1,k)
     &                     * y(j)   * y(j-1) * ( y(j-1) - y(j)   ) )

        fc1 = fc1 + (a / 3.d0) * ( y(j+1)**3 - y(j-1)**3 )
     &            + (b / 2.d0) * ( y(j+1)**2 - y(j-1)**2 )
     &            +  c         * ( y(j+1)    - y(j-1)    )
        end do
c        en(i) = log10(fc1)
       en(i) = log10( fc1 * dsqrt(Re) )
       end do
       write(*,*) ' The results are stored in the file enXX.dat'
 
       en(1) = en(2)
       call der(dendx, en)

       write(nome,'(a,i0.2,a)')'en',k-1,'.dat'
       write(nome2,'(a,i0.2)')'en',k-1
       open (1, file = nome ,status = 'unknown')
       write(1,*) 'VARIABLES="x","energy",der'
       write(1,*) 'ZONE T="',nome2,'", I=',imax - 1
       do i = 2, imax
         x = x0 + dble(i-1) * dx
         write(1,3) x, en(i), dendx(i)
       end do
       close (unit = 1)

      end do
    3 format(1x, 3e21.9e3)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine initval(var)

      ! subroutine to read the values
      implicit none
      include 'par.for'
      character*15 nome
      integer i, j, k, t, shift, my_rank
      real*8 var(imax,jmax,kfour)
      complex*16  ux(ptsx,jmax,kfour),  wx(ptsx,jmax,kfour),
     &            uy(ptsx,jmax,kfour),  wy(ptsx,jmax,kfour),
     &            uz(ptsx,jmax,kfour),  wz(ptsx,jmax,kfour),
     &            th(ptsx,jmax,kfour),
     &           uxt(imax,jmax,kfour), uyt(imax,jmax,kfour),
     &           uzt(imax,jmax,kfour)

      ! Disturbances variables data
      do my_rank = 0, np - 1
        write(nome,'(a,i0.2,a)')'data_',my_rank,'.bin'
        open(2, file = nome, form = 'unformatted')
        read(2) t
        if (my_form .eq. 2) then
          read(2) ux, uy, uz, wx, wy, wz, th
         else
          read(2) ux, uy, uz, wx, wy, wz
        end if
        close (unit = 2)
        shift = my_rank * (ptsx - inter - 1)
        do k = 1, kfour
          do j = 1, jmax
            do i = 1, ptsx
              uxt(i+shift,j,k) = ux(i,j,k)
              uyt(i+shift,j,k) = uy(i,j,k)
              uzt(i+shift,j,k) = uz(i,j,k)
            end do
          end do
        end do
      end do
      
      do k = 1, kfour
        do j = 1, jmax
          do i = 1, imax
c            var(i,j,k) = abs(uxt(i,j,k) * uxt(i,j,k))
            var(i,j,k) = 0.25d0 * ( abs(uxt(i,j,k) * uxt(i,j,k))
     &                            + abs(uyt(i,j,k) * uyt(i,j,k))
     &                            + abs(uzt(i,j,k) * uzt(i,j,k)) )
            if (k .eq. 1) var(i,j,k) = 0.5d0 *
     &                               ( abs(uxt(i,j,k) * uxt(i,j,k))
     &                               + abs(uzt(i,j,k) * uzt(i,j,k)) )
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine rhsx(fc, rhs)
      
      ! RHS for the first derivative calculation in x direction
      implicit none
      include 'par.for'
      integer i
      real*8 rhs(imax), fc(imax)

          rhs(1) = ( - 25.d0 * fc(1) + 48.d0 * fc(2)
     &               - 36.d0 * fc(3) + 16.d0 * fc(4)
     &               -  3.d0 * fc(5) ) / ( 12.d0 * dx )

          rhs(2) = ( - 406.d0 * fc(1) - 300.d0 * fc(2)
     &               + 760.d0 * fc(3) -  80.d0 * fc(4)
     &               +  30.d0 * fc(5) -   4.d0 * fc(6) ) /
     &               ( 120.d0 * dx )
          
          do i = 3, imax - 2
            rhs(i) = (           fc(i+2) - fc(i-2)
     &               + 28.d0 * ( fc(i+1) - fc(i-1) ) ) /
     &               ( 12.d0 * dx )
          end do
          
          rhs(imax-1) = ( - 406.d0 * fc(imax)
     &                    - 300.d0 * fc(imax-1)
     &                    + 760.d0 * fc(imax-2)
     &                    -  80.d0 * fc(imax-3)
     &                    +  30.d0 * fc(imax-4)
     &                    -   4.d0 * fc(imax-5) ) /
     &                  ( - 120.d0 * dx )
          
          rhs(imax) = ( - 74.d0 * fc(imax)
     &                  + 16.d0 * fc(imax-1)
     &                  + 72.d0 * fc(imax-2)
     &                  - 16.d0 * fc(imax-3)
     &                  +  2.d0 * fc(imax-4) ) /
     &                ( - 24.d0 * dx )

      return
      end
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine coef(a, b, c)

      ! mount the LHS of the matrix for the first derivative
      implicit none
      include 'par.for'
      integer i
      real*8 a(imax), b(imax), c(imax)

      a(1)      = 0.d0
      b(1)      = 1.d0
      c(1)      = 0.d0

      a(2)      = 1.d0
      b(2)      = 6.d0
      c(2)      = 2.d0

      do i = 3, imax - 2
        a(i)    = 1.d0
        b(i)    = 3.d0
        c(i)    = 1.d0
      end do

      a(imax-1) = 2.d0
      b(imax-1) = 6.d0
      c(imax-1) = 1.d0

      a(imax)   = 4.d0
      b(imax)   = 1.d0
      c(imax)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine trid(a, b, c, rhs)

      ! solves tridiagonal matrix for the derivatives in x direction
      implicit none
      include 'par.for'
      integer i, k
      real*8 a(imax), b(imax), c(imax), gam(imax), bet,
     &       rhs(imax), u(imax)

      bet  = b(1)
      u(1) = rhs(1) / bet
      do i = 2, imax
        gam(i) = c(i-1) / bet
        bet    = b(i) - a(i) * gam(i)
        u(i)   = ( rhs(i) - a(i) * u(i-1) ) / bet
      end do
      do i = imax - 1, 1, -1
        u(i) = u(i) - gam(i+1) * u(i+1)
      end do
      do i = 1, imax
        rhs(i) = u(i)
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine der(dendx, en)

      ! first derivatives calculation in x direction
      implicit none
      include 'par.for'
      real*8 a(imax), b(imax), c(imax), dendx(imax), en(imax)

      call coef(a, b, c)
      call rhsx(en, dendx)
      call trid(a, b, c, dendx)
      
      return
      end
