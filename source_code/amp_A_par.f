      program amp_A_par 05072012

      implicit none
      include 'par.for'
      character*15 nome, nome2
      integer i, j, k
      real*8 fc1, x, var(imax,jmax,kfour), en(imax), y(jmax), a, b, c,
     &       det

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
     &            + y(j-1) * ( y(j+1)**2 - y(j)**2 )
     &            + y(j) * y(j+1) * ( y(j) - y(j+1) )
            a = (1.d0 / det) * ( var(i,j-1,k) * ( y(j)   - y(j+1) )
     &                         + var(i,j,k)   * ( y(j+1) - y(j-1) )
     &                         + var(i,j+1,k) * ( y(j-1) - y(j)   ) )
            b = (1.d0 / det) *
     &          ( var(i,j-1,k) * ( y(j+1)**2 - y(j)**2   )
     &          + var(i,j,k)   * ( y(j-1)**2 - y(j+1)**2 )
     &          - var(i,j+1,k) * ( y(j-1)**2 - y(j)**2   ) )
            c = (1.d0 / det) * ( var(i,j-1,k)
     &                         * y(j)   * y(j+1) * ( y(j)   - y(j+1) )
     &                         + var(i,j,k)
     &                         * y(j-1) * y(j+1) * ( y(j+1) - y(j-1) )
     &                         + var(i,j+1,k)
     &                         * y(j)   * y(j-1) * ( y(j-1) - y(j)   ) )

            fc1 = fc1 + (a / 3.d0) * ( y(j+1)**3 - y(j-1)**3 )
     &                + (b / 2.d0) * ( y(j+1)**2 - y(j-1)**2 )
     &                +  c         * ( y(j+1)    - y(j-1)    )
          end do
          en(i) = fc1
c         en(i) = log10( fc1 * dsqrt(Re) )
        end do
        write(*,*) ' The results are stored in the file enXX.dat'

        write(nome,'(a,i0.2,a)')'amp',k-1,'.dat'
        write(nome2,'(a,i0.2)')'amp',k-1
        open (1, file = nome ,status = 'unknown')
        write(1,*) 'VARIABLES="x","amplitude"'
        write(1,*) 'ZONE T="',nome2,'", I=',imax - 1
        do i = 2, imax
          x = x0 + dble(i-1) * dx
          write(1,3) x, en(i)
        end do
        close (unit = 1)
      end do
    3 format(1x, 2d17.9)

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
     &            th(ptsx,jmax,kfour), uxt(imax,jmax,kfour),
     &           uyt(imax,jmax,kfour), uzt(imax,jmax,kfour)

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
            var(i,j,k) = abs(uxt(i,j,k) * uxt(i,j,k))
          end do
        end do
      end do

      return
      end
