ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c                  Subroutines transforms binary data in asc                  c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program formated 12032009

      implicit none
      include 'par.for'
      character*50 nome
      integer i, j, k, t, my_rank, shift
      real*8 x, y, z,
     &           uxp(imax,jmax,kphys), wxp(imax,jmax,kphys),
     &           uyp(imax,jmax,kphys), wyp(imax,jmax,kphys),
     &           uzp(imax,jmax,kphys), wzp(imax,jmax,kphys),
     &           thp(imax,jmax,kphys),
     &           uxbt(imax,jmax), uybt(imax,jmax), wzbt(imax,jmax),
     &           uxb(ptsx,jmax), uyb(ptsx,jmax), wzb(ptsx,jmax)
      complex*16 uxt(imax,jmax,kfour), wxt(imax,jmax,kfour),
     &           uyt(imax,jmax,kfour), wyt(imax,jmax,kfour),
     &           uzt(imax,jmax,kfour), wzt(imax,jmax,kfour),
     &           tht(imax,jmax,kfour),
     &            ux(ptsx,jmax,kfour),  wx(ptsx,jmax,kfour),
     &            uy(ptsx,jmax,kfour),  wy(ptsx,jmax,kfour),
     &            uz(ptsx,jmax,kfour),  wz(ptsx,jmax,kfour),
     &           ux2(ptsx,jmax,kfour), wx2(ptsx,jmax,kfour),
     &           uy2(ptsx,jmax,kfour), wy2(ptsx,jmax,kfour),
     &           uz2(ptsx,jmax,kfour), wz2(ptsx,jmax,kfour),
     &            th(ptsx,jmax,kfour)

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
        write(nome,'(a,i0.2,a)')
     &    '../cobem_cs001_semTS/data_',my_rank,'.bin'
        open(2, file = nome, form = 'unformatted')
        read(2) t
        read(2) ux2, uy2, uz2, wx2, wy2, wz2
        close (unit = 2)
        do k = 1, kfour
          do j = 1, jmax
            do i = 1, ptsx
              ux(i,j,k) = ux(i,j,k) - ux2(i,j,k)
              uy(i,j,k) = uy(i,j,k) - uy2(i,j,k)
              uz(i,j,k) = uz(i,j,k) - uz2(i,j,k)
              wx(i,j,k) = wx(i,j,k) - wx2(i,j,k)
              wy(i,j,k) = wy(i,j,k) - wy2(i,j,k)
              wz(i,j,k) = wz(i,j,k) - wz2(i,j,k)
            end do
          end do
        end do
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
              tht(i+shift,j,k) = th(i,j,k)
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
      call f_to_p(thp, tht, imax)

      if (my_form .eq. 2) then
        ! writes data to spacial space to be open by tecplot
        open (3, file = 'spatial.dat',status = 'unknown')
        write(3,*) 'VARIABLES="x","y","z","u","v","w",
     &             "wx","wy","wz","theta"'
        write(3,*) 'ZONE I=',imax,', J=',jmax,', K=',kphys+1,', F=POINT'

        z = -1.d0 * dz
        do j = 1, jmax
          if (stf .eq. 1.d0) then
            ! without stretching
            y = dble(j-1) * dy
           else
            ! with stretching
            y = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)
          end if
          do i = 1, imax
            x = x0 + dble(i-1) * dx
            write(3,6) x, y, z, uxp(i,j,kphys), uyp(i,j,kphys),
     &                 uzp(i,j,kphys), wxp(i,j,kphys), wyp(i,j,kphys),
     &                 wzp(i,j,kphys), thp(i,j,kphys)
          end do
        end do

        do k = 1, kphys
          z = dble(k-1) * dz
          do j = 1, jmax
            if (stf .eq. 1.d0) then
              ! without stretching
              y = dble(j-1) * dy
             else
              ! with stretching
              y = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)
            end if
            do i = 1, imax
              x = x0 + dble(i-1) * dx
            write(3,6) x, y, z, uxp(i,j,k), uyp(i,j,k), uzp(i,j,k),
     &                 wxp(i,j,k), wyp(i,j,k), wzp(i,j,k), thp(i,j,k)
            end do
          end do
        end do
        close (unit = 3)

       else ! my_form = 0 or my_form = 1 or my_form = 4
        ! writes data to spacial space to be open by tecplot
        open (3, file = 'spatial.dat',status = 'unknown')
        write(3,*) 'VARIABLES="x","y","z","u","v","w","wx","wy","wz"'
        write(3,*) 'ZONE I=',imax,', J=',jmax,', K=',kphys+1,',F=POINT'

        z = -1.d0 * dz
        do j = 1, jmax
          if (stf .eq. 1.d0) then
            ! without stretching
            y = dble(j-1) * dy
           else
            ! with stretching
            y = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)  
          end if        
          do i = 1, imax
            x = x0 + dble(i-1) * dx
            write(3,5)x, y, z, uxp(i,j,kphys), uyp(i,j,kphys), 
     &                uzp(i,j,kphys), wxp(i,j,kphys), wyp(i,j,kphys),
     &                wzp(i,j,kphys)
          end do
        end do

        do k = 1, kphys
          z = dble(k-1) * dz
          do j = 1, jmax
            if (stf .eq. 1.d0) then
              ! without stretching
              y = dble(j-1) * dy
             else
              ! with stretching
              y = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)  
            end if        
            do i = 1, imax
              x = x0 + dble(i-1) * dx
            write(3,5) x, y, z, uxp(i,j,k), uyp(i,j,k), uzp(i,j,k),
     &                 wxp(i,j,k), wyp(i,j,k), wzp(i,j,k)
            end do
          end do
        end do
        close (unit = 3)
      end if

    5 format(1x, 3d14.6, 6d17.9)
    6 format(1x, 3d14.6, 7d17.9)

      stop
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c              end subroutines transforms binary data in asc                  c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
