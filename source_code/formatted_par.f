ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c                  subroutines transforms binary data in asc                  c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program formated 20140205

      implicit none
      include 'par.for'
      character*15 nome, name_file
      integer i, j, k, t, my_rank, shift, cont
      real*8 x(imax), y(jmax), z(kphys), z1, xx, yy, zz,
     &       uxbt(imax,jmax), uybt(imax,jmax), wzbt(imax,jmax),
     &       thb(imax,jmax), uxp(imax,jmax,kphys),  
     &       wxp(imax,jmax,kphys), uyp(imax,jmax,kphys),
     &       wyp(imax,jmax,kphys), uzp(imax,jmax,kphys),
     &       wzp(imax,jmax,kphys), thp(imax,jmax,kphys), 
     &       uxb(ptsx,jmax), uyb(ptsx,jmax), wzb(ptsx,jmax),
     &       Uxtt(imax,jmax,kphys), Uytt(imax, jmax, kphys),
     &       Uztt(imax,jmax,kphys), Wxtt(imax,jmax,kphys),
     &       Wytt(imax,jmax,kphys), Wztt(imax,jmax,kphys), 
     &       Thtt(imax,jmax,kphys),  thbt(imax,jmax)

      complex*16 uxt(imax,jmax,kfour),  wxt(imax,jmax,kfour),
     &           uyt(imax,jmax,kfour),  wyt(imax,jmax,kfour),
     &           uzt(imax,jmax,kfour),  wzt(imax,jmax,kfour),
     &           tht(imax,jmax,kfour),
     &           ux(ptsx,jmax,kfour),   wx(ptsx,jmax,kfour),
     &           uy(ptsx,jmax,kfour),   wy(ptsx,jmax,kfour),
     &           uz(ptsx,jmax,kfour),   wz(ptsx,jmax,kfour),
     &           th(ptsx,jmax,kfour)
     
      ! 0 -> Zero Curvature
      ! 1 -> Goertler vortices simulations
      ! 2 -> Goertler vortices simulations with heat transfer
          
      ! Disturbances variables data
        do cont = 0, 12, 4
          write(*,*)'cont=',cont
          do my_rank = 0, np - 1
            if (cont .eq. 0) then
              write(nome,'(a,i0.2,a)')'data_',my_rank,'.bin'
            else
              write(nome,'(a,i0.2,a,i0.2,a)')
     &             'data_',my_rank,'_',cont,'.bin'
            end if
            open(2, file = nome, form = 'unformatted')
            read(2) t
            if (my_form .eq. 1) then
              read(2) ux, uy, uz, wx, wy, wz
            else
              read(2) ux, uy, uz, wx, wy, wz, th
            end if
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
                  if (my_form .eq. 2) then
                    tht(i+shift,j,k) = th(i,j,k)
                  end if
                end do
              end do
            end do
          end do
          write(*,*)"Finishing load file .bin perturbation"
        
          ! reads the boundary layer profile
          open(1,file = 'baseflow2D/basens.bin', form = 'unformatted')
          read(1) uxbt, uybt, wzbt, thbt
          close(unit = 1)

          write(*,*)"Finishing load file .bin baseflow"

          ! fft from Fourier to physical space
          call f_to_p(uxp, uxt, imax)
          write(*,*)"Finishing Fourier Transform u"
          call f_to_p(uyp, uyt, imax)
          write(*,*)"Finishing Fourier Transform v"
          call f_to_p(uzp, uzt, imax)
          write(*,*)"Finishing Fourier Transform w"
          call f_to_p(wxp, wxt, imax)
          write(*,*)"Finishing Fourier Transform wx"
          call f_to_p(wyp, wyt, imax)
          write(*,*)"Finishing Fourier Transform wy"
          call f_to_p(wzp, wzt, imax)
          write(*,*)"Finishing Fourier Transform wz"
          if (my_form .eq. 2) then
            call f_to_p(thp, tht, imax)
            write(*,*)"Finishing Fourier Transform theta"
          end if

          do k = 1, kphys
            do j = 1, jmax
              do i = 1, imax         
                Uxtt(i,j,k) = uxbt(i,j) + uxp(i,j,k)
                Uytt(i,j,k) = uybt(i,j) + uyp(i,j,k)
                Uztt(i,j,k) = uzp(i,j,k) 
                Wxtt(i,j,k) = wxp(i,j,k)
                Wytt(i,j,k) = wyp(i,j,k)
                Wztt(i,j,k) = wzbt(i,j) + wzp(i,j,k)
                if (my_form .eq. 2) then
                  Thtt(i,j,k) = thbt(i,j) + thp(i,j,k)
                end if
              end do
            end do
          end do
          
          write(name_file, '(a, i0.2,a)')'spatial_',cont,'.bin'
          open(1, file = name_file, form = 'unformatted')
          if (my_form .eq. 1) then
            write(1) Uxtt, uyp, uzp 
          else
            write(1) Uxtt, Thtt, uyp, uzp
          end if
          close (unit = 1)

          write(*,*)"Finishing the program"
           
          if (fhz_gv .eq. 0.0) stop

      end do

    5 format(1x, 3d14.6, 3d17.9)
    6 format(1x, 3d14.6, 7d17.9)
    7 format(1x, 3d14.6, 6d17.9)
    8 format(1x, 3d14.6, 9d17.9)

      stop
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c              end subroutines transforms binary data in asc                  c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
