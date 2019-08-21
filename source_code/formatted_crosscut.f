ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c                  subroutines transforms binary data in asc                  c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program formated 20140205

      implicit none
      include 'par.for'
      include 'comm.coef'
      include 'comm.fourier'
      character*30 nome
      integer i, j, k, t, my_rank, shift
      real*8 x(imax), y(jmax), z(kphys), z1, xx, yy, zz,
     &       uxbt(imax,jmax), uybt(imax,jmax), wzbt(imax,jmax),
     &        uxb(ptsx,jmax), uyb(ptsx,jmax), wzb(ptsx,jmax),
     &             uxp(imax,jmax,kphys),  wxp(imax,jmax,kphys),
     &             uyp(imax,jmax,kphys),  wyp(imax,jmax,kphys),
     &             uzp(imax,jmax,kphys),  wzp(imax,jmax,kphys),
     &          duxpdy(imax,jmax,kphys),     duxbdy(ptsx,jmax),
     &               duxbtdy(imax,jmax), duxpdz(imax,jmax,kphys)  
      complex*16   uxt(imax,jmax,kfour),  wxt(imax,jmax,kfour),
     &             uyt(imax,jmax,kfour),  wyt(imax,jmax,kfour),
     &             uzt(imax,jmax,kfour),  wzt(imax,jmax,kfour),
     &              ux(ptsx,jmax,kfour),   wx(ptsx,jmax,kfour),
     &              uy(ptsx,jmax,kfour),   wy(ptsx,jmax,kfour),
     &              uz(ptsx,jmax,kfour),   wz(ptsx,jmax,kfour),
     &           duxdy(ptsx,jmax,kfour), duxtdy(imax,jmax,kfour), 
     &           duxtdz(imax,jmax,kfour) 

      select case (my_form)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! CASE 1 - Goertler vortices simulations !!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(1)
        open(1, file = 'pre_processing/coefs.bin', form = 'unformatted')
        read(1) fp_fd_coef_e
        read(1) sp_fd_coef
        read(1) cp_fd_coef
        read(1) pp_fd_coef
        read(1) lp_fd_coef
        read(1) fp_sd_coef
        read(1) sp_sd_coef
        read(1) cp_sd_coef
        read(1) pp_sd_coef
        read(1) lp_sd_coef
        read(1) sp_poi_coef
        read(1) cp_poi_coef
        read(1) pp_poi_coef
        read(1) lp_poi_coef
        read(1) w_at_w_coef
        read(1) dwydy_coef
        read(1) ! integration in the y-direction, used in baseflow2D
        read(1) ! integration in the y-direction, used in baseflow2D
        read(1) ! integration in the y-direction, used in baseflow2D
        read(1) ! integration in the y-direction, used in baseflow2D
        close(unit = 1)

          call derivs_k
          call create_ctes

          ! Disturbances variables data
          do my_rank = 0, np - 1
            write(nome,'(a,i0.2,a)')'data_',my_rank,'.bin'
            open(2, file = nome, form = 'unformatted')
            read(2) t
            read(2) ux, uy, uz, wx, wy, wz
            close (unit = 2)
            call dery(duxdy, ux)
            shift = my_rank * (ptsx - inter - 1)
            do k = 1, kfour
              do j = 1, jmax
                do i = 1, ptsx
                  uxt(i+shift,j,k) = ux(i,j,k)
                  duxtdy(i+shift,j,k) = duxdy(i,j,k)
                  uyt(i+shift,j,k) = uy(i,j,k)
                  uzt(i+shift,j,k) = uz(i,j,k)
                  wxt(i+shift,j,k) = wx(i,j,k)
                  wyt(i+shift,j,k) = wy(i,j,k)
                  wzt(i+shift,j,k) = wz(i,j,k)
                end do
              end do
            end do
          end do

          do k = 1, kfour
            do j = 1, jmax
              do i = 1, imax
                duxtdz(i,j,k) = v_kb(k)*uxt(i,j,k)
              end do
             end do
           end do
          

          do my_rank = 0, np - 1
            write(nome,'(a,i0.2,a)')'baseflow2D/based_',my_rank,'.bin'
            open(2, file = nome, form = 'unformatted')
            read(2) uxb, uyb, wzb
            close (unit = 2)
            call dery2d(duxbdy, uxb)
            shift = my_rank * (ptsx - inter - 1)
            do j = 1, jmax
              do i = 1, ptsx
                uxbt(i+shift,j) = uxb(i,j)
                duxbtdy(i+shift,j) = duxbdy(i,j)
                uybt(i+shift,j) = uyb(i,j)
                wzbt(i+shift,j) = wzb(i,j)
!               uxbt(i+shift,j) = 0d0
!               uybt(i+shift,j) = 0d0
!               wzbt(i+shift,j) = 0d0
              end do
            end do
          end do

          ! fft from Fourier to physical space
          call f_to_p(uxp, uxt, imax)
          call f_to_p(duxpdy, duxtdy, imax)
          call f_to_p(duxpdz, duxtdz, imax)
          call f_to_p(uyp, uyt, imax)
          call f_to_p(uzp, uzt, imax)
          call f_to_p(wxp, wxt, imax)
          call f_to_p(wyp, wyt, imax)
          call f_to_p(wzp, wzt, imax)

!         ! writes data to spacial space to be open by tecplot
!         open (3, file = 'spatial.dat',status = 'unknown')
!         write(3,*) 'VARIABLES="x","y","z","u","v","w","wx","wy","wz",
!    & "duxdy"'
!         write(3,*) 'ZONE I=',imax,', J=',jmax,', K=',kphys+1,
!    &               ',F=POINT'
!         
!         ! coordinate of z-direction
!         zz = - 1.d0 * dz
!         do j = 1, jmax
!           ! coordinate of y-direction
!           if (stf .eq. 1.d0) then
!             ! without stretching
!             yy = dble(j-1) * dy
!            else
!             ! with stretching
!             yy = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)
!           end if
!           do i = 1, imax
!             ! coordinate of x-direction
!             xx = x0 + dble(i-1) * dx
!             write(3,8) xx, yy, zz, uxp(i,j,kphys)+uxbt(i,j), 
!    &                               uyp(i,j,kphys)+uybt(i,j),
!    &                   uzp(i,j,kphys), wxp(i,j,kphys), wyp(i,j,kphys),
!    &                               wzp(i,j,kphys)+wzbt(i,j),
!    &                            duxpdy(i,j,kphys)+duxbtdy(i,j) 
!           end do
!         end do

!         do k = 1, kphys
!           ! coordinate of z-direction
!           zz = dble(k-1) * dz
!           do j = 1, jmax
!             ! coordinate of y-direction
!             if (stf .eq. 1.d0) then
!               ! without stretching
!               yy = dble(j-1) * dy
!              else
!               ! with stretching
!               yy = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)
!             end if
!             do i = 1, imax
!               ! coordinate of x-direction
!               xx = x0 + dble(i-1) * dx
!             write(3,8) xx, yy, zz, uxp(i,j,k)+uxbt(i,j), 
!    &                               uyp(i,j,k)+uybt(i,j),
!    &                   uzp(i,j,k), wxp(i,j,k), wyp(i,j,k),
!    &                               wzp(i,j,k)+wzbt(i,j),
!    &                            duxpdy(i,j,kphys)+duxbtdy(i,j) 
!             end do
!           end do
!         end do
!         close (unit = 3)

          ! writes data to spacial space to be open by tecplot
          open (3, file = 'spatial.dat',status = 'unknown')
          write(3,*) 'VARIABLES="y","z","u","v","w","wx","wy","wz",
     & "duxdy", "duxdz"'
          write(3,*) 'ZONE J=',jmax,', K=',kphys+1,
     &               ',F=POINT'
          
          ! coordinate of z-direction
          zz = - 1.d0 * dz
          do j = 1, jmax
            ! coordinate of y-direction
            if (stf .eq. 1.d0) then
              ! without stretching
              yy = dble(j-1) * dy
             else
              ! with stretching
              yy = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)
            end if
            do i = 1400, 1400
              ! coordinate of x-direction
              xx = x0 + dble(i-1) * dx
              write(3,9) yy, zz, uxp(i,j,kphys)+uxbt(i,j), 
     &                           uyp(i,j,kphys)+uybt(i,j),
     &                   uzp(i,j,kphys), wxp(i,j,kphys), wyp(i,j,kphys),
     &                           wzp(i,j,kphys)+wzbt(i,j),
     &                           duxpdy(i,j,kphys), 
     &                           duxpdz(i,j,kphys) 
            end do
          end do

          do k = 1, kphys
            ! coordinate of z-direction
            zz = dble(k-1) * dz
            do j = 1, jmax
              ! coordinate of y-direction
              if (stf .eq. 1.d0) then
                ! without stretching
                yy = dble(j-1) * dy
               else
                ! with stretching
                yy = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)
              end if
              do i = 1400, 1400
                ! coordinate of x-direction
                xx = x0 + dble(i-1) * dx
              write(3,9) yy, zz, uxp(i,j,k)+uxbt(i,j), 
     &                           uyp(i,j,k)+uybt(i,j),
     &                   uzp(i,j,k), wxp(i,j,k), wyp(i,j,k),
     &                           wzp(i,j,k)+wzbt(i,j),
     &                           duxpdy(i,j,k), 
     &                           duxpdz(i,j,k) 
              end do
            end do
          end do
          close (unit = 3)


      end select

    5 format(1x, 3d14.6, 3d17.9)
    6 format(1x, 3d14.6, 7d17.9)
    7 format(1x, 3d14.6, 6d17.9)
    8 format(1x, 3d14.6, 7d17.9)
    9 format(1x, 2d14.6, 8d17.9)

      stop
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine create_ctes

      implicit none
      include 'par.for'
      include 'comm.fourier'
      include 'comm.multi'
 
      integer lvl, j, k
      real*8 aux
 
      ! multigrid spatial calculations
 
      ! dy0 at each multigrid level
      do lvl = 1 , msh
        if (stf .ne. 1.d0) then
          v_dy0(lvl) = dy * ( ( stf**(2**(lvl-1)) - 1.d0)
     &               / (stf - 1.d0) )
         else
          v_dy0(lvl) = dy * 2.d0**(lvl-1)
        end if
      end do

      do lvl = 1 , msh
        v_stf(lvl) = stf ** ( 2**(lvl-1) )
        v_dx2(lvl) = (dx * dble(2**(lvl-1)))**2
      end do

      ! dy at each space
      if (stf .ne. 1.d0) then
        do j = 1 , jmax - 1
          v_qdy(j)  = 1.d0 / (v_dy0(1) * v_stf(1)**(j-1))
        end do
       else
        do j = 1 , jmax - 1
          v_qdy(j)  = 1.d0 / (v_dy0(1))
        end do
      end if

      do lvl = 1 , msh
        if (stf .ne. 1.d0) then
          do j = 1 , (jmax - 1) / 2**(lvl-1)
            aux           = (v_dy0(lvl) * v_stf(lvl)**(j-1))
            v_qdy2(j,lvl) = 1.d0 / (aux**2)
          end do
         else
          do j = 1 , (jmax - 1) / 2**(lvl-1)
            aux           = (v_dy0(lvl))
            v_qdy2(j,lvl) = 1.d0 / (aux**2)
          end do
        end if
      end do

      v_ptsx(1) = ptsx
      v_ptsy(1) = jmax
      do lvl = 2 , msh
        v_ptsx(lvl) = (v_ptsx(lvl-1) + 1) / 2
        v_ptsy(lvl) = (v_ptsy(lvl-1) + 1) / 2
      end do

      do k = 1 , kfour
        v_kb(k)   = - im * dble(k-1) * beta
        v_k2b2(k) = - dble(k-1) * dble(k-1) * beta * beta
      end do

      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c              end subroutines transforms binary data in asc                  c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
