      program file_dat

      implicit none
      include '../SB/unsteady/lam36/nssc/fgv15/tt432000/par.for'
      integer i, j, k, pos_x, cont
      character(200) path, path1
      real*8 ux(imax,jmax,kphys), th(imax,jmax,kphys),
     &       uxp(imax,jmax,kphys), uyp(imax,jmax,kphys),
     &       uzp(imax,jmax,kphys),
     &       Y(jmax,kphys), Z(jmax,kphys),
     &       Uyz(jmax,kphys), Tyz(jmax,kphys),
     &       Vp(jmax, kphys), Wp(jmax,kphys)

      do i = 0, 12, 4
        write(path, '(a,i0.2,a)'),
     &   '../SB/unsteady/lam36/nssc/fgv15/tt432000/spatial_',i,'.bin'
        open(unit=8, FILE=path, form='unformatted')
        if (my_form .eq. 1 ) then
          read(8) ux, uyp, uzp 
        else
          read(8) ux, th, uyp, uzp 
        end if
        close(unit=8)

        do k = 1, kphys
          do j = 1, jmax
            Y(j,k) = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)
            Z(j,k) = dble(k-1)*dz
          end do
        end do

        pos_x = 1700 !imax=1465

        write(*,*)'pos_x=', pos_x
        Uyz = ux(pos_x, :, :)
        Tyz = th(pos_x, :, :)
        Vp = uyp(pos_x,:,:)
        Wp = uzp(pos_x,:,:)
        call save_dat(i, Y, Z, Uyz, Tyz, Vp, Wp)

      end do

      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine save_dat(cont, Y, Z, Uyz, Tyz, Vp, Wp)

      implicit none
      include '../SB/unsteady/lam36/nssc/fgv15/tt432000/par.for'
      character(200) nome
      integer i, j, k, cont, dim_z
      real*8 Uyz(jmax,kphys), Y(jmax,kphys), Z(jmax,kphys), 
     &       Tyz(jmax,kphys), Vp(jmax,kphys), Wp(jmax,kphys)

      write(nome, '(a,i0.2,a)'),'lam36_fgv15_perd',cont,'.dat'
      open(unit=1, file=nome, form='formatted')
      do k = 1, kphys
        do j = 1, jmax
          if (my_form .eq. 1) then
            write(1,6) Y(j,k),Z(j,k),Uyz(j,k),
     &                 Vp(j,k),Wp(j,k)
          else
            write(1,6) Y(j,k),Z(j,k),Uyz(j,k),
     &                Tyz(j,k),Vp(j,k),Wp(j,k)
          end if
        end do
      end do
      close(unit=1)

    6  format(F11.8,F12.8,F12.8,F12.8,F13.8,F13.8)
      
      end
