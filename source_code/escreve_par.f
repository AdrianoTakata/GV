cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                         subroutines writes the results                       c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine escreve(t)

      ! write the results in binary form
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      character*15 nome
      integer t

      if (my_form.eq.2) then
          write(nome,'(a,i0.2,a)')'data_',my_rank,'.bin'
          write(*,*) ' The results are stored in the file data_xx.bin' 
          open(1,file = nome, form = 'unformatted')
          write(1) t
          write(1) ux, uy, uz, wx, wy, wz, th
          close (unit = 1)
        else
          write(nome,'(a,i0.2,a)')'data_',my_rank,'.bin'
          write(*,*) ' The results are stored in the file data_xx.bin' 
          open(1,file = nome, form = 'unformatted')
          write(1) t
          write(1) ux, uy, uz, wx, wy, wz
          close (unit = 1)
      end if

      return
      end 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine escreve2(t, fanal)

      ! write the results for Fourier analysis in time
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      character*20 nm
      integer i, j, k, t, fanal, var
      real*8 uxb(ptsx,jmax), uyb(ptsx,jmax), wzb(ptsx,jmax)
      complex*16 uxp(ptsx,jmax,kfour), ux_fta(ptsx/8+1,jmax,kfour)
      common/blas/ uxb, uyb, wzb

      if ( (dble(t-tt) / dble(fanal)) .gt. ((t - tt) / fanal)) return
      var = (t - tt) / fanal + 1

      write(nm,'(a,i0.2,a,i0.2,a)')'pert_',my_rank,'_',var,'.bin'
      write(*,*) ' The results are stored in the files pert_cc_XX' 
      select case (my_form)

        case(0, 4, 5)
          do j = 1, jmax
            do i = 1, ptsx
              uxp(i,j,1) = ux(i,j,1) - uxb(i,j)
            end do
          end do
          do k = 2, kfour
            do j = 1, jmax
              do i = 1, ptsx
                uxp(i,j,k) = ux(i,j,k)
              end do
            end do
          end do
          open(1, file = nm, form = 'unformatted')
          write(1) uxp
          close (unit = 1)

        case default
          do k = 1, kfour
            do j = 1, jmax
              do i = 1, ptsx/8+1
                ux_fta(i,j,k) = ux((i-1)*8+1,j,k)
              end do
            end do
          end do
          open(1, file = nm, form = 'unformatted')
          write(1) ux_fta
          close (unit = 1)

      end select

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine escreve_secondary(t)

      ! write the results for Fourier analysis in time
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      character*20 nm
      integer i, j, k, t
      real*8 uxb(ptsx,jmax), uyb(ptsx,jmax), wzb(ptsx,jmax)
      complex*16 uxp(ptsx,jmax,kfour), ux_fta(ptsx/8+1,jmax,kfour)
      common/blas/ uxb, uyb, wzb

      write(nm,'(a,i0.2,a,i0.5,a)')'pert_',my_rank,'_',t,'.bin'
      write(*,*) ' The results are stored in the files pert_cc_XX'
      select case (my_form)

        case(0, 4, 5)
          do j = 1, jmax
            do i = 1, ptsx
              uxp(i,j,1) = ux(i,j,1) - uxb(i,j)
            end do
          end do
          do k = 2, kfour
            do j = 1, jmax
              do i = 1, ptsx
                uxp(i,j,k) = ux(i,j,k)
              end do
            end do
          end do
          do k = 1, kfour
            do j = 1, jmax
              do i = 1, ptsx/8+1
                ux_fta(i,j,k) = uxp((i-1)*8+1,j,k)
              end do
            end do
          end do
          open(1, file = nm, form = 'unformatted')
          write(1) ux_fta
          close (unit = 1)

        case default
          do k = 1, kfour
            do j = 1, jmax
              do i = 1, ptsx/8+1
                ux_fta(i,j,k) = ux((i-1)*8+1,j,k)
              end do
            end do
          end do
          open(1, file = nm, form = 'unformatted')
          write(1) ux_fta
          close (unit = 1)

      end select

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine escreve_phys(var)

      ! write the results in binary form
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      character*15 nome
      integer i, j, k
      real*8 var(ptsx,jmax,kphys), x, y, z

      write(nome,'(a,i0.2,a)')'var_',my_rank,'.dat'
      open (3, file = nome)
      write(3,*) 'VARIABLES="x","y","z","var"'
      write(3,*) 'ZONE I=',ptsx,', J=',jmax,', K=',kphys+1,', F=POINT'

      z = - dz
      do j = 1, jmax
        if (stf .eq. 1.d0) then
          ! without stretching
          y = dble(j-1) * dy
         else
          ! with stretching
          y = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)
        end if
        do i = 1, ptsx
          x = x0 + dble(i-1+shift) * dx
          write(3,6) x, y, z, var(i,j,kphys)
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
          do i = 1, ptsx
            x = x0 + dble(i-1+shift) * dx
            write(3,6) x, y, z, var(i,j,k)
          end do
        end do
      end do
      close (unit = 3)

    6 format(1x, 3d14.6, 1d17.9)

      return
      end

c writing the file for unsteady perturbation

      subroutine escreve_unsteady(cont,t)

      ! write the results in binary form
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      character*15 nome
      integer t, cont

      if (my_form .eq. 2) then
          write(nome,'(a,i0.2,a,i0.2,a)')'data_',my_rank,'_',cont,'.bin'
          write(*,*) ' The results are stored in the file data_xx.bin' 
          open(1,file = nome, form = 'unformatted')
          write(1) t
          write(1) ux, uy, uz, wx, wy, wz, th
          close (unit = 1)
        else
          write(nome,'(a,i0.2,a,i0.2,a)')'data_',my_rank,'_',cont,'.bin'
          write(*,*) ' The results are stored in the file data_xx.bin' 
          open(1,file = nome, form = 'unformatted')
          write(1) t
          write(1) ux, uy, uz, wx, wy, wz
          close (unit = 1)
      end if

      return
      end 


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                     end subroutines writes the results                       c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
