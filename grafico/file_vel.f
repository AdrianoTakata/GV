      program file_dat

      implicit none
      include '../SB/unsteady/lam18/blasius/par.for'
      integer i, j, k, pos_x
      character(200) path, nome
      real*8 ux(imax,jmax), uy(imax,jmax),
     &       wz(imax,jmax),th(imax,jmax),
     &       Y(imax,jmax), X(imax,jmax)

      write(path, '(a,i0.2,a)'),
     & '../SB/unsteady/lam18/blasius/basens.bin'
      open(unit=8, FILE=path, form='unformatted')
      read(8) ux, uy, wz, th 
      close(unit=8)

      do j = 1, jmax
        do i = 1, imax
          Y(i,j) = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)
          X(i,j) = x0+dble(i-1)*dx
        end do
      end do

      write(nome, '(a)'),'lam18_blasius.dat'
      open(unit=1, file=nome, form='formatted')
      do j = 1, jmax
        do i = 1, imax
          write(1,6) X(i,j), Y(i,j), ux(i,j), th(i,j)
        end do
      end do
      close(unit=1)

    6  format(F11.8, F11.8, F11.8, F11.8)
      
      end
