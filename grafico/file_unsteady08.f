!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      program file_vtk

      implicit none
      include '../marensi/dx202/unsteady/blasius/tt51840/par.for'
      integer i, j, k, pos_x, cont, perd
      character(200) path1, path2, nome
      real*8 Q(imax,jmax,kphys), 
     &       qv(((imax+1)/2)*jmax*(kphys+1)),
     &       ux(imax,jmax,kphys), th(imax,jmax,kphys),
     &       uyp(imax,jmax,kphys), uzp(imax,jmax,kphys),
     &       uxv(((imax+1)/2)*jmax*(kphys+1)),
     &       thv(((imax+1)/2)*jmax*(kphys+1)),
     &       x(((imax+1)/2)*jmax*(kphys+1)),
     &       y(((imax+1)/2)*jmax*(kphys+1)), 
     &       z(((imax+1)/2)*jmax*(kphys+1))

      do perd = 0, 12, 4
        write(*,*)'perd = ', perd
        write(path1, '(a,i0.2,a)')
     &'../marensi/dx202/unsteady/blasius/tt51840/isoq_',perd,'.bin'
        open(unit=8, FILE=path1, form='unformatted')
        read(8) Q
        close(unit=8)

        write(path2, '(a,i0.2,a)')
     &'../marensi/dx202/unsteady/blasius/tt51840/spatial_',perd,'.bin'
        open(unit=8, FILE=path2, form='unformatted')
        if (my_form .eq. 1) then
          read(8) ux, uyp, uzp
        else
          read(8) ux, th, uyp, uzp
        end if
        close(unit=8)   

        cont = 1
        do k = 1, kphys+1
          do j = 1, jmax
            do i = 1, imax, 2
              x(cont) = (x0 + dble(i-1)*dx)/5.0
              y(cont) = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)
              z(cont) = dble(k-1)*dz
              if (k .ne. (kphys+1)) then
                qv(cont) = Q(i,j,k)
                uxv(cont) = ux(i,j,k)
                thv(cont) = th(i,j,k)
              else
                qv(cont) = Q(i,j,1)
                uxv(cont) = ux(i,j,1)
                thv(cont) = th(i,j,1)
              end if
              cont = cont + 1
            end do
          end do
        end do

        write(nome, '(a,i0.2,a)')'lam08_fgv05_perd',perd,'.vtk'
        open(unit=1, file=nome, form='formatted')
        write(1,9),'# vtk DataFile Version 3.0'
        write(1,9),'vtk output'
        write(1,9),'ASCII'
        write(1,9),'DATASET STRUCTURED_GRID'
        write(1,7),'DIMENSIONS',(imax+1)/2, jmax, kphys+1
        write(1,8),'POINTS',((imax+1)/2)*jmax*(kphys+1),'double'
        write(1,9),''

        do i = 1, ((imax+1)/2)*jmax*(kphys+1)
          write(1,5) x(i), y(i), z(i)
        end do
        write(1,10)'CELL_DATA',((imax-1)/2)*(jmax-1)*(kphys)
        write(1,10)'POINT_DATA', (imax/2+1)*jmax*(kphys+1)
        write(1,*)''
        write(1,9)'SCALARS isoq double'
        write(1,9)'LOOKUP_TABLE default'
        do i = 1, ((imax+1)/2)*jmax*(kphys+1)
          write(1,6) qv(i)
        end do

        write(1,10)'CELL_DATA',((imax-1)/2)*(jmax-1)*(kphys)
        write(1,10)'POINT_DATA', (imax/2+1)*jmax*(kphys+1)
        write(1,*)''
        write(1,9)'SCALARS U double'
        write(1,9)'LOOKUP_TABLE default'
        do i = 1, ((imax+1)/2)*jmax*(kphys+1)
          write(1,6) uxv(i)
        end do

        write(1,10)'CELL_DATA',((imax-1)/2)*(jmax-1)*(kphys)
        write(1,10)'POINT_DATA', (imax/2+1)*jmax*(kphys+1)
        write(1,*)''
        write(1,9)'SCALARS U double'
        write(1,9)'LOOKUP_TABLE default'
        do i = 1, ((imax+1)/2)*jmax*(kphys+1)
          write(1,6) thv(i)
        end do

        close(unit=1)
      end do
 
    5  format(F11.8, F11.8, F11.8)
    6  format(F15.8)
    7  format(A, I5, I4, I3)
    8  format(A, I10, A7)
    9  format(A)
   10  format(A, I10)
 
      end
