cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                       begin of pentadiagonal solvers                         c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine ludecomp(a, lu)

      ! decomposes the pentadiagonal a in lu where the l and u 
      ! are in the same pentadiagonal matrix
      ! used in the x direction (domain decomposition)
      implicit none
      include 'par.for'
      integer i, j, k, it, jj, ii, kk
      real*8 a(imax,5), lu(imax,5), soma

      do i = 1, imax - 2
        jj = i - 1
        if (jj.gt.2) jj = 2
        do j = 3, 5
          kk   = 2
          soma = 0.d0
          do k = 1, jj
            soma = soma + lu(i,kk) * lu(i-k,j+k)
            kk   = kk - 1
          end do
          jj = jj - 1
          lu(i,j) = a(i,j) - soma
        end do
        j  = 2
        ii = i - 1
        if (ii.gt.1) ii = 1
        do it = i + 1, i + 2
          soma = 0.d0
          do k = 1, ii
            soma = lu(it,k) * lu(i-1,j+2)
          end do
          lu(it,j) = ( a(it,j) - soma ) / lu(i,3)
          ii = ii - 1
          j  = j  - 1
        end do
      end do
      i = imax - 1
      lu(i,3)   = a(i,3) - lu(i,2) * lu(i-1,4) - lu(i,1) * lu(i-2,5)
      lu(i,4)   = a(i,4) - lu(i,2) * lu(i-1,5)
      lu(i+1,2) = ( a(i+1,2) - lu(i+1,1) * lu(i-1,4) ) / lu(i,3)
      i = imax
      lu(i,3)   = a(i,3) - lu(i,2) * lu(i-1,4) - lu(i,1) * lu(i-2,5)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine lusolver_1(lu, rhs)

      ! solves the pentadiagonal problem using the lu and 
      ! rhs is a simple vector woth ptsx elements
      ! used in the x direction (domain decomposition)
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'mpif.h'
      integer status(MPI_Status_size)
      integer i, i_ini, i_fim
      real*8 lu(imax,5)
      complex*16 aux(3), rhs(ptsx), y(ptsx)

      ! begin of solver Ly = rhs
      if (my_rank.gt.0) then
        call MPI_Recv(aux, 3, MPI_COMPLEX16, my_rank - 1,
     &                159, MPI_COMM_WORLD, status, ierr)
        do i = 1, 3
          y(i) = aux(i)
        end do
        i_ini = 4
       else
        y(1)   = rhs(1)
        y(2)   = rhs(2) - lu(2,2) * y(1)
        i_ini  = 3
      end if

      do i = i_ini, ptsx
        y(i)   = rhs(i) - lu(i+shift,2) * y(i-1)
     &                  - lu(i+shift,1) * y(i-2)
      end do

      if (my_rank.lt.numproc) then
        do i = 1, 3
          aux(i) = y(i+ptsx-inter-1)
        end do
        call MPI_Send(aux, 3, MPI_COMPLEX16, my_rank + 1, 
     &                159, MPI_COMM_WORLD, ierr)
      end if
      ! end of solver Ly = rhs

      ! begin of solver U(rhs) = y
      if (my_rank.lt.numproc) then
        call MPI_Recv(aux, 3, MPI_COMPLEX16, my_rank + 1, 
     &                169, MPI_COMM_WORLD, status, ierr)
        do i = 1, 3
          rhs(i + ptsx - 3) = aux(i)
        end do
        i_fim = ptsx - 3
       else
        i      = ptsx
        rhs(i) = y(i) / lu(i+shift,3)
        i      = ptsx - 1
        rhs(i) = ( y(i) - lu(i+shift,4) * rhs(ptsx) ) /
     &           lu(i+shift,3)
        i_fim  = ptsx - 2
      end if
      
      do i = i_fim, 1, -1
        rhs(i) = ( y(i) - lu(i+shift,4) * rhs(i+1)
     &                  - lu(i+shift,5) * rhs(i+2) ) / lu(i+shift,3)
      end do
      
      if (my_rank.gt.0) then
        do i = 1, 3
          aux(i) = rhs(i + inter - 2)
        end do
        call MPI_Send(aux, 3, MPI_COMPLEX16, my_rank - 1, 
     &                169, MPI_COMM_WORLD, ierr)
      end if
      ! end of solver U(rhs) = y

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine lusolver(lu, rhs)

      ! solves the pentadiagonal problem using the lu and 
      ! rhs is a matrix with ptsx * jmax elements
      ! used in the x direction (domain decomposition)
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'mpif.h'
      integer status(MPI_Status_size)
      integer i, j, i_ini, i_fim
      real*8 lu(imax,5)
      complex*16 aux(3), rhs(ptsx,jmax), y(ptsx,jmax)

      ! begin of solver Ly = rhs
      do j = 2, jmax
        if (my_rank .gt. 0) then
          call MPI_Recv(aux, 3, MPI_COMPLEX16, my_rank - 1, 
     &                  259, MPI_COMM_WORLD, status, ierr)
          do i = 1, 3
            y(i,j) = aux(i)
          end do
          i_ini = 4
         else
          y(1,j) = rhs(1,j)
          y(2,j) = rhs(2,j) - lu(2,2) * y(1,j)
          i_ini  = 3
        end if
        
        do i = i_ini, ptsx
          y(i,j) = rhs(i,j) - lu(i+shift,1) * y(i-2,j)
     &                      - lu(i+shift,2) * y(i-1,j)
        end do
        
        if (my_rank .lt. numproc) then
          do i = 1, 3 
            aux(i) = y(i+ptsx-inter-1,j)
          end do
          call MPI_Send(aux, 3, MPI_COMPLEX16, my_rank + 1, 
     &                  259, MPI_COMM_WORLD, ierr)
        end if
      end do
      ! end of solver Ly = rhs

      ! begin of solver U(rhs) = y
      do j = 2, jmax
        if (my_rank.lt.numproc) then
          call MPI_Recv(aux, 3, MPI_COMPLEX16, my_rank + 1, 
     &                  269, MPI_COMM_WORLD, status, ierr)
          do i = 1, 3
            rhs(ptsx + i - 3,j) = aux(i)
          end do
          i_fim = ptsx - 3
         else
          i        = ptsx
          rhs(i,j) = y(i,j) / lu(i+shift,3)
          i        = ptsx-1
          rhs(i,j) = (y(i,j) - lu(i+shift,4) * rhs(i+1,j)) /
     &                         lu(i+shift,3)
          i_fim    = ptsx - 2
        end if
       
        do i = i_fim, 1, -1
          rhs(i,j) = ( y(i,j) - lu(i+shift,4) * rhs(i+1,j)
     &               - lu(i+shift,5) * rhs(i+2,j) ) / lu(i+shift,3)
        end do
       
        if (my_rank .gt. 0) then
          do i = 1, 3
            aux(i) = rhs(inter + i - 2,j)
          end do
          call MPI_Send(aux, 3, MPI_COMPLEX16, my_rank - 1, 
     &                  269, MPI_COMM_WORLD, ierr)
        end if
      end do
      ! end of solver U(rhs) = y

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine bandy5(a, al, indx)

      ! solve the LHS of one pentadiagonal matrix in y direction
      ! the a is the input and a, al and indx are outputs
      implicit none
      include 'par.for'
      integer m1, indx(jmax), i, j, k, l, mm
      real*8 d, a(jmax,5), al(jmax,5), tiny, dum
      parameter (tiny = 1.d-20)

      m1 = 2
      mm = 5
      l  = m1
      do i = 1, m1
        do j = m1 + 2 - i, mm
          a(i,j-l) = a(i,j)
        end do
        l = l - 1
        do j = mm - l, mm
          a(i,j) = 0.d0
        end do
      end do
      d = 1.d0
      l = m1
      do k = 1, jmax
        dum = a(k,1)
        i   = k
        if (l.lt.jmax) l = l + 1
!       do j = k + 1, l
!         if (dabs(a(j,1)) .gt. dabs(dum)) then
!           dum = a(j,1)
!           i   = j
!         end if
!       end do
        indx(k) = i
        if (dum .eq. 0.d0) a(k,1) = tiny
        if (i .ne. k) then
          d = - d
          do j = 1, mm
            dum    = a(k,j)
            a(k,j) = a(i,j)
            a(i,j) = dum
          end do
        end if
        do i = k + 1, l
          dum       = a(i,1) / a(k,1)
          al(k,i-k) = dum
          do j = 2, mm
            a(i,j-1) = a(i,j) - dum * a(k,j)
          end do
          a(i,mm) = 0.d0
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine banbky5(a, al, indx, rhs)

      ! solve the the pentadiagonal matrix in y direction
      ! the term a and al comes from the subroutine bandy5
      ! the rhs variable is the input and at the end
      ! is the result of the solved problem
      implicit none
      include 'par.for'
      integer indx(jmax), i, k, l
      real*8 a(jmax,5), al(jmax,5)
      complex*16 rhs(jmax), dum

      l = 2
      do k = 1, jmax
        i = indx(k)
        if (i .ne. k) then
          dum    = rhs(k)
          rhs(k) = rhs(i)
          rhs(i) = dum
        end if
        if (l .lt. jmax) l = l + 1
        do i = k + 1, l
          rhs(i) = rhs(i) - al(k,i-k) * rhs(k)
        end do
      end do
      l = 1
      do i = jmax, 1, -1
        dum = rhs(i)
        do k = 2, l
          dum = dum - a(i,k) * rhs(k+i-1)
        end do
        rhs(i) = dum / a(i,1)
        if (l.lt.5) l = l + 1
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine band5_poi(a, n, al, indy)

      ! solve the LHS of the pentadiagonal matrix
      ! n is the size of the matrix
      ! used for uy Poisson solver subroutine
      implicit none
      include 'par.for'
      integer m1, n, indy(jmax), i, j, k, l, mm
      real*8 d, a(jmax,5), al(jmax,5), TINY, dum
      parameter (TINY = 1.e-20)

      m1 = 2
      mm = 5
      l  = m1
      do i = 1, m1
        do j = m1 + 2 - i, mm
          a(i,j-l) = a(i,j)
        end do
        l = l - 1
        do j = mm - l, mm
          a(i,j) = 0.d0
        end do
      end do
      d = 1.d0
      l = m1
      do k = 1, n
        dum = a(k,1)
        i   = k
        if (l.lt.n) l = l + 1
!       do j = k + 1, l
!         if (dabs(a(j,1)) .gt. dabs(dum)) then
!           dum = a(j,1)
!           i   = j
!         end if
!       end do
        indy(k) = i
        if (dum .eq. 0.d0) a(k,1) = TINY
        if (i .ne. k) then
          d = - d
          do j = 1, mm
            dum    = a(k,j)
            a(k,j) = a(i,j)
            a(i,j) = dum
          end do
        end if
        do i = k + 1, l
          dum       = a(i,1) / a(k,1)
          al(k,i-k) = dum
          do j = 2, mm
            a(i,j-1) = a(i,j) - dum * a(k,j)
          end do
          a(i,mm) = 0.d0
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine banbk5_poi(a, n, al, indy, rhs)

      ! solve the the pentadiagonal matrix the term a and al comes from
      ! subroutine band5_poi, the rhs variable is the input and at the
      ! end is the result of matrix
      ! n is the size of the matrix
      ! used for uy Poisson solver subroutine
      implicit none
      include 'par.for'
      integer n, i, k, l, indy(jmax)
      real*8 a(jmax,5), al(jmax,5)
      complex*16 rhs(jmax), dum

      l = 2
      do k = 1, n
        i = indy(k)
        if (i .ne. k) then
          dum    = rhs(k)
          rhs(k) = rhs(i)
          rhs(i) = dum
        end if
        if (l.lt.n) l = l + 1
        do i = k + 1, l
          rhs(i) = rhs(i) - al(k,i-k) * rhs(k)
        end do
      end do
      l = 1
      do i = n, 1, -1
        dum = rhs(i)
        do k = 2, l
          dum = dum - a(i,k) * rhs(k+i-1)
        end do
        rhs(i) = dum / a(i,1)
        if (l.lt.5) l = l + 1
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                      end of pentadiagonal solvers                            c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
