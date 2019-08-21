cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                               loop calculations                              c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine bzone

      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      include 'mpif.h'
      integer i, j, k, bdt
      real*8 ep, fc, bdfc(2,imax),
     &       uxb(ptsx,jmax), uyb(ptsx,jmax), wzb(ptsx,jmax)
      common/blas/ uxb, uyb, wzb
      common/bd/ bdfc

      do k = 1, kfour
        bdt = 1
        if (k .eq. 3) bdt = 2
        do j = 2, jmax
          do i = 1, ptsx
            wx(i,j,k) = wx(i,j,k) * bdfc(bdt,i+shift)
            wy(i,j,k) = wy(i,j,k) * bdfc(bdt,i+shift)
            select case (my_form)
              case(1)!-> Gortler vortices simulation without heat transfer
                wz(i,j,k) = wz(i,j,k) * bdfc(bdt,i+shift)
              case(2)!-> Gortler vortices simulation with heat transfer
                wz(i,j,k) = wz(i,j,k) * bdfc(bdt,i+shift)
                th(i,j,k) = th(i,j,k) * bdfc(bdt,i+shift)
            end select
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine outuy(dwzdx)

      ! outflow uy calculation with 6th order compact appx.
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      include 'comm.coef'
      integer j, k, indx(jmax)
      real*8 a(jmax,5), al(jmax,5), uxb(ptsx,jmax), uyb(ptsx,jmax),
     &       wzb(ptsx,jmax), dyya, dya
      complex*16 dwzdx(ptsx,jmax,kfour), rhs(jmax)
      common/blas/ uxb, uyb, wzb

      call derparx(dwzdx, wz)

      ! uy calculation for the last point in the last computer
      ! these values are used as dirichlet BC at end of domain
      if (my_rank .eq. numproc) then
        do k = 1, kfour

c         call lhsoutuy(a, k)
c         call bandy5(a, al, indx)

c         rhs(1) =  dcmplx(0.d0,0.d0)
c         dyya   =  dyy
c         rhs(2) = - dyya * ( sp_poi_coef(2,1) * dwzdx(ptsx,2,k)
c    &                      + sp_poi_coef(3,1) * dwzdx(ptsx,3,k) )

c         do j = 3, jmax - 2
c           dyya = ( dy * stf**(j-3) )**2
c           rhs(j) = - dyya * ( cp_poi_coef(1,1) * dwzdx(ptsx,j-1,k)
c    &                        + cp_poi_coef(2,1) * dwzdx(ptsx,j,k)
c    &                        + cp_poi_coef(3,1) * dwzdx(ptsx,j+1,k) )
c         end do

c         dyya = ( dy * stf**(jmax-3) )**2
c         rhs(jmax-1) = - dyya * dwzdx(ptsx,jmax-1,k)

c         rhs(jmax) = - dyya * dwzdx(ptsx,jmax,k)
c         dya = dsqrt(dyya)
c         if (k .eq. 1) then
c          rhs(jmax) = - dyya * lp_poi_coef(1,1) * dwzdx(ptsx,jmax,k)
c    &                 + dya  * lp_poi_coef(2,1) * duexmdx(ptsx)
c    &                 - dya  * lp_poi_coef(2,1) * alpha *
c    &                   uyb(ptsx,jmax)
c         end if

c         call banbky5(a, al, indx, rhs)

          do j = 2, jmax
c           uy(ptsx,j,k) = rhs(j)
            uy(ptsx,j,k) = dcmplx(0.d0,0.d0)
          end do
    
        end do
      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine lhsoutuy(a, k)

      ! LHS of the uy calculation
      implicit none
      include 'par.for'
      include 'comm.coef'
      include 'comm.fourier'
      integer j, k
      real*8 a(jmax,5), dyya, dya

      a(1,1) = 0.d0
      a(1,2) = 0.d0
      a(1,3) = 1.d0
      a(1,4) = 0.d0
      a(1,5) = 0.d0

      dyya = dyy / fac_y
      a(2,1) = 0.d0
      a(2,2) = sp_poi_coef(4,1)
      a(2,3) = sp_poi_coef(5,1) + sp_poi_coef(2,1) * v_k2b2(k) * dyya
      a(2,4) = sp_poi_coef(6,1) + sp_poi_coef(3,1) * v_k2b2(k) * dyya
      a(2,5) = sp_poi_coef(7,1)

      do j = 3, jmax - 2
        dyya = ( dy * stf**(j-3) ) ** 2 / fac_y
        a(j,1) = cp_poi_coef(4,1)
        a(j,2) = cp_poi_coef(5,1) + cp_poi_coef(1,1) * v_k2b2(k) * dyya
        a(j,3) = cp_poi_coef(6,1) + cp_poi_coef(2,1) * v_k2b2(k) * dyya
        a(j,4) = cp_poi_coef(7,1) + cp_poi_coef(3,1) * v_k2b2(k) * dyya
        a(j,5) = cp_poi_coef(8,1)
      end do

      dyya = ( dy * stf**(jmax-3) ) ** 2
      a(jmax-1,1) = 0.d0
      a(jmax-1,2) = pp_poi_coef(2,1)
      a(jmax-1,3) = pp_poi_coef(3,1) + v_k2b2(k) * dyya / fac_y
      a(jmax-1,4) = pp_poi_coef(4,1)
      a(jmax-1,5) = 0.d0

      dya = dsqrt(dyya)
      a(jmax,1) = lp_poi_coef(5,1)
      a(jmax,2) = lp_poi_coef(4,1)
      a(jmax,3) = lp_poi_coef(3,1) + v_k2b2(k) * dyya / fac_y
     &          - lp_poi_coef(2,1) * dya * dsqrt( (alpha * alpha
     &          - v_k2b2(k)) / fac_y)
      a(jmax,4) = 0.d0
      a(jmax,5) = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine wx_wall(lapv)

      ! calculate the vorticity in x direction 
      ! at wall using the v-Poisson equation, modified for that.
      ! The variable lapv is calculated here and is also used 
      ! in the subroutine wz_wall
      ! Eq. 2.32 do Stemmer
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      include 'comm.coef'
      include 'comm.fourier'
      integer i, k, i_ini
      real*8 a(imax,5), lu(imax,5), k2b2
      complex*16 d2wydxdy(ptsx), rhs(ptsx),
     &           lapv(ptsx,kfour), d2uydx2(ptsx,kfour)
      common/derw/ d2uydx2

      ! this variable (lapv) is used in wz_wall sub.
      do k = 1, kfour
        do i = 1, ptsx
          lapv(i,k) = ( d2uydx2(i,k) + v_k2b2(k) * uy(i,1,k) ) / fac_y
     &              + ( w_at_w_coef(3) * uy(i,1,k)
     &                + w_at_w_coef(4) * uy(i,2,k)
     &                + w_at_w_coef(5) * uy(i,3,k)
     &                + w_at_w_coef(6) * uy(i,4,k)
     &                + w_at_w_coef(7) * uy(i,5,k)
     &                + w_at_w_coef(8) * uy(i,6,k)
     &                + w_at_w_coef(9) * uy(i,7,k) )
     &              / ( dyy * w_at_w_coef(1) )
        end do
      end do

      i_ini = 1
      if (my_rank .eq. 0) i_ini = 2

      do k = 2, kfour

        k2b2 = dxx * v_k2b2(k) 
        call lhs5(a, k2b2)
        call ludecomp(a, lu)

        call der2wydxdy(d2wydxdy, k)

        rhs(1) = wx(1,1,k)

        rhs(2) = dxx     * ( -         d2wydxdy(1)
     &                       - 10.d0 * d2wydxdy(2)
     &                       -         d2wydxdy(3)
     &         + v_kb(k) * (           lapv(1,k)
     &                       + 10.d0 * lapv(2,k)
     &                       +         lapv(3,k) ) ) / 12.d0

        do i = 3, ptsx - 2
          rhs(i) = 8.d0 * dxx * ( -         d2wydxdy(i-1)
     &                            - 5.5d0 * d2wydxdy(i)
     &                            -         d2wydxdy(i+1)
     &           + v_kb(k)    * (           lapv(i-1,k)
     &                            + 5.5d0 * lapv(i,k)
     &                            +         lapv(i+1,k) ) )
        end do

        rhs(ptsx-1) = dxx     * ( -         d2wydxdy(ptsx)
     &                            - 10.d0 * d2wydxdy(ptsx-1)
     &                            -         d2wydxdy(ptsx-2)
     &              + v_kb(k) * (           lapv(ptsx,k)
     &                            + 10.d0 * lapv(ptsx-1,k)
     &                            +         lapv(ptsx-2,k) ) ) / 12.d0

        rhs(ptsx) = ( d2wydxdy(ptsx) - v_kb(k) * lapv(ptsx,k) )
     &            / ( dble(k - 1) * dble(k - 1) * beta * beta )
        call lusolver_1(lu, rhs)
        
        do i = i_ini, ptsx
          wx(i,1,k) = rhs(i)
        end do

      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine lhs5(a, k2b2)

      ! lhs for wx_wall, wz_wall, poi_ux and poi_uz subroutines
      implicit none
      include 'par.for'
      integer i
      real*8 a(imax,5), k2b2

      a(1,1) = 0.d0
      a(1,2) = 0.d0
      a(1,3) = 1.d0
      a(1,4) = 0.d0
      a(1,5) = 0.d0

      a(2,1) =   0.d0
      a(2,2) =   1.d0 + k2b2 / 12.d0
      a(2,3) = - 2.d0 + k2b2 / 1.2d0
      a(2,4) =   1.d0 + k2b2 / 12.d0
      a(2,5) =   0.d0

      do i = 3, imax - 2
        a(i,1) =     3.d0
        a(i,2) =    48.d0 +  8.d0 * k2b2
        a(i,3) = - 102.d0 + 44.d0 * k2b2
        a(i,4) =    48.d0 +  8.d0 * k2b2
        a(i,5) =     3.d0
      end do

      a(imax-1,1) =   0.d0
      a(imax-1,2) =   1.d0 + k2b2 / 12.d0
      a(imax-1,3) = - 2.d0 + k2b2 / 1.2d0
      a(imax-1,4) =   1.d0 + k2b2 / 12.d0
      a(imax-1,5) =   0.d0

      a(imax,1) = 0.d0
      a(imax,2) = 0.d0
      a(imax,3) = 1.d0
      a(imax,4) = 0.d0
      a(imax,5) = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine der2wydxdy(d2wydxdy, k)

      ! calculates the d2wy/dxdy to be used in wx_wall subroutine
      implicit none
      include 'par.for'
      include 'comm.var'
      include 'comm.coef'
      integer i, k
      complex*16 d2wydxdy(ptsx), dwydy(ptsx)

! alterar esta aproximação com 7 pontos de 7a ordem para calculo do dwydy
      do i = 1, ptsx
        dwydy(i) = ( 1.d0 / dy ) * ( dwydy_coef(2) * wy(i,1,k)
     &                             + dwydy_coef(3) * wy(i,2,k)
     &                             + dwydy_coef(4) * wy(i,3,k)
     &                             + dwydy_coef(5) * wy(i,4,k)
     &                             + dwydy_coef(6) * wy(i,5,k)
     &                             + dwydy_coef(7) * wy(i,6,k)
     &                             + dwydy_coef(8) * wy(i,7,k) )
      end do

      d2wydxdy(1) = ( - 14.7d0 * dwydy(1) + 36.d0 * dwydy(2)
     &                - 45.d0  * dwydy(3) + 40.d0 * dwydy(4)
     &                - 22.5d0 * dwydy(5) + 7.2d0 * dwydy(6)
     &                -          dwydy(7) ) / ( 6.d0 * dx )
      
      d2wydxdy(2) = ( -         dwydy(1) - 7.7d0 * dwydy(2)
     &                + 15.d0 * dwydy(3) - 10.d0 * dwydy(4)
     &                +  5.d0 * dwydy(5) - 1.5d0 * dwydy(6)
     &                + 0.2d0 * dwydy(7) ) / ( 6.d0 * dx )

      d2wydxdy(3) = (   0.2d0 * dwydy(1) - 2.4d0 * dwydy(2)
     &                - 3.5d0 * dwydy(3) +  8.d0 * dwydy(4)
     &                -  3.d0 * dwydy(5) + 0.8d0 * dwydy(6)
     &                - 0.1d0 * dwydy(7) ) / ( 6.d0 * dx )

      do i = 4, ptsx - 3
        d2wydxdy(i) = ( 45.d0 * ( dwydy(i+1) - dwydy(i-1) )
     &                +  9.d0 * ( dwydy(i-2) - dwydy(i+2) )
     &                +           dwydy(i+3) - dwydy(i-3) )
     &              / ( 60.d0 * dx )
      end do

      d2wydxdy(ptsx-2) = (        dwydy(ptsx-4) - 8.d0 * dwydy(ptsx-3)
     &                   + 8.d0 * dwydy(ptsx-1) -        dwydy(ptsx) )
     &                   / ( 12.d0 * dx )
      
      d2wydxdy(ptsx-1) = ( 2.d0 * dwydy(ptsx)   + 3.d0 * dwydy(ptsx-1)
     &                   - 6.d0 * dwydy(ptsx-2) + dwydy(ptsx-3) )
     &                   / ( 6.d0 * dx )
      
      d2wydxdy(ptsx)   = ( 3.d0 * dwydy(ptsx) - 4.d0 * dwydy(ptsx-1)
     &                   +        dwydy(ptsx-2) ) / ( 2.d0 * dx )

      call boundary_exchange_1d(d2wydxdy)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine wz_wall(lapv)

      ! calculate the vorticity in z direction 
      ! at wall using the v-Poisson equation
      ! the variable lapv comes from subroutine wx_wall
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      include 'comm.fourier'
      include 'mpif.h'
      integer status(MPI_status_size)
      integer i, k
      complex*16 rhs(ptsx), lapv(ptsx,kfour), aux(2)

      do k = 1, kfour

        do i = 1, ptsx
          rhs(i) = v_kb(k) * wx(i,1,k) - lapv(i,k)
        end do

        if (my_rank .gt. 0) then
          call MPI_Recv(aux, 2, MPI_COMPLEX16, my_rank - 1, 
     &                  571, MPI_COMM_WORLD, status, ierr)
          wz(1,1,k) = aux(1)
          wz(2,1,k) = aux(2)
         else
          wz(2,1,k) = ( dx * ( 251.d0 * rhs(1)   + 646.d0 * rhs(2)
     &                       - 264.d0 * rhs(3)   + 106.d0 * rhs(4)
     &                       -  19.d0 * rhs(5) ) + 720.d0 * wz(1,1,k) )
     &                / 720.d0
        end if

        do i = 3, ptsx - 2
          wz(i,1,k) = ( dx * (  281.d0 * rhs(i-2) + 2056.d0 * rhs(i-1)
     &                       + 1176.d0 * rhs(i)   -  104.d0 * rhs(i+1)
     &                       +   11.d0 * rhs(i+2) ) / 90.d0
     &                + 11.d0 * wz(i-2,1,k) + 16.d0 * wz(i-1,1,k) )
     &                / 27.d0
        end do

        if (my_rank .lt. numproc) then
          aux(1) = wz(ptsx-inter,1,k)
          aux(2) = wz(ptsx-inter+1,1,k)
          call MPI_Send(aux, 2, MPI_COMPLEX16, my_rank + 1, 
     &                  571, MPI_COMM_WORLD, ierr)
         else
          i = ptsx - 1
          wz(i,1,k) = ( dx * ( 10.d0 * rhs(i-2) + 57.d0 * rhs(i-1)
     &                       + 24.d0 * rhs(i)   -         rhs(i+1) )
     &                + 33.d0 * wz(i-2,1,k) + 24.d0 * wz(i-1,1,k) )
     &                / 57.d0
          i = ptsx
          wz(i,1,k) = ( dx * ( 251.d0 * rhs(i)   + 646.d0 * rhs(i-1)
     &                       - 264.d0 * rhs(i-2) + 106.d0 * rhs(i-3)
     &                       -  19.d0 * rhs(i-4) )
     &                + 720.d0 * wz(i-1,1,k) ) / 720.d0
         end if

      end do

      call boundary_exchange_wwall(wz)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine poi_ux(duydy)

      ! calculate the new ux velocity using the continuity equation for
      ! fundamental mode and Poisson equation for other modes the
      ! argument duydy is calculated here and used in poi_uz sub.
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      include 'comm.fourier'
      include 'comm.fs'
      include 'mpif.h'
      integer status(MPI_status_size)
      integer i, j, k, i_ini
      real*8 a(imax,5), lu(imax,5), k2b2
      complex*16 rhs(ptsx,jmax), aux(2), duydy(ptsx,jmax,kfour),
     &           d2uydxdy(ptsx,jmax,kfour)

      call deryfv(duydy,uy)
      i_ini = 1
      if (my_rank.eq.0) i_ini = 2

      ! for the fundamental mode, the continuity equation is used
      k = 1
      do j = 2, jmax - 1
        if (my_rank .gt. 0) then
          call MPI_Recv(aux, 2, MPI_COMPLEX16, my_rank - 1, 
     &                  153, MPI_COMM_WORLD, status, ierr)
          ux(1,j,k) = aux(1)
          ux(2,j,k) = aux(2)
         else
          ux(2,j,k) = ( - dx * ( 251.d0 * duydy(1,j,k)
     &                         + 646.d0 * duydy(2,j,k)
     &                         - 264.d0 * duydy(3,j,k)
     &                         + 106.d0 * duydy(4,j,k)
     &                         -  19.d0 * duydy(5,j,k) )
     &                + 720.d0 * ux(1,j,k) ) / 720.d0
        end if

        do i = 3, ptsx - 2
          ux(i,j,k) = ( - dx * (  281.d0 * duydy(i-2,j,k)
     &                         + 2056.d0 * duydy(i-1,j,k)
     &                         + 1176.d0 * duydy(i,j,k)
     &                         -  104.d0 * duydy(i+1,j,k)
     &                         +   11.d0 * duydy(i+2,j,k) ) / 90.d0
     &                + 11.d0 * ux(i-2,j,k) + 16.d0 * ux(i-1,j,k) )
     &                / 27.d0
        end do

        if (my_rank .lt. numproc) then
          aux(1) = ux(ptsx-inter,j,k)
          aux(2) = ux(ptsx-inter+1,j,k)
          call MPI_Send(aux, 2, MPI_COMPLEX16, my_rank + 1, 
     &                  153, MPI_COMM_WORLD, ierr)
         else
          i = ptsx - 1
          ux(i,j,k) = ( - dx * ( 10.d0 * duydy(i-2,j,k) 
     &                         + 57.d0 * duydy(i-1,j,k)
     &                         + 24.d0 * duydy(i,j,k)
     &                         -         duydy(i+1,j,k) )
     &                + 33.d0 * ux(i-2,j,k) + 24.d0 * ux(i-1,j,k) )
     &                / 57.d0
          i = ptsx
          ux(i,j,k) = ( - dx * ( 251.d0 * duydy(i,j,k)
     &                         + 646.d0 * duydy(i-1,j,k)
     &                         - 264.d0 * duydy(i-2,j,k)
     &                         + 106.d0 * duydy(i-3,j,k)
     &                         -  19.d0 * duydy(i-4,j,k) )
     &                + 720.d0 * ux(i-1,j,k) ) / 720.d0
        end if
      end do

      call boundary_exchange_1stmode(ux)

      ! for other modes a ux_Poisson equation is used
      call derparx(d2uydxdy, duydy)

      do k = 2, kfour

        k2b2 = dxx * v_k2b2(k) 
        call lhs5(a, k2b2)
        call ludecomp(a, lu)

        do j = 2, jmax

          rhs(1,j) = ux(1,j,k)

          rhs(2,j) = - dxx   * (         d2uydxdy(1,j,k)
     &                         + 10.d0 * d2uydxdy(2,j,k)
     &                         +         d2uydxdy(3,j,k)
     &             + v_kb(k) * (         wy(1,j,k)
     &                         + 10.d0 * wy(2,j,k)
     &                         +         wy(3,j,k) ) ) / 12.d0

          do i = 3, ptsx - 2
            rhs(i,j) = - 8.d0 * dxx * (         d2uydxdy(i-1,j,k)
     &                                + 5.5d0 * d2uydxdy(i,j,k)
     &                                +         d2uydxdy(i+1,j,k)
     &               + v_kb(k)      * (         wy(i-1,j,k)
     &                                + 5.5d0 * wy(i,j,k)
     &                                +         wy(i+1,j,k) ) )
          end do

          rhs(ptsx-1,j) = - dxx   * (         d2uydxdy(ptsx,j,k)
     &                              + 10.d0 * d2uydxdy(ptsx-1,j,k)
     &                              +         d2uydxdy(ptsx-2,j,k)
     &                  + v_kb(k) * (         wy(ptsx,j,k)
     &                              + 10.d0 * wy(ptsx-1,j,k)
     &                              +         wy(ptsx-2,j,k) ) ) / 12.d0

          rhs(ptsx,j) = ( v_kb(k) * wy(ptsx,j,k) + d2uydxdy(ptsx,j,k) )
     &                / ( dble(k - 1) * dble(k - 1) * beta * beta )
        end do

        call lusolver(lu, rhs)

        do j = 2, jmax
          do i = i_ini, ptsx
            ux(i,j,k) = rhs(i,j)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine poi_ux_fi(duydy)

      ! calculate the new ux velocity using 
      ! the continuity equation for fundamental mode and
      ! Poisson equation for other modes
      ! the argument duydy is calculated here and used in poi_uz sub.
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      include 'comm.fourier'
      include 'comm.fs'
      include 'mpif.h'
      integer i, j, k, i_ini
      real*8 a(imax,5), lu(imax,5), k2b2
      complex*16 rhs(ptsx,jmax), duydy(ptsx,jmax,kfour),
     &           d2uydxdy(ptsx,jmax,kfour)

      call deryfv(duydy, uy)
      i_ini = 1
      if (my_rank .eq. 0) i_ini = 2

      ! for all modes a ux_Poisson equation is used
      call derparx(d2uydxdy, duydy)

      do k = 1, kfour

        k2b2 = dxx * v_k2b2(k) 
        call lhs5(a, k2b2)
        call ludecomp(a, lu)

        do j = 2, jmax

          rhs(1,j) = ux(1,j,k)

          rhs(2,j) = - dxx   * (         d2uydxdy(1,j,k)
     &                         + 10.d0 * d2uydxdy(2,j,k)
     &                         +         d2uydxdy(3,j,k)
     &             + v_kb(k) * (         wy(1,j,k)
     &                         + 10.d0 * wy(2,j,k)
     &                         +         wy(3,j,k) ) ) / 12.d0

          do i = 3, ptsx - 2
            rhs(i,j) = - 8.d0 * dxx * (         d2uydxdy(i-1,j,k)
     &                                + 5.5d0 * d2uydxdy(i,j,k)
     &                                +         d2uydxdy(i+1,j,k)
     &               + v_kb(k)      * (         wy(i-1,j,k)
     &                                + 5.5d0 * wy(i,j,k)
     &                                +         wy(i+1,j,k) ) )
          end do

          rhs(ptsx-1,j) = - dxx   * (         d2uydxdy(ptsx,j,k)
     &                              + 10.d0 * d2uydxdy(ptsx-1,j,k)
     &                              +         d2uydxdy(ptsx-2,j,k)
     &                  + v_kb(k) * (         wy(ptsx,j,k)
     &                              + 10.d0 * wy(ptsx-1,j,k)
     &                              +         wy(ptsx-2,j,k) ) ) / 12.d0

          rhs(ptsx,j) = ux(ptsx,j,k)
        end do

        call lusolver(lu, rhs)

        do j = 2, jmax
          do i = i_ini, ptsx
            ux(i,j,k) = rhs(i,j)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine poi_uz(duydy)

      ! calculate the new uz velocity using the uz Poisson equation
      ! the argument duydy comes from poi_ux subroutine
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      include 'comm.fourier'
      include 'mpif.h'
      integer status(MPI_status_size)
      integer i, j, k, i_ini
      real*8 a(imax,5), lu(imax,5), k2b2
      complex*16 rhs(ptsx,jmax), aux(2),
     &           duydy(ptsx,jmax,kfour), dwydx(ptsx,jmax,kfour)

      i_ini = 1
      if (my_rank .eq. 0) i_ini = 2

      ! for the fundamental mode, the vorticity equation is used
      k = 1
      do j = 2, jmax
        if (my_rank .gt. 0) then
          call MPI_Recv(aux, 2, MPI_COMPLEX16, my_rank - 1, 
     &                  253, MPI_COMM_WORLD, status, ierr)
          uz(1,j,k) = aux(1)
          uz(2,j,k) = aux(2)
         else
          uz(2,j,k) = ( dx * ( 251.d0 * wy(1,j,k) + 646.d0 * wy(2,j,k)
     &                       - 264.d0 * wy(3,j,k) + 106.d0 * wy(4,j,k)
     &                       -  19.d0 * wy(5,j,k) )
     &                + 720.d0 * uz(1,j,k) ) / 720.d0
        end if

        do i = 3, ptsx - 2
          uz(i,j,k) = ( dx * (  281.d0 * wy(i-2,j,k)
     &                       + 2056.d0 * wy(i-1,j,k)
     &                       + 1176.d0 * wy(i,j,k)
     &                       -  104.d0 * wy(i+1,j,k)
     &                       +   11.d0 * wy(i+2,j,k) ) / 90.d0
     &                + 11.d0 * uz(i-2,j,k) + 16.d0 * uz(i-1,j,k) )
     &                / 27.d0
        end do

        if (my_rank .lt. numproc) then
          aux(1) = uz(ptsx-inter,j,k)
          aux(2) = uz(ptsx-inter+1,j,k)
          call MPI_Send(aux, 2, MPI_COMPLEX16, my_rank + 1, 
     &                  253, MPI_COMM_WORLD, ierr)
         else
          i = ptsx - 1
          uz(i,j,k) = ( dx * ( 10.d0 * wy(i-2,j,k)
     &                       + 57.d0 * wy(i-1,j,k)
     &                       + 24.d0 * wy(i,j,k)
     &                       -         wy(i+1,j,k) )
     &                       + 33.d0 * uz(i-2,j,k)
     &                       + 24.d0 * uz(i-1,j,k) ) / 57.d0
          i = ptsx
          uz(i,j,k) = ( dx * ( 251.d0 * wy(i,j,k)
     &                       + 646.d0 * wy(i-1,j,k)
     &                       - 264.d0 * wy(i-2,j,k)
     &                       + 106.d0 * wy(i-3,j,k)
     &                       -  19.d0 * wy(i-4,j,k) )
     &                       + 720.d0 * uz(i-1,j,k) ) / 720.d0
        end if
      end do

      call boundary_exchange_1stmode(uz)

      ! for other modes the Poisson equation is used
      call derparx(dwydx, wy)

      do k = 2, kfour

        k2b2 = dxx * v_k2b2(k)
        call lhs5(a, k2b2)
        call ludecomp(a, lu)

        do j = 2, jmax
          rhs(1,j) = uz(1,j,k)

          rhs(2,j) = dxx     * (         dwydx(1,j,k)
     &                         + 10.d0 * dwydx(2,j,k)
     &                         +         dwydx(3,j,k)
     &             - v_kb(k) * (         duydy(1,j,k)
     &                         + 10.d0 * duydy(2,j,k)
     &                         +         duydy(3,j,k) ) ) / 12.d0

          do i = 3, ptsx - 2
            rhs(i,j) = 8.d0 * dxx * (         dwydx(i-1,j,k)
     &                              + 5.5d0 * dwydx(i,j,k)
     &                              +         dwydx(i+1,j,k)
     &               - v_kb(k)    * (         duydy(i-1,j,k)
     &                              + 5.5d0 * duydy(i,j,k)
     &                              +         duydy(i+1,j,k) ) )
          end do

          rhs(ptsx-1,j) = dxx     * (         dwydx(ptsx,j,k)
     &                              + 10.d0 * dwydx(ptsx-1,j,k)
     &                              +         dwydx(ptsx-2,j,k)
     &                  - v_kb(k) * (         duydy(ptsx,j,k)
     &                              + 10.d0 * duydy(ptsx-1,j,k)
     &                              +         duydy(ptsx-2,j,k) ) )
     &                  / 12.d0

          rhs(ptsx,j) = ( - dwydx(ptsx,j,k)
     &                + v_kb(k) * duydy(ptsx,j,k) )
     &                / ( dble(k - 1) * dble(k - 1) * beta * beta )
        end do

        call lusolver(lu, rhs)

        do j = 2, jmax
          do i = i_ini, ptsx
            uz(i,j,k) = rhs(i,j)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine bzonew

      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      include 'mpif.h'
      integer i, k, bdt
      real*8 bdfc(2,imax), uxb(ptsx,jmax), uyb(ptsx,jmax),
     &       wzb(ptsx,jmax)
      common/blas/ uxb, uyb, wzb
      common/bd/ bdfc

      do k = 1, kfour
        bdt = 1
        if (k .eq. 3) bdt = 2
        do i = 1, ptsx
          wx(i,1,k) = wx(i,1,k) * bdfc(bdt,i+shift)
          wy(i,1,k) = wy(i,1,k) * bdfc(bdt,i+shift)
          wz(i,1,k) = wz(i,1,k) * bdfc(bdt,i+shift)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine boundary_exchange_1d(var)

      implicit none
      include 'par.for'
      include 'comm.par'
      include 'mpif.h'
      integer i
      integer status(MPI_Status_size)
      complex*16 aux(3), var(ptsx)

      if (my_rank .lt. numproc) then
        ! send the penultimate columns
        do i = 1, 3
          aux(i) = var(i + ptsx - inter - 1)
        end do
        call MPI_Send(aux, 3, MPI_COMPLEX16, my_rank + 1, 
     &                131, MPI_COMM_WORLD, ierr)
        ! receive the last columns
        call MPI_Recv(aux, 3, MPI_COMPLEX16, my_rank + 1,
     &                141, MPI_COMM_WORLD, status, ierr)
        do i = 1, 3
          var(i + ptsx - 3) = aux(i)
        end do
      end if

      if (my_rank .gt. 0) then
        ! receive the first columns
        call MPI_Recv(aux, 3, MPI_COMPLEX16, my_rank - 1,
     &                131, MPI_COMM_WORLD, status, ierr)
        do i = 1, 3
          var(i) = aux(i)
        end do
        ! send the second columns
        do i = 1, 3
          aux(i) = var(i + inter - 2)
        end do
        call MPI_Send(aux, 3, MPI_COMPLEX16, my_rank - 1,
     &                141, MPI_COMM_WORLD, ierr)
      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine boundary_exchange_wwall(var)

      implicit none
      include 'par.for'
      include 'comm.par'
      include 'mpif.h'
      integer i, k
      integer status(MPI_Status_size)
      complex*16 aux(3), var(ptsx,jmax,kfour)

      do k = 1, kfour
        if (my_rank .lt. numproc) then
          ! receive the last columns
          call MPI_Recv(aux, 3, MPI_COMPLEX16, my_rank + 1,
     &                  140, MPI_COMM_WORLD, status, ierr)
          do i = 1, 3
            var(i + ptsx - 3, 1, k) = aux(i)
          end do
        end if
        if (my_rank .gt. 0) then
          ! send the second columns
          do i = 1, 3
            aux(i) = var(i + inter - 2, 1, k)
          end do
          call MPI_Send(aux, 3, MPI_COMPLEX16, my_rank - 1,
     &                  140, MPI_COMM_WORLD, ierr)
        end if
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine boundary_exchange_1stmode(var)

      implicit none
      include 'par.for'
      include 'comm.par'
      include 'mpif.h'
      integer i, j
      integer status(MPI_Status_size)
      complex*16 aux(3*jmax), var(ptsx,jmax,kfour)

      if (my_rank .lt. numproc) then
          ! receive the last columns
        call MPI_Recv(aux, 3*jmax, MPI_COMPLEX16, my_rank + 1,
     &                240, MPI_COMM_WORLD, status, ierr)
        do i = 1, 3
          do j = 1, jmax
            var(i + ptsx - 3,j,1) = aux(j+(i-1)*jmax)
          end do
        end do
      end if
      if (my_rank .gt. 0) then
        ! send the second columns
        do i = 1, 3
          do j = 1, jmax
            aux(j+(i-1)*jmax) = var(i + inter - 2,j,1)
          end do
        end do
        call MPI_Send(aux, 3*jmax, MPI_COMPLEX16, my_rank - 1,

     &                240, MPI_COMM_WORLD, ierr)
      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                        end of loop calculations                              c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
