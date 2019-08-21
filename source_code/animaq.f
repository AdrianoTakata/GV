cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                 program that calculates Q for visualization                  c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program isoq 19012009

      implicit none
      include 'par.for'
      include 'comm.fourier'
      integer i, j, k, var
      real*8    duxdxp(imax,jmax,kphys), duxdyp(imax,jmax,kphys),
     &          duydxp(imax,jmax,kphys), duydyp(imax,jmax,kphys),
     &          duzdxp(imax,jmax,kphys), duzdyp(imax,jmax,kphys),
     &          duzdzp(imax,jmax,kphys), duxdzp(imax,jmax,kphys),
     &          duydzp(imax,jmax,kphys),      Q(imax,jmax,kphys),
     &            thtp(imax,jmax,kphys)
      complex*16 duxdx(imax,jmax,kfour),  duxdy(imax,jmax,kfour),
     &           duydx(imax,jmax,kfour),  duydy(imax,jmax,kfour),
     &           duzdx(imax,jmax,kfour),  duzdy(imax,jmax,kfour),
     &           duzdz(imax,jmax,kfour),  duxdz(imax,jmax,kfour),
     &           duydz(imax,jmax,kfour),
     &             uxt(imax,jmax,kfour),    uyt(imax,jmax,kfour),
     &             uzt(imax,jmax,kfour),    tht(imax,jmax,kfour)
      common/iso/ uxt, uyt, uzt, tht

      call derivs_kt

      do var = 1, 16
        call initval(var)

        call derxt(duxdx, uxt)
        call derxt(duydx, uyt)
        call derxt(duzdx, uzt)
        call deryt(duxdy, uxt)
        call deryfvt(duydy, uyt)
        call deryt(duzdy, uzt)

        do k = 1, kfour
          do j = 1, jmax
            do i = 1, imax
              duxdz(i,j,k) = v_kb(k) * uxt(i,j,k)
              duydz(i,j,k) = v_kb(k) * uyt(i,j,k)
              duzdz(i,j,k) = v_kb(k) * uzt(i,j,k)
            end do
          end do
        end do

        call f_to_p(duxdxp, duxdx, imax)
        call f_to_p(duydxp, duydx, imax)
        call f_to_p(duzdxp, duzdx, imax)
        call f_to_p(duxdyp, duxdy, imax)
        call f_to_p(duydyp, duydy, imax)
        call f_to_p(duzdyp, duzdy, imax)
        call f_to_p(duxdzp, duxdz, imax)
        call f_to_p(duydzp, duydz, imax)
        call f_to_p(duzdzp, duzdz, imax)
        call f_to_p(thtp,   tht,   imax)

        do k = 1, kphys
          do j = 1, jmax
            do i = 1, imax
                Q(i,j,k) = - 0.5d0 * (  duxdxp(i,j,k) * duxdxp(i,j,k)
     &                     + duydyp(i,j,k) * duydyp(i,j,k)
     &                     + duzdzp(i,j,k) * duzdzp(i,j,k)
     &                     + 2.d0 * ( duxdyp(i,j,k) * duydxp(i,j,k)
     &                     + duxdzp(i,j,k) * duzdxp(i,j,k)
     &                     + duydzp(i,j,k) * duzdyp(i,j,k) ) )
            end do
          end do
        end do

        call escreveq(Q, thtp, var)
      end do

      stop
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine initval(var)

      implicit none
      include 'par.for'
      character*40 nome
      integer i, j, t, my_rank, k, shift, var
      real*8 uxb(ptsx,jmax), uyb(ptsx,jmax), wzb(ptsx,jmax),
     &       thb(ptsx,jmax)
      complex*16  ux(ptsx,jmax,kfour),  uy(ptsx,jmax,kfour),
     &            uz(ptsx,jmax,kfour),  th(ptsx,jmax,kfour),
     &           uxt(imax,jmax,kfour), uyt(imax,jmax,kfour),
     &           uzt(imax,jmax,kfour), tht(imax,jmax,kfour)
      common/iso/ uxt, uyt, uzt, tht

      ! variables data
      t = (var - 1) * 80 + 12800

      select case (my_form)

        case(0)
          do my_rank = 0, np - 1
            write(nome,'(a,i0.2,a,i0.5,a)')
     &           'dataiso/dataiso_',my_rank,'_',t,'.bin'
            open(2, file = nome, form = 'unformatted')
            read(2) ux, uy, uz
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

        case(1)
          do my_rank = 0, np - 1
            write(nome,'(a,i0.2,a)')'baseflow2D/based_',my_rank,'.bin'
            open (1, file = nome, form = 'unformatted')
            read(1) uxb, uyb, wzb
            close(unit = 1)
            write(nome,'(a,i0.2,a,i0.5,a)')
     &           'dataiso/dataiso_',my_rank,'_',t,'.bin'
            open(2, file = nome, form = 'unformatted')
            read(2) ux, uy, uz
            close (unit = 2)
            shift = my_rank * (ptsx - inter - 1)
            do j = 1, jmax
              do i = 1, ptsx
                uxt(i+shift,j,1) = ux(i,j,1) + uxb(i,j)
                uyt(i+shift,j,1) = uy(i,j,1) + uyb(i,j)
                uzt(i+shift,j,1) = uz(i,j,1)
              end do
            end do
            do k = 2, kfour
              do j = 1, jmax
                do i = 1, ptsx
                  uxt(i+shift,j,k) = ux(i,j,k)
                  uyt(i+shift,j,k) = uy(i,j,k)
                  uzt(i+shift,j,k) = uz(i,j,k)
                end do
              end do
            end do
          end do

        case(2)
          do my_rank = 0, np - 1
            write(nome,'(a,i0.2,a)')'baseflow2D/based_',my_rank,'.bin'
            open (1, file = nome,form = 'unformatted')
            read(1) uxb, uyb, wzb, thb
            close(unit = 1)
            write(nome,'(a,i0.2,a,i0.5,a)')
     &           'dataiso_',my_rank,'_',t,'.bin'
            open(2, file = nome, form = 'unformatted')
            read(2) ux, uy, uz, th
            close (unit = 2)
            shift = my_rank * (ptsx - inter - 1)
            do j = 1, jmax
              do i = 1, ptsx
                uxt(i+shift,j,1) = ux(i,j,1) + uxb(i,j)
                uyt(i+shift,j,1) = uy(i,j,1) + uyb(i,j)
                uzt(i+shift,j,1) = uz(i,j,1)
                tht(i+shift,j,1) = th(i,j,1) + thb(i,j)
              end do
            end do
            do k = 2, kfour
              do j = 1, jmax
                do i = 1, ptsx
                  uxt(i+shift,j,k) = ux(i,j,k)
                  uyt(i+shift,j,k) = uy(i,j,k)
                  uzt(i+shift,j,k) = uz(i,j,k)
                  tht(i+shift,j,k) = th(i,j,k)
                end do
              end do
            end do
          end do
      end select

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine escreveq(Q, thtp, var)

      ! write the results in Fourier modes
      implicit none
      include 'par.for'
      character*30 nome
      integer i, j, k, var
      real*8 x, y, z, Q(imax,jmax,kphys), thtp(imax,jmax,kphys)

      ! writes data to spacial space to be open by tecplot
      write(nome,'(a,i0.2,a)')'isoq_',var,'.dat'
      open (3, file = nome,status = 'unknown')
      write(3,*) 'VARIABLES="x","y","z","Q","th"'
      write(3,*) 'ZONE I=',imax/2+1,',J=',jmax,',K=',kphys+1,',F=POINT'

      z = - dz
      do j = 1, jmax
        y = dble(j-1) * dy
        do i = 1, imax, 2
          x = x0 + dble(i-1) * dx
          if (abs(Q(i,j,kphys)) .lt. 1e-15) Q(i,j,kphys) = 0.d0
          write(3,5) x, y, z, Q(i,j,kphys), thtp(i,j,kphys)
        end do
      end do

      do k = 1, kphys
        z = dble(k-1) * dz
        do j = 1, jmax
          y = dble(j-1) * dy
          do i = 1, imax, 2
            x = x0 + dble(i-1) * dx
            if (abs(Q(i,j,k)) .lt. 1e-15) Q(i,j,k) = 0.d0
            write(3,5) x, y, z, Q(i,j,k), thtp(i,j,k)
          end do
        end do
      end do
      close (unit = 3)

    5 format(1x, 3d14.6, 2d17.9)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine derivs_kt

      ! calculate the tri-diagonal LHS for the derivative calculation
      implicit none
      include 'par.for'
      include 'comm.coef'
      include 'comm.multi'
      real*8  a1x(imax),  b1x(imax),  c1x(imax),
     &        a1y(jmax),  b1y(jmax),  c1y(jmax),
     &       a1fv(jmax), b1fv(jmax), c1fv(jmax)
      common/der1x/  a1x,  b1x,  c1x
      common/der1y/  a1y,  b1y,  c1y
      common/der1fv/ a1fv, b1fv, c1fv

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
      read(1) ! integration in the y direction, used in baseflow2D
      read(1) ! integration in the y direction, used in baseflow2D
      read(1) ! integration in the y direction, used in baseflow2D
      read(1) ! integration in the y direction, used in baseflow2D
      close(unit = 1)

      call create_ctes

      call coefx(a1x, b1x, c1x)
      call coefy(a1y, b1y, c1y)
      call coeffvt(a1fv, b1fv, c1fv)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine derxt(ddx, fc)

      ! first derivative calculation in x direction
      implicit none
      include 'par.for'
      integer i, j, k
      real*8 a1x(imax), b1x(imax), c1x(imax)
      complex*16 rhs(imax), u(imax),
     &           fc(imax,jmax,kfour), ddx(imax,jmax,kfour)
      common/der1x/ a1x, b1x, c1x

      do k = 1, kfour
        do j = 1, jmax
          call rhsxt(fc, rhs, j, k)
          call tridxt(rhs, a1x, b1x, c1x, u)
          do i = 1, imax
            ddx(i,j,k) = u(i)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine deryt(ddy, fc)

      ! first derivative calculation in y direction
      implicit none
      include 'par.for'
      integer i, j, k
      real*8 a1y(jmax), b1y(jmax), c1y(jmax)
      complex*16 rhs(jmax), u(jmax),
     &           fc(imax,jmax,kfour), ddy(imax,jmax,kfour)
      common/der1y/ a1y, b1y, c1y

      do k = 1, kfour
        do i = 1, imax
          call rhsyt(fc, rhs, i, k)
          call tridyt(rhs, a1y, b1y, c1y, u)
          do j = 1, jmax
            ddy(i,j,k) = u(j)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine deryfvt(ddy, fc)

      ! first derivatives calculation in y direction for uy
      implicit none
      include 'par.for'
      integer i, j, k
      real*8 a1fv(jmax), b1fv(jmax), c1fv(jmax)
      complex*16 rhs(jmax), u(jmax),
     &           fc(imax,jmax,kfour), ddy(imax,jmax,kfour)
      common/der1fv/ a1fv, b1fv, c1fv

      do k = 1, kfour
        do i = 1, imax
          call rhsyfvt(fc, rhs, i, k)
          call tridyt(rhs, a1fv, b1fv, c1fv, u)
          do j = 1, jmax
            ddy(i,j,k) = u(j)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsxt(fc, rhs, j, k)

      ! RHS for the first derivative calculation in x direction
      implicit none
      include 'par.for'
      integer i, j, k
      complex*16 rhs(imax), fc(imax,jmax,kfour)

      rhs(1) = ( - 74.d0 * fc(1,j,k) + 16.d0 * fc(2,j,k)
     &           + 72.d0 * fc(3,j,k) - 16.d0 * fc(4,j,k)
     &           +  2.d0 * fc(5,j,k) ) / ( 24.d0 * dx )

      rhs(2) = ( - 406.d0 * fc(1,j,k) - 300.d0 * fc(2,j,k)
     &           + 760.d0 * fc(3,j,k) -  80.d0 * fc(4,j,k)
     &           +  30.d0 * fc(5,j,k) -   4.d0 * fc(6,j,k) ) /
     &           ( 120.d0 * dx )

      do i = 3, imax - 2
        rhs(i) = (           fc(i+2,j,k) - fc(i-2,j,k)
     &           + 28.d0 * ( fc(i+1,j,k) - fc(i-1,j,k) ) ) /
     &           ( 12.d0 * dx )
      end do

      rhs(imax-1) = ( - 406.d0 * fc(imax,j,k)
     &                - 300.d0 * fc(imax-1,j,k)
     &                + 760.d0 * fc(imax-2,j,k)
     &                -  80.d0 * fc(imax-3,j,k)
     &                +  30.d0 * fc(imax-4,j,k)
     &                -   4.d0 * fc(imax-5,j,k) ) /
     &              ( - 120.d0 * dx )

      rhs(imax) = ( - 74.d0 * fc(imax,j,k)   + 16.d0 * fc(imax-1,j,k)
     &              + 72.d0 * fc(imax-2,j,k) - 16.d0 * fc(imax-3,j,k)
     &              +  2.d0 * fc(imax-4,j,k) ) / ( - 24.d0 * dx )

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsyt(fc, rhs, i, k)

      ! RHS for the first derivative calculation in y direction
      implicit none
      include 'par.for'
      include 'comm.coef'
      include 'comm.multi'
      integer i, j, k
      complex*16 rhs(jmax), fc(imax,jmax,kfour)

      rhs(1) = v_qdy(1) * ( fp_fd_coef_e(2) * fc(i,1,k)
     &                    + fp_fd_coef_e(3) * fc(i,2,k)
     &                    + fp_fd_coef_e(4) * fc(i,3,k)
     &                    + fp_fd_coef_e(5) * fc(i,4,k)
     &                    + fp_fd_coef_e(6) * fc(i,5,k) )

      rhs(2) = v_qdy(1) * ( sp_fd_coef(4) * fc(i,1,k)
     &                    + sp_fd_coef(5) * fc(i,2,k)
     &                    + sp_fd_coef(6) * fc(i,3,k)
     &                    + sp_fd_coef(7) * fc(i,4,k)
     &                    + sp_fd_coef(8) * fc(i,5,k)
     &                    + sp_fd_coef(9) * fc(i,6,k) )

      do j = 3, jmax - 2
        rhs(j) = v_qdy(j-2) * ( cp_fd_coef(4) * fc(i,j-2,k)
     &                        + cp_fd_coef(5) * fc(i,j-1,k)
     &                        + cp_fd_coef(6) * fc(i,j,k)
     &                        + cp_fd_coef(7) * fc(i,j+1,k)
     &                        + cp_fd_coef(8) * fc(i,j+2,k) )
      end do

      rhs(jmax-1) = v_qdy(jmax-5)
     &            * ( pp_fd_coef(4) * fc(i,jmax,k)
     &              + pp_fd_coef(5) * fc(i,jmax-1,k)
     &              + pp_fd_coef(6) * fc(i,jmax-2,k)
     &              + pp_fd_coef(7) * fc(i,jmax-3,k)
     &              + pp_fd_coef(8) * fc(i,jmax-4,k)
     &              + pp_fd_coef(9) * fc(i,jmax-5,k) )

      rhs(jmax) = v_qdy(jmax-4) * ( lp_fd_coef(3) * fc(i,jmax,k)
     &                            + lp_fd_coef(4) * fc(i,jmax-1,k)
     &                            + lp_fd_coef(5) * fc(i,jmax-2,k)
     &                            + lp_fd_coef(6) * fc(i,jmax-3,k)
     &                            + lp_fd_coef(7) * fc(i,jmax-4,k) )

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsyfvt(fc, rhs, i, k)

      ! RHS for the first derivative calculation in y direction
      implicit none
      include 'par.for'
      include 'comm.coef'
      include 'comm.multi'
      integer i, j, k
      complex*16 rhs(jmax), fc(imax,jmax,kfour)

      rhs(1) = dcmplx(0.d0,0.d0)

      rhs(2) = v_qdy(1) * ( sp_fd_coef(4) * fc(i,1,k)
     &                    + sp_fd_coef(5) * fc(i,2,k)
     &                    + sp_fd_coef(6) * fc(i,3,k)
     &                    + sp_fd_coef(7) * fc(i,4,k)
     &                    + sp_fd_coef(8) * fc(i,5,k)
     &                    + sp_fd_coef(9) * fc(i,6,k) )

      do j = 3, jmax - 2
        rhs(j) = v_qdy(j-2) * ( cp_fd_coef(4) * fc(i,j-2,k)
     &                        + cp_fd_coef(5) * fc(i,j-1,k)
     &                        + cp_fd_coef(6) * fc(i,j,k)
     &                        + cp_fd_coef(7) * fc(i,j+1,k)
     &                        + cp_fd_coef(8) * fc(i,j+2,k) )
      end do

      rhs(jmax-1) = v_qdy(jmax-5) * ( pp_fd_coef(4) * fc(i,jmax,k)
     &                              + pp_fd_coef(5) * fc(i,jmax-1,k)
     &                              + pp_fd_coef(6) * fc(i,jmax-2,k)
     &                              + pp_fd_coef(7) * fc(i,jmax-3,k)
     &                              + pp_fd_coef(8) * fc(i,jmax-4,k)
     &                              + pp_fd_coef(9) * fc(i,jmax-5,k) )

      rhs(jmax) = v_qdy(jmax-4) * ( lp_fd_coef(3) * fc(i,jmax,k)
     &                            + lp_fd_coef(4) * fc(i,jmax-1,k)
     &                            + lp_fd_coef(5) * fc(i,jmax-2,k)
     &                            + lp_fd_coef(6) * fc(i,jmax-3,k)
     &                            + lp_fd_coef(7) * fc(i,jmax-4,k) )

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tridxt(rr, a, b, c, u)

      ! solves the tridiagonal matrix for derivatives in x direction
      implicit none
      include 'par.for'
      integer i
      real*8 a(imax), b(imax), c(imax), gam(imax), bet
      complex*16 rr(imax), u(imax)

      bet  = b(1)
      u(1) = rr(1) / bet
      do i = 2, imax
        gam(i) = c(i-1) / bet
        bet    = b(i) - a(i) * gam(i)
        u(i)   = ( rr(i) - a(i) * u(i-1) ) / bet
      end do
      do i = imax - 1, 1, -1
        u(i) = u(i) - gam(i+1) * u(i+1)
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tridyt(rr, a, b, c, u)

      ! solves the tridiagonal matrix for derivatives in y direction
      implicit none
      include 'par.for'
      integer j
      real*8 a(jmax), b(jmax), c(jmax), gam(jmax), bet
      complex*16 rr(jmax), u(jmax)

      bet  = b(1)
      u(1) = rr(1) / bet
      do j = 2, jmax
        gam(j) = c(j-1) / bet
        bet    = b(j) - a(j) * gam(j)
        u(j)   = ( rr(j) - a(j) * u(j-1) ) / bet
      end do
      do j = jmax - 1, 1, -1
        u(j) = u(j) - gam(j+1) * u(j+1)
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine coefx(a, b, c)

      ! mount the LHS of the matrix for the first derivative
      implicit none
      include 'par.for'
      integer i
      real*8 a(imax), b(imax), c(imax)

      a(1)      = 0.d0
      b(1)      = 1.d0
      c(1)      = 4.d0

      a(2)      = 1.d0
      b(2)      = 6.d0
      c(2)      = 2.d0

      do i = 3, imax - 2
        a(i)    = 1.d0
        b(i)    = 3.d0
        c(i)    = 1.d0
      end do

      a(imax-1) = 2.d0
      b(imax-1) = 6.d0
      c(imax-1) = 1.d0

      a(imax)   = 4.d0
      b(imax)   = 1.d0
      c(imax)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine coefy(a, b, c)

      ! mount the LHS of the matrix for the first derivative
      implicit none
      include 'par.for'
      include 'comm.coef'
      integer j
      real*8 a(jmax), b(jmax), c(jmax)

      a(1)      = 0.d0
      b(1)      = fp_fd_coef_e(1)
      c(1)      = 0.d0

      a(2)      = sp_fd_coef(1)
      b(2)      = sp_fd_coef(2)
      c(2)      = sp_fd_coef(3)

      do j = 3, jmax - 2
        a(j)    = cp_fd_coef(1)
        b(j)    = cp_fd_coef(2)
        c(j)    = cp_fd_coef(3)
      end do

      a(jmax-1) = pp_fd_coef(3)
      b(jmax-1) = pp_fd_coef(2)
      c(jmax-1) = pp_fd_coef(1)

      a(jmax)   = lp_fd_coef(2)
      b(jmax)   = lp_fd_coef(1)
      c(jmax)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine coeffvt(a, b, c)

      implicit none
      include 'par.for'
      include 'comm.coef'
      integer j
      real*8 a(jmax), b(jmax), c(jmax)

      a(1)      = 0.d0
      b(1)      = 1.d0
      c(1)      = 0.d0

      a(2)      = sp_fd_coef(1)
      b(2)      = sp_fd_coef(2)
      c(2)      = sp_fd_coef(3)

      do j = 3, jmax - 2
        a(j)    = cp_fd_coef(1)
        b(j)    = cp_fd_coef(2)
        c(j)    = cp_fd_coef(3)
      end do

      a(jmax-1) = pp_fd_coef(3)
      b(jmax-1) = pp_fd_coef(2)
      c(jmax-1) = pp_fd_coef(1)

      a(jmax)   = lp_fd_coef(2)
      b(jmax)   = lp_fd_coef(1)
      c(jmax)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine create_ctes

      implicit none
      include 'par.for'
      include 'comm.multi'
      include 'comm.fourier'
 
      integer lvl, j, k
      real*8 aux
 
      ! Multigrid spatial calculations
 
      ! dy0 at each multigrid level
      do lvl = 1 , msh
        if (stf .ne. 1.d0) then
         v_dy0(lvl) = dy * (( stf**(2**(lvl-1)) - 1.d0) / (stf - 1.d0))
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
          do j = 1 , ( jmax - 1 ) / 2**(lvl-1)
            aux           = (v_dy0(lvl) * v_stf(lvl)**(j-1))
            v_qdy2(j,lvl) = 1.d0 / (aux**2)
          end do
         else
          do j = 1 , ( jmax - 1 ) / 2**(lvl-1)
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
        v_k2b2(k) = - dble(k-1) * dble(k-1) * beta * beta
        v_kb(k)   = - im * dble(k-1) * beta
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c              end program that calculates Q for visualization                 c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
