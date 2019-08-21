cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                         derivative calculations                              c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine derivs_k

      ! mounts the tri-diagonal LHS for derivative calculations
      implicit none
      include 'par.for'
      real*8  a1x(ptsx),  b1x(ptsx),  c1x(ptsx),
     &        a2x(ptsx),  b2x(ptsx),  c2x(ptsx),
     &        a1y(jmax),  b1y(jmax),  c1y(jmax),
     &        a2y(jmax),  b2y(jmax),  c2y(jmax),
     &       a1fv(jmax), b1fv(jmax), c1fv(jmax)
      common/der1x/  a1x,  b1x,  c1x
      common/der2x/  a2x,  b2x,  c2x
      common/der1y/  a1y,  b1y,  c1y
      common/der2y/  a2y,  b2y,  c2y
      common/der1fv/ a1fv, b1fv, c1fv

      call coefx(a1x, b1x, c1x)
      call coeffx(a2x, b2x, c2x)
      call coefy(a1y, b1y, c1y)
      call coefy_fv(a1fv, b1fv, c1fv)
      call coeffy(a2y, b2y, c2y)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine derparxue(ddx, fc)

      ! first derivatives calculation in x direction (parallel)
      implicit none
      include 'par.for'
      real*8 a1x(ptsx), b1x(ptsx), c1x(ptsx), fc(ptsx), ddx(ptsx)
      common/der1x/ a1x, b1x, c1x

      call rhsxue(fc, ddx)
      call tridparxue(a1x, b1x, c1x, ddx)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine derparxr(ddx, fc)

      ! first derivatives calculation in x direction
      implicit none
      include 'par.for'
      real*8 a1x(ptsx), b1x(ptsx), c1x(ptsx),
     &       fc(ptsx,jmax), ddx(ptsx,jmax)
      common/der1x/ a1x, b1x, c1x

      call rhsxr(fc, ddx)
      call tridparxr(a1x, b1x, c1x, ddx)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine derparxr3d(ddx, fc)

      ! first derivatives calculation in x direction
      implicit none
      include 'par.for'
      real*8 a1x(ptsx), b1x(ptsx), c1x(ptsx),
     &       fc(ptsx,jmax,kphys), ddx(ptsx,jmax,kphys)
      common/der1x/ a1x, b1x, c1x

      call rhsxr3d(fc, ddx)
      call tridparxr3d(a1x, b1x, c1x, ddx)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine derparx(ddx, fc)

      ! first derivatives calculation in x direction
      implicit none
      include 'par.for'
      real*8 a1x(ptsx), b1x(ptsx), c1x(ptsx)
      complex*16 fc(ptsx,jmax,kfour), ddx(ptsx,jmax,kfour)
      common/der1x/ a1x, b1x, c1x

      call rhsx(fc, ddx)
c     call tridparx(a1x, b1x, c1x, ddx)
      call tridseqx(a1x, b1x, c1x, ddx)
      call boundary_exchange_derivs(ddx)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine derparxx(d2dx2, fc)

      ! second derivative calculations in x direction
      implicit none
      include 'par.for'
      real*8 a2x(ptsx), b2x(ptsx), c2x(ptsx)
      complex*16 fc(ptsx,jmax,kfour), d2dx2(ptsx,jmax,kfour)
      common/der2x/ a2x, b2x, c2x

      call rhsxx(fc, d2dx2)
c     call tridparx(a2x, b2x, c2x, d2dx2)
      call tridseqx(a2x, b2x, c2x, d2dx2)
      call boundary_exchange_derivs(d2dx2)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dery(ddy, fc)

      ! first derivative calculation in y direction
      implicit none
      include 'par.for'
      real*8 a1y(jmax), b1y(jmax), c1y(jmax)
      complex*16 fc(ptsx,jmax,kfour), ddy(ptsx,jmax,kfour)
      common/der1y/ a1y, b1y, c1y

      call rhsy(fc, ddy)
      call tridy(a1y, b1y, c1y, ddy)
      
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine deryr3d(ddy, fc)

      ! first derivative calculation in y direction
      implicit none
      include 'par.for'
      real*8 a1y(jmax), b1y(jmax), c1y(jmax),
     &       fc(ptsx,jmax,kphys), ddy(ptsx,jmax,kphys)
      common/der1y/ a1y, b1y, c1y

      call rhsyr3d(fc, ddy)
      call tridyr3d(a1y, b1y, c1y, ddy)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine deryfv(ddy, fc)

      ! first derivative calculation in y direction for uy
      implicit none
      include 'par.for'
      real*8 a1fv(jmax), b1fv(jmax), c1fv(jmax)
      complex*16 fc(ptsx,jmax,kfour), ddy(ptsx,jmax,kfour)
      common/der1fv/ a1fv, b1fv, c1fv

      call rhsyfv(fc, ddy)
      call tridy(a1fv, b1fv, c1fv, ddy)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine deryy(d2dy2, fc)

      ! second derivative calculation in y direction
      implicit none
      include 'par.for'
      real*8 a2y(jmax), b2y(jmax), c2y(jmax)
      complex*16 fc(ptsx,jmax,kfour), d2dy2(ptsx,jmax,kfour)
      common/der2y/ a2y, b2y, c2y

      call rhsyy(fc, d2dy2)
      call tridy(a2y, b2y, c2y, d2dy2)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsxue(fc, rhs)
      
      ! RHS for the first derivative calculation in x direction
      implicit none
      include 'par.for'
      integer i
      real*8 rhs(ptsx), fc(ptsx)

      rhs(1) = ( - 25.d0 * fc(1) + 48.d0 * fc(2)
     &           - 36.d0 * fc(3) + 16.d0 * fc(4)
     &           -  3.d0 * fc(5) ) / ( 12.d0 * dx )
      
      rhs(2) = ( - 406.d0 * fc(1) - 300.d0 * fc(2)
     &           + 760.d0 * fc(3) -  80.d0 * fc(4)
     &           +  30.d0 * fc(5) -   4.d0 * fc(6) )
     &       /   ( 120.d0 * dx )
      
      do i = 3, ptsx - 2
        rhs(i) = (           fc(i+2) - fc(i-2)
     &           + 28.d0 * ( fc(i+1) - fc(i-1) ) )
     &         / ( 12.d0 * dx )
      end do
      
      rhs(ptsx-1) = ( - 406.d0 * fc(ptsx)
     &                - 300.d0 * fc(ptsx-1)
     &                + 760.d0 * fc(ptsx-2)
     &                -  80.d0 * fc(ptsx-3)
     &                +  30.d0 * fc(ptsx-4)
     &                -   4.d0 * fc(ptsx-5) )
     &            / ( - 120.d0 * dx )
      
      rhs(ptsx) = ( - 74.d0 * fc(ptsx)
     &              + 16.d0 * fc(ptsx-1)
     &              + 72.d0 * fc(ptsx-2)
     &              - 16.d0 * fc(ptsx-3)
     &              +  2.d0 * fc(ptsx-4) )
     &          / ( - 24.d0 * dx )

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsxr(fc, rhs)
      
      ! RHS for the first derivative calculation in x direction
      implicit none
      include 'par.for'
      integer i, j
      real*8 rhs(ptsx,jmax), fc(ptsx,jmax)

      do j = 1, jmax
        rhs(1,j) = ( - 25.d0 * fc(1,j) + 48.d0 * fc(2,j) 
     &               - 36.d0 * fc(3,j) + 16.d0 * fc(4,j) 
     &               -  3.d0 * fc(5,j) ) / ( 12.d0 * dx )
        
        rhs(2,j) = ( - 406.d0 * fc(1,j) - 300.d0 * fc(2,j)
     &               + 760.d0 * fc(3,j) -  80.d0 * fc(4,j)
     &               +  30.d0 * fc(5,j) -   4.d0 * fc(6,j) )
     &           / (   120.d0 * dx )
        
        do i = 3, ptsx - 2
          rhs(i,j) = (           fc(i+2,j) - fc(i-2,j)
     &               + 28.d0 * ( fc(i+1,j) - fc(i-1,j) ) )
     &             / ( 12.d0 * dx )
        end do
        
        rhs(ptsx-1,j) = ( - 406.d0 * fc(ptsx,j)
     &                    - 300.d0 * fc(ptsx-1,j)
     &                    + 760.d0 * fc(ptsx-2,j)
     &                    -  80.d0 * fc(ptsx-3,j)
     &                    +  30.d0 * fc(ptsx-4,j)
     &                    -   4.d0 * fc(ptsx-5,j) )
     &                / ( - 120.d0 * dx )
        
        rhs(ptsx,j) = ( - 74.d0 * fc(ptsx,j)   + 16.d0 * fc(ptsx-1,j)
     &                  + 72.d0 * fc(ptsx-2,j) - 16.d0 * fc(ptsx-3,j)
     &                  +  2.d0 * fc(ptsx-4,j) ) / ( - 24.d0 * dx )
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsx(fc, rhs)
      
      ! RHS for the first derivative calculation in x direction
      implicit none
      include 'par.for'
      integer i, j, k
      complex*16 rhs(ptsx,jmax,kfour), fc(ptsx,jmax,kfour)

      do k = 1, kfour
        do j = 1, jmax
          rhs(1,j,k) = ( - 25.d0 * fc(1,j,k) + 48.d0 * fc(2,j,k)
     &                   - 36.d0 * fc(3,j,k) + 16.d0 * fc(4,j,k)
     &                   -  3.d0 * fc(5,j,k) ) / ( 12.d0 * dx )
          
          rhs(2,j,k) = ( - 406.d0 * fc(1,j,k) - 300.d0 * fc(2,j,k)
     &                   + 760.d0 * fc(3,j,k) -  80.d0 * fc(4,j,k)
     &                   +  30.d0 * fc(5,j,k) -   4.d0 * fc(6,j,k) )
     &               / (   120.d0 * dx )
          
          do i = 3, ptsx - 2
            rhs(i,j,k) = (           fc(i+2,j,k) - fc(i-2,j,k)
     &                   + 28.d0 * ( fc(i+1,j,k) - fc(i-1,j,k) ) )
     &                 / ( 12.d0 * dx )
          end do
          
          rhs(ptsx-1,j,k) = ( - 406.d0 * fc(ptsx,j,k)
     &                        - 300.d0 * fc(ptsx-1,j,k)
     &                        + 760.d0 * fc(ptsx-2,j,k)
     &                        -  80.d0 * fc(ptsx-3,j,k)
     &                        +  30.d0 * fc(ptsx-4,j,k)
     &                        -   4.d0 * fc(ptsx-5,j,k) )
     &                    / ( - 120.d0 * dx )
          
          rhs(ptsx,j,k) = ( - 74.d0 * fc(ptsx,j,k)
     &                      + 16.d0 * fc(ptsx-1,j,k)
     &                      + 72.d0 * fc(ptsx-2,j,k)
     &                      - 16.d0 * fc(ptsx-3,j,k)
     &                      +  2.d0 * fc(ptsx-4,j,k) )
     &                  / ( - 24.d0 * dx )
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsxr3d(fc, rhs)
      
      ! RHS for the first derivative calculation in x direction
      implicit none
      include 'par.for'
      integer i, j, k
      real*8 rhs(ptsx,jmax,kphys), fc(ptsx,jmax,kphys)

      do k = 1, kphys
        do j = 1, jmax
          rhs(1,j,k) = ( - 25.d0 * fc(1,j,k) + 48.d0 * fc(2,j,k) 
     &                   - 36.d0 * fc(3,j,k) + 16.d0 * fc(4,j,k) 
     &                   -  3.d0 * fc(5,j,k) ) / ( 12.d0 * dx )
          
          rhs(2,j,k) = ( - 406.d0 * fc(1,j,k) - 300.d0 * fc(2,j,k)
     &                   + 760.d0 * fc(3,j,k) -  80.d0 * fc(4,j,k)
     &                   +  30.d0 * fc(5,j,k) -   4.d0 * fc(6,j,k) )
     &               / (   120.d0 * dx )
          
          do i = 3, ptsx - 2
            rhs(i,j,k) = (           fc(i+2,j,k) - fc(i-2,j,k)
     &                   + 28.d0 * ( fc(i+1,j,k) - fc(i-1,j,k) ) )
     &                 / ( 12.d0 * dx )
          end do
          
          rhs(ptsx-1,j,k) = ( - 406.d0 * fc(ptsx,j,k)
     &                        - 300.d0 * fc(ptsx-1,j,k)
     &                        + 760.d0 * fc(ptsx-2,j,k)
     &                        -  80.d0 * fc(ptsx-3,j,k)
     &                        +  30.d0 * fc(ptsx-4,j,k)
     &                        -   4.d0 * fc(ptsx-5,j,k) )
     &                    / ( - 120.d0 * dx )
          
          rhs(ptsx,j,k) = ( - 74.d0 * fc(ptsx,j,k)
     &                      + 16.d0 * fc(ptsx-1,j,k)
     &                      + 72.d0 * fc(ptsx-2,j,k)
     &                      - 16.d0 * fc(ptsx-3,j,k)
     &                      +  2.d0 * fc(ptsx-4,j,k) )
     &                  / ( - 24.d0 * dx )
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsxx(fc, rhs)
      
      ! RHS for the second derivative calculation in x direction
      implicit none
      include 'par.for'
      integer i, j, k
      complex*16 rhs(ptsx,jmax,kfour), fc(ptsx,jmax,kfour)

      do k = 1, kfour
        do j = 1, jmax
          rhs(1,j,k) = (  9775.d0 * fc(1,j,k) - 20285.d0 * fc(2,j,k)
     &                 + 11170.d0 * fc(3,j,k) -   550.d0 * fc(4,j,k)
     &                 -   145.d0 * fc(5,j,k) +    35.d0 * fc(6,j,k) )
     &               / ( 60.d0 * dxx )
          
          rhs(2,j,k) = ( 4834.d0 * fc(1,j,k) - 8424.d0 * fc(2,j,k)
     &                 + 1890.d0 * fc(3,j,k) + 2320.d0 * fc(4,j,k)
     &                 -  810.d0 * fc(5,j,k) +  216.d0 * fc(6,j,k)
     &                 -   26.d0 * fc(7,j,k) ) / ( 360.d0 * dxx )
          
          do i = 3, ptsx - 2
            rhs(i,j,k) = (   3.d0 * ( fc(i+2,j,k) + fc(i-2,j,k) )
     &                   +  48.d0 * ( fc(i+1,j,k) + fc(i-1,j,k) )
     &                   - 102.d0 *   fc(i,j,k) ) / ( 4.d0 * dxx )
          end do
          
         rhs(ptsx-1,j,k) = ( 4834.d0 * fc(ptsx,j,k)
     &                     - 8424.d0 * fc(ptsx-1,j,k)
     &                     + 1890.d0 * fc(ptsx-2,j,k)
     &                     + 2320.d0 * fc(ptsx-3,j,k)
     &                     -  810.d0 * fc(ptsx-4,j,k)
     &                     +  216.d0 * fc(ptsx-5,j,k)
     &                     -   26.d0 * fc(ptsx-6,j,k) )
     &                   / (  360.d0 * dxx )
          
          rhs(ptsx,j,k) = (  9775.d0 * fc(ptsx,j,k)
     &                    - 20285.d0 * fc(ptsx-1,j,k)
     &                    + 11170.d0 * fc(ptsx-2,j,k)
     &                    -   550.d0 * fc(ptsx-3,j,k)
     &                    -   145.d0 * fc(ptsx-4,j,k)
     &                    +    35.d0 * fc(ptsx-5,j,k) )
     &                  / (    60.d0 * dxx )
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsy(fc,rhs)

      ! RHS for the first derivative calculation in y direction
      implicit none
      include 'par.for'
      include 'comm.coef'
      include 'comm.multi'
      integer i, j, k
      complex*16 rhs(ptsx,jmax,kfour), fc(ptsx,jmax,kfour)

      do k = 1, kfour
        do i = 1, ptsx
          rhs(i,1,k) = v_qdy(1) * ( fp_fd_coef_e(2) * fc(i,1,k)
     &                            + fp_fd_coef_e(3) * fc(i,2,k)
     &                            + fp_fd_coef_e(4) * fc(i,3,k)
     &                            + fp_fd_coef_e(5) * fc(i,4,k)
     &                            + fp_fd_coef_e(6) * fc(i,5,k) )
          
          rhs(i,2,k) = v_qdy(1) * ( sp_fd_coef(4) * fc(i,1,k)
     &                            + sp_fd_coef(5) * fc(i,2,k)
     &                            + sp_fd_coef(6) * fc(i,3,k)
     &                            + sp_fd_coef(7) * fc(i,4,k)
     &                            + sp_fd_coef(8) * fc(i,5,k)
     &                            + sp_fd_coef(9) * fc(i,6,k) )

         do j = 3, jmax - 2
           rhs(i,j,k) = v_qdy(j-2) * ( cp_fd_coef(4) * fc(i,j-2,k)
     &                               + cp_fd_coef(5) * fc(i,j-1,k)
     &                               + cp_fd_coef(6) * fc(i,j,k)
     &                               + cp_fd_coef(7) * fc(i,j+1,k)
     &                               + cp_fd_coef(8) * fc(i,j+2,k) )
         end do
         
         rhs(i,jmax-1,k) = v_qdy(jmax-5)
     &                   * ( pp_fd_coef(4) * fc(i,jmax,k)
     &                     + pp_fd_coef(5) * fc(i,jmax-1,k)
     &                     + pp_fd_coef(6) * fc(i,jmax-2,k)
     &                     + pp_fd_coef(7) * fc(i,jmax-3,k)
     &                     + pp_fd_coef(8) * fc(i,jmax-4,k)
     &                     + pp_fd_coef(9) * fc(i,jmax-5,k) )
         
         rhs(i,jmax,k) = v_qdy(jmax-4)
     &                 * ( lp_fd_coef(3) * fc(i,jmax,k)
     &                   + lp_fd_coef(4) * fc(i,jmax-1,k)
     &                   + lp_fd_coef(5) * fc(i,jmax-2,k)
     &                   + lp_fd_coef(6) * fc(i,jmax-3,k)
     &                   + lp_fd_coef(7) * fc(i,jmax-4,k) )
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsy_i0(fc,rhs)

      ! RHS for the first derivative calculation in y direction at x = 1
      ! usado no CASE 4
      implicit none
      include 'par.for'
      include 'comm.coef'
      include 'comm.multi'
      integer j, k
      complex*16 rhs(jmax,kfour), fc(jmax,kfour)

      do k = 1, kfour
          rhs(1,k) =   v_qdy(1) * ( fp_fd_coef_e(2) * fc(1,k)
     &                            + fp_fd_coef_e(3) * fc(2,k)
     &                            + fp_fd_coef_e(4) * fc(3,k)
     &                            + fp_fd_coef_e(5) * fc(4,k)
     &                            + fp_fd_coef_e(6) * fc(5,k) )
          
          rhs(2,k) =   v_qdy(1) * ( sp_fd_coef(4) * fc(1,k)
     &                            + sp_fd_coef(5) * fc(2,k)
     &                            + sp_fd_coef(6) * fc(3,k)
     &                            + sp_fd_coef(7) * fc(4,k)
     &                            + sp_fd_coef(8) * fc(5,k)
     &                            + sp_fd_coef(9) * fc(6,k) )

         do j = 3, jmax - 2
           rhs(j,k) =   v_qdy(j-2) * ( cp_fd_coef(4) * fc(j-2,k)
     &                               + cp_fd_coef(5) * fc(j-1,k)
     &                               + cp_fd_coef(6) * fc(j,k)
     &                               + cp_fd_coef(7) * fc(j+1,k)
     &                               + cp_fd_coef(8) * fc(j+2,k) )
         end do
         
         rhs(jmax-1,k) = v_qdy(jmax-5) *
     &                     ( pp_fd_coef(4) * fc(jmax,k)
     &                     + pp_fd_coef(5) * fc(jmax-1,k)
     &                     + pp_fd_coef(6) * fc(jmax-2,k)
     &                     + pp_fd_coef(7) * fc(jmax-3,k)
     &                     + pp_fd_coef(8) * fc(jmax-4,k)
     &                     + pp_fd_coef(9) * fc(jmax-5,k) )
         
         rhs(jmax,k) = v_qdy(jmax-4) *
     &                   ( lp_fd_coef(3) * fc(jmax,k)
     &                   + lp_fd_coef(4) * fc(jmax-1,k)
     &                   + lp_fd_coef(5) * fc(jmax-2,k)
     &                   + lp_fd_coef(6) * fc(jmax-3,k)
     &                   + lp_fd_coef(7) * fc(jmax-4,k) )
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsyr3d(fc, rhs)

      ! RHS for the first derivative calculation in y direction
      implicit none
      include 'par.for'
      include 'comm.coef'
      include 'comm.multi'
      integer i, j, k
      real*8 rhs(ptsx,jmax,kphys), fc(ptsx,jmax,kphys)

      do k = 1, kphys
        do i = 1, ptsx
          rhs(i,1,k) = v_qdy(1) * ( fp_fd_coef_e(2) * fc(i,1,k)
     &                            + fp_fd_coef_e(3) * fc(i,2,k)
     &                            + fp_fd_coef_e(4) * fc(i,3,k)
     &                            + fp_fd_coef_e(5) * fc(i,4,k)
     &                            + fp_fd_coef_e(6) * fc(i,5,k) )
          
          rhs(i,2,k) = v_qdy(1) * ( sp_fd_coef(4) * fc(i,1,k)
     &                            + sp_fd_coef(5) * fc(i,2,k)
     &                            + sp_fd_coef(6) * fc(i,3,k)
     &                            + sp_fd_coef(7) * fc(i,4,k)
     &                            + sp_fd_coef(8) * fc(i,5,k)
     &                            + sp_fd_coef(9) * fc(i,6,k) )

         do j = 3, jmax - 2
           rhs(i,j,k) = v_qdy(j-2) * ( cp_fd_coef(4) * fc(i,j-2,k)
     &                               + cp_fd_coef(5) * fc(i,j-1,k)
     &                               + cp_fd_coef(6) * fc(i,j,k)  
     &                               + cp_fd_coef(7) * fc(i,j+1,k)
     &                               + cp_fd_coef(8) * fc(i,j+2,k) )
         end do
         
         rhs(i,jmax-1,k) = v_qdy(jmax-5)
     &                   * ( pp_fd_coef(4) * fc(i,jmax,k)
     &                     + pp_fd_coef(5) * fc(i,jmax-1,k)
     &                     + pp_fd_coef(6) * fc(i,jmax-2,k)
     &                     + pp_fd_coef(7) * fc(i,jmax-3,k)
     &                     + pp_fd_coef(8) * fc(i,jmax-4,k)
     &                     + pp_fd_coef(9) * fc(i,jmax-5,k) )
         
         rhs(i,jmax,k) = v_qdy(jmax-4)
     &                 * ( lp_fd_coef(3) * fc(i,jmax,k)
     &                   + lp_fd_coef(4) * fc(i,jmax-1,k)
     &                   + lp_fd_coef(5) * fc(i,jmax-2,k)
     &                   + lp_fd_coef(6) * fc(i,jmax-3,k)
     &                   + lp_fd_coef(7) * fc(i,jmax-4,k) )
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsyfv(fc, rhs)

      ! RHS for the first derivative calculation in y direction
      implicit none
      include 'par.for'
      include 'comm.coef'
      include 'comm.multi'
      integer i, j, k
      complex*16 rhs(ptsx,jmax,kfour), fc(ptsx,jmax,kfour)

      do k = 1, kfour
        do i = 1, ptsx
          rhs(i,1,k) = dcmplx(0.d0,0.d0)
          
          rhs(i,2,k) = v_qdy(1) * ( sp_fd_coef(4) * fc(i,1,k)
     &                            + sp_fd_coef(5) * fc(i,2,k)
     &                            + sp_fd_coef(6) * fc(i,3,k)
     &                            + sp_fd_coef(7) * fc(i,4,k)
     &                            + sp_fd_coef(8) * fc(i,5,k)
     &                            + sp_fd_coef(9) * fc(i,6,k) )

         do j = 3, jmax - 2
           rhs(i,j,k) = v_qdy(j-2) * ( cp_fd_coef(4) * fc(i,j-2,k)
     &                               + cp_fd_coef(5) * fc(i,j-1,k)
     &                               + cp_fd_coef(6) * fc(i,j,k)  
     &                               + cp_fd_coef(7) * fc(i,j+1,k)
     &                               + cp_fd_coef(8) * fc(i,j+2,k) )
         end do
         
         rhs(i,jmax-1,k) = v_qdy(jmax-5)
     &                   * ( pp_fd_coef(4) * fc(i,jmax,k)
     &                     + pp_fd_coef(5) * fc(i,jmax-1,k)
     &                     + pp_fd_coef(6) * fc(i,jmax-2,k)
     &                     + pp_fd_coef(7) * fc(i,jmax-3,k)
     &                     + pp_fd_coef(8) * fc(i,jmax-4,k)
     &                     + pp_fd_coef(9) * fc(i,jmax-5,k) )
         
         rhs(i,jmax,k) = v_qdy(jmax-4)
     &                 * ( lp_fd_coef(3) * fc(i,jmax,k)
     &                   + lp_fd_coef(4) * fc(i,jmax-1,k)
     &                   + lp_fd_coef(5) * fc(i,jmax-2,k)
     &                   + lp_fd_coef(6) * fc(i,jmax-3,k)
     &                   + lp_fd_coef(7) * fc(i,jmax-4,k) )
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsyy(fc, rhs)

      ! RHS for the second derivative calculation in y direction
      implicit none
      include 'par.for'
      include 'comm.coef'
      include 'comm.multi'
      integer i, j, k
      complex*16 rhs(ptsx,jmax,kfour), fc(ptsx,jmax,kfour)

      do k = 1, kfour
        do i = 1, ptsx

          rhs(i,1,k) = v_qdy2(1,1)
     &               * ( fp_sd_coef(3) * fc(i,1,k)
     &                 + fp_sd_coef(4) * fc(i,2,k)
     &                 + fp_sd_coef(5) * fc(i,3,k)
     &                 + fp_sd_coef(6) * fc(i,4,k)
     &                 + fp_sd_coef(7) * fc(i,5,k)
     &                 + fp_sd_coef(8) * fc(i,6,k) )
          
          rhs(i,2,k) = v_qdy2(1,1)
     &               * ( sp_sd_coef(4)  * fc(i,1,k)
     &                 + sp_sd_coef(5)  * fc(i,2,k)
     &                 + sp_sd_coef(6)  * fc(i,3,k)
     &                 + sp_sd_coef(7)  * fc(i,4,k)
     &                 + sp_sd_coef(8)  * fc(i,5,k)
     &                 + sp_sd_coef(9)  * fc(i,6,k)
     &                 + sp_sd_coef(10) * fc(i,7,k) )
          
         do j = 3, jmax - 2
           rhs(i,j,k) = v_qdy2(j-2,1)
     &                * ( cp_sd_coef(4) * fc(i,j-2,k)
     &                  + cp_sd_coef(5) * fc(i,j-1,k)
     &                  + cp_sd_coef(6) * fc(i,j,k)
     &                  + cp_sd_coef(7) * fc(i,j+1,k)
     &                  + cp_sd_coef(8) * fc(i,j+2,k) )
         end do
          
        rhs(i,jmax-1,k) = v_qdy2(jmax-6,1)
     &                  * ( pp_sd_coef(4)  * fc(i,jmax,k)
     &                    + pp_sd_coef(5)  * fc(i,jmax-1,k)
     &                    + pp_sd_coef(6)  * fc(i,jmax-2,k)
     &                    + pp_sd_coef(7)  * fc(i,jmax-3,k)
     &                    + pp_sd_coef(8)  * fc(i,jmax-4,k)
     &                    + pp_sd_coef(9)  * fc(i,jmax-5,k)
     &                    + pp_sd_coef(10) * fc(i,jmax-6,k) )
         
        rhs(i,jmax,k) = v_qdy2(jmax-5,1)
     &                * ( lp_sd_coef(3) * fc(i,jmax,k)
     &                  + lp_sd_coef(4) * fc(i,jmax-1,k)
     &                  + lp_sd_coef(5) * fc(i,jmax-2,k)
     &                  + lp_sd_coef(6) * fc(i,jmax-3,k)
     &                  + lp_sd_coef(7) * fc(i,jmax-4,k)
     &                  + lp_sd_coef(8) * fc(i,jmax-5,k) )
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tridparxue(a, b, c, rhs)
 
      ! solves tridiagonal matrix for the derivatives in x direction
      ! domain decomposition in x direction -> parallel subroutine
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'mpif.h'
      integer status(MPI_status_size)
      integer i, i_ini, i_fim
      real*8 a(ptsx), b(ptsx), c(ptsx), aux(2*inter), rhs(ptsx),
     &       u(ptsx), gam(ptsx), bet
 
      if (my_rank .eq. 0) then
        bet    = b(1)
        u(1)   = rhs(1) / bet
        gam(2) = c(1) / bet
        bet    = b(2) - a(2) * gam(2)
        u(2)   = ( rhs(2) - a(2) * u(1) ) / bet
       else
        call MPI_Recv(aux, 2*inter, mpi_double_precision, my_rank-1,
     &                57, MPI_COMM_WORLD, status, ierr)
        do i = 1, inter
          u(i)   = aux(i)
          gam(i) = aux(i+inter)
        end do
      end if
      i_ini = inter
      if (my_rank .eq. 0) i_ini = 3
      i_fim = ptsx - 2
      if (my_rank .eq. numproc) i_fim = ptsx
 
      bet = b(i_ini-1) - a(i_ini-1) * gam(i_ini-1)
      do i = i_ini, i_fim
        gam(i) = c(i-1) / bet
        bet    = b(i) - a(i) * gam(i)
        u(i)   = ( rhs(i) - a(i) * u(i-1) ) / bet
      end do
 
      if (my_rank .lt. numproc) then
        do i = 1, inter
          aux(i)       = u(ptsx-inter+i-1)
          aux(i+inter) = gam(ptsx-inter+i-1)
        end do
        call MPI_Send(aux, 2*inter, mpi_double_precision, my_rank+1,
     &                57, MPI_COMM_WORLD, ierr)
      end if

      if (my_rank .lt. numproc) then
        call MPI_Recv(aux, 2*inter, mpi_double_precision, my_rank+1,
     &                67, MPI_COMM_WORLD, status, ierr)
        do i = 1, inter
          u(ptsx-inter+i) = aux(i)
          gam(ptsx-inter+i) = aux(i+inter)
        end do
      end if
 
      i_fim = ptsx - inter
      if (my_rank .eq. numproc) i_fim = ptsx - 1
      do i = i_fim, 1, -1
        u(i) = u(i) - gam(i+1) * u(i+1)
      end do

      if (my_rank .gt. 0) then
        do i = 1, inter
          aux(i)       = u(i+1)
          aux(i+inter) = gam(i+1)
        end do
        call MPI_Send(aux, 2*inter, mpi_double_precision, my_rank-1,
     &                67, MPI_COMM_WORLD, ierr)
      end if
      do i = 1, ptsx
        rhs(i) = u(i)
      end do

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tridparxr(a, b, c, rhs)
      
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'mpif.h'
      integer status(MPI_status_size)
      integer i, j, i_ini, i_fim
      real*8 a(ptsx), b(ptsx), c(ptsx), aux(2*inter), rhs(ptsx,jmax),
     &       u(ptsx,jmax), gam(ptsx), bet
      
      do j = 1, jmax
        if (my_rank .eq. 0) then
          bet    = b(1)
          u(1,j) = rhs(1,j) / bet
          gam(2) = c(1) / bet
          bet    = b(2) - a(2) * gam(2)
          u(2,j) = (rhs(2,j) - a(2) * u(1,j)) / bet
         else
          call MPI_Recv(aux, 2*inter, MPI_double_precision, my_rank-1,
     &                  57, MPI_COMM_WORLD, status, ierr)
          do i = 1, inter
            u(i,j) = aux(i)
            gam(i) = aux(i+inter)
          end do
        end if
        i_ini = inter
        if (my_rank .eq. 0) i_ini = 3
        i_fim = ptsx - 2
        if (my_rank .eq. numproc) i_fim = ptsx
        
        bet = b(i_ini-1) - a(i_ini-1) * gam(i_ini-1)
        do i = i_ini, i_fim
          gam(i) = c(i-1) / bet
          bet    = b(i) - a(i) * gam(i)
          u(i,j) = (rhs(i,j) - a(i) * u(i-1,j)) / bet
        end do
        
        if (my_rank .lt. numproc) then
          do i = 1, inter
            aux(i)       = u(ptsx-inter+i-1,j)
            aux(i+inter) = gam(ptsx-inter+i-1)
          end do
          call MPI_Send(aux, 2*inter, MPI_double_precision, my_rank+1,
     &                  57, MPI_COMM_WORLD, ierr)
        end if
      end do
      
      do j = 1, jmax
        if (my_rank .lt. numproc) then
          call MPI_Recv(aux, 2*inter, MPI_double_precision, my_rank+1, 
     &                  67, MPI_COMM_WORLD, status, ierr)
          do i = 1, inter
            u(ptsx-inter+i,j) = aux(i)
            gam(ptsx-inter+i) = aux(i+inter)
          end do
        end if
        
        i_fim = ptsx - inter
        if (my_rank .eq. numproc) i_fim = ptsx - 1
        do i = i_fim, 1, -1
          u(i,j) = u(i,j) - gam(i+1)*u(i+1,j)
        end do
        
        if (my_rank .gt. 0) then
          do i = 1, inter
            aux(i)       = u(i+1,j)
            aux(i+inter) = gam(i+1)
          end do
          call MPI_Send(aux, 2*inter, MPI_double_precision, my_rank-1,
     &                  67, MPI_COMM_WORLD, ierr)
        end if
        do i = 1, ptsx
          rhs(i,j) = u(i,j)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tridparx(a, b, c, rhs)
      
      ! solves tridiagonal matrix for the derivatives in x direction
      ! domain decomposition in x direction -> parallel subroutine
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'mpif.h'
      integer status(MPI_status_size)
      integer i, j, k, i_ini, i_fim
      real*8 a(ptsx), b(ptsx), c(ptsx)
      complex*16 aux(2*inter), rhs(ptsx,jmax,kfour), u(ptsx,jmax),
     &           gam(ptsx), bet
      
      do k = 1, kfour
        do j = 1, jmax
          if (my_rank .eq. 0) then
            bet    = b(1)
            u(1,j) = rhs(1,j,k) / bet
            gam(2) = c(1) / bet
            bet    = b(2) - a(2) * gam(2)
            u(2,j) = ( rhs(2,j,k) - a(2) * u(1,j) ) / bet
           else
            call MPI_Recv(aux, 2*inter, MPI_complex16, my_rank-1,
     &                    57, MPI_COMM_WORLD, status, ierr)
            do i = 1, inter
              u(i,j) = aux(i)
              gam(i) = aux(i+inter)
            end do
          end if
          i_ini = inter
          if (my_rank .eq. 0) i_ini = 3
          i_fim = ptsx - 2
          if (my_rank .eq. numproc) i_fim = ptsx
          
          bet = b(i_ini-1) - a(i_ini-1) * gam(i_ini-1)
          do i = i_ini, i_fim
            gam(i) = c(i-1) / bet
            bet    = b(i) - a(i) * gam(i)
            u(i,j) = ( rhs(i,j,k) - a(i) * u(i-1,j) ) / bet
          end do
          
          if (my_rank .lt. numproc) then
            do i = 1, inter
              aux(i)       = u(ptsx-inter+i-1,j)
              aux(i+inter) = gam(ptsx-inter+i-1)
            end do
            call MPI_Send(aux, 2*inter, MPI_complex16, my_rank+1,
     &                    57, MPI_COMM_WORLD, ierr)
          end if
        end do
        
        do j = 1, jmax
          if (my_rank .lt. numproc) then
            call MPI_Recv(aux, 2*inter, MPI_complex16, my_rank+1, 
     &                    67, MPI_COMM_WORLD, status, ierr)
            do i = 1, inter
              u(ptsx-inter+i,j) = aux(i)
              gam(ptsx-inter+i) = aux(i+inter)
            end do
          end if
          
          i_fim = ptsx - inter
          if (my_rank .eq. numproc) i_fim = ptsx - 1
          do i = i_fim, 1, -1
            u(i,j) = u(i,j) - gam(i+1) * u(i+1,j)
          end do
          
          if (my_rank .gt. 0) then
            do i = 1, inter
              aux(i)       = u(i+1,j)
              aux(i+inter) = gam(i+1)
            end do
            call MPI_Send(aux, 2*inter, MPI_complex16, my_rank-1,
     &                    67, MPI_COMM_WORLD, ierr)
          end if
          do i = 1, ptsx
            rhs(i,j,k) = u(i,j)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tridparxr3d(a, b, c, rhs)
      
      ! solves tridiagonal matrix for the derivatives in x direction
      ! domain decomposition in x direction -> parallel subroutine
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'mpif.h'
      integer status(MPI_status_size)
      integer i, j, k, i_ini, i_fim
      real*8 a(ptsx), b(ptsx), c(ptsx), aux(2*inter),
     &       rhs(ptsx,jmax,kphys), u(ptsx,jmax), gam(ptsx), bet
      
      do k = 1, kphys
        do j = 1, jmax
          if (my_rank .eq. 0) then
            bet    = b(1)
            u(1,j) = rhs(1,j,k) / bet
            gam(2) = c(1) / bet
            bet    = b(2) - a(2) * gam(2)
            u(2,j) = ( rhs(2,j,k) - a(2) * u(1,j) ) / bet
           else
            call MPI_Recv(aux,2*inter,MPI_double_precision,my_rank-1,
     &                    57, MPI_COMM_WORLD, status, ierr)
            do i = 1, inter
              u(i,j) = aux(i)
              gam(i) = aux(i+inter)
            end do
          end if
          i_ini = inter
          if (my_rank .eq. 0) i_ini = 3
          i_fim = ptsx - 2
          if (my_rank .eq. numproc) i_fim = ptsx
          
          bet = b(i_ini-1) - a(i_ini-1) * gam(i_ini-1)
          do i = i_ini, i_fim
            gam(i) = c(i-1) / bet
            bet    = b(i) - a(i) * gam(i)
            u(i,j) = ( rhs(i,j,k) - a(i) * u(i-1,j) ) / bet
          end do
          
          if (my_rank .lt. numproc) then
            do i = 1, inter
              aux(i)       = u(ptsx-inter+i-1,j)
              aux(i+inter) = gam(ptsx-inter+i-1)
            end do
            call MPI_Send(aux,2*inter,MPI_double_precision,my_rank+1,
     &                    57, MPI_COMM_WORLD, ierr)
          end if
        end do
        
        do j = 1, jmax
          if (my_rank .lt. numproc) then
            call MPI_Recv(aux,2*inter,MPI_double_precision,my_rank+1, 
     &                    67, MPI_COMM_WORLD, status, ierr)
            do i = 1, inter
              u(ptsx-inter+i,j) = aux(i)
              gam(ptsx-inter+i) = aux(i+inter)
            end do
          end if
          
          i_fim = ptsx - inter
          if (my_rank .eq. numproc) i_fim = ptsx - 1
          do i = i_fim, 1, -1
            u(i,j) = u(i,j) - gam(i+1) * u(i+1,j)
          end do
          
          if (my_rank .gt. 0) then
            do i = 1, inter
              aux(i)       = u(i+1,j)
              aux(i+inter) = gam(i+1)
            end do
            call MPI_Send(aux,2*inter,MPI_double_precision,my_rank-1,
     &                    67, MPI_COMM_WORLD, ierr)
          end if
          do i = 1, ptsx
            rhs(i,j,k) = u(i,j)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tridy(a, b, c, rhs)

      ! solves tridiagonal matrix for the derivatives in y direction
      implicit none
      include 'par.for'
      integer i, j, k
      real*8 a(jmax), b(jmax), c(jmax), gam(jmax), bet
      complex*16 rhs(ptsx,jmax,kfour), u(jmax)

      do k = 1, kfour
        do i = 1, ptsx
          bet  = b(1)
          u(1) = rhs(i,1,k) / bet
          do j = 2, jmax
            gam(j) = c(j-1) / bet
            bet    = b(j) - a(j) * gam(j)
            u(j)   = ( rhs(i,j,k) - a(j) * u(j-1) ) / bet
          end do
          do j = jmax - 1, 1, -1
            u(j) = u(j) - gam(j+1) * u(j+1)
          end do
          do j = 1, jmax
            rhs(i,j,k) = u(j)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tridy_i0(a, b, c, rhs)

      ! solves tridiagonal matrix for the derivatives in y direction at x = 1
      ! usado no CASE 4
      implicit none
      include 'par.for'
      integer j, k
      real*8 a(jmax), b(jmax), c(jmax), gam(jmax), bet
      complex*16 rhs(jmax,kfour), u(jmax)

      do k = 1, kfour
        bet  = b(1)
        u(1) = rhs(1,k) / bet
        do j = 2, jmax
          gam(j) = c(j-1) / bet
          bet    = b(j) - a(j) * gam(j)
          u(j)   = ( rhs(j,k) - a(j) * u(j-1) ) / bet
        end do
        do j = jmax - 1, 1, -1
          u(j) = u(j) - gam(j+1) * u(j+1)
        end do
        do j = 1, jmax
          rhs(j,k) = u(j)
        end do
      end do

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tridyr3d(a, b, c, rhs)

      ! solves tridiagonal matrix for the derivatives in y direction
      implicit none
      include 'par.for'
      integer i, j, k
      real*8 a(jmax), b(jmax), c(jmax), gam(jmax), bet,
     &       rhs(ptsx,jmax,kphys), u(jmax)

      do k = 1, kphys
        do i = 1, ptsx
          bet  = b(1)
          u(1) = rhs(i,1,k) / bet
          do j = 2, jmax
            gam(j) = c(j-1) / bet
            bet    = b(j) - a(j) * gam(j)
            u(j)   = ( rhs(i,j,k) - a(j) * u(j-1) ) / bet
          end do
          do j = jmax - 1, 1, -1
            u(j) = u(j) - gam(j+1) * u(j+1)
          end do
          do j = 1, jmax
            rhs(i,j,k) = u(j)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tridseqx(a, b, c, rhs)

      ! solves tridiagonal matrix for the derivatives in x direction
      implicit none
      include 'par.for'
      integer i, j, k
      real*8 a(ptsx), b(ptsx), c(ptsx), gam(ptsx), bet
      complex*16 rhs(ptsx,jmax,kfour), u(ptsx)

      do k = 1, kfour
        do j = 1, jmax
          bet  = b(1)
          u(1) = rhs(1,j,k) / bet
          do i = 2, ptsx
            gam(i) = c(i-1) / bet
            bet    = b(i) - a(i) * gam(i)
            u(i)   = ( rhs(i,j,k) - a(i) * u(i-1) ) / bet
          end do
          do i = ptsx - 1, 1, -1
            u(i) = u(i) - gam(i+1) * u(i+1)
          end do
          do i = 1, ptsx
            rhs(i,j,k) = u(i)
          end do
        end do
      end do

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

      subroutine coefy_fv(a, b, c)

      ! mount the LHS of the matrix for the first derivative of uy
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

      subroutine coeffy(a, b, c)

      ! mount the LHS of the matrix for the second derivative
      implicit none
      include 'par.for'
      include 'comm.coef'
      integer j
      real*8 a(jmax), b(jmax), c(jmax)

      a(1)      = 0.d0
      b(1)      = fp_sd_coef(1)
      c(1)      = fp_sd_coef(2)

      a(2)      = sp_sd_coef(1)
      b(2)      = sp_sd_coef(2)
      c(2)      = sp_sd_coef(3)

      do j = 3, jmax - 2
        a(j)    = cp_sd_coef(1)
        b(j)    = cp_sd_coef(2)
        c(j)    = cp_sd_coef(3)
      end do

      a(jmax-1) = pp_sd_coef(3)
      b(jmax-1) = pp_sd_coef(2)
      c(jmax-1) = pp_sd_coef(1)

      a(jmax)   = lp_sd_coef(2)
      b(jmax)   = lp_sd_coef(1)
      c(jmax)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine coefx(a, b, c)

      ! mount the LHS of the matrix for the first derivative
      implicit none
      include 'par.for'
      integer i
      real*8 a(ptsx), b(ptsx), c(ptsx)

      a(1)      = 0.d0
      b(1)      = 1.d0
      c(1)      = 0.d0

      a(2)      = 1.d0
      b(2)      = 6.d0
      c(2)      = 2.d0

      do i = 3, ptsx - 2
        a(i)    = 1.d0
        b(i)    = 3.d0
        c(i)    = 1.d0
      end do

      a(ptsx-1) = 2.d0
      b(ptsx-1) = 6.d0
      c(ptsx-1) = 1.d0

      a(ptsx)   = 4.d0
      b(ptsx)   = 1.d0
      c(ptsx)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine coeffx(a, b, c)

      ! mount the LHS of the matrix for the second derivative
      implicit none
      include 'par.for'
      integer i
      real*8 a(ptsx), b(ptsx), c(ptsx)

      a(1)      = 0.d0
      b(1)      = 13.d0
      c(1)      = 137.d0

      a(2)      = 1.d0
      b(2)      = 12.d0
      c(2)      = 3.d0

      do i = 3, ptsx - 2
        a(i)    = 2.d0
        b(i)    = 11.d0
        c(i)    = 2.d0
      end do

      a(ptsx-1) = 3.d0
      b(ptsx-1) = 12.d0
      c(ptsx-1) = 1.d0

      a(ptsx)   = 137.d0
      b(ptsx)   = 13.d0
      c(ptsx)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine boundary_exchange_derivs(var)

      ! exchange values of the boundaries
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'mpif.h'
      integer status(MPI_Status_size)
      integer i, j, k, tam
      complex*16 aux(meshdx*jmax*kfour), var(ptsx,jmax,kfour)

      ! variable used to calculate the number of columns needed 
      ! to go forward or backward (interm)

      tam = jmax * meshdx

      if (my_rank .lt. numproc) then
        ! Sending the near right boundary columns to node + 1
       do k = 1, kfour 
        do i = 1, meshdx
          do j = 1, jmax
            aux(j+(i-1)*jmax+(k-1)*tam) = var(i+ptsx-24-1,j,k)
          end do
        end do
       end do 
        call MPI_Send(aux, meshdx*jmax*kfour, MPI_COMPLEX16, my_rank+1,
     &                10, MPI_COMM_WORLD, ierr)
        ! Receiving the new right boundary columns from node + 1
        call MPI_Recv(aux, meshdx*jmax*kfour, MPI_COMPLEX16, my_rank+1,
     &                20, MPI_COMM_WORLD, status, ierr)
       do k = 1, kfour 
        do i = 1, meshdx
          do j = 1, jmax
            var(i+ptsx-meshdx,j,k) = aux(j+(i-1)*jmax+(k-1)*tam)
          end do
        end do
       end do 
       
       end if

      if (my_rank .gt. 0) then
        ! Receiving the new left boundary columns from node - 1
        call MPI_Recv(aux, meshdx*jmax*kfour, MPI_COMPLEX16, my_rank-1,
     &                10, MPI_COMM_WORLD, status, ierr)
       do k = 1, kfour 
        do i = 1, meshdx
          do j = 1, jmax
            var(i,j,k) = aux(j+(i-1)*jmax+(k-1)*tam)
          end do
        end do
       end do 
        ! Sending the new right boundary columns to node - 1
       do k = 1, kfour 
        do i = 1, meshdx
          do j = 1, jmax
            aux(j+(i-1)*jmax+(k-1)*tam) = var(i+24-meshdx+1,j,k)
          end do
        end do
       end do 
        call MPI_Send(aux, meshdx*jmax*kfour, MPI_COMPLEX16, my_rank-1,
     &                20, MPI_COMM_WORLD, ierr)
      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                      end of derivative calculations                          c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
