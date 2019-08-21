ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c                  calculates de Reynolds statistics                          c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine statistics_reynolds(l,wzp) 

      implicit none
      include 'par.for'
      character*30 nome
      integer i, j, k, l, my_rank, shift
      real*8  x, y, var(6), utau, soma1,
     &        ubt(ptsx,jmax,kphys), vbt(ptsx,jmax,kphys), 
     &        wbt(ptsx,jmax,kphys),
     &        uubt(ptsx,jmax,kphys), vvbt(ptsx,jmax,kphys),
     &        wwbt(ptsx,jmax,kphys),
     &        uvbt(ptsx,jmax,kphys), uwbt(ptsx,jmax,kphys), 
     &        vwbt(ptsx,jmax,kphys),
     &        wzbt(ptsx,jmax,kphys),
     &        ub(imax,jmax,kphys), vb(imax,jmax,kphys), 
     &        wb(imax,jmax,kphys),
     &        uub(imax,jmax,kphys), vvb(imax,jmax,kphys),
     &        wwb(imax,jmax,kphys),
     &        uvb(imax,jmax,kphys), uwb(imax,jmax,kphys), 
     &        vwb(imax,jmax,kphys),
     &        wzbb(imax,jmax,kphys),
     &        uulb(imax,jmax,kphys), vvlb(imax,jmax,kphys),
     &        wwlb(imax,jmax,kphys),
     &        uvlb(imax,jmax,kphys), uwlb(imax,jmax,kphys),
     &        vwlb(imax,jmax,kphys), 
     &        dudy(imax, jmax, kphys),
     &        wzp(imax, jmax, kphys),
     &        x_ad, x_dim, delta
      common/media_rey/ ub, vb, wb
      common/calculo_utau/ utau, delta


cc    reading media u, v, w, uu, vv, ww, uv, uw, vw, wz 
      do my_rank = 0, np - 1
        write(nome,'(a,i0.2,a)')'media_esp_temp_',my_rank,'.bin'
        open(2, file = nome, form = 'unformatted')
        read(2) ubt, vbt, wbt, uubt, vvbt, wwbt, uvbt, uwbt, vwbt, wzbt
        close (unit = 2)

        shift = my_rank * (ptsx - inter - 1)
        do k = 1, kphys
          do j = 1, jmax
            do i = 1, ptsx
              ub(i+shift,j,k)  = ubt(i,j,k)
              vb(i+shift,j,k)  = vbt(i,j,k)
              wb(i+shift,j,k)  = wbt(i,j,k)
              uub(i+shift,j,k) = uubt(i,j,k)
              vvb(i+shift,j,k) = vvbt(i,j,k)
              wwb(i+shift,j,k) = wwbt(i,j,k)
              uvb(i+shift,j,k) = uvbt(i,j,k)
              uwb(i+shift,j,k) = uwbt(i,j,k)
              vwb(i+shift,j,k) = vwbt(i,j,k)
              wzbb(i+shift,j,k) = wzbt(i,j,k)
            end do
          end do
        end do
      end do


      i = l !ponto x em que sao calculados utau, delta e estatisticas


cc    calculate utau x fixo e k = 1, kphys (ADIMENSIONAL)
      call derivs_k_phys
      call dery_phys(dudy, ub)
 
      !soma1  = dsqrt((1.0d0/Re)*dabs(wzp(i,1,1)))
      !soma1  = dsqrt((1.0d0/Re)*dabs(dudy(i,1,1)))
      soma1 = dsqrt((1.0d0/Re)*dabs(wzbb(i,1,1)))
      do k = 1, kphys
        !soma1 = soma1 + dsqrt((1.0d0/Re)*(wzp(i,1,k)))
        !soma1 = soma1 + dsqrt((1.0d0/Re)*(dudy(i,1,k)))
        soma1 = soma1 + dsqrt((1.0d0/Re)*(wzbb(i,1,k)))
      enddo
      utau = soma1/dble(kphys+1) !media utau na parede somente entre i0 e i3
      write(*,*) 'utau=', utau


cc    turbulent boundary layer thickness at x = i0 (DIMENSIONAL)
      x_ad  = 1.d0                                       ! non dimensional x at inflow
      x_dim = (x_ad + dble(i-1)*dx) * L_1                ! dimesnional x no ponto i0
      delta = 0.37d0 * x_dim / (U_1*x_dim/N_1)**0.2      ! espessura CL turb em x = i0


cc    calculate the Reynolds stress (ADIMENSIONAL)

      ! resolved stresses (variance of f_bar)
      do k = 1, kphys
        do j = 1, jmax
          do i = 1, imax
            uulb(i,j,k) = uub(i,j,k)-(ub(i,j,k)*ub(i,j,k))
            vvlb(i,j,k) = vvb(i,j,k)-(vb(i,j,k)*vb(i,j,k))
            wwlb(i,j,k) = wwb(i,j,k)-(wb(i,j,k)*wb(i,j,k))
            
            uvlb(i,j,k) = uvb(i,j,k)-(ub(i,j,k)*vb(i,j,k))
            uwlb(i,j,k) = uwb(i,j,k)-(ub(i,j,k)*wb(i,j,k))
            vwlb(i,j,k) = vwb(i,j,k)-(vb(i,j,k)*wb(i,j,k))
          enddo
        enddo
      enddo  

      write(nome,'(a,i0.2,a)')'statistics.dat'
      open (44,file = nome, status = 'unknown')
      write(44,*) 'VARIABLES= "y","y/delta","uu","vv","ww",
     &                         "uv","uw","vw"'
      write(44,*) 'ZONE J=',jmax, ',F=POINT'
      i = l
      do j = 1, jmax
        if (stf .eq. 1.d0) then
          y = dble(j-1) * dy
          else
           y = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)
        end if
c       do i = 1, imax
          x = x0 + dble(i-1) * dx
          do k = 1, 6
            var(k) = 0.d0
          enddo
          var(1) = uulb(i,j,1)
          var(2) = vvlb(i,j,1)
          var(3) = wwlb(i,j,1)
          var(4) = uvlb(i,j,1)
          var(5) = uwlb(i,j,1)
          var(6) = vwlb(i,j,1)
          do k = 1, kphys
            var(1) = var(1) + uulb(i,j,k)
            var(2) = var(2) + vvlb(i,j,k)
            var(3) = var(3) + wwlb(i,j,k)
            var(4) = var(4) + uvlb(i,j,k)
            var(5) = var(5) + uwlb(i,j,k)
            var(6) = var(6) + vwlb(i,j,k)
          enddo
          do k = 1, 6
            var(k) = var(k) / dble(kphys+1)
          enddo
          write(44,6) y,y/(delta/L_1),dsqrt(var(1))/utau,
     &   dsqrt(var(2))/utau, dsqrt(var(3))/utau,
     &   var(4)/(utau**2), var(5)/(utau**2), var(6)/(utau**2)
c       enddo
      enddo
      close (unit = 44)

    6 format(1x, 2d14.6, 6d17.9)
     
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c                  calculates de u+ versus y+                                 c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine mean_vel(l) 

      implicit none
      include 'par.for'
      character*30 nome
      integer i, j, k, l 
      real*8  x, y, utau, yplus, soma, soma1, soma2, soma3, 
     &        aux(imax,kphys),
     &        ub(imax,jmax,kphys), vb(imax,jmax,kphys), 
     &        wb(imax,jmax,kphys),
     &        x_ad, x_dim, delta
      common/media_rey/ ub, vb, wb
      common/calculo_utau/ utau, delta


      i = l !ponto x em que as médias são calculadas


      write(nome,'(a,i0.2,a)')'umean.dat'
      open (22,file = nome, status = 'unknown')
      write(22,*) 'VARIABLES= "y+","y/delta","umean+","vmean+","wmean+"'
      write(22,*) 'ZONE J=',jmax-1, ',F=POINT'
      do j = 2, jmax
        if (stf .eq. 1.d0) then
          y = dble(j-1) * dy
          else
           y = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)
        end if
        yplus = y * utau * Re
        soma1 = ub(i,j,1)
        soma2 = vb(i,j,1)
        soma3 = wb(i,j,1)
        do k = 1, kphys
          soma1 = soma1 + ub(i,j,k)
          soma2 = soma2 + vb(i,j,k)
          soma3 = soma3 + wb(i,j,k)
        end do
        soma1 = (soma1 / dble(kphys+1))/utau
        soma2 = (soma2 / dble(kphys+1))/utau
        soma3 = (soma3 / dble(kphys+1))/utau
        write(22,*) yplus, y/(delta/L_1), soma1, soma2, soma3

      end do
      close(unit = 22)
 
      return
      end 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c    calculates the integral properties: delta1 (displacement thickness)      c 
c                                        delat2 (momnetum thickness)          c
c                                        H12 (sahpe factor)                   c
c                                               (dimensionais)                c       
c                                                                             c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine integ_properties 

      implicit none
      include 'par.for'
      character*30 nome
      integer i, j, k, t, my_rank, shift, ii
      real*8 x, y(jmax), a, b, c, det, fc1, fc2, umax,
     &       delta1(i3-i0+1), delta2(i3-i0+1), H12(i3-i0+1),
     &       var1(i3-i0+1,jmax), var2(i3-i0+1,jmax),
     &       ub(imax,jmax,kphys), vb(imax,jmax,kphys), 
     &       wb(imax,jmax,kphys), 
     &       uxp(imax,jmax,kphys), uyp(imax,jmax,kphys), 
     &       uzp(imax,jmax,kphys)
      common/media_rey/ ub, vb, wb


      do j = 1, jmax
        if (stf .ne. 1.d0) then 
         y(j) = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)
        else
         y(j) = dy * dble(j-1)
        end if        
      end do

      k = 1
      do i = i0, i3
        ii = (i - i0) + 1
        umax = ub(i,jmax,k) !! ver se para escoamento turbulento é isso mesmo
        do j = 1, jmax
          var1(ii,j) = 1.d0 - ub(i,j,k) / umax
          var2(ii,j) = ub(i,j,k) / umax * (1.d0 - ub(i,j,k) / umax)
        end do
      end do

      do i = i0, i3
        ii = (i - i0) + 1
        fc1 = 0.d0
        fc2 = 0.d0
        do j = 2, jmax - 1, 2
          det = y(j-1)**2 * ( y(j) - y(j+1) )
     &        + y(j-1) * ( y(j+1)**2 - y(j)**2 )
     &        + y(j) * y(j+1) * ( y(j) - y(j+1) )
          a   = (1.d0 / det) * ( var1(ii,j-1) * (y(j)   - y(j+1))
     &                         + var1(ii,j)   * (y(j+1) - y(j-1))
     &                         + var1(ii,j+1) * (y(j-1) - y(j)))
          b   = (1.d0 / det) * ( var1(ii,j-1) * (y(j+1)**2 - y(j)**2)
     &                         + var1(ii,j)   * (y(j-1)**2 - y(j+1)**2)
     &                         - var1(ii,j+1) * (y(j-1)**2 - y(j)**2))
          c   = (1.d0 / det) * ( var1(ii,j-1) *
     &                           y(j)  * y(j+1) * (y(j)   - y(j+1))
     &                         + var1(ii,j)   *
     &                           y(j-1)* y(j+1) * (y(j+1) - y(j-1))
     &                         + var1(ii,j+1) *
     &                           y(j)  * y(j-1) * (y(j-1) - y(j)))
          fc1 = fc1 + (a / 3.d0) * (y(j+1)**3 - y(j-1)**3)
     &              + (b / 2.d0) * (y(j+1)**2 - y(j-1)**2)
     &              +  c         * (y(j+1)    - y(j-1))

          a   = (1.d0 / det) * ( var2(ii,j-1) * (y(j)   - y(j+1))
     &                         + var2(ii,j)   * (y(j+1) - y(j-1))
     &                         + var2(ii,j+1) * (y(j-1) - y(j)))
          b   = (1.d0 / det) * ( var2(ii,j-1) * (y(j+1)**2 - y(j)**2)
     &                         + var2(ii,j)   * (y(j-1)**2 - y(j+1)**2)
     &                         - var2(ii,j+1) * (y(j-1)**2 - y(j)**2))
          c   = (1.d0 / det) * ( var2(ii,j-1)
     &                         * y(j)  * y(j+1) * (y(j)   - y(j+1))
     &                         + var2(ii,j)
     &                         * y(j-1)* y(j+1) * (y(j+1) - y(j-1))
     &                         + var2(ii,j+1)
     &                         * y(j)  * y(j-1) * (y(j-1) - y(j)))
          fc2 = fc2 + (a / 3.d0) * (y(j+1)**3 - y(j-1)**3)
     &              + (b / 2.d0) * (y(j+1)**2 - y(j-1)**2)
     &              +  c         * (y(j+1)    - y(j-1))
        end do

        delta1(ii) = (fc1 / dsqrt(fac_y))*L_1    ! dimensionalizando
        delta2(ii) = (fc2 / dsqrt(fac_y))*L_1    ! dimensionalizando
        H12(ii) = (delta1(ii) / delta2(ii))*L_1  ! dimensionalizando
      enddo 
 
      open (66, file = 'integ_properties.dat',status = 'unknown')
      write(66,*) 'VARIABLES="x","delta1","delta2","H12"'
      write(66,*) 'ZONE I=',i3-i0+1,',  F=POINT'
      do i = i0, i3
        ii = (i - i0) + 1
        x  = x0 + dble(i-1) * dx
        write(66,5) x, delta1(ii), delta2(ii), H12(ii)
      end do
      close (unit = 1)
 
    5 format(1x, 4d25.17)
      return

      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine derivs_k_phys

      ! mounts the tri-diagonal LHS for derivative calculations
      implicit none
      include 'par.for'
      real*8  a1y(jmax),  b1y(jmax),  c1y(jmax)
      common/der1y/  a1y,  b1y,  c1y

      call coefy_phys(a1y, b1y, c1y)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine coefy_phys(a, b, c)

      ! mount the LHS of the matrix for the first derivative
      implicit none
      include 'par.for'
      include 'comm.coef'
      integer j
      real*8 a(jmax), b(jmax), c(jmax)

      ! reads the derivative and Poisson coefficients
      open(1, file = 'pre_processing/coefs.bin',form = 'unformatted')
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

      subroutine dery_phys(ddy, fc)

      ! first derivative calculation in y direction
      implicit none
      include 'par.for'
      real*8 a1y(jmax), b1y(jmax), c1y(jmax)
      real*8 fc(imax,jmax,kphys), ddy(imax,jmax,kphys)
      common/der1y/ a1y, b1y, c1y

      call rhsy_phys(fc, ddy)
      call tridy_phys(a1y, b1y, c1y, ddy)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsy_phys(fc,rhs)

      ! RHS for the first derivative calculation in y direction
      implicit none
      include 'par.for'
      include 'comm.coef'
      integer i, j, k, lvl
      real*8 rhs(imax,jmax,kphys), fc(imax,jmax,kphys), v_qdy(jmax),
     &       v_dy0(msh), v_stf(msh)

      ! dy0 at each multigrid level
      do lvl = 1 , msh
        v_stf(lvl) = stf ** ( 2**(lvl-1) )
        if (stf .ne. 1.d0) then
          v_dy0(lvl) = dy * ( ( stf**(2**(lvl-1)) - 1.d0) /
     &                 (stf - 1.d0) )
         else
          v_dy0(lvl) = dy * 2.d0**(lvl-1)
        end if
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

      do k = 1, kphys
        do i = 1, imax
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
         
         rhs(i,jmax-1,k) = v_qdy(jmax-5) *
     &                     ( pp_fd_coef(4) * fc(i,jmax,k)
     &                     + pp_fd_coef(5) * fc(i,jmax-1,k)
     &                     + pp_fd_coef(6) * fc(i,jmax-2,k)
     &                     + pp_fd_coef(7) * fc(i,jmax-3,k)
     &                     + pp_fd_coef(8) * fc(i,jmax-4,k)
     &                     + pp_fd_coef(9) * fc(i,jmax-5,k) )
         
         rhs(i,jmax,k) = v_qdy(jmax-4) *
     &                   ( lp_fd_coef(3) * fc(i,jmax,k)
     &                   + lp_fd_coef(4) * fc(i,jmax-1,k)
     &                   + lp_fd_coef(5) * fc(i,jmax-2,k)
     &                   + lp_fd_coef(6) * fc(i,jmax-3,k)
     &                   + lp_fd_coef(7) * fc(i,jmax-4,k) )
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tridy_phys(a, b, c, rhs)

      ! solves tridiagonal matrix for the derivatives in y direction
      implicit none
      include 'par.for'
      integer i, j, k
      real*8 a(jmax), b(jmax), c(jmax), gam(jmax), bet
      real*8 rhs(imax,jmax,kphys), u(jmax)

      do k = 1, kphys
        do i = 1, imax
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
