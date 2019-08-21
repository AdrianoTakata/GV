cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                         nonlinear terms calculation                          c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine lterms_gv(a, b, c, d)

      implicit none
      include 'par.for'
      include 'comm.var'
      integer i, j, k
      real*8 uxb(ptsx,jmax), uyb(ptsx,jmax), wzb(ptsx,jmax)
      complex*16 a(ptsx,jmax,kfour),   b(ptsx,jmax,kfour),
     &           c(ptsx,jmax,kfour),   d(ptsx,jmax,kfour)
      common/blas/ uxb, uyb, wzb

      do k = 1, kfour
        do j = 1, jmax
          do i = 1, ptsx
            a(i,j,k) =   uyb(i,j) * wx(i,j,k) - wy(i,j,k) * uxb(i,j)
            b(i,j,k) =   uxb(i,j) * wz(i,j,k) + ux(i,j,k) * wzb(i,j)
            c(i,j,k) = - uyb(i,j) * wz(i,j,k) - uy(i,j,k) * wzb(i,j)
            d(i,j,k) =   2.d0 * uxb(i,j) * ux(i,j,k)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine nlterms_gv(a, b, c, d)

      ! calculate the non linear terms of transport
      ! equations, in terms are ux, uy, uz, wx, wy and wz and  
      ! out terms are a, b, c and d
      !		a = uy * wx - ux * wy
      !		b = ux * wz - uz * wx
      !		c = uz * wy - uy * wz
      !		d = (Uxb + ux)^2
      implicit none
      include 'par.for'
      include 'comm.var'
      integer i, j, k
      real*8 uxb(ptsx,jmax), uyb(ptsx,jmax), wzb(ptsx,jmax),
     &          uxp(ptsx,jmax,kphys), wxp(ptsx,jmax,kphys),
     &          uyp(ptsx,jmax,kphys), wyp(ptsx,jmax,kphys),
     &          uzp(ptsx,jmax,kphys), wzp(ptsx,jmax,kphys),
     &           ap(ptsx,jmax,kphys),  bp(ptsx,jmax,kphys),
     &           cp(ptsx,jmax,kphys),  dp(ptsx,jmax,kphys)
      complex*16  a(ptsx,jmax,kfour),   b(ptsx,jmax,kfour),
     &            c(ptsx,jmax,kfour),   d(ptsx,jmax,kfour)
      common/blas/ uxb, uyb, wzb

      ! fft transforms from Fourier to physical space
      call f_to_p(uxp, ux, ptsx)
      call f_to_p(uyp, uy, ptsx)
      call f_to_p(uzp, uz, ptsx)
      call f_to_p(wxp, wx, ptsx)
      call f_to_p(wyp, wy, ptsx)
      call f_to_p(wzp, wz, ptsx)

      ! nonlinear terms calculation in physical space
      do k = 1, kphys
        do j = 1, jmax
          do i = 1, ptsx
           ap(i,j,k) = uyp(i,j,k) * wxp(i,j,k) - uxp(i,j,k) * wyp(i,j,k)
           bp(i,j,k) = uxp(i,j,k) * wzp(i,j,k) - uzp(i,j,k) * wxp(i,j,k)
           cp(i,j,k) = uzp(i,j,k) * wyp(i,j,k) - uyp(i,j,k) * wzp(i,j,k)
           dp(i,j,k) = uxp(i,j,k) * uxp(i,j,k)
          end do
        end do
      end do

      ! fft back from physical to Fourier space
      call p_to_f(ap, a, ptsx)
      call p_to_f(bp, b, ptsx)
      call p_to_f(cp, c, ptsx)
      call p_to_f(dp, d, ptsx)

      ! adding the base flow to the nonlinear product
      do k = 1, kfour
        do j = 1, jmax
          do i = 1, ptsx
            a(i,j,k) =   uyb(i,j) * wx(i,j,k) - uxb(i,j) * wy(i,j,k)
     &               +   a(i,j,k)
            b(i,j,k) =   uxb(i,j) * wz(i,j,k) + ux(i,j,k) * wzb(i,j)
     &               +   b(i,j,k)
            c(i,j,k) = - uyb(i,j) * wz(i,j,k) - uy(i,j,k) * wzb(i,j)
     &               +   c(i,j,k)
            d(i,j,k) =   2.d0 * uxb(i,j) * ux(i,j,k) + d(i,j,k)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine lterms_th(a, b, c, d, uth, vth, wth)

      implicit none
      include 'par.for'
      include 'comm.var'
      integer i, j, k
      real*8 uxb(ptsx,jmax), uyb(ptsx,jmax), wzb(ptsx,jmax),
     &       thb(ptsx,jmax)
      complex*16   a(ptsx,jmax,kfour),   b(ptsx,jmax,kfour),
     &             c(ptsx,jmax,kfour),   d(ptsx,jmax,kfour),
     &           uth(ptsx,jmax,kfour), vth(ptsx,jmax,kfour),
     &           wth(ptsx,jmax,kfour)
      common/blas/ uxb, uyb, wzb
      common/blast/ thb

      do k = 1, kfour
        do j = 1, jmax
          do i = 1, ptsx
            a(i,j,k)   =   uyb(i,j) * wx(i,j,k) - wy(i,j,k) * uxb(i,j)
            b(i,j,k)   =   uxb(i,j) * wz(i,j,k) + ux(i,j,k) * wzb(i,j)
            c(i,j,k)   = - uyb(i,j) * wz(i,j,k) - uy(i,j,k) * wzb(i,j)
            d(i,j,k)   =   2.d0 * uxb(i,j) * ux(i,j,k)
            uth(i,j,k) =   uxb(i,j) * th(i,j,k) + ux(i,j,k) * thb(i,j)
            vth(i,j,k) =   uyb(i,j) * th(i,j,k) + uy(i,j,k) * thb(i,j)
            wth(i,j,k) =   uz(i,j,k) * thb(i,j)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine nlterms_th(a, b, c, d, uth, vth, wth)

      ! calculate the non linear terms of transport
      ! equations, in terms are ux, uy, uz, wx, wy and wz and  
      ! out terms are a, b, c and d
      !		a = uy * wx - ux * wy
      !		b = ux * wz - uz * wx
      !		c = uz * wy - uy * wz
      !		d = (Uxb + ux)^2
      implicit none
      include 'par.for'
      include 'comm.var'
      integer i, j, k
      real*8    uxb(ptsx,jmax), uyb(ptsx,jmax), wzb(ptsx,jmax),
     &          thb(ptsx,jmax),
     &          uxp(ptsx,jmax,kphys), wxp(ptsx,jmax,kphys),
     &          uyp(ptsx,jmax,kphys), wyp(ptsx,jmax,kphys),
     &          uzp(ptsx,jmax,kphys), wzp(ptsx,jmax,kphys),
     &         uthp(ptsx,jmax,kphys),vthp(ptsx,jmax,kphys),
     &         wthp(ptsx,jmax,kphys), thp(ptsx,jmax,kphys),
     &           ap(ptsx,jmax,kphys),  bp(ptsx,jmax,kphys),
     &           cp(ptsx,jmax,kphys),  dp(ptsx,jmax,kphys)
      complex*16 uth(ptsx,jmax,kfour), vth(ptsx,jmax,kfour),
     &           wth(ptsx,jmax,kfour),
     &             a(ptsx,jmax,kfour),   b(ptsx,jmax,kfour),
     &             c(ptsx,jmax,kfour),   d(ptsx,jmax,kfour)
      common/blas/ uxb, uyb, wzb
      common/blast/ thb

      ! fft transforms from Fourier to physical space
      call f_to_p(uxp, ux, ptsx)
      call f_to_p(uyp, uy, ptsx)
      call f_to_p(uzp, uz, ptsx)
      call f_to_p(wxp, wx, ptsx)
      call f_to_p(wyp, wy, ptsx)
      call f_to_p(wzp, wz, ptsx)
      call f_to_p(thp, th, ptsx)

      ! nonlinear terms calculation in physical space
      do k = 1, kphys
        do j = 1, jmax
          do i = 1, ptsx
           ap(i,j,k)   = uyp(i,j,k) * wxp(i,j,k)
     &                 - uxp(i,j,k) * wyp(i,j,k)
           bp(i,j,k)   = uxp(i,j,k) * wzp(i,j,k)
     &                 - uzp(i,j,k) * wxp(i,j,k)
           cp(i,j,k)   = uzp(i,j,k) * wyp(i,j,k)
     &                 - uyp(i,j,k) * wzp(i,j,k)
           dp(i,j,k)   = uxp(i,j,k) * uxp(i,j,k)
           uthp(i,j,k) = uxp(i,j,k) * thp(i,j,k)
           vthp(i,j,k) = uyp(i,j,k) * thp(i,j,k)
           wthp(i,j,k) = uzp(i,j,k) * thp(i,j,k)
          end do
        end do
      end do

      ! fft back from physical to Fourier space
      call p_to_f(ap, a, ptsx)
      call p_to_f(bp, b, ptsx)
      call p_to_f(cp, c, ptsx)
      call p_to_f(dp, d, ptsx)
      call p_to_f(uthp, uth, ptsx)
      call p_to_f(vthp, vth, ptsx)
      call p_to_f(wthp, wth, ptsx)

      ! adding the base flow to the nonlinear product
      do k = 1, kfour
        do j = 1, jmax
          do i = 1, ptsx
            a(i,j,k)   =   uyb(i,j) * wx(i,j,k) - uxb(i,j) * wy(i,j,k)
     &                 +   a(i,j,k)
            b(i,j,k)   =   uxb(i,j) * wz(i,j,k) + ux(i,j,k) * wzb(i,j)
     &                 +   b(i,j,k)
            c(i,j,k)   = - uyb(i,j) * wz(i,j,k) - uy(i,j,k) * wzb(i,j)
     &                 +   c(i,j,k)
            d(i,j,k)   =   2.d0 * uxb(i,j) * ux(i,j,k) + d(i,j,k)
            uth(i,j,k) =   uxb(i,j) * th(i,j,k) + ux(i,j,k) * thb(i,j)
     &                 +   uth(i,j,k)
            vth(i,j,k) =   uyb(i,j) * th(i,j,k) + uy(i,j,k) * thb(i,j)
     &                 +   vth(i,j,k)
            wth(i,j,k) =   uz(i,j,k) * thb(i,j) + wth(i,j,k)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                     end of nonlinear terms calculation                       c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
