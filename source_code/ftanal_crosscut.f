ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c                       Fourier analysis in time prg                          c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program fourier analisys_crosscut 08122009

      implicit none
      include 'par.for'
      character*20 nome, nome2, nome3, nome4
      integer i, j, k, f, N, M, jj, kk, var, cross
      parameter ( N = 64, M = 6 )! 2^M = N
      real*8 a(N), b(N), c(N), d(N), 
     &       mode00(imax/8+1,jmax),  Amp_max_00(imax/8+1),
     &       alpha_x(imax/8+1,N), alpha_xd2(imax/8+1,N), 
     &       Amp_max(imax/8+1,N), log_Amp_max(imax/8+1,N),
     &       Amp_maxd2(imax/8+1,N), log_Amp_maxd2(imax/8+1,N),
     &       Amp(imax/8+1,jmax,kphys,N), phase(imax/8+1,jmax,kphys,N),
     &       Ampd2(imax/8+1,jmax,kphys,N),
     &       phased2(imax/8+1,jmax,kphys,N),
     &       ux_p(imax/8+1,jmax,kphys,N), d2ux_p(imax/8+1,jmax,kphys,N),
     &       x(3), y(3), coef(3), yy, ymax, xad, z

      call initval(ux_p, d2ux_p, mode00, N)

      ! max amplitude of mode (0,0)
      do i = 1, imax/8 + 1
        yy = 1.d-30
        jj = 3
        do j = 1, jmax
          if (mode00(i,j) .gt. yy) then 
            jj = j
            yy = dabs(mode00(i,j))
          end if
        end do
        x(1) = dble(jj-2) * dy
        x(2) = dble(jj-1) * dy
        x(3) = dble(jj)   * dy
        y(1) = dabs(mode00(i,jj-1))
        y(2) = dabs(mode00(i,jj))
        y(3) = dabs(mode00(i,jj+1))
        call polcoe(x, y, coef)
        ymax             = - 0.5d0 * coef(2) / coef(3)
        Amp_max_00(i)    =   coef(3) * ymax * ymax
     &                   +   coef(2) * ymax + coef(1)
        log_Amp_max(i,1) =   dlog(Amp_max_00(i))
      end do
      do i = 5, imax/8 + 1 - 3
        alpha_x(i,1) = ( -           log_Amp_max(i-3,1)
     &                   +           log_Amp_max(i+3,1)
     &                   +  9.d0 * ( log_Amp_max(i-2,1)
     &                             - log_Amp_max(i+2,1) )
     &                   - 45.d0 * ( log_Amp_max(i-1,1)
     &                             - log_Amp_max(i+1,1) ) ) /
     &                   ( 60.d0 * dx * 4.d0 )
      end do
      open (unit = 1, file = 'mode0_0.dat', status = 'unknown')
      write(1,*) 'VARIABLES="x","U_max","alpha_x",
     &"U_maxd2","alpha_xd2"'
      write(1,*) 'ZONE T="mode(0,0)", I=',imax/8
      do i = 2, imax/8 + 1
        xad = x0 + dble(i-1) * dx * 4.d0
        write(1,8) xad, Amp_max_00(i), alpha_x(i,1), 
     &             Amp_max_00(i), alpha_x(i,1)
      end do
      close (unit = 1)
      ! end of max amplitude of mode (0,0)

      ! max amplitude of all modes in planes z,y (omega_t=0 to N/2*omega)
      do i = 1, imax/8 + 1
        do k = 1, kphys
          do j = 1, jmax
            do f = 1, N
              a(f) = ux_p(i,j,k,f)
              b(f) = 0.d0
              c(f) = d2ux_p(i,j,k,f)
              d(f) = 0.d0
            end do
            call fft1d ( a, b, .FALSE., N, M)
            call fft1d ( c, d, .FALSE., N, M)
            Amp(i,j,k,1)   = dsqrt(a(1) * a(1) + b(1) * b(1))
            phase(i,j,k,1) = datan2(b(1), a(1))
            do f = 2, N / 2
              Amp(i,j,k,f)     = dsqrt(a(f) * a(f) + b(f) * b(f))
              phase(i,j,k,f)   = datan2(b(f), a(f))
              Ampd2(i,j,k,f)   = dsqrt(c(f) * c(f) + d(f) * d(f)) /
     &                           ( (dble(f-1) * omega)**2 )
              phased2(i,j,k,f) = datan2(d(f), c(f))
            end do
          end do
        end do
      end do
      ! max amplitude of each frequency for the function
      do f = 1, N / 4 + 1
        do i = 2, imax/8 + 1
          yy = 1.d-30
          kk = 1
          jj = 3
          do k = 1, kphys
            do j = 3, jmax - 1
              if (Amp(i,j,k,f) .gt. yy) then
                jj = j
                kk = k
                yy = dabs(Amp(i,j,k,f))
              end if
            end do
          end do
          x(1) = dble(jj-2) * dy
          x(2) = dble(jj-1) * dy
          x(3) = dble(jj)   * dy
          y(1) = dabs(Amp(i,jj-1,kk,f))
          y(2) = dabs(Amp(i,jj,kk,f))
          y(3) = dabs(Amp(i,jj+1,kk,f))
          call polcoe(x, y, coef)
          ymax             = - 0.5d0 * coef(2) / coef(3)
          Amp_max(i,f)     =   coef(3) * ymax * ymax
     &                     +   coef(2) * ymax + coef(1)
          log_Amp_max(i,f) =   dlog(Amp_max(i,f))
        end do
      end do
      ! end of max amplitude of each frequency for the function
      ! max amplitude of each frequency for the second derivative of the function
      do f = 2, N / 4 + 1
        do i = 2, imax/8 + 1
          yy = 1.d-30
          kk = 1
          jj = 3
          do k = 1, kphys
            do j = 3, jmax - 1
              if (Ampd2(i,j,k,f) .gt. yy) then 
                jj = j
                kk = k
                yy = dabs(Ampd2(i,j,k,f))
              end if
            end do
          end do
          x(1) = dble(jj-2) * dy
          x(2) = dble(jj-1) * dy
          x(3) = dble(jj)   * dy
          y(1) = dabs(Ampd2(i,jj-1,kk,f))
          y(2) = dabs(Ampd2(i,jj,kk,f))
          y(3) = dabs(Ampd2(i,jj+1,kk,f))
          call polcoe(x, y, coef)
          ymax               = - 0.5d0 * coef(2) / coef(3)
          Amp_maxd2(i,f)     =   coef(3) * ymax * ymax
     &                       +   coef(2) * ymax + coef(1)
          log_Amp_maxd2(i,f) =   dlog(Amp_maxd2(i,f))
        end do
      end do
      ! end of max amplitude of each frequency for the second derivative of the function

      ! alpha calculation (growth rate)
      call derfax(alpha_x, log_Amp_max, N)
      call derfaxd2(alpha_xd2, log_Amp_maxd2, N)
      ! end of alpha calculation (growth rate)

      ! write the results
      do f = 1, N / 4 + 1

        write(nome,'(a,i0.2,a)')'mode',f-1,'.dat'
        write(nome2,'(a,i0.2,a)')'mode(',f-1,')'
        open (1, file = nome, status = 'unknown')
        write(1,*) 'VARIABLES="x","U_max","alpha_x",
     &              "U_maxd2","alpha_xd2"'
        write(1,*) 'ZONE T="',nome2,'", I=',imax / 8
        do i = 2, imax / 8 + 1
          xad = x0 + dble(i-1) * dx * 4.d0
          write(1,8) xad, Amp_max(i,f), alpha_x(i,f), 
     &               Amp_maxd2(i,f), alpha_xd2(i,f)
        end do
        close (unit = 1)

        do cross = 1, 10
          if (cross .eq. 1)  i =  40
          if (cross .eq. 2)  i =  54
          if (cross .eq. 3)  i = 108
          if (cross .eq. 4)  i = 144
          if (cross .eq. 5)  i = 180
          if (cross .eq. 6)  i = 196
          if (cross .eq. 7)  i = 217
          if (cross .eq. 8)  i = 234
          if (cross .eq. 9)  i = 252
          if (cross .eq. 10) i = 270
          var = 1000.d0 * ( ( 4.d0 * (i - 1) * dx ) + x0 ) * L_1
          write(nome3,'(a,i0.2,a,i0.4)')'cc_mode_',f-1,'_',var
          write(nome4,'(a,i0.2,a,i0.4)')'cc_mode_d2',f-1,'_',var
          open (2, file = nome3, status = 'unknown')
          write(2,*) 'VARIABLES="z","y","ux","phase"'
          write(2,*) 'ZONE I=',kphys+1,', J=',jmax,', F=POINT'
          open (3, file = nome4, status = 'unknown')
          write(3,*) 'VARIABLES="z","y","ux","phase"'
          write(3,*) 'ZONE I=',kphys+1,', J=',jmax,', F=POINT'

          do j = 1, jmax
            yy =   dble(j-1) * dy
            z  = - 1.d0 * 2.d0 * pi / ( dble(kphys) * beta )
            write(2,5) z, yy, Amp(i,j,kphys,f), phase(i,j,kphys,f)
            write(3,5) z, yy, Ampd2(i,j,kphys,f), phased2(i,j,kphys,f)
            do k = 1, kphys
              z = dble(k-1) * 2.d0 * pi / (dble(kphys) * beta)
              write(2,5) z, yy, Amp(i,j,k,f), phase(i,j,k,f)
              write(3,5) z, yy, Ampd2(i,j,k,f), phased2(i,j,k,f)
            end do
          end do
          close (unit = 2)
          close (unit = 3)
        end do

      end do
      ! write the results

    5 format(1x, 2d14.6, 2d17.9)
    8 format(1x, 1d14.6, 4d17.9)

      stop
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine initval(ux_p, d2ux_p, mode00, N)

      ! subroutine to read the values
      implicit none
      include 'par.for'
      character*17 nm, nmp1, nmp2
      integer i, j, k, f, N, shift, my_rank, t
      real*8        ux_p(imax/8+1,jmax,kphys,N),
     &            d2ux_p(imax/8+1,jmax,kphys,N),
     &            mode00(imax/8+1,jmax)
      complex*16 u_tanal(ptsx/8+1,jmax,kfour),
     &         u_tanalp1(ptsx/8+1,jmax,kfour),
     &         u_tanalp2(ptsx/8+1,jmax,kfour),
     &              ux_f(imax/8+1,jmax,kfour,N),
     &            d2ux_f(imax/8+1,jmax,kfour,N)

      do f = 1, N
        t = 12800 + (f-1) * 64
        do my_rank = 0, 7
          write(nm,'(a,i0.2,a,i0.5,a)')'pert_',my_rank,'_',t,'.bin'
          write(nmp1,'(a,i0.2,a,i0.5,a)')'pert_',my_rank,'_',t+2,'.bin'
          write(nmp2,'(a,i0.2,a,i0.5,a)')'pert_',my_rank,'_',t+4,'.bin'
          open(unit = 1, file = nm,   form = 'unformatted')
          open(unit = 2, file = nmp1, form = 'unformatted')
          open(unit = 3, file = nmp2, form = 'unformatted')
          read(1) u_tanal
          read(2) u_tanalp1
          read(3) u_tanalp2
          close(unit = 1)
          close(unit = 2)
          close(unit = 3)
c          shift = my_rank * (ptsx / 4 - inter / 4)
          shift = my_rank * ((ptsx-1) / 8 - inter / 8)
          do k = 1, kfour
            do j = 1, jmax
              do i = 1, ptsx/8 + 1
                ux_f(i+shift,j,k,f)   = u_tanal(i,j,k)
                d2ux_f(i+shift,j,k,f) = (        u_tanal(i,j,k)
     &                                  - 2.d0 * u_tanalp1(i,j,k)
     &                                  +        u_tanalp2(i,j,k) ) /
     &                                  ( 4.d0 * dt * dt )
              end do
            end do
          end do
        end do
      end do

      do i = 1, imax/8 + 1
        do j = 1, jmax
          mode00(i,j) = 0.d0
          do f = 1, N
            mode00(i,j) = mode00(i,j) + dreal(ux_f(i,j,1,f))
            ux_f(i,j,1,f) = dcmplx(0.d0,0.d0)
          end do
          mode00(i,j) = mode00(i,j) / dble(N)
        end do
      end do
        
      call f_to_p_cc(ux_p, ux_f, N)
      call f_to_p_cc(d2ux_p, d2ux_f, N)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine polcoe(x, y, coef)

      implicit none
      integer i, j, k, n
      real*8 x(3), y(3), coef(3), b, ff, phi, s(3)

      n = 3
      do i = 1, n
        s(i)    = 0.d0
        coef(i) = 0.d0
      end do
      s(n) = - x(1)
      
      do i = 2, n
        do j = n+1-i, n-1
          s(j) = s(j) - x(i) * s(j+1)
        end do
        s(n) = s(n) - x(i)
      end do
      do j = 1, n
        phi = dble(n)
        do k = n-1, 1, -1
          phi = dble(k) * s(k+1) + x(j) * phi
        end do
        ff = y(j) / phi
        b = 1.d0
        do k = n, 1, -1
          coef(k) = coef(k) + b * ff
          b = s(k) + x(j) * b
        end do
      end do

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine derfax(ddx, fc, N)

      ! first derivatives calculation in x direction
      implicit none
      include 'par.for'
      integer i, f, N
      real*8 fc(imax/8+1,N), ddx(imax/8+1,N)

      do f = 1, N / 2
        do i = 4, imax/2 + 1 - 3
          ddx(i,f) = ( -           fc(i-3,f) + fc(i+3,f)
     &                 +  9.d0 * ( fc(i-2,f) - fc(i+2,f) )
     &                 - 45.d0 * ( fc(i-1,f) - fc(i+1,f) ) ) /
     &                 ( 60.d0 * dx * 4.d0 )
        end do
      end do
      
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine derfaxd2(ddx, fc, N)

      ! first derivatives calculation in x direction
      implicit none
      include 'par.for'
      integer i, f, N
      real*8 fc(imax/8+1,N), ddx(imax/8+1,N)

      do f = 2, N / 2
        do i = 4, imax/2 + 1 - 3
          ddx(i,f) = ( -           fc(i-3,f) + fc(i+3,f)
     &                 +  9.d0 * ( fc(i-2,f) - fc(i+2,f) )
     &                 - 45.d0 * ( fc(i-1,f) - fc(i+1,f) ) ) /
     &                 ( 60.d0 * dx * 4.d0 )
        end do
      end do
      
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine fft1d(a, b, INV, N, M)

      implicit none
      include 'par.for'
      logical INV
      integer j, nby2, l, k, me, mm, i, lpk, M, N 
      real*8 a(n), b(n), tr, ti, sr, si, ninv,
     &       ur(n/2), ui(n/2), wr, wi, sign, t, pibyk

c************************************************************************
c
c   ******************************************************************
c   *                                                                *
c   *                   Subroutine fft1d                             *
c   *                                                                *
c   ******************************************************************
c   Written by G. A. Fenton, Jan. 6, 1988.
c
c                 Single Precision Version
c
c    Computes the discrete Fourier transform of an observed series X(t)
c    which can be expressed as
c
c             X(t_i) = A(t_i) + iB(t_i)
c
c    where A(t_i) and B(t_i) are the real and imaginary parts of X(t) at the
c    time t_i. The vectors A and B are overwritten with the calculated Fourier
c    coefficients A(k) and B(k) corresponding to the frequencies 2*pi*k/(N*s)
c    where s is the time increment. This routine can also perform the inverse
c    Fourier transform by setting the flag INV.
c    If INV = .TRUE. then the original series X(t) is reconstructed from
c    the Fourier coefficients A(k) and B(k).
c    Input and output variables are;
c
c      A   real vector of length N. On input A contains the real part of the
c          observed series. On output A contains the Fourier coefficients A(k)
c          corresponding to the frequencies 2*pi*k/(N*s) where s is the time
c          increment.
c
c      B   real vector of length N. On input B contains the imaginary
c          components of the observed series. On ouput, B contains the Fourier
c          coefficients B(k).
c
c      m   the value of m is such that N = 2**m
c
c      N   the number of data points entered. N must be a power of 2.
c
c    INV   logical flag which is false if the Fourier coefficients are to be
c          calculated, true if the inverse transform is required.
c
c   ur,ui  temporary double precision vectors of length N/2
c
c   Note: this routine is modified from a routine presented by D.E. Newland
c         (originally written by J.W. Cooley et al.)
c
c
c............................................................................

      sign = 1.d0
      if ( .not. INV ) then
        sign = - 1.d0
        ninv =   1.d0 / dble(n)
        do j = 1,N
          a(j) = a(j) * ninv
          b(j) = b(j) * ninv
        end do
      end if

c..........................................................reorder vector
      nby2 = n / 2
      j = 1

      do 30 l = 1, n-1
        if (l .lt. j) then
          t    = a(j)
          a(j) = a(l)
          a(l) = t
          t    = b(j)
          b(j) = b(l)
          b(l) = t
        end if
        k = nby2
  20    if (k .ge. j) go to 30
        j = j - k
        k = k / 2
        go to 20
  30  j = j + k
      
c..........................................Calculate Fourier coefficients
      me = 1
      do mm = 1, m
        k     = me
        pibyk = pi / dble(k)
        me    = 2 * me
        wr    = dcos(pibyk)
        wi    = sign * dsin(pibyk)

        ur(1) = 1.d0
        ui(1) = 0.d0
        do i = 2, k
          ur(i) = ur(i-1) * wr - ui(i-1) * wi
          ui(i) = ur(i-1) * wi + ui(i-1) * wr
        end do

        do j = 1,k
          sr = ur(j)
          si = ui(j)
          do l = j, n, me
            lpk    = l + k
            tr     = a(lpk) * sr - b(lpk) * si
            ti     = b(lpk) * sr + a(lpk) * si
            a(lpk) = a(l) - tr
            b(lpk) = b(l) - ti
            a(l)   = a(l) + tr
            b(l)   = b(l) + ti
          end do
        end do
      end do

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine f_to_p_cc(datap, dataf, N)

      ! transform from Fourier space to Physical space
      ! dataf -> Fourier space (in)
      ! datap -> Physical space (out)
      implicit none
      include 'par.for'
      integer p, p1, p2, p3, p4, n2p3, i, j, k, f, N
      real*8 c1, h1i, h1r, h2i, h2r, wis, wrs, theta,
     &       wi, wpi, wpr, wr, wtemp, datap(imax/8+1,jmax,kphys,N)
      complex*16 dataf(imax/8+1,jmax,kfour,N)

      theta = - 2.d0 * pi / dble(kphys)
      c1    =   0.5d0
      wpr   = - 2.0d0 * dsin(0.5d0 * theta)**2
      wpi   =   dsin(theta)
      
      do f = 1, N
        do j = 1, jmax
          do i = 1, imax/8 + 1
            do k = 0, kfour - 1
              datap(i,j,2*k+1,f) = dreal(dataf(i,j,k+1,f))
              datap(i,j,2*k+2,f) = dimag(dataf(i,j,k+1,f))
            end do
            do k = 2 * kfour + 1, kphys ! if 2*kfour > 2/3 kphys = alias
              datap(i,j,k,f) = 0.d0
            end do
            datap(i,j,1,f) = datap(i,j,1,f) * 2.d0
            wr             = 1.0d0 + wpr
            wi             = wpi
            n2p3           = kphys + 3
            do p = 2, kphys/4
              p1              =   2 * p - 1
              p2              =   p1 + 1
              p3              =   n2p3 - p2
              p4              =   p3 + 1
              wrs             =   wr
              wis             =   wi
              h1r             =   c1 * ( datap(i,j,p1,f)
     &                                 + datap(i,j,p3,f) )
              h1i             =   c1 * ( datap(i,j,p2,f)
     &                                 - datap(i,j,p4,f) )
              h2r             = - c1 * ( datap(i,j,p2,f)
     &                                 + datap(i,j,p4,f) )
              h2i             =   c1 * ( datap(i,j,p1,f)
     &                                 - datap(i,j,p3,f) )
              datap(i,j,p1,f) =   h1r + wrs * h2r - wis * h2i
              datap(i,j,p2,f) =   h1i + wrs * h2i + wis * h2r
              datap(i,j,p3,f) =   h1r - wrs * h2r + wis * h2i
              datap(i,j,p4,f) = - h1i + wrs * h2i + wis * h2r
              wtemp           =   wr
              wr              =   wr * wpr - wi    * wpi + wr
              wi              =   wi * wpr + wtemp * wpi + wi
            end do
            h1r            = datap(i,j,1,f)
            datap(i,j,1,f) = c1 * ( h1r + datap(i,j,2,f) )
            datap(i,j,2,f) = c1 * ( h1r - datap(i,j,2,f) )
            call four1_cc(datap, -1, i, j, f, N)
          end do
        end do
      end do

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine four1_cc(dataff, isig, ii, jj, f, N)

      implicit none
      include 'par.for'
      integer isig, i, istep, z, m, mmax, ii, jj, f, N
      real*8 tempi, tempr, dataff(imax/8+1,jmax,kphys,N),
     &       theta, wi, wpi, wpr, wr, wtemp

      z = 1
      do i = 1, kphys, 2 
        if (z .gt. i) then
          tempr               = dataff(ii,jj,z,f)
          tempi               = dataff(ii,jj,z+1,f)
          dataff(ii,jj,z,f)   = dataff(ii,jj,i,f)
          dataff(ii,jj,z+1,f) = dataff(ii,jj,i+1,f)
          dataff(ii,jj,i,f)   = tempr
          dataff(ii,jj,i+1,f) = tempi
        end if
        m = kphys / 2
    1   if ((m .ge. 2) .and. (z .gt. m)) then
          z = z - m
          m = m / 2
          goto 1
        end if
        z = z + m
      end do 
      mmax = 2
    2 if (kphys .gt. mmax) then
        istep =   2 * mmax
        theta =   2.d0 * pi / dble(isig * mmax)
        wpr   = - 2.d0 * dsin(0.5d0 * theta)**2
        wpi   =   dsin(theta)
        wr    =   1.d0
        wi    =   0.d0
        do m = 1, mmax, 2
          do i = m, kphys, istep
            z                   = i + mmax
            tempr               = wr * dataff(ii,jj,z,f)
     &                          - wi * dataff(ii,jj,z+1,f)
            tempi               = wr * dataff(ii,jj,z+1,f)
     &                          + wi * dataff(ii,jj,z,f)
            dataff(ii,jj,z,f)   = dataff(ii,jj,i,f)   - tempr
            dataff(ii,jj,z+1,f) = dataff(ii,jj,i+1,f) - tempi
            dataff(ii,jj,i,f)   = dataff(ii,jj,i,f)   + tempr
            dataff(ii,jj,i+1,f) = dataff(ii,jj,i+1,f) + tempi
          end do 
          wtemp = wr
          wr    = wr * wpr - wi    * wpi + wr
          wi    = wi * wpr + wtemp * wpi + wi
        end do 
        mmax = istep
        goto 2
      end if

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c                   end of Fourier analysis in time prg                       c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
