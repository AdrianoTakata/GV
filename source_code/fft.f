ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                             c
c                       Fast Fourier Transformation                           c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine f_to_p(datap, dataf, nx)

      ! transform from Fourier space to Physical space
      ! dataf -> Fourier space (in)
      ! datap -> Physical space (out)
      implicit none
      include 'par.for'
      integer p, p1, p2, p3, p4, n2p3, i, j, k, nx
      real*8 c1, wis, wrs, theta, wi, wpi, wpr, wr, wtemp,
     &       h1i(nx,jmax), h1r(nx,jmax), h2i(nx,jmax), h2r(nx,jmax),
     &       datap(nx,jmax,kphys)
      complex*16 dataf(nx,jmax,kfour)

      theta = - 2.d0 * pi / dble(kphys)
      c1    =   0.5d0
      wpr   = - 2.d0 * dsin(0.5d0 * theta)**2
      wpi   =   dsin(theta)

      do k = 0, kfour-1
        datap(:,:,2*k+1) = dreal(dataf(:,:,k+1))
        datap(:,:,2*k+2) = dimag(dataf(:,:,k+1))
      end do
      datap(:,:,1) = 2.d0 * datap(:,:,1)
      datap(:,:,2) = 0.d0
      do k = 2 * kfour + 1, kphys ! if 2*kfour > 2/3 kphys = alias
        datap(:,:,k) = 0.d0
      end do
      wr         = 1.d0 + wpr
      wi         = wpi
      n2p3       = kphys + 3
      do p = 2, kphys/4
        p1            =   2 * p - 1
        p2            =   p1 + 1
        p3            =   n2p3 - p2
        p4            =   p3 + 1
        wrs           =   wr
        wis           =   wi
        h1r(:,:)      =   c1 * (datap(:,:,p1) + datap(:,:,p3))
        h1i(:,:)      =   c1 * (datap(:,:,p2) - datap(:,:,p4))
        h2r(:,:)      = - c1 * (datap(:,:,p2) + datap(:,:,p4))
        h2i(:,:)      =   c1 * (datap(:,:,p1) - datap(:,:,p3))
        datap(:,:,p1) =   h1r(:,:) + wrs * h2r(:,:) - wis * h2i(:,:)
        datap(:,:,p2) =   h1i(:,:) + wrs * h2i(:,:) + wis * h2r(:,:)
        datap(:,:,p3) =   h1r(:,:) - wrs * h2r(:,:) + wis * h2i(:,:)
        datap(:,:,p4) = - h1i(:,:) + wrs * h2i(:,:) + wis * h2r(:,:)
        wtemp         =   wr
        wr            =   wr * wpr - wi    * wpi + wr
        wi            =   wi * wpr + wtemp * wpi + wi
      end do
      h1r(:,:)     = datap(:,:,1)
      datap(:,:,1) = c1 * (h1r(:,:) + datap(:,:,2))
      datap(:,:,2) = c1 * (h1r(:,:) - datap(:,:,2))
      call four1(datap,-1,nx)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine p_to_f(datap, dataf, nx)

      ! transform from Physical space to Fourier space
      ! datap -> Physical space data (in)
      ! dataf -> Fourier space data (out)
      implicit none
      include 'par.for'
      integer p, p1, p2, p3, p4, n2p3, i, j, k, nx
      real*8 c1, wis, wrs, theta, wi, wpi, wpr, wr, wtemp,
     &       h1i(nx,jmax), h1r(nx,jmax), h2i(nx,jmax), h2r(nx,jmax),
     &       datap(nx,jmax,kphys)
      complex*16 dataf(nx,jmax,kfour)

      theta =   2.d0 * pi / dble(kphys)
      c1    =   0.5d0
      wpr   = - 2.d0 * dsin(0.5d0 * theta)**2
      wpi   =   dsin(theta)

      call four1(datap,+1,nx)
      wr   = 1.d0 + wpr
      wi   = wpi
      n2p3 = kphys + 3
      do p = 2, kphys/4
        p1            =   2 * p - 1
        p2            =   p1 + 1
        p3            =   n2p3 - p2
        p4            =   p3 + 1
        wrs           =   wr
        wis           =   wi
        h1r(:,:)      =   c1 * (datap(:,:,p1) + datap(:,:,p3))
        h1i(:,:)      =   c1 * (datap(:,:,p2) - datap(:,:,p4))
        h2r(:,:)      =   c1 * (datap(:,:,p2) + datap(:,:,p4))
        h2i(:,:)      = - c1 * (datap(:,:,p1) - datap(:,:,p3))
        datap(:,:,p1) =   h1r(:,:) + wrs * h2r(:,:) - wis * h2i(:,:)
        datap(:,:,p2) =   h1i(:,:) + wrs * h2i(:,:) + wis * h2r(:,:)
        datap(:,:,p3) =   h1r(:,:) - wrs * h2r(:,:) + wis * h2i(:,:)
        datap(:,:,p4) = - h1i(:,:) + wrs * h2i(:,:) + wis * h2r(:,:)
        wtemp         =   wr
        wr            =   wr * wpr - wi    * wpi + wr
        wi            =   wi * wpr + wtemp * wpi + wi
      end do
      h1r(:,:)     = datap(:,:,1)
      datap(:,:,1) = ( h1r(:,:) + datap(:,:,2) ) / dble(kphys)
      datap(:,:,2) =   h1r(:,:) - datap(:,:,2)
      do k = 2, kphys
        datap(:,:,k) = 2.d0 * datap(:,:,k) / dble(kphys)
      end do
      do k = 0, kfour - 1
       do j = 1, jmax
        do i = 1, nx
         if (dabs(datap(i,j,2*k+1)).lt.1d-14) datap(i,j,2*k+1) = 0.d0
         if (dabs(datap(i,j,2*k+2)).lt.1d-14) datap(i,j,2*k+2) = 0.d0
        end do
       end do
       dataf(:,:,k+1) = dcmplx(datap(:,:,2*k+1),datap(:,:,2*k+2))
      end do
      dataf(:,:,1) = dcmplx(dreal(dataf(:,:,1)),0.d0)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      subroutine four1(dataff, isig, nx)

      implicit none
      include 'par.for'
      integer isig, k, istep, z, mm, mmax, nx
      real*8 dataff(nx,jmax,kphys), tempi(nx,jmax),tempr(nx,jmax),
     &       theta, wi, wpi, wpr, wr, wtemp

      z = 1
      do k = 1, kphys, 2 
        if (z .gt. k) then
          tempr           = dataff(:,:,z)
          tempi           = dataff(:,:,z+1)
          dataff(:,:,z)   = dataff(:,:,k)
          dataff(:,:,z+1) = dataff(:,:,k+1)
          dataff(:,:,k)   = tempr
          dataff(:,:,k+1) = tempi
        end if
        mm = kphys / 2
        do while ((mm .ge. 2). and. (z .gt. mm))
          z  = z - mm
          mm = mm / 2
        end do
        z = z + mm
      end do
      mmax = 2
      do while (kphys .gt. mmax)
        istep =   2 * mmax
        theta =   2.d0 * pi / dble(isig * mmax)
        wpr   = - 2.d0 * dsin(0.5d0 * theta)**2
        wpi   =   dsin(theta)
        wr    =   1.d0
        wi    =   0.d0
        do mm = 1, mmax, 2
          do k = mm, kphys, istep
            z               = k + mmax
            tempr           = wr * dataff(:,:,z)
     &                      - wi * dataff(:,:,z+1)
            tempi           = wr * dataff(:,:,z+1)
     &                      + wi * dataff(:,:,z)
            dataff(:,:,z)   = dataff(:,:,k)   - tempr
            dataff(:,:,z+1) = dataff(:,:,k+1) - tempi
            dataff(:,:,k)   = dataff(:,:,k)   + tempr
            dataff(:,:,k+1) = dataff(:,:,k+1) + tempi
          end do
          wtemp = wr
          wr    = wr * wpr - wi    * wpi + wr
          wi    = wi * wpr + wtemp * wpi + wi
        end do
        mmax = istep
      end do
      
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine f_to_p_1d(datap, dataf)

      ! transform from Fourier space to Physical space
      ! dataf -> Fourier space (in)
      ! datap -> Physical space (out)
      implicit none
      include 'par.for'
      integer p, p1, p2, p3, p4, n2p3, k
      real*8 c1, h1i, h1r, h2i, h2r, wis, wrs, theta,
     &       wi, wpi, wpr, wr, wtemp, datap(kphys)
      complex*16 dataf(kfour)

      theta = - 2.d0 * pi / dble(kphys)
      c1    =   0.5d0
      wpr   = - 2.d0 * dsin(0.5d0 * theta)**2
      wpi   =   dsin(theta)

      do k = 0, kfour-1
        datap(2*k+1) = dreal(dataf(k+1))
        datap(2*k+2) = dimag(dataf(k+1))
      end do
      datap(1) = 2.d0 * datap(1)
      datap(2) = 0.d0
      do k = 2 * kfour + 1, kphys ! if 2*kfour > 2/3 kphys = alias
        datap(k) = 0.d0
      end do
      wr           = 1.0d0 + wpr
      wi           = wpi
      n2p3         = kphys + 3
      do p = 2, kphys/4
        p1            =   2 * p - 1
        p2            =   p1 + 1
        p3            =   n2p3 - p2
        p4            =   p3 + 1
        wrs           =   wr
        wis           =   wi
        h1r           =   c1 * (datap(p1) + datap(p3))
        h1i           =   c1 * (datap(p2) - datap(p4))
        h2r           = - c1 * (datap(p2) + datap(p4))
        h2i           =   c1 * (datap(p1) - datap(p3))
        datap(p1) =   h1r + wrs * h2r - wis * h2i
        datap(p2) =   h1i + wrs * h2i + wis * h2r
        datap(p3) =   h1r - wrs * h2r + wis * h2i
        datap(p4) = - h1i + wrs * h2i + wis * h2r
        wtemp         =   wr
        wr            =   wr * wpr - wi    * wpi + wr
        wi            =   wi * wpr + wtemp * wpi + wi
      end do
      h1r          = datap(1)
      datap(1) = c1 * (h1r + datap(2))
      datap(2) = c1 * (h1r - datap(2))
      call four1_1d(datap, -1)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine four1_1d(dataff, isig)

      implicit none
      include 'par.for'
      integer isig, i, istep, z, mm, mmax
      real*8 tempi, tempr, dataff(kphys),
     &       theta, wi, wpi, wpr, wr, wtemp

      z = 1
      do i = 1, kphys, 2
        if (z .gt. i) then
          tempr             = dataff(z)
          tempi             = dataff(z+1)
          dataff(z)   = dataff(i)
          dataff(z+1) = dataff(i+1)
          dataff(i)   = tempr
          dataff(i+1) = tempi
        end if
        mm = kphys / 2
    1   if ((mm .ge. 2) .and. (z .gt. mm)) then
          z  = z - mm
          mm = mm / 2
          goto 1
        end if
        z = z + mm
      end do
      mmax = 2
    2 if (kphys .gt. mmax) then
        istep =   2 * mmax
        theta =   2.d0 * pi / dble(isig * mmax)
        wpr   = - 2.d0 * dsin(0.5d0 * theta)**2
        wpi   =   dsin(theta)
        wr    =   1.d0
        wi    =   0.d0
        do mm = 1, mmax, 2
          do i = mm, kphys, istep
            z                 = i + mmax
            tempr             = wr * dataff(z)
     &                        - wi * dataff(z+1)
            tempi             = wr * dataff(z+1)
     &                        + wi * dataff(z)
            dataff(z)   = dataff(i)   - tempr
            dataff(z+1) = dataff(i+1) - tempi
            dataff(i)   = dataff(i)   + tempr
            dataff(i+1) = dataff(i+1) + tempi
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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine p_to_f_2d(datap, dataf)

      ! transform from Physical space to Fourier space
      ! datap -> Physical space data (in)
      ! dataf -> fourier space data (out)
      implicit none
      include 'par.for'
      integer p, p1, p2, p3, p4, n2p3, j, k
      real*8 c1, h1i, h1r, h2i, h2r, wis, wrs, theta,
     &       wi, wpi, wpr, wr, wtemp, datap(jmax,kphys)
      complex*16 dataf(jmax,kfour)

      theta =   2.d0 * pi / dble(kphys)
      c1    =   0.5d0
      wpr   = - 2.d0 * dsin(0.5d0 * theta)**2
      wpi   =   dsin(theta)

      do j = 1, jmax
        call four1_2d(datap,+1,j)
        wr   = 1.0d0 + wpr
        wi   = wpi
        n2p3 = kphys + 3
        do p = 2, kphys/4
          p1          =   2 * p - 1
          p2          =   p1 + 1
          p3          =   n2p3 - p2
          p4          =   p3 + 1
          wrs         =   wr
          wis         =   wi
          h1r         =   c1 * (datap(j,p1) + datap(j,p3))
          h1i         =   c1 * (datap(j,p2) - datap(j,p4))
          h2r         =   c1 * (datap(j,p2) + datap(j,p4))
          h2i         = - c1 * (datap(j,p1) - datap(j,p3))
          datap(j,p1) =   h1r + wrs * h2r - wis * h2i
          datap(j,p2) =   h1i + wrs * h2i + wis * h2r
          datap(j,p3) =   h1r - wrs * h2r + wis * h2i
          datap(j,p4) = - h1i + wrs * h2i + wis * h2r
          wtemp       =   wr
          wr          =   wr * wpr - wi    * wpi + wr
          wi          =   wi * wpr + wtemp * wpi + wi
        end do
        h1r        = datap(j,1)
        datap(j,1) = (h1r + datap(j,2)) / dble(kphys)
        datap(j,2) =  h1r - datap(j,2)
        do k = 2, kphys
          datap(j,k) = 2.d0 * datap(j,k) / dble(kphys)
        end do
        if (dabs(datap(j,1)) .lt. 1d-14) datap(j,1) = 0.d0
        dataf(j,1) = dcmplx(datap(j,1),0.d0)
        do k = 1, kfour - 1
          if (dabs(datap(j,2*k+1)) .lt. 1d-14)
     &       datap(j,2*k+1) = 0.d0
          if (dabs(datap(j,2*k+2)) .lt. 1d-14)
     &       datap(j,2*k+2) = 0.d0
          dataf(j,k+1) = dcmplx(datap(j,2*k+1),datap(j,2*k+2))
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine four1_2d(dataff, isig, jj)

      implicit none
      include 'par.for'
      integer isig, i, istep, z, mm, mmax, jj
      real*8 tempi, tempr, dataff(jmax,kphys),
     &       theta, wi, wpi, wpr, wr, wtemp

      z = 1
      do i = 1, kphys, 2
        if (z.gt.i) then
          tempr          = dataff(jj,z)
          tempi          = dataff(jj,z+1)
          dataff(jj,z)   = dataff(jj,i)
          dataff(jj,z+1) = dataff(jj,i+1)
          dataff(jj,i)   = tempr
          dataff(jj,i+1) = tempi
        end if
        mm = kphys / 2
    1   if ((mm .ge. 2) .and. (z .gt. mm)) then
          z = z - mm
          mm = mm / 2
          goto 1
        end if
        z = z + mm
      end do 
      mmax = 2
    2 if (kphys .gt. mmax) then
        istep =   2 * mmax
        theta =   2.d0 * pi / dble(isig * mmax)
        wpr   = - 2.d0 * dsin(0.5d0 * theta)**2
        wpi   =   dsin(theta)
        wr    =   1.d0
        wi    =   0.d0
        do mm = 1, mmax, 2
          do i = mm, kphys, istep
            z              = i + mmax
            tempr          = wr * dataff(jj,z)   - wi * dataff(jj,z+1)
            tempi          = wr * dataff(jj,z+1) + wi * dataff(jj,z)
            dataff(jj,z)   = dataff(jj,i)   - tempr
            dataff(jj,z+1) = dataff(jj,i+1) - tempi
            dataff(jj,i)   = dataff(jj,i)   + tempr
            dataff(jj,i+1) = dataff(jj,i+1) + tempi
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
c                       Fast Fourier Transformation                           c
c                                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
