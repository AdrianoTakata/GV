cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                   subroutines transforms binary data in asc                  c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program analysis 29052009

      implicit none
      include 'par.for'
      character*30 nome, name_file
      integer i, j, k, t, my_rank, shift, var, cont
      real*8  x, y, z,
     &        uxbs(ptsx,jmax), uybs(ptsx,jmax), wzbs(ptsx,jmax),
     &        thbs(ptsx,jmax),
     &        uxb(imax,jmax), uyb(imax,jmax), wzb(imax,jmax),
     &        thb(imax,jmax),
     &        uxp(imax,jmax,kphys),    wxp(imax,jmax,kphys),
     &        uyp(imax,jmax,kphys),    wyp(imax,jmax,kphys),
     &        uzp(imax,jmax,kphys),    wzp(imax,jmax,kphys),
     &        thp(imax,jmax,kphys),
     &        duxbdy(imax,jmax),       duxbdz(imax,jmax),
     &        dthbdy(imax,jmax),       dthbdz(imax,jmax),
     &        dthmdy(imax),            duxmdy(imax),
     &        duxdyp(imax,jmax,kphys), dthdyp(imax,jmax,kphys)
      complex*16 uxt(imax,jmax,kfour),    wxt(imax,jmax,kfour),
     &           uyt(imax,jmax,kfour),    wyt(imax,jmax,kfour),
     &           uzt(imax,jmax,kfour),    wzt(imax,jmax,kfour),
     &           tht(imax,jmax,kfour),   
     &            ux(ptsx,jmax,kfour),     wx(ptsx,jmax,kfour),
     &            uy(ptsx,jmax,kfour),     wy(ptsx,jmax,kfour),
     &            uz(ptsx,jmax,kfour),     wz(ptsx,jmax,kfour),
     &            th(ptsx,jmax,kfour),  duxdy(imax,jmax,kfour),
     &         dthdy(imax,jmax,kfour)

      call derivs_kt

      ! Boundary layer data
      do my_rank = 0, np - 1
        write(nome,'(a,i0.2,a)')'baseflow2D/based_',my_rank,'.bin'
        open(1, file = nome, form = 'unformatted')
        read(1) uxbs, uybs, wzbs, thbs
        close(unit = 1)
        shift = my_rank * (ptsx - inter - 1)
        do j = 1, jmax
          do i = 1, ptsx
            uxb(i+shift,j) = uxbs(i,j)
            uyb(i+shift,j) = uybs(i,j)
            wzb(i+shift,j) = wzbs(i,j)
            thb(i+shift,j) = thbs(i,j)
          end do
        end do
      end do

      ! Disturbances variables data
      do cont = 0, 15
        write(*,*)'cont=', cont
        do my_rank = 0, np - 1
          if (cont .eq. 0) then
            write(nome,'(a,i0.2,a)')'data_',my_rank,'.bin'
          else
            write(nome,'(a,i0.2,a,i0.2,a)')
     &               'data_',my_rank,'_',cont,'.bin'
          end if
          open(2, file = nome, form = 'unformatted')
          read(2) t
          read(2) ux, uy, uz, wx, wy, wz, th
          close (unit = 2)
          shift = my_rank * (ptsx - inter - 1)
          do k = 1, kfour
            do j = 1, jmax
              do i = 1, ptsx
                uxt(i+shift,j,k) = ux(i,j,k)
                uyt(i+shift,j,k) = uy(i,j,k)
                uzt(i+shift,j,k) = uz(i,j,k)
                wxt(i+shift,j,k) = wx(i,j,k)
                wyt(i+shift,j,k) = wy(i,j,k)
                wzt(i+shift,j,k) = wz(i,j,k)
                tht(i+shift,j,k) = th(i,j,k)
              end do
            end do
          end do
        end do
  
        call f_to_p(uxp, uxt, imax)
        call f_to_p(uyp, uyt, imax)
        call f_to_p(uzp, uzt, imax)
        call f_to_p(wxp, wxt, imax)
        call f_to_p(wyp, wyt, imax)
        call f_to_p(wzp, wzt, imax)
        call f_to_p(thp, tht, imax)


      ! writes data to spacial space to be open by tecplot
c     do i = 201, 738, 67
!     i = 157
!       var = 100 * ( (i - 1) * dx + x0 )
!       write(nome,'(a,i0.4,a)')'crosscutxeq',var,'.dat'
!       open (3, file = nome, status = 'unknown')
!       write(3,*) 'VARIABLES=,"z","y","uxb","uyb","wzb","thb",
!    &    "uxp","uyp","uzp","thp","qq"'
!       write(3,*) 'ZONE I=',kphys+1,', J=',jmax,', F=POINT'
!       
!       do j = 1, jmax
!         if (stf .eq. 1.d0) then
!           y = dble(j-1) * dy
!          else
!           y = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)  
!         end if
!         z = - 1.d0 * 2.d0 * pi / ( dble(kphys) * beta )
!         write(3,5) z, y, uxb(i,j), uyb(i,j), wzb(i,j), thb(i,j),
!    &               uxp(i,j,kphys), uyp(i,j,kphys), uzp(i,j,kphys),
!    &               thp(i,j,kphys), 0.d0
!       
!         do k = 1, kphys
!           z = dble(k-1) * 2.d0 * pi / ( dble(kphys) * beta )
!           write(3,5) z, y, uxb(i,j), uyb(i,j), wzb(i,j), thb(i,j),
!    &                 uxp(i,j,k), uyp(i,j,k), uzp(i,j,k), thp(i,j,k), 0.d0
!         end do
!       end do
!       close (unit = 3)
!     i = 214
!       var = 100 * ( (i - 1) * dx + x0 )
!       write(nome,'(a,i0.4,a)')'crosscutxeq',var,'.dat'
!       open (3, file = nome,status = 'unknown')
!       write(3,*) 'VARIABLES=,"z","y","uxb","uyb","wzb","thb",
!    &    "uxp","uyp","uzp","thp","qq"'
!       write(3,*) 'ZONE I=',kphys+1,', J=',jmax,', F=POINT'
!       
!       do j = 1, jmax
!         if (stf .eq. 1.d0) then
!           y = dble(j-1) * dy
!          else
!           y = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)  
!         end if
!         z = - 1.d0 * 2.d0 * pi / ( dble(kphys) * beta )
!         write(3,5) z, y, uxb(i,j), uyb(i,j), wzb(i,j), thb(i,j),
!    &               uxp(i,j,kphys), uyp(i,j,kphys), uzp(i,j,kphys),
!    &               thp(i,j,kphys), 0.d0
!       
!         do k = 1, kphys
!           z = dble(k-1) * 2.d0 * pi / ( dble(kphys) * beta )
!           write(3,5) z, y, uxb(i,j), uyb(i,j), wzb(i,j), thb(i,j),
!    &                 uxp(i,j,k), uyp(i,j,k), uzp(i,j,k), thp(i,j,k), 0.d0
!         end do
!       end do
!       close (unit = 3)
!     i = 428
!       var = 100 * ( (i - 1) * dx + x0 )
!       write(nome,'(a,i0.4,a)')'crosscutxeq',var,'.dat'
!       open (3, file = nome,status = 'unknown')
!       write(3,*) 'VARIABLES=,"z","y","uxb","uyb","wzb","thb",
!    &    "uxp","uyp","uzp","thp","qq"'
!       write(3,*) 'ZONE I=',kphys+1,', J=',jmax,', F=POINT'
!       
!       do j = 1, jmax
!         if (stf .eq. 1.d0) then
!           y = dble(j-1) * dy
!          else
!           y = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)  
!         end if
!         z = - 1.d0 * 2.d0 * pi / ( dble(kphys) * beta )
!         write(3,5) z, y, uxb(i,j), uyb(i,j), wzb(i,j), thb(i,j),
!    &               uxp(i,j,kphys), uyp(i,j,kphys), uzp(i,j,kphys),
!    &               thp(i,j,kphys), 0.d0
!       
!         do k = 1, kphys
!           z = dble(k-1) * 2.d0 * pi / ( dble(kphys) * beta )
!           write(3,5) z, y, uxb(i,j), uyb(i,j), wzb(i,j), thb(i,j),
!    &                 uxp(i,j,k), uyp(i,j,k), uzp(i,j,k), thp(i,j,k), 0.d0
!         end do
!       end do
!       close (unit = 3)
!     i = 571
!       var = 100 * ( (i - 1) * dx + x0 )
!       write(nome,'(a,i0.4,a)')'crosscutxeq',var,'.dat'
!       open (3, file = nome,status = 'unknown')
!       write(3,*) 'VARIABLES=,"z","y","uxb","uyb","wzb","thb",
!    &    "uxp","uyp","uzp","thp","qq"'
!       write(3,*) 'ZONE I=',kphys+1,', J=',jmax,', F=POINT'
!       
!       do j = 1, jmax
!         if (stf .eq. 1.d0) then
!           y = dble(j-1) * dy
!          else
!           y = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)  
!         end if
!         z = - 1.d0 * 2.d0 * pi / ( dble(kphys) * beta )
!         write(3,5) z, y, uxb(i,j), uyb(i,j), wzb(i,j), thb(i,j),
!    &               uxp(i,j,kphys), uyp(i,j,kphys), uzp(i,j,kphys),
!    &               thp(i,j,kphys), 0.d0
!       
!         do k = 1, kphys
!           z = dble(k-1) * 2.d0 * pi / ( dble(kphys) * beta )
!           write(3,5) z, y, uxb(i,j), uyb(i,j), wzb(i,j), thb(i,j),
!    &                 uxp(i,j,k), uyp(i,j,k), uzp(i,j,k), thp(i,j,k), 0.d0
!         end do
!       end do
!       close (unit = 3)
!     i = 714
!       var = 100 * ( (i - 1) * dx + x0 )
!       write(nome,'(a,i0.4,a)')'crosscutxeq',var,'.dat'
!       open (3, file = nome,status = 'unknown')
!       write(3,*) 'VARIABLES=,"z","y","uxb","uyb","wzb","thb",
!    &    "uxp","uyp","uzp","thp","qq"'
!       write(3,*) 'ZONE I=',kphys+1,', J=',jmax,', F=POINT'
!       
!       do j = 1, jmax
!         if (stf .eq. 1.d0) then
!           y = dble(j-1) * dy
!          else
!           y = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)  
!         end if
!         z = - 1.d0 * 2.d0 * pi / ( dble(kphys) * beta )
!         write(3,5) z, y, uxb(i,j), uyb(i,j), wzb(i,j), thb(i,j),
!    &               uxp(i,j,kphys), uyp(i,j,kphys), uzp(i,j,kphys),
!    &               thp(i,j,kphys), 0.d0
!       
!         do k = 1, kphys
!           z = dble(k-1) * 2.d0 * pi / ( dble(kphys) * beta )
!           write(3,5) z, y, uxb(i,j), uyb(i,j), wzb(i,j), thb(i,j),
!    &                 uxp(i,j,k), uyp(i,j,k), uzp(i,j,k), thp(i,j,k), 0.d0
!         end do
!       end do
!       close (unit = 3)
!     i = 778
!       var = 100 * ( (i - 1) * dx + x0 )
!       write(nome,'(a,i0.4,a)')'crosscutxeq',var,'.dat'
!       open (3, file = nome,status = 'unknown')
!       write(3,*) 'VARIABLES=,"z","y","uxb","uyb","wzb","thb",
!    &    "uxp","uyp","uzp","thp","qq"'
!       write(3,*) 'ZONE I=',kphys+1,', J=',jmax,', F=POINT'
!       
!       do j = 1, jmax
!         if (stf .eq. 1.d0) then
!           y = dble(j-1) * dy
!          else
!           y = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)  
!         end if
!         z = - 1.d0 * 2.d0 * pi / ( dble(kphys) * beta )
!         write(3,5) z, y, uxb(i,j), uyb(i,j), wzb(i,j), thb(i,j),
!    &               uxp(i,j,kphys), uyp(i,j,kphys), uzp(i,j,kphys),
!    &               thp(i,j,kphys), 0.d0
!       
!         do k = 1, kphys
!           z = dble(k-1) * 2.d0 * pi / ( dble(kphys) * beta )
!           write(3,5) z, y, uxb(i,j), uyb(i,j), wzb(i,j), thb(i,j),
!    &                 uxp(i,j,k), uyp(i,j,k), uzp(i,j,k), thp(i,j,k), 0.d0
!         end do
!       end do
!       close (unit = 3)
!     i = 864
!       var = 100 * ( (i - 1) * dx + x0 )
!       write(nome,'(a,i0.4,a)')'crosscutxeq',var,'.dat'
!       open (3, file = nome,status = 'unknown')
!       write(3,*) 'VARIABLES=,"z","y","uxb","uyb","wzb","thb",
!    &    "uxp","uyp","uzp","thp","qq"'
!       write(3,*) 'ZONE I=',kphys+1,', J=',jmax,', F=POINT'
!       
!       do j = 1, jmax
!         if (stf .eq. 1.d0) then
!           y = dble(j-1) * dy
!          else
!           y = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)  
!         end if
!         z = - 1.d0 * 2.d0 * pi / ( dble(kphys) * beta )
!         write(3,5) z, y, uxb(i,j), uyb(i,j), wzb(i,j), thb(i,j),
!    &               uxp(i,j,kphys), uyp(i,j,kphys), uzp(i,j,kphys),
!    &               thp(i,j,kphys), 0.d0
!       
!         do k = 1, kphys
!           z = dble(k-1) * 2.d0 * pi / ( dble(kphys) * beta )
!           write(3,5) z, y, uxb(i,j), uyb(i,j), wzb(i,j), thb(i,j),
!    &                 uxp(i,j,k), uyp(i,j,k), uzp(i,j,k), thp(i,j,k), 0.d0
!         end do
!       end do
!       close (unit = 3)
!     i = 931
!       var = 100 * ( (i - 1) * dx + x0 )
!       write(nome,'(a,i0.4,a)')'crosscutxeq',var,'.dat'
!       open (3, file = nome,status = 'unknown')
!       write(3,*) 'VARIABLES=,"z","y","uxb","uyb","wzb","thb",
!    &    "uxp","uyp","uzp","thp","qq"'
!       write(3,*) 'ZONE I=',kphys+1,', J=',jmax,', F=POINT'
!       
!       do j = 1, jmax
!         if (stf .eq. 1.d0) then
!           y = dble(j-1) * dy
!          else
!           y = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)  
!         end if
!         z = - 1.d0 * 2.d0 * pi / ( dble(kphys) * beta )
!         write(3,5) z, y, uxb(i,j), uyb(i,j), wzb(i,j), thb(i,j),
!    &               uxp(i,j,kphys), uyp(i,j,kphys), uzp(i,j,kphys),
!    &               thp(i,j,kphys), 0.d0
!       
!         do k = 1, kphys
!           z = dble(k-1) * 2.d0 * pi / ( dble(kphys) * beta )
!           write(3,5) z, y, uxb(i,j), uyb(i,j), wzb(i,j), thb(i,j),
!    &                 uxp(i,j,k), uyp(i,j,k), uzp(i,j,k), thp(i,j,k), 0.d0
!         end do
!       end do
!       close (unit = 3)
!     i = 1005
!       var = 100 * ( (i - 1) * dx + x0 )
!       write(nome,'(a,i0.4,a)')'crosscutxeq',var,'.dat'
!       open (3, file = nome,status = 'unknown')
!       write(3,*) 'VARIABLES=,"z","y","uxb","uyb","wzb","thb",
!    &    "uxp","uyp","uzp","thp","qq"'
!       write(3,*) 'ZONE I=',kphys+1,', J=',jmax,', F=POINT'
!       
!       do j = 1, jmax
!         if (stf .eq. 1.d0) then
!           y = dble(j-1) * dy
!          else
!           y = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)  
!         end if
!         z = - 1.d0 * 2.d0 * pi / ( dble(kphys) * beta )
!         write(3,5) z, y, uxb(i,j), uyb(i,j), wzb(i,j), thb(i,j),
!    &               uxp(i,j,kphys), uyp(i,j,kphys), uzp(i,j,kphys),
!    &               thp(i,j,kphys), 0.d0
!       
!         do k = 1, kphys
!           z = dble(k-1) * 2.d0 * pi / ( dble(kphys) * beta )
!           write(3,5) z, y, uxb(i,j), uyb(i,j), wzb(i,j), thb(i,j),
!    &                 uxp(i,j,k), uyp(i,j,k), uzp(i,j,k), thp(i,j,k), 0.d0
!         end do
!       end do
!       close (unit = 3)
c     end do

    5 format(1x, 2d14.6, 9d17.9)

        call deryt(duxdy, uxt)
        call deryb(duxbdy, uxb)
        call deryt(dthdy, tht)
        call deryb(dthbdy, thb)
  
        call f_to_p(duxdyp, duxdy, imax)
        call f_to_p(dthdyp, dthdy, imax)
  
        do i = 1, imax
          duxmdy(i) = 0.d0
          dthmdy(i) = 0.d0
        end do
        do i = 1, imax
          do k = 1, kphys
            duxmdy(i) = duxmdy(i) + duxdyp(i,1,k)
            dthmdy(i) = dthmdy(i) + dthdyp(i,1,k)
          end do
          duxmdy(i) = duxmdy(i) / dble(kphys) + duxbdy(i,1)
          dthmdy(i) = dthmdy(i) / dble(kphys) + dthbdy(i,1)
        end do

        write(name_file, '(a, i0.2,a)')'heat_coeffs_',cont,'.dat'
        open (4, file = name_file, status = 'unknown')
        write(4,*) 'VARIABLES=,"x","z","dudy","dthdy","dubdy","dthbdy",
     &  "dumdy","dthmdy"'
        write(4,*) 'ZONE I=',imax,', K=',kphys+1,', F=POINT'

        z = - 1.d0 * 2.d0 * pi / ( dble(kphys) * beta )
        do i = 1, imax
          x = x0 + dble(i-1) * dx
          write(4,6) x, z, duxdyp(i,1,kphys)+duxbdy(i,1), 
     &               dthdyp(i,1,kphys)+dthbdy(i,1),
     &               duxbdy(i,1), dthbdy(i,1),duxmdy(i), dthmdy(i)
        end do

        do k = 1, kphys
          z = dble(k-1) * 2.d0 * pi / ( dble(kphys) * beta )
          do i = 1, imax
            x = x0 + dble(i-1) * dx
            write(4,6) x, z, duxdyp(i,1,k)+duxbdy(i,1), 
     &                 dthdyp(i,1,k)+dthbdy(i,1),
     &                 duxbdy(i,1), dthbdy(i,1),duxmdy(i), dthmdy(i)
          end do
        end do
      end do
    6 format(1x, 2d14.6, 6d17.9)

      stop
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine derivs_kt

      ! calculate the tri-diagonal LHS for the derivative calculation
      implicit none
      include 'par.for'
      real*8  a1y(jmax), b1y(jmax), c1y(jmax)
      common/der1y/ a1y, b1y, c1y

      call coeft(a1y, b1y, c1y, jmax)

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

      subroutine deryb(ddy, fc)

      ! first derivative calculation in y direction
      implicit none
      include 'par.for'
      integer i, j, k
      real*8 a1y(jmax), b1y(jmax), c1y(jmax), rhs(jmax), u(jmax),
     &       fc(imax,jmax), ddy(imax,jmax)
      common/der1y/ a1y, b1y, c1y

      do i = 1, imax
        call rhsyb(fc, rhs, i)
        call tridyb(rhs, a1y, b1y, c1y, u)
        do j = 1, jmax
          ddy(i,j) = u(j)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsyt(fc, rhs, i, k)

      ! RHS for the first derivative calculation in y direction
      implicit none
      include 'par.for'
      integer i, j, k
      complex*16 rhs(jmax), fc(imax,jmax,kfour)

      rhs(1) = ( - 74.d0 * fc(i,1,k) + 16.d0 * fc(i,2,k)
     &           + 72.d0 * fc(i,3,k) - 16.d0 * fc(i,4,k)
     &           +  2.d0 * fc(i,5,k) ) / ( 24.d0 * dy )

      rhs(2) = ( - 406.d0 * fc(i,1,k) - 300.d0 * fc(i,2,k)
     &           + 760.d0 * fc(i,3,k) -  80.d0 * fc(i,4,k)
     &           +  30.d0 * fc(i,5,k) -   4.d0 * fc(i,6,k) ) /
     &           ( 120.d0 * dy )

      do j = 3, jmax - 2
        rhs(j) = (           fc(i,j+2,k) - fc(i,j-2,k)
     &           + 28.d0 * ( fc(i,j+1,k) - fc(i,j-1,k) ) ) /
     &           ( 12.d0 * dy )
      end do

      rhs(jmax-1) = ( - 406.d0 * fc(i,jmax,k)
     &                - 300.d0 * fc(i,jmax-1,k)
     &                + 760.d0 * fc(i,jmax-2,k)
     &                -  80.d0 * fc(i,jmax-3,k)
     &                +  30.d0 * fc(i,jmax-4,k)
     &                -   4.d0 * fc(i,jmax-5,k) ) /
     &              ( - 120.d0 * dy )

      rhs(jmax) = ( - 74.d0 * fc(i,jmax,k)   + 16.d0 * fc(i,jmax-1,k)
     &              + 72.d0 * fc(i,jmax-2,k) - 16.d0 * fc(i,jmax-3,k)
     &              +  2.d0 * fc(i,jmax-4,k) ) / ( - 24.d0 * dy )

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsyb(fc, rhs, i)

      ! RHS for the first derivative calculation in y direction
      implicit none
      include 'par.for'
      integer i, j, k
      real*8 rhs(jmax), fc(imax,jmax)

      rhs(1) = ( - 74.d0 * fc(i,1) + 16.d0 * fc(i,2)
     &           + 72.d0 * fc(i,3) - 16.d0 * fc(i,4)
     &           +  2.d0 * fc(i,5) ) / ( 24.d0 * dy )

      rhs(2) = ( - 406.d0 * fc(i,1) - 300.d0 * fc(i,2)
     &           + 760.d0 * fc(i,3) -  80.d0 * fc(i,4)
     &           +  30.d0 * fc(i,5) -   4.d0 * fc(i,6) ) /
     &           ( 120.d0 * dy )

      do j = 3, jmax - 2
        rhs(j) = (           fc(i,j+2) - fc(i,j-2)
     &           + 28.d0 * ( fc(i,j+1) - fc(i,j-1) ) ) /
     &           ( 12.d0 * dy )
      end do

      rhs(jmax-1) = ( - 406.d0 * fc(i,jmax)   - 300.d0 * fc(i,jmax-1)
     &                + 760.d0 * fc(i,jmax-2) -  80.d0 * fc(i,jmax-3)
     &                +  30.d0 * fc(i,jmax-4) -   4.d0 * fc(i,jmax-5) )
     &            / ( - 120.d0 * dy )

      rhs(jmax) = ( - 74.d0 * fc(i,jmax)   + 16.d0 * fc(i,jmax-1)
     &              + 72.d0 * fc(i,jmax-2) - 16.d0 * fc(i,jmax-3)
     &              +  2.d0 * fc(i,jmax-4) ) / ( - 24.d0 * dy )

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tridyt(rhs, a, b, c, u)

      ! solves the tridiagonal matrix for the derivatives in y direction
      implicit none
      include 'par.for'
      integer j
      real*8 a(jmax), b(jmax), c(jmax), gam(jmax), bet
      complex*16 rhs(jmax), u(jmax)

      bet  = b(1)
      u(1) = rhs(1) / bet
      do j = 2, jmax
        gam(j) = c(j-1) / bet
        bet    = b(j) - a(j) * gam(j)
        u(j)   = ( rhs(j) - a(j) * u(j-1) ) / bet
      end do
      do j = jmax - 1, 1, -1
        u(j) = u(j) - gam(j+1) * u(j+1)
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tridyb(rhs, a, b, c, u)

      ! solves the tridiagonal matrix for the derivatives in y direction
      implicit none
      include 'par.for'
      integer j
      real*8 a(jmax), b(jmax), c(jmax), gam(jmax), bet,
     &       rhs(jmax), u(jmax)

      bet  = b(1)
      u(1) = rhs(1) / bet
      do j = 2, jmax
        gam(j) = c(j-1) / bet
        bet    = b(j) - a(j) * gam(j)
        u(j)   = ( rhs(j) - a(j) * u(j-1) ) / bet
      end do
      do j = jmax - 1, 1, -1
        u(j) = u(j) - gam(j+1) * u(j+1)
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine coeft(a, b, c, lmax)

      ! mount the LHS of the matrix for the first derivative
      implicit none
      integer l, lmax
      real*8 a(lmax), b(lmax), c(lmax)
      
      a(1)      = 0.d0
      b(1)      = 1.d0
      c(1)      = 4.d0

      a(2)      = 1.d0
      b(2)      = 6.d0
      c(2)      = 2.d0

      do l = 3, lmax - 2
        a(l)    = 1.d0
        b(l)    = 3.d0
        c(l)    = 1.d0
      end do

      a(lmax-1) = 2.d0
      b(lmax-1) = 6.d0
      c(lmax-1) = 1.d0

      a(lmax)   = 4.d0
      b(lmax)   = 1.d0
      c(lmax)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                 end subroutines transforms binary data in asc                c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
