      integer imax, jmax, kfour, kphys, tt, i0, i1, i2, i3, i4, i5,
     &        stpp, start, msh, ptsx, sonda, media, read_media, stencil,
     &        meshpx, meshdx, tt_base, my_form, np, inter, type_cur,
     &        type_sec
      real*8 Re, dx, dxx, dyy, dy, dt, x0, omega, alpha, beta, dt_base,
     &       dz, stf, alphaf, af, bf, cf, df, Pr, U_1, Radius, Go_base,
     &       L_1, N_1, pi, fac_y, C_s, C_w, delta_les, lambda_z, amp_gv,
     &       omega_gv, fhz_gv, Go2Re_base, fhz_sec, amp_sec
      complex*16 im

c     my_form = 0 -> Zero Curvature 
c     my_form = 1 -> Gortler vortices simulations without heat transfer
c     my_form = 2 -> Gortler vortices simulations with heat transfer
      parameter ( my_form = 2 )

c     type_cur = 0 -> varg = kc (Constant)
c     type_cur = 1 -> varg = - kc*tanh(3*(x-8))
c     type_cur = 2 -> varg = kc * (1/2)*(1-tanh(3(x-10.24)))
c     type_cur = 3 -> varg = kc * 2*sin((5*x+3)pi/48)
      parameter( type_cur = 0 )

c     type_sec = 0 -> Varicoso
c     type_sec = 1 -> Sinuoso
      parameter (type_sec = 0)

c     start the program from t = 0 (start=0) or from a given time t (start=1)
      parameter ( start = 1 )

      ! pi number
      parameter ( pi = 3.141592653589793d0 )

      ! imaginary number (sqrt(-1))
      parameter ( im = (0.d0,1.d0) )

      ! Characteristic velocity (dimensional [m / s])
      parameter ( U_1 = 5.0d0 )

      ! lenght scale parameter L (dimensional [m])
      parameter ( L_1 = 0.1d0 )

      ! dynamic viscosity (dimensional [m^2 / s])
      parameter ( N_1 = 1.5094795344339218d-005 )

      ! Curvature wall Radius (dimensional [m] )
      parameter ( Radius = 3.2d0 )

      ! spanwise wavelength (dimensional [m])
      parameter ( lambda_z = 0.009d0 )

      ! value of beta (spanwise wavenumber) (2 * pi / lambda_z (non_dimensional) )
      parameter ( beta = (2.d0 * pi) / (lambda_z / L_1) )

      ! nondimensionalization parameter in the y-direction (Re or 1.d0)
      parameter ( fac_y = 1.d0 )
c     parameter ( fac_y = Re )

      ! Reynolds Number, Prandtl Number and initial x of the domain
      ! Görtler Number and curvature parameter
      parameter ( Re = U_1*L_1/N_1)
      parameter ( Pr = 0.72d0, x0 = 1.d0 )
      parameter ( Go_base = Re**0.25d0*dsqrt(L_1/Radius) )

      ! Curvature term in baseflow Go_base
      parameter ( Go2Re_base =  0.0d0 ) 
c      parameter ( Go2Re_base = (Go_base*Go_base)/dsqrt(Re) )

      ! dimensional frequency for unsteady disturbance F_gv(Hz)
      parameter ( fhz_gv = 0.0d0 )
      parameter ( fhz_sec = 0.0d0 )

      ! alpha (streamwise wavelength) and omega (frequency) omega = F(Hz)*2*pi*L/U
c      parameter ( alpha = 22.6d0 )
      parameter ( alpha = 0.45723605d0 )
      parameter ( omega = fhz_sec * 2.d0 * pi * L_1 / U_1 )
      parameter ( omega_gv = (fhz_gv * 2.0d0 * pi * L_1) / U_1 )

      ! amplitude of the disturbance
      parameter ( amp_gv = 3.9d-3 )
      parameter ( amp_sec = 1.0d-5 )

      ! number of points in y-direction, delta y (dy / sqrt(re * x0))
      parameter ( jmax = 257,  dy = 5.0d-4 )
      parameter ( dyy = dy * dy, stf = 1.01d0 )

      ! number of processing elements
      parameter ( np = 10 )

c     number of points in x direction, delta x, dx^2 and ptsx
      parameter ( imax = 1305, dx = 1.0d-2 )
      parameter ( dxx = dx * dx, ptsx = (imax + (np - 1) * 25) / np )

      ! initial and end points of inflow damping zone, and disturbance strip
      parameter ( i0 = 5, i1 = 10, i2 = 20, i3 = 52 )

      ! initial and end points of outflow damping zone in x-direction
      parameter ( i4 = imax-90, i5 = imax-45 )

      ! steps per period, number of time steps and time step  (2 * pi / omega / stpp)
      parameter ( stpp = 256, tt = 400 * stpp )
      parameter ( dt = 1.0d-3 )
c      parameter ( dt = (2.d0 * pi) / omega_gv / stpp )
c      parameter ( dt = (2.0d0 * pi) / omega / stpp)

      ! number of meshes used in the multigrid solver
      parameter ( msh = 4 )

      ! number of modes in Fourier and physical space, and delta z
      parameter ( kfour = 21, kphys = 64 )
      parameter ( dz = (2.d0 * pi) / (kphys * beta) )

      ! for baseflow use these parameters
      parameter ( dt_base = 5.0d-2 * dx )
      parameter ( tt_base = 1000 *imax ) 

      ! filter constants (Lele C.2.5)
      parameter ( alphaf = 0.48d0 )
      parameter ( af = (   11.d0 + 10.d0 * alphaf) / 16.d0 )
      parameter ( bf = (   15.d0 + 34.d0 * alphaf) / 64.d0 ) !/2
      parameter ( cf = ( -  3.d0 +  6.d0 * alphaf) / 32.d0 ) !/2
      parameter ( df = (    1.d0 -  2.d0 * alphaf) / 64.d0 ) !/2

c     parameters used in parallelized Poisson subroutines
      parameter ( stencil = 5 )
      parameter ( inter   = 2**( msh - 1 ) * ( stencil - 2 ) )
      parameter ( meshdx  = 12, meshpx  = ( stencil - 1 ) / 2 )
