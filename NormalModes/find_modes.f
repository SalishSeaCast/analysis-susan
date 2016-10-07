      program find_modes

      parameter (max_no_depths = 5000)

      parameter (npts = max_no_depths)

      ! these parameters only influence the size of a, not the mode shape
      parameter (g = 9.8)
      parameter (rho0 = 1025.)
      parameter (small = 1e-17)
      parameter (f = 1e-4)

      real*8 zin(max_no_depths), den(max_no_depths), N(max_no_depths)
      real*8 N2

      real*8 p(npts), z, aN, dN, dp(npts)

      open (11,file='density.dat')   ! input file
      open (12,file="assign4.dat")   ! output file

      read (11,*) zin(1), den(1)
      read (11,*) zin(2), den(2)
      N2 = g * (den(2) - den(1)) / (rho0 * (zin(2) - zin(1)))
      if (N2 .lt. small) then
         N(2) = small
      else
         N(2) = sqrt(N2)
      endif
      do i = 3, max_no_depths
         read (11, *, END=999) zin(i),den(i)
         N2 = g*(den(i)-den(i-1))/(rho0 * (zin(i) - zin(i-1)))
         if (N2.lt.small) then
            N(i) = small
         else
            N(i) = sqrt(N2)
         endif
      enddo
 999  id = i - 1
      depth = zin(id)
      N(id+1) = N(id)
      N(1) = N(2)
      dz = depth/(id-1)

  100 write (*,*) 'Guess the lengthscale a?'
      write (*,*) "If you enter -1, Pi profile from last a value"
      write (*,*) " entered will be printed"
      read (*,*) a
      af = a*f
      if (a.gt.0) then
         p(id+1) = 1.
         dp(id+1) = 0.

         pmax = 1.
         do i=id+1, 2, -1
c     find the indexes of the N values above and below this depth
            Nbelow = i+1
            Nabove = i
c     linearly interpolate N
            aN = 0.5 * (N(Nabove) + N(Nbelow))
            dN = (N(Nabove) - N(Nbelow))/dz
            dp(i-1) = dp(i) + dz * (2 / aN * dN *dp(i)
     >              - p(i) * (aN * aN / af**2))
            p(i-1) = p(i) + dz * dp(i)
            if (abs(p(i-1)).gt.pmax) pmax = abs(p(i-1))
         enddo

         write (*,*) "Your a value", a,' m',
     >        " gives a dpi/dz at the surface of ",dp(1)/pmax

         goto 100
      else
         do i=1,npts
            z = dz*(i-1.)
            write (12,*) z, p(i)
         enddo
      endif

      end




