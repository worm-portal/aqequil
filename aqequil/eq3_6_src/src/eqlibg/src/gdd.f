      subroutine gdd(dgpit,fxi,gpit,ipbtmx,napmax,napt,palpha)
c
c     This subroutine computes the functions g(xi) and their ionic
c     strength derivatives (xi = alpha(i) * I**1/2). These functions
c     are used in Pitzer's equations. Note: alpha(i) = 0 implies
c     that g(xi) and its derivatives are zero.
c
c     This subroutine is called by:
c
c       EQLIBG/gcoeff.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       fxi    = the ionic strength (the 2nd-order electrostatic
c                  moment function I)
c       napt   = number of alpha coefficient pairs in the palpha array
c
c     Principal output:
c
c       gpit   = array of values for the g(x) function
c       dgpit  = array of values for the ionic strength derivatives
c                  of the g(x) function
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ipbtmx,napmax
c
      integer napt
c
      real(8) dgpit(2,ipbtmx,napmax),gpit(ipbtmx,napmax),
     $ palpha(ipbtmx,napmax)
      real(8) fxi
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,nap
c
      real(8) alphai,dxdi,emx,gpx,gppx,x,xsq
c
c-----------------------------------------------------------------------
c
      do nap = 1,napt
        do i = 1,ipbtmx
          gpit(i,nap) = 0.
          dgpit(1,i,nap) = 0.
          dgpit(2,i,nap) = 0.
          alphai = palpha(i,nap)
          x = alphai * sqrt(fxi)
c
          if (x .ne. 0.) then
            xsq = x*x
            emx = exp(-x)
c
c           Calculate g(x).
c
            gpit(i,nap) = (2./xsq) * (1. - emx*(1. + x))
c
c           Calculate the derivatives with respect to x.
c
            gpx = (-2./(xsq*x)) * (2. - emx*(2. + x*(2. + x)))
            gppx = (-2./(xsq*xsq))
     $      * (-6. + emx*(6. + x*(6. + x*(3. + x))))
c
c           Calculate the derivatives with respect to ionic strength.
c
            dxdi = alphai/(2.*sqrt(fxi))
            dgpit(1,i,nap) = dxdi * gpx
            dgpit(2,i,nap) = dxdi * (gppx*dxdi - 0.5*gpx/fxi)
          endif
        enddo
      enddo
c
      end
