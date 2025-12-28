      subroutine rscaly(avxmax,avy,avys,eps100,nmax)
c
c     This subroutine rescales the elements in the avys array. The
c     results are placed in the avy array. This subroutine is normally
c     used in conjunction with EQLIBU/scalx1.f and EQLIBU/polfit.f
c     to fitting interpolating polynomials to data in which avx
c     contains values of the independent variable, avxs contains
c     values of the independent variable scaled to the interval
c     (-1,1), avys contains values of the dependent variable as
c     fit with the independent variable scaled, and avy contains
c     values of the dependent variable corrected (rescaled) for
c     scaling of the in dependent variable.
c
c     This subroutine is called by:
c
c       EQPT/intrp.f
c       Any
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       avys   = array whose contents are to be rescaled
c       avxmax = the max norm of avx
c       n      = the number of elements in avys
c
c     Output:
c
c       avy    = array whose contents are rescaled
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nmax
c
      real*8 avy(nmax),avys(nmax)
      real*8 avxmax,eps100
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i
c
      real*8 ax,ay
c
c-----------------------------------------------------------------------
c
c     Note- nmax is expected to be relatively small. Therefore, loop
c     unrolling or calling of subroutines which employ loop unrolling
c     is not used here.
c
c     Test avxmax.
c
      ax = avxmax - 1.0
      if (abs(ax) .le. eps100) then
        do i = 1,nmax
          avy(i) = avys(i)
        enddo
        go to 999
      endif
c
c     Rescale the array.
c
      ay = 1.
      avy(1) = avys(1)
      do i = 2,nmax
        ay = ay/avxmax
        avy(i) = avys(i)*ay
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
