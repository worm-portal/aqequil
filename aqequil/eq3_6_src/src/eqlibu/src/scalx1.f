      subroutine scalx1(avx,avxmax,avxs,ier,nmax)
c
c     This subroutine scales the elements in the avx array to the
c     interval (-1,1). The results are placed in the avxs array. This
c     subroutine is normally used in conjunction with EQLIBU/rscaly.f
c     and EQLIBU/polfit.f to fitting interpolating polynomials to
c     data in which avx contains values of the independent variable.
c
c     This subroutine is called by:
c
c       EQPT/eqpt.f
c       Any
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       avx    = array whose contents are to be scaled
c       avxmax = the max norm of avx
c       n      = the number of elements in avx
c
c     Output:
c
c       avxs   = array whose contents are to be scaled
c       ier    = error flag
c                  = 0  okay
c                  = 1  all the elements of avx are zero
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ier,nmax
c
      real*8 avx(nmax),avxs(nmax)
      real*8 avxmax
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i
c
      real*8 ax
c
c-----------------------------------------------------------------------
c
c     Initialize error flag.
c
      ier = 0
c
c     Note- nmax is expected to be relatively small. Therefore, loop
c     unrolling or calling of subroutines which employ loop unrolling
c     is not used here.
c
c     Find the max norm of avx.
c
      avxmax = 0.
      do i = 1,nmax
        ax = abs(avx(i))
        if (ax .gt. avxmax) avxmax = ax
      enddo
c
      if (avxmax .le. 0.) then
        ier = 1
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Scale the array.
c
      do i = 1,nmax
        avxs(i) = avx(i)/avxmax
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
