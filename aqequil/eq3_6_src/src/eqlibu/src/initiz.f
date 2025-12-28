      subroutine initiz(iarray,nmax)
c
c     This subroutine initializes the integer array iarray to zero over
c     the first nmax positions. Normally, nmax would be the dimension
c     of the 1D iarray. However, nmax could be less than the true
c     dimension. Also, iarray could actually have more than one
c     dimension. For example, if the array is really ib, which has
c     dimensions of i1 and i2, this subroutine could be used by calling
c     it in the following manner: call initiz(ib,(i1*i2)). Use this
c     subroutine with caution if nmax is not the product of the true
c     (declared) dimensions.
c
c     To initialize an integer array to a non-zero value, use
c     EQLIBU/initav.f instead.
c
c     NOTE: It may be more efficient to just initialize the array
c     to zero in a DO loop in the calling subroutine.
c
c     This subroutine is called by:
c
c       Any
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       iarray = the array
c       nmax   = dimension or pseudo-dimension of iarray
c
c     Output:
c
c       iarray = iarray, with the first nmax positions set to zero
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
      integer iarray(nmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i
c
c-----------------------------------------------------------------------
c
c     Caution: efficiency may be best served by not unrolling the
c     following loop.
c
      do i = 1,nmax
        iarray(i) = 0
      enddo
c
      end
