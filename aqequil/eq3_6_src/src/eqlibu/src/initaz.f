      subroutine initaz(array,nmax)
c
c     This subroutine initializes the real*8 array "array" to zero over
c     the first nmax positions. Normally, nmax would be the dimension
c     of the 1D "array". However, nmax could be less than the true
c     dimension. Also, array could actually have more than one
c     dimension. For example, if the array is really b, which has
c     dimensions of i1 and i2, this subroutine could be used by calling
c     it in the following manner: call initaz(b,(i1*i2)). Use this
c     subroutine with caution if the two arrays do not have identical
c     dimensioning or nmax is not the product of the true (declared)
c     dimensions.
c
c     To initialize a real*8 array to a non-zero value, use
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
c       array  = the array
c       nmax   = dimension or pseudo-dimension of array
c
c     Output:
c
c       array  = array, with the first nmax positions set to zero
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
      real*8 array(nmax)
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
        array(i) = 0.0
      enddo
c
      end
