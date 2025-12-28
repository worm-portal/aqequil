      subroutine initcb(uarray,nmax)
c
c     This subroutine initializes the character array uarray to blanks
c     over the first nmax positions. Normally, nmax would be the
c     dimension of the 1D uarray. However, nmax could be less than the
c     true dimension. Also, uarray could actually have more than one
c     dimension. For example, if the array is really ub, having
c     dimensions of i1 and i2, this subroutine could be used by calling
c     it in the following manner: call initcb(ub,i1*i2)). Use this
c     subroutine with caution if nmax is not the product of the true
c     (declared) dimensions.
c
c     To initialize a character array to a non-blank string, use
c     EQLIBU/initcv.f instead.
c
c     This subroutine is called by:
c
c       Any
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       uarray = the array
c       nmax   = dimension or pseudo-dimension of uarray
c
c     Output:
c
c       uarray = uarray, with the first nmax positions set to zero
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
      character*(*) uarray(nmax)
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
        uarray(i) = ' '
      enddo
c
      end
