      subroutine initcv(uarray,nmax,uvalue)
c
c     This subroutine initializes the character array uarray to the
c     string contained in uvalue over the first nmax positions.
c     Normally, nmax would be the dimension of the 1D uarray. However,
c     nmax could be less than the true dimension. Also, uarray could
c     actually have more than one dimension. For example, if the array
c     is really ua, which has dimensions of i1 and i2, this subroutine
c     could be used by calling it in the following manner:
c     call initcv((ua,(i1*i2),uv). Use this subroutine with caution if
c     nmax is not the product of the true (declared) dimensions.
c
c     To initialize a character array to blanks, use EQLIBU/initcb.f
c     instead.
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
c       uvalue = the assigned value
c
c     Output:
c
c       uarray = uarray, with the first nmax positions set to the
c                  string contained in uv
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
      character*(*) uarray(nmax),uvalue
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ileft
c
      character*80 uv
c
c-----------------------------------------------------------------------
c
c     Note the assignment of the value to a local variable.
c     Note that the loop is unrolled.
c
      uv = uvalue
      ileft = (nmax/8)*8
c
      do i = 1,ileft,8
        uarray(i) = uv
        uarray(i + 1) = uv
        uarray(i + 2) = uv
        uarray(i + 3) = uv
        uarray(i + 4) = uv
        uarray(i + 5) = uv
        uarray(i + 6) = uv
        uarray(i + 7) = uv
      enddo
c
      do i = ileft + 1,nmax
        uarray(i) = uv
      enddo
c
      end
