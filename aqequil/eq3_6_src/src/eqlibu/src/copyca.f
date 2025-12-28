      subroutine copyca(uarr1,uarr2,nmax)
c
c     This subroutine copies the first nmax elements of the character
c     array uarr1 into the character array uarr2. Normally, nmax would
c     be the dimension of both arrays. However, nmax could be less
c     than the true dimension. Also, the arrays could actually have
c     more than one dimension. For example, if uarr1 and uarr2
c     are really a and b, respectively, having common dimensions of
c     i1 and i2, this subroutine could be used by calling it in the
c     following manner: call copyca(a,b,(i1*i2)). Use this subroutine
c     with caution if the two arrays do not have identical dimensioning,
c     nmax is not the product of the actual (declared) dimensions, or
c     the character lengths of the two arrays are not identical.
c
c     This subroutine is called by:
c
c       Any
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       uarr1  = the first array
c       nmax   = dimension or pseudo-dimension of uarr1 and uarr2
c
c     Output:
c
c       uarr2  = the second array
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
      character*(*) uarr1(nmax),uarr2(nmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ileft
c
c-----------------------------------------------------------------------
c
c     Note that the loop is unrolled.
c
      ileft = (nmax/8)*8
c
      do i = 1,ileft,8
        uarr2(i) = uarr1(i)
        uarr2(i + 1) = uarr1(i + 1)
        uarr2(i + 2) = uarr1(i + 2)
        uarr2(i + 3) = uarr1(i + 3)
        uarr2(i + 4) = uarr1(i + 4)
        uarr2(i + 5) = uarr1(i + 5)
        uarr2(i + 6) = uarr1(i + 6)
        uarr2(i + 7) = uarr1(i + 7)
      enddo
c
      do i = ileft + 1,nmax
        uarr2(i) = uarr1(i)
      enddo
c
      end
