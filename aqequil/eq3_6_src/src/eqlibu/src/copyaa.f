      subroutine copyaa(array1,array2,nmax)
c
c     This subroutine copies the first nmax elements of the real*8
c     array array1 into the real*8 array array2. Normally, nmax would
c     be the dimension of both arrays. However, nmax could be less
c     than the true dimension. Also, the arrays could actually have
c     more than one dimension. For example, if array1 and array2
c     are really a and b, respectively, having common dimensions of
c     i1 and i2, this subroutine could be used by calling it in the
c     following manner: call copyaa(a,b,(i1*i2)). Use this subroutine
c     with caution if the two arrays do not have identical dimensioning
c     or nmax is not the product of the true (declared) dimensions.
c
c     This subroutine is called by:
c
c       Any
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       array1 = the first array
c       nmax   = dimension or pseudo-dimension of array1 and array2
c
c     Output:
c
c       array2 = the second array
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
      real*8 array1(nmax),array2(nmax)
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
        array2(i) = array1(i)
        array2(i + 1) = array1(i + 1)
        array2(i + 2) = array1(i + 2)
        array2(i + 3) = array1(i + 3)
        array2(i + 4) = array1(i + 4)
        array2(i + 5) = array1(i + 5)
        array2(i + 6) = array1(i + 6)
        array2(i + 7) = array1(i + 7)
      enddo
c
      do i = ileft + 1,nmax
        array2(i) = array1(i)
      enddo
c
      end
