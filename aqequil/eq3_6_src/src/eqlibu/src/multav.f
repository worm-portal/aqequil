      subroutine multav(array1,nmax,avalue,array2)
c
c     This subroutine multiplies the first nmax elements of the real*8
c     array array1 by the real*8 value avalue. The result is placed
c     in the real*8 array array2. Normally, nmax would be the dimension
c     of both arrays. However, nmax could be less than the true
c     dimension. Also, the arrays could actually have more than one
c     dimension. For example, if array1 and array2 are really a and b,
c     respectively, having common dimensions of i1 and i2, this
c     subroutine could be used by calling it in the following manner:
c     call multav(a,(i1*i2),avalue,b). Use this subroutine with caution
c     if the two arrays do not have identical dimensioning or nmax is
c     not the product of the true (declared) dimensions.
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
c       avalue = the value by which to multipy array1
c
c     Output:
c
c       array2 = the second array, the product of avalue and array1
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
      real*8 avalue
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ileft
c
      real*8 av
c
c-----------------------------------------------------------------------
c
c     Note the assignment of the value to a local variable.
c     Note also that the loop is unrolled.
c
      av = avalue
      ileft = (nmax/8)*8
c
      do i = 1,ileft,8
        array2(i) = av*array1(i)
        array2(i + 1) = av*array1(i + 1)
        array2(i + 2) = av*array1(i + 2)
        array2(i + 3) = av*array1(i + 3)
        array2(i + 4) = av*array1(i + 4)
        array2(i + 5) = av*array1(i + 5)
        array2(i + 6) = av*array1(i + 6)
        array2(i + 7) = av*array1(i + 7)
      enddo
c
      do i = ileft + 1,nmax
        array2(i) = av*array1(i)
      enddo
c
      end
