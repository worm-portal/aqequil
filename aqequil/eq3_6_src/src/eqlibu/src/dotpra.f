      real*8 function dotpra(array1,array2,nmax)
c
c     This subroutine computes the dot product of two real*8 arrays
c     array1 and array2.
c
c     This subroutine is nearly identical to EQLIBU/ddot.f, which has
c     the same function, but a different calling sequence.
c
c     This subroutine is called by:
c
c       EQLIBU/sgeco.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       array1 =  real*8 vector with nmax elements
c       array2 =  real*8 vector with nmax elements
c
c     Output:
c
c       dotpra = the dot product
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
      real*8 dx
c
c-----------------------------------------------------------------------
c
c     Note the use of a local variable (ax) within the loop.
c     Note also that the loop is unrolled.
c
      dx = 0.
      ileft = (nmax/8)*8
c
      do i = 1,ileft,8
        dx = dx + array1(i)*array2(i) + array1(i + 1)*array2(i + 1)
     $  + array1(i + 2)*array2(i + 2) + array1(i + 3)*array2(i + 3)
     $  + array1(i + 4)*array2(i + 4) + array1(i + 5)*array2(i + 5)
     $  + array1(i + 6)*array2(i + 6) + array1(i + 7)*array2(i + 7)
      enddo
c
      do i = ileft + 1,nmax
        dx = dx + array1(i)*array2(i)
      enddo
c
      dotpra = dx
c
      end
