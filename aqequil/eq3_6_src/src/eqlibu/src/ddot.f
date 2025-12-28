      real*8 function ddot(nmax,array1,incx,array2,incy)
c
c     This subroutine computes the dot product (ddot) of the two real*8
c     arrays "array1" and "array2". This subroutine is an adaptation of
c     the 1979 Linpack subroutine of the same name. The size and order
c     in the original calling sequence has been preserved:
c     (n,dx,incx,dy,incy) = (nmax,array1,incx,array2,incy).
c     The increments incx and incy are not used here (each is taken
c     as having a value of 1; input of any other value constitutes
c     an error). This is a pseudo-Linpack BLAS (Basic Linear Algebra
c     Subsystem) subroutine.
c
c     This subroutine is nearly identical to EQLIBU/dotpra.f, which has
c     the same function, but a different calling sequence.
c
c     This subroutine is called by:
c
c       Any
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       nmax   = dimension or pseudo-dimension of array
c       array1 = the first array
c       incx   = the x increment (must have a value of 1)
c       array2 = the second array
c       incy   = the y increment (must have a value of 1)
c
c     Output:
c
c       ddot   = the dot product
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
      integer incx,incy
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
c     Trap illegal values of incx and inc6.
c
      if (incx .ne. 1) then
        write (6,1000) incx
 1000   format(/' * Error - (EQLIBU/ddot) The argument incx has a',
     $  /7x,'value of ',i5,'. This version of ddot only allows this',
     $  /7x,'argument to have a value of 1. This is a programming',
     $  /7x,'error or a linking error.')
        stop
      endif
c
      if (incy .ne. 1) then
        write (6,1010) incy
 1010   format(/' * Error - (EQLIBU/ddot) The argument incy has a',
     $  /7x,'value of ',i5,'. This version of ddot only allows this',
     $  /7x,'argument to have a value of 1. This is a programming',
     $  /7x,'error or a linking error.')
        stop
      endif
c
c     Note the use of a local variable (dx) within the loop.
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
      ddot = dx
c
      end
