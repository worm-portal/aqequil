      subroutine daxpy(nmax,scfact,array1,incx,array2,incy)
c
c     This subroutine multiplies the real*8 array "array1" by the
c     real*8 constant scfact and adds that to the real*8 array "array2".
c     The result is returned in "array2". The relevant equation is
c     (y) = a*(x) + (y). This subroutine is an adaptation of the 1979
c     Linpack subroutine of the same name. The size and order in the
c     original calling sequence has been preserved:
c     (n,da,dx,incx,dy,incy) == (nmax,scfact,array1,incx,array2,incy).
c     The increments incx and incy are not used here (each is taken
c     as having a value of 1; input of any other value constitutes
c     an error). This is a pseudo-Linpack BLAS (Basic Linear Algebra
c     Subsystem) subroutine.
c
c     This subroutine is called by:
c
c       EQLIBU/dgefa.f
c       EQLIBU/dgesl.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       nmax   = dimension or pseudo-dimension of array
c       scfact = the scale factor
c       array1 = the first array (x)
c       array2 = the second array (y)
c       incx   = the x increment (must have a value of 1)
c       incy   = the y increment (must have a value of 1)
c
c     Output:
c
c       array2 = the original array multiplied by the scale factor.
c
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
      real*8 scfact
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ileft
c
      real*8 scx
c
c-----------------------------------------------------------------------
c
c     Trap illegal values of incx and incy.
c
      if (incx .ne. 1) then
        write (6,1000) incx
 1000   format(/' * Error - (EQLIBU/daxpy) The argument incx has a',
     $  /7x,'value of ',i5,'. This version of daxpy only allows this',
     $  /7x,'argument to have a value of 1. This is a programming',
     $  /7x,'error or a linking error.')
        stop
      endif
c
      if (incy .ne. 1) then
        write (6,1010) incy
 1010   format(/' * Error - (EQLIBU/daxpy) The argument incy has a',
     $  /7x,'value of ',i5,'. This version of daxpy only allows this',
     $  /7x,'argument to have a value of 1. This is a programming',
     $  /7x,'error or a linking error.')
        stop
      endif
c
c     Note that the loop is unrolled.
c
      scx = scfact
      ileft = (nmax/8)*8
c
      do i = 1,ileft,8
        array2(i) = array2(i) + scx*array1(i)
        array2(i + 1) = array2(i + 1) + scx*array1(i + 1)
        array2(i + 2) = array2(i + 2) + scx*array1(i + 2)
        array2(i + 3) = array2(i + 3) + scx*array1(i + 3)
        array2(i + 4) = array2(i + 4) + scx*array1(i + 4)
        array2(i + 5) = array2(i + 5) + scx*array1(i + 5)
        array2(i + 6) = array2(i + 6) + scx*array1(i + 6)
        array2(i + 7) = array2(i + 7) + scx*array1(i + 7)
      enddo
c
      do i = ileft + 1,nmax
        array2(i) = array2(i) + scx*array1(i)
      enddo
c
      end
