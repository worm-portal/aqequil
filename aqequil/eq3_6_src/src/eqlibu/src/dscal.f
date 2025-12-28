      subroutine dscal(nmax,scfact,array,incx)
c
c     This subroutine multiplies the real*8 array "array" by the real*8
c     constant scfact. The result is returned in "array". This
c     subroutine is an adaptation of the 1979 Linpack subroutine of the
c     same name. The size and order in the original calling sequence has
c     been preserved: (n,da,dx,incx) == (nmax,scfact,array,incx).
c     The increment incx is not used here (it is taken as having
c     a value of 1; input of any other value constitutes an error).
c     This is a pseudo-Linpack BLAS (Basic Linear Algebra Subsystem)
c     subroutine.
c
c     This subroutine is called by:
c
c       EQLIBU/dgefa.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       nmax   = dimension or pseudo-dimension of array
c       scfact = the scale factor
c       array  = the array
c       incx   = the increment (must have a value of 1)
c
c     Output:
c
c       array  = the original array multiplied by the scale factor.
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
      integer incx
c
      real*8 array(nmax)
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
c     Trap illegal value of incx.
c
      if (incx .ne. 1) then
        write (6,1000) incx
 1000   format(/' * Error - (EQLIBU/dscal) The argument incx has a',
     $  /7x,'value of ',i5,'. This version of dscal only allows this',
     $  /7x,'argument to have a value of 1. This is a programming',
     $  /7x,'error or a linking error.')
        stop
      endif
c
c     Note that the loop is unrolled.
c
      ileft = (nmax/8)*8
      scx = scfact
c
      do i = 1,ileft,8
        array(i) = scx*array(i)
        array(i + 1) = scx*array(i + 1)
        array(i + 2) = scx*array(i + 2)
        array(i + 3) = scx*array(i + 3)
        array(i + 4) = scx*array(i + 4)
        array(i + 5) = scx*array(i + 5)
        array(i + 6) = scx*array(i + 6)
        array(i + 7) = scx*array(i + 7)
      enddo
c
      do i = ileft + 1,nmax
        array(i) = scx*array(i)
      enddo
c
      end
