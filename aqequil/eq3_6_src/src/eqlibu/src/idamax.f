      integer function idamax(nmax,array,incx)
c
c     This subroutine finds the index of the element corresponding to
c     the max norm of the real*8 array "array". If there is more than
c     one such element, the index returned corresponds to the first one.
c     This subroutine is an adaptation of the 1979 Linpack subroutine of
c     the same name. The size and order in the original calling sequence
c     has been preserved: (n,dx,incx) == (nmax,array,incx).
c     The increment incx is not used here (it is taken as having
c     a value of 1; input of any other value constitutes an error).
c     This is a pseudo-Linpack BLAS (Basic Linear Algebra Subsystem)
c     subroutine.
c
c     This subroutine is nearly identical to EQLIBU/iarmnx.f, for which
c     the calling sequence arguments are reversed.
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
c       array  = the array
c       incx   = the increment (must have a value of 1)
c
c     Output:
c
c       idamax = index of the element corresponding to the max norm
c                  (zero if nmax .le. 0)
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
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ileft,ix,j
c
      real*8 aax,arrmxn
c
c-----------------------------------------------------------------------
c
c     Trap illegal value of incx.
c
      if (incx .ne. 1) then
        write (6,1000) incx
 1000   format(/' * Error - (EQLIBU/idamax) The argument incx has a',
     $  /7x,'value of ',i5,'. This version of idamax only allows this',
     $  /7x,'argument to have a value of 1. This is a programming',
     $  /7x,'error or a linking error.')
        stop
      endif
c
      idamax = 0
      arrmxn = 0.
      if (nmax .le. 0) go to 999
c
c     Note the use of a local variable (ix) within the loop.
c     Note also that the loop is unrolled.
c
      ix = 0
      ileft = (nmax/8)*8
c
      do i = 1,ileft,8
        aax = abs(array(i))
        if (aax .gt. arrmxn) then
          ix = i
          arrmxn = aax
        endif
        j = i + 1
        aax = abs(array(j))
        if (aax .gt. arrmxn) then
          ix = j
          arrmxn = aax
        endif
        j = i + 2
        aax = abs(array(j))
        if (aax .gt. arrmxn) then
          ix = j
          arrmxn = aax
        endif
        j = i + 3
        aax = abs(array(j))
        if (aax .gt. arrmxn) then
          ix = j
          arrmxn = aax
        endif
        j = i + 4
        aax = abs(array(j))
        if (aax .gt. arrmxn) then
          ix = j
          arrmxn = aax
        endif
        j = i + 5
        aax = abs(array(j))
        if (aax .gt. arrmxn) then
          ix = j
          arrmxn = aax
        endif
        j = i + 6
        aax = abs(array(j))
        if (aax .gt. arrmxn) then
          ix = j
          arrmxn = aax
        endif
        j = i + 7
        aax = abs(array(j))
        if (aax .gt. arrmxn) then
          ix = j
          arrmxn = aax
        endif
      enddo
c
      do i = ileft + 1,nmax
        aax = abs(array(i))
        if (aax .gt. arrmxn) then
          ix = i
          arrmxn = aax
        endif
      enddo
c
      idamax = ix
c
  999 continue
      end
