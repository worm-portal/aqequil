      real*8 function arrmxn(array,nmax)
c
c     This subroutine finds the max norm of the real*8 array "array".
c     Note that array could actually have more than one dimension. For
c     example, if the array is really "a", which has dimensions of i1
c     and i2, this subroutine could be used by calling it in the
c     following manner: call arrmxn(a,(i1*i2)). Use this subroutine
c     with caution if nmax is not the product of the true (declared)
c     dimensions.
c
c     To find the index of the first element corresponding to the max
c     norm, call EQLIBU/iarmxn.f instead. Then just use the index to
c     get the max norm if this is also desired.
c
c     This subroutine is called by:
c
c       Any
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       array  = the array
c       nmax   = the dimension or pseudo-dimension of the array
c
c     Output:
c
c       arrmxn = the max norm of the array
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
      real*8 array(nmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ileft
c
      real*8 aax,ax
c
c-----------------------------------------------------------------------
c
      arrmxn = 0.
      if (nmax .le. 0) go to 999
c
c     Note the use of a local variable (ax) within the loop.
c     Note also that the loop is unrolled.
c
      ax = 0.
      ileft = (nmax/8)*8
c
      do i = 1,ileft,8
        aax = abs(array(i))
        if (aax .gt. ax) ax = aax
        aax = abs(array(i + 1))
        if (aax .gt. ax) ax = aax
        aax = abs(array(i + 2))
        if (aax .gt. ax) ax = aax
        aax = abs(array(i + 3))
        if (aax .gt. ax) ax = aax
        aax = abs(array(i + 4))
        if (aax .gt. ax) ax = aax
        aax = abs(array(i + 5))
        if (aax .gt. ax) ax = aax
        aax = abs(array(i + 6))
        if (aax .gt. ax) ax = aax
        aax = abs(array(i + 7))
        if (aax .gt. ax) ax = aax
      enddo
c
      do i = ileft + 1,nmax
        aax = abs(array(i))
        if (aax .gt. ax) ax = aax
      enddo
c
      arrmxn = ax
c
  999 continue
      end
