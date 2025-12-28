      integer function iarmxn(array,nmax)
c
c     This subroutine finds the index of the element corresponding to
c     the max norm of the real*8 array "array". If there is more than
c     one such element, the index returned corresponds to the first one.
c     Note that array could actually have more than one dimension. For
c     example, if the array is really "a", which has dimensions of i1
c     and i2, this subroutine could be used by calling it in the
c     following manner: call iarmxn(a,(i1*i2)). Use this subroutine
c     with caution if nmax is not the product of the true (declared)
c     dimensions.
c
c     To just find the max norm itself, use EQLIBU/arrmxn.f instead.
c
c     This subroutine is nearly identical to EQLIBU/idamax.f, for which
c     the calling sequence arguments are reversed.
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
c       nmax   = dimension or pseudo-dimension of the array
c
c     Output:
c
c       iarmxn = index of the element corresponding to the max norm
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
      iarmxn = 0
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
      iarmxn = ix
c
  999 continue
      end
