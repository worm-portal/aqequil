      integer function jiarmn(iarray,nmax)
c
c     This subroutine finds the minimum element of the integer array
c     iarray. Note that iarray could actually have more than one
c     dimension. For example, if the array is really "ia", which has
c     dimensions of i1 and i2, this subroutine could be used by calling
c     it in the following manner: call jiarmn(ia,(i1*i2)). Use this
c     subroutine with caution if nmax is not the product of the true
c     (declared) dimensions.
c
c     To find the corresponding maximum value, use EQLIBU/jiarmx.f
c     instead.
c
c     This subroutine is called by:
c
c       Any
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       iarray = the array
c       nmax   = dimension or pseudo-dimension of the array
c
c
c     Output:
c
c       jiarmn = the minimum element of the array (zero if nmax .le. 0)
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
      integer iarray(nmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ileft,ix
c
c-----------------------------------------------------------------------
c
      jiarmn = 0
      if (nmax .le. 0) go to 999
c
c     Note the use of a local variable (ix) within the loop.
c     Note also that the loop is unrolled.
c
      ix = iarray(1)
      ileft = (nmax/8)*8
c
      do i = 1,ileft,8
        ix = min(ix,iarray(i),iarray(i + 1),iarray(i + 2),iarray(i + 3),
     $       iarray(i + 4),iarray(i + 5),iarray(i + 6),iarray(i + 7))
      enddo
c
      do i = ileft + 1,nmax
        ix = min(ix,iarray(i))
      enddo
c
      jiarmn = ix
c
  999 continue
      end
