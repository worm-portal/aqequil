      subroutine initii(iarray,nmax)
c
c     This subroutine initializes each element of the integer array
c     iarray its index value over the first nmax positions. Normally,
c     nmax would be the dimension of the 1D iarray. However, nmax could
c     be less than the true dimension. Also, iarray could actually have
c     more than one dimension. For example, if the array is really ib,
c     which has dimensions of i1 and i2, this subroutine could be used
c     by calling it in the following manner: call initii(ib,(i1*i2)).
c     Use this subroutine with caution if nmax is not the product of the
c     true (declared) dimensions.
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
c       nmax   = dimension or pseudo-dimension of iarray
c
c     Output:
c
c       iarray = iarray, with the first nmax positions set to the
c                  corresponding index values
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
      integer i,ileft
c
c-----------------------------------------------------------------------
c
c     Note that the loop is unrolled.
c
      ileft = (nmax/8)*8
c
      do i = 1,ileft,8
        iarray(i) = i
        iarray(i + 1) = i + 1
        iarray(i + 2) = i + 2
        iarray(i + 3) = i + 3
        iarray(i + 4) = i + 4
        iarray(i + 5) = i + 5
        iarray(i + 6) = i + 6
        iarray(i + 7) = i + 7
      enddo
c
      do i = ileft + 1,nmax
        iarray(i) = i
      enddo
c
      end
