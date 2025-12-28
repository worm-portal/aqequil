      subroutine initiv(iarray,nmax,ivalue)
c
c     This subroutine initializes the integer array iarray to ivalue
c     over the first nmax positions. Normally, nmax would be the
c     dimension of the 1D iarray. However, nmax could be less than the
c     true dimension. Also, iarray could actually have more than one
c     dimension. For example, if the array is really ia, which has
c     dimensions of i1 and i2, this subroutine could be used by calling
c     it in the following manner: call initiv(ia,(i1*i2),ivalue).
c     Use this subroutine with caution if nmax is not the product of the
c     true (declared) dimensions.
c
c     To initialize an integer array to a zero value, use
c     EQLIBU/initiz.f instead.
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
c       nmax   = dimension or pseudo-dimension of array
c       ivalue = the initialization value
c
c     Output:
c
c       iarray = iarray, with the first nmax positions set to the
c                initialization value
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
      integer iarray(nmax),ivalue
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ileft,iv
c
c-----------------------------------------------------------------------
c
c     Note the assignment of the value to a local variable.
c     Note that the loop is unrolled.
c
      iv = ivalue
      ileft = (nmax/8)*8
c
      do i = 1,ileft,8
        iarray(i) = iv
        iarray(i + 1) = iv
        iarray(i + 2) = iv
        iarray(i + 3) = iv
        iarray(i + 4) = iv
        iarray(i + 5) = iv
        iarray(i + 6) = iv
        iarray(i + 7) = iv
      enddo
c
      do i = ileft + 1,nmax
        iarray(i) = iv
      enddo
c
      end
