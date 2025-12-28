      subroutine initav(array,nmax,value)
c
c     This subroutine initializes the real*8 array "array" to the value
c     of "value" over the first nmax positions. Normally, nmax would
c     be the dimension of the 1D "array". However, nmax could be less
c     than the true dimension. Also, array could actually have more
c     than one dimension. For example, if the array is really "a",
c     having dimensions of i1 and i2, this subroutine could be used by
c     calling it in the following manner: call initar(a,(i1*i2),value).
c     Use this subroutine with caution if the two arrays do not have
c     identical dimensioning or nmax is not the product of the true
c     (declared) dimensions.
c
c     To initialize a real*8 array to a zero value, use EQLIBU/initaz.f
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
c       array  = the array
c       nmax   = dimension or pseudo-dimension of array
c       value  = the initialization value
c
c     Output:
c
c       array  = array, with the first nmax positions set to the
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
      real*8 array(nmax),value
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ileft
c
      real*8 vx
c
c-----------------------------------------------------------------------
c
c     Note the assignment of the value to a local variable.
c     Note that the loop is unrolled.
c
      vx = value
      ileft = (nmax/8)*8
c
      do i = 1,ileft,8
        array(i) = vx
        array(i + 1) = vx
        array(i + 2) = vx
        array(i + 3) = vx
        array(i + 4) = vx
        array(i + 5) = vx
        array(i + 6) = vx
        array(i + 7) = vx
      enddo
c
      do i = ileft + 1,nmax
        array(i) = vx
      enddo
c
      end
