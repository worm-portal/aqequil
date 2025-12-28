      subroutine initlv(qarray,nmax,qvalue)
c
c     This subroutine initializes the logical array qarray to qvalue
c     over the first nmax positions. Normally, nmax would be the
c     dimension of the 1D qarray. However, nmax could be less than the
c     true dimension. Also, qarray could actually have more than one
c     dimension. For example, if the array is really qa, having
c     dimensions of i1 and i2, this subroutine could be used by calling
c     it in the following manner: call initlv(qa,(i1*i2),qvalue).
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
c       qarray = the array
c       nmax   = dimension or pseudo-dimension of array
c       qvalue = the initialization value (.true. or .false.)
c
c     Output:
c
c       qarray = qarray, with the first nmax positions set to the
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
      logical qarray(nmax),qvalue
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ileft
c
      logical qv
c
c-----------------------------------------------------------------------
c
c     Note the assignment of the value to a local variable.
c     Note that the loop is unrolled.
c
      qv = qvalue
      ileft = (nmax/8)*8
c
      do i = 1,ileft,8
        qarray(i) = qv
        qarray(i + 1) = qv
        qarray(i + 2) = qv
        qarray(i + 3) = qv
        qarray(i + 4) = qv
        qarray(i + 5) = qv
        qarray(i + 6) = qv
        qarray(i + 7) = qv
      enddo
c
      do i = ileft + 1,nmax
        qarray(i) = qv
      enddo
c
      end
