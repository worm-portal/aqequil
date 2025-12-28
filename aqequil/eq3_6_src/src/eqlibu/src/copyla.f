      subroutine copyla(qarr1,qarr2,nmax)
c
c     This subroutine copies the first nmax elements of the logical
c     array qarr1 into the logical array qarr2. Normally, nmax would
c     be the dimension of both arrays. However, nmax could be less
c     than the true dimension. Also, the arrays could actually have
c     more than one dimension. For example, if qarr1 and qarr2
c     are really qa and qb, respectively, having common dimensions of
c     i1 and i2, this subroutine could be used by calling it in the
c     following manner: call copyla(qa,qb,(i1*i2)). Use this subroutine
c     with caution if the two arrays do not have identical dimensioning
c     or nmax is not the product of the actual (declared) dimensions.
c
c     This subroutine is called by:
c
c       Any
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       qarr1  = the first array
c       nmax   = dimension or pseudo-dimension of qarr1 and qarr2
c
c     Output:
c
c       qarr2  = the second array
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
      logical qarr1(nmax),qarr2(nmax)
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
        qarr2(i) = qarr1(i)
        qarr2(i + 1) = qarr1(i + 1)
        qarr2(i + 2) = qarr1(i + 2)
        qarr2(i + 3) = qarr1(i + 3)
        qarr2(i + 4) = qarr1(i + 4)
        qarr2(i + 5) = qarr1(i + 5)
        qarr2(i + 6) = qarr1(i + 6)
        qarr2(i + 7) = qarr1(i + 7)
      enddo
c
      do i = ileft + 1,nmax
        qarr2(i) = qarr1(i)
      enddo
c
      end
