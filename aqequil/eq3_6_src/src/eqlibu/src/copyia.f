      subroutine copyia(iarr1,iarr2,nmax)
c
c     This subroutine copies the first nmax elements of the integer
c     array iarr1 into the integer array iarr2. Normally, nmax would
c     be the dimension of both arrays. However, nmax could be less
c     than the true dimension. Also, the arrays could actually have
c     more than one dimension. For example, if iarr1 and iarr2
c     are really ia and ib, respectively, having common dimensions of
c     i1 and i2, this subroutine could be used by calling it in the
c     following manner: call copyia(ia,ib,(i1*i2)). Use this subroutine
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
c       iarr1  = the first array
c       nmax   = dimension or pseudo-dimension of iarr1 and iarr2
c
c     Output:
c
c       iarr2  = the second array
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
      integer iarr1(nmax),iarr2(nmax)
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
        iarr2(i) = iarr1(i)
        iarr2(i + 1) = iarr1(i + 1)
        iarr2(i + 2) = iarr1(i + 2)
        iarr2(i + 3) = iarr1(i + 3)
        iarr2(i + 4) = iarr1(i + 4)
        iarr2(i + 5) = iarr1(i + 5)
        iarr2(i + 6) = iarr1(i + 6)
        iarr2(i + 7) = iarr1(i + 7)
      enddo
c
      do i = ileft + 1,nmax
        iarr2(i) = iarr1(i)
      enddo
c
      end
