      subroutine xderiv(akmat0,dxx0,fdxx0,nord,nordmx,nrd1mx)
c
c     This subroutine computes estimates of the derivatives (dxx0)
c     for some quantity (such as pH) from the corresponding finite
c     differences (fdxx0). Note that (dxx0) = (akmat0)(fdxx0).
c
c     Compare with:
c
c       EQ6/bderiv.f
c       EQ6/pderiv.f
c       EQ6/rderiv.f
c       EQ6/zderiv.f
c
c     This subroutine is called by:
c
c       EQ6/path.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c
c     Principal output:
c
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer nordmx,nrd1mx
c
      integer nord
c
      real*8 akmat0(nrd1mx,nrd1mx),dxx0(nordmx),fdxx0(nordmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer k,n
c
      real*8 dx
c
c-----------------------------------------------------------------------
c
      do n = 1,nord
        dx = 0.
        do k = n,nord
          dx = dx + fdxx0(k)*akmat0(n,k)
        enddo
        dxx0(n) = dx
      enddo
c
      end
