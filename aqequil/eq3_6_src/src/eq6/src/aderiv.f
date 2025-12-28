      subroutine aderiv(akmat0,daffp0,fdaf0,nord,nordmx,npt,
     $ nptmax,nrd1mx)
c
c     This subroutine computes estimates of the derivatives of the phase
c     affinities (affp) from the corresponding finite differences.
c     Note that (daffp0) = (akmat0)(fdaf0).
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
      integer nordmx,nptmax,nrd1mx
c
      integer nord,npt
c
      real*8 akmat0(nrd1mx,nrd1mx),daffp0(nordmx,nptmax),
     $ fdaf0(nordmx,nptmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer k,n,np
c
      real*8 dx
c
c-----------------------------------------------------------------------
c
      do np = 1,npt
        do n = 1,nord
          dx = 0.
          do k = n,nord
            dx = dx + fdaf0(k,np)*akmat0(n,k)
          enddo
          daffp0(n,np) = dx
        enddo
      enddo
c
      end
