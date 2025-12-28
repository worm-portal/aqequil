      subroutine pderiv(akmat0,demop0,fdpe0,nord,nordmx,npet,
     $ npetmx,nrd1mx)
c
c     This subroutine computes estimates of the derivatives of the
c     numbers of moles of the phases present in the equilibrium
c     system. These derivatives (demop0) are computed from the
c     corresponding finite differences (fdpe0). The relation is:
c     (demop0) = (akmat0)(fdpe0).
c
c     Compare with:
c
c       EQ6/aderiv.f
c       EQ6/bderiv.f
c       EQ6/rderiv.f
c       EQ6/sderiv.f
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
      integer nordmx,npetmx,nrd1mx
c
      integer nord,npet
c
      real*8 akmat0(nrd1mx,nrd1mx),demop0(nordmx,npetmx),
     $ fdpe0(nordmx,npetmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer k,n,npe
c
      real*8 dx
c
c-----------------------------------------------------------------------
c
      do npe = 1,npet
        do n = 1,nord
          dx = 0.
          do k = n,nord
            dx = dx + fdpe0(k,npe)*akmat0(n,k)
          enddo
          demop0(n,npe) = dx
        enddo
      enddo
c
      end
