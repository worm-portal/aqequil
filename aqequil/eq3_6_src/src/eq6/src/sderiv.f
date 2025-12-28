      subroutine sderiv(akmat0,demos0,fdse0,nord,nordmx,nrd1mx,
     $ nset,nsetmx)
c
c     This subroutine computes estimates of the derivatives of the
c     numbers of moles of selected species in the equilibrium system.
c     These species include H2O(l) for the aqueous solution and all
c     species of the non-aqueous phases. These derivatives (demos0)
c     are computed from the corresponding finite differences (fdse0).
c     The relation is: (demos0) = (akmat0)(fdse0).
c
c     Compare with:
c
c       EQ6/aderiv.f
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
      integer nordmx,nrd1mx,nsetmx
c
      integer nord,nset
c
      real*8 akmat0(nrd1mx,nrd1mx),demos0(nordmx,nsetmx),
     $ fdse0(nordmx,nsetmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer k,n,nse
c
      real*8 dx
c
c-----------------------------------------------------------------------
c
      do nse = 1,nset
        do n = 1,nord
          dx = 0.
          do k = n,nord
            dx = dx + fdse0(k,nse)*akmat0(n,k)
          enddo
          demos0(n,nse) = dx
        enddo
      enddo
c
      end
