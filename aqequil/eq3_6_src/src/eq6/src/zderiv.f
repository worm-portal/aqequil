      subroutine zderiv(akmat0,dzvc0,fdzv0,kdim,kmax,nord,nrd1mx)
c
c     This subroutine computes estimates of the derivatives of the
c     master algebraic variables (z vector) from the corresponding
c     finite differences. Note that (dzvc0) = (akmat0)(fdzv0).
c
c     Compare with:
c       EQ6/aderiv.f
c       EQ6/bderiv.f
c       EQ6/pderiv.f
c       EQ6/rderiv.f
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
      integer kmax,nrd1mx
c
      integer kdim,nord
c
      real*8 akmat0(nrd1mx,nrd1mx),dzvc0(nrd1mx,kmax),fdzv0(nrd1mx,kmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer k,kcol,n
c
      real*8 dx
c
c-----------------------------------------------------------------------
c
      do kcol = 1,kdim
        do n = 1,nord
          dx = 0.
          do k = n,nord
            dx = dx + fdzv0(k,kcol)*akmat0(n,k)
          enddo
          dzvc0(n,kcol) = dx
        enddo
      enddo
c
      end
