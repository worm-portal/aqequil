      subroutine bderiv(akmat0,dafrc0,fdar0,jreac,nord,nordmx,nrct,
     $ nrctmx,nrd1mx)
c
c     This subroutine computes estimates of the derivatives of the
c     affinities of reactants (afrc1) from the corresponding finite
c     differences. Note that (dafrc0) = (akmat0)(fdar0).
c
c     Compare with:
c
c       EQ6/aderiv.f
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
      integer nordmx,nrctmx,nrd1mx
c
      integer nord,nrct
c
      integer jreac(nrctmx)
c
      real*8 akmat0(nrd1mx,nrd1mx),dafrc0(nordmx,nrctmx),
     $ fdar0(nordmx,nrctmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer k,n,nrc
c
      real*8 dx
c
c-----------------------------------------------------------------------
c
c     Recall the following jreac reactant status flag
c     conventions:
c
c       jreac =  0: set to react
c       jreac =  1: exhausted
c       jreac = -1: saturated, but the remaining reactant mass
c                     continues to react irreversibly
c       jreac =  2: saturated; the status of any remaining reactant
c                     mass is changed to that of a product phase
c
c     For jreac(nrc) = -1 or 2, the affinity is fixed at zero.
c     However, the actual calculated affinity values may be non-zero
c     owing to convergence tolerances. In EQ6/stepfd.f, the
c     corresponding finite differences are set to zero. In the present
c     subroutine, the corresponding derivatives are set to zero.
c
      do nrc = 1,nrct
        if (jreac(nrc).eq.-1 .or. jreac(nrc).eq.2) then
          do n = 1,nord
            dafrc0(n,nrc) = 0.
          enddo
        else
          do n = 1,nord
            dx = 0.
            do k = n,nord
              dx = dx + fdar0(k,nrc)*akmat0(n,k)
            enddo
            dafrc0(n,nrc) = dx
          enddo
        endif
      enddo
c
      end
