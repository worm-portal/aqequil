      subroutine gslam(dgpit,dpslm,gpit,ipbtmx,nalpha,napmax,nslt,
     $ nsltmx,pslamn,pslm)
c
c     This subroutine computes the S-lambda coefficients and their
c     ionic strength derivatives. These coefficients are used
c     in Pitzer's equations.
c
c     This subroutine is called by:
c
c       EQLIBG/gcoeff.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       gpit   = array of values for the g(x) function
c       dgpit  = array of values of the ionic strength derivatives
c                  of g(x)
c       nslt   = number of S-lambda coeffcient parameter sets
c       pslamn = array of S-lambda(n) coefficient parameters
c       nalpha = array that gives the index in the palpha array
c                  corresponding to a given set of S-lambda coefficient
c                  parameters.
c
c     Principal output:
c
c       pslm   = array of the S-lambda functions
c       dpslm  = array of ionic strength derivatives of the
c                  S-lambda functions
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ipbtmx,napmax,nsltmx
c
      integer nalpha(nsltmx)
      integer nslt
c
      real(8) dgpit(2,ipbtmx,napmax),dpslm(2,nsltmx),
     $ gpit(ipbtmx,napmax),pslamn(0:ipbtmx,nsltmx),pslm(nsltmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,k,nap,nsl
c
      real(8) dpx,px
c
c----------------------------------------------------------------------
c
      do nsl = 1,nslt
c
        nap = nalpha(nsl)
c
        px = pslamn(0,nsl)
        do i = 1,ipbtmx
          px = px + gpit(i,nap)*pslamn(i,nsl)
        enddo
        pslm(nsl) = px
c
        do k = 1,2
          dpx = 0.
          do i = 1,ipbtmx
            dpx = dpx + dgpit(k,i,nap)*pslamn(i,nsl)
          enddo
          dpslm(k,nsl) = dpx
        enddo
      enddo
c
      end
