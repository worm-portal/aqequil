      subroutine calk(alkc,conc,nstmax,ntfx,ntfxmx,ntfxt,tfx)
c
c     This subroutine calculates the alkalinity (eq/kg H2O). A sorted
c     summation is not done here, because relatively few species
c     contribute to alkalinity. Also, the structure of the titration
c     factor arrays (tfx, ntfx) is not amenable to an efficient
c     calculation of that sort.
c
c     This subroutine is called by:
c
c       EQLIB/betas.f
c       EQLIB/prtalk.f
c       EQ3NR/scripx.f
c       EQ6/scripz.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       conc   = array of species concentrations
c       ntfx   = array of species indices corresponding to the
c                  alkalinity factors in the tfx array
c       ntfxt  = the number of species contributing to alkalinity
c       tfx    = array of alkalinity factors
c
c     Principal output:
c
c       alkc   = the calculated alkalinity (eq/kg H20)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nstmax,ntfxmx
c
      integer ntfx(ntfxmx)
      integer ntfxt
c
      real*8 conc(nstmax),tfx(ntfxmx)
      real*8 alkc
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ileft,n
c
      real*8 ax
c
c-----------------------------------------------------------------------
c
c     Note that the loop is unrolled.
c
      ax = 0.
      ileft = (ntfxt/8)*8
c
      do n = 1,ileft,8
        ax = ax + tfx(n)*conc(ntfx(n))
     $          + tfx(n + 1)*conc(ntfx(n + 1))
     $          + tfx(n + 2)*conc(ntfx(n + 2))
     $          + tfx(n + 3)*conc(ntfx(n + 3))
     $          + tfx(n + 4)*conc(ntfx(n + 4))
     $          + tfx(n + 5)*conc(ntfx(n + 5))
     $          + tfx(n + 6)*conc(ntfx(n + 6))
     $          + tfx(n + 7)*conc(ntfx(n + 7))
      enddo
c
      do n = ileft + 1,ntfxt
        ax = ax + tfx(n)*conc(ntfx(n))
      enddo
c
      alkc = ax
c
      end
