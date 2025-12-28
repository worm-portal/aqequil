      subroutine lamgex(acflgc,cgexj,jern1,jern2,jetmax,jgext,net,
     $ netmax,nstmax,xbarlg)
c
c     This subroutine computes the activity coefficients of generic ion
c     exchanger species. The exchanger phase is taken to be ideal
c     in the site-mixing sense. The activity (a) and activity coefficent
c     (lambda) are related to the mole fraction (x) as follows:
c
c       log a = N log x + N log lambda
c
c     where N is the site stoichiometric factor. This ensures that the
c     activity in the ideal case is equal to x**N, a treatment that
c     generally works better in hybrid Newton-Raphson iteration than
c     a treatment based on:
c
c       log a = log x + log lambda
c
c     whcih, in the ideal case, requires that:
c
c       log lambda = (N - 1) log x
c
c     The treatment adopted here is also preferable in the the activity
c     coefficient is reserved for treating actual non-ideality, rather
c     than as a correction used to relate two different definitions of
c     ideality.
c
c     This subroutine is called by:
c
c       EQLIB/ngcadv.f
c       EQ3NR/arrset.f
c       EQ6/exivar.f
c       EQ6/raff.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       cgexj  = array of site stoichiometry factors for site of generic
c                  ion exchange phases
c       jern1  = array giving the start of the range in the species
c                  list corresponding to species in the je-th site
c                  of the ne-th exchanger.
c       jern2  = array giving the end of the range in the species
c                  list corresponding to species in the je-th site
c                  of the ne-th exchanger.
c       jgext  = array giving the number of exchange sites in each of
c                  the ion exchange phases
c       xbarlg = array of mole fractions of species
c
c      Principal output:
c
c       acflgc = array of log activity coefficient values
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer jetmax,netmax,nstmax
c
      integer jern1(jetmax,netmax),jern2(jetmax,netmax),jgext(netmax)
c
      integer net
c
      real*8 acflgc(nstmax),cgexj(jetmax,netmax),xbarlg(nstmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer je,ne,nr1,nr2,ns
c
c-----------------------------------------------------------------------
c
c     Compute log lambda for generic exchanger species.
c
      do ne = 1,net
        do je = 1,jgext(ne)
          nr1 = jern1(je,ne)
          nr2 = jern2(je,ne)
          do ns = nr1,nr2
            acflgc(ns) = 0.
          enddo
        enddo
      enddo
c
      end
