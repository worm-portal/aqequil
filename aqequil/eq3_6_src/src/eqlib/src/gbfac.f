      subroutine gbfac(beta,bfac,efac,iindx1,kbt,kmax,nbt,nbtmax,nfac)
c
c     This subroutine calculates the bfac array, which is used in making
c     continued fraction corrections. It resolves conflicts when the
c     same aqueous species dominates more than one mass balance (the
c     species dominating a given mass balance is determined by
c     EQLIB/fdomsp). The continued fraction algorithm can be applied
c     to the master species associated with only one mass balance in
c     such a set, otherwise, oscillatory behavior will occur. In each
c     set of mass balances with a common dominant species, this
c     subroutine finds the mass balance with the greatest residual and
c     completes the calculation of its bfac factor by doing the
c     appropriate exponentiation. It sets bfac to unity for the others
c     in the set.
c
c     This subroutine is called by:
c
c       EQ3NR/arrset.f
c       EQ6/optmzr.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       nfac   = array of indices of dominant aqueous species
c       beta   = array of normalized Newton-Raphson residual functions
c       efac   = array of corresponding reciprocal stoichiometric
c                  weights
c       kbt    = the number of active basis species
c       nbt    = the number of active basis species
c
c     Principal output:
c
c       bfac   = array (in terms of matrix indexing) of the quantity
c                  used in making a continued fraction correction
c                  (e.g., conc (new) = conc (old) / bfac ); in general,
c                  bfac = ( beta + 1. )**efac, but bfac may be set to
c                  unity for a species in resolving a conflict
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer kmax,nbtmax
c
      integer iindx1(kmax),nfac(nbtmax)
      integer kbt,nbt
c
      real*8 beta(kmax),bfac(nbtmax),efac(nbtmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer krow,nb,nb1,nb2
c
      real*8 bpx1
c
c-----------------------------------------------------------------------
c
c     Calculate the bfac correction factors: bfac = (beta + 1)**efac.
c
      do krow = 1,kbt
        nb = iindx1(krow)
        bpx1 = beta(krow) + 1.
        if (bpx1 .le. 0.) then
c
c         Protect against a singularity in case (beta + 1) is not
c         a positive number. The value set in such a case allows a
c         generous enough one-step increase in the value of the
c         corrected concentration, number of moles, etc. (20 orders
c         of magnitude). A more restrictive limit may be applied
c         later when the actual correction is made.
c
          bfac(nb) = 1.e-20
        else
c
c         Calculate bfac from the usual formula.
c
          bfac(nb) = bpx1**efac(nb)
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Eliminate conflicts.
c
      do nb1 = 1,nbt - 1
        if (nfac(nb1) .gt. 0) then
          do nb2 = nb1 + 1,nbt
            if (nfac(nb2) .eq. nfac(nb1)) then
              if (bfac(nb1) .gt. bfac(nb2)) then
                nfac(nb2) = 0
                bfac(nb2) = 1.
              else
                nfac(nb1) = 0
                bfac(nb1) = 1.
                go to 15
              endif
            endif
          enddo
        endif
   15   continue
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
