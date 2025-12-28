      subroutine gsgsm(conc,dpslm,na,natmax,nsltmx,nstmax,nsxi,nsxmax,
     $ nsxx,pslm,slsum,slsump,uspec)
c
c     This subroutine computes the following first-order sums used
c     in Pitzer's equations:
c
c       SUM(j) S-lambda(ij)*m(j)
c       SUM(j) S-lambda'(ij)*m(j)
c
c     This subroutine is called by:
c
c       EQLIBG/gcoeff.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       conc   = array of species concentrations
c       na     = aqueous species index (i refers to ns, where ns is
c                  the corresponding species index)
c       nsxi   = range pointer array into the nsxx array:
c                  nsxi(1,na) and nsxi(2,na) are the first and last
c                  values of the second subscript (nsx) of the nsxx
c                  array for entries pertaining to the na-th
c                  aqueous species (ns-th species)
c       nsxx   = pointer array:
c                  nsxx(1,nsx) = the species index of the second
c                  species in the nsl-th pair, where
c                  nsxx(2,nsx) = nsl
c       pslm   = array of S-lambda functions
c       dpslm  = array of ionic strength derivatives of the
c                  S-lambda functions
c
c     Principal input:
c
c       slsum  = the sum: SUM(j) S-lambda(ij)*m(j)
c       slsump = the sum: SUM(j) S-lambda'(ij)*m(j)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer natmax,nsltmx,nstmax,nsxmax
c
      integer nsxi(2,natmax),nsxx(2,nsxmax)
      integer na
c
      character(len=48) uspec(nstmax)
c
      real(8) conc(nstmax),dpslm(2,nsltmx),pslm(nsltmx)
      real(8) slsum,slsump
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ixf,ixl,nsl,nsx,ns2
c
c-----------------------------------------------------------------------
c
      slsum = 0.
      slsump = 0.
      ixf = nsxi(1,na)
      ixl = nsxi(2,na)
c
      do nsx = ixf,ixl
        ns2 = nsxx(1,nsx)
        nsl = nsxx(2,nsx)
        slsum = slsum + pslm(nsl)*conc(ns2)
        slsump = slsump + dpslm(1,nsl)*conc(ns2)
      enddo
c
      end
