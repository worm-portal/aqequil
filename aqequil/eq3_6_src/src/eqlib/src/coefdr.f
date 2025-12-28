      real*8 function coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns,nstmax)
c
c     This subroutine finds the coefficient of the nse-th species
c     in the reaction for the destruction of the ns-th species.
c
c     This subroutine is called by:
c
c       EQLIB/elim.f
c       EQLIB/ncmpex.f
c       EQLIB/switch.f
c       EQLIB/swtchb.f
c       EQLIB/swtchk.f
c       EQ3NR/arrset.f
c       EQ3NR/arrsim.f
c       EQ3NR/dawfix.f
c       EQ3NR/balcon.f
c       EQ3NR/betas.f
c       EQ3NR/eq3nr.f
c       EQ3NR/matrix.f
c       EQ3NR/scripx.f
c       EQ6/balcmz.f
c       EQ6/jgibbs.f
c       EQ6/matrxz.f
c       EQ6/mincsp.f
c       EQ6/scanlm.f
c       EQ6/tstrdx.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       cdrs   = array of reaction coefficients
c       ndrs   = array of indices of the species corresponding to the
c                  coefficients in the cdrs array
c       ndrsr  = array giving the range in the cdrs and ndrs arrays
c                  containing the reaction for a given species
c
c     Principal output:
c
c       coefdr = the reaction coefficient of the nse-th species in
c                  the reaction for the ns-th species
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ndrsmx,nstmax
c
      integer ndrs(ndrsmx),ndrsr(2,nstmax)
      integer nse,ns
c
      real*8 cdrs(ndrsmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n,n1,n2
c
c-----------------------------------------------------------------------
c
      coefdr = 0.
c
      n1 = ndrsr(1,ns)
      n2 = ndrsr(2,ns)
      do n = n1,n2
        if (nse .eq. ndrs(n)) then
          coefdr = cdrs(n)
          go to 999
        endif
      enddo
c
  999 continue
      end
