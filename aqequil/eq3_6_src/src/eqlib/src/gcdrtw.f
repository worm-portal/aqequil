      subroutine gcdrtw(cdrs,cdrtw,narn1,narn2,ndrs,ndrsmx,ndrsr,
     $ nelect,no2gaq,nst,nstmax)
c
c     This subroutine computes the cdrtw array. Each element of this
c     array contains the sum of the reaction coefficients of the
c     aqueous solute species in the corresponding reaction. The sum
c     represents the number of times that the number of moles of
c     solvent water is implicitly represented in the reaction via
c     solute molalities. Note that generic ion exchanger species do
c     not count in these calculations, even though molalities can
c     be calculated for such species. That is because their
c     thermodynamic activities are defined by mole fractions.
c
c     Do not confuse this array with the cdrw array computed by
c     EQLIB/gcdrw.f. The present array is really only used by
c     EQ6. It presently exists in EQ3NR only for the sake of
c     consistency.
c
c     This subroutine is called by:
c
c       EQLIB/absswa.f
c       EQ3NR/eq3nr.f
c       EQ6/path.f
c       EQ6/absswb.f
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
c     Calling sequence variable declarations.
c
      integer ndrsmx,nstmax
c
      integer ndrs(ndrsmx),ndrsr(2,nstmax)
c
      integer narn1,narn2,nelect,no2gaq,nst
c
      real*8 cdrs(ndrsmx),cdrtw(nstmax)
c
c----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n,nn,nr1,nr2,ns
c
      real*8 cx
c
c----------------------------------------------------------------------
c
      do ns = 1,nst
        nr1 = ndrsr(1,ns)
        nr2 = ndrsr(2,ns)
        cx = 0.
        do n = nr1,nr2
          nn = ndrs(n)
          if (nn.gt.narn1 .and. nn.le.narn2) then
            if (nn.ne.no2gaq .and. nn.ne.nelect) cx = cx + cdrs(n)
          endif
        enddo
        cdrtw(ns) = cx
      enddo
c
      end
