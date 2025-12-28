      subroutine gcdrw(cdrs,cdrw,narn1,ndrs,ndrsmx,ndrsr,nst,nstmax)
c
c     This subroutine computes the cdrw array. Each element of this
c     array contains the reaction coefficient for liquid water in the
c     corresponding reaction.
c
c     Do not confuse this array with the cdrtw array computed by
c     EQLIB/gcdrtw.f.
c
c     This subroutine is called by:
c
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
      integer narn1,nst
c
      real*8 cdrs(ndrsmx),cdrw(nstmax)
c
c----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n,nn,nr1,nr2,ns
c
c----------------------------------------------------------------------
c
      do ns = 1,nst
        cdrw(ns) = 0.
        nr1 = ndrsr(1,ns)
        nr2 = ndrsr(2,ns)
        do n = nr1,nr2
          nn = ndrs(n)
          if (nn .eq. narn1) then
            cdrw(ns) = cdrs(n)
            go to 100
          endif
        enddo
  100   continue
      enddo
c
      end
