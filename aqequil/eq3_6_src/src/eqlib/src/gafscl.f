      subroutine gafscl(cdrsd,cscale,ndrsmx,ndrsrd,nst,nstmax)
c
c     This subroutine calculates the cscale array of affinity scaling
c     factors. Affinity scaling is used to help make decisions on
c     which phases are the best choices to precipitate when there
c     are multiple supersaturations.
c
c     This subroutine is called by:
c
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       cdrsd  = array of reaction coefficients ('d' set)
c
c     Principal output:
c
c       cscale = array of affinity scaling factors
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
      integer ndrsrd(2,nstmax)
      integer nst
c
      real*8 cdrsd(ndrsmx),cscale(nstmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n,ns,nr1,nr2
c
      real*8 cx
c
c-----------------------------------------------------------------------
c
      do ns = 1,nst
        cx = 0.
        nr1 = ndrsrd(1,ns)
        nr2 = ndrsrd(2,ns)
        do n = nr1 + 1,nr2
          cx = cx + abs(cdrsd(n))
        enddo
        if (cx .le. 0.) cx = 1.0
        cscale(ns) = cx
      enddo
c
      end
