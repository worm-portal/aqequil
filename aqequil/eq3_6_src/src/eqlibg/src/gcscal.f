      subroutine gcscal(acflgc,delacf,narn1,narn2,nref,nstmax,zchar)
c
c     This subroutine rescales the activity coefficients of aqeuous
c     ionic species to make them consistent with a given pH scale,
c     which has been used to define the correction parameter "delacf".
c
c     This subroutine is called by:
c
c       EQLIBG/gcoeff.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       acflgc = calculated activity coefficient array  (unscaled)
c       delacf = log gamma of the reference ion (new scale) -
c                  log gamma of the reference ion (old scale)
c       nref   = the species index of the reference ion
c       nstmax = the maximum number of species
c       zchar  = electrical charge array
c
c     Output:
c
c       acflgc = calculated activity coefficient array (scaled)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nstmax
c
      integer narn1,narn2,nref
c
      real*8 acflgc(nstmax),zchar(nstmax)
      real*8 delacf
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ileft,ns,nval
c
      real*8 zref
c
c-----------------------------------------------------------------------
c
      zref = zchar(nref)
c
c     Note that the loop can run from narn1 to narn2 without
c     causing a problem in the case of water because water has no
c     electrical charge. Note also that the loop is unrolled.
c
      nval = narn2 - narn1 + 1
      ileft = (nval/8)*8 + narn1 - 1
c
      do ns = narn1,ileft,8
        acflgc(ns) = acflgc(ns) + (zchar(ns)/zref)*delacf
        acflgc(ns + 1) = acflgc(ns + 1) + (zchar(ns + 1)/zref)*delacf
        acflgc(ns + 2) = acflgc(ns + 2) + (zchar(ns + 2)/zref)*delacf
        acflgc(ns + 3) = acflgc(ns + 3) + (zchar(ns + 3)/zref)*delacf
        acflgc(ns + 4) = acflgc(ns + 4) + (zchar(ns + 4)/zref)*delacf
        acflgc(ns + 5) = acflgc(ns + 5) + (zchar(ns + 5)/zref)*delacf
        acflgc(ns + 6) = acflgc(ns + 6) + (zchar(ns + 6)/zref)*delacf
        acflgc(ns + 7) = acflgc(ns + 7) + (zchar(ns + 7)/zref)*delacf
      enddo
c
      do ns = ileft + 1,narn2
        acflgc(ns) = acflgc(ns) + (zchar(ns)/zref)*delacf
      enddo
c
      end
