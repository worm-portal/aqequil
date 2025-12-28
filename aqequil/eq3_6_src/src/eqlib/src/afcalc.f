      subroutine afcalc(actlg,af,afcnst,cdrs,jflag,jsflag,ndrs,ndrsmx,
     $ ndrsr,ns,nstmax,si,xlks)
c
c     This subroutine computes the affinity function:
c
c       A = 2.303 RT log Q/K
c
c     and the saturation index:
c
c       SI = RT log Q/K
c
c     for the ns-th reaction. Note that this affinity is the affinity
c     for the forward direction (A(+)), hence it has the same sign as
c     the saturation index (SI).
c
c     EQLIB/gafsir.f makes the same calculations for all reactions.
c
c     This subroutine is called by:
c
c       EQLIB/betas.f
c       EQLIB/hpsat.f
c       EQ6/raff.f
c       EQ6/satchk.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       actlg  = array of log activities of species
c       afcnst = the factor 2.303 RT
c       cdrs   = array of reaction coefficients
c       jflag  = flag array denoting whether species are in the active
c                  basis set (jflag is not 30) or out of it (jflag = 30)
c       ndrs   = array of indices of species whose reaction coefficients
c                  are in the cdrs array
c       ndrsr  = pointer array for the cdrs/ndrs arrays, denoting the
c                  range of entries corresponding to a given reaction
c       ns     = index of the species for whose associated reaction
c                  the affinity and saturation index are to be
c                  calculated
c       xlks   = array of equilibrium constants
c
c     Principal output:
c
c       af     = affinity
c       si     = saturation index (SI)
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
      integer jflag(nstmax),jsflag(nstmax),ndrs(ndrsmx),ndrsr(2,nstmax)
      integer ns
c
      real*8 actlg(nstmax),cdrs(ndrsmx),xlks(nstmax)
      real*8 af,afcnst,si
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n,nr1,nr2,nse,nt
c
c-----------------------------------------------------------------------
c
      if (jsflag(ns) .eq. 2) then
c
c       Species is suppressed.
c
        si = -9999999.
        af = -9999999.
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (jflag(ns) .eq. 30) then
c
c       Species is active, but is required to satisfy equilibrium.
c
        si = 0.
        af = 0.
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      nr1 = ndrsr(1,ns)
      nr2 = ndrsr(2,ns)
      nt = nr2 - nr1 + 1
      if (nt .lt. 2) then
c
c       Species has no reaction. Take it to be in equilibrium with
c       itself.
c
        si = 0.
        af = 0.
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (xlks(ns) .ge. 9999999.) then
c
c       The equilibrium constant is unknown (set to pseudo-infinity).
c       Set the SI and affinity to negative pseudo-infinity.
c
        si = -9999999.
        af = -9999999.
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the SI and affinity in the normal way.
c
      si = -xlks(ns)
      do n = nr1,nr2
        nse = ndrs(n)
        if (nse .eq. 0) go to 200
        si = si + cdrs(n)*actlg(nse)
      enddo
      af = afcnst*si
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  200 continue
c
c     The species has no valid reaction. It is a detached auxiliary
c     basis species. Take it to be in equilibrium with itself.
c
      si = 0.
      af = 0.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
