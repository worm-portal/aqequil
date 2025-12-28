      subroutine cdrssd(adhfs,adhfsd,advfs,advfsd,axhfs,axhfsd,
     $ axlks,axlksd,axvfs,axvfsd,cdrs,cdrsd,ipch,ipchmx,ipcv,ipcvmx,
     $ narxmx,nbasp,nbaspd,nbtmax,ndrs,ndrsd,ndrsmx,ndrsr,ndrsrd,
     $ nstmax,ntprmx)
c
c     This subroutine copies the reactions and reaction properties as
c     they are currently written into the 'd' set.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       adhfs  = standard array of coefficients for computing
c                  pressure derivatives of enthalpy functions as a
c                  function of temperature
c       advfs  = standard array of coefficients for computing
c                  pressure derivatives of volume functions as a
c                  function of temperature
c       axhfs  = standard array of coefficients for computing
c                  enthalpy functions as a function of temperature
c       axlks  = standard array of coefficients for computing
c                  equilibrium constants as a function of temperature
c       axvfs  = standard array of coefficients for computing
c                  volume functions as a function of temperature
c       cdrs   = standard array of reaction coefficients
c       narxmx = number of coefficient elements of axlks per species
c                  per temperature range
c       nbtmax = maximum number of basis species
c       ndrs   = standard  array of indices of species corresponding to
c                  reaction coefficients
c       ndrsr  = standard pointer array for indices of species appearing
c                  in reactions
c       nstmax = maximum number of species
c       ntprmx = maximum number of temperature ranges
c
c     Principal output:
c
c       adhfsd = 'data file' array of coefficients for computing
c                  pressure derivatives of enthalpy functions as a
c                  function of temperature
c       advfsd = 'data file' array of coefficients for computing
c                  pressure derivatives of volume functions as a
c                  function of temperature
c       axhfsd = 'data file' array of coefficients for computing
c                  enthalpy functions as a function of temperature
c       axlksd = 'data file' array of coefficients for computing
c                  equilibrium constants as a function of temperature
c       axvfsd = 'data file' array of coefficients for computing
c                  volume functions as a function of temperature
c       cdrsd  = 'data file' array of reaction coefficients
c       ndrsd  = 'data file'  array of indices of species corresponding
c                  to reaction coefficients
c       ndrsrd = 'data file' pointer array for indices of species
c                  appearing in reactions
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ipchmx,ipcvmx,narxmx,nbtmax,ndrsmx,nstmax,ntprmx
c
      integer nbasp(nbtmax),nbaspd(nbtmax),ndrs(ndrsmx),ndrsd(ndrsmx),
     $ ndrsr(2,nstmax),ndrsrd(2,nstmax)
      integer ipch,ipcv
c
      real*8 adhfs(narxmx,ntprmx,ipchmx,nstmax),
     $ adhfsd(narxmx,ntprmx,ipchmx,nstmax),
     $ advfs(narxmx,ntprmx,ipcvmx,nstmax),
     $ advfsd(narxmx,ntprmx,ipcvmx,nstmax),
     $ axhfs(narxmx,ntprmx,nstmax),axhfsd(narxmx,ntprmx,nstmax),
     $ axlks(narxmx,ntprmx,nstmax),axlksd(narxmx,ntprmx,nstmax),
     $ axvfs(narxmx,ntprmx,nstmax),axvfsd(narxmx,ntprmx,nstmax),
     $ cdrs(ndrsmx),cdrsd(ndrsmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jhmax,jvmax,nmax
c
c-----------------------------------------------------------------------
c
      nmax = narxmx*ntprmx*nstmax
      jhmax = narxmx*ntprmx*ipchmx*nstmax
      jvmax = narxmx*ntprmx*ipcvmx*nstmax
c
      call copyaa(axlks,axlksd,nmax)
c
      if (ipch .ge. 0) then
        call copyaa(axhfs,axhfsd,nmax)
        if (ipch .ge. 1) call copyaa(adhfs,adhfsd,jhmax)
      endif
c
      if (ipcv .ge. 0) then
        call copyaa(axvfs,axvfsd,nmax)
        if (ipcv .ge. 1) call copyaa(advfs,advfsd,jvmax)
      endif
c
      nmax = 2*nstmax
      call copyia(ndrsr,ndrsrd,nmax)
c
      call copyaa(cdrs,cdrsd,ndrsmx)
      call copyia(ndrs,ndrsd,ndrsmx)
c
      call copyia(nbasp,nbaspd,nbtmax)
c
      end
