      subroutine cdrscx(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,axlks,
     $ axlksx,axvfs,axvfsx,cdrs,cdrsx,ipch,ipchmx,ipcv,ipcvmx,
     $ narxmx,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,nstmax,ntprmx)
c
c     This subroutine copies the scratch arrays for reactions and
c     reaction properties into the standard arrays.
c
c     This subroutine is called by:
c
c       EQLIB/elim.f
c       EQLIB/mdrgex.f
c       EQLIB/switch.f
c       EQLIB/swtchb.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       adhfsx = scratch array of coefficients for computing
c                  pressure derivatives of enthalpy functions as a
c                  function of temperature
c       advfsx = scratch array of coefficients for computing
c                  pressure derivatives of volume functions as a
c                  function of temperature
c       axhfsx = scratch array of coefficients for computing
c                  enthalpy functions as a function of temperature
c       axlksx = scratch array of coefficients for computing
c                  equilibrium constants as a function of temperature
c       axvfsx = scratch array of coefficients for computing
c                  volume functions as a function of temperature
c       cdrsx  = scratch array of reaction coefficients
c       narxmx = number of coefficient elements of axlks per species
c                  per temperature range
c       ndrsx  = scratch array of indices of species corresponding to
c                  reaction coefficients
c       ndrsrx = scratch pointer array for indices of species appearing
c                  in reactions
c       nstmax = maximum number of species
c       ntprmx = maximum number of temperature ranges
c
c     Principal output:
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
c       ndrs   = standard  array of indices of species corresponding to
c                  reaction coefficients
c       ndrsr  = standard pointer array for indices of species
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
      integer ipchmx,ipcvmx,narxmx,ndrsmx,nstmax,ntprmx
c
      integer ndrs(ndrsmx),ndrsx(ndrsmx),ndrsr(2,nstmax),
     $ ndrsrx(2,nstmax)
      integer ipch,ipcv
c
      real*8 adhfs(narxmx,ntprmx,ipchmx,nstmax),
     $ adhfsx(narxmx,ntprmx,ipchmx,nstmax),
     $ advfs(narxmx,ntprmx,ipcvmx,nstmax),
     $ advfsx(narxmx,ntprmx,ipcvmx,nstmax),
     $ axhfs(narxmx,ntprmx,nstmax),axhfsx(narxmx,ntprmx,nstmax),
     $ axlks(narxmx,ntprmx,nstmax),axlksx(narxmx,ntprmx,nstmax),
     $ axvfs(narxmx,ntprmx,nstmax),axvfsx(narxmx,ntprmx,nstmax),
     $ cdrs(ndrsmx),cdrsx(ndrsmx)
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
      call copyaa(axlksx,axlks,nmax)
c
      if (ipch .ge. 0) then
        call copyaa(axhfsx,axhfs,nmax)
        if (ipch .ge. 1) call copyaa(adhfsx,adhfs,jhmax)
      endif
c
      if (ipcv .ge. 0) then
        call copyaa(axvfsx,axvfs,nmax)
        if (ipcv .ge. 1) call copyaa(advfsx,advfs,jvmax)
      endif
c
      nmax = 2*nstmax
      call copyia(ndrsrx,ndrsr,nmax)
c
      call copyaa(cdrsx,cdrs,ndrsmx)
      call copyia(ndrsx,ndrs,ndrsmx)
c
      end
