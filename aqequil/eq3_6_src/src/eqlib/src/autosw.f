      subroutine autosw(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,
     $ axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,eps100,ibswx,iindx1,
     $ ipch,ipchmx,ipcv,ipcvmx,jflag,jsflag,kbt,kmax,narn1,narxmx,
     $ nbasp,nbaspd,nbaspx,nbt,nbtmax,nbw,ndrs,ndrsmx,ndrsx,ndrsr,
     $ ndrsrx,noutpt,nst,nstmax,ntprmx,nttyo,qbassw,uspec,uzvec1)
c
c     This subroutine executes automatic basis switching for the purpose
c     of reducing mass balance residuals. EQLIB/fbassw.f finds
c     candidates for basis switching, and EQLIB/gabswx.f resolves any
c     conflicts.
c
c     This subroutine is called by:
c
c       EQLIB/absswa.f
c       EQ6/absswb.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       axlks  = array of polynomial coefficients for computing log K
c                  values (altered by this subroutine)
c       ibswx  = array defining switches to be made
c       nbasp  = array defining the species in the basis set
c                  (altered by this subroutine)
c       nbaspd = array defining the species in the data file basis set
c
c     Principal output:
c
c       axlks  = array of polynomial coefficients for computing log K
c                  values
c       cdrs   = array of reaction coefficients
c       nbasp  = array defining the species in the basis set
c       ndrs   = array of species indices corresponding to the reaction
c                  coefficients in the cdrs array
c       ndrsr  = array defining the range in the cdrs and ndrs arrays
c                  corresponding to the reaction for a given species
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ipchmx,ipcvmx,kmax,narxmx,nbtmax,ndrsmx,nstmax,ntprmx
c
      integer noutpt,nttyo
c
      integer ibswx(nbtmax),iindx1(kmax),jflag(nstmax),jsflag(nstmax),
     $ nbasp(nbtmax),nbaspd(nbtmax),nbaspx(nbtmax),ndrs(ndrsmx),
     $ ndrsx(ndrsmx),ndrsr(2,nstmax),ndrsrx(2,nstmax)
c
      integer ipch,ipcv,kbt,narn1,nbt,nbw,nst
c
      logical qbassw
c
      character*48 uspec(nstmax),uzvec1(kmax)
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
      real*8 eps100
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer nb,nsd,ns1,ns2,krow
c
      logical qbswok
c
c-----------------------------------------------------------------------
c
c     Save the nbasp array.
c
      call copyia(nbasp,nbaspx,nbt)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do krow = 1,kbt
        nb = iindx1(krow)
        ns2 = ibswx(nb)
        if (ns2 .gt. 0) then
          ns1 = nbaspx(nb)
          nsd = nbaspd(nb)
c
c         If the ns1-th species has not previously been switched with
c         the nsd-th species, go make the currently requested switch.
c         If it has been switched, first undo that switch. Take care
c         not to undo it twice.
c
          if (nsd .ne. ns2) then
            if (ns1 .ne. nsd) then
c
c             First undo the existing switch.
c
              nbasp(nb) = nsd
c
c             Calling sequence substitutions:
c               nsd for ns2
c
              call switch(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,
     $        axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,eps100,ipch,ipchmx,
     $        ipcv,ipcvmx,jflag,jsflag,narn1,narxmx,nbasp,nbaspd,nbaspx,
     $        nb,nbt,nbtmax,nbw,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,noutpt,
     $        nsd,nst,nstmax,ntprmx,nttyo,qbassw,qbswok,uspec)
            endif
          endif
c
          nbasp(nb) = ns2
c
          call switch(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,
     $    axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,eps100,ipch,ipchmx,
     $    ipcv,ipcvmx,jflag,jsflag,narn1,narxmx,nbasp,nbaspd,nbaspx,
     $    nb,nbt,nbtmax,nbw,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,noutpt,
     $    ns2,nst,nstmax,ntprmx,nttyo,qbassw,qbswok,uspec)
c
c         Update the names in the uzvec1 array.
c
          uzvec1(krow) = uspec(ns2)
        endif
      enddo
c
      end
