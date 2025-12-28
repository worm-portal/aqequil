      subroutine mdrgex(adhfs,adhfsd,adhfsx,advfs,advfsd,advfsx,
     $ axhfs,axhfsd,axhfsx,axlks,axlksd,axlksx,axvfs,axvfsd,axvfsx,
     $ cdrs,cdrsd,cdrsx,ipch,ipchmx,ipcv,ipcvmx,narxmx,nbasp,nbaspd,
     $ nbtmax,ndrs,ndrsd,ndrsx,ndrsmx,ndrsr,ndrsrd,ndrsrx,nern1,
     $ nern2,noutpt,nst,nstmax,ntprmx,nttyo)
c
c     This subroutine folds the reactions and reaction properties for
c     generic ion exchangers into the 'd' set. This is generally done
c     after those reactions and properties have been manipulated.
c
c     Note that the modified arrays must first be written into a set
c     of scratch arrays. This is because the ion exchanger sections
c     must generally go in the middle of the new arrays, and their
c     lengths are generally such that the positioning of the start of
c     the terminal "unmodified" sections is modified.
c
c     This subroutine is called by:
c
c       EQLIB/chsgex.f
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
      integer noutpt,nttyo
c
      integer nbasp(nbtmax),nbaspd(nbtmax),ndrs(ndrsmx),ndrsd(ndrsmx),
     $ ndrsx(ndrsmx),ndrsr(2,nstmax),ndrsrd(2,nstmax),ndrsrx(2,nstmax)
c
      integer ipch,ipcv,nern1,nern2,nst
c
      real*8 adhfs(narxmx,ntprmx,ipchmx,nstmax),
     $ adhfsd(narxmx,ntprmx,ipchmx,nstmax),
     $ adhfsx(narxmx,ntprmx,ipchmx,nstmax),
     $ advfs(narxmx,ntprmx,ipcvmx,nstmax),
     $ advfsd(narxmx,ntprmx,ipcvmx,nstmax),
     $ advfsx(narxmx,ntprmx,ipcvmx,nstmax),
     $ axhfs(narxmx,ntprmx,nstmax),axhfsd(narxmx,ntprmx,nstmax),
     $ axhfsx(narxmx,ntprmx,nstmax),axlks(narxmx,ntprmx,nstmax),
     $ axlksd(narxmx,ntprmx,nstmax),axlksx(narxmx,ntprmx,nstmax),
     $ axvfs(narxmx,ntprmx,nstmax),axvfsd(narxmx,ntprmx,nstmax),
     $ axvfsx(narxmx,ntprmx,nstmax),cdrs(ndrsmx),cdrsd(ndrsmx),
     $ cdrsx(ndrsmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ipc,j,jhmax,jvmax,j2,j3,n,nmax,nb,nn,nr1,nr2,ns,nt
c
      integer ilnobl
c
      character*8 ux8a,ux8b
c
c-----------------------------------------------------------------------
c
      nmax = narxmx*ntprmx*(nern1 -1)
      jhmax = narxmx*ntprmx*ipchmx*(nern1 -1)
      jvmax = narxmx*ntprmx*ipcvmx*(nern1 -1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Coefficients for equilibrium constants.
c
      if (nern1 .gt. 1) call copyaa(axlksd,axlksx,nmax)
c
      do ns = nern1,nern2
        do j = 1,ntprmx
          do i = 1,narxmx
            axlksx(i,j,ns) = axlks(i,j,ns)
          enddo
        enddo
      enddo
c
      if (nern2 .lt. nst) then
        do ns = nern2 + 1,nst
          do j = 1,ntprmx
            do i = 1,narxmx
              axlksx(i,j,ns) = axlksd(i,j,ns)
            enddo
          enddo
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Coefficients related to enthalpy.
c
      if (ipch .ge. 0) then
c
c       Coefficients for enthalpy functions.
c
        if (nern1 .gt. 1) call copyaa(axhfsd,axhfsx,nmax)
c
        do ns = nern1,nern2
          do j = 1,ntprmx
            do i = 1,narxmx
              axhfsx(i,j,ns) = axhfs(i,j,ns)
            enddo
          enddo
        enddo
c
        if (nern2 .lt. nst) then
          do ns = nern2 + 1,nst
            do j = 1,ntprmx
              do i = 1,narxmx
                axhfsx(i,j,ns) = axhfsd(i,j,ns)
              enddo
            enddo
          enddo
        endif
c
        if (ipch .ge. 1) then
c
c         Coefficients for derivatives of enthalpy functions.
c
          if (nern1 .gt. 1) call copyaa(adhfsd,adhfsx,jhmax)
c
          do ns = nern1,nern2
            do ipc = 1,ipchmx
              do j = 1,ntprmx
                do i = 1,narxmx
                  adhfsx(i,j,ipc,ns) = adhfs(i,j,ipc,ns)
                enddo
              enddo
            enddo
          enddo
c
          if (nern2 .lt. nst) then
            do ns = nern2 + 1,nst
              do ipc = 1,ipchmx
                do j = 1,ntprmx
                  do i = 1,narxmx
                    adhfsx(i,j,ipc,ns) = adhfsd(i,j,ipc,ns)
                  enddo
                enddo
              enddo
            enddo
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Coefficients related to volume.
c
      if (ipcv .ge. 0) then
c
c         Coefficients for volume functions.
c
        if (nern1 .gt. 1) call copyaa(axvfsd,axvfsx,nmax)
c
        do ns = nern1,nern2
          do j = 1,ntprmx
            do i = 1,narxmx
              axvfsx(i,j,ns) = axvfs(i,j,ns)
            enddo
          enddo
        enddo
c
        if (nern2 .lt. nst) then
          do ns = nern2 + 1,nst
            do j = 1,ntprmx
              do i = 1,narxmx
                axvfsx(i,j,ns) = axvfsd(i,j,ns)
              enddo
            enddo
          enddo
        endif
c
        if (ipcv .ge. 1) then
c
c         Coefficients for derivatives of volume functions.
c
          if (nern1 .gt. 1) call copyaa(advfsd,advfsx,jvmax)
c
          do ns = nern1,nern2
            do ipc = 1,ipcvmx
              do j = 1,ntprmx
                do i = 1,narxmx
                  advfsx(i,j,ipc,ns) = advfs(i,j,ipc,ns)
                enddo
              enddo
            enddo
          enddo
c
          if (nern2 .lt. nst) then
            do ns = nern2 + 1,nst
              do ipc = 1,ipcvmx
                do j = 1,ntprmx
                  do i = 1,narxmx
                    advfsx(i,j,ipc,ns) = advfsd(i,j,ipc,ns)
                  enddo
                enddo
              enddo
            enddo
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Index range pointers for reactions.
c
      if (nern1 .gt. 1) then
        nmax = 2*(nern1 - 1)
        call copyia(ndrsrd,ndrsrx,nmax)
      endif
c
      do ns = nern1,nern2
        ndrsrx(1,ns) = ndrsr(1,ns)
        ndrsrx(1,ns) = ndrsr(1,ns)
      enddo
c
      if (nern2 .lt. nst) then
        do ns = nern2 + 1,nst
          ndrsrx(1,ns) = ndrsrd(1,ns)
          ndrsrx(1,ns) = ndrsrd(1,ns)
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Reaction coefficients and indices of associated species.
c
      nt = 0
      if (nern1 .gt. 1) nt = ndrsrd(2,nern1 - 1)
      if (nern2 .ge. nern1)
     $  nt = nt + ndrsr(2,nern2) - ndrsr(1,nern1) + 1
      if (nst .gt. nern2)
     $ nt = nt + ndrsrd(2,nst) - ndrsrd(1,nern2 + 1) + 1
c
      if (nt .gt. ndrsmx) then
        write (ux8a,'(i5)') ndrsmx
        write (ux8b,'(i5)') nt
        j2 = ilnobl(ux8a)
        j3 = ilnobl(ux8b)
        write (noutpt,1000) ux8a(1:j2),ux8b(1:j3)
        write (nttyo,1000)  ux8a(1:j2),ux8b(1:j3)
 1000   format(/' * Error - (EQLIB/mdrgex.f) Have insufficient array',
     $  ' space to expand the',/7x,'reaction coefficient arrays to',
     $  ' accommodate changes to the',/7x,'ion exchanger section',
     $  " of the 'd' set of reaction data. Increase",/7x,'the',
     $  ' dimensioning parameter ndrspa from ',a,' to at least ',
     $  a,'.')
        stop
      endif
c
      if (nern1 .gt. 1) then
        nmax = ndrsrd(2,nern1 - 1)
        call copyaa(cdrsd,cdrsx,nmax)
        call copyia(ndrsd,ndrsx,nmax)
      endif
c
      if (nern1 .gt. 1) then
        nn = ndrsrd(2,nern1 - 1)
      else
        nn = 0
      endif
      do ns = nern1,nern2
        nr1 = ndrsr(1,ns)
        nr2 = ndrsr(2,ns)
        do n = nr1,nr2
          nn = nn + 1
          ndrsx(nn) = ndrs(n)
          cdrsx(nn) = cdrs(n)
        enddo
      enddo
c
      nn = ndrsrx(2,nern2)
      do ns = nern2 + 1,nst
        nr1 = ndrsrd(1,ns)
        nr2 = ndrsrd(2,ns)
        do n = nr1,nr2
          nn = nn + 1
          ndrsx(nn) = ndrsd(n)
          cdrsx(nn) = cdrsd(n)
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do nb = 1,nbtmax
        ns = nbasp(nb)
        if (ns.ge.nern1 .and. ns.le.nern2) nbaspd(nb) = nbasp(nb)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Copy the scratch arrays into the corresponding 'd' set arrays.
c
c     Copy the new reactions in the 'x' set into the standard arrays
c     in the 'd' set.
c
c     Calling sequence substitutions:
c       adhfsd for adhfs
c       advfsd for advfs
c       axhfsd for axhfs
c       axlksd for axlks
c       axvfsd for avhfs
c       cdrsd for cdrs
c       ndrsd for ndrs
c       ndrsrd for ndrsr
c
      call cdrscx(adhfsd,adhfsx,advfsd,advfsx,axhfsd,axhfsx,axlksd,
     $ axlksx,axvfsd,axvfsx,cdrsd,cdrsx,ipch,ipchmx,ipcv,ipcvmx,
     $ narxmx,ndrsd,ndrsmx,ndrsx,ndrsrd,ndrsrx,nstmax,ntprmx)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
