      subroutine switch(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,
     $ axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,eps100,ipch,ipchmx,
     $ ipcv,ipcvmx,jflag,jsflag,narn1,narxmx,nbasp,nbaspd,nbaspx,
     $ nb,nbt,nbtmax,nbw,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,noutpt,
     $ ns2,nst,nstmax,ntprmx,nttyo,qbassw,qbswok,uspec)
c
c     This subroutine executes an ordinary basis switch. The ns2-th
c     species (not in the current active basis set) is exchanged into
c     that set for the ns1-th species. All reactions are rewritten where
c     necessary for consistency. Here nbaspx is the unmodified copy of
c     the nbasp array. It contains a record of the basis set prior to
c     basis switching.
c
c     Basis switching is subject to the following rules:
c
c        1. Neither species may be suppressed.
c
c        2. The species switched in must have reaction that gives the
c           species to be switched out as a product.
c
c        3. The species to be switched into the active basis set must
c           not already occupy a position of its own in that set (e.g.,
c           be an active auxiliary basis species). To switch an active
c           auxiliary basis species into the strict basis set (e.g.,
c           make a special basis switch), use EQLIB/swtchb.f.
c
c     The species indices are not interchanged by the switch. If both
c     species are basis species, however, their basis indices are
c     interchanged.
c
c     This subroutine is called by:
c
c       EQLIB/autosw.f
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c       EQ6/eqcalc.f
c       EQ6/eqphas.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c
c     Principal output:
c
c       qbassw = logical flag:
c                  = .true. if the basis set after the current switch
c                    is not identical to the data file basis set
c       qbswok = logical flag:
c                  = .true. if the current switch was completed okay
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
      integer jflag(nstmax),jsflag(nstmax),nbasp(nbtmax),
     $ nbaspd(nbtmax),nbaspx(nbtmax),ndrs(ndrsmx),ndrsx(ndrsmx),
     $ ndrsr(2,nstmax),ndrsrx(2,nstmax)
c
      integer ipch,ipcv,narn1,nb,nbt,nbw,nst,ns2
c
      logical qbassw
c
      character(len=48) uspec(nstmax)
c
      real(8) adhfs(narxmx,ntprmx,ipchmx,nstmax),
     $ adhfsx(narxmx,ntprmx,ipchmx,nstmax),
     $ advfs(narxmx,ntprmx,ipcvmx,nstmax),
     $ advfsx(narxmx,ntprmx,ipcvmx,nstmax),
     $ axhfs(narxmx,ntprmx,nstmax),axhfsx(narxmx,ntprmx,nstmax),
     $ axlks(narxmx,ntprmx,nstmax),axlksx(narxmx,ntprmx,nstmax),
     $ axvfs(narxmx,ntprmx,nstmax),axvfsx(narxmx,ntprmx,nstmax),
     $ cdrs(ndrsmx),cdrsx(ndrsmx)
c
      real(8) eps100
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ipc,j,jlen,jlen1,jlen2,jfl,j2,n,nbb,nb2,nmax,nrf1,nrf2,
     $ nrl1,nrl2,nr1,nr2,ns,nse,nsi,ns1,ntf,nx
c
      integer ilnobl,nbasis
c
      logical qbswok
c
      character(len=56) uspn56,usp156,usp256
      character(len=8) ux8
c
      real(8) axx,cx,cxes,cxe1,cxe2,cx1s,cx11,cx12,stofac
c
      real(8) coefdr
c
c-----------------------------------------------------------------------
c
c     Zero scratch arrays.
c
      nmax = narxmx*ntprmx*ipchmx*nstmax
      call initaz(adhfsx,nmax)
c
      nmax = narxmx*ntprmx*ipcvmx*nstmax
      call initaz(advfsx,nmax)
c
      nmax = narxmx*ntprmx*nstmax
      call initaz(axhfsx,nmax)
      call initaz(axlksx,nmax)
      call initaz(axvfsx,nmax)
c
      nmax = 2*nstmax
      call initiz(ndrsrx,nmax)
c
      nmax = ndrsmx
      call initaz(cdrsx,nmax)
      call initiz(ndrsx,nmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      qbswok = .false.
c
      ns1 = nbaspx(nb)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check to see if the switch is okay.
c
      call swtchk(cdrs,jflag,jsflag,nbaspx,nbt,nbtmax,ndrs,ndrsmx,
     $ ndrsr,noutpt,ns1,ns2,nstmax,nttyo,uspec)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calling sequence substitutions:
c       jlen1 for jlen
c       uspec(ns1) for unam48
c       usp156 for uspn56
c
      call fmspnx(jlen1,uspec(ns1),usp156)
c
c     Calling sequence substitutions:
c       jlen2 for jlen
c       uspec(ns2) for unam48
c       usp256 for uspn56
c
      call fmspnx(jlen2,uspec(ns2),usp256)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check whether or not the species being switched in is not already
c     in the active basis set (associated with some other mass balance).
c     This is necessary to ensure that the switch to be made is an
c     ordinary basis switch (the only type allowed in the present
c     subroutine).
c
c     Calling sequence substitutions:
c       nbaspx for nbasp
c       ns2 for ns
c
      nb2 = nbasis(nbaspx,nbt,nbtmax,ns2)
c
      if (nb2.gt.0 .and. jflag(ns2).ne.30) then
        write (noutpt,1000) usp156(1:jlen1),usp256(1:jlen2)
        write (nttyo,1000) usp156(1:jlen1),usp256(1:jlen2)
 1000   format(/' * Error - (EQLIB/switch) Programming error trap:',
     $  " Can't replace",/7x,'the species ', a,' in the basis set with',
     $  ' ',a,',',/7x,'doing an ordinary basis switch, because the',
     $  ' former is already',/7x,'in the active basis set.')
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1010) usp156(1:jlen1),usp256(1:jlen2)
      write (nttyo,1010) usp156(1:jlen1),usp256(1:jlen2)
 1010 format(/' Making an ordinary basis switch: replacing ',a,
     $ /1x,'in the active basis set with ',a,'.',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nb .eq. nbw) then
c
c       Switching solvent water out of the basis set. Write a warning.
c
        write (noutpt,1020) usp156(1:jlen1)
        write (nttyo,1020) usp156(1:jlen1)
 1020   format(/' * Warning - (EQLIB/switch) Making an ordinary',
     $  ' basis switch',/7x,'which moves ',a,' out of the active',
     $  ' basis set.')
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1040)
 1040 format(/'  Before the switch:',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print the existing reaction for the ns1-th species, if any.
c
c     Calling sequence substitutions:
c       noutpt for nf
c       ns1 for ns
c
      call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns1,nstmax,uspec)
c
c     Print the existing linking reaction.
c
c     Calling sequence substitutions:
c       noutpt for nf
c       ns2 for ns
c
      call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns2,nstmax,uspec)
c
      write (noutpt,1050)
 1050 format(1x)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calling sequence substitutions:
c       ns1 for nse
c       ns2 for ns
c
      cx12 = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns1,ns2,nstmax)
      nrl1 = ndrsr(1,ns2)
      nrl2 = ndrsr(2,ns2)
c
      nrf1 = ndrsr(1,ns1)
      nrf2 = ndrsr(2,ns1)
      cx11 = cdrs(nrf1)
c
      write (ux8,'(i5)') ndrsmx
      call lejust(ux8)
      j2 = ilnobl(ux8)
c
      nx = 0
      do ns = 1,nst
        nr1 = ndrsr(1,ns)
        nr2 = ndrsr(2,ns)
        ndrsrx(1,ns) = nx + 1
        if (ns .eq. ns1) then
c
c         Invert the linking reaction. Put the coefficient for the
c         ns1-th species first in the range.
c
          nx = nx + 1
          if (nx .gt. ndrsmx) then
            write (noutpt,1070) ux8(1:j2),usp256(1:jlen2),
     $      usp156(1:jlen1),usp156(1:jlen1)
            write (nttyo,1070) ux8(1:j2),usp256(1:jlen2),
     $      usp156(1:jlen1),usp156(1:jlen1)
 1070       format(/' * Error - (EQLIB/switch) The maximum ',a,
     $      ' entries in the cdrs and ndrs',/7x,'arrays has been',
     $      ' exceeded in trying to switch ',a,/7x,'into the active',
     $      ' basis set for ',a,' while writing',/7x,'the reaction',
     $      ' for ',a,'. Increase the value of the',/7x,
     $      ' dimensioning parameter ndrspa.')
            stop
          endif
          cdrsx(nx) = -cx12
          ndrsx(nx) = ns
          do n = nrl1,nrl2
            if (ndrs(n) .ne. ns) then
              nx = nx + 1
              if (nx .gt. ndrsmx) then
                write (noutpt,1070) ux8(1:j2),usp256(1:jlen2),
     $          usp156(1:jlen1),usp156(1:jlen1)
                write (nttyo,1070) ux8(1:j2),usp256(1:jlen2),
     $          usp156(1:jlen1),usp156(1:jlen1)
                stop
              endif
              cdrsx(nx) = -cdrs(n)
              ndrsx(nx) = ndrs(n)
            endif
          enddo
c
c         Log K coefficients.
c
          do j = 1,ntprmx
            do i = 1,narxmx
              axlksx(i,j,ns) = -axlks(i,j,ns2)
            enddo
          enddo
c
          if (ipch .ge. 0) then
c
c           Enthalpy function coefficients.
c
            do j = 1,ntprmx
              do i = 1,narxmx
                axhfsx(i,j,ns) = -axhfs(i,j,ns2)
              enddo
            enddo
            do ipc = 1,ipch
              do j = 1,ntprmx
                do i = 1,narxmx
                  adhfsx(i,j,ipc,ns) = -adhfs(i,j,ipc,ns2)
                enddo
              enddo
            enddo
          endif
c
          if (ipcv .ge. 0) then
c
c           Volume function coefficients.
c
            do j = 1,ntprmx
              do i = 1,narxmx
                axvfsx(i,j,ns) = -axvfs(i,j,ns2)
              enddo
            enddo
            do ipc = 1,ipcv
              do j = 1,ntprmx
                do i = 1,narxmx
                  advfsx(i,j,ipc,ns) = -advfs(i,j,ipc,ns2)
                enddo
              enddo
            enddo
          endif
c
        elseif (ns .eq. ns2) then
c
          ntf = nrf2 - nrf1 + 1
          if (ntf .lt. 2) then
c
c           If the ns1-th species was a strict basis species, make
c           the ns2-th species into one. Write a null reaction.
c
            nx = nx + 1
            if (nx .gt. ndrsmx) then
              write (noutpt,1070) ux8(1:j2),usp256(1:jlen2),
     $        usp156(1:jlen1),usp256(1:jlen2)
              write (nttyo,1070) ux8(1:j2),usp256(1:jlen2),
     $        usp156(1:jlen1),usp256(1:jlen2)
              stop
            endif
            cdrsx(nx) = 0.
            ndrsx(nx) = 0
c
c           Log K coefficients. Log K = 0 for a strict basis species.
c
            do j = 1,ntprmx
              do i = 1,narxmx
                axlksx(i,j,ns) = 0.
              enddo
            enddo
c
            if (ipch .ge. 0) then
c
c             Enthalpy function coefficients. This function is the
c             partial molar enthalpy of formation for a strict
c             basis species.
c
              do j = 1,ntprmx
                do i = 1,narxmx
                  axx = axhfs(i,j,ns)
                  do n = nr1 + 1,nr2
                    nsi = ndrs(n)
                    cx = cdrs(n)
                    axx = axx - cx*axhfs(i,j,nsi)
                  enddo
                  cx = cdrs(nr1)
                  axhfsx(i,j,ns) = axx/cx
                enddo
              enddo
              do ipc = 1,ipch
                do j = 1,ntprmx
                  do i = 1,narxmx
                    axx = adhfs(i,j,ipc,ns)
                    do n = nr1 + 1,nr2
                      nsi = ndrs(n)
                      cx = cdrs(n)
                      axx = axx - cx*adhfs(i,j,ipc,nsi)
                    enddo
                    cx = cdrs(nr1)
                    adhfsx(i,j,ipc,ns) = axx/cx
                  enddo
                enddo
              enddo
            endif
c
            if (ipcv .ge. 0) then
c
c             Volume function coefficients. This function is the
c             partial molar volume for a strict basis species.
c
              do j = 1,ntprmx
                do i = 1,narxmx
                  axx = axvfs(i,j,ns)
                  do n = nr1 + 1,nr2
                    nsi = ndrs(n)
                    cx = cdrs(n)
                    axx = axx - cx*axvfs(i,j,nsi)
                  enddo
                  cx = cdrs(nr1)
                  axvfsx(i,j,ns) = axx/cx
                enddo
              enddo
              do ipc = 1,ipcv
                do j = 1,ntprmx
                  do i = 1,narxmx
                    axx = advfs(i,j,ipc,ns)
                    do n = nr1 + 1,nr2
                      nsi = ndrs(n)
                      cx = cdrs(n)
                      axx = axx - cx*advfs(i,j,ipc,nsi)
                    enddo
                    cx = cdrs(nr1)
                    advfsx(i,j,ipc,ns) = axx/cx
                  enddo
                enddo
              enddo
            endif
c
          else
c
c           If the ns1-th species was an auxiliary basis species, make
c           the ns2-th species into one. Write a reaction which is the
c           original reaction for the breakdown of the the ns1-th
c           species, but using the original reaction for the ns1-th
c           species instead of the linking reaction to eliminate the
c           ns1-th species.
c
            stofac = -cx12/cx11
c
c           Go through the species of the existing reaction for the
c           ns2-th species (this is the linking reaction).
c
            do n = nr1,nr2
              nse = ndrs(n)
c
c             Calling sequence substitutions:
c               ns1 for ns
c
              cxe1 = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns1,nstmax)
              cx = cdrs(n) + stofac*cxe1
              if (abs(cx) .gt. eps100) then
                nx = nx + 1
                if (nx .gt. ndrsmx) then
                  write (noutpt,1070) ux8(1:j2),usp256(1:jlen2),
     $            usp156(1:jlen1),usp256(1:jlen2)
                  write (nttyo,1070) ux8(1:j2),usp256(1:jlen2),
     $            usp156(1:jlen1),usp256(1:jlen2)
                  stop
                endif
                cdrsx(nx) = cx
                ndrsx(nx) = nse
              endif
            enddo
c
c           Add to the new reaction any species in the original
c           reaction for the breakdown of the ns1-th species.
c
            do n = nrf1,nrf2
              nse = ndrs(n)
              cxes = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns,nstmax)
              if (cxes .eq. 0.) then
                cx = stofac*cdrs(n)
                nx = nx + 1
                if (nx .gt. ndrsmx) then
                  write (noutpt,1070) ux8(1:j2),usp256(1:jlen2),
     $            usp156(1:jlen1),usp256(1:jlen2)
                  write (nttyo,1070) ux8(1:j2),usp256(1:jlen2),
     $            usp156(1:jlen1),usp256(1:jlen2)
                  stop
                endif
                cdrsx(nx) = cx
                ndrsx(nx) = nse
              endif
            enddo
c
c           Log K coefficients.
c
            do j = 1,ntprmx
              do i = 1,narxmx
                axlksx(i,j,ns) = axlks(i,j,ns) + stofac*axlks(i,j,ns1)
              enddo
            enddo
c
            if (ipch .ge. 0) then
c
c             Enthalpy function coefficients.
c
              do j = 1,ntprmx
                do i = 1,narxmx
                  axhfsx(i,j,ns) = axhfs(i,j,ns) + stofac*axhfs(i,j,ns1)
                enddo
              enddo
              do ipc = 1,ipch
                do j = 1,ntprmx
                  do i = 1,narxmx
                    adhfsx(i,j,ipc,ns) = adhfs(i,j,ipc,ns)
     $              + stofac*adhfs(i,j,ipc,ns1)
                  enddo
                enddo
              enddo
            endif
c
            if (ipcv .ge. 0) then
c
c             Volume function coefficients.
c
              do j = 1,ntprmx
                do i = 1,narxmx
                  axvfsx(i,j,ns) = axvfs(i,j,ns)
     $            + stofac*axvfs(i,j,ns1)
                enddo
              enddo
              do ipc = 1,ipcv
                do j = 1,ntprmx
                  do i = 1,narxmx
                    advfsx(i,j,ipc,ns) = advfs(i,j,ipc,ns)
     $              + stofac*advfs(i,j,ipc,ns1)
                  enddo
                enddo
              enddo
            endif
          endif
c
        else
c
c         Do reactions for all other species.
c
c         Calling sequence substitutions:
c           uspec(ns) for unam48
c
          call fmspnx(jlen,uspec(ns),uspn56)
c
c         Calling sequence substitutions:
c           ns1 for nse
c
          cx1s = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns1,ns,nstmax)
          if (cx1s .eq. 0.) then
c
c           Have a reaction which is not connected. Copy as is.
c
            do n = nr1,nr2
              nx = nx + 1
              if (nx .gt. ndrsmx) then
              write (noutpt,1070) ux8(1:j2),usp256(1:jlen2),
     $          usp156(1:jlen1),uspn56(1:jlen)
                write (nttyo,1070) ux8(1:j2),usp256(1:jlen2),
     $          usp156(1:jlen1),uspn56(1:jlen)
                stop
              endif
              cdrsx(nx) = cdrs(n)
              ndrsx(nx) = ndrs(n)
            enddo
c
c           Log K coefficients.
c
            do j = 1,ntprmx
              do i = 1,narxmx
                axlksx(i,j,ns) = axlks(i,j,ns)
              enddo
            enddo
c
            if (ipch .ge. 0) then
c
c             Enthalpy function coefficients.
c
              do j = 1,ntprmx
                do i = 1,narxmx
                  axhfsx(i,j,ns) = axhfs(i,j,ns)
                enddo
              enddo
              do ipc = 1,ipch
                do j = 1,ntprmx
                  do i = 1,narxmx
                    adhfsx(i,j,ipc,ns) = adhfs(i,j,ipc,ns)
                  enddo
                enddo
              enddo
            endif
c
            if (ipcv .ge. 0) then
c
c             Volume function coefficients.
c
              do j = 1,ntprmx
                do i = 1,narxmx
                  axvfsx(i,j,ns) = axvfs(i,j,ns)
                enddo
              enddo
              do ipc = 1,ipcv
                do j = 1,ntprmx
                  do i = 1,narxmx
                    advfsx(i,j,ipc,ns) = advfs(i,j,ipc,ns)
                  enddo
                enddo
              enddo
            endif
c
          else
c
c           Have a reaction which is connected.
c
            stofac = -cx1s/cx12
c
c           Go through the species of the existing reaction for the
c           ns-th species.
c
            do n = nr1,nr2
              nse = ndrs(n)
c
c             Calling sequence substitutions:
c               ns2 for ns
c
              cxe2 = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns2,nstmax)
              cx = cdrs(n) + stofac*cxe2
              if (abs(cx) .gt. eps100) then
                nx = nx + 1
                if (nx .gt. ndrsmx) then
                  write (noutpt,1070) ux8(1:j2),usp256(1:jlen2),
     $            usp156(1:jlen1),uspn56(1:jlen)
                  write (nttyo,1070) ux8(1:j2),usp256(1:jlen2),
     $            usp156(1:jlen1),uspn56(1:jlen)
                  stop
                endif
                cdrsx(nx) = cx
                ndrsx(nx) = nse
              endif
            enddo
c
c           Add to the new reaction any species in the linking reaction
c           that were not in the original reaction for the breakdown of
c           the ns-th species.
c
            do n = nrl1,nrl2
              nse = ndrs(n)
              cxes = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns,nstmax)
              if (cxes .eq. 0.) then
                cx = stofac*cdrs(n)
                nx = nx + 1
                if (nx .gt. ndrsmx) then
                  write (noutpt,1070) ux8(1:j2),usp256(1:jlen2),
     $            usp156(1:jlen1),uspn56(1:jlen)
                  write (nttyo,1070) ux8(1:j2),usp256(1:jlen2),
     $            usp156(1:jlen1),uspn56(1:jlen)
                  stop
                endif
                cdrsx(nx) = cx
                ndrsx(nx) = nse
              endif
            enddo
c
c           Log K coefficients.
c
            do j = 1,ntprmx
              if (axlks(1,j,ns) .lt. 9999999.) then
                do i = 1,narxmx
                  axlksx(i,j,ns) = axlks(i,j,ns) + stofac*axlks(i,j,ns2)
                enddo
              endif
            enddo
c
            if (ipch .ge. 0) then
c
c             Enthalpy function coefficients.
c
              do j = 1,ntprmx
                if (axhfs(1,j,ns) .lt. 9999999.) then
                  do i = 1,narxmx
                    axhfsx(i,j,ns) = axhfs(i,j,ns)
     $              + stofac*axhfs(i,j,ns2)
                  enddo
                endif
              enddo
              do ipc = 1,ipch
                do j = 1,ntprmx
                  do i = 1,narxmx
                    adhfsx(i,j,ipc,ns) = adhfs(i,j,ipc,ns)
     $              + stofac*adhfs(i,j,ipc,ns2)
                  enddo
                enddo
              enddo
            endif
c
            if (ipcv .ge. 0) then
c
c             Enthalpy function coefficients.
c
              do j = 1,ntprmx
                do i = 1,narxmx
                  if (axvfs(1,j,ns) .lt. 9999999.) then
                    axvfsx(i,j,ns) = axvfs(i,j,ns)
     $              + stofac*axvfs(i,j,ns2)
                  endif
                enddo
              enddo
              do ipc = 1,ipcv
                do j = 1,ntprmx
                  do i = 1,narxmx
                    advfsx(i,j,ipc,ns) = advfs(i,j,ipc,ns)
     $              + stofac*advfs(i,j,ipc,ns2)
                  enddo
                enddo
              enddo
            endif
c
          endif
c
        endif
        ndrsrx(2,ns) = nx
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Copy the new reactions from the 'x' arrays into the standard
c     arrays.
c
      call cdrscx(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,axlks,
     $ axlksx,axvfs,axvfsx,cdrs,cdrsx,ipch,ipchmx,ipcv,ipcvmx,
     $ narxmx,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,nstmax,ntprmx)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1090)
 1090 format(/'  After the switch:',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print the new linking reaction.
c
c     Calling sequence substitutions:
c       noutpt for nf
c       ns1 for ns
c
      call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns1,nstmax,uspec)
c
c     Print the new reaction for the ns2-th species.
c
c     Calling sequence substitutions:
c       noutpt for nf
c       ns2 for ns
c
      call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns2,nstmax,uspec)
c
      write (noutpt,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Exchange jflag values.
c
      jfl = jflag(ns1)
      jflag(ns1) = jflag(ns2)
      jflag(ns2) = jfl
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Clear the content in the nbasp and nbaspx arrays which constitutes
c     the instruction to make the current switch. Change the nbaspx
c     entry to clear the instruction for the current switch.
c
      nbaspx(nb) = ns2
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     If necessary, reset the basis index of water (nbw).
c
      if (ns1 .eq. narn1) then
        nbw = 0
      elseif (ns2 .eq. narn1) then
        nbw = nb
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set the flag indicating that the switch completed successfully.
c
      qbswok = .true.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check the value of the qbassw flag. Set it to .true. if
c     the basis set is not identical to the corresponding data file
c     basis set.
c
      qbassw = .false.
      do nbb = 1,nbt
        ns1 = nbaspd(nbb)
        ns2 = nbasp(nbb)
        if (ns2 .ne. ns1) go to 200
      enddo
      qbassw = .true.
  200 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
