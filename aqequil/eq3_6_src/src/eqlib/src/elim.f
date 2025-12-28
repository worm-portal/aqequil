      subroutine elim(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,
     $ axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,eps100,ipch,ipchmx,
     $ ipcv,ipcvmx,jsflag,narxmx,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,
     $ nse,nst,nstmax,ntprmx,noutpt,nttyo,uspec)
c
c     This subroutine rewrites reaction equations so that the auxiliary
c     basis species with index nse and jflag = 30 is eliminated from
c     all reactions except the one linking it with its corresponding
c     strict basis variable. The polynomial coefficients for computing
c     the equilibrium coefficients are recomputed accordingly.
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
      integer ipchmx,ipcvmx,narxmx,ndrsmx,nstmax,ntprmx
c
      integer jsflag(nstmax),ndrs(ndrsmx),ndrsx(ndrsmx),
     $ ndrsr(2,nstmax),ndrsrx(2,nstmax)
      integer ipch,ipcv,nse,nst,noutpt,nttyo
c
      character*48 uspec(nstmax)
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
      integer i,ipc,j,jlen,n,nnx,nre1,nre2,nr1,nr2,ns,nss,nt
c
      character*56 uspn56
c
      real*8 cx,cxe,cxee,cxse,cxss,stofac
c
      real*8 coefdr
c
c-----------------------------------------------------------------------
c
c
      nt = ndrsr(2,nse) - ndrsr(1,nse) + 1
      if (nt .lt. 2) then
c
c       Calling sequence substitutions:
c         uspec(nse) for unam48
c
        call fmspnx(jlen,uspec(nse),uspn56)
        write (noutpt,1000) uspn56(1:jlen)
        write (nttyo,1000) uspn56(1:jlen)
 1000   format(/' * Error - (EQLIB/elim) The species ',a,
     $  /7x,'is in the strict basis and therefore can not be',
     $  ' eliminated',/7x,'from the working basis set.')
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calling sequence substitutions:
c       nse for ns
c
      cxee = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,nse,nstmax)
      nre1 = ndrsr(1,nse)
      nre2 = ndrsr(2,nse)
c
      nnx = 0
      do ns = 1,nst
        nr1 = ndrsr(1,ns)
        nr2 = ndrsr(2,ns)
        ndrsrx(1,ns) = nnx + 1
        cxe = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns,nstmax)
        if (ns.eq.nse .or. cxe.eq.0. .or. jsflag(ns).ge.2) then
c
c         Have found a reaction that is not to be changed.
c
          do n = nr1,nr2
            nnx = nnx + 1
            if (nnx .gt. ndrsmx) then
c
c             Calling sequence substitutions:
c               uspec(nse) for unam48
c
              call fmspnx(jlen,uspec(nse),uspn56)
              write (noutpt,1005) ndrsmx,uspn56(1:jlen),ndrsmx
              write (nttyo,1005) ndrsmx,uspn56(1:jlen),ndrsmx
 1005         format(/' * Error - (EQLIB/elim) The maximum ',i7,
     $        ' entries in the',/7x,'cdrs and ndrs arrays has been',
     $        ' exceeded in trying to eliminate',/7x,'the species ',
     $        a,' from the working basis set.',/7x,'Increase the',
     $        ' dimensioning  parameter ndrspa from its current',
     $        /7x,'value of ',i6,'.')
              stop
            endif
            cdrsx(nnx) = cdrs(n)
            ndrsx(nnx) = ndrs(n)
          enddo
          ndrsrx(2,ns) = nnx
c
c         Log K coefficients.
c
          do j = 1,ntprmx
            do i = 1,narxmx
              axlksx(i,j,ns) = axlks(i,j,ns)
            enddo
          enddo
c
          if (ipch .ge. 0) then
c
c           Enthalpy function coefficients.
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
c           Volume function coefficients.
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
          stofac = cxe/cxee
c
c         Do species appearing in the existing reaction for the ns-th
c         species.
c
          do n = nr1,nr2
            nss = ndrs(n)
c
c           Calling sequence substitutions:
c             nss for nse
c             nse for ns
c
            cxse = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nss,nse,nstmax)
            cx = cdrs(n) - stofac*cxse
            if (abs(cx) .le. eps100) cx = 0.
            if (cx .ne. 0.) then
              nnx = nnx + 1
              if (nnx .gt. ndrsmx) then
c
c               Calling sequence substitutions:
c                 uspec(nse) for unam48
c
                call fmspnx(jlen,uspec(nse),uspn56)
                write (noutpt,1005) ndrsmx,uspn56(1:jlen),ndrsmx
                write (nttyo,1005) ndrsmx,uspn56(1:jlen),ndrsmx
                stop
              endif
              cdrsx(nnx) = cx
              ndrsx(nnx) = nss
            endif
          enddo
c
c         Do species appearing in the existing reaction for the nse-th
c         species but not in that for the ns-th species.
c
          do n = nre1,nre2
            nss = ndrs(n)
c
c           Calling sequence substitutions:
c             nss for nse
c             nse for ns
c
            cxse = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nss,nse,nstmax)
c
c           Calling sequence substitutions:
c             nss for nse
c
            cxss = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nss,ns,nstmax)
            if (cxss .eq. 0.) then
              cx = -stofac*cxse
              if (abs(cx) .le. eps100) cx = 0.
              if (cx .ne. 0.) then
                nnx = nnx + 1
                if (nnx .gt. ndrsmx) then
c
c                 Calling sequence substitutions:
c                   uspec(nse) for unam48
c
                  call fmspnx(jlen,uspec(nse),uspn56)
                  write (noutpt,1005) ndrsmx,uspn56(1:jlen),ndrsmx
                  write (nttyo,1005) ndrsmx,uspn56(1:jlen),ndrsmx
                  stop
                endif
                cdrsx(nnx) = cx
                ndrsx(nnx) = nss
              endif
            endif
          enddo
c
c         If the new reaction has no entries, put in a null entry.
c
          if (nnx .lt. ndrsrx(1,ns)) then
            nnx = nnx + 1
            if (nnx .gt. ndrsmx) then
c
c             Calling sequence substitutions:
c               uspec(nse) for unam48
c
              call fmspnx(jlen,uspec(nse),uspn56)
              write (noutpt,1005) ndrsmx,uspn56(1:jlen),ndrsmx
              write (nttyo,1005) ndrsmx,uspn56(1:jlen),ndrsmx
              stop
            endif
            cdrsx(nnx) = 0.
            ndrsx(nnx) = 0
          endif
          ndrsrx(2,ns) = nnx
c
c         Log K coefficients.
c
          do j = 1,ntprmx
            if (axlks(1,j,ns) .lt. 9999999.) then
              do i = 1,narxmx
                axlksx(i,j,ns) = axlks(i,j,ns) - stofac*axlks(i,j,nse)
              enddo
            endif
          enddo
c
          if (ipch .ge. 0) then
c
c           Enthalpy function coefficients.
c
            do j = 1,ntprmx
              if (axhfs(1,j,ns) .lt. 9999999.) then
                do i = 1,narxmx
                  axhfsx(i,j,ns) = axhfs(i,j,ns) - stofac*axhfs(i,j,nse)
                enddo
              endif
            enddo
            do ipc = 1,ipch
              do j = 1,ntprmx
                do i = 1,narxmx
                  adhfsx(i,j,ipc,ns) = adhfs(i,j,ipc,ns)
     $            - stofac*adhfs(i,j,ipc,nse)
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
              if (axvfs(1,j,ns) .lt. 9999999.) then
                do i = 1,narxmx
                  axvfsx(i,j,ns) = axvfs(i,j,ns) - stofac*axvfs(i,j,nse)
                enddo
              endif
            enddo
            do ipc = 1,ipcv
              do j = 1,ntprmx
                do i = 1,narxmx
                  advfsx(i,j,ipc,ns) = advfs(i,j,ipc,ns)
     $            - stofac*advfs(i,j,ipc,nse)
                enddo
              enddo
            enddo
          endif
c
        endif
c
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Copy the new reactions into the standard arrays.
c
      call cdrscx(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,axlks,
     $ axlksx,axvfs,axvfsx,cdrs,cdrsx,ipch,ipchmx,ipcv,ipcvmx,
     $ narxmx,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,nstmax,ntprmx)
c
      end
