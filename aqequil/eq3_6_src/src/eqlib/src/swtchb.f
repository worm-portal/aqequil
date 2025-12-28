      subroutine swtchb(adhfsd,adhfsx,advfsd,advfsx,axhfsd,axhfsx,
     $ axlksd,axlksx,axvfsd,axvfsx,cdrsd,cdrsx,ipch,ipchmx,ipcv,
     $ ipcvmx,narxmx,nbaspd,nbtmax,nbw,nb1,nb2,ndrsd,ndrsmx,ndrsx,
     $ ndrsrd,ndrsrx,noutpt,ns1,ns2,nsta,nstmax,ntprmx,nttyo,uspeca)
c
c     This subroutine performs a special kind of basis switch. It
c     interchanges the status of a strict basis species with an
c     auxiliary basis species. The effect of this action is equivalent
c     to modifying the data file so that a different basis species is
c     the strict basis species corresponding to a chemical element.
c
c     Here nb1 is the basis index of the species originally in the
c     strict set, and nb2 is the basis index of the species originally
c     in the auxiliary basis set. Here also ns1 and ns2 are the
c     corresponding species indices. The species indices are not
c     interchanged by the switch. The basis indices, as represented by
c     the contents of the nbaspd array, are interchanged.
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
      integer ipchmx,ipcvmx,narxmx,nbtmax,ndrsmx,nstmax,ntprmx
c
      integer noutpt,nttyo
c
      integer nbaspd(nbtmax),ndrsd(ndrsmx),ndrsx(ndrsmx),
     $ ndrsrd(2,nstmax),ndrsrx(2,nstmax)
      integer ipch,ipcv,nbw,nb1,nb2,ns1,ns2,nsta
c
      character*48 uspeca(nstmax)
c
      real*8 adhfsd(narxmx,ntprmx,ipchmx,nstmax),
     $ adhfsx(narxmx,ntprmx,ipchmx,nstmax),
     $ advfsd(narxmx,ntprmx,ipcvmx,nstmax),
     $ advfsx(narxmx,ntprmx,ipcvmx,nstmax),
     $ axhfsd(narxmx,ntprmx,nstmax),axhfsx(narxmx,ntprmx,nstmax),
     $ axlksd(narxmx,ntprmx,nstmax),axlksx(narxmx,ntprmx,nstmax),
     $ axvfsd(narxmx,ntprmx,nstmax),axvfsx(narxmx,ntprmx,nstmax),
     $ cdrsd(ndrsmx),cdrsx(ndrsmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ipc,j,jlen,jlen1,jlen2,n,nerr,nrl1,nrl2,nr1,nr2,ns,nsi,
     $ nt1,nt2,nx
c
      character*56 uspn56,usp156,usp256
c
      real*8 axx,cx,cx12
c
      real*8 coefdr
c
c-----------------------------------------------------------------------
c
      nerr = 0
c
c     Calling sequence substitutions:
c       jlen1 for jlen
c       uspeca(ns1) for unam48
c       usp156 for uspn56
c
      call fmspnx(jlen1,uspeca(ns1),usp156)
c
c     Calling sequence substitutions:
c       jlen2 for jlen
c       uspeca(ns2) for unam48
c       usp256 for uspn56
c
      call fmspnx(jlen2,uspeca(ns2),usp256)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nb1 .eq. nb2) then
        write (noutpt,1004) usp156(1:jlen1)
        write (nttyo,1004) usp156(1:jlen1)
 1004   format(/" * Error - (EQLIB/swtchb) Can't make a special basis",
     $  ' switch',/7x,'replacing ',a,' with itself.')
        nerr = nerr + 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nb2 .lt. nb1) then
        write (noutpt,1006) usp156(1:jlen1),usp256(1:jlen2)
        write (nttyo,1006) usp156(1:jlen1),usp256(1:jlen2)
 1006   format(/" * Error - (EQLIB/swtchb) Can't make a special",
     $  ' switch',/7x,'replacing ',a,' with ',a,' because these',
     $  ' species are',/7x,'not in the proper hierarchical order.',
     $  ' The former must appear before',/7x,'the latter in the list',
     $  ' of basis species at the time the special switch',
     $  /7x,'is executed. Check the ordering on the data file and any',
     $  /7x,'changes made by previously executed special basis',
     $  ' switches.')
        nerr = nerr + 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make sure that the first species is in the strict basis set.
c
      nt1 = ndrsrd(2,ns1) - ndrsrd(1,ns1) + 1
      if (nt1 .ge. 2) then
        write (noutpt,1010) usp156(1:jlen1),usp256(1:jlen2)
        write (nttyo,1010) usp156(1:jlen1),usp256(1:jlen2)
 1010   format(/" * Error - (EQLIB/swtchb) Can't make a special",
     $  ' switch',/7x,'replacing ',a,' with ',a,' because the',
     $  ' former species',/7x,'is not currently in the strict',
     $  ' basis set.',/7x,'Check the status of this species on',
     $  ' the data file and any',/7x,'changes made by previously',
     $  ' executed special basis switches.')
        nerr = nerr + 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check the linking reaction.
c
      nt2 = ndrsrd(2,ns2) - ndrsrd(1,ns2) + 1
      if (nt2 .lt. 2) then
c
c       The second species must not be in the strict basis, because
c       there is then no possibility of a linking reaction.
c
        write (noutpt,1030) usp156(1:jlen1),usp256(1:jlen2)
        write (nttyo,1030) usp156(1:jlen1),usp256(1:jlen2)
 1030   format(/" * Error - (EQLIB/swtchb) Can't make a special",
     $  ' switch',/7x,'replacing ',a,' with ',a,' because the',
     $  ' latter species',/7x,'is currently in the strict basis',
     $  ' set, precluding.',/7x,'the possibility of it having a',
     $  ' reaction linking it to the',/7x,'former. Check the status',
     $  ' of the latter species on the',/7x,'data file and any',
     $  ' changes made by previously executed',/7x,'special',
     $  /7x,'basis switches.')
        nerr = nerr + 1
      else
c
c       Make sure that the first species appears as a product in
c       the reaction belonging to the second species.
c
c       Calling sequence substitutions:
c         cdrsd for cdrs
c         ndrsd for ndrs
c         ndrsrd for ndrsr
c         ns1 for nse
c         ns2 for ns
c
        cx12 = coefdr(cdrsd,ndrsd,ndrsmx,ndrsrd,ns1,ns2,nstmax)
        if (cx12 .eq. 0.) then
          write (noutpt,1040) usp156(1:jlen1),usp256(1:jlen2)
          write (nttyo,1040) usp156(1:jlen1),usp256(1:jlen2)
 1040     format(/" * Error - (EQLIB/swtchb) Can't make a special",
     $    ' switch',/7x,'replacing ',a,' with ',a,' because the',
     $    ' former species',/7x,'does not appear in the reaction',
     $    ' for the former species. Check',/7x,'the reaction',
     $    ' on the data file and any changes made by',
     $    /7x,'previously executed special basis switches.')
          ns = ns2
c
c         Calling sequence substitutions:
c           cdrsd for cdrs
c           ndrsd for ndrs
c           ndrsrd for ndrsr
c           noutpt for nf
c           uspeca for uspec
c
          call prreac(cdrsd,ndrsd,ndrsmx,ndrsrd,noutpt,ns,nstmax,uspeca)
          nerr = nerr + 1
        endif
      endif
c
      if (nerr .gt. 0) stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1050) usp156(1:jlen1),usp256(1:jlen2)
      write (nttyo,1050) usp156(1:jlen1),usp256(1:jlen2)
 1050 format(/' Making a special basis switch: replacing ',a,
     $ ' with ',a,'.'/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nb1 .eq. nbw) then
        write (noutpt,1052) usp156(1:jlen1)
        write (nttyo,1052) usp156(1:jlen1)
 1052   format(/' * Warning - (EQLIB/swtchb) Making a special basis',
     $  ' switch',/7x,'replacing ',a,' with another species.')
      endif
c
      if (nb2 .eq. nbw) then
        write (noutpt,1054) usp256(1:jlen2)
        write (nttyo,1054) usp256(1:jlen2)
 1054   format(/' * Warning - (EQLIB/swtchb) Making a special basis',
     $  ' switch',/7x,'replacing a species with ',a,'.')
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print the existing linking reaction.
c
c     Calling sequence substitutions:
c       cdrsd for cdrs
c       ndrsd for ndrs
c       ndrsrd for ndrsr
c       noutpt for nf
c       ns2 for ns
c       uspeca for uspec
c
      call prreac(cdrsd,ndrsd,ndrsmx,ndrsrd,noutpt,ns2,nstmax,uspeca)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      nrl1 = ndrsrd(1,ns2)
      nrl2 = ndrsrd(2,ns2)
c
      nx = 0
      do ns = 1,nsta
        nr1 = ndrsrd(1,ns)
        nr2 = ndrsrd(2,ns)
        ndrsrx(1,ns) = nx + 1
c
c       Calling sequence substitutions:
c         uspeca(ns) for unam48
c
        call fmspnx(jlen,uspeca(ns),uspn56)
c
        if (ns .eq. ns1) then
c
c         Invert the linking reaction. Put the coefficient for the
c         ns1-th species first in the range.
c
          nx = nx + 1
          if (nx .gt. ndrsmx) then
            write (noutpt,1060) ndrsmx,usp156(1:jlen1),usp256(1:jlen2),
     $      uspn56(1:jlen),ndrsmx
            write (nttyo,1060) ndrsmx,usp156(1:jlen1),usp256(1:jlen2),
     $      uspn56(1:jlen),ndrsmx
 1060       format(/' * Error - (EQLIB/swtchb) The maximum ',i7,
     $      ' entries in the',/7x,'cdrsd and ndrsd arrays has been',
     $      ' exceeded in trying to',/7x,'make a special basis switch',
     $      ' replacing ',a,' with',/7x,a,' while writing the reaction',
     $      ' for ',a,'. Increase',/7x,'the dimensioning parameter',
     $      ' ndrspa from its present value of ',i7,'.')
            stop
          endif
          cdrsx(nx) = -cx12
          ndrsx(nx) = ns
          do n = nrl1,nrl2
            if (ndrsd(n) .ne. ns) then
              nx = nx + 1
              if (nx .gt. ndrsmx) then
                write (noutpt,1060) ndrsmx,usp156(1:jlen1),
     $          usp256(1:jlen2),uspn56(1:jlen),ndrsmx
                write (nttyo,1060) ndrsmx,usp156(1:jlen1),
     $          usp256(1:jlen2),uspn56(1:jlen),ndrsmx
                stop
              endif
              cdrsx(nx) = -cdrsd(n)
              ndrsx(nx) = ndrsd(n)
            endif
          enddo
          ndrsrx(2,ns) = nx
c
c         Log K coefficients.
c
          do j = 1,ntprmx
            do i = 1,narxmx
              axlksx(i,j,ns) = -axlksd(i,j,ns2)
            enddo
          enddo
c
          if (ipch .ge. 0) then
c
c           Enthalpy function coefficients.
c
            do j = 1,ntprmx
              do i = 1,narxmx
                axhfsx(i,j,ns) = -axhfsd(i,j,ns2)
              enddo
            enddo
            do ipc = 1,ipch
              do j = 1,ntprmx
                do i = 1,narxmx
                  adhfsx(i,j,ipc,ns) = -adhfsd(i,j,ipc,ns2)
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
                axvfsx(i,j,ns) = -axvfsd(i,j,ns2)
              enddo
            enddo
            do ipc = 1,ipcv
              do j = 1,ntprmx
                do i = 1,narxmx
                  advfsx(i,j,ipc,ns) = -advfsd(i,j,ipc,ns2)
                enddo
              enddo
            enddo
          endif
c
        elseif (ns .eq. ns2) then
c
c         Write a null reaction.
c
          nx = nx + 1
          if (nx .gt. ndrsmx) then
            write (noutpt,1060) ndrsmx,usp156(1:jlen1),
     $      usp256(1:jlen2),uspn56(1:jlen),ndrsmx
            write (nttyo,1060) ndrsmx,usp156(1:jlen1),
     $      usp256(1:jlen2),uspn56(1:jlen),ndrsmx
            stop
          endif
          cdrsx(nx) = 0.
          ndrsx(nx) = 0
          ndrsrx(2,ns) = nx
c
c         Log K coefficients. Log K = 0 for a strict basis species.
c
          do j = 1,ntprmx
            do i = 1,narxmx
              axlksx(i,j,ns) = 0.
            enddo
          enddo
c
          if (ipch .ge. 0) then
c
c           Enthalpy function coefficients. This function is the
c           partial molar enthalpy of formation for a strict
c           basis species.
c
            do j = 1,ntprmx
              do i = 1,narxmx
                axx = axhfsd(i,j,ns)
                do n = nr1 + 1,nr2
                  nsi = ndrsd(n)
                  cx = cdrsd(n)
                  axx = axx - cx*axhfsd(i,j,nsi)
                enddo
                cx = cdrsd(nr1)
                axhfsx(i,j,ns) = axx/cx
              enddo
            enddo
            do ipc = 1,ipch
              do j = 1,ntprmx
                do i = 1,narxmx
                  axx = adhfsd(i,j,ipc,ns)
                  do n = nr1 + 1,nr2
                    nsi = ndrsd(n)
                    cx = cdrsd(n)
                    axx = axx - cx*adhfsd(i,j,ipc,nsi)
                  enddo
                  cx = cdrsd(nr1)
                  adhfsx(i,j,ipc,ns) = axx/cx
                enddo
              enddo
            enddo
          endif
c
          if (ipcv .ge. 0) then
c
c           Volume function coefficients. This function is the
c           partial molar volume for a strict basis species.
c
            do j = 1,ntprmx
              do i = 1,narxmx
                axx = axvfsd(i,j,ns)
                do n = nr1 + 1,nr2
                  nsi = ndrsd(n)
                  cx = cdrsd(n)
                  axx = axx - cx*axvfsd(i,j,nsi)
                enddo
                cx = cdrsd(nr1)
                axvfsx(i,j,ns) = axx/cx
              enddo
            enddo
            do ipc = 1,ipcv
              do j = 1,ntprmx
                do i = 1,narxmx
                  axx = advfsd(i,j,ipc,ns)
                  do n = nr1 + 1,nr2
                    nsi = ndrsd(n)
                    cx = cdrsd(n)
                    axx = axx - cx*advfsd(i,j,ipc,nsi)
                  enddo
                  cx = cdrsd(nr1)
                  advfsx(i,j,ipc,ns) = axx/cx
                enddo
              enddo
            enddo
          endif
c
        else
c
c         Copy as is reactions for all other species.
c
          do n = nr1,nr2
            nx = nx + 1
            if (nx .gt. ndrsmx) then
              write (noutpt,1060) ndrsmx,usp156(1:jlen1),
     $        usp256(1:jlen2),uspn56(1:jlen),ndrsmx
              write (nttyo,1060) ndrsmx,usp156(1:jlen1),
     $        usp256(1:jlen2),uspn56(1:jlen),ndrsmx
              stop
            endif
            cdrsx(nx) = cdrsd(n)
            ndrsx(nx) = ndrsd(n)
          enddo
          ndrsrx(2,ns) = nx
c
c         Log K coefficients.
c
          do j = 1,ntprmx
            do i = 1,narxmx
              axlksx(i,j,ns) = axlksd(i,j,ns)
            enddo
          enddo
c
          if (ipch .ge. 0) then
c
c           Enthalpy function coefficients.
c
            do j = 1,ntprmx
              do i = 1,narxmx
                axhfsx(i,j,ns) = axhfsd(i,j,ns)
              enddo
            enddo
            do ipc = 1,ipch
              do j = 1,ntprmx
                do i = 1,narxmx
                  adhfsx(i,j,ipc,ns) = adhfsd(i,j,ipc,ns)
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
                axvfsx(i,j,ns) = axvfsd(i,j,ns)
              enddo
            enddo
            do ipc = 1,ipcv
              do j = 1,ntprmx
                do i = 1,narxmx
                  advfsx(i,j,ipc,ns) = advfsd(i,j,ipc,ns)
                enddo
              enddo
            enddo
          endif
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
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
c     Print the new linking reaction.
c
c     Calling sequence substitutions:
c       cdrsd for cdrs
c       ndrsd for ndrs
c       ndrsrd for ndrsr
c       ns1 for ns
c       noutpt for nf
c       uspeca for uspec
c
      call prreac(cdrsd,ndrsd,ndrsmx,ndrsrd,noutpt,ns1,nstmax,uspeca)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Interchange the basis indices.
c
      nbaspd(nb1) = ns2
      nbaspd(nb2) = ns1
c
      if (nb1 .eq. nbw) then
        nbw = nb2
      elseif (nb2 .eq. nbw) then
        nbw = nb1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
