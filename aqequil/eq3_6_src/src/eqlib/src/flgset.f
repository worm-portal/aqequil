      subroutine flgset(axlksd,iopt,jflag,jpflag,jsflag,kxmod,narn1a,
     $ narn2a,narxmx,nbaspd,nbtd,nbtmax,ncmpra,ncta,ndrsd,ndrsmx,
     $ ndrsrd,noptmx,noutpt,npta,nptmax,nrdxsp,nsta,nstmax,ntpr,
     $ ntprmx,nttyo,nxmdmx,nxmod,uphasa,uptypa,uspeca,uxmod)
c
c     This subroutine sets up the status arrays jpflag and jsflag.
c     These flags denote the statuses, respectively, of phases
c     and species. The relevant values and their meanings are
c     as follows:
c
c        jpflag:
c          = 0   The species is present or potentially present
c          = 1   The physical presence of the phase in the model
c                  is suppressed; associated variables, such as a
c                  reaction affinity, may be calculated
c          = 2   The phase is completely ignored by the code
c
c        jsflag:
c          = 0   The species is present or potentially present
c          = 1   The physical presence of the species in the model
c                  is suppressed; associated variables, such as a
c                  reaction affinity, may be calculated
c          = 2   The species is completely ignored by the code
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
c       jpflag = phase status flag array
c       jsflag = species status flag array
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer narxmx,nbtmax,ndrsmx,noptmx,nptmax,nstmax,ntprmx,nxmdmx
c
      integer noutpt,nttyo
c
      integer iopt(noptmx),jflag(nstmax),jpflag(nptmax),
     $ jsflag(nstmax),kxmod(nxmdmx),nbaspd(nbtmax),ncmpra(2,nptmax),
     $ ndrsd(ndrsmx),ndrsrd(2,nstmax)
c
      integer narn1a,narn2a,nbtd,ncta,npta,nrdxsp,nsta,ntpr,nxmod
c
      character(len=48) uspeca(nstmax),uxmod(nxmdmx)
      character(len=24) uphasa(nptmax),uptypa(nptmax)
c
      real(8) axlksd(narxmx,ntprmx,nstmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jlen,n,nb,ncount,nn,np,nr1,nr2,ns,nse,nss,nt
c
      integer nbasis
c
      logical qheadr
c
      character(len=56) uspn56
      character(len=24) uptsld,uptliq,ux24,uy24
c
c-----------------------------------------------------------------------
c
      data uptsld /'Solid                   '/
      data uptliq /'Liquid                  '/
c
c-----------------------------------------------------------------------
c
c     Note: the following statements don't really do anything except
c     cause the compiler not to complain that uptliq is not used.
c     are not used.
c
      ux24 = uptliq
      uy24 = ux24
      uptliq = uy24
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Zero the jsflag and jpflag arrays.
c
      do ns = 1,nsta
        jsflag(ns) = 0
      enddo
      do np = 1,npta
        jpflag(np) = 0
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set jsflag to 2 if jflag = -1 for a species. If not, if the
c     species is an active dependent species, set jsflag to 2 if
c     jsflag = -1 for any other species appearing in its reaction.
c
      do ns = 1,nsta
        nr1 = ndrsrd(1,ns)
        nr2 = ndrsrd(2,ns)
        nt = nr2 - nr1 + 1
        if (jflag(ns) .eq. -1) then
          jsflag(ns) = 2
        elseif (ns.lt.narn1a .or. ns.gt.narn2a
     $    .or. jflag(ns) .eq. 30) then
          if (nt .ge. 2) then
            do n = nr1,nr2
              nss = ndrsd(n)
              if (jflag(nss) .eq. -1) then
                jsflag(ns) = 2
                go to 130
              endif
            enddo
          endif
        endif
  130   continue
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Now provide an exception to the above. Scan auxiliary basis
c     species with jflag .le. 1 and ensure that any basis species that
c     they are linked to do not have jsflag = 2. If so, reset jsflag
c     to 1. Do this recursively, to ensure that no species is skipped.
c     The idea here is to keep these strict basis species in the
c     compressed species set produced by EQLIB\cmpdat.f, so that the
c     one-to-one association with a chemical element is preserved.
c
      nn = max(ncta + 1,nrdxsp + 1)
      qheadr = .true.
  140 ncount = 0
      do nb = nn,nbtd
        ns = nbaspd(nb)
        if (jsflag(ns) .le. 1) then
          nr1 = ndrsrd(1,ns)
          nr2 = ndrsrd(2,ns)
          do n = nr1,nr2
            nse = ndrsd(n)
            if (jsflag(nse) .ge. 2) then
              jsflag(nse) = 1
              ncount = ncount + 1
              if (qheadr) then
                write (noutpt,1000)
                write (nttyo,1000)
 1000           format(/' Preserving the following basis species',
     $          ' in the model',/' to maintain linkage to strict',
     $          ' basis species associated',/' one-to-one with',
     $          ' active chemical elements:',/)
              endif
              qheadr = .false.
c
c             Calling sequence substitutions:
c               uspeca(nse) for unam48
c
              call fmspnx(jlen,uspeca(nse),uspn56)
c
              write (noutpt,1010) uspn56(1:jlen)
              write (nttyo,1010) uspn56(1:jlen)
 1010         format('   ',a)
            endif
          enddo
        endif
      enddo
      if (ncount .gt. 0) go to 140
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set jsflag to 2 if there exist no log K data in the current
c     temperature range, unless the species happens to be in the
c     basis set.
c
      do ns = 1,nsta
        if (axlksd(1,ntpr,ns) .ge. 9999999.) then
c
c         Calling sequence substitutions:
c           nbaspd for nbasp
c           nbtd for nbt
c
          nb = nbasis(nbaspd,nbtd,nbtmax,ns)
          if (nb .eq. 0) jsflag(ns) = 2
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Execute the nxmod suppression options.
c
      call supprs(kxmod,jpflag,jsflag,ncmpra,noutpt,npta,nptmax,
     $ nsta,nstmax,nttyo,nxmdmx,nxmod,uphasa,uspeca,uxmod)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Suppress any phase with zero species.
c
      do np = 1,npta
        nt = ncmpra(2,np) - ncmpra(1,np) + 1
        if (nt .le. 0) jpflag(np) = 2
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Exercise the option to hard suppress all solid solutions.
c
      if (iopt(4) .le. 0) then
        do np = 1,npta
          if (uptypa(np)(1:24) .eq. uptsld(1:24)) then
            nr1 = ncmpra(1,np)
            nr2 = ncmpra(2,np)
            nt = nr2 - nr1 + 1
            if (nt .ge. 2) then
              jpflag(np) = 2
              do ns = nr1,nr2
                jsflag(ns) = 2
              enddo
            endif
          endif
        enddo
      endif
c
      end
