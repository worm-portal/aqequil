      subroutine supprs(kxmod,jpflag,jsflag,ncmpra,noutpt,npta,nptmax,
     $ nsta,nstmax,nttyo,nxmdmx,nxmod,uphasa,uspeca,uxmod)
c
c     This subroutine suppresses phases/species as directed by what is
c     on the input file. Here uxmod is the name of the associated
c     species and kxmod = -1. EQLIB/alters.f handles the log K alter
c     function (kxmod = 0, 1, or 2). Suppression of a phase results in
c     suppression of all of its component species.
c
c     The string uxmod is 48 characters in length. This can contain
c     a phase-specific species name in which the first 24 characters
c     contain the species name proper and the second 24 characters
c     contain the phase name. If the second 24 characters are blank,
c     then every species or phase whose name matches what is in the
c     first 24 characters will be suppressed. Note that it is not
c     possible to specify a name in the second 24 character field
c     with nothing in the first such field, as the contents of the
c     uxmod field read from the input file is automatically left-
c     adjusted when its contents are read. To specify a phase name
c     only, put the string 'Phase_name_only:' in the first 24
c     characters, where a species name proper would ordinarily
c     appear, and put the phase name in the second 24 characters.
c     In practice, this should not be necessary, as phase names
c     should be unique from species names, save for the case of
c     a pure phase, for which the species name and the phase
c     name are identical.
c
c     This subroutine is called by:
c
c       EQLIB/flgset.f
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
      integer nptmax,nstmax,nxmdmx
c
      integer noutpt,nttyo
c
      integer kxmod(nxmdmx),jsflag(nstmax),jpflag(nptmax),
     $ ncmpra(2,nptmax)
      integer npta,nsta,nxmod
c
      character(len=24) uphasa(nptmax)
      character(len=48) uspeca(nstmax)
      character(len=48) uxmod(nxmdmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer, parameter :: nlpmax = 25,nlsmax = 25
c
      integer jlen,jlenx,j2,n,nchar,nerr,nhitp,nhitpl,nhits,
     $ nhitsl,nlist,nn,np,nr1,nr2,ns
c
      integer ilnobl
c
      logical qponly
c
      character(len=24), dimension(:), allocatable :: uusupp
      character(len=56), dimension(:), allocatable :: uusups
c
      character(len=56) uspn56,ux56
      character(len=48) unam48,ux48
      character(len=24) ublk24
      character(len=8) ufix,ux8
c
c-----------------------------------------------------------------------
c
      data ublk24 /'                        '/
      data ufix   /'fix     '/
c
c-----------------------------------------------------------------------
c
      nerr = 0
c
      if (nxmod .le. 0) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Allocate arrays for lists of user-suppressed phases and species.
c
      ALLOCATE (uusupp(nlpmax))
      ALLOCATE (uusups(nlsmax))
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      nhitpl = 0
      nhitsl = 0
c
      do n = 1,nxmod
        if (kxmod(n) .ge. 0) go to 150
        unam48 = uxmod(n)
        nchar = 48
c
        qponly = .false.
        if (unam48(1:16) .eq. 'Phase_name_only:') then
          qponly = .true.
          ux48 = unam48(17:48)
          call lejust(ux48)
          j2 = ilnobl(ux48)
          if (j2 .gt. 24) then
            write (noutpt,1010) ux48(1:j2)
            write (nttyo,1010) ux48(1:j2)
 1010       format(/' * Error - (EQLIB/supprs) Have encountered ',
     $      ' an input file directive',/7x,'to suppress the phase',
     $      ' "',a,'".',/7x,'This name exceeds the allowed',
     $      ' 24 characters.')
            nerr = nerr + 1
            go to 150
          endif
          unam48(1:24) = ux48(1:24)
          unam48(25:48) = ublk24(1:24)
        endif
c
        j2 = ilnobl(unam48(1:24))
        call fmspnm(jlen,unam48,uspn56)
c
        if (unam48(25:48) .eq. ublk24(1:24)) nchar = nchar - 24
        if (unam48(1:24) .eq. ublk24(1:24)) nchar = nchar - 24
c
        if (nchar .eq. 48) then
          if (unam48(1:24) .eq. unam48(25:48)) then
c
c           Have specified a pure phase species like Albite (Albite).
c           Set up to suppress this as the pure phase.
c
            unam48(25:48) = ublk24(1:24)
            qponly = .true.
            nchar = 24
          endif
        endif
c
        if (nchar .eq. 0) then
          write (noutpt,1020)
          write (nttyo,1020)
 1020     format(/' * Error - (EQLIB/supprs) Have encountered a blank',
     $    ' uxmod input for an',/7x,'nxmod suppress option.')
          nerr = nerr + 1
          go to 150
        endif
c
        nhitp = 0
c
        if (nchar .eq. 24) then
c
c         Look in the list of phases.
c
          do np = 1,npta
            if (unam48(1:24) .eq. uphasa(np)) then
              nhitp = nhitp + 1
              nhitpl = nhitpl + 1
              if (nhitp .le. nlpmax) uusupp(nhitpl) = uphasa(np)
              if (jpflag(np) .le. 0) jpflag(np) = 1
              nr1 = ncmpra(1,np)
              nr2 = ncmpra(2,np)
              do ns = nr1,nr2
                if (jsflag(ns) .le. 0) jsflag(ns) = 1
              enddo
              if (qponly) go to 150
              go to 110
            endif
          enddo
        endif
c
        if (qponly) then
c
c         Was only looking to suppress a phase, and did not find it.
c
          write (noutpt,1050) unam48(1:j2)
          write (nttyo,1050) unam48(1:j2)
 1050     format(/" * Warning - (EQLIB/supprs) Can't find the phase ",
     $    '"',a,'",',/7x,'which is specified in an nxmod suppress',
     $    ' option. Check to make sure',/7x,'that a phase of this',
     $    ' name appears on the supporting data file.')
          go to 150
        endif
c
  110   continue
c
c       Look in the list of species.
c
        nhits = 0
c
        do ns = 1,nsta
          if (unam48(1:nchar) .eq. uspeca(ns)(1:nchar)) then
            nhits = nhits + 1
            if (jsflag(ns) .lt. 1) jsflag(ns) = 1
            if (uspeca(ns)(1:24).ne.uspeca(ns)(25:48)) then
              call fmspnm(jlenx,uspeca(ns),ux56)
              nhitsl = nhitsl + 1
              if (nhitsl .le. nlsmax) uusups(nhitsl) = ux56
            endif
          endif
        enddo
        if (nhits .gt. 0) go to 150
c
c       The species to be suppressed was not found.
c
        if (nchar.eq.48) then
          write (noutpt,1120) uspn56(1:jlen)
          write (nttyo,1120) uspn56(1:jlen)
 1120     format(/" * Warning - (EQLIB/supprs) Can't find the",
     $    ' species "',a,'"',/7x,'which is specified in an nxmod',
     $    ' suppress option. Check to make sure',/7x,'that a',
     $    ' species of this name appears on the supporting',
     $    ' data file.')
        elseif (nhitp .le. 0) then
c
c         The name could also have referred to a phase, but no
c         such phase was found.
c
          write (noutpt,1130) uspn56(1:jlen)
          write (nttyo,1130) uspn56(1:jlen)
 1130     format(/" * Warning - (EQLIB/supprs) Can't find an",
     $    ' entity "',a,'"',/7x,'which is specified in an nxmod',
     $    ' suppress option. Check to make sure',/7x,'that a',
     $    ' phase or species of this name appears on the supporting',
     $    ' data file.')
        endif
c
c       See if an attempt was made to suppress a fictive fugacity
c       fixing phase.
c
        if (unam48(1:8) .eq. ufix(1:8)) then
          write (noutpt,1150) unam48(1:j2)
          write (nttyo,1150) unam48(1:j2)
 1150     format(/' * Note - (EQLIB/supprs) The phase "',a,
     $    '" specified',/7x,'in an nxmod suppress option is a',
     $    ' fictive fugacity-fixing phase.',/7x,"Such a phase can't",
     $    ' be suppressed.')
        endif
c
  150   continue
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write lists of user-suppressed phases and species.
c
      nlist = min(nhitpl,nlpmax)
      if (nlist .gt. 0) then
        if (nlist .eq. 1) then
          write (noutpt,1200)
          write (nttyo,1200)
 1200     format(/' The following phase has been user-suppressed:',/)
        else
          write (noutpt,1210)
          write (nttyo,1210)
 1210     format(/' The following phases have been user-suppressed:',/)
        endif
c
        do n = 1,nlist
          j2 = ilnobl(uusupp(n))
          write (noutpt,1220) uusupp(n)(1:j2)
          write (nttyo,1220) uusupp(n)(1:j2)
 1220     format(4x,a)
        enddo
c
        if (nhitpl .gt. nlist) then
          nn = nhitpl - nlist
          write (ux8,'(i8)') nn
          call lejust(ux8)
          j2 = ilnobl(ux8)
          write (noutpt,1230) ux8(1:j2)
          write (nttyo,1230) ux8(1:j2)
 1230     format(6x,'plus ',a,' others')
        endif
      endif
c
      nlist = min(nhitsl,nlsmax)
      if (nlist .gt. 0) then
        if (nlist .eq. 1) then
          write (noutpt,1250)
          write (nttyo,1250)
 1250     format(/' The following species has been user-suppressed:',/)
        else
          write (noutpt,1260)
          write (nttyo,1260)
 1260     format(/' The following species have been user-suppressed:',/)
        endif
c
        do n = 1,nlist
          j2 = ilnobl(uusups(n))
          write (noutpt,1220) uusups(n)(1:j2)
          write (nttyo,1220) uusups(n)(1:j2)
        enddo
c
        if (nhitsl .gt. nlist) then
          nn = nhitsl - nlist
          write (ux8,'(i8)') nn
          call lejust(ux8)
          j2 = ilnobl(ux8)
          write (noutpt,1230) ux8(1:j2)
          write (nttyo,1230) ux8(1:j2)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nerr .gt. 0) stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Deallocate arrays for lists of user-suppressed phases and species.
c
      DEALLOCATE (uusupp)
      DEALLOCATE (uusups)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
