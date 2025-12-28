      subroutine chkinz(ier,imchmx,imech,iopt,jcode,kmax,kxt,nelect,
     $ noptmx,noutpt,no2gaq,nrct,nrctmx,nrk,nstmax,nttyo,rkb,ureac,
     $ uspeca,uzveci,zvclgi)
c
c     This subroutine checks the code input for various kinds of
c     errors and inconsistencies. Here ier accumulates the
c     number of errors caught.
c
c     This subroutine is called by:
c
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
      integer imchmx,kmax,noptmx,nrctmx,nstmax
c
      integer noutpt,nttyo
c
      integer imech(2,nrctmx),iopt(noptmx),jcode(nrctmx),nrk(2,nrctmx)
c
      integer ier,kxt,nelect,no2gaq,nrct
c
      character*48 uspeca(nstmax),uzveci(kmax)
      character*24 ureac(nrctmx)
c
      real*8 rkb(imchmx,2,nrctmx),zvclgi(kmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,jlen,j2,kcol,nrc
c
      integer ilnobl
c
      character*56 uspn56
c
c-----------------------------------------------------------------------
c
      ier = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Rate law data.
c
c       nrk(1,nrc)= forward rate law code:
c
c         -1 = Use the backward rate law form (legal only if
c              nrk(2,nrc) = 2)
c          0 = Illegal value
c          1 = Relative rate
c          2 = Transition state theory net rate
c          3 = Linear rate
c
c       For the case nrk(1,nrc) = 2, the reactant must be
c       either a pure mineral (jcode = 0) or a solid solution
c       (jcode = 1).
c
c       nrk(2,nrc)= backward rate law code:
c
c         -1 = Use the forward rate law form (legal only if
c              nrk(1,nrc) = 2)
c          0 = No rate law specified; the reaction may be controlled
c              by partial equilibrium.
c          1 = Relative rate
c          2 = Transition state theory net rate
c          3 = Linear rate
c
c       For the case nrk(2,nrc) = 2, the reactant must be
c       either a pure mineral (jcode = 0) or a solid solution
c       (jcode = 1)
c
      do nrc = 1,nrct
        j2 = ilnobl(ureac(nrc))
c
        if (iopt(2) .le. 0) then
c
c         Check for rate law choices that can't be used without
c         an explicit time frame.
c
          if (nrk(2,nrc) .ge. 2) then
            write (noutpt,1000) ureac(nrc)(1:j2),nrk(2,nrc)
            write (nttyo,1000) ureac(nrc)(1:j2),nrk(2,nrc)
 1000       format(/' * Error - (EQ6/chkinz) The forward rate law code',
     $      /7x,'for ',a,' is ',i2,', which is not valid',
     $      /7x,'unless iopt(2) = 1 (the model has a time frame).')
            ier = ier + 1
          endif
c
          if (nrk(2,nrc) .ge. 2) then
            write (noutpt,1010) ureac(nrc)(1:j2),nrk(2,nrc)
            write (nttyo,1010) ureac(nrc)(1:j2),nrk(2,nrc)
 1010       format(/' * Error - (EQ6/chkinz) The backward rate law',
     $      ' code',/7x,'for ',a,' is ',i2,', which is not valid',
     $      /7x,'unless iopt(2) = 1 (the model has a time frame).')
            ier = ier + 1
          endif
        endif
c
c       Check for invalid forward and backward rate law combinations.
c       The following combinations are invalid:
c
c         nrk(1,nrc) = 1 and nrk(2,nrc) = 1
c         nrk(1,nrc) = 3 and nrk(2,nrc) = 3
c         nrk(1,nrc) = -1 and nrk(2,nrc) is not 2
c         nrk(2,nrc) = 0
c
c
        if (nrk(1,nrc).eq.1 .and. nrk(2,nrc).eq.1) then
          write (noutpt,1020) ureac(nrc)(1:j2)
          write (nttyo,1020) ureac(nrc)(1:j2)
 1020     format(' * Error - (EQ6/chkinz) The forward and backward',
     $    /7x,'rate law codes for ',a,' are both 1, which is not a',
     $    /7x,'valid combination.')
          ier = ier + 1
        endif
c
        if (nrk(1,nrc).eq.3 .and. nrk(2,nrc).eq.3) then
          write (noutpt, 1030) ureac(nrc)(1:j2)
          write (nttyo, 1030) ureac(nrc)(1:j2)
 1030     format(' * Error - (EQ6/chkinz) The forward and backward',
     $    /7x,'rate law codes for ',a,' are both 3, which is not a',
     $    /7x,'valid combination.')
          ier = ier + 1
        endif
c
        if (nrk(1,nrc) .eq. -1) then
          if (nrk(2,nrc) .ne. 2) then
            write (noutpt, 1040) ureac(nrc)(1:j2)
            write (nttyo, 1040) ureac(nrc)(1:j2)
 1040       format(' * Error - (EQ6/chkinz) The forward rate law code',
     $      /7x,'for ',a,' is -1, which is not a valid value  unless',
     $      /7x,'the backward rate law code is 2.')
            ier = ier + 1
          endif
        endif
c
        if (nrk(1,nrc) .eq. 0) then
          write (noutpt,1050) ureac(nrc)(1:j2)
          write (nttyo,1050) ureac(nrc)(1:j2)
 1050     format(' * Error - (EQ6/chkinz) The forward rate law code',
     $    /7x,'for ',a,' is 0, which is not a valid value.')
          ier = ier + 1
        endif
c
c       Check the parameters for the forward rate laws.
c
        if (nrk(1,nrc) .eq. 2) then
c
c         Transition state theory.
c
          do i = 1,imech(1,nrc)
            if (rkb(i,1,nrc) .lt. 0.) then
              write (noutpt,1060) i,rkb(i,1,nrc),ureac(nrc)(1:j2),
     $        nrk(1,nrc)
              write (nttyo,1060) i,rkb(i,1,nrc),ureac(nrc)(1:j2),
     $        nrk(1,nrc)
 1060         format (/' * Error - (EQ6/chkinz) Forward rate constant',
     $        /7x,i2,' has a value of ',1pe11.4,' for reactant'
     $        /7x,a,'. A negative rate constant is not valid for a',
     $        /7x,'forward rate code of ',i2,'.')
              ier = ier + 1
            endif
          enddo
c
          if (jcode(nrc) .gt. 1) then
            write (noutpt,1070) nrk(1,nrc),ureac(nrc)(1:j2)
            write (nttyo,1070) nrk(1,nrc),ureac(nrc)(1:j2)
 1070       format(/' * Error - (EQ6/chkinz) Have specified a',
     $      /7x,'forward rate law code of ',i2,' for reactant',
     $      /7x,a,'. This rate law incorporates an affinity',
     $      /7x,'dependence and can only be used for a reactant which',
     $      /7x,'is a mineral or solid solution.')
            ier = ier + 1
          endif
        elseif (nrk(1,nrc) .eq. 3) then
c
c         Linear rate law
c
          if (rkb(1,1,nrc) .lt. 0.) then
            write (noutpt,1060) i,rkb(1,1,nrc),ureac(nrc)(1:j2),
     $      nrk(1,nrc)
            write (nttyo,1060) i,rkb(1,1,nrc),ureac(nrc)(1:j2),
     $      nrk(1,nrc)
            ier = ier + 1
          endif
        endif
c
c       Check the parameters for the backward rate laws.
c
        if (nrk(2,nrc) .eq. -1) then
c
c         Use the forward rate law for the backward case.
c
          if (nrk(1,nrc) .ne. 2) then
            write (noutpt,1080) ureac(nrc)(1:j2),nrk(1,nrc)
            write (nttyo,1080) ureac(nrc)(1:j2),nrk(1,nrc)
 1080       format(' * Error - (EQ6/chkinz) The backward rate law',
     $      /7x,'code for ',a,' is  -1, which is invalid because',
     $      /7x,'the forward rate law code is ',i2,'.')
            ier = ier + 1
          endif
        elseif (nrk(2,nrc) .eq. 2) then
c
c         Transition state rate law.
c
          do i = 1,imech(2,nrc)
            if (rkb(i,2,nrc) .lt. 0.) then
              write (noutpt,1090) i,rkb(i,2,nrc),ureac(nrc)(1:j2),
     $        nrk(2,nrc)
              write (nttyo,1090) i,rkb(i,2,nrc),ureac(nrc)(1:j2),
     $        nrk(2,nrc)
 1090         format (/' * Error - (EQ6/chkinz) Backward rate constant',
     $        /7x,i2,' has a value of ',1pe11.4,' for reactant'
     $        /7x,a,'. A negative rate constant is not valid for a',
     $        /7x,'forward rate code of ',i2,'.')
              write (noutpt,1060) ureac(nrc)(1:j2)
              write (nttyo,1060) ureac(nrc)(1:j2)
              ier = ier + 1
            endif
          enddo
c
          if (jcode(nrc) .gt. 1) then
            write (noutpt,1100) nrk(2,nrc),ureac(nrc)(1:j2)
            write (nttyo,1100) nrk(2,nrc),ureac(nrc)(1:j2)
 1100       format(/' * Error - (EQ6/chkinz) Have specified a',
     $      /7x,'backward rate law code of ',i2,' for reactant',
     $      /7x,a,'. This rate law incorporates an affinity',
     $      /7x,'dependence and can only be used for a reactant which',
     $      /7x,'is a mineral or solid solution.')
            write (noutpt,1070)
            write (nttyo,1070)
            ier = ier + 1
          endif
c
        elseif (nrk(2,nrc) .eq. 3) then
c
c         Linear  rate law.
c
          if (rkb(1,2,nrc) .lt. 0.) then
            write (noutpt,1090) i,rkb(1,2,nrc),ureac(nrc)(1:j2),
     $      nrk(2,nrc)
            write (nttyo,1090) i,rkb(1,2,nrc),ureac(nrc)(1:j2),
     $      nrk(2,nrc)
            ier = ier + 1
          endif
c
          if (jcode(nrc) .gt. 1) then
            if (nrk(1,nrc).eq.3 .or. nrk(1,nrc).eq.1) then
              write (noutpt,1110) ureac(nrc)(1:j2),nrk(1,nrc)
              write (nttyo,1110) ureac(nrc)(1:j2),nrk(1,nrc)
 1110         format(' * Error - (EQ6/chkinz) The forward rate law',
     $        /7x,'code for reactant ',a,' may not be 3 if the',
     $        /7x,'backward rate law code is',i2,'. This is permitted',
     $        /7x,'only for a reactant which is a pure mineral or',
     $        /7x,'a solid solution..')
              ier = ier + 1
            endif
          endif
c
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Matrix variable values.
c
c        zvclgi = logarithmic basis variable, starting estimate
c           or pick-up value
c
      do kcol = 1,kxt
        if (zvclgi(kcol) .le. -99999.) then
          if (uzveci(kcol)(1:8).ne.uspeca(no2gaq)(1:8) .and.
     $      uzveci(kcol)(1:8).ne.uspeca(nelect)(1:8)) then
            call fmspnx(jlen,uzveci(kcol),uspn56)
            write (noutpt,1200) uspn56(1:jlen)
            write (nttyo,1200) uspn56(1:jlen)
 1200       format(/' * Error - (EQ6/chkinz) The basis species ',a,
     $      /7x,'has a log number of moles (zvclgi) value less than',
     $      ' or equal to -99999.',/7x,'on the input file.')
            ier = ier + 1
          endif
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
