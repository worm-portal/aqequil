      subroutine rd3d7(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,
     $ jflagb,jxmod,kxmod,ncompb,ninpts,nodbmx,nopgmx,noprmx,noptmx,
     $ nsq,nsqmax,ntitl,ntitmx,nttyo,nxmdmx,nxmod,nxtb,nxtmax,qend,
     $ qrderr,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,ubasis,
     $ uebal,umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,xbarb,
     $ uxmd24,xlkmod)
c
c     This subroutine reads the EQ3NR input file in menu-style ("D")
c     format for versions 7.0-7.2. It thus encompasses two version
c     levels.
c
c     This subroutine is called by:
c
c       XCON3/xcon3.f
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer iktmax,nodbmx,nopgmx,noprmx,noptmx,nsqmax,ntitmx,nxmdmx,
     $ nxtmax
c
      integer iodb(nodbmx),iopg(nopgmx),iopr(noprmx),iopt(noptmx),
     $ jflagb(nsqmax),jxmod(nxmdmx),kxmod(nxmdmx),ncompb(nxtmax)
c
      integer itermx,ninpts,nsq,ntitl,nttyo,nxmod,nxtb
c
      logical qend,qrderr
c
      character*80 utitl(ntitmx)
      character*24 ubasis(nsqmax),umemb(iktmax,nxtmax),uphas1(nsqmax),
     $ uphas2(nsqmax),usolb(nxtmax),uspecb(nsqmax),uxmd24(nxmdmx)
      character*24 uebal,uredox
c
      real*8 cspb(nsqmax),xbarb(iktmax,nxtmax),xlkmod(nxmdmx)
c
      real*8 fep,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat
c
c-----------------------------------------------------------------------
c
      include 'xcon3/x3op7.h'
c
c-----------------------------------------------------------------------
c
c     Local parameter declarations.
c
c       nfldpa = maximum number of fields per line
c       nlchpa = character length of a line
c
      integer nfldpa,nlchpa
c
      parameter (nfldpa = 8,nlchpa = 80)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,idesc,idescx,iktb,ivar,j,jdesc,j1,j2,j3,n,ncount,
     $ nfldmx,nfldt,nfldtx,nlchmx,nmark,nmarks
      integer ilnobl
c
      logical qrdxcp
c
      character*(nlchpa) ufield(nfldpa)
      character*(nlchpa) uline1,uline2,ulscr
      character*(nlchpa) uheadr,uheadx
c
      character*80 ustr
c
      real*8 var
c
c-----------------------------------------------------------------------
c
      nfldmx = nfldpa
      nlchmx = nlchpa
c
      qrderr = .false.
c
c-----------------------------------------------------------------------
c
c     Title.
c
      qend = .false.
      read (ninpts,1000,end=100,err=990) uline1
 1000 format(a80)
      call parsln(nfldmx,nfldt,nlchmx,ufield,uline1,ulscr)
      ustr = ufield(1)
      if (ustr(1:8) .ne. '--------') then
        j2 = ilnobl(uline1)
        write (nttyo,1010) uline1(1:j2)
 1010   format(/' * Error - (XCON3/rd3d7) The first line of a "D"',
     $  /7x,'format input file must begin with "|--------".',
     $  /7x,'The first line read from the old input file begins with-',
     $  /7x,'"',a,'".')
        go to 990
      endif
      go to 105
c
  100 qend = .true.
      go to 999
c
  105 do 110 n = 1,ntitmx + 1
        read (ninpts,1000,err=990) uline1
        call parsln(nfldmx,nfldt,nlchmx,ufield,uline1,ulscr)
        ustr = ufield(1)
        if (ustr(1:8) .eq. '--------') go to 120
        utitl(n) = ufield(1)
  110 continue
c
      write (nttyo,1015) ntitmx
 1015 format(/' * Error - (XCON3/rd3d7) Have too many lines in the',
     $ /7x,'main title. The code is only dimensioned for ',i4,
     $ /7x,'lines. Reduce the size of the title or increase the',
     $ /7x,'dimensioning parameter ntitpa.')
      go to 990
c
  120 ntitl = n - 1
c
c-----------------------------------------------------------------------
c
c     Temperature and density.
c
      tempc = 0.
      rho = 0.
c
      uheadx = 'Temperature (C)'
      nfldtx = 4
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      tempc = var
      ustr = ufield(4)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      rho = var
c
c-----------------------------------------------------------------------
c
c     Total dissolved salts.
c
      tdspkg = 0.
      tdspl = 0.
c
      uheadx = 'Total dissolved salts'
      nfldtx = 5
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
      nmark = 0
      ncount = 0
c
      do 125 n = 3,5
        j = index(ufield(n),'*')
        if (j .gt. 0) then
          ncount = ncount + 1
          if (ncount .eq. 1) nmark = n
        endif
  125 continue
c
      if (ncount .eq. 0) then
        j2 = ilnobl(uline1)
        write (nttyo,1017) uline1(1:j2)
 1017   format(/' * Warning - (XCON3/rd3d7) None of the options'
     $  /7x,'is marked with an asterisk on the line beginning with',
     $  /7x,'"',a,'".')
      endif
c
      if (ncount .gt. 1) then
        j2 = ilnobl(uline1)
        write (nttyo,1018) uline1(1:j2)
 1018   format(/' * Warning - (XCON3/rd3d7) More than one of the'
     $  /7x,'options is marked with an asterisk on the line',
     $  ' beginning with',
     $  /7x,'"',a,'".')
      endif
c
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
c
      if (nmark .eq. 3) then
        tdspkg = var
      elseif (nmark .eq. 4) then
        tdspl = var
      endif
c
c-----------------------------------------------------------------------
c
c     Electrical balancing species.
c
      uebal = ' '
c
      uheadx = 'Electrical balancing'
      nfldtx = 4
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
      nmark = 0
      ncount = 0
c
      do 127 n = 3,4
        j = index(ufield(n),'*')
        if (j .gt. 0) then
          ncount = ncount + 1
          if (ncount .eq. 1) nmark = n
        endif
  127 continue
c
      if (ncount .gt. 1) write (nttyo,1018) uline1
c
      uebal = ufield(2)
      if (nmark .eq. 3) uebal = 'pick1.'
      if (nmark .eq. 4) uebal = ' '
c
c-----------------------------------------------------------------------
c
c     Basis species and associated constraints. If a solid solution
c     end-member is part of a constraint, read the solid solution
c     composition here. Other solid solution compositions are read
c     in the following block.
c
      uheadx = 'SPECIES'
      nfldtx = 4
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
      nxtb = 0
      uredox = ' '
      fep = 0.
      qrdxcp = .false.
c
      nsq = 0
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
c
  130 ustr = ufield(1)
      if (ustr(1:8) .eq. '--------') go to 180
c
      nfldtx = 4
      if (nfldt .ne. nfldtx) then
        j2 = ilnobl(uline1)
        write (nttyo,1030) nfldt,nfldtx,uline1(1:j2)
 1030   format(/' * Warning - (XCON3/rd3d7) Found ',i2,' fields',
     $  /7x,'where ',i2,' were expected on the line beginning with-',
     $  /7x,'"',a,'".')
      endif
c
      uheadr = ufield(1)
      call locase(uheadr)
c
      ustr = ufield(3)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
c
      ustr = ufield(4)
      call locase(ustr)
c
      if (uheadr(1:8) .eq. 'redox') then
        fep = var
        if (ustr(1:6) .eq. 'logfo2') then
          iopt(1) = 0
        elseif (ustr(1:2) .eq. 'eh') then
          iopt(1) = -1
        elseif (ustr(1:2) .eq. 'pe') then
          iopt(1) = -2
        elseif (ustr(1:12) .eq. 'redox couple') then
          iopt(1) = 1
          qrdxcp = .true.
        else
          j2 = ilnobl(ustr)
          write (nttyo,1040) ustr(1:j2)
 1040     format(/" * Error - (XCON3/rd3d7) Can't identify the",
     $    /7x,'following redox option string- "',a,'". This must',
     $    /7x,'be one of "LogfO2", "Eh", "pe", or "redox couple".')
          go to 990
        endif
c
        nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        go to 130
      endif
c
      nsq = nsq + 1
c
      if (nsq .gt. nsqmax) then
        write (nttyo,1050) nsqmax
 1050   format(/' * Error - (XCON3/rd3d7) Have too many basis',
     $  /7x,'species. The code is only dimensioned for ',i3,
     $  /7x,'such species. Reduce the number of such species',
     $  /7x,'or increase the dimensioning parameter nsqpar.')
        go to 990
      endif
c
      uspecb(nsq) = ufield(1)
      if (qrdxcp) uredox = ufield(1)
      qrdxcp = .false.
      if (uspecb(nsq)(1:5).eq.'o2(g)' .or.
     $    uspecb(nsq)(1:5).eq.'O2(g)') then
        iopt(1) = -3
        uredox = ' '
      endif
      cspb(nsq) = var
c
      do 140 n = -1,njf7pa
        uheadx = ujflg7(n)
        call locase(uheadx)
        if (ustr(1:16) .eq. uheadx(1:16)) then
          jflagb(nsq) = n
          go to 150
        endif
  140 continue
c
      j2 = ilnobl(ustr)
      write (nttyo,1060) ustr(1:j2)
 1060 format(/" * Error - (XCON3/rd3d7) Can't identify the",
     $ /7x,'following jflag option string- "',a,'". This should',
     $ /7x,'be one of the strings defined in the ujflg7 array.')
      go to 990
c
  150 ustr = ufield(2)
      if (jflagb(nsq).ge.17 .and. jflagb(nsq).le.21) then
        uphas1(nsq) = ustr
      else
        ubasis(nsq) = ustr
      endif
c
      if (jflagb(nsq) .eq. 20) then
        nxtb = nxtb + 1
c
        if (nxtb .gt. nxtmax) then
          write (nttyo,1310) nxtmax
 1310     format(/' * Error - (XCON3/rd3d7) Have too many solid',
     $    /7x,'solutions present. The code is only dimensioned',
     $    /7x,'for ',i3,' solid solutions. Reduce the number of such',
     $    /7x,'phases or increase the dimensioning parameter nxtpar.')
          go to 990
        endif
c
        usolb(nxtb) = ustr
        iktb = 0
  160   nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
c
        if (ufield(1)(1:8) .eq. '        ') then
          nfldtx = 4
          if (nfldt .ne. nfldtx) write (nttyo,1030) nfldt,nfldtx,uline1
          iktb = iktb + 1
          if (iktb .eq. 1) uphas2(nsq) = ufield(2)
          umemb(iktb,nxtb) = ufield(2)
          ustr = ufield(3)
          call chreal(nttyo,qrderr,ustr,var)
          if (qrderr) go to 999
          xbarb(iktb,nxtb) = var
          go to 160
        else
          ncompb(nxtb) = iktb
          go to 130
        endif
      endif
c
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      go to 130
c
  180 continue
c
c-----------------------------------------------------------------------
c
c     Mole fractions of solid solutions.
c
      uheadx = 'input solid solutions'
      nfldtx = 1
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
c
      ustr = ufield(1)
      call locase(ustr)
      if (ustr(1:8) .eq. '--------') go to 260
c
      nfldtx = 4
      if (nfldt .ne. nfldtx) write (nttyo,1030) nfldt,nfldtx,uline1
c
      if (ustr(1:4) .eq. 'none') then
        nfldtx = 1
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        ustr = ufield(1)
        call locase(ustr)
        if (ustr(1:8) .ne. '--------') then
          j2 = ilnobl(uline1)
          write (nttyo,1305) uline1(1:j2)
 1305     format(/' * Error - (XCON3/rd3d7) Found the line beginning',
     $    ' with-',/7x,'"',a,'"',
     $    /7x,'where a dashed separator line was expected.')
          go to 990
        endif
        go to 260
      endif
c
  220 nxtb = nxtb + 1
c
      if (nxtb .gt. nxtmax) then
        write (nttyo,1310) nxtmax
        go to 990
      endif
c
      usolb(nxtb) = ufield(1)
      iktb = 0
c
  230 nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      call locase(ustr)
      if (ustr(1:8) .ne. '        ') then
        ncompb(nxtb) = iktb
        if (ustr(1:8) .eq. '--------') go to 260
        go to 220
      endif
      nfldtx = 4
      if (nfldt .ne. nfldtx) write (nttyo,1030) nfldt,nfldtx,uline1
c
      ustr = ufield(1)
      if (ustr(1:8) .ne. '        ') then
        ncompb(nxtb) = iktb
        go to 220
      endif
c
      iktb = iktb + 1
c
      if (iktb .gt. iktmax) then
        j2 = ilnobl(usolb(nxtb))
        write (nttyo,1330) usolb(nxtb)(1:j2),iktmax
 1330   format(/' * Error - (XCON3/rd3d7) Solid solution',
     $  /7x,'"',a,'" has too many end-members present.',
     $  /7x,'This code is only dimensioned for ',i3,' end-members',
     $  /7x,'per solid solution. Reduce the number of end-members',
     $  /7x,'or increase the dimensioning parameter iktpar.')
        go to 990
      endif
c
      umemb(iktb,nxtb) = ufield(2)
      ustr = ufield(3)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      xbarb(iktb,nxtb) = var
c
      go to 230
c
  260 continue
c
c-----------------------------------------------------------------------
c
c     Nxmod options.
c
      uheadx = 'suppressed species'
      nfldtx = 1
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
c
      n = 0
      ustr = ufield(1)
      call locase(ustr)
      if (ustr(1:8) .eq. '--------') go to 570
c
      if (ustr(1:4) .eq. 'none') go to 560
c
      nfldtx = 4
      if (nfldt .ne. nfldtx) write (nttyo,1030) nfldt,nfldtx,uline1
c
  550 n = n + 1
c
      if (n .gt. nxmdmx) then
        write (nttyo,1132) nxmdmx
 1132   format(/' * Error - (XCON3/rd3d7) Have too many nxmod',
     $  /7x,'alter/suppress options. The code is only dimensioned',
     $  /7x,'for ',i3,' such options. Reduce the number of such',
     $  ' options',/7x,'or increase the dimensioning parameter',
     $  ' nxmdpa.')
        go to 990
      endif
c
      uxmd24(n) = ufield(1)
c
      ustr = ufield(2)
      call locase(ustr)
      if (ustr(1:7) .eq. 'aqueous') then
        jxmod(n) = 0
      elseif (ustr(1:7) .eq. 'mineral') then
        jxmod(n) = 1
      elseif (ustr(1:3) .eq. 'gas') then
        jxmod(n) = 2
      elseif (ustr(1:14) .eq. 'solid solution') then
        jxmod(n) = 3
      else
        j2 = ilnobl(ustr)
        write (nttyo,1135) ustr(1:j2)
 1135   format(/" * Error - (XCON3/rd3d7) Can't identify the",
     $  /7x,'following alter/suppress species type string- "',a,'".',
     $  /7x,'This must be one of "aqueous", "mineral", "gas",',
     $  /7x,' or "solid solution".')
        go to 990
      endif
c
      ustr = ufield(3)
      call locase(ustr)
      if (ustr(1:8) .eq. 'suppress') then
        kxmod(n) = -1
      elseif (ustr(1:7) .eq. 'replace') then
        kxmod(n) = 0
      elseif (ustr(1:8) .eq. 'augmentk') then
        kxmod(n) = 1
      elseif (ustr(1:8) .eq. 'augmentg') then
        kxmod(n) = 2
      else
        j2 = ilnobl(ustr)
        write (nttyo,1140) ustr(1:j2)
 1140   format(/" * Error - (XCON3/rd3d7) Can't identify the",
     $  /7x,'following alter/suppress option string- "',a,'".',
     $  /7x,'This must be one of "suppress", "replace", "augmentk",',
     $  /7x,' or "augmentg".')
        go to 990
      endif
c
      ustr = ufield(4)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      xlkmod(n) = var
c
  560 nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      call locase(ustr)
      if (ustr(1:8) .eq. '--------') go to 570
      nfldtx = 4
      if (nfldt .ne. nfldtx) write (nttyo,1030) nfldt,nfldtx,uline1
c
      go to 550
c
  570 nxmod = n
c
c-----------------------------------------------------------------------
c
c     Options. These are iopt, iopg, and iopr option switches,
c     skipping ones which may be classified as development options.
c     The development options are read below, in a somewhat different
c     manner. The iodb option switches are read still further below,
c     in a manner similar to that employed here.
c
c     Note: iopt(1) = iopt1, etc.
c
      nfldtx = 1
      uheadx = 'options'
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
      i = 0
      nfldtx = 1
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      if (ustr(1:8) .eq. '--------') go to 420
c
      ulscr = ufield(1)
      if (ulscr(1:2) .ne. '- ') then
        j = ilnobl(uline1(1:70))
        write (nttyo,1250) uline1(1:j)
 1250   format(/' * Error - (XCON3/rd3d7) Found the line beginning',
     $  ' with-',/7x,'"',a,'"',
     $  /7x,'where an iopt, iopg, iopr, or iodb option switch header',
     $  /7x,'was expected.')
        go to 990
      endif
c
c     Have found an option header.
c
  340 i = i + 1
      nmark = 0
      nmarks = 0
c
c     Put the option header string minus the preceding "- " string
c     into the variable uheadr.
c
      uheadr = ulscr(3:nlchmx)
      call locase(uheadr)
c
c     Remove the succeeding " -" string, if any.
c
      j = ilnobl(uheadr)
      if (uheadr(j - 1:j) .eq. ' -') uheadr(j:j) = ' '
c
c     Identify the corresponding option.
c
      do 370 idescx = 1,nop3pa
        uheadx = uopt3(idescx)
        call locase(uheadx)
        if (uheadr(1:40) .eq. uheadx(1:40)) then
          if (idescx .ne. i) then
            j2 = ilnobl(uopt3(idescx))
            j3 = ilnobl(uvar3(idescx))
            write (nttyo,1260) uopt3(idescx)(1:j2),
     $      uvar3(idescx)(1:j3),index3(idescx),i,idescx
 1260       format(/' * Warning - (XCON3/rd3d7) The input for the',
     $      /7x,'option whose identifying string which begins with-',
     $      /7x,'"',a,'"',/7x,'(',a,'(',i2,')) is out of order,',
     $      ' in place ',i3,' instead of place ',i3,'.')
          endif
          go to 380
        endif
  370 continue
c
      j2 = ilnobl(uheadr)
      write (nttyo,1270) uheadr(1:j2)
 1270 format(/" * Error - (XCON3/rd3d7) Can't identify the following",
     $ /7x,'iopt, iopg, iopr, or iodb option header-',
     $ /7x,'"',a,'".')
      go to 990
c
  380 continue
c
c     Read the first line after an option  header.
c
      jdesc = 0
  390 jdesc = jdesc + 1
      nfldtx = 1
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ulscr = ufield(1)
      if (ulscr(1:8) .eq. '--------') then
        if (nmarks .gt. 1) then
          j2 = ilnobl(uopt3(idescx))
          j3 = ilnobl(uvar3(idescx))
          write (nttyo,1280) uopt3(idescx)(1:j2),
     $    uvar3(idescx)(1:j3),index3(idescx)
 1280     format(/' * Error - (XCON3/rd3d7) More than one choice',
     $    /7x,'was marked for the following option-',
     $    /7x,'"',a,'"',/9x,a,'(',i2,').')
          go to 990
        endif
        go to 420
      endif
c
c     Is the current line a new option header?
c
      if (ulscr(1:2) .eq. '- ') then
        if (nmarks .gt. 1) then
          j2 = ilnobl(uopt3(idescx))
          j3 = ilnobl(uvar3(idescx))
          write (nttyo,1280) uopt3(idescx)(1:j2),
     $    uvar3(idescx)(1:j3),index3(idescx)
          go to 990
        endif
        go to 340
      endif
c
      j = 1
      nmark = 0
      if (ulscr(1:1) .eq. '*') then
        j = 2
        nmark = 1
        nmarks = nmarks + 1
      endif
      uheadr = ulscr(j:nlchmx)
      call locase(uheadr)
      call lejust(uheadr)
c
      do 400 idesc = 1,nod3pa
        if (iopti3(idesc) .eq. idescx) then
          uheadx = udesc3(idesc)
          call locase(uheadx)
          if (uheadr(1:40) .eq. uheadx(1:40)) go to 410
        endif
  400 continue
c
      j = ilnobl(uheadr)
      j2 = ilnobl(uopt3(idescx))
      j3 = ilnobl(uvar3(idescx))
      write (nttyo,1290) uheadr(1:j),uopt3(idescx)(1:j2),
     $ uvar3(idescx)(1:j3),index3(idescx)
 1290 format(/" * Error - (XCON3/rd3d7) Can't identify the following",
     $ /7x,'option string-',/7x,'"',a,'".',
     $ /7x,'It was given for the following option-',
     $ /7x,'"',a,'"',/9x,a,'(',i2,').')
      go to 990
c
  410 if (nmark .eq. 1) then
        ustr = uvar3(idescx)
        call locase(ustr)
        if (ustr(1:4) .eq. 'iopt') then
          iopt(index3(idescx)) = ivalu3(idesc)
        elseif (ustr(1:4) .eq. 'iopg') then
          iopg(index3(idescx)) = ivalu3(idesc)
        elseif (ustr(1:4) .eq. 'iopr') then
          iopr(index3(idescx)) = ivalu3(idesc)
        elseif (ustr(1:4) .eq. 'iodb') then
          iodb(index3(idescx)) = ivalu3(idesc)
        else
          j = ilnobl(uheadr)
          j3 = ilnobl(uvar3(idescx))
          write (nttyo,1300) uvar3(idescx)(1:j3),uheadr(1:j)
 1300     format(/" * Error - (XCON3/rd3d7) Don't recognize the",
     $    /7x,'option type string "',a,'", which was given for',
     $    /7x,'the following option-',/7x,'"',a,'".')
          go to 990
        endif
      endif
c
      go to 390
c
  420 continue
c
c-----------------------------------------------------------------------
c
c     Iodb options.
c
c     Note: iodb(1) = iodb, etc.
c
      nfldtx = 1
      uheadx = 'debugging switches'
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      uheadr = ufield(1)
      call locase(uheadr)
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .eq. 0) then
        j2 = ilnobl(uheadx)
        j3 = ilnobl(uline1)
        write (nttyo,1020) uheadx(1:j2),uline1(1:j3)
 1020   format(/' * Error - (XCON3/rd3d7) Was expecting to find the',
     $  /7x,'header beginning with-',/7x,'"',a,'"',
     $  /7x,'on the line beginning with-',/7x,'"',a,'".')
        go to 990
      endif
c
      i = 0
  440 i = i + 1
      nfldtx = 1
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ulscr = ufield(1)
      if (ulscr(1:8) .eq. '--------') go to 470
c
      j = index(ulscr,' ')
      j2 = j - 1
      if (j2 .eq. 0) j2 = 1
      ustr = ulscr(1:j2)
      j1 = j + 1
      if (j1 .gt. nlchmx) j1 = nlchmx
      uheadr = ulscr(j1:nlchmx)
      call locase(uheadr)
      call lejust(uheadr)
      do 450 idescx = 1,ndb3pa
        uheadx = udebug(idescx)
        call locase(uheadx)
        if (uheadr(1:40) .eq. uheadx(1:40)) go to 460
  450 continue
c
      j2 = ilnobl(uheadr)
      write (nttyo,1110) uheadr(1:j2)
 1110 format(/" * Error - (XCON3/rd3d7) Can't identify the following",
     $ /7x,'iodb option header-',
     $ /7x,'"',a,'".')
      go to 990
c
  460 call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 999
      iodb(idbugi(idescx)) = ivar
      go to 440
c
  470 continue
c
c-----------------------------------------------------------------------
c
c     Development options. There are none.
c
      uheadx = 'development options'
      nfldtx = 1
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      uheadr = ufield(1)
      call locase(uheadr)
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .eq. 0) then
        write (nttyo,1020) uheadx,uline1
        go to 990
      endif
c
      nfldtx = 1
      uheadx = 'none'
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      uheadr = ufield(1)
      call locase(uheadr)
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .eq. 0) then
        write (nttyo,1020) uheadx,uline1
        go to 990
      endif
c
c-----------------------------------------------------------------------
c
c     Tolerances.
c
      uheadx = 'tolerances'
      nfldtx = 1
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      uheadr = ufield(1)
      call locase(uheadr)
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .eq. 0) then
        write (nttyo,1020) uheadx,uline1
        go to 990
      endif
c
      i = 0
  510 i = i + 1
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      if (ustr(1:8) .eq. '--------') go to 520
c
      nfldtx = 3
      if (nfldt .ne. nfldtx) write (nttyo,1030) nfldt,nfldtx,uline1
c
      uheadr = ufield(1)
      call locase(uheadr)
c
      do 515 idescx = 1,nto3pa
        uheadx = utol3(idescx)
        call lejust(uheadx)
        call locase(uheadx)
        if (uheadr(1:32) .eq. uheadx(1:32)) go to 517
  515 continue
c
      j2 = ilnobl(uheadr)
      write (nttyo,1130) uheadr(1:j2)
 1130 format(/" * Error - (XCON3/rd3d7) Can't identify the following",
     $ /7x,'tolerance parameter header-',
     $ /7x,'"',a,'".')
      go to 990
c
  517 ustr = ufield(2)
      if (idescx .eq. 1) then
c
c       Tolbt.
c
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        tolbt = var
      elseif (idescx .eq. 2) then
c
c       Toldl.
c
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        toldl = var
      elseif (idescx .eq. 3) then
c
c       Tolsat.
c
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        tolsat = var
      elseif (idescx .eq. 4) then
c
c       Itermx.
c
        call chrint(ivar,nttyo,qrderr,ustr)
        if (qrderr) go to 999
        itermx = ivar
      endif
c
      go to 510
c
  520 continue
c-----------------------------------------------------------------------
c
      go to 999
c
  990 qrderr = .true.
c
  999 continue
      end
