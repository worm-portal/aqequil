      subroutine rd3d8b(cgexj,cgexpi,covali,ehi,egexsi,fo2lgi,iebal3,
     $ ietmax,iodb,iopg,iopr,iopt,irdxc3,itdsf3,itermx,jgext,jetmax,
     $ jflgi,jgexti,jpres3,kxmod,mwtges,nbti,nbtmax,ncmpri,net,neti,
     $ netmax,ngexti,ninpts,ngexrt,nobswt,nodbmx,nopgmx,noprmx,
     $ noptmx,noutpt,nprob,nsbswt,ntitl,ntitmx,nttyo,nxicmx,nxmdmx,
     $ nxmod,nxti,nxtimx,pei,press,qend,qrderr,rho,scamas,tdspkg,
     $ tdspl,tempc,tgexp,tolbt,toldl,tolspf,ucospi,uebal,ugexj,
     $ ugexji,ugexmo,ugexp,ugexpi,ugexr,ugexsi,umemi,uobsw,uredox,
     $ usbsw,usoli,uspeci,utitl,uhfgex,uvfgex,uxkgex,uxmod,xbari,
     $ xgexsi,xhfgex,xlkgex,xvfgex,xlkmod,zgexj)
c
c     This subroutine reads the EQ3NR input file in menu-style ("D")
c     format for version 8.0.
c
c     This subroutine is a near-clone of EQ3NR/rd3ind.f. However, the
c     present subroutine embodies only a pure read function (it does
c     only minimal checking of what is read, to ensure that what
c     follows is readable). EQ3NR/rd3ind.f differs in that it also
c     writes an instant echo of what is read to the EQ3NR output file.
c
c     The calling sequence of this subroutine is identical to that of
c     EQ3NR/rd3inw.f, EQ3NR/rd3ind.f, and XCON3/rd3w8.f.
c
c     This subroutine is called by:
c
c       XCON3/xcon3.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ninpts = unit number of the stripped input file
c       noutpt = unit number of the output file
c       nttyo  = unit number of the screen file
c
c     Principal output:
c
c       qrderr = flag denoting a read error
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ietmax,jetmax,nbtmax,netmax,nodbmx,nopgmx,noprmx,noptmx,
     $ ntitmx,nxmdmx,nxicmx,nxtimx
c
      integer ninpts,noutpt,nttyo
c
      integer iodb(nodbmx),iopg(nopgmx),iopr(noprmx),iopt(noptmx),
     $ jgext(netmax),jgexti(netmax),jflgi(nbtmax),kxmod(nxmdmx),
     $ ncmpri(2,nxtimx),ngexrt(jetmax,netmax),ngexti(jetmax,netmax)
c
      integer iebal3,irdxc3,itdsf3,itermx,jpres3,nbti,net,neti,nobswt,
     $ nprob,nsbswt,ntitl,nxmod,nxti
c
      logical qend,qrderr
c
      character*80 utitl(ntitmx)
      character*56 ugexr(ietmax,jetmax,netmax)
      character*48 ucospi(nbtmax),uobsw(2,nbtmax),usbsw(2,nbtmax),
     $ uspeci(nbtmax),uxmod(nxmdmx)
      character*24 ugexmo(netmax),ugexp(netmax),ugexpi(netmax),
     $ ugexsi(ietmax,jetmax,netmax),umemi(nxicmx),usoli(nxtimx)
      character*24 uebal,uredox
      character*8 ugexj(jetmax,netmax),ugexji(jetmax,netmax),
     $ uhfgex(ietmax,jetmax,netmax),uvfgex(ietmax,jetmax,netmax),
     $ uxkgex(ietmax,jetmax,netmax)
c
      real*8 cgexj(jetmax,netmax),cgexpi(netmax),covali(nbtmax),
     $ egexsi(ietmax,jetmax,netmax),mwtges(netmax),tgexp(netmax),
     $ xbari(nxicmx),xhfgex(ietmax,jetmax,netmax),
     $ xlkgex(ietmax,jetmax,netmax),xvfgex(ietmax,jetmax,netmax),
     $ xlkmod(nxmdmx),zgexj(jetmax,netmax),xgexsi(ietmax,jetmax,netmax)
c
      real*8 ehi,fo2lgi,pei,press,rho,scamas,tdspkg,tdspl,tempc,tolbt,
     $ toldl,tolspf
c
c-----------------------------------------------------------------------
c
      include 'eqlib/eqlj8.h'
      include 'eqlib/eqlo8.h'
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
      integer icount,iei,ii,ival,ivar,j,je,jei,jj,jlast,j2,j3,j4,j5,
     $ k1,k2,k3,k4,kxmd,n,nbi,ne,nei,nfldmx,nfldt,nfldtx,nlchmx,nmark,
     $ nn,nnn,nobsw,nsbsw,nxi,nxic
c
      integer ilnobl
c
      logical qgexef,qmark,qnone,qnonei,qnonep,qnones,qnoner,
     $ qokay,qstop
c
      character*(nlchpa) ufield(nfldpa)
      character*(nlchpa) uline1,uline2,ulscr
      character*(nlchpa) uheadr,uheadx
c
      character*80 ustr
      character*48 ux48
      character*24 ustr24,ustrn
      character*1 ux1
c
      real*8 var
c
c-----------------------------------------------------------------------
c
c     The following is a bit of nonsense so compiler warnings will
c     not be generated saying that nprob is not used.
c
      n = nprob
      nprob = n
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      nfldmx = nfldpa
      nlchmx = nlchpa
c
      qrderr = .false.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check some dimensioning parameters.
c
      qstop = .false.
c
      if (ndbxpa .ne. nodbmx) then
        write (nttyo,3000) ndbxpa,nodbmx
 3000   format(/' * Error - (XCON3/rd3d8b) The dimensioning parameter',
     $  ' for the',/7x,'number of iodb debugging print option switches',
     $  ' with string definitions',/7x,'(ndbxpa) has a value of ',
     $  i3,', but the dimensioning',/7x,'parameter for the number',
     $  ' of such switches (nodbpa) has a',/7x,'value of ',i3,'.')
        qstop = .true.
      endif
c
      if (npgxpa .ne. nopgmx) then
        write (nttyo,3010) npgxpa,nopgmx
 3010   format(/' * Error - (XCON3/rd3d8b) The dimensioning parameter',
     $  ' for the',/7x,'number of iopg activity coefficient option',
     $  ' switches with string definitions',/7x,'(npgxpa) has a value',
     $  ' of ' ,i3,', but the dimensioning',/7x,'parameter for the',
     $  ' number of such switches (nopgpa) has a',/7x,'value of ',i3,
     $  '.')
        qstop = .true.
      endif
c
      if (nprxpa .ne. noprmx) then
        write (nttyo,3020) nprxpa,noprmx
 3020   format(/' * Error - (XCON3/rd3d8b) The dimensioning parameter',
     $  ' for the',/7x,'number of iopr print option switches',
     $  ' with string definitions',/7x,'(nprxpa) has a value of ',
     $  i3,', but the dimensioning',/7x,'parameter for the number',
     $  ' of such switches (noprpa) has a',/7x,'value of ',i3,'.')
        qstop = .true.
      endif
c
      if (nptxpa .ne. noptmx) then
        write (nttyo,3030) nptxpa,noptmx
 3030   format(/' * Error - (XCON3/rd3d8b) The dimensioning parameter',
     $  ' for the',/7x,'number of iopt model option switches',
     $  ' with string definitions',/7x,'(nptxpa) has a value of ',
     $  i3,', but the dimensioning',/7x,'parameter for the number',
     $  ' of such switches (noptpa) has a',/7x,'value of ',i3,'.')
        qstop = .true.
      endif
c
      if (qstop) stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Title.
c
      qend = .false.
c
c     Read the first line of the file. This must be a separator line
c     ("|----- ... -----|").
c
      read (ninpts,1000,end=100,err=990) uline1
 1000 format(a80)
      if (uline1(1:8) .ne. '|-------') then
        j2 = ilnobl(uline1)
        j2 = min(j2,70)
        write (nttyo,1010) uline1(1:j2)
 1010   format(/' * Error - (XCON3/rd3d8b) The first line of a "D"',
     $  ' format input file',/7x,'should be a separator line;',
     $  ' therefore, it should begin with',/7x,'"|-------". The first',
     $  ' line begins instead with',/7x,'"',a,'".')
        go to 990
      endif
      go to 105
c
  100 qend = .true.
      go to 999
c
  105 continue
c
c     Read the block title ("Title") from a two-line header.
c
      uheadx = 'Title'
      nfldtx = 2
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c     Now read the title itself.
c
      n = 0
      do nn = 1,ntitmx + 1
        read (ninpts,1000,err=990) uline1
        call parsln(nfldmx,nfldt,nlchmx,ufield,uline1,ulscr)
        ustr = ufield(1)
c
c       A separator line terminates the this block. It is not part
c       of the title itself.
c
        if (ustr(1:8) .eq. '--------') go to 120
c
        n = n + 1
c
        if (n .gt. ntitmx) then
          write (nttyo,1015) ntitmx
 1015     format(/' * Error - (XCON3/rd3d8b) Have too many lines in',
     $    /7x,'the title. The code is only dimensioned for ',i4,
     $    /7x,'lines. Reduce the size of the title or increase the',
     $    /7x,'dimensioning parameter ntitpa.')
          go to 990
        endif
c
        utitl(n) = ufield(1)
      enddo
c
  120 ntitl = n
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Special basis switches.
c
c     Read the block title from a two-line header.
c
      uheadx = 'Special Basis Switches'
      nfldtx = 2
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
      nsbsw = 0
      do n = 1,nbtmax + 1
c
c       Read a line. If the block has not been completely read,
c       this contains the name of a species to "Replace", and a
c       sub-block for that species follows. Otherwise, this line is
c       the first line of the next block.
c
        nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        ustr = ufield(1)
        uheadx = 'Replace'
        call locase(ustr)
        call locase(uheadx)
        j2 = ilnobl(ustr)
        j3 = ilnobl(uheadx)
c
        if (ustr(1:j2) .ne. uheadx(1:j3)) then
c
c         Back up.
c
          backspace ninpts
          go to 125
        endif
c
        ustr = ufield(2)
        ustrn = ustr
        call locase(ustrn)
        if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') ustr = 'None'
        qnone = ustr(1:5).eq.'None '
c
        if (.not.qnone) then
          nsbsw = nsbsw + 1
          usbsw(1,nsbsw) = ufield(2)
        endif
c
c       Read the name of the "with" species from a two-line header.
c
        uheadx = 'with'
        nfldtx = 3
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
        if (qrderr) go to 999
        if (.not.qnone) usbsw(2,nsbsw) = ufield(2)
      enddo
  125 nsbswt = nsbsw
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Temperature.
c
      tempc = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Temperature (C)'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      tempc = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Pressure.
c
      jpres3 = 0
      press = 0.
      icount = 0
c
c     Read a one-line header.
c
      uheadx = 'Pressure option (jpres3):'
      nfldtx = 1
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
c
c     Read the first option from a data line.
c
      nfldtx = 1
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      ustr24 = ustr(5:29)
      uheadx = '( 0) Data file reference'
      call locase(ustr24)
      call locase(uheadx)
      j2 = ilnobl(ustr24)
      j3 = ilnobl(uheadx)
      if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
 1020   format(/' * Error - (XCON3/rd3d8b) Was expecting to find a',
     $  ' line beginning with',/7x,'"',a,'", instead found one',
     $  ' beginning with',/7x,'"',a,'".')
        qrderr = .true.
        go to 999
      endif
      if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (nttyo,1022) ustr(1:j2)
 1022   format(/' * Error - (XCON3/rd3d8b) Was expecting to find an',
     $  /7x,'option check box "[ ]" in the line beginning with',
     $  /7x,'"',a,'".')
        qrderr = .true.
        go to 999
      endif
      ux1 = ustr(2:2)
      if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        jpres3 = 0
        press = 0.
        icount = icount + 1
      endif
c
c     Read the second option from a data line.
c
      nfldtx = 1
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      ustr24 = ustr(5:29)
      uheadx = '( 1) 1.013-bar/steam-sat'
      call locase(ustr24)
      call locase(uheadx)
      j2 = ilnobl(ustr24)
      j3 = ilnobl(uheadx)
      if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
      endif
      if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
      endif
      ux1 = ustr(2:2)
      if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        jpres3 = 1
        press = 0.
        icount = icount + 1
      endif
c
c     Read the third (last) option from a two-line combination
c     (a data line plus a separator line).
c
      nfldtx = 3
      call rdd2l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      ustr24 = ustr(5:29)
      uheadx = '( 2) Value (bars)'
      call locase(ustr24)
      call locase(uheadx)
      j2 = ilnobl(ustr24)
      j3 = ilnobl(uheadx)
      if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        go to 999
      endif
      if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
      endif
      ux1 = ustr(2:2)
      if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        jpres3 = 2
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        press = var
        icount = icount + 1
      endif
c
      if (icount .eq. 0) then
        write (nttyo,1024)
 1024   format(/' * Note - (XCON3/rd3d8b) No option was selected for',
     $  ' the pressure.',/7x,'The pressure will be set to be in',
     $  ' accord with the',/7x,'data file reference pressure curve.')
        jpres3 = 0
        press = 0.
      elseif (icount .gt. 1) then
        write (nttyo,1026)
 1026   format(/' * Warning - (XCON3/rd3d8b) Multiple options were',
     $  ' selected for',/7x,'the pressure. The pressure will be set',
     $  ' to be in accord with the',/7x,'data file reference pressure',
     $  ' curve.')
        jpres3 = 0
        press = 0.
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Density.
c
      rho = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Density'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      rho = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Total dissolved solutes.
c
      itdsf3 = 0
      tdspkg = 0.
      tdspl = 0.
      icount = 0
c
c     Read the data from a two-line header.
c
      uheadx = 'Total dissolved solutes option (itdsf3):'
      nfldtx = 1
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
c
c     Read the first option from a data line.
c
      nfldtx = 3
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      ustr24 = ustr(5:29)
      uheadx = '( 0) Value (mg/kg.sol)'
      call locase(ustr24)
      call locase(uheadx)
      j2 = ilnobl(ustr24)
      j3 = ilnobl(uheadx)
      if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
      endif
      if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
      endif
      ux1 = ustr(2:2)
      if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        itdsf3 = 0
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        tdspkg = var
        icount = icount + 1
      endif
c
c     Read the second (last) option from a two-line combination
c     (a data line plus a separator line).
c
      nfldtx = 3
      call rdd2l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      ustr24 = ustr(5:29)
      uheadx = '( 1) Value (mg/L)'
      call locase(ustr24)
      call locase(uheadx)
      j2 = ilnobl(ustr24)
      j3 = ilnobl(uheadx)
      if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        go to 999
      endif
      if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
      endif
      ux1 = ustr(2:2)
      if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        itdsf3 = 1
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        tdspl = var
        icount = icount + 1
      endif
c
      if (icount .eq. 0) then
        write (nttyo,1034)
 1034   format(/' * Note - (XCON3/rd3d8b) No option was selected',
     $  ' for the total',/7x,'dissolved solutes (TDS). The TDS will',
     $  ' be set to 0 mg/kg.sol.')
        itdsf3 = 0
        tdspl = 0.
        tdspkg = 0.
      elseif (icount .gt. 1) then
        write (nttyo,1036) tdspkg
 1036   format(/' * Warning - (XCON3/rd3d8b) Multiple options were',
     $  ' selected for the',/7x,'total dissolved solutes (TDS).',
     $  ' The TDS will be set to',/7x,g12.5,' mg/kg.sol.')
        itdsf3 = 0
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Electrical balancing species.
c
      iebal3 = 0
      uebal = ' '
      icount = 0
c
c     Read a one-line header.
c
      uheadx = 'Electrical balancing option (iebal3):'
      nfldtx = 1
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
c
c     Read the first option from a data line.
c
      nfldtx = 1
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      ustr24 = ustr(5:29)
      uheadx = '( 0) No balancing is don'
      call locase(ustr24)
      call locase(uheadx)
      j2 = ilnobl(ustr24)
      j3 = ilnobl(uheadx)
      if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
      endif
      if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
      endif
      ux1 = ustr(2:2)
      if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        iebal3 = 0
        uebal = 'None'
        icount = icount + 1
      endif
c
c     Read the second (last) option from a two-line combination
c     (a data line plus a separator line).
c
      nfldtx = 3
      call rdd2l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      ustr24 = ustr(5:29)
      uheadx = '( 1) Balance on species'
      call locase(ustr24)
      call locase(uheadx)
      j2 = ilnobl(ustr24)
      j3 = ilnobl(uheadx)
      if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
      endif
      if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
      endif
      ux1 = ustr(2:2)
      if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        iebal3 = 1
        ustr = ufield(2)
        uebal = ustr
        icount = icount + 1
      endif
c
      if (icount .eq. 0) then
        write (nttyo,1044)
 1044   format(/' * Note - (XCON3/rd3d8b) No electrical balancing',
     $  ' option was',/7x,'selected. No balancing will be done.')
        iebal3 = 0
        uebal = 'None'
      elseif (icount .gt. 1) then
        write (nttyo,1046)
 1046   format(/' * Warning - (XCON3/rd3d8b) Multiple electrical',
     $  ' balancing options.',/7x,'were selected. No balancing',
     $  ' will be done.')
        iebal3 = 0
        uebal = 'None'
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Default redox constraint.
c
      irdxc3 = 0
      fo2lgi = 0.
      ehi = 0.
      pei = 0.
      uredox = ' '
      icount = 0
c
c     Read a one-line header.
c
      uheadx = 'Default redox constraint (irdxc3):'
      nfldtx = 1
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
c
c     Read the first option from a data line.
c
      nfldtx = 1
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      ustr24 = ustr(5:29)
      uheadx = '(-3) Use O2(g) line in t'
      call locase(ustr24)
      call locase(uheadx)
      j2 = ilnobl(ustr24)
      j3 = ilnobl(uheadx)
      if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
      endif
      if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
      endif
      ux1 = ustr(2:2)
      if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        irdxc3 = -3
        icount = icount + 1
      endif
c
c     Read the second option from a data line.
c
      nfldtx = 3
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      ustr24 = ustr(5:29)
      uheadx = '(-2) pe (pe units)'
      call locase(ustr24)
      call locase(uheadx)
      j2 = ilnobl(ustr24)
      j3 = ilnobl(uheadx)
      if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
      endif
      if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
      endif
      ux1 = ustr(2:2)
      if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        irdxc3 = -2
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        pei = var
        icount = icount + 1
      endif
c
c     Read the third option from a data line.
c
      nfldtx = 3
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      ustr24 = ustr(5:29)
      uheadx = '(-1) Eh (volts)'
      call locase(ustr24)
      call locase(uheadx)
      j2 = ilnobl(ustr24)
      j3 = ilnobl(uheadx)
      if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
      endif
      if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
      endif
      ux1 = ustr(2:2)
      if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        irdxc3 = -1
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        ehi = var
        icount = icount + 1
      endif
c
c     Read the fourth option from a data line.
c
      nfldtx = 3
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      ustr24 = ustr(5:29)
      uheadx = '( 0) Log fO2 (log bars)'
      call locase(ustr24)
      call locase(uheadx)
      j2 = ilnobl(ustr24)
      j3 = ilnobl(uheadx)
      if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
      endif
      if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
      endif
      ux1 = ustr(2:2)
      if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        irdxc3 = 0
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        fo2lgi = var
        icount = icount + 1
      endif
c
c     Read the fifth (last) option from a two-line combination
c     (a data line plus a separator line).
c
      nfldtx = 3
      call rdd2l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      ustr24 = ustr(5:29)
      uheadx = '( 1) Couple (aux. sp.)'
      call locase(ustr24)
      call locase(uheadx)
      j2 = ilnobl(ustr24)
      j3 = ilnobl(uheadx)
      if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
        qrderr = .true.
        go to 999
      endif
      if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (nttyo,1022) ustr(1:j2)
        qrderr = .true.
        go to 999
      endif
      ux1 = ustr(2:2)
      if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        irdxc3 = 1
        ustr = ufield(2)
        uredox = ustr
        icount = icount + 1
      endif
c
      if (icount .eq. 0) then
        write (nttyo,1047) fo2lgi
 1047   format(/' * Note - (XCON3/rd3d8b) No default redox constraint',
     $  ' option was',/7x,'selected. A default condition of log fO2= ',
     $  g12.5,' (log bars) will be used.')
        irdxc3 = 0
      elseif (icount .gt. 1) then
        write (nttyo,1048)
 1048   format(/' * Warning - (XCON3/rd3d8b) Multiple default redox',
     $  ' constraint options.',/7x,'were selected. The last option',
     $  ' that was selected will be used.')
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Aqueous basis species and associated constraints.
c
c     Read the first part of the block title from a one-line header.
c
      uheadx = 'Aqueous Basis Species/Constraint Species'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
c
c     Read the second part of the block title from a two-line header.
c
      uheadx = '(uspeci(n)/ucospi(n))'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
      nbi = 0
c
c     Read the first data line.
c
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
c
c     The label below marks a return point for processing subsequent
c     lines in the current data block.
c
  130 ustr = ufield(1)
c
c     There are no data remaining in the current block if a
c     separator line has been encountered.
c
      if (ustr(1:8) .eq. '--------') go to 180
c
      nfldtx = 3
      if (nfldt .ne. nfldtx) then
        j2 = ilnobl(uline1)
        j2 = min(j2,70)
        write (nttyo,1030) nfldt,nfldtx,uline1(1:j2)
 1030   format(/' * Warning - (XCON3/rd3d8b) Found ',i2,' fields',
     $  ' where ',i2,/7x,'were expected on the line which begins with',
     $  /7x,'"',a,'".')
      endif
      ux48 = ufield(1)
c
      nbi = nbi + 1
c
      if (nbi .gt. nbtmax) then
        j2 = ilnobl(ux48)
        write (nttyo,1050) nbtmax,ux48(1:j2)
 1050   format(/' * Error - (XCON3/rd3d8b) The number of basis',
     $  ' species read',/7x,'from the input file exceeded the',
     $  ' maximum of ',i3,' while',/7x,'trying to read data for',
     $  ' the species ',a,'.',/7x,'Increase the dimensioning',
     $  ' parameter nbtpar.')
        go to 990
      endif
c
      uspeci(nbi) = ux48
c
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      covali(nbi) = var
c
      ustr = ufield(3)
      call locase(ustr)
c
      do n = -1,njfxpa
        uheadx = ujf3(n)
        call locase(uheadx)
        if (ustr(1:16) .eq. uheadx(1:16)) then
          jflgi(nbi) = n
          go to 150
        endif
      enddo
c
      j2 = ilnobl(ustr)
      write (nttyo,1060) ustr(1:j2)
 1060 format(/" * Error - (XCON3/rd3d8b) Don't recognize the",
     $ ' jflag option string',/7x,' "',a,'". This should',
     $ ' be one of the strings defined in the',/7x,'ujf3 array.')
      go to 990
c
  150 if (jflgi(nbi).eq.17 .or. jflgi(nbi).eq.18 .or.
     $ jflgi(nbi).eq.25) then
c
c       Have an option that requires a second line of data to
c       complete the constraint.
c
        nfldtx = 3
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        ustr24 = ufield(1)
        uheadx = '->'
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)
        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
          write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
          qrderr = .true.
          go to 999
        endif
        ustr = ufield(2)
        ucospi(nbi)(1:48) = ustr
      endif
c
c     Read the next data line. Go back to process it.
c
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      go to 130
c
  180 continue
      nbti = nbi
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Ion exchanger creation. This section reads the directives for
c     creating ion exchange phases and species and their associated
c     intrinsic (including thermodynamic) properties. The data is
c     essentially that which might be read from a supporting data file.
c     Scenario-specific data (e.g., a specified amount of exchanger or
c     a specific composition for an exchanger phase) are not included
c     in this section.
c
      net = 0
      qnonep = .false.
c
c     Read a two-line header for the block.
c
      uheadx = 'Create Ion Exchangers'
      nfldtx = 2
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c     Loop on exchanger phases.
c
      ne = 0
      do nn = 1,netmax + 1
c
c       Read a line. If the block has not been completely read,
c       this contains the name of an exchanger phase, and a sub-block
c       for that phase follows. Otherwise, this line is the first line
c       of the next block.
c
        nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        ustr = ufield(1)
        uheadx = 'Exchanger phase'
        call locase(ustr)
        call locase(uheadx)
        j2 = ilnobl(ustr)
        j3 = ilnobl(uheadx)
c
        if (ustr(1:j2) .ne. uheadx(1:j3)) then
c
c         Back up.
c
          backspace ninpts
          go to 240
        endif
c
        ustr = ufield(2)
        ustrn = ustr
        call locase(ustrn)
        if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') ustr = 'None'
        qnonep = ustr(1:5).eq.'None '
c
        if (.not.qnonep) then
          net = net + 1
          ne = ne + 1
c
          if (ne .gt. netmax) then
            write (nttyo,1620) netmax,ne
 1620       format(/' * Error - (XCON3/rd3d8b) Have exceeded the',
     $      ' maximum number of ',i3,/7x,'generic ion exchange phases',
     $      ' while reading the data to create',/7x,'such phases.',
     $      ' Increase the dimensioning parameter netpar',/7x,'to at',
     $      ' least ',i3,'.')
            go to 990
          endif
c
          ugexp(ne) = ustr
          jgext(ne) = 0
        endif
c
c       Read the separator line following the data line containing
c       the name of an exchanger phase.
c
        nfldtx = 1
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        ustr24 = ufield(1)
        uheadx = '--------'
        if (ustr24(1:8) .ne. uheadx(1:8)) then
          j2 = ilnobl(uline1)
          j2 = min(j2,70)
          write (nttyo,1070) uline1(1:j2)
 1070     format(/' * Error - (XCON3/rd3d8b) Read the following line',
     $    ' where a separator line',/7x,'("|------- ... -------|")',
     $    ' was expected:',/7x,'"',a,'".')
          qrderr = .true.
          go to 999
        endif
c
c       Read the molecular weight of the bare exchange ligand (Z)
c       from a two-line header.
c
        uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $   nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
        if (qrderr) go to 999
        ustr24 = ufield(2)
        uheadx = 'Mol. Wt. (Z)'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)
        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
          write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
          qrderr = .true.
          go to 999
        endif
        ustr = ufield(3)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        if (.not.qnonep) mwtges(ne) = var
c
c       Read the exchange model name from a two-line header.
c
        uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $   nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
        if (qrderr) go to 999
        ustr24 = ufield(2)
        uheadx = 'Exchange model'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)
        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
          write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
          qrderr = .true.
          go to 999
        endif
        if (.not.qnonep) ugexmo(ne) = ufield(3)
c
c       Read the reference temperature (C) for the thermodynamic
c       data from a two-line header.
c
        uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $   nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
        if (qrderr) go to 999
        ustr = ufield(3)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        ustr24 = ufield(2)
        uheadx = 'Ref. Temp. (C)'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)
        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
          write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
          qrderr = .true.
          go to 999
        endif
        if (.not.qnonep) tgexp(ne) = var
c
c       Loop on exchange sites.
c
        je = 0
        do jj = 1,jetmax + 1
c
c         Read a line. If the sub-block (for the current exchanger
c         phase) has not been completely read, this contains the name
c         of an exchange site, and a sub-sub-block for that site
c         follows. Otherwise, this line is the first line of the next
c         sub-block (for the next exchanger phase).
c
          nfldtx = 0
          call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uline1,ulscr)
          if (qrderr) go to 999
          ustr = ufield(1)
          uheadx = '->'
          call locase(ustr)
          call locase(uheadx)
          j2 = ilnobl(ustr)
          j3 = ilnobl(uheadx)
          if (ustr(1:j2) .eq. uheadx(1:j3)) then
            ustr24 = ufield(2)
            uheadx = 'Exchange site'
            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)
            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
              write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
              qrderr = .true.
              go to 999
            endif
c
            ustr = ufield(3)
            ustrn = ustr
            call locase(ustrn)
            if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ')
     $      ustr = 'None'
            qnones = ustr(1:5).eq.'None '
c
            qokay = .not.qnones .and. .not.qnonep
            if (qokay) then
c
c             Have found another exchange site sub-sub-block.
c
              jgext(ne) = jgext(ne) + 1
              je = je + 1
c
              if (je .gt. jetmax) then
                j2 = ilnobl(ugexp(ne))
                write (nttyo,1710) jetmax,ugexp(ne)(1:j2),je
 1710           format(/' * Error - (XCON3/rd3d8b) Have exceeded the',
     $          ' maximum number of ',i3,/7x,'exchange sites on a',
     $          ' generic ion exchange phase while reading',/7x,'the',
     $          ' data to create ',a,'. Increase the',
     $          /7x,'dimensioning parameter jetpar to at least ',i3,'.')
                go to 990
              endif
c
              ugexj(je,ne) = ufield(3)
              ngexrt(je,ne) = 0
            endif
          else
c
c           Back up.
c
            backspace ninpts
c
            ustr = ufield(1)
            uheadx = 'Exchanger phase'
            call locase(ustr)
            call locase(uheadx)
            j2 = ilnobl(ustr)
            j3 = ilnobl(uheadx)
            if (ustr(1:j2) .eq. uheadx(1:j3)) then
c
c             Have found another exchanger phase.
c
              go to 230
            endif
c
c           Have found the end of the current block.
c
            go to 240
          endif
c
c         Read the separator line following the data line containing
c         the name of an exchange site.
c
          nfldtx = 1
          call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uline1,ulscr)
          if (qrderr) go to 999
          ustr24 = ufield(1)
          uheadx = '--------'
          if (ustr24(1:8) .ne. uheadx(1:8)) then
            j2 = ilnobl(uline1)
            j2 = min(j2,70)
            write (nttyo,1070) uline1(1:j2)
            qrderr = .true.
            go to 999
          endif
c
c         Read the site stoichiometric number (moles of site per mole
c         of bare exchange ligand) from a two-line header.
c
          uheadx = '--->'
          nfldtx = 4
          call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
          if (qrderr) go to 999
          ustr24 = ufield(2)
          uheadx = 'Stoich. number'
          call locase(ustr24)
          call locase(uheadx)
          j2 = ilnobl(ustr24)
          j3 = ilnobl(uheadx)
          if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
            qrderr = .true.
            go to 999
          endif
          ustr = ufield(3)
          call chreal(nttyo,qrderr,ustr,var)
          if (qrderr) go to 999
          if (qokay) cgexj(je,ne) = var
c
c         Read the intrinsic electrical charge of the site (charge
c         number for one mole of site) from a two-line header.
c
          uheadx = '--->'
          nfldtx = 4
          call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
          if (qrderr) go to 999
          ustr24 = ufield(2)
          uheadx = 'Electr. charge'
          call locase(ustr24)
          call locase(uheadx)
          j2 = ilnobl(ustr24)
          j3 = ilnobl(uheadx)
          if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
            qrderr = .true.
            go to 999
          endif
          ustr = ufield(3)
          call chreal(nttyo,qrderr,ustr,var)
          if (qrderr) go to 999
          if (qokay) zgexj(je,ne) = var
c
c         Loop on exchange reactions.
c
          n = 0
          do nnn = 1,ietmax + 1
c
c           Read a line. If the sub-sub-block (for the current
c           exchange site) has not been completely read, this contains
c           a string denoting an exchange reaction in condensed format,
c           and a sub-sub-sub-block for that reaction follows.
c           Otherwise, this line is the first line of the next
c           sub-sub-block (for the next site) or the first line of the
c           next sub-block for the next exchanger phase).
c
            nfldtx = 0
            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $      nttyo,qrderr,ufield,uline1,ulscr)
            if (qrderr) go to 999
            ustr = ufield(1)
            uheadx = '--->'
            call locase(ustr)
            call locase(uheadx)
            j2 = ilnobl(ustr)
            j3 = ilnobl(uheadx)
            if (ustr(1:j2) .eq. uheadx(1:j3)) then
              ustr24 = ufield(2)
              uheadx = 'Reaction'
              call locase(ustr24)
              call locase(uheadx)
              j2 = ilnobl(ustr24)
              j3 = ilnobl(uheadx)
              if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
                qrderr = .true.
                go to 999
              endif
c
              ustr = ufield(3)
              ustrn = ustr
              call locase(ustrn)
              if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ')
     $        ustr = 'None'
              qnoner = ustr(1:5).eq.'None '
c
              qokay = .not.qnoner .and. .not.qnones .and. .not.qnonep
              if (qokay) then
c
c               Have found another exchange reaction sub-sub-sub-block.
c
                ngexrt(je,ne) = ngexrt(je,ne) + 1
                n = n + 1
c
                if (n .gt. ietmax) then
                  j2 = ilnobl(ugexp(ne))
                  j3 = ilnobl(ugexj(ne,ne))
                  write (nttyo,1820) netmax,ugexj(je,ne)(1:j3),
     $            ugexp(ne)(1:j2),n
 1820             format(/' * Error - (XCON3/rd3d8b) Have exceeded the',
     $            ' maximum number of ',i3,/7x,'reactions for a site',
     $            ' belonging to a generic ion exchange',/7x,'phase',
     $            ' while reading the data for site ',a,' of exchange'
     $            ' phase',/7x,a,'. Increase the dimensioning',
     $            ' parameter',/7x,'ietpar to at least ',i3,'.')
                  go to 990
                endif
c
                ugexr(n,je,ne) = ufield(3)
              endif
            else
c
c             Back up.
c
              backspace ninpts
c
              ustr = ufield(1)
              uheadx = '->'
              call locase(ustr)
              call locase(uheadx)
              j2 = ilnobl(ustr)
              j3 = ilnobl(uheadx)
              if (ustr(1:j2) .eq. uheadx(1:j3)) then
                ustr = ufield(2)
                uheadx = 'Exchange site'
                call locase(ustr)
                call locase(uheadx)
                j2 = ilnobl(ustr)
                j3 = ilnobl(uheadx)
                if (ustr(1:j2) .eq. uheadx(1:j3)) then
c
c                 Have found another exchanger site.
c
                  go to 220
                endif
              endif
              ustr = ufield(1)
              uheadx = 'Exchanger phase'
              call locase(ustr)
              call locase(uheadx)
              j2 = ilnobl(ustr)
              j3 = ilnobl(uheadx)
              if (ustr(1:j2) .eq. uheadx(1:j3)) then
c
c               Have found another exchanger phase.
c
                go to 230
              endif
c
c             Have found the end of the current block.
c
              go to 240
            endif
c
c           Read the separator line following the data line containing
c           the string containing an exchange reaction in condensed
c           format.
c
            nfldtx = 1
            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $      nttyo,qrderr,ufield,uline1,ulscr)
            if (qrderr) go to 999
            ustr24 = ufield(1)
            uheadx = '--------'
            if (ustr24(1:8) .ne. uheadx(1:8)) then
              j2 = ilnobl(uline1)
              j2 = min(j2,70)
              write (nttyo,1070) uline1(1:j2)
              qrderr = .true.
              go to 999
            endif
c
c           Read the table header for the thermodynamic data for the
c           current exchange reaction from a two-line header.
c
            uheadx = '----->'
            nfldtx = 5
            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $      nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
            if (qrderr) go to 999
            ustr24 = ufield(2)
            uheadx = 'Parameter'
            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)
            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
              write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
              qrderr = .true.
              go to 999
            endif
c
c           Read the log K/Delta G data for the reaction from a
c           one-line header.
c
            uheadx = '----->'
            nfldtx = 5
            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $      nttyo,qrderr,ufield,uline1,ulscr)
            if (qrderr) go to 999
            ustr24 = ufield(2)
            uheadx = 'K func.'
            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)
            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
              write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
              qrderr = .true.
              go to 999
            endif
            ustr = ufield(3)
            call chreal(nttyo,qrderr,ustr,var)
            if (qrderr) go to 999
            if (qokay) then
              xlkgex(n,je,ne) = var
              uxkgex(n,je,ne) = ufield(4)
            endif
c
c           Read the Delta H data for the reaction from a one-line
c           header.
c
            uheadx = '----->'
            nfldtx = 5
            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $      nttyo,qrderr,ufield,uline1,ulscr)
            if (qrderr) go to 999
            ustr24 = ufield(2)
            uheadx = 'DelH0r'
            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)
            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
              write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
              qrderr = .true.
              go to 999
            endif
            ustr = ufield(3)
            call chreal(nttyo,qrderr,ustr,var)
            if (qrderr) go to 999
            if (qokay) then
              xhfgex(n,je,ne) = var
              uhfgex(n,je,ne) = ufield(4)
            endif
c
c           Read the Delta V data for the reaction from a two-line
c           header. The separator line completes the sub-sub-sub-block.
c
            uheadx = '----->'
            nfldtx = 5
            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $      nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
            if (qrderr) go to 999
            ustr24 = ufield(2)
            uheadx = 'DelV0r'
            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)
            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
              write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
              qrderr = .true.
              go to 999
            endif
            ustr = ufield(3)
            call chreal(nttyo,qrderr,ustr,var)
            if (qrderr) go to 999
            if (qokay) then
              xvfgex(n,je,ne) = var
              uvfgex(n,je,ne) = ufield(4)
            endif
c
          enddo
  220     continue
        enddo
  230   continue
      enddo
  240 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Specified generic ion exchanger compositions.
c
      neti = 0
      qnonep = .false.
c
c     Read a two-line header for the block.
c
      uheadx = 'Ion Exchanger Compositions'
      nfldtx = 2
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c     Loop on exchanger phases.
c
      nei = 0
      do nn = 1,netmax + 1
c
c       Read a line. If the block has not been completely read,
c       this contains the name of an exchanger phase, and a sub-block
c       for that phase follows. Otherwise, this line is the first line
c       of the next block.
c
        nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        ustr = ufield(1)
        uheadx = 'Exchanger phase'
        call locase(ustr)
        call locase(uheadx)
        j2 = ilnobl(ustr)
        j3 = ilnobl(uheadx)
c
        if (ustr(1:j2) .ne. uheadx(1:j3)) then
c
c         Back up.
c
          backspace ninpts
          go to 340
        endif
c
        ustr = ufield(2)
        ustrn = ustr
        call locase(ustrn)
        if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') ustr = 'None'
        qnonep = ustr(1:5).eq.'None '
c
        if (qnonep) then
          qgexef = .true.
        else
          neti = neti + 1
          nei = nei + 1
c
          if (nei .gt. netmax) then
            write (nttyo,1930) netmax,nei
 1930       format(/' * Error - (XCON3/rd3d8b) Have exceeded the',
     $      ' maximum number of ',i3,/7x,'generic ion exchange phases',
     $      ' while reading the data for concentrations',7x,'and',
     $      ' compositions of such phases. Increase the dimensioning',
     $      /7x,'parameter netpar to at least ',i3,'.')
            go to 990
          endif
c
          ugexpi(nei) = ustr
          j3 = ilnobl(ugexpi(nei))
          jgexti(nei) = 0
c
c         Find the corresponding ne index.
c
          do ne = 1,net
            j2 = ilnobl(ugexp(ne))
            if (ugexp(nei)(1:j3) .eq. ugexp(ne)(1:j2)) go to 250
          enddo
c
          write (nttyo,1932) ugexpi(nei)(1:j3)
 1932     format(/' * Error - (XCON3/rd3d8b) Data are present on the',
     $    ' input file for the',/7x,'generic ion exchange phase "',a,
     $    '", but this phase',/7x,"hasn't been previously defined.",
     $    " Can't determine the requisite model,",/7x,"therefore can't",
     $    ' finish reading the current data for this phase.')
          qrderr = .true.
          stop
c
  250     j2 = ilnobl(ugexmo(ne))
          if (ugexmo(ne)(1:j2) .eq. 'Gapon' .or.
     $      ugexmo(ne)(1:6) .eq. 'Gapon-' .or.
     $      ugexmo(ne)(1:j2) .eq. 'Vanselow' .or.
     $      ugexmo(ne)(1:9) .eq. 'Vanselow-') then
c
c           Input composition is described in terms of equivalent
c           fractions on the sites.
c
            qgexef = .true.
          else
c
c           Input composition is described in terms of mole fractions
c           on the sites.
c
            qgexef = .false.
          endif
        endif
c
c       Read the separator line following the data line containing
c       the name of an exchanger phase.
c
        nfldtx = 1
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        ustr24 = ufield(1)
        uheadx = '--------'
        if (ustr24(1:8) .ne. uheadx(1:8)) then
          j2 = ilnobl(uline1)
          j2 = min(j2,70)
          write (nttyo,1070) uline1(1:j2)
          qrderr = .true.
          go to 999
        endif
c
c       Read the concentration of the exchanger phase (moles/kg.H2O)
c       from a two-line header.
c
        uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $   nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
        if (qrderr) go to 999
        ustr24 = ufield(2)
        uheadx = 'Moles/kg.H2O'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)
        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
          write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
          qrderr = .true.
          go to 999
        endif
        ustr = ufield(3)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        if (.not.qnonep) cgexpi(nei) = var
c
c       Loop on exchange sites.
c
        jei = 0
        do jj = 1,jetmax + 1
c
c         Read a line. If the sub-block (for the current exchanger
c         phase) has not been completely read, this contains the name
c         of an exchange site, and a sub-sub-block for that site
c         follows. Otherwise, this line is the first line of the next
c         sub-block (for the next exchanger phase).
c
          nfldtx = 0
          call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uline1,ulscr)
          ustr = ufield(1)
          uheadx = '->'
          call locase(ustr)
          call locase(uheadx)
          j2 = ilnobl(ustr)
          j3 = ilnobl(uheadx)
          if (ustr(1:j2) .eq. uheadx(1:j3)) then
            ustr24 = ufield(2)
            uheadx = 'Exchange site'
            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)
            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
              write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
              qrderr = .true.
              go to 999
            endif
c
            ustr = ufield(3)
            ustrn = ustr
            call locase(ustrn)
            if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ')
     $      ustr = 'None'
            qnones = ustr(1:5).eq.'None '
c
            qokay = .not.qnones .and. .not.qnonep
            if (qokay) then
c
c             Have found another exchange site sub-sub-block.
c
              jgexti(ne) = jgexti(ne) + 1
              jei = jei + 1
c
              if (jei .gt. jetmax) then
                j2 = ilnobl(ugexpi(nei))
                write (nttyo,2000) jetmax,ugexpi(nei)(1:j2),jei
 2000           format(/' * Error - (XCON3/rd3d8b) Have exceeded the',
     $          ' maximum number of ',i3,/7x,'exchange sites on a',
     $          ' generic ion exchange phase while reading',/7x,'the',
     $          ' data for the concentration and composition of',
     $          /7x,a,'. Increase the dimensioning parameter jetpar',
     $          /7x,'to at least ',i3,'.')
                go to 990
              endif
c
              ugexji(jei,nei) = ufield(3)
              ngexti(jei,nei) = 0
            endif
          else
c
c           Back up.
c
            backspace ninpts
c
            ustr = ufield(1)
            uheadx = 'Exchanger phase'
            call locase(ustr)
            call locase(uheadx)
            j2 = ilnobl(ustr)
            j3 = ilnobl(uheadx)
            if (ustr(1:j2) .eq. uheadx(1:j3)) then
c
c             Have found another exchanger phase.
c
              go to 330
            endif
c
c           Have found the end of the current block.
c
            go to 340
          endif
c
c         Read the separator line following the data line containing
c         the name of an exchange site.
c
          nfldtx = 1
          call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uline1,ulscr)
          if (qrderr) go to 999
          ustr24 = ufield(1)
          uheadx = '--------'
          if (ustr24(1:8) .ne. uheadx(1:8)) then
            j2 = ilnobl(uline1)
            j2 = min(j2,70)
            write (nttyo,1070) uline1(1:j2)
            qrderr = .true.
            go to 999
          endif
c
c         Read a table header for the exchange site composition from
c         from a two-line header.
c
          uheadx = '--->'
          nfldtx = 4
          call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
          if (qrderr) go to 999
          ustr24 = ufield(2)
          uheadx = 'Exchange species'
          call locase(ustr24)
          call locase(uheadx)
          j2 = ilnobl(ustr24)
          j3 = ilnobl(uheadx)
          if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
            qrderr = .true.
            go to 999
          endif
c
c         Check the measure of concentration: equivalent fraction
c         versus mole fraction.
c
          ustr24 = ufield(3)
c
          if (qnonep) then
            if (ustr24(1:10) .eq. 'Mole frac.') then
              ustr24 = 'Eq. frac.'
            endif
          endif
c
          if (qgexef) then
            uheadx = 'Eq. frac.'
          else
            uheadx = 'Mole frac.'
          endif
          j2 = ilnobl(ustr24)
          j3 = ilnobl(uheadx)
          if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (nttyo,2010) uheadx(1:j3),ustr24(1:j2)
 2010       format(/' * Error - (XCON3/rd3d8b) Was expecting to find',
     $      ' a line containing',/7x,'"',a,'", instead found one',
     $      ' containing "',a,'".')
            if (.not.qnonep) then
              j4 = ilnobl(ugexmo(ne))
              j5 = ilnobl(ugexpi(nei))
              write (nttyo,2012) ugexmo(ne)(1:j4),ugexpi(nei)(1:j5)
 2012         format(7x,'The former is required for the ',a,
     $        ' exchange model,',/7x,'which was specified for',
     $        ' the generic ion exchange phase',/7x,a,'.')
            endif
            qrderr = .true.
            go to 999
          endif
c
c         Loop on exchange species.
c
          iei = 0
          do ii = 1,ietmax + 1
c
c           Read a line. If the sub-sub-block (for the current
c           exchange site) has not been completely read, this contains
c           the name of the exchange species and its mole fraction.
c           Otherwise, this line is the first line of the next
c           sub-sub-block (for the next site) or the first line of the
c           next sub-block for the next exchanger phase).
c
            nfldtx = 0
            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $      nttyo,qrderr,ufield,uline1,ulscr)
            if (qrderr) go to 999
            ustr = ufield(1)
            uheadx = '--->'
            call locase(ustr)
            call locase(uheadx)
            j2 = ilnobl(ustr)
            j3 = ilnobl(uheadx)
c
            if (ustr(1:j2) .eq. uheadx(1:j3)) then
              ustr = ufield(2)
              ustrn = ustr
              call locase(ustrn)
              if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ')
     $        ustr = 'None'
              qnonei = ustr(1:5).eq.'None '
c
              qokay = .not.qnonei .and. .not.qnones .and. .not.qnonep
              if (qokay) then
c
c               Have found another exchange species line.
c
                ngexti(jei,nei) = ngexti(jei,nei) + 1
                iei = iei + 1
c
                if (iei .gt. ietmax) then
                  write (nttyo,2040) ietmax,ugexji(jei,nei),
     $            ugexpi(nei),iei
 2040             format(/' * Error - (XCON3/rd3d8b) Have exceeded the',
     $            ' maximum number of ',i3,/7x,'species on an exchange',
     $            ' site of a generic ion exchange phase',/7x,'while',
     $            ' reading',/7x,'the data for the composition of',
     $            ' site ',a,' of',/7x,'exchange phase ',a,'. Increase',
     $            ' the dimensioning parameter ietpar to at least ',
     $            i3,'.')
                  go to 990
                endif
c
                ugexsi(iei,jei,nei) = ufield(2)
                ustr = ufield(3)
                call chreal(nttyo,qrderr,ustr,var)
                if (qrderr) go to 999
c
c               Check the array name noted: egexsi versus xgexsi.
c
                ustr24 = ' '
                k2 = index(uline1,'egexsi')
                k3 = index(uline1,'xgexsi')
                if (k2 .gt. 0) ustr24 = 'egexsi'
                if (k3 .gt. 0) ustr24 = 'xgexsi'
                if (qgexef) then
                  uheadx = 'egexsi'
                else
                  uheadx = 'xgexsi'
                endif
                j2 = ilnobl(ustr24)
                j3 = ilnobl(uheadx)
                j4 = ilnobl(ugexmo(ne))
                j5 = ilnobl(ugexpi(nei))
                if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                  write (nttyo,2014) uheadx(1:j3),ustr24(1:j2),
     $            ugexmo(ne)(1:j4),ugexpi(nei)(1:j5)
 2014             format(/' * Warning - (XCON3/rd3d8b) Was expecting',
     $            ' to find a line containing',/7x,'"',a,'", instead',
     $            ' found one containing "',a,'".',/7x,'The former',
     $            ' matches the ',a,' exchange model, which was',
     $            /7x,'specified for the generic ion exchange phase ',
     $            a,'.')
                endif
c
                if (qgexef) then
                  egexsi(iei,jei,nei) = var
                else
                  xgexsi(iei,jei,nei) = var
                endif
              endif
            else
              uheadx = '--------'
              if (ustr(1:8) .eq. uheadx(1:8)) then
c
c               Have found the end of the sub-sub-block for the
c               composition of the current exchange site
c
                go to 320
              else
c
c               Have unrecognized input.
c
                j3 = ilnobl(ugexji(jei,nei))
                j4 = ilnobl(ugexpi(nei))
                write (nttyo,1360) ustr(1:j2),ugexji(jei,nei)(1:j3),
     $          ugexpi(nei)(1:j4)
 1360           format(/' * Error - (XCON3/rd3d8b) Have unrecognized',
     $          ' input on the line',/7x,'beginning with "',a,'".',
     $          ' This line should contain the',/7x,'name of an',
     $          ' exchanger species and corresponding mole fraction',
     $          /7x,'for site ',a,' of exchanger phase ',a,', else it',
     $          /7x,'should be a separator line terminating a sequence',
     $          ' of such lines.')
c
              endif
            endif
c
          enddo
  320     continue
        enddo
  330   continue
      enddo
  340 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Specified solid solution compositions.
c
      nxti = 0
      nxic = 0
      qnonep = .false.
c
c     Read a two-line header for the block.
c
      uheadx = 'Solid Solution Compositions'
      nfldtx = 2
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c     Loop on solid solution phases.
c
      nxi = 0
      do nn = 1,nxtimx + 1
c
c       Read a line. If the block has not been completely read,
c       this contains the name of a solid solution, and a sub-block
c       for that phase follows. Otherwise, this line is the first line
c       of the next block.
c
        nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        ustr = ufield(1)
        uheadx = 'Solid Solution'
        call locase(ustr)
        call locase(uheadx)
        j2 = ilnobl(ustr)
        j3 = ilnobl(uheadx)
c
        if (ustr(1:j2) .ne. uheadx(1:j3)) then
c
c         Back up.
c
          backspace ninpts
          go to 440
        endif
c
        ustr = ufield(2)
        ustrn = ustr
        call locase(ustrn)
        if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') ustr = 'None'
        qnonep = ustr(1:5).eq.'None '
c
        if (.not.qnonep) then
          nxti = nxti + 1
          nxi = nxi + 1
c
          if (nxi .gt. nxtimx) then
            write (nttyo,1430) nxtimx,ustr(1:j2)
 1430       format(/' * Error - (XCON3/rd3d8b) Have exceeded the',
     $      ' maximum',/7x,'number of ',i5,' solid solutions for which',
     $      /7x,'compositions may be specified on the input file. This',
     $      /7x,'occurred while reading data for ',a,'.',
     $      /7x,'Increase the dimensioning parameter nxtipa.')
            go to 990
          endif
c
          usoli(nxi) = ustr
          ncmpri(1,nxi) = nxic + 1
        endif
c
c       Read the separator line following the data line containing
c       the name of a solid solution.
c
        nfldtx = 1
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        ustr24 = ufield(1)
        uheadx = '--------'
        if (ustr24(1:8) .ne. uheadx(1:8)) then
          j2 = ilnobl(uline1)
          j2 = min(j2,70)
          write (nttyo,1070) uline1(1:j2)
          qrderr = .true.
          go to 999
        endif
c
c       Read a table header for the solid solution composition from
c       from a two-line header.
c
        uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
        if (qrderr) go to 999
        ustr24 = ufield(2)
        uheadx = 'Component'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)
        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
          write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
          qrderr = .true.
          go to 999
        endif
c
c       Loop on component species.
c
        do nnn = 1,nxicmx + 1
c
c         Read a line. If the sub-block (for the current
c         solid solution) has not been completely read, this contains
c         the name of the component species and its mole fraction.
c         Otherwise, this line is the first line of the next
c         sub-block (for the next solid solution).
c
          nfldtx = 0
          call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uline1,ulscr)
          ustr = ufield(1)
          uheadx = '->'
          call locase(ustr)
          call locase(uheadx)
          j2 = ilnobl(ustr)
          j3 = ilnobl(uheadx)
          if (ustr(1:j2) .eq. uheadx(1:j3)) then
            ustr = ufield(2)
            ustrn = ustr
            call locase(ustrn)
            if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ')
     $      ustr = 'None'
            qnonei = ustr(1:5).eq.'None '
c
            qokay = .not.qnonep .and. .not.qnonei
            if (qokay) then
c
c             Have found another component species line.
c
              nxic = nxic + 1
c
              if (nxic .gt. nxicmx) then
                j2 = ilnobl(ustr)
                j3 = ilnobl(usoli(nxi))
                write (nttyo,1460) nxicmx,ustr(1:j2),usoli(nxi)(1:j3)
 1460           format(/' * Error - (XCON3/rd3d8b) Have exceeded the',
     $          ' maximum number',/7x,'of ',i5,' solid solution',
     $          ' species for which mole fractions',/7x,'are specified',
     $          ' on the input file. This occurred while reading',
     $          /7x,'data for the species ',a,'(',a,').',
     $          /7x,'Increase the dimensioning parameter nxipar.')
                go to 990
              endif
c
              umemi(nxic) = ufield(2)
              ustr = ufield(3)
              call chreal(nttyo,qrderr,ustr,var)
              if (qrderr) go to 999
              xbari(nxic) = var
            endif
          else
            uheadx = '--------'
            if (ustr(1:8) .eq. uheadx(1:8)) then
c
c             Have found the end of the sub-block for the
c             composition of the current solid solution.
c
              go to 420
            else
c
c             Have unrecognized input.
c
              j3 = ilnobl(usoli(nxi))
              write (nttyo,1470) ustr(1:j2),usoli(nxi)(1:j3)
 1470         format(/' * Error - (XCON3/rd3d8b) Have unrecognized',
     $        ' input on the line',/7x,'beginning with "',a,'".',
     $        ' This line should contain the',/7x,'name of a',
     $        ' component species and corresponding mole fraction',
     $        /7x,'for solid solution ',a,', else it should be a',
     $        ' separator line terminating a sequence of such lines.')
c
            endif
          endif
c
        enddo
  420   if (.not.qnonep) ncmpri(2,nxi) = nxic
      enddo
  440 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Nxmod options.
c
      uheadx = 'Alter/Suppress options'
      nfldtx = 2
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c     Read the first part of a table header from a one-line header.
c
      uheadx = 'Species'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
c
c     Read the second part of the table header from a two-line header.
c
      uheadx = '(uxmod(n))'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
      n = 0
c
c     Loop on the number of alter/suppress options.
c
      do nn = 1,nxmdmx + 1
c
c       Read a line. This contains an alter/suppress option, else it is
c       a separator line marking the end of the current block.
c
        nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
c
        ustr = ufield(1)
        if (ustr(1:8) .eq. '--------') go to 570
c
        call locase(ustr)
        ustrn = ustr
        call locase(ustrn)
        if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') ustr = 'None'
        if (ustr(1:5) .eq. 'None ') go to 560
c
        n = n + 1
c
        if (n .gt. nxmdmx) then
          write (nttyo,1500) nxmdmx
 1500     format(/' * Error - (XCON3/rd3d8b) Have too many nxmod',
     $    /7x,'alter/suppress options. The code is only dimensioned',
     $    /7x,'for ',i3,' such options. Reduce the number of such',
     $    /7x,'options or increase the dimensioning parameter nxmdpa.')
          go to 990
        endif
c
        uxmod(n) = ufield(1)
c
        ustr = ufield(2)
        call locase(ustr)
        do kxmd = -1,2
          uheadx = ukxm(kxmd)
          call locase(uheadx)
          if (ustr(1:16) .eq. uheadx(1:16)) then
            kxmod(n) = kxmd
            go to 152
          endif
        enddo
c
        j2 = ilnobl(ustr)
        write (nttyo,1510) ustr(1:j2)
 1510   format(/" * Error - (XCON3/rd3d8b) Don't recognize the",
     $  ' alter/suppress option',/7x,'string "',a,'". This should',
     $  ' be one of the strings',/7x,'defined in the ukxm array.',
     $  ' The valid strings are:',/)
        do kxmd = -1,2
          j3 = ilnobl(ukxm(kxmd))
          write (nttyo,1514) ukxm(kxmd)(1:j3)
 1514     format(9x,a)
        enddo
        go to 990
c
  152   ustr = ufield(3)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        xlkmod(n) = var
c
  560   continue
      enddo
c
  570 nxmod = n
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iopt Model Option Switches.
c     Note: iopt(1) = iopt1, etc.
c
      uheadx = 'Iopt Model Option Switches ("( 0)" marks default choices
     $)'
      nfldtx = 1
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c     Loop on option switches.
c
      do nn = 1,noptmx
c
c       Read the option title string from a one-line header.
c
        nfldtx = 1
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        uheadx = ufield(1)
        j2 = ilnobl(uheadx)
c
c       Check for the end of the block.
c
        if (uheadx(1:4) .ne. 'iopt') then
c
c         Back up.
c
          backspace ninpts
          go to 650
        endif
c
        k1 = index(uheadx,'(')
        k2 = index(uheadx,')')
        k3 = index(uheadx,'- ')
c
        if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or. k3.le.0 .or.
     $  k3.le.k2) then
          write (nttyo,1600) uheadx(1:j2)
 1600     format(/' * Error - (XCON3/rd3d8b) The iopt option switch',
     $    ' title line',/7x,'"',a,'"',/7x,'read from the',
     $    " input file isn't in the required format.")
          go to 990
        endif
c
c       Get the index of the option.
c
        ustr = uheadx(k1 + 1:k2 - 1)
        call chrint(ivar,nttyo,qrderr,ustr)
        if (qrderr) go to 999
        n = ivar
c
        if (n.lt.1 .or. n.gt.nptxpa) then
          write (nttyo,1610) uheadx(1:j2),n
 1610     format(/' * Error - (XCON3/rd3d8b) The iopt option switch',
     $    ' title line',/7x,'"',a,'"',/7x,'read from the',
     $    ' input file references an option switch index',
     $    /7x,'of ',i3,', which is out of range.')
          go to 990
        endif
c
c       Check the full title string.
c
        ustr = uopttx(n)
        j3 = ilnobl(ustr)
        if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
          write (nttyo,1612) uheadx(1:j2),ustr(1:j3)
 1612     format(/' * Error - (XCON3/rd3d8b) The iopt option switch',
     $    ' title string',/7x,'"',a,'"',
     $    /7x,"read from the input file doesn't contain the",
     $    ' matching defined string',/7x,'"',a,'".')
          go to 990
        endif
c
        iopt(n) = 0
        nmark = 0
c
c       Loop on option choices.
c
        do jj = 1,jptxpa + 1
c
c         Read a line. This contains an option choice line, else it
c         is a separator line marking the end of the option choice
c         lines for the current option.
c
          nfldtx = 1
          call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uline1,ulscr)
          if (qrderr) go to 999
c
          uheadx = ufield(1)
          j2 = ilnobl(uheadx)
          if (uheadx(1:8) .eq. '--------') go to 620
c
          if (jj .gt. jptxpa) then
            j3 = ilnobl(uopttx(n))
            write (nttyo,1630) uopttx(n)(1:j3),jptxpa
 1630       format(/' * Error - (XCON3/rd3d8b) Have too many option',
     $      /7x,'choice lines for the iopt option switch whose title',
     $      ' string is',/7x,'"',a,'".',/7x,' The code is only',
     $      ' dimensioned for ',i3,' such lines. Reduce the',
     $      /7x,'number of such lines or increase the dimensioning',
     $      ' parameter jptxpa.')
            go to 990
          endif
c
          k1 = index(uheadx,'[')
          k2 = index(uheadx,']')
          k3 = index(uheadx,'(')
          k4 = index(uheadx,')')
c
          if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or.
     $    k3.le.0 .or. k4.le.0 .or. k4.lt.(k3 + 2)) then
            write (nttyo,1640) uheadx(1:j2)
 1640       format(/' * Error - (XCON3/rd3d8b) The following iopt',
     $      /7x,'option switch choice line read from the input file',
     $      /7x,"isn't in the required format:",/7x,'"',a,'"')
            go to 990
          endif
c
c         Get the index of the option choice.
c
          ustr = uheadx(k3 + 1:k4 - 1)
          call chrint(ivar,nttyo,qrderr,ustr)
          if (qrderr) go to 999
          ival = ivar
c
          do j = 1,jptxpa
            ii = ioptox(j,n)
            if (ival .eq. ii) go to 600
          enddo
c
          write (nttyo,1650) uheadx(1:j2)
 1650     format(/' * Error - (XCON3/rd3d8b) The iopt option switch',
     $    /7x,'choice line',/7x,'"',a,'"',
     $    /7x,'read from the input file references an out-of-range',
     $    ' option',/7x,'choice index.')
          go to 990
c
  600     continue
c
c         Check the full option choice string.
c
          ustr = uoptox(j,n)
          j3 = ilnobl(ustr)
          j4 = ilnobl(uopttx(n))
          if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
            write (nttyo,1660) uheadx(1:j2),ustr(1:j3),uopttx(n)(1:j4)
 1660       format(/' * Error - (XCON3/rd3d8b) The iopt option switch',
     $      ' choice line',/7x,'"',a,'"',
     $      /7x," read from the input file doesn't contain the",
     $      ' matching defined string',/7x,'"',a,'".',
     $      /7x,'This line belongs to the option switch whose title',
     $      ' string is',/7x,'"',a,'".')
            go to 990
          endif
c
          k1 = index(uheadx,'[')
          k2 = index(uheadx,']')
          ustr24 = uheadx(k1 + 1:k2 -1)
          call lejust(ustr24)
          j3 = ilnobl(ustr24)
          qmark = .false.
          if (j3 .gt. 0) then
            if (index(ustr24,'*') .ge. 1) then
              qmark = .true.
            elseif (index(ustr24,'x') .ge. 1) then
              qmark = .true.
            elseif (index(ustr24,'X') .ge. 1) then
              qmark = .true.
            else
              j4 = ilnobl(uopttx(n))
              write (nttyo,1670) ustr24(1:j3),uheadx(1:j2),
     $        uopttx(n)(1:j4)
 1670         format(/" * Error - (XCON3/rd3d8b) Don't recognize the",
     $        ' string "',a,'"',/7x,'that appears on the iopt',
     $        ' option switch choice line',/7x,'"',a,'"',
     $        /7x,'read from the input file. An option choice should',
     $        ' be chosen by',/7x,'placing a "*", "x", or "X" in the',
     $        ' checkbox ("[ ]"). This',/7x,'choice line belongs to',
     $        ' the option switch whose title string is',/7x,'"',a,'".')
              go to 990
            endif
          endif
c
          if (qmark) then
            nmark = nmark + 1
            iopt(n) = ival
            jlast = j
          endif
c
        enddo
  620   continue
c
        if (nmark .le. 0) then
          j2 = ilnobl(uopttx(n))
          write (nttyo,1680) uopttx(n)(1:j2)
 1680     format(/' * Warning - (XCON3/rd3d8b) No option choice was',
     $    ' checked on the input file',/7x,'for the iopt option',
     $    ' switch whose title string is',/7x,'"',a,'".')
c
          do j = 1,jptxpa
            ival = ioptox(j,n)
            if (ival .eq. 0) go to 630
          enddo
c
          write (nttyo,1690)
 1690     format(/7x,'A default value of 0 has been applied, but no',
     $    /7x,'matching option choice string is defined.')
          go to 640
c
  630     j3 = ilnobl(uoptox(j,n))
          write (nttyo,1692) uoptox(j,n)(1:j3)
 1692     format(/7x,'A default value of 0 has been applied. The',
     $    ' matching string is',/7x,'"',a,'".')
  640     continue
        endif
c
        if (nmark .gt. 1) then
          j2 = ilnobl(uopttx(n))
          j = jlast
          j3 = ilnobl(uoptox(j,n))
           write (nttyo,1694) uopttx(n)(1:j2),uoptox(j,n)(1:j3)
 1694     format(/' * Warning - (XCON3/rd3d8b) More than one option',
     $    ' choice was checked',/7x,'on the input file for the iopt',
     $    ' option switch whose title string is',/7x,'"',a,'".',
     $    /7x,'The last choice checked will be used. The',
     $    ' matching string is',/7x,'"',a,'".')
        endif
c
      enddo
  650 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iopg Print Option Switches.
c     Note: iopg(1) = iopg1, etc.
c
      uheadx = 'Iopg Activity Coefficient Option Switches ("( 0)" marks
     $default choices)'
      nfldtx = 1
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c     Loop on option switches.
c
      do nn = 1,nopgmx
c
c       Read the option title string from a one-line header.
c
        nfldtx = 1
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        uheadx = ufield(1)
        j2 = ilnobl(uheadx)
c
c       Check for the end of the block.
c
        if (uheadx(1:4) .ne. 'iopg') then
c
c         Back up.
c
          backspace ninpts
          go to 750
        endif
c
        k1 = index(uheadx,'(')
        k2 = index(uheadx,')')
        k3 = index(uheadx,'- ')
c
        if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or. k3.le.0 .or.
     $  k3.le.k2) then
          write (nttyo,1700) uheadx(1:j2)
 1700     format(/' * Error - (XCON3/rd3d8b) The iopg option switch',
     $    ' title line',/7x,'"',a,'"',/7x,'read from the',
     $    " input file isn't in the required format.")
          go to 990
        endif
c
c       Get the index of the option.
c
        ustr = uheadx(k1 + 1:k2 - 1)
        call chrint(ivar,nttyo,qrderr,ustr)
        if (qrderr) go to 999
        n = ivar
c
        if (n.lt.1 .or. n.gt.npgxpa) then
          write (nttyo,1712) uheadx(1:j2),n
 1712     format(/' * Error - (XCON3/rd3d8b) The iopg option switch',
     $    ' title line',/7x,'"',a,'"',/7x,'read from the',
     $    ' input file references an option switch index',
     $    /7x,'of ',i3,', which is out of range.')
          go to 990
        endif
c
c       Check the full title string.
c
        ustr = uopgtx(n)
        j3 = ilnobl(ustr)
        if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
          write (nttyo,1720) uheadx(1:j2),ustr(1:j3)
 1720     format(/' * Error - (XCON3/rd3d8b) The iopg option switch',
     $    ' title string',/7x,'"',a,'"',
     $    /7x,"read from the input file doesn't contain the",
     $    ' matching defined string',/7x,'"',a,'".')
          go to 990
        endif
c
        iopg(n) = 0
        nmark = 0
c
c       Loop on option choices.
c
        do jj = 1,jpgxpa + 1
c
c         Read a line. This contains an option choice line, else it
c         is a separator line marking the end of the option choice
c         lines for the current option.
c
          nfldtx = 1
          call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uline1,ulscr)
          if (qrderr) go to 999
c
          uheadx = ufield(1)
          j2 = ilnobl(uheadx)
          if (uheadx(1:8) .eq. '--------') go to 720
c
          if (jj .gt. jpgxpa) then
            j3 = ilnobl(uopgtx(n))
            write (nttyo,1730) uopgtx(n)(1:j3),jpgxpa
 1730       format(/' * Error - (XCON3/rd3d8b) Have too many option',
     $      /7x,'choice lines for the iopg option switch whose title',
     $      ' string is',/7x,'"',a,'".',/7x,' The code is only',
     $      ' dimensioned for ',i3,' such lines. Reduce the',
     $      /7x,'number of such lines or increase the dimensioning',
     $      ' parameter jpgxpa.')
            go to 990
          endif
c
          k1 = index(uheadx,'[')
          k2 = index(uheadx,']')
          k3 = index(uheadx,'(')
          k4 = index(uheadx,')')
c
          if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or.
     $    k3.le.0 .or. k4.le.0 .or. k4.lt.(k3 + 2)) then
            write (nttyo,1740) uheadx(1:j2)
 1740       format(/' * Error - (XCON3/rd3d8b) The following iopg',
     $      /7x,'option switch choice line read from the input file',
     $      /7x,"isn't in the required format:",/7x,'"',a,'"')
            go to 990
          endif
c
c         Get the index of the option choice.
c
          ustr = uheadx(k3 + 1:k4 - 1)
          call chrint(ivar,nttyo,qrderr,ustr)
          if (qrderr) go to 999
          ival = ivar
c
          do j = 1,jpgxpa
            ii = iopgox(j,n)
            if (ival .eq. ii) go to 700
          enddo
c
          write (nttyo,1750) uheadx(1:j2)
 1750     format(/' * Error - (XCON3/rd3d8b) The iopg option switch',
     $    /7x,'choice line',/7x,'"',a,'"',
     $    /7x,'read from the input file references an out-of-range',
     $    ' option',/7x,'choice index.')
          go to 990
c
  700     continue
c
c         Check the full option choice string.
c
          ustr = uopgox(j,n)
          j3 = ilnobl(ustr)
          j4 = ilnobl(uopgtx(n))
          if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
            write (nttyo,1760) uheadx(1:j2),ustr(1:j3),uopgtx(n)(1:j4)
 1760       format(/' * Error - (XCON3/rd3d8b) The iopg option switch',
     $      ' choice line',/7x,'"',a,'"',
     $      /7x," read from the input file doesn't contain the",
     $      ' matching defined string',/7x,'"',a,'".',
     $      /7x,'This line belongs to the option switch whose title',
     $      ' string is',/7x,'"',a,'".')
            go to 990
          endif
c
          k1 = index(uheadx,'[')
          k2 = index(uheadx,']')
          ustr24 = uheadx(k1 + 1:k2 -1)
          call lejust(ustr24)
          j3 = ilnobl(ustr24)
          qmark = .false.
          if (j3 .gt. 0) then
            if (index(ustr24,'*') .ge. 1) then
              qmark = .true.
            elseif (index(ustr24,'x') .ge. 1) then
              qmark = .true.
            elseif (index(ustr24,'X') .ge. 1) then
              qmark = .true.
            else
              j4 = ilnobl(uopgtx(n))
              write (nttyo,1770) ustr24(1:j3),uheadx(1:j2),
     $        uopgtx(n)(1:j4)
 1770         format(/" * Error - (XCON3/rd3d8b) Don't recognize the",
     $        ' string "',a,'"',/7x,'that appears on the iopg',
     $        ' option switch choice line',/7x,'"',a,'"',
     $        /7x,'read from the input file. An option choice should',
     $        ' be chosen by',/7x,'placing a "*", "x", or "X" in the',
     $        ' checkbox ("[ ]"). This',/7x,'choice line belongs to',
     $        ' the option switch whose title string is',/7x,'"',a,'".')
              go to 990
            endif
          endif
c
          if (qmark) then
            nmark = nmark + 1
            iopg(n) = ival
            jlast = j
          endif
c
        enddo
  720   continue
c
        if (nmark .le. 0) then
          j2 = ilnobl(uopgtx(n))
          write (nttyo,1780) uopgtx(n)(1:j2)
 1780     format(/' * Warning - (XCON3/rd3d8b) No option choice was',
     $    ' checked on the input file',/7x,'for the iopg option',
     $    ' switch whose title string is',/7x,'"',a,'".')
c
          do j = 1,jpgxpa
            ival = iopgox(j,n)
            if (ival .eq. 0) go to 730
          enddo
c
          write (nttyo,1790)
 1790     format(/7x,'A default value of 0 has been applied, but no',
     $    /7x,'matching option choice string is defined.')
          go to 740
c
  730     j3 = ilnobl(uopgox(j,n))
          write (nttyo,1792) uopgox(j,n)(1:j3)
 1792     format(/7x,'A default value of 0 has been applied. The',
     $    ' matching string is',/7x,'"',a,'".')
  740     continue
        endif
c
        if (nmark .gt. 1) then
          j2 = ilnobl(uopgtx(n))
          j = jlast
          j3 = ilnobl(uopgox(j,n))
           write (nttyo,1794) uopgtx(n)(1:j2),uopgox(j,n)(1:j3)
 1794     format(/' * Warning - (XCON3/rd3d8b) More than one option',
     $    ' choice was checked',/7x,'on the input file for the iopg',
     $    ' option switch whose title string is',/7x,'"',a,'".',
     $    /7x,'The last choice checked will be used. The',
     $    ' matching string is',/7x,'"',a,'".')
        endif
c
      enddo
  750 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iopr Print Option Switches.
c     Note: iopr(1) = iopt1, etc.
c
      uheadx = 'Iopr Print Option Switches ("( 0)" marks default choices
     $)'
      nfldtx = 1
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c     Loop on option switches.
c
      do nn = 1,noprmx
c
c       Read the option title string from a one-line header.
c
        nfldtx = 1
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        uheadx = ufield(1)
        j2 = ilnobl(uheadx)
c
c       Check for the end of the block.
c
        if (uheadx(1:4) .ne. 'iopr') then
c
c         Back up.
c
          backspace ninpts
          go to 850
        endif
c
        k1 = index(uheadx,'(')
        k2 = index(uheadx,')')
        k3 = index(uheadx,'- ')
c
        if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or. k3.le.0 .or.
     $  k3.le.k2) then
          write (nttyo,1800) uheadx(1:j2)
 1800     format(/' * Error - (XCON3/rd3d8b) The iopr option switch',
     $    ' title line',/7x,'"',a,'"',/7x,'read from the',
     $    " input file isn't in the required format.")
          go to 990
        endif
c
c       Get the index of the option.
c
        ustr = uheadx(k1 + 1:k2 - 1)
        call chrint(ivar,nttyo,qrderr,ustr)
        if (qrderr) go to 999
        n = ivar
c
        if (n.lt.1 .or. n.gt.nprxpa) then
          write (nttyo,1810) uheadx(1:j2),n
 1810     format(/' * Error - (XCON3/rd3d8b) The iopr option switch',
     $    ' title line',/7x,'"',a,'"',/7x,'read from the',
     $    ' input file references an option switch index',
     $    /7x,'of ',i3,', which is out of range.')
          go to 990
        endif
c
c       Check the full title string.
c
        ustr = uoprtx(n)
        j3 = ilnobl(ustr)
        if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
          write (nttyo,1812) uheadx(1:j2),ustr(1:j3)
 1812     format(/' * Error - (XCON3/rd3d8b) The iopr option switch',
     $    ' title string',/7x,'"',a,'"',
     $    /7x,"read from the input file doesn't contain the",
     $    ' matching defined string',/7x,'"',a,'".')
          go to 990
        endif
c
        iopr(n) = 0
        nmark = 0
c
c       Loop on option choices.
c
        do jj = 1,jprxpa + 1
c
c         Read a line. This contains an option choice line, else it
c         is a separator line marking the end of the option choice
c         lines for the current option.
c
          nfldtx = 1
          call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uline1,ulscr)
          if (qrderr) go to 999
c
          uheadx = ufield(1)
          j2 = ilnobl(uheadx)
          if (uheadx(1:8) .eq. '--------') go to 820
c
          if (jj .gt. jprxpa) then
            j3 = ilnobl(uoprtx(n))
            write (nttyo,1830) uoprtx(n)(1:j3),jprxpa
 1830       format(/' * Error - (XCON3/rd3d8b) Have too many option',
     $      /7x,'choice lines for the iopr option switch whose title',
     $      ' string is',/7x,'"',a,'".',/7x,' The code is only',
     $      ' dimensioned for ',i3,' such lines. Reduce the',
     $      /7x,'number of such lines or increase the dimensioning',
     $      ' parameter jprxpa.')
            go to 990
          endif
c
          k1 = index(uheadx,'[')
          k2 = index(uheadx,']')
          k3 = index(uheadx,'(')
          k4 = index(uheadx,')')
c
          if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or.
     $    k3.le.0 .or. k4.le.0 .or. k4.lt.(k3 + 2)) then
            write (nttyo,1840) uheadx(1:j2)
 1840       format(/' * Error - (XCON3/rd3d8b) The following iopr',
     $      /7x,'option switch choice line read from the input file',
     $      /7x,"isn't in the required format:",/7x,'"',a,'"')
            go to 990
          endif
c
c         Get the index of the option choice.
c
          ustr = uheadx(k3 + 1:k4 - 1)
          call chrint(ivar,nttyo,qrderr,ustr)
          if (qrderr) go to 999
          ival = ivar
c
          do j = 1,jprxpa
            ii = ioprox(j,n)
            if (ival .eq. ii) go to 800
          enddo
c
          write (nttyo,1850) uheadx(1:j2)
 1850     format(/' * Error - (XCON3/rd3d8b) The iopr option switch',
     $    /7x,'choice line',/7x,'"',a,'"',
     $    /7x,'read from the input file references an out-of-range',
     $    ' option',/7x,'choice index.')
          go to 990
c
  800     continue
c
c         Check the full option choice string.
c
          ustr = uoprox(j,n)
          j3 = ilnobl(ustr)
          j4 = ilnobl(uoprtx(n))
          if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
            write (nttyo,1860) uheadx(1:j2),ustr(1:j3),uoprtx(n)(1:j4)
 1860       format(/' * Error - (XCON3/rd3d8b) The iopr option switch',
     $      ' choice line',/7x,'"',a,'"',
     $      /7x," read from the input file doesn't contain the",
     $      ' matching defined string',/7x,'"',a,'".',
     $      /7x,'This line belongs to the option switch whose title',
     $      ' string is',/7x,'"',a,'".')
            go to 990
          endif
c
          k1 = index(uheadx,'[')
          k2 = index(uheadx,']')
          ustr24 = uheadx(k1 + 1:k2 -1)
          call lejust(ustr24)
          j3 = ilnobl(ustr24)
          qmark = .false.
          if (j3 .gt. 0) then
            if (index(ustr24,'*') .ge. 1) then
              qmark = .true.
            elseif (index(ustr24,'x') .ge. 1) then
              qmark = .true.
            elseif (index(ustr24,'X') .ge. 1) then
              qmark = .true.
            else
              j4 = ilnobl(uoprtx(n))
              write (nttyo,1870) ustr24(1:j3),uheadx(1:j2),
     $        uoprtx(n)(1:j4)
 1870         format(/" * Error - (XCON3/rd3d8b) Don't recognize the",
     $        ' string "',a,'"',/7x,'that appears on the iopr',
     $        ' option switch choice line',/7x,'"',a,'"',
     $        /7x,'read from the input file. An option choice should',
     $        ' be chosen by',/7x,'placing a "*", "x", or "X" in the',
     $        ' checkbox ("[ ]"). This',/7x,'choice line belongs to',
     $        ' the option switch whose title string is',/7x,'"',a,'".')
              go to 990
            endif
          endif
c
          if (qmark) then
            nmark = nmark + 1
            iopr(n) = ival
            jlast = j
          endif
c
        enddo
  820   continue
c
        if (nmark .le. 0) then
          j2 = ilnobl(uoprtx(n))
          write (nttyo,1880) uoprtx(n)(1:j2)
 1880     format(/' * Warning - (XCON3/rd3d8b) No option choice was',
     $    ' checked on the input file',/7x,'for the iopr option',
     $    ' switch whose title string is',/7x,'"',a,'".')
c
          do j = 1,jprxpa
            ival = ioprox(j,n)
            if (ival .eq. 0) go to 830
          enddo
c
          write (nttyo,1890)
 1890     format(/7x,'A default value of 0 has been applied, but no',
     $    /7x,'matching option choice string is defined.')
          go to 840
c
  830     j3 = ilnobl(uoprox(j,n))
          write (nttyo,1892) uoprox(j,n)(1:j3)
 1892     format(/7x,'A default value of 0 has been applied. The',
     $    ' matching string is',/7x,'"',a,'".')
  840     continue
        endif
c
        if (nmark .gt. 1) then
          j2 = ilnobl(uoprtx(n))
          j = jlast
          j3 = ilnobl(uoprox(j,n))
           write (nttyo,1894) uoprtx(n)(1:j2),uoprox(j,n)(1:j3)
 1894     format(/' * Warning - (XCON3/rd3d8b) More than one option',
     $    ' choice was checked',/7x,'on the input file for the iopr',
     $    ' option switch whose title string is',/7x,'"',a,'".',
     $    /7x,'The last choice checked will be used. The',
     $    ' matching string is',/7x,'"',a,'".')
        endif
c
      enddo
  850 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iodb debugging print option switches.
c     Note: iodb(1) = iodb1, etc.
c
      uheadx = 'Iodb Debugging Print Option Switches ("( 0)" marks defau
     $lt choices)'
      nfldtx = 1
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c     Loop on option switches.
c
      do nn = 1,nodbmx
c
c       Read the option title string from a one-line header.
c
        nfldtx = 1
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        uheadx = ufield(1)
        j2 = ilnobl(uheadx)
c
c       Check for the end of the block.
c
        if (uheadx(1:4) .ne. 'iodb') then
c
c         Back up.
c
          backspace ninpts
          go to 950
        endif
c
        k1 = index(uheadx,'(')
        k2 = index(uheadx,')')
        k3 = index(uheadx,'- ')
c
        if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or. k3.le.0 .or.
     $  k3.le.k2) then
          write (nttyo,1900) uheadx(1:j2)
 1900     format(/' * Error - (XCON3/rd3d8b) The iodb option switch',
     $    ' title line',/7x,'"',a,'"',/7x,'read from the',
     $    " input file isn't in the required format.")
          go to 990
        endif
c
c       Get the index of the option.
c
        ustr = uheadx(k1 + 1:k2 - 1)
        call chrint(ivar,nttyo,qrderr,ustr)
        if (qrderr) go to 999
        n = ivar
c
        if (n.lt.1 .or. n.gt.ndbxpa) then
          write (nttyo,1910) uheadx(1:j2),n
 1910     format(/' * Error - (XCON3/rd3d8b) The iodb option switch',
     $    ' title line',/7x,'"',a,'"',/7x,'read from the',
     $    ' input file references an option switch index',
     $    /7x,'of ',i3,', which is out of range.')
          go to 990
        endif
c
c       Check the full title string.
c
        ustr = uodbtx(n)
        j3 = ilnobl(ustr)
        if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
          write (nttyo,1920) uheadx(1:j2),ustr(1:j3)
 1920     format(/' * Error - (XCON3/rd3d8b) The iodb option switch',
     $    ' title string',/7x,'"',a,'"',
     $    /7x,"read from the input file doesn't contain the",
     $    ' matching defined string',/7x,'"',a,'".')
          go to 990
        endif
c
        iodb(n) = 0
        nmark = 0
c
c       Loop on option choices.
c
        do jj = 1,jdbxpa + 1
c
c         Read a line. This contains an option choice line, else it
c         is a separator line marking the end of the option choice
c         lines for the current option.
c
          nfldtx = 1
          call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uline1,ulscr)
          if (qrderr) go to 999
c
          uheadx = ufield(1)
          j2 = ilnobl(uheadx)
          if (uheadx(1:8) .eq. '--------') go to 920
c
          if (jj .gt. jdbxpa) then
            j3 = ilnobl(uodbtx(n))
            write (nttyo,1922) uodbtx(n)(1:j3),jdbxpa
 1922       format(/' * Error - (XCON3/rd3d8b) Have too many option',
     $      /7x,'choice lines for the iodb option switch whose title',
     $      ' string is',/7x,'"',a,'".',/7x,' The code is only',
     $      ' dimensioned for ',i3,' such lines. Reduce the',
     $      /7x,'number of such lines or increase the dimensioning',
     $      ' parameter jdbxpa.')
            go to 990
          endif
c
          k1 = index(uheadx,'[')
          k2 = index(uheadx,']')
          k3 = index(uheadx,'(')
          k4 = index(uheadx,')')
c
          if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2) .or.
     $    k3.le.0 .or. k4.le.0 .or. k4.lt.(k3 + 2)) then
            write (nttyo,1940) uheadx(1:j2)
 1940       format(/' * Error - (XCON3/rd3d8b) The following iodb',
     $      /7x,'option switch choice line read from the input file',
     $      /7x,"isn't in the required format:",/7x,'"',a,'"')
            go to 990
          endif
c
c         Get the index of the option choice.
c
          ustr = uheadx(k3 + 1:k4 - 1)
          call chrint(ivar,nttyo,qrderr,ustr)
          if (qrderr) go to 999
          ival = ivar
c
          do j = 1,jdbxpa
            ii = iodbox(j,n)
            if (ival .eq. ii) go to 900
          enddo
c
          write (nttyo,1950) uheadx(1:j2)
 1950     format(/' * Error - (XCON3/rd3d8b) The iodb option switch',
     $    /7x,'choice line',/7x,'"',a,'"',
     $    /7x,'read from the input file references an out-of-range',
     $    ' option',/7x,'choice index.')
          go to 990
c
  900     continue
c
c         Check the full option choice string.
c
          ustr = uodbox(j,n)
          j3 = ilnobl(ustr)
          j4 = ilnobl(uodbtx(n))
          if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
            write (nttyo,1960) uheadx(1:j2),ustr(1:j3),uodbtx(n)(1:j4)
 1960       format(/' * Error - (XCON3/rd3d8b) The iodb option switch',
     $      ' choice line',/7x,'"',a,'"',
     $      /7x," read from the input file doesn't contain the",
     $      ' matching defined string',/7x,'"',a,'".',
     $      /7x,'This line belongs to the option switch whose title',
     $      ' string is',/7x,'"',a,'".')
            go to 990
          endif
c
          k1 = index(uheadx,'[')
          k2 = index(uheadx,']')
          ustr24 = uheadx(k1 + 1:k2 -1)
          call lejust(ustr24)
          j3 = ilnobl(ustr24)
          qmark = .false.
          if (j3 .gt. 0) then
            if (index(ustr24,'*') .ge. 1) then
              qmark = .true.
            elseif (index(ustr24,'x') .ge. 1) then
              qmark = .true.
            elseif (index(ustr24,'X') .ge. 1) then
              qmark = .true.
            else
              j4 = ilnobl(uodbtx(n))
              write (nttyo,1970) ustr24(1:j3),uheadx(1:j2),
     $        uodbtx(n)(1:j4)
 1970         format(/" * Error - (XCON3/rd3d8b) Don't recognize the",
     $        ' string "',a,'"',/7x,'that appears on the iodb',
     $        ' option switch choice line',/7x,'"',a,'"',
     $        /7x,'read from the input file. An option choice should',
     $        ' be chosen by',/7x,'placing a "*", "x", or "X" in the',
     $        ' checkbox ("[ ]"). This',/7x,'choice line belongs to',
     $        ' the option switch whose title string is',/7x,'"',a,'".')
              go to 990
            endif
          endif
c
          if (qmark) then
            nmark = nmark + 1
            iodb(n) = ival
            jlast = j
          endif
c
        enddo
  920   continue
c
        if (nmark .le. 0) then
          j2 = ilnobl(uodbtx(n))
          write (nttyo,1980) uodbtx(n)(1:j2)
 1980     format(/' * Warning - (XCON3/rd3d8b) No option choice was',
     $    ' checked on the input file',/7x,'for the iodb option',
     $    ' switch whose title string is',/7x,'"',a,'".')
c
          do j = 1,jdbxpa
            ival = iodbox(j,n)
            if (ival .eq. 0) go to 930
          enddo
c
          write (nttyo,1990)
 1990     format(/7x,'A default value of 0 has been applied, but no',
     $    /7x,'matching option choice string is defined.')
          go to 940
c
  930     j3 = ilnobl(uodbox(j,n))
          write (nttyo,1992) uodbox(j,n)(1:j3)
 1992     format(/7x,'A default value of 0 has been applied. The',
     $    ' matching string is',/7x,'"',a,'".')
  940     continue
        endif
c
        if (nmark .gt. 1) then
          j2 = ilnobl(uodbtx(n))
          j = jlast
          j3 = ilnobl(uodbox(j,n))
           write (nttyo,1994) uodbtx(n)(1:j2),uodbox(j,n)(1:j3)
 1994     format(/' * Warning - (XCON3/rd3d8b) More than one option',
     $    ' choice was checked',/7x,'on the input file for the iodb',
     $    ' option switch whose title string is',/7x,'"',a,'".',
     $    /7x,'The last choice checked will be used. The',
     $    ' matching string is',/7x,'"',a,'".')
        endif
c
      enddo
  950 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Numerical parameters.
c
c     Read the block title from a two-line header.
c
      uheadx = 'Numerical parameters'
      nfldtx = 1
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c     Read the beta convergence tolerance (tolbt) from a one-line
c     header.
c
      uheadx = 'Beta convergence tolerance'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      tolbt = var
c
c     Read the del convergence tolerance (toldl) from a one-line header.
c
      uheadx = 'Del convergence tolerance'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      toldl = var
c
c     Read the maximum number of Newton-Raphson iterations (itermx)
c     from a two-line header (the separator line is the end of the
c     current block).
c
      uheadx = 'Max. Number of N-R Iterations'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 999
      itermx = ivar
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Ordinary basis switches.
c
c     Read the block title from a two-line header.
c
      uheadx = 'Ordinary Basis Switches'
      nfldtx = 2
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
      nobsw = 0
      do n = 1,nbtmax + 1
c
c       Read a line. If the block has not been completely read,
c       this contains the name of a species to "Replace", and a
c       sub-block for that species follows. Otherwise, this line is
c       the first line of the next block.
c
        nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        ustr = ufield(1)
        uheadx = 'Replace'
        call locase(ustr)
        call locase(uheadx)
        j2 = ilnobl(ustr)
        j3 = ilnobl(uheadx)
c
        if (ustr(1:j2) .ne. uheadx(1:j3)) then
c
c         Back up.
c
          backspace ninpts
          go to 525
        endif
c
        ustr = ufield(2)
        ustrn = ustr
        call locase(ustrn)
        if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') ustr = 'None'
        qnone = ustr(1:5).eq.'None '
c
        if (.not.qnone) then
          nobsw = nobsw + 1
          uobsw(1,nobsw) = ufield(2)
        endif
c
c       Read the name of the "with" species from a two-line header.
c
        uheadx = 'with'
        nfldtx = 3
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
        if (qrderr) go to 999
        if (.not.qnone) uobsw(2,nobsw) = ufield(2)
      enddo
  525 nobswt = nobsw
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Saturation flag tolerance (tolspf).
c
      tolspf = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Sat. flag tolerance'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      tolspf = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Scale factor for the mass of aqueous solution to write
c     on the pickup file.
c
      scamas = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Aq. Phase Scale Factor'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      scamas = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read a two-line header marking the end of the current problem
c     input.
c
      uheadx = 'End of problem'
      nfldtx = 1
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  990 qrderr = .true.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
