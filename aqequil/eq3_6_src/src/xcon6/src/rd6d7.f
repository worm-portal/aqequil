      subroutine rd6d7(cdac,cesrb,csigma,dlzmx1,dlzmx2,dlzidp,
     $ dzprlg,dzprnt,eact,electr,fk,iact,iktbt,iktmax,imchmx,
     $ imech,iodb,iopg,iopr,iopt,ioscan,itermx,jcode,jreac,jtemp,
     $ jxmod,kct,kdim,kmax,kmt,kprs,ksq,ksppmx,kstpmx,kxmod,kxt,
     $ hact,modr,moffg,morr,mprs,mteaqb,mteb,nctmax,ndact,ndctmx,
     $ nesrbt,nffg,nffgmx,ninpts,nmodl1,nmodl2,nodbmx,nopgmx,
     $ noprmx,noptmx,nordlm,npslmx,nprmn,nprmx,nprsmx,nrct,nrctmx,
     $ nrk,nsk,nsrtmx,nsscmx,nsslmx,ntitl1,ntitl2,ntitmx,ntrymx,
     $ nttkmx,nttyo,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrtmx,
     $ qend,qrderr,rk0,rxbarb,sscrew,sk,tempci,tempc0,timemx,tolbt,
     $ toldl,tolsat,tolsst,tolx,trk0,tstrt,ttk,udac,uelemb,uendb,
     $ uesrb,uffg,undms,unrms,uprs,ureac,utitl1,utitl2,uxct16,uxmd24,
     $ uxopex,uxopt,vreac,xlkffg,xlkmod,zimax,zistrt,zkfac,zklogl,
     $ zklogu,zvclgi)
c
c     This subroutine reads the EQ6 input file in menu-style ("D")
c     format for versions 7.0-7.2. It thus encompasses two version
c     levels, '7.0' and '7.2'. The input file in this format is not
c     identical for these two version levels. However, the line
c     parsing capability can handle the differences, which consist
c     only of differences in field sizes.
c
c     This subroutine is called by:
c
c       XCON6/xcon6.f
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer iktmax,imchmx,kmax,nctmax,ndctmx,nffgmx,nodbmx,
     $ nopgmx,noprmx,noptmx,nprsmx,nrctmx,nsrtmx,nsscmx,ntitmx,
     $ nttkmx,nxmdmx,nxopmx,nxpemx,nxrtmx
c
      integer iact(imchmx,2,nrctmx),iktbt(nxrtmx),imech(2,nrctmx),
     $ iodb(nodbmx),iopg(nopgmx),iopr(noprmx),iopt(noptmx),
     $ jcode(nrctmx),jreac(nrctmx),jxmod(nxmdmx),kxmod(nxmdmx),
     $ ndact(imchmx,2,nrctmx),nesrbt(nsrtmx),nrk(2,nrctmx),nsk(nrctmx)
c
      integer ioscan,itermx,jtemp,kct,kdim,ksq,kmt,kprs,ksppmx,kstpmx,
     $ kxt,nffg,ninpts,nmodl1,nmodl2,nordlm,npslmx,nprmn,nprmx,nrct,
     $ nsslmx,ntitl1,ntitl2,ntrymx,nttyo,nxmod,nxopex,nxopt
c
      logical qend,qrderr
c
      character*80 utitl1(ntitmx),utitl2(ntitmx)
      character*48 uprs(nprsmx)
      character*24 udac(ndctmx,imchmx,2,nrctmx),uendb(iktmax,nxrtmx),
     $ uffg(nffgmx),undms(kmax),unrms(kmax),ureac(nrctmx),
     $ uxmd24(nxmdmx),uxopex(nxpemx)
      character*16 uxct16(nxopmx)
      character*8 uelemb(nctmax),uesrb(nctmax,nsrtmx),uxopt(nxopmx)
c
      real*8 cdac(ndctmx,imchmx,2,nrctmx),cesrb(nctmax,nsrtmx),
     $ csigma(imchmx,2,nrctmx),eact(imchmx,2,nrctmx),fk(nrctmx),
     $ hact(imchmx,2,nrctmx),modr(nrctmx),moffg(nffgmx),morr(nrctmx),
     $ mprs(nprsmx),mteaqb(nctmax),mteb(nctmax),rk0(imchmx,2,nrctmx),
     $ rxbarb(iktmax,nxrtmx),sscrew(nsscmx),sk(nrctmx),ttk(nttkmx),
     $ trk0(imchmx,2,nrctmx),vreac(nrctmx),xlkffg(nffgmx),
     $ xlkmod(nxmdmx),zvclgi(kmax)
c
      real*8 dlzmx1,dlzmx2,dlzidp,dzprlg,dzprnt,electr,tempci,tempc0,
     $ timemx,tolbt,toldl,tolsat,tolsst,tolx,tstrt,zkfac,zklogl,zklogu,
     $ zimax,zistrt
c
c-----------------------------------------------------------------------
c
      include 'xcon6/x6op7.h'
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
      integer i,im,idesc,idescx,iktb,ivar,j,jdesc,j1,j2,ksb,
     $ n,ncb,nfldmx,nfldt,nfldtx,nlchmx,nmark,nmarks,nrc,
     $ nsrt,nxrt
c
      integer ilnobl
c
      character*(nlchpa) ufield(nfldpa)
      character*(nlchpa) uline1,uline2,ulscr
      character*(nlchpa) uheadr,uheadx
c
      character*80 ustr
      character*24 ux
      character*1 ux1
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
        write (nttyo,1010) uline1
 1010   format(/' * Error - (XCON6/rd6d7) The first line of a "D"',
     $  /7x,'format input file must begin with "|--------".',
     $  /7x,'The first line read from the old input file begins with-',
     $  /7x,'"',a70,'".')
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
        utitl1(n) = ufield(1)
  110 continue
c
      write (nttyo,1015) ntitmx
 1015 format(/' * Error - (XCON6/rd6d7) Have too many lines in the',
     $ /7x,'main title. The code is only dimensioned for ',i4,
     $ /7x,'lines. Reduce the size of the title or increase the',
     $ /7x,'dimensioning parameter ntitpa.')
      go to 990
c
  120 ntitl1 = n - 1
c
c-----------------------------------------------------------------------
c
c     Nmodl option switches.
c
      nmodl2 = 0
      nmodl1 = 0
c
      uheadx = 'calculational mode'
      nfldtx = 4
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      call gmarko(nfldmx,nfldt,nmark,nttyo,ufield,uline1)
      if (nmark .gt. 0) nmodl2 = nmark - 2
c
      uheadx = 'model type'
      nfldtx = 4
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      call gmarko(nfldmx,nfldt,nmark,nttyo,ufield,uline1)
      if (nmark .gt. 0) nmodl1 = nmark - 1
c
c-----------------------------------------------------------------------
c
c     Temperature parameters.
c     Note: ttk(1) = tk1, etc.
c
      jtemp = 0
c
      uheadx = 'temperature model'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      call gmarko(nfldmx,nfldt,nmark,nttyo,ufield,uline1)
      if (nmark .gt. 0) jtemp = nmark - 2
c
      uheadx = 'tstart(c)'
      nfldtx = 8
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      tempc0 = var
      do 130 i = 1,3
        ustr = ufield(2*i + 2)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 990
        ttk(i) = var
  130 continue
c
c-----------------------------------------------------------------------
c
c     Zi and time parameters.
c
c     Note: cplim does not appear on the "D" format input file.
c
      uheadx = 'starting value of zi'
      nfldtx = 4
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      zistrt = var
      ustr = ufield(4)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      zimax = var
c
      uheadx = 'starting time (sec)'
      nfldtx = 4
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      tstrt = var
      ustr = ufield(4)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      timemx = var
c
      uheadx = 'max. steps'
      nfldtx = 4
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 990
      kstpmx = ivar
      ustr = ufield(4)
      call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 990
      ksppmx = ivar
c
c-----------------------------------------------------------------------
c
c     Print interval parameters.
c
      uheadx = 'linear print interval'
      nfldtx = 4
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      dzprnt = var
      ustr = ufield(4)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      dzprlg = var
c
c-----------------------------------------------------------------------
c
c     Plot interval parameters.
c
c     Note: dzplot, dzpllg, and ksplmx currently do not appear on the
c     "D" format input file.
c
c-----------------------------------------------------------------------
c
c     Ifile.
c
c     Note: ifile does not appear on the "D" format input file.
c
c-----------------------------------------------------------------------
c
c     Nxopt mineral subset selection suppression options.
c
      nxopt = 0
c
      uheadx = 'suppress mineral phases'
      nfldtx = 1
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
      uheadx = 'phases w/ elements'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
c
  140 do 150 n = 2,3
        ustr = ufield(n)
        call locase(ustr)
        if (ustr(1:8) .ne. '        ') then
          nxopt = nxopt + 1
c
          if (nxopt .gt. nxopmx) then
            write (nttyo,1025) nxopmx
 1025       format(/' * Error - (XCON6/rd6d7) Have too many mineral',
     $      /7x,'subset-selection suppression options. The code is',
     $      /7x,'only dimensioned for ',i3,' such options. Reduce the',
     $      /7x,'number of options or increase the dimensioning',
     $      /7x,'parameter nxoppa.')
            go to 990
          endif
c
          if (ustr(1:8) .eq. 'all     ') then
            uxopt(nxopt) = 'all'
            uxct16(nxopt) = '        '
          else
            uxopt(nxopt) = 'alwith'
            uxct16(nxopt) = ustr
          endif
        endif
  150 continue
c
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      if (ustr(1:8) .eq. '--------') go to 175
      uheadr = ufield(1)
      call locase(uheadr)
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .gt. 0) then
        nfldtx = 3
        if (nfldt .ne. nfldtx) then
          write (nttyo,1030) nfldt,nfldtx,uline1
 1030     format(/' * Warning - (XCON6/rd6d7) Found ',i2,' fields',
     $    /7x,'where ',i2,' were expected on the line beginning with-',
     $    /7x,'"',a70,'".')
        endif
        go to 140
      endif
c
      nxopex = 0
c
      uheadx = 'phases except'
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .eq. 0) then
        write (nttyo,1020) uheadx,uline1
 1020   format(/' * Error - (XCON6/rd6d7) Was expecting to find the',
     $  /7x,'header beginning with-',/7x,'"',a70,'"',
     $  /7x,'on the line beginning with-',/7x,'"',a70,'".')
        go to 990
      endif
c
  160 do 170 n = 2,3
        ustr = ufield(n)
        if (ustr(1:8) .ne. '        ') then
          nxopex = nxopex + 1
c
          if (nxopex .gt. nxpemx) then
            write (nttyo,1027) nxpemx
 1027       format(/' * Error - (XCON6/rd6d7) Have too many',
     $      /7x,'exceptions specified to the mineral subset-selection',
     $      /7x,'suppression options. The code is only dimensioned',
     $      /7x,'for ',i3,'exceptions. Reduce the number of exceptions',
     $      /7x,'or increase the dimensioning parameter nxoppa.')
            go to 990
          endif
c
          uxopex(nxopex) = ustr
        endif
  170 continue
c
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      uheadr = ufield(1)
      call locase(uheadr)
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .gt. 0) then
        nfldtx = 3
        if (nfldt .ne. nfldtx) write (nttyo,1030) nfldt,nfldtx,uline1
        go to 160
      endif
c
      ustr = ufield(1)
      if (ustr(1:8) .ne. '--------') then
        write (nttyo,1040) uline1
 1040   format(/' * Error - (XCON6/rd6d7) Found the line beginning',
     $  ' with-',/7x,'"',a70,'"',
     $  /7x,'where a dashed separator line was expected.')
        go to 990
      endif
c
  175 continue
c
c-----------------------------------------------------------------------
c
c     Nffg options.
c
      nffg = 0
c
      uheadx = 'fixed fugacity phases- species, '
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
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
c
      ustr = ufield(1)
      call locase(ustr)
      if (ustr(1:4) .eq. 'none') go to 190
c
      nfldtx = 3
      if (nfldt .ne. nfldtx) write (nttyo,1030) nfldt,nfldtx,uline1
c
  180 nffg = nffg + 1
c
      if (nffg .gt. nffgmx) then
        write (nttyo,1035) nffgmx
 1035   format(/' * Error - (XCON6/rd6d7) Have too many gases whose',
     $  /7x,'fugacities are to be fixed. The code is only dimensioned',
     $  /7x,'for ',i4,' such gases. Reduce the number of gases or',
     $  /7x,'increase the dimensioning parameter nffgpa.')
        go to 990
      endif
c
      uffg(nffg) = ufield(1)
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      moffg(nffg) = var
      ustr = ufield(3)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      xlkffg(nffg) = var
c
  190 nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      if (ustr(1:8) .eq. '--------') go to 200
      nfldtx = 3
      if (nfldt .ne. nfldtx) write (nttyo,1030) nfldt,nfldtx,uline1
      go to 180
c
c     Read a second dashed separator line.
c
  200 nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      if (ustr(1:8) .ne. '--------') then
        write (nttyo,1040) uline1
        go to 990
      endif
c
c-----------------------------------------------------------------------
c
c     Reactants.
c
      nrct = 0
      nsrt = 0
      nxrt = 0
c
      uheadx = 'reactants'
      nfldtx = 1
      call rdd2l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,uline2,ulscr)
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
  210 nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
c
      ustr = ufield(1)
      if (ustr(1:8).eq.'--------') go to 330
c
      uheadr = ufield(1)
      call locase(uheadr)
c
      uheadx = 'REACTANT '
      call locase(uheadx)
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .eq. 0) then
        write (nttyo,1020) uheadx,uline1
        go to 990
      endif
c
      nfldtx = 4
      if (nfldt .ne. nfldtx) write (nttyo,1030) nfldt,nfldtx,uline1
c
      ustr = ufield(2)
      call locase(ustr)
      if (ustr(1:4) .eq. 'none') then
        nfldtx = 1
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        ustr = ufield(1)
        if (ustr(1:8) .ne. '--------') then
          write (nttyo,1040) uline1
          go to 990
        endif
        go to 330
      endif
c
  215 nrc = nrc + 1
c
      if (nrc .gt. nrctmx) then
        write (nttyo,1037) nrctmx
 1037   format(/' * Error - (XCON6/rd6d7) Have too many reactants',
     $  /7x,'The code is only dimensioned for ',i4,' reactants.',
     $  /7x,'Reduce the number of reactants or increase the',
     $  /7x,'dimensioning parameter nrctpa.')
        go to 990
      endif
c
      nrct = nrc
      ureac(nrc) = ufield(2)
      ustr = ufield(4)
      call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 990
      jreac(nrc) = ivar
c
      uheadx = 'moles remaining'
      nfldtx = 4
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      morr(nrc) = var
      ustr = ufield(4)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      modr(nrc) = var
c
      uheadx = 'reactant type'
      nfldtx = 4
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      if (ustr(1:7) .eq. 'mineral') then
        jcode(nrc) = 0
      elseif (ustr(1:14) .eq. 'solid solution') then
        jcode(nrc) = 1
      elseif (ustr(1:7) .eq. 'special') then
        jcode(nrc) = 2
      elseif (ustr(1:7) .eq. 'aqueous') then
        jcode(nrc) = 3
      elseif (ustr(1:3) .eq. 'gas') then
        jcode(nrc) = 4
      else
        jcode(nrc) = 0
      endif
      ustr = ufield(4)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      sk(nrc) = var
c
      uheadx = 'surface type'
      nfldtx = 4
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 990
      nsk(nrc) = ivar
      ustr = ufield(4)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      fk(nrc) = var
c
      if (jcode(nrc) .eq. 1) then
        nxrt = nxrt + 1
c
        if (nxrt .gt. nxrtmx) then
          write (nttyo,1038) nxrtmx
 1038     format(/' * Error - (XCON6/rd6d7) Have too many solid',
     $    /7x,'solution reactants. The code is only dimensioned',
     $    /7x,'for ',i4,' such reactants. Reduce the number of such',
     $    /7x,'reactants or increase the dimensioning parameter',
     $    ' nxrtpa.')
          go to 990
        endif
c
        uheadx = 'end-member'
        nfldtx = 4
        call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
        if (qrderr) go to 999
c
        iktb = 0
  220   iktb = iktb + 1
c
        if (iktb .gt. iktmax) then
          write (nttyo,1039) ureac(nrc),iktmax
 1039     format(/' * Error - (XCON6/rd6d7) Have too many end-members',
     $    /7x,'in the solid solution reactant "',a24,'".',
     $    /7x,'The code is only dimensioned for ',i4,' end-members per',
     $    /7x,'solid solution. Reduce the number of end-members or',
     $    /7x,'increase the dimensioning parameter iktpar.')
          go to 990
        endif
c
        uendb(iktb,nxrt) = ufield(2)
        ustr = ufield(4)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 990
        rxbarb(iktb,nxrt) = var
c
        nfldtx = 4
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        uheadr = ufield(1)
        call locase(uheadr)
        j2 = ilnobl(uheadx)
        j = index(uheadr,uheadx(1:j2))
        if (j .gt. 0) go to 220
        iktbt(nxrt) = iktb
      else
        uheadx = 'end-member'
        nfldtx = 4
  230   call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        uheadr = ufield(1)
        call locase(uheadr)
        j2 = ilnobl(uheadx)
        j = index(uheadr,uheadx(1:j2))
        if (j .gt. 0) go to 230
      endif
c
      if (jcode(nrc) .eq. 2) then
        nsrt = nsrt + 1
c
        if (nsrt .gt. nsrtmx) then
          write (nttyo,1042) nsrtmx
 1042     format(/' * Error - (XCON6/rd6d7) Have too many special',
     $    /7x,'reactants. The code is only dimensioned for ',i4,
     $    /7x,'such reactants. Reduce the number of such reactants',
     $    /7x,'or increase the dimensioning parameter nsrtpa.')
          go to 990
        endif
c
        nrct = nrc
        uheadx = 'volume'
        j2 = ilnobl(uheadx)
        j = index(uheadr,uheadx(1:j2))
        if (j .eq. 0) then
          write (nttyo,1020) uheadx,uline1
          go to 990
        endif
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 990
        vreac(nsrt) = var
c
        uheadx = 'element'
        nfldtx = 4
        call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
        if (qrderr) go to 999
c
        ncb = 0
  240   ncb = ncb + 1
c
        if (ncb .gt. nctmax) then
          write (nttyo,1044) ureac(nrc),nctmax
 1044     format(/' * Error - (XCON6/rd6d7) Have too many chemical',
     $    /7x,'elements in the special  reactant "',a24,'".',
     $    /7x,'The code is only dimensioned for ',i4,' elements.',
     $    /7x,'Reduce the number of elements or increase the',
     $    /7x,'dimensioning parameter nctpar.')
          go to 990
        endif
c
        uesrb(ncb,nsrt) = ufield(2)
        ustr = ufield(4)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 990
        cesrb(ncb,nsrt) = var
c
        nfldtx = 4
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        uheadr = ufield(1)
        call locase(uheadr)
        j2 = ilnobl(uheadx)
        j = index(uheadr,uheadx(1:j2))
        if (j .gt. 0) go to 240
        nesrbt(nsrt) = ncb
      else
        uheadx = 'volume'
        j2 = ilnobl(uheadx)
        j = index(uheadr,uheadx(1:j2))
c
        uheadx = 'element'
        nfldtx = 4
  250   call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        uheadr = ufield(1)
        call locase(uheadr)
        j2 = ilnobl(uheadx)
        j = index(uheadr,uheadx(1:j2))
        if (j .gt. 0) go to 250
      endif
c
      uheadx = 'DISSOLUTION LAW'
      call locase(uheadx)
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .eq. 0) then
        write (nttyo,1020) uheadx,uline1
        go to 990
      endif
      ustr = ufield(2)
      call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 990
      nrk(1,nrc) = ivar
c
c     Rate law parameters, forward direction.
c
      nfldtx = 4
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
c
      if (nrk(1,nrc) .le. 0) then
        imech(1,nrc) = 0
        go to 290
      endif
c
      im = 0
  260 im = im + 1
c
      uheadr = ufield(1)
      call locase(uheadr)
c
      uheadx = 'PRECIPITATION LAW'
      call locase(uheadx)
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .gt. 0) then
        im = im - 1
        go to 290
      endif
c
      if (nrk(1,nrc) .le. 0) then
        write (nttyo,1020) uheadx,uline1
        go to 990
      endif
c
      uheadx = 'rate constant rk'
      write (ux1,'(i1)') im
      uheadx(17:17) = ux1(1:1)
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .eq. 0) then
        write (nttyo,1020) uheadx,uline1
        go to 990
      endif
c
      if (im .gt. imchmx) then
        write (nttyo,1045) ureac(nrc),imchmx
 1045   format(/' * Error - (XCON6/rd6d7) Have too many rate',
     $  /7x,'constants in the forward direction rate law for',
     $  /7x,'reactant "',a24,'". The code is only',
     $  /7x,'dimensioned for ',i2,' rate constants per rate law.',
     $  /7x,'Reduce the number of rate constants or increase the',
     $  /7x,'dimensioning parameter imchpa.')
        go to 990
      endif
c
      i = im
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      rk0(i,1,nrc) = var
      if (nrk(1,nrc) .eq. 2) then
        ustr = ufield(4)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 990
        csigma(i,1,nrc) = var
      endif
c
      if (nrk(1,nrc) .eq. 1) go to 280
c
      n = 0
      uheadx = 'aqueous species'
      nfldtx = 4
  270 call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      uheadr = ufield(1)
      call locase(uheadr)
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .gt. 0) then
        n = n + 1
c
        if (n .gt. ndctmx) then
          write (nttyo,1046) ureac(nrc),i,ndctmx
 1046     format(/' * Error - (XCON6/rd6d7) Have too many species',
     $    /7x,'in the activity product in term ',i2,' of the',
     $    /7x,'forward direction rate law for reactant',
     $    /7x,'"',a24,'". The code is only dimensioned for ',i3,
     $    /7x,'such species. Reduce the number of such species',
     $    /7x,'or increase the dimensioning parameter ndctpa.')
        go to 990
        endif
c
        udac(n,i,1,nrc) = ufield(2)
        ustr = ufield(4)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 990
        cdac(n,i,1,nrc) = var
        go to 270
      endif
      ndact(i,1,nrc) = n
c
      uheadx = 'PRECIPITATION LAW'
      call locase(uheadx)
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .gt. 0) go to 290
c
      uheadx = 'temperature (c)'
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .eq. 0) then
        write (nttyo,1020) uheadx,uline1
        go to 990
      endif
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      trk0(i,1,nrc) = var
c
      iact(i,1,nrc) = 0
c
      nfldtx = 4
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      uheadr = ufield(1)
      call locase(uheadr)
c
      uheadx = 'PRECIPITATION LAW'
      call locase(uheadx)
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .gt. 0) go to 290
c
      uheadx = 'rate constant rk'
      call locase(uheadx)
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .gt. 0) go to 260
c
      uheadx = 'act. energy-kcal'
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .eq. 0) then
        write (nttyo,1020) uheadx,uline1
        go to 990
      endif
c
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      eact(i,1,nrc) = var
      if (var .gt. 0.) iact(i,1,nrc) = 1
      ustr = ufield(4)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      hact(i,1,nrc) = var
      if (var .gt. 0.) iact(i,1,nrc) = 2
c
  280 nfldtx = 4
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      go to 260
c
  290 imech(1,nrc) = im
c
c     Rate law parameters, backward direction.
c
      ustr = ufield(2)
      call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 990
      nrk(2,nrc) = ivar
c
      if (nrk(2,nrc) .le. 0) then
        imech(2,nrc) = 0
        go to 210
      endif
c
      im = 0
 300  im = im + 1
c
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      if (ustr(1:8) .eq. '--------') then
        imech(2,nrc) = 0
        go to 330
      endif
c
      nfldtx = 4
      if (nfldt .ne. nfldtx) write (nttyo,1030) nfldt,nfldtx,uline1
c
      uheadr = ufield(1)
      call locase(uheadr)
c
      uheadx = 'REACTANT'
      call locase(uheadx)
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .gt. 0) then
        im = im - 1
        imech(2,nrc) = im
        go to 215
      endif
c
      uheadx = 'rate constant rk'
      write (ux1,'(i1)') im
      uheadx(17:17) = ux1(1:1)
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .eq. 0) then
        write (nttyo,1020) uheadx,uline1
        go to 990
      endif
c
      if (im .gt. imchmx) then
        write (nttyo,1047) ureac(nrc),imchmx
 1047   format(/' * Error - (XCON6/rd6d7) Have too many rate',
     $  /7x,'constants in the backward direction rate law for',
     $  /7x,'reactant "',a24,'". The code is only',
     $  /7x,'dimensioned for ',i2,' rate constants per rate law.',
     $  /7x,'Reduce the number of rate constants or increase the',
     $  /7x,'dimensioning parameter imchpa.')
        go to 990
      endif
c
      i = im
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      rk0(i,2,nrc) = var
      if (nrk(2,nrc) .eq. 2) then
        ustr = ufield(4)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 990
        csigma(i,2,nrc) = var
      endif
c
      if (nrk(2,nrc) .eq. 1) go to 210
c
      n = 0
      uheadx = 'aqueous species'
      nfldtx = 4
  310 call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      uheadr = ufield(1)
      call locase(uheadr)
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .gt. 0) then
        n = n + 1
c
        if (n .gt. ndctmx) then
          write (nttyo,1048) ureac(nrc),i,ndctmx
 1048     format(/' * Error - (XCON6/rd6d7) Have too many species',
     $    /7x,'in the activity product in term ',i2,' of the',
     $    /7x,'forward direction rate law for reactant',
     $    /7x,'"',a24,'". The code is only dimensioned for ',i3,
     $    /7x,'such species. Reduce the number of such species',
     $    /7x,'or increase the dimensioning parameter ndctpa.')
          go to 990
        endif
c
        udac(n,i,2,nrc) = ufield(2)
        ustr = ufield(4)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 990
        cdac(n,i,2,nrc) = var
        go to 310
      endif
      ndact(i,2,nrc) = n
c
      ustr = ufield(1)
      if (ustr(1:8) .eq. '--------') then
        imech(2,nrc) = im
        go to 330
      endif
c
      uheadx = 'temperature (c)'
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .eq. 0) then
        write (nttyo,1020) uheadx,uline1
        go to 990
      endif
c
      nfldtx = 4
      if (nfldt .ne. nfldtx) write (nttyo,1030) nfldt,nfldtx,uline1
c
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      trk0(i,2,nrc) = var
c
      iact(i,2,nrc) = 0
c
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
c
      ustr = ufield(1)
      if (ustr(1:8) .eq. '--------') then
        imech(2,nrc) = im
        go to 330
      endif
c
      uheadr = ufield(1)
      call locase(uheadr)
c
      uheadx = 'REACTANT'
      call locase(uheadx)
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .gt. 0) then
        imech(2,nrc) = im
        go to 215
      endif
c
      uheadx = 'rate constant rk'
      call locase(uheadx)
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .gt. 0) go to 300
c
      uheadx = 'act. energy-kcal'
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .eq. 0) then
        write (nttyo,1020) uheadx,uline1
        go to 990
      endif
c
      nfldtx = 4
      if (nfldt .ne. nfldtx) write (nttyo,1030) nfldt,nfldtx,uline1
c
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      eact(i,2,nrc) = var
      if (var .gt. 0.) iact(i,2,nrc) = 1
      ustr = ufield(4)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      hact(i,2,nrc) = var
      if (var .gt. 0.) iact(i,2,nrc) = 2
c
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      go to 300
c
  330 continue
c
c-----------------------------------------------------------------------
c
c     Options. These are iopt, iopr, and iodb option switches,
c     skipping ones which may be classified as development options.
c     The development options are read below, in a somewhat different
c     manner. The iopg option switches are read still further below,
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
        write (nttyo,1050) uline1
 1050   format(/' * Error - (XCON6/rd6d7) Found the line beginning',
     $  ' with-',/7x,'"',a70,'"',
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
      do 370 idescx = 1,nop6pa
        uheadx = uopt6(idescx)
        call locase(uheadx)
        if (uheadr(1:40) .eq. uheadx(1:40)) then
          if (idescx .ne. i) then
            write (nttyo,1060) uopt6(idescx),uvar6(idescx),
     $      index6(idescx),i,idescx
 1060       format(/' * Warning - (XCON6/rd6d7) The input for the',
     $      /7x,'option whose identifying string which begins with-',
     $      /7x,'"',a70,'"',
     $      /7x,'(',a4,'(',i2,')) is out of order, in place ',i3,
     $      ' instead of place ',i3,'.')
          endif
          go to 380
        endif
  370 continue
c
      write (nttyo,1070) uheadr
 1070 format(/" * Error - (XCON6/rd6d7) Can't identify the following",
     $ /7x,'iopt, iopg, iopr, or iodb option header-',
     $ /7x,'"',a70,'".')
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
          write (nttyo,1080) uopt6(idescx),uvar6(idescx),index6(idescx)
 1080     format(/' * Error - (XCON6/rd6d7) More than one choice',
     $    /7x,'was marked for the following option-',
     $    /7x,'"',a70,'"',/9x,a4,'(',i2,').')
          go to 990
        endif
        go to 420
      endif
c
c     Is the current line a new option header?
c
      if (ulscr(1:2) .eq. '- ') then
        if (nmarks .gt. 1) then
          write (nttyo,1080) uopt6(idescx),uvar6(idescx),index6(idescx)
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
      do 400 idesc = 1,nod6pa
        if (iopti6(idesc) .eq. idescx) then
          uheadx = udesc6(idesc)
          call locase(uheadx)
          if (uheadr(1:40) .eq. uheadx(1:40)) go to 410
        endif
  400 continue
c
      write (nttyo,1090) uheadr,uopt6(idescx),uvar6(idescx),
     $ index6(idescx)
 1090 format(/" * Error - (XCON6/rd6d7) Can't identify the following",
     $ /7x,'option string-',/7x,'"',a70,'".',
     $ /7x,'It was given for the following option-',
     $ /7x,'"',a70,'"',/9x,a4,'(',i2,').')
      go to 990
c
  410 if (nmark .eq. 1) then
        ustr = uvar6(idescx)
        call locase(ustr)
        if (ustr(1:4) .eq. 'iopt') then
          iopt(index6(idescx)) = ivalu6(idesc)
        elseif (ustr(1:4) .eq. 'iopg') then
          iopg(index6(idescx)) = ivalu6(idesc)
        elseif (ustr(1:4) .eq. 'iopr') then
          iopr(index6(idescx)) = ivalu6(idesc)
        elseif (ustr(1:4) .eq. 'iodb') then
          iodb(index6(idescx)) = ivalu6(idesc)
        else
          write (nttyo,1100) uheadr,uvar6(idescx)
 1100     format(/" * Error - (XCON6/rd6d7) Don't recognize the",
     $    /7x,'option type string "',a4,'", which was given for',
     $    /7x,'the following option-',/7x,'"',a70,'".')
          go to 990
        endif
      endif
c
      go to 390
c
  420 continue
c
c     Get the time frame flag (iopt(1)).
c
      iopt(1) = 0
      if (nrct. gt. 0) then
        do 430 nrc = 1,nrct
          if (nrk(1,nrc) .ge. 2) then
            iopt(1) = 1
            go to 435
          endif
          if (nrk(2,nrc) .ge. 2) then
            iopt(1) = 1
            go to 435
          endif
  430   continue
  435   continue
      endif
c
c-----------------------------------------------------------------------
c
c     Development options.
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
      do 450 idescx = 1,ndv6pa
        uheadx = udevl6(idescx)
        call locase(uheadx)
        if (uheadr(1:40) .eq. uheadx(1:40)) go to 460
  450 continue
c
      write (nttyo,1110) uheadr
 1110 format(/" * Error - (XCON6/rd6d7) Can't identify the following",
     $ /7x,'development option header-',
     $ /7x,'"',a70,'".')
      go to 990
c
  460 call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 990
      ustr = udv6vr(idescx)
      call locase(ustr)
      if (ustr(1:4) .eq. 'iopt') then
        iopt(index6(idescx)) = ivar
      elseif (ustr(1:4) .eq. 'iopr') then
        iopr(index6(idescx)) = ivar
      elseif (ustr(1:4) .eq. 'iodb') then
        iodb(index6(idescx)) = ivar
      elseif (ustr(1:4) .eq. 'iopg') then
        iopg(index6(idescx)) = ivar
      else
        write (nttyo,1120) uheadr,udv6vr(idescx)
 1120   format(/" * Error - (XCON6/rd6d7) Don't recognize the",
     $  /7x,'option type string "',a4,'", which was  given for',
     $  /7x,'the following option-',/7x,'"',a70,'".')
        go to 990
      endif
      go to 440
c
  470 continue
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
      do 515 idescx = 1,nto6pa
        uheadx = utol6(idescx)
        call locase(uheadx)
        if (uheadr(1:32) .eq. uheadx(1:32)) go to 517
  515 continue
c
      write (nttyo,1130) uheadr
 1130 format(/" * Error - (XCON6/rd6d7) Can't identify the following",
     $ /7x,'tolerance parameter header-',
     $ /7x,'"',a70,'".')
      go to 990
c
  517 ustr = ufield(2)
      if (idescx .eq. 1) then
c
c       Itermx.
c
        call chrint(ivar,nttyo,qrderr,ustr)
        if (qrderr) go to 999
        itermx = ivar
      elseif (idescx .eq. 2) then
c
c       Dlzidp.
c
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        dlzidp = var
      elseif (idescx .eq. 3) then
c
c       Tolbt.
c
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        tolbt = var
      elseif (idescx .eq. 4) then
c
c       Toldl.
c
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        toldl = var
      elseif (idescx .eq. 5) then
c
c       Tolx.
c
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        tolx = var
      elseif (idescx .eq. 6) then
c
c       Tolsat.
c
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        tolsat = var
      elseif (idescx .eq. 7) then
c
c       Tolsst.
c
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        tolsst = var
      elseif (idescx .eq. 8) then
c
c       Sscrew1.
c
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        sscrew(1) = var
      elseif (idescx .eq. 9) then
c
c       Sscrew2.
c
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        sscrew(2) = var
      elseif (idescx .eq. 10) then
c
c       Sscrew3.
c
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        sscrew(3) = var
      elseif (idescx .eq. 11) then
c
c       Sscrew4.
c
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        sscrew(4) = var
      elseif (idescx .eq. 12) then
c
c       Sscrew5.
c
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        sscrew(5) = var
      elseif (idescx .eq. 13) then
c
c       Sscrew6.
c
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        sscrew(6) = var
      elseif (idescx .eq. 14) then
c
c       Zklogu.
c
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        zklogu = var
      elseif (idescx .eq. 15) then
c
c       Zklogl.
c
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        zklogl = var
      elseif (idescx .eq. 16) then
c
c       Zkfac.
c
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        zkfac = var
      elseif (idescx .eq. 17) then
c
c       Dlzmx1.
c
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        dlzmx1 = var
      elseif (idescx .eq. 18) then
c
c       Dlzmx2.
c
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        dlzmx2 = var
      elseif (idescx .eq. 19) then
c
c       Nordlm.
c
        call chrint(ivar,nttyo,qrderr,ustr)
        if (qrderr) go to 999
        nordlm = ivar
      elseif (idescx .eq. 20) then
c
c       Ntrymx.
c
        call chrint(ivar,nttyo,qrderr,ustr)
        if (qrderr) go to 999
        ntrymx = ivar
      elseif (idescx .eq. 21) then
c
c       Npslmx.
c
        call chrint(ivar,nttyo,qrderr,ustr)
        if (qrderr) go to 999
        npslmx = ivar
      elseif (idescx .eq. 22) then
c
c       Nsslmx.
c
        call chrint(ivar,nttyo,qrderr,ustr)
        if (qrderr) go to 999
        nsslmx = ivar
      elseif (idescx .eq. 23) then
c
c       Ioscan.
c
        call chrint(ivar,nttyo,qrderr,ustr)
        if (qrderr) go to 999
        ioscan = ivar
      endif
c
      go to 510
c
  520 continue
c
c-----------------------------------------------------------------------
c
c     Process the bottom half of the current output file.
c
c     Title.
c
      n = 1
      do 530 n = 1,ntitmx + 1
        read (ninpts,1000,err=990) uline1
        call parsln(nfldmx,nfldt,nlchmx,ufield,uline1,ulscr)
        ustr = ufield(1)
        if (ustr(1:8) .eq. '--------') go to 540
        utitl2(n) = ufield(1)
  530 continue
c
        write (nttyo,1115) ntitmx
 1115   format(/' * Error - (XCON6/rd6d7) Have too many lines in the',
     $  /7x,'secondary title. The code is only dimensioned for ',i4,
     $  /7x,'lines. Reduce the size of the title or increase the',
     $  /7x,'dimensioning parameter ntitpa.')
        go to 990
c
  540 ntitl2 = n - 1
c
c-----------------------------------------------------------------------
c
c     Original temperature.
c
      uheadx = 'temperature (c)'
      nfldtx = 2
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      tempci = var
c
c-----------------------------------------------------------------------
c
c     Electrical imbalance.
c
      uheadx = 'electrical imbalance'
      nfldtx = 2
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      electr = var
c
c-----------------------------------------------------------------------
c
c     Number of aqueous basis species.
c
      uheadx = 'number of aqueous master species'
      nfldtx = 2
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 990
      ksb = ivar
c
      ksq = ksb
      kct = ksb - 1
c
c-----------------------------------------------------------------------
c
c     Postion of last pure mineral.
c
      uheadx = 'position of last pure mineral'
      nfldtx = 2
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 990
      kmt = ivar
c
c-----------------------------------------------------------------------
c
c     Postion of last solid solution.
c
      uheadx = 'position of last solid solution'
      nfldtx = 2
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 990
      kxt = ivar
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
      if (ustr(1:8) .eq. 'none') go to 560
c
      nfldtx = 4
      if (nfldt .ne. nfldtx) write (nttyo,1030) nfldt,nfldtx,uline1
c
  550 n = n + 1
c
      if (n .gt. nxmdmx) then
        write (nttyo,1132) nxmdmx
 1132   format(/' * Error - (XCON6/rd6d7) Have too many nxmod',
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
 1135   format(/" * Error - (XCON6/rd6d7) Can't identify the",
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
 1140   format(/" * Error - (XCON6/rd6d7) Can't identify the",
     $  /7x,'following alter/suppress option string- "',a,'".',
     $  /7x,'This must be one of "suppress", "replace", "augmentk",',
     $  /7x,' or "augmentg".')
        go to 990
      endif
c
      ustr = ufield(4)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
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
c     Iopg options.
c     Note: iopg(1) = iopg1, etc.
c
      uheadx = 'iopg options'
      nfldtx = 1
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
      i = 0
      nfldtx = 1
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ulscr = ufield(1)
      if (ulscr(1:8) .eq. '--------') go to 720
c
      if (ulscr(1:2) .ne. '- ') then
        write (nttyo,1050) uline1
        go to 990
      endif
c
c     Have found an option header.
c
  640 i = i + 1
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
      do 670 idescx = 1,nop6pa
        uheadx = uopt6(idescx)
        call locase(uheadx)
        if (uheadr(1:40) .eq. uheadx(1:40)) then
c
c         This warning test is commented out.  The warning doesn't
c         really apply as the version 7 options are structured
c         (the strings and such for the iopg options are mixed in
c         with those of the iopt, iopr, and iodb options).
c
c         if (idescx .ne. i) then
c           write (nttyo,1060) uopt6(idescx),uvar6(idescx),
c    $      index6(idescx),i,idescx
c         endif
          go to 680
        endif
  670 continue
c
      write (nttyo,1070) uheadr
      go to 990
c
  680 continue
c
c     Read the first line after an option  header.
c
      jdesc = 0
  690 jdesc = jdesc + 1
      nfldtx = 1
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ulscr = ufield(1)
      if (ulscr(1:8) .eq. '--------') then
        if (nmarks .gt. 1) then
          write (nttyo,1080) uopt6(idescx),uvar6(idescx),index6(idescx)
        endif
        go to 720
      endif
c
c     Is the current line a new option header?
c
      if (ulscr(1:2) .eq. '- ') then
        if (nmarks .gt. 1) then
          write (nttyo,1080) uopt6(idescx),uvar6(idescx),index6(idescx)
          go to 990
        endif
        go to 640
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
      do 700 idesc = 1,nod6pa
        if (iopti6(idesc) .eq. idescx) then
          uheadx = udesc6(idesc)
          call locase(uheadx)
          if (uheadr(1:40) .eq. uheadx(1:40)) go to 710
        endif
  700 continue
c
      write (nttyo,1090) uheadr,uopt6(idescx),uvar6(idescx),
     $ index6(idescx)
      go to 990
c
  710 if (nmark .eq. 1) then
        ustr = uvar6(idescx)
        call locase(ustr)
        if (ustr(1:4) .eq. 'iopt') then
          iopt(index6(idescx)) = ivalu6(idesc)
        elseif (ustr(1:4) .eq. 'iopg') then
          iopg(index6(idescx)) = ivalu6(idesc)
        elseif (ustr(1:4) .eq. 'iopr') then
          iopr(index6(idescx)) = ivalu6(idesc)
        elseif (ustr(1:4) .eq. 'iodb') then
          iodb(index6(idescx)) = ivalu6(idesc)
        else
          write (nttyo,1100) uheadr,uvar6(idescx)
          go to 990
        endif
      endif
c
      go to 690
c
  720 continue
c
c-----------------------------------------------------------------------
c
c     Balance totals.
c
      uheadx = 'elements, moles and moles aqueous'
      nfldtx = 1
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
      i = 0
  730 nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ulscr = ufield(1)
      if (ulscr(1:8) .eq. '--------') go to 750
c
      nfldtx = 3
      if (nfldt .ne. nfldtx) write (nttyo,1030) nfldt,nfldtx,uline1
c
      i = i + 1
c
      if (i .gt. nctmax) then
        write (nttyo,1134) nctmax
 1134   format(/' * Error - (XCON6/rd6d7) Have too many chemical',
     $  /7x,'elements present. The code is only dimensioned',
     $  /7x,'for ',i3,' elements. Reduce the number of elements',
     $  /7x,'or increase the dimensioning parameter nctpar.')
        go to 990
      endif
c
      uelemb(i) = ufield(1)
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      mteb(i) = var
      ustr = ufield(3)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      mteaqb(i) = var
      go to 730
c
  750 continue
c
c-----------------------------------------------------------------------
c
c     Basis variable data.
c
      uheadx = 'master species and logarithmic basis variables'
      nfldtx = 1
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
      i = 0
  760 nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ulscr = ufield(1)
      if (ulscr(1:8) .eq. '--------') go to 770
c
      nfldtx = 3
      if (nfldt .ne. nfldtx) write (nttyo,1030) nfldt,nfldtx,uline1
c
      i = i + 1
c
      if (i .gt. kmax) then
        write (nttyo,1137) kmax
 1137   format(/' * Error - (XCON6/rd6d7) Have too many master',
     $  /7x,'variables. The code is only dimensioned for ',i3,
     $  /7x,'master variables. Reduce the number of such variables',
     $  /7x,'or increase the dimensioning parameter kpar.')
        go to 990
      endif
c
      unrms(i) = ufield(1)
      undms(i) = ufield(2)
      ustr = ufield(3)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      zvclgi(i) = var
      go to 760
c
  770 kdim = i
c
c-----------------------------------------------------------------------
c
c     Physically removed system data.
c
      uheadx = 'physically removed subsystem'
      nfldtx = 1
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
      nprmn = 0
      nprmx = 0
      kprs = 0
c
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      if (ustr(1:8) .eq. '--------') go to 820
c
      if (ustr(1:4) .eq. 'none') then
        nfldtx = 1
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        ustr = ufield(1)
        if (ustr(1:8) .ne. '--------') then
          write (nttyo,1040) uline1
          go to 990
        endif
        go to 820
      endif
c
      nfldtx = 3
      if (nfldt .ne. nfldtx) write (nttyo,1030) nfldt,nfldtx,uline1
      ulscr = ufield(1)
      call locase(ulscr)
c
c     Pure minerals.
c
  780 ustr = ufield(1)
      if (ustr(1:8) .eq. '--------') go to 820
      nfldtx = 3
      if (nfldt .ne. nfldtx) write (nttyo,1030) nfldt,nfldtx,uline1
      ulscr = ufield(1)
      if (ulscr(1:8) .ne. '        ') go to 790
      nprmn = nprmn + 1
c
      if (nprmn .gt. nprsmx) then
        write (nttyo,1150) nprsmx
 1150   format(/' * Error - (XCON6/rd6d7) Have too many mineral',
     $  /7x,'species in the physically removed system. The code is',
     $  /7x,'only dimensioned for ',i3,' such species. Reduce the',
     $  /7x,'number of such species or increase the dimensioning',
     $  /7x,'parameter nprspa.')
        go to 990
      endif
c
      uprs(nprmn)(1:24) = ufield(2)
      uprs(nprmn)(25:48) = ufield(2)
      ustr = ufield(3)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      mprs(nprmn) = var
c
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      go to 780
c
c     Solid solutions.
c
  790 nprmx = nprmn
c
  800 ux = ufield(1)
c
  810 nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      if (ustr(1:8) .eq. '--------') go to 820
c
      nfldtx = 3
      if (nfldt .ne. nfldtx) write (nttyo,1030) nfldt,nfldtx,uline1
      ulscr = ufield(1)
      if (ulscr(1:8) .ne. '        ') go to 800
c
      nprmx = nprmx + 1
c
      if (nprmn .gt. nprsmx) then
        write (nttyo,1150) nprsmx
        go to 990
      endif
c
      uprs(nprmx)(1:24) = ux
      uprs(nprmx)(25:48) = ufield(2)
      ustr = ufield(3)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 990
      mprs(nprmx) = var
      go to 810
c
  820 n = nprmn + nprmx
      if (n .gt. 0) kprs = 1
c
c-----------------------------------------------------------------------
c
      go to 999
c
  990 qrderr = .true.
c
  999 continue
      end
