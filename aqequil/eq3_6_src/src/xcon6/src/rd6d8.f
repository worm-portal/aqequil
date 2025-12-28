      subroutine rd6d8(awmaxi,awmini,cbsri,cdac,cesri,cgexj,
     $ csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,
     $ dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,
     $ dlxplo,dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,
     $ iact,ibsrti,iesrti,ietmax,iktmax,imchmx,imech,iodb,iopg,
     $ iopr,iopt,igerti,itermx,ixrti,jcode,jgerti,jetmax,jflgi,
     $ jgext,jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,
     $ ksppmx,kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,
     $ mtbaqi,mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,
     $ nert,nertmx,net,netmax,nffg,nffgmx,ngexrt,ninpts,nobswt,
     $ nodbmx,nopgmx,noprmx,noptmx,nordmx,noutpt,nprob,nprpmx,
     $ nprpti,nprsmx,nprsti,nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,
     $ nsrt,nsrtmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,nttyo,nxmdmx,
     $ nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrt,nxrtmx,o2maxi,o2mini,
     $ phmaxi,phmini,pressb,pressi,ptk,qend,qgexsh,qrderr,rkb,
     $ rxbari,sfcar,ssfcar,tempcb,tempci,tgexp,timmxi,tistti,tolbt,
     $ toldl,tolsat,tolxsf,trkb,ttk,ubmtbi,ubsri,ucxri,udac,uesri,
     $ uffg,ugerji,ugermo,ugersi,ugexj,ugexmo,ugexp,ugexr,uhfgex,
     $ uobsw,uprphi,uprspi,ureac,usbsw,utitl1,utitl2,uvfgex,uxcat,
     $ uxkgex,uxmod,uxopex,uxopt,uzveci,vreac,xgersi,xhfgex,ximaxi,
     $ xistti,xlkffg,xlkgex,xlkmod,xvfgex,zgexj,zvclgi)
c
c     This subroutine reads the EQ6 input file in menu-style ("D")
c     format for version 8.0.
c
c     This subroutine is a near-clone of EQ6/rd6ind.f. However, the
c     present subroutine embodies only a pure read function (it does
c     only minimal checking of what is read to ensure that what
c     follows is readable). EQ6/rd6ind.f differs in that it also
c     writes an instant echo of what is read to the EQ6 output file.
c
c     The calling sequence of this subroutine is identical to that of
c     EQ6/rd6ind.f, EQ6/rd6inw.f, and XCON6/rd6w8.f.
c
c     This subroutine is called by:
c
c       XCON6/xcon6.f
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
      integer ietmax,iktmax,imchmx,jetmax,kmax,nbtmax,nbt1mx,nctmax,
     $ ndctmx,nertmx,netmax,nffgmx,nodbmx,nopgmx,noprmx,noptmx,nordmx,
     $ nprpmx,nprsmx,nptkmx,nrctmx,nsrtmx,ntitmx,nttkmx,nxmdmx,nxopmx,
     $ nxpemx,nxrtmx
c
      integer ninpts,noutpt,nttyo
c
      integer iact(imchmx,2,nrctmx),ibsrti(nsrtmx),
     $ igerti(jetmax,nertmx),iesrti(nsrtmx),imech(2,nrctmx),
     $ iodb(nodbmx),iopg(nopgmx),iopr(noprmx),iopt(noptmx),
     $ ixrti(nxrtmx),jcode(nrctmx),jgerti(nertmx),jflgi(nbtmax),
     $ jgext(netmax),jreac(nrctmx),kxmod(nxmdmx),ndact(imchmx,2,nrctmx),
     $ ngexrt(jetmax,netmax),nrk(2,nrctmx),nsk(nrctmx)
c
      integer itermx,jpress,jtemp,kbt,kct,kdim,kmt,kprs,ksplmx,ksppmx,
     $ kstpmx,kxt,nbti,nert,net,nffg,nobswt,nprob,nprpti,nprsti,nrct,
     $ nsbswt,nsrt,ntitl1,ntitl2,ntrymx,nxmod,nxopex,nxopt,nxrt
c
      logical qend,qgexsh,qrderr
c
      character*80 utitl1(ntitmx),utitl2(ntitmx)
      character*56 ugexr(ietmax,jetmax,netmax)
      character*48 ubmtbi(nbtmax),uobsw(2,nbtmax),uprspi(nprsmx),
     $ usbsw(2,nbtmax),uxmod(nxmdmx),uzveci(kmax)
      character*24 ubsri(nbt1mx,nsrtmx),ucxri(iktmax,nxrtmx),
     $ udac(ndctmx,imchmx,2,nrctmx),uffg(nffgmx),ugermo(nertmx),
     $ ugersi(ietmax,jetmax,nertmx),ugexmo(netmax),ugexp(netmax),
     $ uprphi(nprpmx),ureac(nrctmx),uxcat(nxopmx),uxopex(nxpemx)
      character*8 uesri(nctmax,nsrtmx),ugerji(jetmax,nertmx),
     $ ugexj(jetmax,netmax),uhfgex(ietmax,jetmax,netmax),
     $ uvfgex(ietmax,jetmax,netmax),uxkgex(ietmax,jetmax,netmax),
     $ uxopt(nxopmx)
c
      real*8 cbsri(nbt1mx,nsrtmx),cdac(ndctmx,imchmx,2,nrctmx),
     $ cesri(nctmax,nsrtmx),cgexj(jetmax,netmax),
     $ csigma(imchmx,2,nrctmx),eact(imchmx,2,nrctmx),
     $ egersi(ietmax,jetmax,nertmx),fkrc(nrctmx),hact(imchmx,2,nrctmx),
     $ modr(nrctmx),moffg(nffgmx),morr(nrctmx),mprphi(nprpmx),
     $ mprspi(nprsmx),mtbaqi(nbtmax),mtbi(nbtmax),
     $ mwtges(netmax),ptk(nptkmx),rkb(imchmx,2,nrctmx),
     $ rxbari(iktmax,nxrtmx),sfcar(nrctmx),ssfcar(nrctmx),
     $ tgexp(netmax),trkb(imchmx,2,nrctmx),ttk(nttkmx),vreac(nrctmx),
     $ xgersi(ietmax,jetmax,nertmx),xhfgex(ietmax,jetmax,netmax),
     $ xlkffg(nffgmx),xlkgex(ietmax,jetmax,netmax),xlkmod(nxmdmx),
     $ xvfgex(ietmax,jetmax,netmax),zvclgi(kmax),zgexj(jetmax,netmax)
c
      real*8 awmaxi,awmini,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,
     $ dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,
     $ dlxplo,dlxprl,dlxprn,ehmaxi,ehmini,electr,o2maxi,o2mini,phmaxi,
     $ phmini,pressb,pressi,tempcb,tempci,timmxi,tistti,tolbt,toldl,
     $ tolsat,tolxsf,ximaxi,xistti
c
c-----------------------------------------------------------------------
c
      include 'eqlib/eqlo8.h'
      include 'eqlib/eqlk8.h'
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
      integer i,icount,ii,iki,imh,ival,ivar,j,jcox,jd,je,jei,jj,jlast,
     $ jrex,jnrk,j2,j3,j4,krow,k1,k2,k3,k4,kxmd,n,nbi,ndt,ne,nei,ner,
     $ nfi,nfldmx,nfldt,nfldtx,nlchmx,nmark,nn,nnn,nobsw,npi,nsbsw,nsr,
     $ nsi,nxi,nxic,nrc,nxr,nci
c
      integer ilnobl
c
      logical qgexbf,qgexbs,qgexrd,qmark,qnone,qnonep,qnones,qnoner,
     $ qokay,qstop
c
      character*(nlchpa) ufield(nfldpa)
      character*(nlchpa) uline1,uline2,ulscr
      character*(nlchpa) uheadr,uheadx
c
      character*80 ustr
      character*48 ux48
      character*24 ustr24,ustrn
      character*8 ux8
      character*1 ux1
c
      real*8 var
c
c-----------------------------------------------------------------------
c
c     The following is nonsense to avoid compiler "unused variable"
c     warnings. Here noutpt and nprob are not actually used. They are
c     included in the calling sequence to allow it to match that of
c     EQ6/rd6inw.f.
c
      noutpt = nttyo
      i = nprob
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
 3000   format(/' * Error - (XCON6/rd6d8) The dimensioning parameter',
     $  ' for the',/7x,'number of iodb debugging print option switches',
     $  ' with string definitions',/7x,'(ndbxpa) has a value of ',
     $  i3,', but the dimensioning',/7x,'parameter for the number',
     $  ' of such switches (nodbpa) has a',/7x,'value of ',i3,'.')
        qstop = .true.
      endif
c
      if (npgxpa .ne. nopgmx) then
        write (nttyo,3010) npgxpa,nopgmx
 3010   format(/' * Error - (XCON6/rd6d8) The dimensioning parameter',
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
 3020   format(/' * Error - (XCON6/rd6d8) The dimensioning parameter',
     $  ' for the',/7x,'number of iopr print option switches',
     $  ' with string definitions',/7x,'(nprxpa) has a value of ',
     $  i3,', but the dimensioning',/7x,'parameter for the number',
     $  ' of such switches (noprpa) has a',/7x,'value of ',i3,'.')
        qstop = .true.
      endif
c
      if (nptxpa .ne. noptmx) then
        write (nttyo,3030) nptxpa,noptmx
 3030   format(/' * Error - (XCON6/rd6d8) The dimensioning parameter',
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
 1010   format(/' * Error - (XCON6/rd6d8) The first line of a "D"',
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
c     Read the block title ("Main Title") from a two-line header.
c
      uheadx = 'Main Title'
      nfldtx = 2
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c     Now read the main title itself.
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
 1015     format(/' * Error - (XCON6/rd6d8) Have too many lines in',
     $    /7x,'the main title. The code is only dimensioned for ',i4,
     $    /7x,'lines. Reduce the size of the main title or increase',
     $    /7x,'the dimensioning parameter ntitpa.')
          go to 990
        endif
c
        utitl1(n) = ufield(1)
      enddo
c
  120 ntitl1 = n
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Temperature.
c
      jtemp = -1
      icount = 0
      tempcb = 0.
      ttk(1) = 0.
      ttk(2) = 0.
c
c     Read a one-line header.
c
      uheadx = 'Temperature option (jtemp):'
      nfldtx = 1
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
c
c     Read the first option (constant temperature) from a one-line
c     header.
c
      nfldtx = 1
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      ustr24 = ustr(5:29)
      uheadx = '( 0) Constant temperatur'
      call locase(ustr24)
      call locase(uheadx)
      j2 = ilnobl(ustr24)
      j3 = ilnobl(uheadx)
      if (ustr24(1:j2) .ne. uheadx(1:j3)) then
        write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
 1020   format(/' * Error - (XCON6/rd6d8) Was expecting to find a',
     $  ' line beginning with',/7x,'"',a,'", instead found one',
     $  ' beginning with',/7x,'"',a,'".')
        qrderr = .true.
        go to 999
      endif
      if (ustr(1:1).ne.'[' .or. ustr(3:3).ne.']') then
        j2 = ilnobl(ustr)
        j2 = min(j2,70)
        write (nttyo,1022) ustr(1:j2)
 1022   format(/' * Error - (XCON6/rd6d8) Was expecting to find an',
     $  /7x,'option check box "[ ]" in the line beginning with',
     $  /7x,'"',a,'".')
        qrderr = .true.
        go to 999
      endif
      ux1 = ustr(2:2)
      if (ux1.eq.'*' .or. ux1.eq.'x' .or. ux1.eq.'X') then
        jtemp = 0
        icount = icount + 1
      endif
c
c     Read the associated temperature ("tempcb", C) from a one-line
c     header.
c
      uheadx = 'Value (C)'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      if (jtemp .eq. 0) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        tempcb = var
      endif
c
c     Read the second option (linear tracking in Xi) from a one-line
c     header.
c
      nfldtx = 1
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      ustr24 = ustr(5:29)
      uheadx = '( 1) Linear tracking in '
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
        jtemp = 1
        icount = icount + 1
      endif
c
c     Read the associated base temperature ("tempcb", C) from a one-line
c     header.
c
      uheadx = 'Base Value (C)'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      if (jtemp .eq. 1) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        tempcb = var
      endif
c
c     Read the associated derivative ("ttk(1)", dT/dXi, C/mol) from a
c     one-line header.
c
      uheadx = 'Derivative'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      if (jtemp .eq. 1) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        ttk(1) = var
      endif
c
c     Read the third option (linear tracking in time) from a one-line
c     header.
c
      nfldtx = 1
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      ustr24 = ustr(5:29)
      uheadx = '( 2) Linear tracking in '
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
        jtemp = 2
        icount = icount + 1
      endif
c
c     Read the associated base temperature ("tempcb", C) from a one-line
c     header.
c
      uheadx = 'Base Value (C)'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      if (jtemp .eq. 2) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        tempcb = var
      endif
c
c     Read the associated derivative ("ttk(1)", dT/dt, C/sec) from a
c     one-line header.
c
      uheadx = 'Derivative'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      if (jtemp .eq. 2) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        ttk(1) = var
      endif
c
c     Read the fourth option (fluid mixing tracking) option from a
c     one-line header.
c
      nfldtx = 1
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      ustr24 = ustr(5:29)
      uheadx = '( 3) Fluid mixing tracki'
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
        jtemp = 3
        icount = icount + 1
      endif
c
c     Read the temperature of fluid 1 ("tempcb", C) from a one-line
c     header.
c
      uheadx = 'T of fluid 1 (C)'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      if (jtemp .eq. 3) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        tempcb = var
      endif
c
c     Read the temperature of fluid 2 ("ttk(2)", C) from a one-line
c     header.
c
      uheadx = 'T of fluid 2 (C)'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      if (jtemp .eq. 3) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        ttk(2) = var
      endif
c
c     Read the mass ratio factor ("ttk(1)") from a two-line header.
c     This this the mass ratio of fluid2/fluid1 at Xi = 1.
c
      uheadx = 'Mass ratio factor'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      if (jtemp .eq. 3) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        ttk(1) = var
      endif
c
      if (icount.eq.0 .or. jtemp.lt.0) then
        write (nttyo,1023) tempcb
 1023   format(/' * Warning - (XCON6/rd6d8) No option was selected for',
     $  /7x,'the temperature. The temperature will be set at a fixed',
     $  ' value',/7x,'of ',f7.2,'C.')
        jtemp = 0
        ttk(1) = 0.
        ttk(2) = 0.
      elseif (icount .gt. 1) then
        write (nttyo,1024) tempcb
 1024   format(/' * Warning - (XCON6/rd6d8) Multiple options were',
     $  ' selected for',/7x,'the temperature. The temperature will be',
     $  ' set at a fixed value',/7x,'of ',f7.2,'C.')
        jtemp = 0
        ttk(1) = 0.
        ttk(2) = 0.
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Pressure.
c
      jpress = -1
      icount = 0
      pressb = 0.
      ptk(1) = 0.
      ptk(2) = 0.
c
c     Read a one-line header.
c
      uheadx = 'Pressure option (jpress):'
      nfldtx = 1
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
c
c     Read the first option (follow the data file reference pressure
c     curve) from a one-line header.
c
      nfldtx = 1
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      ustr24 = ustr(5:29)
      uheadx = '( 0) Follow the data fil'
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
        jpress = 0
        pressb = 0.
        icount = icount + 1
      endif
c
c     Read the second option (follow the 1.013-bar/steam-saturation
c     curve) from a one-line header.
c
      nfldtx = 1
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      ustr24 = ustr(5:29)
      uheadx = '( 1) Follow the 1.013-ba'
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
        jpress = 1
        pressb = 0.
        icount = icount + 1
      endif
c
c     Read the third option (specified constant pressure) from a
c     one-line header.
c
      nfldtx = 1
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      ustr24 = ustr(5:29)
      uheadx = '( 2) Constant pressure: '
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
        jpress = 2
        icount = icount + 1
      endif
c
c     Read the associated pressure ("pressb", bars) from a one-line
c     header.
c
      uheadx = 'Value (bars)'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      if (jpress .eq. 2) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        pressb = var
      endif
c
c     Read the fourth option (linear tracking in Xi) from a one-line
c     header.
c
      nfldtx = 1
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      ustr24 = ustr(5:29)
      uheadx = '( 3) Linear tracking in '
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
        jpress = 3
        icount = icount + 1
      endif
c
c     Read the base pressure ("pressb", bars) from a one-line header.
c
      uheadx = 'Base Value (bars)'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      if (jpress .eq. 3) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        pressb = var
      endif
c
c     Read the associated derivative ("ptk(1)", dP/dXi, bars/mol) from
c     a one-line header.
c
      uheadx = 'Derivative'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      if (jpress.eq. 3) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        ptk(1) = var
      endif
c
c     Read the fifth option (linear tracking in time) from a one-line
c     header.
c
      nfldtx = 1
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(1)
      ustr24 = ustr(5:29)
      uheadx = '( 4) Linear tracking in '
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
        jpress = 4
        icount = icount + 1
      endif
c
c     Read the base pressure ("pressb", bars) from a one-line header.
c
      uheadx = 'Base Value (bars)'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      if (jpress .eq. 4) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        pressb = var
      endif
c
c     Read the associated derivative ("ptk(1)", dP/dt, bars/sec) from
c     a two-line header.
c
      uheadx = 'Derivative'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      if (jpress .eq. 4) then
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        ptk(1) = var
      endif
c
      if (icount.eq.0 .or. jpress.lt.0) then
        write (nttyo,1030)
 1030   format(/' * Warning - (XCON6/rd6d8) No option was selected for',
     $  ' the pressure.',/7x,'The pressure will be set to be in',
     $  ' accord with the',/7x,'data file reference pressure curve.')
        jpress = 0
        pressb = 0.
        ptk(1) = 0.
        ptk(2) = 0.
      elseif (icount .gt. 1) then
        write (nttyo,1026)
 1026   format(/' * Warning - (XCON6/rd6d8) Multiple options were',
     $  ' selected for',/7x,'the pressure. The pressure will be set',
     $  ' to be in accord with the',/7x,'data file reference pressure',
     $  ' curve.')
        jpress = 0
        pressb = 0.
        ptk(1) = 0.
        ptk(2) = 0.
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Reactants.
c
      nrct = 0
      nsrt = 0
      nxrt = 0
      nert = 0
c
c     Read a two-line header for the block.
c
      uheadx = 'Reactants (Irreversible Reactions)'
      nfldtx = 2
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c     Loop on Reactants.
c
      nrc = 0
      nsr = 0
      nxr = 0
      ner = 0
      do nn = 1,nrctmx + 1
c
c       Read a line. If the block has not been completely read, this
c       contains the name of a reactant. Otherwise, this line is the
c       first line of the block following the reactants super-block.
c
        nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        ustr = ufield(1)
        uheadx = 'Reactant'
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
c
c         Done reading the reactants super-block.
c
          go to 240
        else
c
c         Back up.
c
          backspace ninpts
c
c         Re-read the reactant line with the test on the number
c         of fields.
c
          nfldtx = 3
          call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uline1,ulscr)
          if (qrderr) go to 999
        endif
c
        ustr = ufield(2)
        ustrn = ustr
        call locase(ustrn)
        if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') ustr = 'None'
        if (ustr(1:5) .eq. 'None ') then
          nrc = 0
        else
          nrc = nrc + 1
        endif
c
        if (nrc .gt. nrctmx) then
          write (nttyo,1520) nrctmx
 1520     format(/' * Error - (XCON6/rd6d8) Have too many reactants',
     $    /7x,'The code is only dimensioned for ',i4,' reactants.',
     $    /7x,'Reduce the number of reactants or increase the',
     $    /7x,'dimensioning parameter nrctpa.')
          go to 990
        endif
c
        if (nrc .gt. 0) ureac(nrc) = ustr
c
c       Read the separator line following the line containing the name
c       of a reactant.
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
c       Read the type of reactant from a two-line header.
c
        uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
        if (qrderr) go to 999
        ustr24 = ufield(2)
        uheadx = 'Type'
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
        call locase(ustr)
c
        if (nrc .eq. 0) go to 150
        do jcox = 0,5
          uheadx = urcjco(jcox)
          call locase(uheadx)
          if (ustr(1:16) .eq. uheadx(1:16)) then
            jcode(nrc) = jcox
            go to 150
          endif
        enddo
c
        j2 = ilnobl(ustr)
        write (nttyo,1060) ustr(1:j2)
 1060   format(/" * Error - (XCON6/rd6d8) Don't recognize the",
     $  ' jcode option string',/7x,' "',a,'". This should',
     $  ' be one of the strings',/7x,'defined in the urcjco array.',
     $  ' The valid strings are:',/)
        do jcox = 0,5
          j3 = ilnobl(urcjco(jcox))
          write (nttyo,1062) urcjco(jcox)(1:j3)
 1062     format(9x,a)
        enddo
        go to 990
c
c       Read the reactant status from a two-line header.
c
  150   uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
        if (qrderr) go to 999
        ustr24 = ufield(2)
        uheadx = 'Status'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)
        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
          write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
          qrderr = .true.
          go to 999
        endif
        if (nrc .eq. 0) go to 160
        ustr = ufield(3)
        call locase(ustr)
c
        do jrex = -1,3
          uheadx = urcjre(jrex)
          call locase(uheadx)
          if (ustr(1:16) .eq. uheadx(1:16)) then
            jreac(nrc) = jrex
            go to 160
          endif
        enddo
c
        j2 = ilnobl(ustr)
        write (nttyo,1065) ustr(1:j2)
 1065   format(/" * Error - (XCON6/rd6d8) Don't recognize the",
     $  ' jreac option string',/7x,' "',a,'". This should',
     $  ' be one of the strings',/7x,'defined in the urcjre array.',
     $  ' The valid strings are:',/)
        do jrex = -1,3
          j3 = ilnobl(urcjre(jrex))
          write (nttyo,1067) urcjre(jrex)(1:j3)
 1067     format(9x,a)
        enddo
        go to 990
c
c       Read the amount of reactant remaining (moles)
c       from a two-line header.
c
  160   uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
        if (qrderr) go to 999
        ustr24 = ufield(2)
        uheadx = 'Amount remaining (moles)'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)
        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
          write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
          qrderr = .true.
          go to 999
        endif
        if (nrc .gt. 0) then
          ustr = ufield(3)
          call chreal(nttyo,qrderr,ustr,var)
          if (qrderr) go to 999
          morr(nrc) = var
        endif
c
c       Read the amount of destroyed remaining (moles)
c       from a two-line header.
c
        uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
        if (qrderr) go to 999
        ustr24 = ufield(2)
        uheadx = 'Amount destroyed (moles)'
        call locase(ustr24)
        call locase(uheadx)
        j2 = ilnobl(ustr24)
        j3 = ilnobl(uheadx)
        if (ustr24(1:j2) .ne. uheadx(1:j3)) then
          write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
          qrderr = .true.
          go to 999
        endif
        if (nrc .gt. 0) then
          ustr = ufield(3)
          call chreal(nttyo,qrderr,ustr,var)
          if (qrderr) go to 999
          modr(nrc) = var
        endif
c
c       Skip jcode (reactant type) tests if no reactant is present.
c
        if (nrc .eq. 0) go to 235
c
        if (jcode(nrc) .eq. 1) then
c
c         Have a solid solution.
c
          nxr = nxr + 1
c
          if (nxr .gt. nxrtmx) then
            write (nttyo,1550) nxrtmx
 1550       format(/' * Error - (XCON6/rd6d8) Have too many solid',
     $      /7x,'solution reactants. The code is only dimensioned',
     $      /7x,'for ',i4,' such reactants. Reduce the number of',
     $      /7x,'such reactants or increase the dimensioning',
     $      ' parameter nxrtpa.')
            go to 990
          endif
c
c         Composition sub-block header.
c
          uheadx = '->'
          nfldtx = 2
          call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
          if (qrderr) go to 999
          ustr24 = ufield(2)
          uheadx = 'Composition'
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
c         Composition table header.
c
          uheadx = '--->'
          nfldtx = 4
          call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
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
c         Loop on end-member components.
c
          iki = 0
          do jj = 1,iktmax + 1
c
c           Read a line. If the sub-block (for the components) has not
c           been completely read, this contains the name of a component
c           and its mole fraction. Otherwise, this line is the separator
c           line before the next sub-block (surface area).
c
            nfldtx = 0
            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $      nttyo,qrderr,ufield,uline1,ulscr)
            if (qrderr) go to 999
            ustr = ufield(1)
            if (ustr(1:8) .eq. '--------') go to 210
c
c           Have found another component.
c
            iki = iki + 1
c
            if (iki .gt. iktmax) then
              j2 = ilnobl(ureac(nrc))
              write (nttyo,1590) ureac(nrc)(1:j2),iktmax
 1590         format(/' * Error - (XCON6/rd6d8) Have too many',
     $        ' end-members',/7x,'in the solid solution reactant',
     $        a,'.',/7x,'The code is only dimensioned for ',
     $        i4,' end-members per',/7x,'solid solution. Reduce',
     $        ' the number of end-members or',
     $        /7x,'increase the dimensioning parameter iktpar.')
              go to 990
            endif
c
            ustr = ufield(3)
            call chreal(nttyo,qrderr,ustr,var)
            if (qrderr) go to 999
            rxbari(iki,nxr) = var
            ucxri(iki,nxr) = ufield(2)
          enddo
  210     ixrti(nxr) = iki
        endif
c
        if (jcode(nrc) .eq. 2) then
c
c         Have a special reactant. Read its molar volume, composition,
c         and reaction.
c
          nsr = nsr + 1
c
          if (nsr .gt. nsrtmx) then
            write (nttyo,1600) nsrtmx
 1600       format(/' * Error - (XCON6/rd6d8) Have too many special',
     $      /7x,'reactants. The code is only dimensioned for ',i4,
     $      /7x,'such reactants. Reduce the number of such reactants',
     $      /7x,'or increase the dimensioning parameter nsrtpa.')
            go to 990
          endif
c
c         Read the molar volume (cm3/mol) from a two-line header.
c
          uheadx = '->'
          nfldtx = 4
          call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
          if (qrderr) go to 999
          ustr24 = ufield(2)
          uheadx = 'Molar volume (cm3/mol)'
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
          vreac(nrc) = var
c
c         Composition sub-block header.
c
          uheadx = '->'
          nfldtx = 2
          call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
          if (qrderr) go to 999
          ustr24 = ufield(2)
          uheadx = 'Composition'
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
c         Composition table header.
c
          uheadx = '--->'
          nfldtx = 4
          call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
          if (qrderr) go to 999
          ustr24 = ufield(2)
          uheadx = 'Element'
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
c         Loop on elements.
c
          nci = 0
          do jj = 1,nctmax + 1
c
c           Read a line. If the sub-block (for the elements) has not
c           been completely read, this contains the name of an element
c           and  the stoichiometric number. Otherwise, this line is the
c           separator line before the next sub-block (reaction).
c
            nfldtx = 0
            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $      nttyo,qrderr,ufield,uline1,ulscr)
            if (qrderr) go to 999
            ustr = ufield(1)
            if (ustr(1:8) .eq. '--------') go to 220
c
c           Have found another element.
c
            nci = nci + 1
c
            if (nci .gt. nctmax) then
              j2 = ilnobl(ureac(nrc))
              write (nttyo,1660) ureac(nrc)(1:j2),nctmax
 1660         format(/' * Error - (XCON6/rd6d8) Have too many',
     $        '  chemical',/7x,'elements in the special reactant ',a,
     $        '.',/7x,'The code is only dimensioned for ',i4,
     $        ' elements.',/7x,'Reduce the number of elements or',
     $        ' increase the',/7x,'dimensioning parameter nctpar.')
              go to 990
            endif
c
            ustr = ufield(3)
            call chreal(nttyo,qrderr,ustr,var)
            if (qrderr) go to 999
            cesri(nci,nsr) = var
            uesri(nci,nsr) = ufield(2)
          enddo
  220     iesrti(nsr) = nci
c
c         Reaction header.
c
          uheadx = '->'
          nfldtx = 2
          call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
          if (qrderr) go to 999
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
c         Species in the reaction.
c
          uheadx = '--->'
          nfldtx = 4
          call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
          if (qrderr) go to 999
          ustr24 = ufield(2)
          uheadx = 'Species'
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
c         Loop on species.
c
          nbi = 0
          do n = 1,nbt1mx+ 1
c
c           Read a line. If the sub-block (for the species) has not
c           been completely read, this contains the name of a species
c           and  the reaction coefficient. Otherwise, this line is the
c           separator line before the next sub-block (reaction).
c
            nfldtx = 0
            call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $      nttyo,qrderr,ufield,uline1,ulscr)
            if (qrderr) go to 999
            ustr = ufield(1)
            if (ustr(1:8) .eq. '--------') go to 230
c
c           Have found another species.
c
            nbi = nbi + 1
c
            if (nbi .gt. nbt1mx) then
              j2 = ilnobl(ureac(nrc))
              write (nttyo,1700) ureac(nrc)(1:j2),nbtmax
 1700         format(/' * Error - (XCON6/rd6d8) Have too many basis',
     $        ' basis species in the',/7x,'in the reaction for the',
     $        ' special reactant ',a,'.',/7x,'The code is only',
     $        ' dimensioned for ',i4,' basis species. Increase',
     $        /7x,'the dimensioning parameter nbtpar.')
              go to 990
            endif
c
            ustr = ufield(3)
            call chreal(nttyo,qrderr,ustr,var)
            if (qrderr) go to 999
            cbsri(n,nsr) = var
            ubsri(n,nsr) = ufield(2)
c
          enddo
  230     ibsrti(nsr) = nbi
        endif
c
        if (jcode(nrc) .eq. 5) then
c
c         Have a generic ion exchanger.
c
          ner = ner + 1
c
          if (ner .gt. nertmx) then
            write (nttyo,1710) nertmx
 1710       format(/' * Error - (XCON6/rd6d8) Have too many generic',
     $      ' ion exchanger',/7x,'reactants. The code is only',
     $      ' dimensioned for ',i4,' such reactants.',/7x,'Reduce',
     $      ' the number of such reactants or increase the',
     $      /7x,'dimensioning parameter nertpa.')
            go to 990
          endif
c
c         Read the exchange model name for the exchanger reactant
c         from a two-line header.
c
          uheadx = '->'
          nfldtx = 4
          call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
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
          ugermo(ner) = ufield(3)
c
c         Composition sub-block header.
c
          uheadx = '->'
          nfldtx = 2
          call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
          if (qrderr) go to 999
          ustr24 = ufield(2)
          uheadx = 'Composition'
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
c         Loop on exchange sites.
c
          jei = 0
          do jj = 1,jetmax + 1
c
c           Read a line. If the sub-block (for the composition of the
c           current exchanger reactant) has not been completely read,
c           this contains the name of an exchange site, and a sub-sub-
c           block for that site follows. Otherwise, this line is the
c           first line following the composition sub-block for the
c           current exchanger reactant.
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
c             Have found another exchange site sub-sub-block.
c
              jei = jei + 1
c
              if (jei .gt. jetmax) then
                write (nttyo,1716) ureac(nrc)(1:j2),jetmax
 1716           format(/' * Error - (XCON6/rd6d8) Have too many',
     $          ' exchange sites',/7x,'in the generic ion exchanger',
     $          ' reactant ',a,'.',/7x,'The code is only dimensioned',
     $          ' for ',i4,' exchange sites per',/7x,'generic ion',
     $          ' exchanger. Reduce the number of exchange sites or',
     $          /7x,'increase the dimensioning parameter jetpar.')
                go to 990
              endif
c
              ugerji(jei,ner) = ufield(3)
            else
c
c             Should have found the first line after the sub-block for
c             the composition of the current exchanger reactant.
c
c             Back up.
c
              backspace ninpts
c
              go to 234
            endif
c
c           Read the separator line following the line containing
c           the name of an exchange site.
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
c           Read the sub-sub-sub-block for the composition of the
c           current site. This consists of a table header, followed
c           by a table.
c
c           Composition table header.
c
            uheadx = '----->'
            nfldtx = 4
            call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $      nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
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
c           Loop on components for the current site.
c
            nei = 0
            do nnn = 1,netmax + 1
c
c             Read a line. If the sub-sub-sub-block (for the components
c             for the current site) has not been completely read, this
c             contains the name of a component and its mole fraction on
c             the current site. Otherwise, this line is the separator
c             line before the next sub-sub-block (for the next site
c             of the exchanger reactant) or the sub-block (surface
c             area) following the composition sub-block.
c
              nfldtx = 0
              call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uline1,ulscr)
              if (qrderr) go to 999
              ustr = ufield(1)
              if (ustr(1:8) .eq. '--------') go to 232
c
c             Have found another component.
c
              nei = nei + 1
c
              if (nei .gt. netmax) then
                j2 = ilnobl(ureac(nrc))
                j3 = ilnobl(ugerji(jei,ner))
                write (nttyo,1724) ureac(nrc)(1:j2),
     $          ugerji(jei,ner)(1:j3),netmax
 1724           format(/' * Error - (XCON6/rd6d8) Have too many',
     $          ' species on exchange',/7x,'site ',a,' of the generic',
     $          ' ion exchanger reactant',/7x,a,'. The code is only',
     $          ' dimensioned for species per',/7x,'exchange site.',
     $          ' Reduce the number of end-members or increase',
     $          /7x,'the dimensioning parameter netpar.')
                go to 990
              endif
c
              ustr = ufield(3)
              call chreal(nttyo,qrderr,ustr,var)
              if (qrderr) go to 999
              ugersi(nei,jei,ner) = ufield(2)
              egersi(nei,jei,ner) = var
              xgersi(nei,jei,ner) = var
            enddo
  232       igerti(jei,ner) = nei
          enddo
  234     jgerti(ner) = jei
        endif
c
c       Surface area option.
c
        icount = 0
        nsk(nrc) = 0.
        sfcar(nrc) = 0.
        ssfcar(nrc) = 0.
c
c       Read a one-line header.
c
        uheadx = '->'
        nfldtx = 2
        call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
        if (qrderr) go to 999
        ustr24 = ufield(2)
        uheadx = 'Surface area option (nsk'
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
c       Read the first option (constant surface area) from a one-line
c       header.
c
        nfldtx = 2
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        ustr = ufield(2)
        ustr24 = ustr(5:29)
        uheadx = '( 0) Constant surface ar'
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
          nsk(nrc) = 0
          icount = icount + 1
        endif
c
c       Read the associated constant surface area value (sfcar(n),
c       cm2) from a one-line header.
c
        uheadx = '->'
        nfldtx = 4
        call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
        if (qrderr) go to 999
        ustr24 = ufield(2)
        uheadx = 'Value (cm2)'
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
        if (nsk(nrc) .le. 0) sfcar(nrc) = var
c
c       Read the second option (constant specific surface area)
c       from a one-line header.
c
        nfldtx = 2
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        ustr = ufield(2)
        ustr24 = ustr(5:29)
        uheadx = '( 1) Constant specific s'
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
          nsk(nrc) = 1
          icount = icount + 1
        endif
c
c       Read the associated constant specific surface area value
c       (ssfcar(n), cm2/g) from a one-line header.
c
        uheadx = '->'
        nfldtx = 4
        call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
        if (qrderr) go to 999
        ustr24 = ufield(2)
        uheadx = 'Value (cm2/g)'
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
        if (nsk(nrc) .eq. 1) ssfcar(nrc) = var
c
c       Read the third option (n**2/3 growth law: current surface
c       area) from a two-line header.
c
        nfldtx = 2
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        ustr = ufield(2)
        ustr24 = ustr(5:29)
        uheadx = '( 2) n**2/3 growth law-'
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
          nsk(nrc) = 2
          icount = icount + 1
        endif
c
c       Read the associated current surface area value (sfcar(n),
c       cm2) from a two-line header.
c
        uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
        if (qrderr) go to 999
        ustr24 = ufield(2)
        uheadx = 'Value (cm2)'
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
        if (nsk(nrc) .eq. 2) sfcar(nrc) = var
c
        if (icount.eq.0) then
          j2 = ilnobl(ureac(nrc))
          write (nttyo,1033) ureac(nrc)(1:j2)
 1033     format(/' * Warning - (XCON6/rd6d8) No option was selected',
     $    ' for the',/7x,'surface area of the reactant ',a,'. The',
     $    /7x,'surface area and specific surface area will each be set',
     $    /7x,'to zero.')
          sfcar(nrc) = 0.0
          ssfcar(nrc) = 0.0
        elseif (icount .gt. 1) then
          j2 = ilnobl(ureac(nrc))
          write (nttyo,1035) ureac(nrc)(1:j2),sfcar(nrc)
 1035     format(/' * Warning - (XCON6/rd6d8) Multiple options were',
     $    ' selected for',/7x,'the surface area of the reactant ',a,
     $    '. The',/7x,'surface area will be set to a fixed value of ',
     $    1pe10.3,' cm3.')
          ssfcar(nrc) = 0.0
        endif
c
c       Read the ratio of active to total surface area (f) from a
c       two-line header.
c
        uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
        if (qrderr) go to 999
        ustr24 = ufield(2)
        uheadx = 'Surface area factor'
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
        fkrc(nrc)= var
c
c       Loop over forward and backward directions. The forward
c       direction corresponds to the disappearance (e.g., dissolution,
c       dissociation) of the associated species. The backward direction
c       corresponds to its formation (e.g., precipitation, association).
c
        do jd = 1,2
          if (jd .eq. 1) then
            ux8 = 'forward'
          else
            ux8 = 'backward'
          endif
c
c         Read the defining rate law string from a two-line header.
c
          uheadx = '->'
          nfldtx = 4
          call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
          if (qrderr) go to 999
          ustr24 = ufield(2)
          if (jd .eq. 1) then
            uheadx = 'Forward rate law'
          else
            uheadx = 'Backward rate law'
          endif
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
          call locase(ustr)
c
          do jnrk = -1,3
            uheadx = urcnrk(jnrk,jd)
            call locase(uheadx)
            if (ustr(1:16) .eq. uheadx(1:16)) then
              nrk(jd,nrc) = jnrk
              go to 170
            endif
          enddo
c
  165     j2 = ilnobl(ustr)
          write (nttyo,1075) ustr(1:j2)
 1075     format(/" * Error - (XCON6/rd6d8) Don't recognize the",
     $    ' rate law option string',/7x,' "',a,'". This should',
     $    ' be one of the strings',/7x,'defined in the urcnrk array.',
     $    ' The valid strings are:',/)
          do jnrk = -1,3
            j3 = ilnobl(urcnrk(jnrk,jd))
            write (nttyo,1077) urcnrk(jnrk,jd)(1:j3)
 1077       format(9x,a)
          enddo
          go to 990
c
  170     continue
c
          if (jd .eq. 1) then
c
c           Trap illegal option for the forward direction.
c
            if (nrk(jd,nrc) .eq. 0) go to 165
          endif
c
c         Trap options requiring no direct input (e.g., use of the
c         rate law for the opposite direction, precipitation according
c         to instantaneous partial equilibrium).
c
          if (nrk(jd,nrc) .le. 0) go to 195
c
c         Continue with options requiring direct input.
c
          if (nrk(jd,nrc) .eq. 1) then
c
c           Arbitrary kinetics (relative rates, indifferent to time).
c
            imech(jd,nrc) = 3
c
            if (imech(jd,nrc) .gt. imchmx) then
              j2 = ilnobl(ux8)
              j3 = ilnobl(ureac(nrc))
              write (nttyo,1960) ux8(1:j2),ureac(nrc)(1:j3),imchmx
 1960         format(/' * Error - (XCON6/rd6d8) Have too many rate',
     $        /7x,'constants or corresponding mechanisms in the ',a,
     $        /7x,'rate law for reactant ',a,'. The code is only',
     $        /7x,'dimensioned for ',i2,' rate constants per rate law.',
     $        /7x,'Reduce the number of rate constants or increase the',
     $        /7x,'dimensioning parameter imchpa.')
              go to 990
            endif
c
            do i = 1,3
              uheadx = '--->'
              nfldtx = 4
              call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
              if (qrderr) go to 999
              ustr24 = ufield(2)
              uheadx = urcrel(i,1)
              call locase(ustr24)
              call locase(uheadx)
              j2 = ilnobl(ustr24)
              j3 = ilnobl(uheadx)
c
              if (ustr24(1:j2) .ne. uheadx(1:j3)) then
                write (nttyo,1080) ustr24(1:j2),uheadx(1:j3)
 1080           format(/" * Error - (XCON6/rd6d8) Have an incorrect",
     $          ' relative rate law header',/7x,'"',a,'".',
     $          ' Was expecting the string "',a,'".')
                go to 990
              endif
c
              ustr = ufield(3)
              call chreal(nttyo,qrderr,ustr,var)
              if (qrderr) go to 999
              rkb(i,jd,nrc) = var
            enddo
c
          elseif (nrk(jd,nrc) .eq. 2) then
c
c           Transition state theory (TST) rate law. Up to imchmx
c           parallel mechanisms are allowed.
c
c           Loop on mechanisms.
c
            imh = 0
            do ii = 1,imchmx + 1
c
c             Read a line. If the sub-block for a mechanism has not
c             been completely read, this contains the name of a
c             mechanism, and a sub-sub-block follows. Otherwise, if
c             jd = 1, this line is the first line of the next sub-block
c             (backward rate law), or, if jd = 2, the first line of the
c             block for the next reactant, if any, else the first line
c             of the block following the reactants super-block.
c
              uheadx = '--->'
              nfldtx = 0
              call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uline1,ulscr)
              if (qrderr) go to 999
              ustr = ufield(2)
              uheadx = 'Mechanism'
              call locase(ustr)
              call locase(uheadx)
c
              if (ustr(1:9) .ne. uheadx(1:9)) then
c
c               Done reading mechanisms.
c
                if (imh .le. 0) then
                  j2 = ilnobl(ux8)
                  j3 = ilnobl(ureac(nrc))
                  write (nttyo,1958) ux8(1:j2),ureac(nrc)(1:j3)
 1958             format(/' * Error - (XCON6/rd6d8) Have no',
     $            ' "mechanisms" specified in the ',a,
     $            /7x,'rate law for reactant ',a,'.')
                  go to 990
                endif
c
c               Set the number of mechanisms for the last
c               direction and reactant processed.
c
                imech(jd,nrc) = imh
c
c               Back up.
c
                backspace ninpts
                if (jd .eq. 1) then
c
c                 Go read the backward rate law data.
c
                  go to 195
                else
c
c                 Go read the data for the next reactant, if any.
c
                  go to 235
                endif
              endif
c
c             Read the separator line following the mechanism header.
c
              nfldtx = 1
              call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uline1,ulscr)
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
              imh = imh + 1
c
              if (imh .gt. imchmx) then
                j2 = ilnobl(ux8)
                j3 = ilnobl(ureac(nrc))
                write (nttyo,1960) ux8(1:j2),ureac(nrc)(1:j3),imchmx
                go to 990
              endif
c
c             Read the sigma (stoichiometric correction) factor from
c             a two-line header. In some treatments, m appears instead
c             (sigma = 1/m). Here sigma is represented by the "csigma"
c             variable.
c
              uheadx = '----->'
              nfldtx = 4
              call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
              if (qrderr) go to 999
              ustr24 = ufield(2)
              if (jd .eq. 1) then
                uheadx = 'sigma(i,+,n)'
              elseif (jd .eq. 2) then
                uheadx = 'sigma(i,-,n)'
              endif
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
              csigma(imh,jd,nrc) = var
c
c             Read the rate constant (k, "rkb", mol/cm2/sec) at the
c             base or reference temperature ("trkb", see below) from a
c             two-line header.
c
              uheadx = '----->'
              nfldtx = 4
              call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
              if (qrderr) go to 999
              ustr24 = ufield(2)
              if (jd .eq. 1) then
                uheadx = 'k(i,+,n) (mol/cm2/sec)'
              elseif (jd .eq. 2) then
                uheadx = 'k(i,-,n) (mol/cm2/sec)'
              endif
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
              rkb(imh,jd,nrc) = var
c
c             Read the base (or reference) temperature ("trkb", C)
c             from a two-line header.
c
              uheadx = '----->'
              nfldtx = 4
              call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
              if (qrderr) go to 999
              ustr24 = ufield(2)
              uheadx = 'Ref. Temperature (C)'
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
              trkb(imh,jd,nrc) = var
c
c             Temperature dependence option.
c
              icount = 0
              iact(imh,jd,nrc) = 0
c
c             Read a one-line header.
c
              uheadx = '----->'
              nfldtx = 2
              call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uline1,ulscr)
              if (qrderr) go to 999
              ustr24 = ufield(2)
              uheadx = 'Temperature dependence o'
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
c             Read the first option (no temperature dependence) from
c             a one-line header.
c
              uheadx = '----->'
              nfldtx = 2
              call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uline1,ulscr)
              if (qrderr) go to 999
              ustr = ufield(2)
              ustr24 = ustr(5:29)
              uheadx = '( 0) No temperature depe'
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
                iact(imh,jd,nrc) = 0
                icount = icount + 1
              endif
c
c             Read the second option (constant activation energy)
c             from a one-line header.
c
              uheadx = '----->'
              nfldtx = 2
              call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uline1,ulscr)
              if (qrderr) go to 999
              ustr = ufield(2)
              ustr24 = ustr(5:29)
              uheadx = '( 1) Constant activation'
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
                iact(imh,jd,nrc) = 1
                icount = icount + 1
              endif
c
c             Read the associated activation energy from a one-line
c             header.
c
              uheadx = '----->'
              nfldtx = 4
              call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uline1,ulscr)
              if (qrderr) go to 999
              ustr24 = ufield(2)
              uheadx = 'Value (kcal/mol)'
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
              eact(imh,jd,nrc) = var
c
c             Read the third option (constant activation enthalpy)
c             from a one-line header.
c
              uheadx = '----->'
              nfldtx = 2
              call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uline1,ulscr)
              if (qrderr) go to 999
              ustr = ufield(2)
              ustr24 = ustr(5:29)
              uheadx = '( 2) Constant activation'
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
                iact(imh,jd,nrc) = 2
                icount = icount + 1
              endif
c
c             Read the associated activation enthalpy from a two-line
c             header.
c
              uheadx = '----->'
              nfldtx = 4
              call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
              if (qrderr) go to 999
              ustr24 = ufield(2)
              uheadx = 'Value (kcal/mol)'
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
              hact(imh,jd,nrc) = var
c
              if (icount .ne. 1) then
                j2 = ilnobl(ux8)
                j3 = ilnobl(ureac(nrc))
                if (icount .le. 0) then
                  write (nttyo,1118) ux8(1:j2),imh,ureac(nrc)(1:j3)
 1118             format(/' * Warning - (XCON6/rd6d8) No option was',
     $            ' selected for',/7x,'treating the temperature',
     $            ' dependence of the ',a,' rate constant',/7x,'for',
     $            ' term ',i2,' of the rate equation for the reactant',
     $            ' ',a,'.',/7x,'The temperature dependence will be',
     $            ' set to zero.')
                else
                  write (nttyo,1120) ux8(1:j2),imh,ureac(nrc)(1:j3)
 1120             format(/' * Warning - (XCON6/rd6d8) Multiple options',
     $            ' were selected for',/7x,'treating the temperature',
     $            ' dependence of the ',a,' rate constant',/7x,'for',
     $            ' term ',i2,' of the rate equation for the reactant',
     $            ' ',a,'.',/7x,'The temperature dependence will be',
     $            ' set to zero.')
                endif
                iact(imh,jd,nrc) = 0
                eact(imh,jd,nrc) = 0.
                hact(imh,jd,nrc) = 0.
              endif
c
c             Read the species and associated stoichiometric numbers,
c             if any, appearing in the kinetic activity product.
c
c             Read the kinetic activity product header from a
c             two-line header.
c
              uheadx = '----->'
              nfldtx = 2
              call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
              if (qrderr) go to 999
              ustr24 = ufield(2)
              uheadx = 'Kinetic activity product'
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
c             Read the first line of the table header for the
c             kinetic activity product.
c
              uheadx = '------->'
              nfldtx = 3
              call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uline1,ulscr)
              if (qrderr) go to 999
              ustr24 = ufield(2)
              uheadx = 'Species'
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
c             Read the second line of the table header for the
c             kinetic activity product.
c
              uheadx = '------->'
              nfldtx = 3
              call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
              if (qrderr) go to 999
              ustr24 = ufield(2)
              if (jd .eq. 1) then
                uheadx = '(udac(j,i,1,n))'
              else
                uheadx = '(udac(j,i,2,n))'
              endif
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
c             Loop on species appearing in the kinetic activity
c             product.
c
              ndt = 0
              do jj = 1,ndctmx + 1
c
c               Read a line. If the sub-sub-block (for the current
c               mechanism) has not been completely read, this contains
c               the name of a species. Otherwise, this line is the
c               first line of the next sub-block (backward rate law).
c
                uheadx = '------->'
                nfldtx = 0
                call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $          nttyo,qrderr,ufield,uline1,ulscr)
                if (qrderr) go to 999
                ustr24 = ufield(2)
                uheadx = 'None'
                call locase(ustr24)
                call locase(uheadx)
                j2 = ilnobl(ustr24)
                j3 = ilnobl(uheadx)
                if (ustr24(1:j2) .eq. uheadx(1:j3)) then
c
c                 Have no species appearing in the kinetic activity
c                 product.
c
                  udac(1,imh,jd,nrc) = ufield(2)
                  cdac(1,imh,jd,nrc) = 0.
                  uheadx = '--------'
c
c                 Read a separator line.
c
                  nfldtx = 1
                  call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $            nttyo,qrderr,ufield,uline1,ulscr)
                  if (qrderr) go to 999
                  go to 175
                endif
c
                ndt = ndt + 1
c
                if (ndt .gt. ndctmx) then
                  j2 = ilnobl(ux8)
                  j3 = ilnobl(ureac(nrc))
                  write (nttyo,1130) i,ux8(1:j2),ureac(nrc)(1:j3),
     $            ndctmx
 1130             format(/' * Error - (XCON6/rd6d8) Have too many',
     $            /7x,'species in the activity product in term ',i2,
     $            /7x,'of the ',a,' direction rate law for reactant',
     $            /7x,a,'. The code is only dimensioned for ',i3,
     $            /7x,'such species. Reduce the number of such species',
     $            /7x,'or increase the dimensioning parameter ndctpa.')
                  go to 990
                endif
c
                udac(ndt,imh,jd,nrc) = ufield(2)
                ustr = ufield(3)
                call chreal(nttyo,qrderr,ustr,var)
                if (qrderr) go to 999
                cdac(ndt,imh,jd,nrc) = var
c
c               Read a line. If the sub-sub-sub-block (for the current
c               activity product) has been completely read, this line
c               is a separator line. Otherwise, it contains data for
c               another species in the current kinetic activity product.
c
                nfldtx = 0
                call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $          nttyo,qrderr,ufield,uline1,ulscr)
                if (qrderr) go to 999
                ustr = ufield(1)
                uheadx = '--------'
                if (ustr(1:8) .eq. uheadx(1:8)) go to 175
c
              enddo
c
  175         ndact(imh,jd,nrc) = ndt
c
            enddo
c
          elseif (nrk(jd,nrc) .eq. 3) then
c
c           Linear rate law.
c
c           Loop on mechanisms. However, only one "mechanism" is
c           currently allowed.
c
            imh = 0
            do ii = 1,imchmx + 1
c
c             Read a line. If the sub-block for a mechanism has not
c             been completely read, this contains the name of a
c             mechanism, and a sub-sub-block follows. Otherwise, if
c             jd = 1, this line is the first line of the next sub-block
c             (backward rate law), or, if jd = 2, the first line of the
c             block for the next reactant, if any, else the first line
c             of the block following the reactants super-block.
c
              uheadx = '--->'
              nfldtx = 0
              call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uline1,ulscr)
              if (qrderr) go to 999
              ustr = ufield(2)
              uheadx = 'Mechanism'
              call locase(ustr)
              call locase(uheadx)
c
              if (ustr(1:9) .ne. uheadx(1:9)) then
c
c               Done reading mechanisms.
c
                if (imh .le. 0) then
                  j2 = ilnobl(ux8)
                  j3 = ilnobl(ureac(nrc))
                  write (nttyo,1958) ux8(1:j2),ureac(nrc)(1:j3)
                  go to 990
                endif
c
c               Set the number of mechanisms for the last
c               direction and reactant processed.
c
                imech(jd,nrc) = imh
c
c               Back up.
c
                backspace ninpts
                if (jd .eq. 1) then
c
c                 Go read the backward rate law data.
c
                  go to 195
                else
c
c                 Go read the data for the next reactant, if any.
c
                  go to 235
                endif
              endif
c
c             Read the separator line following the mechanism header.
c
              nfldtx = 1
              call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uline1,ulscr)
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
              imh = imh + 1
c
              if (imh .gt. imchmx) then
                j2 = ilnobl(ux8)
                j3 = ilnobl(ureac(nrc))
                write (nttyo,1960) ux8(1:j2),ureac(nrc)(1:j3),imchmx
                go to 990
              endif
c
              if (imh .gt. 1) then
                j2 = ilnobl(ux8)
                j3 = ilnobl(ureac(nrc))
                write (nttyo,1964) ux8(1:j2),ureac(nrc)(1:j3)
 1964           format(/' * Error - (XCON6/rd6d8) Have more than the',
     $          ' allowed one "mechanism"',/7x,'in the ',a,' rate',
     $          'law for reactant ',a,'.')
                go to 990
              endif
c
c             Read the rate constant (k, "rkb", mol/cm2/sec) at the
c             base or reference temperature ("trkb", see below) from a
c             two-line header.
c
              uheadx = '----->'
              nfldtx = 4
              call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
              if (qrderr) go to 999
              ustr24 = ufield(2)
              if (jd .eq. 1) then
                uheadx = 'k(i,+,n) (mol/cm2/sec)'
              elseif (jd .eq. 2) then
                uheadx = 'k(i,-,n) (mol/cm2/sec)'
              endif
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
              rkb(imh,jd,nrc) = var
c
c             Read the base (or reference) temperature ("trkb", C)
c             from a two-line header.
c
              uheadx = '----->'
              nfldtx = 4
              call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
              if (qrderr) go to 999
              ustr24 = ufield(2)
              uheadx = 'Ref. Temperature (C)'
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
              trkb(imh,jd,nrc) = var
c
c             Temperature dependence option.
c
              icount = 0
              iact(imh,jd,nrc) = 0
c
c             Read a one-line header.
c
              uheadx = '----->'
              nfldtx = 2
              call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uline1,ulscr)
              if (qrderr) go to 999
              ustr24 = ufield(2)
              uheadx = 'Temperature dependence o'
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
c             Read the first option (no temperature dependence) from
c             a one-line header.
c
              uheadx = '----->'
              nfldtx = 2
              call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uline1,ulscr)
              if (qrderr) go to 999
              ustr = ufield(2)
              ustr24 = ustr(5:29)
              uheadx = '( 0) No temperature depe'
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
                iact(imh,jd,nrc) = 0
                icount = icount + 1
              endif
c
c             Read the second option (constant activation energy)
c             from a one-line header.
c
              uheadx = '----->'
              nfldtx = 2
              call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uline1,ulscr)
              if (qrderr) go to 999
              ustr = ufield(2)
              ustr24 = ustr(5:29)
              uheadx = '( 1) Constant activation'
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
                iact(imh,jd,nrc) = 1
                icount = icount + 1
              endif
c
c             Read the associated activation energy from a one-line
c             header.
c
              uheadx = '----->'
              nfldtx = 4
              call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uline1,ulscr)
              if (qrderr) go to 999
              ustr24 = ufield(2)
              uheadx = 'Value (kcal/mol)'
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
              eact(imh,jd,nrc) = var
c
c             Read the third option (constant activation enthalpy)
c             from a one-line header.
c
              uheadx = '----->'
              nfldtx = 2
              call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uline1,ulscr)
              if (qrderr) go to 999
              ustr = ufield(2)
              ustr24 = ustr(5:29)
              uheadx = '( 2) Constant activation'
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
                iact(imh,jd,nrc) = 2
                icount = icount + 1
              endif
c
c             Read the associated activation enthalpy from a two-line
c             header.
c
              uheadx = '----->'
              nfldtx = 4
              call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $        nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
              if (qrderr) go to 999
              ustr24 = ufield(2)
              uheadx = 'Value (kcal/mol)'
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
              hact(imh,jd,nrc) = var
c
              if (icount .ne. 1) then
                j3 = ilnobl(ux8)
                j2 = ilnobl(ureac(nrc))
                if (icount .le. 0) then
                  write (nttyo,1118) ux8(1:j3),imh,ureac(nrc)(1:j2)
                else
                  write (nttyo,1120) ux8(1:j3),imh,ureac(nrc)(1:j2)
                endif
                iact(imh,jd,nrc) = 0
                eact(imh,jd,nrc) = 0.
                hact(imh,jd,nrc) = 0.
              endif
            enddo
c
          else
c
c           Error, unknown rate law code.
c
            j3 = ilnobl(ux8)
            j2 = ilnobl(ureac(nrc))
            write (nttyo,1952) ux8(1:j3),nrk(jd,nrc),ureac(nrc)(1:j2)
 1952       format(/' * Error - (XCON6/rd6d8) Programming error trap:',
     $      ' The ',a,/7x,'rate law code has an unrecognized value',
     $      ' of ',i2,' for reactant',/7x,a,'.')
            go to 990
          endif
c
  195     continue
        enddo
c
c       End of the loop on the forward and backward directions of the
c       irreversible reaction corresponding to a reactant.
c
  235   continue
      enddo
c
c     End of the loop on reactants.
c
  240 nrct = nrc
      nsrt = nsr
      nxrt = nxr
      nert = ner
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Starting, minimum, and maximum values of key run parameters.
c     Read superblock header.
c
c     Read the data from a two-line header.
c
      uheadx =
     $ 'Starting, minimum, and maximum values of key run parameters.'
      nfldtx = 1
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Starting value of Xi.
c
      xistti = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Starting Xi value'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      xistti = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum value of Xi.
c
      ximaxi = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Maximum Xi value'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      ximaxi = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Starting value of time.
c
      tistti = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Starting time (seconds)'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      tistti = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum value of time.
c
      timmxi = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Maximum time (seconds)'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      timmxi = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Minimum value of pH.
c
      phmini = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Minimum value of pH'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      phmini = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum value of pH.
c
      phmaxi = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Maximum value of pH'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      phmaxi = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Minimum value of Eh (v).
c
      ehmini = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Minimum value of Eh (v)'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      ehmini = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum value of Eh (v).
c
      ehmaxi = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Maximum value of Eh (v)'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      ehmaxi = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Minimum value of log fO2.
c
      o2mini = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Minimum value of log fO2'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      o2mini = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum value of log fO2.
c
      o2maxi = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Maximum value of log fO2'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      o2maxi = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Minimum value of aw.
c
      awmini = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Minimum value of aw'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      awmini = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum value of aw.
c
      awmaxi = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Maximum value of aw'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      awmaxi = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum number of steps.
c
      kstpmx = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Maximum number of steps'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 999
      kstpmx = ivar
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print interval parameters. Read superblock header.
c
c     Read the data from a two-line header.
c
      uheadx = 'Print interval parameters.'
      nfldtx = 1
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Xi print interval.
c
      dlxprn = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Xi print interval'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      dlxprn = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Log Xi print interval.
c
      dlxprl = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Log Xi print interval'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      dlxprl = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Time print interval.
c
      dltprn = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Time print interval'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      dltprn = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Log time print interval.
c
      dltprl = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Log time print interval'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      dltprl = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     pH print interval.
c
      dlhprn = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'pH print interval'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      dlhprn = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Eh (v) print interval.
c
      dleprn = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Eh (v) print interval'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      dleprn = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Log fO2 print interval.
c
      dloprn = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Log fO2 print interval'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      dloprn = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Activity of water print interval.
c
      dlaprn = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'aw print interval'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      dlaprn = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Steps print interval.
c
      ksppmx = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Steps print interval'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 999
      ksppmx = ivar
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Plot interval parameters. Read superblock header.
c
c     Read the data from a two-line header.
c
      uheadx = 'Plot interval parameters.'
      nfldtx = 1
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Xi plot interval.
c
      dlxplo = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Xi plot interval'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      dlxplo = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Log Xi plot interval.
c
      dlxpll = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Log Xi plot interval'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      dlxpll = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Time plot interval.
c
      dltplo = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Time plot interval'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      dltplo = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Log time plot interval.
c
      dltpll = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Log time plot interval'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      dltpll = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     pH plot interval.
c
      dlhplo = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'pH plot interval'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      dlhplo = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Eh (v) plot interval.
c
      dleplo = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Eh (v) plot interval'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      dleplo = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Log fO2 plot interval.
c
      dloplo = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Log fO2 plot interval'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      dloplo = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Activity of water plot interval.
c
      dlaplo = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'aw plot interval'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      dlaplo = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Steps plot interval.
c
      ksplmx = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Steps plot interval'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 999
      ksplmx = ivar
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
          write (nttyo,1607) uheadx(1:j2)
 1607     format(/' * Error - (XCON6/rd6d8) The iopt option switch',
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
 1610     format(/' * Error - (XCON6/rd6d8) The iopt option switch',
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
 1612     format(/' * Error - (XCON6/rd6d8) The iopt option switch',
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
 1630       format(/' * Error - (XCON6/rd6d8) Have too many option',
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
 1640       format(/' * Error - (XCON6/rd6d8) The following iopt',
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
 1650     format(/' * Error - (XCON6/rd6d8) The iopt option switch',
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
          if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
            j4 = ilnobl(uopttx(n))
            write (nttyo,1662) uheadx(1:j2),ustr(1:j3),uopttx(n)(1:j4)
 1662       format(/' * Error - (XCON6/rd6d8) The iopt option switch',
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
 1670         format(/" * Error - (XCON6/rd6d8) Don't recognize the",
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
 1680     format(/' * Warning - (XCON6/rd6d8) No option choice was',
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
 1694     format(/' * Warning - (XCON6/rd6d8) More than one option',
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
c     Iopr Print Option Switches.
c     Note: iopr(1) = iopr1, etc.
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
 1800     format(/' * Error - (XCON6/rd6d8) The iopr option switch',
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
 1810     format(/' * Error - (XCON6/rd6d8) The iopr option switch',
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
 1812     format(/' * Error - (XCON6/rd6d8) The iopr option switch',
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
 1830       format(/' * Error - (XCON6/rd6d8) Have too many option',
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
 1840       format(/' * Error - (XCON6/rd6d8) The following iopr',
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
 1850     format(/' * Error - (XCON6/rd6d8) The iopr option switch',
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
          if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
            j4 = ilnobl(uoprtx(n))
            write (nttyo,1860) uheadx(1:j2),ustr(1:j3),uoprtx(n)(1:j4)
 1860       format(/' * Error - (XCON6/rd6d8) The iopr option switch',
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
 1870         format(/" * Error - (XCON6/rd6d8) Don't recognize the",
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
 1880     format(/' * Warning - (XCON6/rd6d8) No option choice was',
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
 1894     format(/' * Warning - (XCON6/rd6d8) More than one option',
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
        nfldtx = 0
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
 1900     format(/' * Error - (XCON6/rd6d8) The iodb option switch',
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
 1910     format(/' * Error - (XCON6/rd6d8) The iodb option switch',
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
 1920     format(/' * Error - (XCON6/rd6d8) The iodb option switch',
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
 1922       format(/' * Error - (XCON6/rd6d8) Have too many option',
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
 1940       format(/' * Error - (XCON6/rd6d8) The following iodb',
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
 1950     format(/' * Error - (XCON6/rd6d8) The iodb option switch',
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
          if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
            j4 = ilnobl(uodbtx(n))
            write (nttyo,1962) uheadx(1:j2),ustr(1:j3),uodbtx(n)(1:j4)
 1962       format(/' * Error - (XCON6/rd6d8) The iodb option switch',
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
 1970         format(/" * Error - (XCON6/rd6d8) Don't recognize the",
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
 1980     format(/' * Warning - (XCON6/rd6d8) No option choice was',
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
 1994     format(/' * Warning - (XCON6/rd6d8) More than one option',
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
c     Mineral sub-set selection suppression options.
c
c     Read the block title from a two-line header.
c
      uheadx = 'Mineral Sub-Set Selection Suppression Options'
      nfldtx = 2
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c     Read the option title from a two-line header.
c
      uheadx = 'Option'
      nfldtx = 0
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr24 = ufield(1)
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
      nxi = 0
c
c     Read the first line.
c
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
c
c     The label below marks a return point for processing subsequent
c     lines in the current data block.
c
  330 ustr = ufield(1)
c
c     There are no data remaining in the current block if a
c     separator line has been encountered.
c
      if (ustr(1:8) .eq. '--------') go to 380
c
      ustrn = ustr
      call locase(ustrn)
      if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') ustr = 'None'
      qnone = ustr(1:5).eq.'None '
c
      if (.not.qnone)  then
        nxi = nxi + 1
c
        if (nxi .gt. nxopmx) then
          write (nttyo,1320) nxopmx
 1320     format(/' * Error - (XCON6/rd6d8) Have too many mineral',
     $    /7x,'subset-selection suppression options. The code is',
     $    /7x,'only dimensioned for ',i3,' such options. Reduce the',
     $    /7x,'number of options or increase the dimensioning',
     $    /7x,'parameter nxoppa.')
          go to 990
        endif
c
        ustrn = ustr
        call locase(ustrn)
c
        do n = 1,4
          uheadx = uxopti(n)
          call locase(uheadx)
          if (ustrn(1:16) .eq. uheadx(1:16)) then
            uxopt(nxi) = uxopti(n)
            go to 350
          endif
        enddo
c
        j2 = ilnobl(ustr)
        write (nttyo,3060) ustr(1:j2)
 3060   format(/" * Error - (XCON6/rd6d8) Don't recognize the",
     $  ' uxopt option string',/7x,'"',a,'". This should',
     $  ' be one of the strings',/7x,'defined in the uxopti array.',
     $  ' The valid strings are:',/)
        do n = 1,4
          j3 = ilnobl(uxopti(n))
          write (nttyo,3062) uxopti(n)(1:j3)
 3062     format(9x,a)
        enddo
        go to 990
c
  350   uxcat(nxi) = ufield(2)
      endif
c
c     Read the next line. Go back to process it.
c
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      go to 330
c
  380 nxopt = nxi
c
c     Exceptions to the mineral sub-set selection suppression options.
c
c     Read the block title from a two-line header.
c
      uheadx = 'Exceptions to the Mineral Sub-Set Selection Suppression'
      nfldtx = 2
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c     Read the mineral title from a two-line header.
c
      uheadx = 'Mineral'
      nfldtx = 0
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr24 = ufield(1)
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
      nxic = 0
c
c     Read the first line.
c
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
c
c     The label below marks a return point for processing subsequent
c     lines in the current data block.
c
  410 ustr = ufield(1)
c
c     There are no data remaining in the current block if a
c     separator line has been encountered.
c
      if (ustr(1:8) .eq. '--------') go to 420
c
      ustrn = ustr
      call locase(ustrn)
      if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') ustr = 'None'
      qnone = ustr(1:5).eq.'None '
c
      if (.not.qnone)  then
        nxic = nxic + 1
c
        if (nxic .gt. nxpemx) then
          write (nttyo,1370) nxpemx
 1370     format(/' * Error - (XCON6/rd6d8) Have too many',
     $    /7x,'exceptions specified to the mineral subset-selection',
     $    /7x,'suppression options. The code is only dimensioned',
     $    /7x,'for ',i3,'exceptions. Reduce the number of exceptions',
     $    /7x,'or increase the dimensioning parameter nxpepa.')
          go to 990
        endif
c
        uxopex(nxic) = ufield(1)
      endif
c
c     Read the next line. Go back to process it.
c
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      go to 410
c
  420 nxopex = nxic
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Fixed fugacity options.
c
c     Read the block title from a two-line header.
c
      uheadx = 'Fixed Fugacity Options'
      nfldtx = 2
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c     Read the second part of the block title from two lines.
c
      uheadx = 'Gas'
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      ustr24 = ufield(1)
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
      uheadx = '(uffg(n))'
      nfldtx = 0
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr24 = ufield(1)
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
      nfi = 0
c
c     Read the first line.
c
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
c
c     The label below marks a return point for processing subsequent
c     lines in the current data block.
c
  430 ustr = ufield(1)
c
c     There are no data remaining in the current block if a
c     separator line has been encountered.
c
      if (ustr(1:8) .eq. '--------') go to 440
c
      ustrn = ustr
      call locase(ustrn)
      if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') ustr = 'None'
      qnone = ustr(1:5).eq.'None '
c
      if (.not.qnone) then
        nfi = nfi + 1
c
        if (nfi .gt. nffgmx) then
          write (nttyo,1420) nffgmx
 1420     format(/' * Error - (XCON6/rd6d8) Have too many gases whose',
     $    ' fugacities are to be fixed.',/7x,'The code is only',
     $    ' dimensioned for ',i4,' such gases. Reduce the number',
     $    /7x,'of gases or increase the dimensioning parameter nffgpa.')
          go to 990
        endif
c
        uffg(nfi) = ufield(1)
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        moffg(nfi) = var
        ustr = ufield(3)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        xlkffg(nfi) = var
      endif
c
c     Read the next line. Go back to process it.
c
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      go to 430
c
  440 nffg = nfi
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
c     Read the maximum finite-difference order (nordmx)
c     from a one-line header.
c
      uheadx = 'Max. finite-difference order'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 999
      nordmx = ivar
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
c     from a one-line header.
c
      uheadx = 'Max. No. of N-R iterations'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 999
      itermx = ivar
c
c     Read the search/find convergence tolerance (tolxsf) from
c     a one-line header.
c
      uheadx = 'Search/find convergence tolerance'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      tolxsf = var
c
c     Read the saturation tolerance (tolsat) from a one-line header.
c
      uheadx = 'Saturation tolerance'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      tolsat = var
c
c     Read the maximum number of phase assemblage tries (ntrymx) from
c     a one-line header.
c
      uheadx = 'Max. No. of Phase Assemblage Tries'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 999
      ntrymx = ivar
c
c     Read the zero order step size (in Xi) (dlxmx0) from
c     a one-line header.
c
      uheadx = 'Zero order step size (in Xi)'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      dlxmx0 = var
c
c     Read the maximum interval in Xi between PRS transfers
c     from a two-line header (the separator line is the end of the
c     current block).
c
      uheadx = 'Max. interval in Xi between PRS transfers'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      dlxdmp = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Secondary title.
c
c     Read the block title ("Secondary Title") from a two-line header.
c
      uheadx = 'Secondary Title'
      nfldtx = 2
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c     Now read the secondary title itself.
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
        if (ustr(1:8) .eq. '--------') go to 520
c
        n = n + 1
c
        if (n .gt. ntitmx) then
          write (nttyo,5015) ntitmx
 5015     format(/' * Error - (XCON6/rd6d8) Have too many lines in',
     $    /7x,'the secondary title. The code is only dimensioned for ',
     $    i4,/7x,'lines. Reduce the size of the secondary title or ',
     $    'increase',/7x,'the dimensioning parameter ntitpa.')
          go to 990
        endif
c
        utitl2(n) = ufield(1)
      enddo
c
  520 ntitl2 = n
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
  525 nsbswt = nsbsw
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Original temperature.
c
      tempci = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Original temperature (C)'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      tempci = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Original pressure.
c
      pressi = 0.
c
c     Read the data from a two-line header.
c
      uheadx = 'Original pressure (bars)'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      pressi = var
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Ion exchanger creation.
c
      net = 0
c
c     Read a two-line header for the block.
c
      uheadx = 'Create Ion Exchangers'
      nfldtx = 2
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c     Read the advisory from one line: one or more exchanger creation
c     blocks follow or not. This advisory will be confirmed by
c     examining the input that actually follows. Here qgexbf is
c     the advisory flag. It is .true. if the advisory indicates
c     that one or more exchanger creation blocks follow. Here
c     also qgexbs is .true. if qgexbf has been set.
c
      qgexbs = .false.
      nfldtx = 1
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
c
      ustr = ufield(1)
      call lejust(ustr)
      uheadx = 'Advisory:'
      call locase(ustr)
      call locase(uheadx)
      if (ustr(1:9) .ne. uheadx(1:9)) then
        qrderr = .true.
        j3 = ilnobl(uline1(1:70))
        write (nttyo,1210) uheadx(1:9),uline1(1:j3)
 1210   format(/' * Error - (XCON6/rd6d8) Was expecting to find the',
     $  /7x,'string "',a,'" at the start of the line beginning with',
     $  /7x,'"',a,'".')
        go to 999
      endif
c
      uheadx = 'Advisory: at least one exchanger creation block'
     $ // ' follows on this file.'
      call locase(ustr)
      call locase(uheadx)
      j2 = ilnobl(ustr)
      j3 = ilnobl(uheadx)
      if (ustr(1:j2) .eq. uheadx(1:j3)) then
        qgexbf = .true.
        qgexbs = .true.
      else
        uheadx =
     $  'Advisory: no exchanger creation blocks follow on this file.'
        call locase(uheadx)
        j3 = ilnobl(uheadx)
        if (ustr(1:j2) .eq. uheadx(1:j3)) then
          qgexbf = .false.
          qgexbs = .true.
        else
          qgexbf = .false.
          qgexbs = .false.
          j3 = ilnobl(uline1(1:70))
          write (nttyo,1220) uline1(1:j3)
 1220     format(/' * Warning - (XCON6/rd6d8) Could not interpret the',
     $    ' advisory on whether',/7x,'or not one or more exchanger',
     $    ' blocks follow. The problem is with',/7x,'the line',
     $    ' beginning with',/7x,'"',a,'".',/7x,'The presence of such',
     $    ' blocks will be determined directly.')
        endif
      endif
c
c     Read the first line of the qgexsh option (on processing, show
c     at least one exchanger creation block on a "D" format input or
c     pickup file).
c
      uheadx = 'Option: on further processing (writing a pickup file'
     $  // ' or running XCON6 on the'
      nfldtx = 1
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
c
c     Read the second line of the qgexsh option.
c
      uheadx = 'present file), force the inclusion of at least one'
     $ // ' such block (qgexsh):'
      nfldtx = 1
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
c
c     Read the third and last line of the qgexsh option, plus the
c     following separator line.
c
      nfldtx = 1
      call rdd2l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
      k1 = index(uline1,'[')
      k2 = index(uline1,']')
      j3 = ilnobl(uline1)
c
      if (k1.le.0 .or. k2.le.0 .or. k2.lt.(k1 + 2)) then
        write (nttyo,1230) uline1(1:j3)
 1230   format(/' * Error - (XCON6/rd6d8) Was expecting to find an',
     $  /7x,'option check box "[ ]" in the line beginning with',
     $  /7x,'"',a,'"',/7x,'The expected line contains the box for',
     $  ' the qgexsh flag.')
        go to 990
      endif
c
c     Check the full option choice string.
c
      ustr = '(.true.)'
      if (index(uline1(1:j3),ustr(1:8)) .le. 0) then
        write (nttyo,1240) uline1(1:j2),ustr(1:8)
 1240   format(/' * Error - (XCON6/rd6d8) The qgexsh flag line',
     $  /7x,'"',a,'"',/7x,"read from the input file doesn't contain",
     $  ' the expected string',/7x,'"',a,'".')
        go to 990
      endif
c
      ustr24 = uline1(k1 + 1:k2 -1)
      call lejust(ustr24)
      j2 = ilnobl(ustr24)
      qmark = .false.
      if (j2 .gt. 0) then
        if (index(ustr24,'*') .ge. 1) then
          qmark = .true.
        elseif (index(ustr24,'x') .ge. 1) then
          qmark = .true.
        elseif (index(ustr24,'X') .ge. 1) then
          qmark = .true.
        else
          write (nttyo,1250) ustr24(1:j2),uline1(1:j3)
 1250     format(/" * Error - (XCON6/rd6d8) Don't recognize the",
     $    ' string "',a,'"',/7x,'that appears on the qgexsh',
     $    ' option switch choice line',/7x,'"',a,'"',
     $    /7x,'read from the input file. An option choice should',
     $    ' be chosen by',/7x,'placing a "*", "x", or "X" in the',
     $    ' checkbox ("[ ]").')
          go to 990
        endif
      endif
c
      qgexsh = qmark
c
c     Check to see if an exchanger block actually follows. This is a
c     test of the qgexbf flag. Read a line. If it contains the string
c     'Exchanger phase' in the first field, an exchanger block is
c     present.
c
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
c
      ustr = ufield(1)
      uheadx = 'Exchanger phase'
      call locase(ustr)
      call locase(uheadx)
      j2 = ilnobl(ustr)
      j3 = ilnobl(uheadx)
c
      if (ustr(1:j2) .eq. uheadx(1:j3)) then
c
c       An exchanger block is present.
c
        qgexrd = .true.
        if (qgexbs .and. .not.qgexbf) then
          write (nttyo,1260)
 1260     format(/' * Note - (XCON6/rd6d8) The advisory on the',
     $    ' presence of',/7x,'ion exchanger blocks was negative,',
     $    ' but one or more blocks is',/7x,'actually present.')
        endif
      else
c
c       An exchanger block is not present.
c
        qgexrd = .false.
        if (qgexbs .and. qgexbf) then
          write (nttyo,1270)
 1270     format(/' * Note - (XCON6/rd6d8) The advisory on the',
     $    ' presence of',/7x,'ion exchanger blocks was positive,',
     $    ' but no blocks are actually',/7x,'present.')
        endif
      endif
c
c     Back up.
c
      backspace ninpts
c
      if (.not.qgexrd) go to 575
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
          go to 570
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
 1620       format(/' * Error - (XCON6/rd6d8) Have exceeded the',
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
c       Read the separator line following the line containing
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
 1070     format(/' * Error - (XCON6/rd6d8) Read the following line',
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
     $  nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
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
     $  nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
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
     $  nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
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
                write (nttyo,1290) jetmax,ugexp(ne)(1:j2),je
 1290           format(/' * Error - (XCON6/rd6d8) Have exceeded the',
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
              go to 560
            endif
c
c           Have found the end of the current block.
c
            go to 570
          endif
c
c         Read the separator line following the line containing
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
                  j3 = ilnobl(ugexj(je,ne))
                  write (nttyo,1820) netmax,ugexj(je,ne)(1:j3),
     $            ugexp(ne)(1:j2),n
 1820             format(/' * Error - (XCON6/rd6d8) Have exceeded the',
     $            ' maximum number of ',i3,/7x,'reactions for a site',
     $            ' belonging to a generic ion exchange',/7x,'phase',
     $            ' while reading the data for site ',a,' of exchange',
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
                  go to 550
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
                go to 560
              endif
c
c             Have found the end of the current block.
c
              go to 570
            endif
c
c           Read the separator line following the line containing
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
  550     continue
        enddo
  560   continue
      enddo
  570 continue
  575 continue
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
        if (ustr(1:8) .eq. '--------') go to 590
c
        call locase(ustr)
        ustrn = ustr
        call locase(ustrn)
        if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') ustr = 'None'
        if (ustr(1:5) .eq. 'None ') go to 580
c
        n = n + 1
c
        if (n .gt. nxmdmx) then
          write (nttyo,1500) nxmdmx
 1500     format(/' * Error - (XCON6/rd6d8) Have too many nxmod',
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
 1510   format(/" * Error - (XCON6/rd6d8) Don't recognize the",
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
  580   continue
      enddo
c
  590 nxmod = n
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
          write (nttyo,1705) uheadx(1:j2)
 1705     format(/' * Error - (XCON6/rd6d8) The iopg option switch',
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
 1712     format(/' * Error - (XCON6/rd6d8) The iopg option switch',
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
 1720     format(/' * Error - (XCON6/rd6d8) The iopg option switch',
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
 1730       format(/' * Error - (XCON6/rd6d8) Have too many option',
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
 1740       format(/' * Error - (XCON6/rd6d8) The following iopg',
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
 1750     format(/' * Error - (XCON6/rd6d8) The iopg option switch',
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
          if (index(uheadx(1:j2),ustr(1:j3)) .le. 0) then
            j4 = ilnobl(uopgtx(n))
            write (nttyo,1760) uheadx(1:j2),ustr(1:j3),uopgtx(n)(1:j4)
 1760       format(/' * Error - (XCON6/rd6d8) The iopg option switch',
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
 1770         format(/" * Error - (XCON6/rd6d8) Don't recognize the",
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
 1780     format(/' * Warning - (XCON6/rd6d8) No option choice was',
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
 1794     format(/' * Warning - (XCON6/rd6d8) More than one option',
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
c     Matrix Index Limits.
c
c     Read the block title from a two-line header.
c
      uheadx = 'Matrix Index Limits'
      nfldtx = 1
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c     Read the number of chemical elements (kct) from a one-line header.
c
      uheadx = 'No. of chem. elements'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 999
      kct = ivar
c
      if (kct .gt. nctmax) then
        write (nttyo,3240) nctmax
 3240   format(/' * Error - (XCON6/rd6d8) Have too many chemical',
     $  /7x,'elements present. The code is only dimensioned',
     $  /7x,'for ',i3,' elements. Reduce the number of elements',
     $  /7x,'or increase the dimensioning parameter nctpar.')
        go to 990
      endif
c
c     Read the number of basis species (kbt) from a one-line header.
c
      uheadx = 'No. of basis species'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 999
      kbt = ivar
c
      if (kbt .gt. nbtmax) then
        write (nttyo,3250) nbtmax
 3250   format(/' * Error - (XCON6/rd6d8) Have too many basis',
     $  /7x,'species present. The code is only dimensioned',
     $  /7x,'for ',i3,' basis species. Reduce the number of elements',
     $  /7x,'or increase the dimensioning parameter nctpar.')
        go to 990
      endif
c
c     Read the matrix index of the last pure min. (kmt) from a
c     one-line header.
c
      uheadx = 'Index of last pure min.'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 999
      kmt = ivar
c
c     Read the matrix index of the last solid solution component
c     from a one-line header.
c
      uheadx = 'Index of last sol-sol.'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 999
      kxt = ivar
c
c     Read the matrix size (kdim) from a one-line header.
c
      uheadx = 'Matrix size'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 999
      kdim = ivar
c
      if (kdim .gt. kmax) then
        write (nttyo,3260) kmax
 3260   format(/' * Error - (XCON6/rd6d8) Have too many matrix',
     $  /7x,'variables. The code is only dimensioned for ',i3,
     $  /7x,'matrix variables. Reduce the number of such variables',
     $  /7x,'or increase the dimensioning parameter kpar.')
        go to 990
      endif
c
c     Read the PRS data flag (kprs) from a two-line header
c     (the separator line is the end of the current block).
c
      uheadx = 'PRS data flag'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr = ufield(2)
      call chrint(ivar,nttyo,qrderr,ustr)
      if (qrderr) go to 999
      kprs = ivar
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Species for which mass balances are defined.
c
c     Read the first part of the block title from a one-line header.
c
      uheadx = 'Mass Balance Species (Matrix Row Variables)'
      nfldtx = 3
      call rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
      if (qrderr) go to 999
c
c     Read the second part of the block title from a two-line header.
c
      uheadx = '(ubmtbi(n))'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
      nbi = 0
c
c     Read the first line.
c
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
c
c     The label below marks a return point for processing subsequent
c     lines in the current data block.
c
  910 ustr = ufield(1)
c
c     There are no data remaining in the current block if a
c     separator line has been encountered.
c
      if (ustr(1:8) .eq. '--------') go to 915
c
c     Check for no species for which mass balances are defined.
c
      ustrn = ustr
      call locase(ustrn)
      if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ') ustr = 'None'
      if (ustr(1:5) .eq. 'None ') then
        nbi = 0
        nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        ustr = ufield(1)
        if (ustr(1:8) .eq. '--------') go to 915
      endif
c
      nfldtx = 3
      if (nfldt .ne. nfldtx) then
        j2 = ilnobl(uline1)
        j2 = min(j2,70)
        write (nttyo,5010) nfldt,nfldtx,uline1(1:j2)
 5010   format(/' * Error - (XCON6/rd6d8) Found ',i2,' fields',
     $  ' where ',i2,/7x,'were expected on the line which begins with',
     $  /7x,'"',a,'".')
        go to 990
      endif
      ux48 = ufield(1)
c
      nbi = nbi + 1
c
      if (nbi .gt. nbtmax) then
        j2 = ilnobl(ux48)
        write (nttyo,5020) nbtmax,ux48(1:j2)
 5020   format(/' * Error - (XCON6/rd6d8) The number of mass balance',
     $  ' species read',/7x,'from the input file exceeded the',
     $  ' maximum of ',i3,' while',/7x,'trying to read data for',
     $  ' the species ',a,'.',/7x,'Increase the dimensioning',
     $  ' parameter nbtpar.')
        go to 990
      endif
c
      ubmtbi(nbi) = ux48
c
      ustr = ufield(2)
      ustrn = ustr
      call locase(ustrn)
      if (ustrn(1:16) .eq. 'moles           ') then
        jflgi(nbi) = 0
      elseif (ustrn(1:16) .eq. 'make non-basis  ') then
        jflgi(nbi) = 30
      else
        j2 = ilnobl(ustr)
        write (nttyo,5030) ustr(1:j2)
 5030   format(/" * Error - (XCON6/rd6d8) Can't identify the",
     $  /7x,'following Mass-Balance-Species Units/constraints string- ',
     $  a,'.',
     $  /7x,'This must be one of "Moles", or "Make non-basis".')
        go to 990
      endif
c
c     Read the next line. Go back to process it.
c
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      go to 910
c
  915 continue
      nbti = nbi
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Mass Balance Totals (moles).
c
c     Read the first part of the block title from a two-line header.
c
      uheadx = 'Mass Balance Totals (moles)'
      nfldtx = 1
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c     Read the next part of the block title from a one-line header.
c
      uheadx = 'Basis species (info. only)'
      nfldtx = 3
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
c
c     Read the second part of the block title from a two-line header.
c
      uheadx = '(ubmtbi(n))'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c     Read the lines.
c
      do nbi = 1,nbti
        nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
c
        ustr = ufield(1)
        ux48 = ubmtbi(nbi)
        j2 = ilnobl(ustr)
        j3 = ilnobl(ux48)
        if (ustr(1:32) .ne. ux48(1:32)) then
          write (nttyo,5040) ux48(1:j3),ustr(1:j2)
 5040     format(/' * Error - (XCON6/rd6d8) Found that the mass',
     $    ' balance species name',/7x,'"',a,'".',/7x,"doesn't match",
     $    ' the mass balance totals species name tag',/7x,'"',a,'"',
     $    /7x,'(which has a maximum length of only 32 characters).')
          qrderr = .true.
          go to 999
        endif
c
        ustr = ufield(2)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        mtbi(nbi) = var
        ustr = ufield(3)
        call chreal(nttyo,qrderr,ustr,var)
        if (qrderr) go to 999
        mtbaqi(nbi) = var
c
      enddo
c
c     Read the electrical imbalance from a two-line header.
c
      uheadx = 'Electrical imbalance'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
      ustr24 = ufield(1)
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
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      electr = var
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
          go to 925
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
  925 nobswt = nobsw
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Matrix column variables and values.
c
c     Read the first part of the block title from a two-line header.
c
      uheadx = 'Matrix Column Variables and Values'
      nfldtx = 1
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c     Read the second part of the block title from a two-line header.
c
      uheadx = 'Basis species (uzveci(n))'
      nfldtx = 3
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
      krow = 0
c
c     Read the first line.
c
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
c
c     The label below marks a return point for processing subsequent
c     lines in the current data block.
c
  960 ustr = ufield(1)
c
c     There are no data remaining in the current block if a
c     separator line has been encountered.
c
      if (ustr(1:8) .eq. '--------') go to 965
c
      nfldtx = 3
      if (nfldt .ne. nfldtx) then
        j2 = ilnobl(uline1)
        j2 = min(j2,70)
        write (nttyo,5010) nfldt,nfldtx,uline1(1:j2)
        go to 990
      endif
      ux48 = ufield(1)
c
      krow = krow + 1
c
      if (krow .gt. kmax) then
        j2 = ilnobl(ux48)
        write (nttyo,5050) kmax,ux48(1:j2)
 5050   format(/' * Error - (XCON6/rd6d8) The number of matrix',
     $  ' column variables',/7x,'and values read from the input file',
     $  ' exceeded the maximum',/7x,'of ',i3,' while trying to read',
     $  ' data for the basis species',/7x,a,'. Increase the',
     $  ' dimensioning parameter kmax.')
        go to 990
      endif
c
      uzveci(krow) = ux48
      ustr = ufield(2)
      call chreal(nttyo,qrderr,ustr,var)
      if (qrderr) go to 999
      zvclgi(krow) = var
c
c     Read the next line. Go back to process it.
c
      nfldtx = 0
      call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
      if (qrderr) go to 999
      go to 960
c
  965 continue
      kdim = krow
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Phases and species in the PRS.
c
c     Read a two-line header for the block.
c
      uheadx = 'Phases and Species in the PRS'
      nfldtx = 1
      call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
      if (qrderr) go to 999
c
c     Loop on phases.
c
      nsi = 0
      npi = 0
      do nn = 1,nprpmx + 1
c
c       Read a line. If the block has not been completely read,
c       this contains the name of a phase, and a sub-block
c       for that phase follows. Otherwise, this line is the first line
c       of the next block.
c
        nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        ustr = ufield(1)
        uheadx = 'Phase'
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
        if (.not.qnonep) then
          npi = npi + 1
c
          if (npi .gt. nprpmx) then
            write (nttyo,3640) nprpmx
 3640       format(/' * Error - (XCON6/rd6d8) Have too many phases',
     $      /7x,'in the physically removed system (PRS) on the',
     $      /7x,'input file. The code is only dimensioned for ',i3,
     $      /7x,'such phases. Reduce the number of such phases',
     $      /7x,'or increase the dimensioning parameter nprppa.')
            go to 990
          endif
c
          uprphi(npi )= ustr
        endif
c
c       Read the separator line following the line containing
c       the phase name.
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
c       Read the number of moles from a two-line header.
c
        uheadx = '->'
        nfldtx = 4
        call rdd2lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uheadr,uheadx,uline1,uline2,ulscr)
        if (qrderr) go to 999
        ustr24 = ufield(2)
        uheadx = 'No. of Moles'
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
        if (.not.qnonep) mprphi(npi) = var
c
c       Read the first title line.
c
        nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        ustr = ufield(1)
        uheadx = '--->'
        call locase(ustr)
        call locase(uheadx)
        j2 = ilnobl(ustr)
        j3 = ilnobl(uheadx)
        if (ustr(1:j2) .eq. uheadx(1:j3)) then
          ustr24 = ufield(2)
          uheadx = 'Species'
          call locase(ustr24)
          call locase(uheadx)
          j2 = ilnobl(ustr24)
          j3 = ilnobl(uheadx)
          if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
            qrderr = .true.
            go to 999
          endif
        endif
c
c       Read the second title line.
c
        nfldtx = 0
        call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $  nttyo,qrderr,ufield,uline1,ulscr)
        if (qrderr) go to 999
        ustr = ufield(1)
        uheadx = '--->'
        call locase(ustr)
        call locase(uheadx)
        j2 = ilnobl(ustr)
        j3 = ilnobl(uheadx)
        if (ustr(1:j2) .eq. uheadx(1:j3)) then
          ustr24 = ufield(2)
          uheadx = '(uprspi(i,n))'
          call locase(ustr24)
          call locase(uheadx)
          j2 = ilnobl(ustr24)
          j3 = ilnobl(uheadx)
          if (ustr24(1:j2) .ne. uheadx(1:j3)) then
            write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
            qrderr = .true.
            go to 999
          endif
        endif
c
c       Read the separator line following the second title line.
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
c       Loop on component species.
c
        iki = 0
        do jj = 1,iktmax + 1
c
c         Read a line.
c
          nfldtx = 0
          call rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $    nttyo,qrderr,ufield,uline1,ulscr)
          if (qrderr) go to 999
c
          ustr = ufield(1)
          uheadx = '--------'
          if (ustr(1:8) .ne. uheadx(1:8)) then
            ustr24 = ufield(1)
            uheadx = '--->'
            call locase(ustr24)
            call locase(uheadx)
            j2 = ilnobl(ustr24)
            j3 = ilnobl(uheadx)
            if (ustr24(1:j2) .ne. uheadx(1:j3)) then
              write (nttyo,1020) uheadx(1:j3),ustr24(1:j2)
              qrderr = .true.
              go to 999
             endif
            ustr = ufield(2)
            ustrn = ustr
            call locase(ustrn)
            if (ustrn(1:5).eq.'none ' .or. ustr(1:1).eq.' ')
     $      ustr = 'None'
            qnones = ustr(1:5).eq.'None '
c
            qokay = .not.qnones .and. .not.qnonep
            if (qokay) then
c
c             Have found another species.
c
              iki = iki + 1
c
              if (iki .gt. iktmax) then
                j2 = ilnobl(uprphi(npi))
                write (nttyo,3680) uprphi(npi)(1:j2),iktmax
 3680           format(/' * Error - (XCON6/rd6d8) Have too many',
     $          ' end-members',/7x,'in the PRS solid solution ',a,'.',
     $          /7x,'The code is only dimensioned for ',
     $          i4,' end-members per',/7x,'solid solution. Reduce',
     $          ' the number of end-members or',
     $          /7x,'increase the dimensioning parameter iktpar.')
                go to 990
              endif
c
              nsi = nsi + 1
c
              if (nsi .gt. nprsmx) then
                write (nttyo,2000) nprsmx,nsi
 2000           format(/' * Error - (XCON6/rd6d8) Have exceeded the',
     $          ' maximum number of ',i3,/7x,'species while reading',
     $          ' the data for phases and species in the PRS',/7x,
     $          ' Increase the dimensioning parameter',
     $          ' nprsmx to at least ',i3,'.')
                go to 990
              endif
c
              uprspi(nsi) = ufield(2)
              ustr = ufield(3)
              call chreal(nttyo,qrderr,ustr,var)
              if (qrderr) go to 999
              mprspi(nsi) = var
            endif
          else
c
c           Have found end of species.
c
            go to 335
          endif
c
        enddo
  335   continue
      enddo
  340 continue
      nprpti = npi
      nprsti = nsi
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
