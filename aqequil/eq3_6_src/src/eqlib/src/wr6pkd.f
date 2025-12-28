      subroutine wr6pkd(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,
     $ dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,
     $ dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,
     $ dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,iact,
     $ ibsrti,iesrti,ietmax,igerti,iktmax,imchmx,imech,iodb,iopg,
     $ iopr,iopt,itermx,ixrti,jcode,jetmax,jflgi,jgerti,jgext,
     $ jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,
     $ kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,
     $ mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,
     $ nertmx,net,netmax,newin,nffg,nffgmx,ngexrt,nobswt,nodbmx,
     $ nopgmx,noprmx,noptmx,nordmx,nprpmx,nprpti,nprsmx,nprsti,
     $ nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,
     $ ntitmx,ntrymx,nttkmx,nxmdmx,nxmod,nxopex,nxopmx,nxopt,
     $ nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,pressb,
     $ pressi,ptk,qgexsh,rkb,rxbari,sfcar,ssfcar,tempcb,tempci,
     $ tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,trkb,ttk,
     $ ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,ugersi,
     $ ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,ureac,
     $ usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,uxopt,
     $ uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,
     $ xlkmod,xvfgex,zgexj,zvclgi)
c
c     This subroutine writes the EQ6 pickup file in menu-style ("D")
c     format for version 8.0.
c
c     This subroutine is a near-clone of XCON6/wr6d8.f.
c
c     The calling sequence of this subroutine is identical to that of
c     EQLIB/wr6pkw.f, XCON6/wr6d8.f, and XCON6.wr6w8.f.
c
c     The calling sequence of this subroutine is identical to that of
c     EQ6/rd6ind.f, EQ6/rd6inw.f, XCON6/rd6d8.f, and XCON6/rd6w8.f,
c     except that newin is added and ninpts, nprob, noutpt, qend,
c     and qrderr are deleted.
c
c     To correctly perform the write function in EQ6, certain variables
c     and arrays must first be updated by calling EQ6/setpk6.f.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ6/path.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c
c     Principal output:
c
c       None
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
      integer newin
c
      integer iact(imchmx,2,nrctmx),ibsrti(nsrtmx),iesrti(nsrtmx),
     $ igerti(jetmax,nertmx),imech(2,nrctmx),iodb(nodbmx),iopg(nopgmx),
     $ iopr(noprmx),iopt(noptmx),ixrti(nxrtmx),jcode(nrctmx),
     $ jflgi(nbtmax),jgerti(nertmx),jgext(netmax),jreac(nrctmx),
     $ kxmod(nxmdmx),ndact(imchmx,2,nrctmx),ngexrt(jetmax,netmax),
     $ nrk(2,nrctmx),nsk(nrctmx)
c
      integer itermx,jpress,jtemp,kbt,kct,kdim,kmt,kprs,ksplmx,ksppmx,
     $ kstpmx,kxt,nbti,nffg,nert,net,nobswt,nprpti,nprsti,nrct,nsbswt,
     $ nsrt,ntitl1,ntitl2,ntrymx,nxmod,nxopex,nxopt,nxrt
c
      logical qgexsh
c
      character(len=80) utitl1(ntitmx),utitl2(ntitmx)
      character(len=56) ugexr(ietmax,jetmax,netmax)
      character(len=48) ubmtbi(nbtmax),uobsw(2,nbtmax),uprspi(nprsmx),
     $ usbsw(2,nbtmax),uxmod(nxmdmx),uzveci(kmax)
      character(len=24) ubsri(nbt1mx,nsrtmx),ucxri(iktmax,nxrtmx),
     $ udac(ndctmx,imchmx,2,nrctmx),uffg(nffgmx),ugermo(nertmx),
     $ ugersi(ietmax,jetmax,nertmx),ugexmo(netmax),ugexp(netmax),
     $ uprphi(nprpmx),ureac(nrctmx),uxcat(nxopmx),uxopex(nxpemx)
      character(len=8) uesri(nctmax,nsrtmx),ugerji(jetmax,nertmx),
     $ ugexj(jetmax,netmax),uhfgex(ietmax,jetmax,netmax),
     $ uvfgex(ietmax,jetmax,netmax),uxkgex(ietmax,jetmax,netmax),
     $ uxopt(nxopmx)
c
      real(8) cbsri(nbt1mx,nsrtmx),cdac(ndctmx,imchmx,2,nrctmx),
     $ cesri(nctmax,nsrtmx),cgexj(jetmax,netmax),
     $ csigma(imchmx,2,nrctmx),eact(imchmx,2,nrctmx),
     $ egersi(ietmax,jetmax,nertmx),fkrc(nrctmx),
     $ hact(imchmx,2,nrctmx),modr(nrctmx),moffg(nffgmx),morr(nrctmx),
     $ mprphi(nprpmx),mprspi(nprsmx),mtbaqi(nbtmax),mtbi(nbtmax),
     $ mwtges(netmax),ptk(nptkmx),rkb(imchmx,2,nrctmx),
     $ rxbari(iktmax,nxrtmx),sfcar(nrctmx),ssfcar(nrctmx),
     $ tgexp(netmax),trkb(imchmx,2,nrctmx),ttk(nttkmx),vreac(nrctmx),
     $ xgersi(ietmax,jetmax,nertmx),xhfgex(ietmax,jetmax,netmax),
     $ xlkffg(nffgmx),xlkgex(ietmax,jetmax,netmax),xlkmod(nxmdmx),
     $ xvfgex(ietmax,jetmax,netmax),zvclgi(kmax),zgexj(jetmax,netmax)
c
      real(8) awmaxi,awmini,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,
     $ dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,
     $ dlxplo,dlxprl,dlxprn,ehmaxi,ehmini,electr,o2maxi,o2mini,phmaxi,
     $ phmini,pressb,pressi,tempcb,tempci,timmxi,tistti,tolbt,toldl,
     $ tolsat,tolxsf,ximaxi,xistti
c
c-----------------------------------------------------------------------
c
      include 'eqlib/eqlk8.h'
      include 'eqlib/eqlo8.h'
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,iei,ij,iki,j,jcox,jd,je,jei,jfl,jfli,jj,jnrk,jrex,j2,
     $ j3,j4,j5,j6,k,krow,kxmd,n,nbi,nci,ne,ner,npi,nrc,nr1,nsi,nsr,nxr
c
      integer nttyo
c
      integer ilnobl
c
      logical qgexbf,qstop,qx
c
      character(len=80) ux80
      character(len=48) ux48,ux48a,ux48b
      character(len=24) ux24
      character(len=16) ux16
      character(len=8) ux1(5)
      character(len=8) uxn,uxo,uxs,ux8
c
      real(8) px3,px4,px5,ptk14,ptk15,sx1,sx2,sx3,tx1,tx2,tx3,tx4,
     $ ttk12,ttk22,ttk14,ttk24,xx
c
c-----------------------------------------------------------------------
c
      data nttyo  /6/
c
c-----------------------------------------------------------------------
c
c     Check some dimensioning parameters.
c
      qstop = .false.
c
      if (ndbxpa .ne. nodbmx) then
        write (nttyo,1000) ndbxpa,nodbmx
 1000   format(/' * Error - (EQLIB/wr6pkd) The dimensioning parameter',
     $  ' for the',/7x,'number of iodb debugging print option switches',
     $  ' with string definitions',/7x,'(ndbxpa) has a value of ',
     $  i3,', but the dimensioning',/7x,'parameter for the number',
     $  ' of such switches (nodbpa) has a',/7x,'value of ',i3,'.')
        qstop = .true.
      endif
c
      if (npgxpa .ne. nopgmx) then
        write (nttyo,1010) npgxpa,nopgmx
 1010   format(/' * Error - (EQLIB/wr6pkd) The dimensioning parameter',
     $  ' for the',/7x,'number of iopg activity coefficient option',
     $  ' switches with string definitions',/7x,'(npgxpa) has a value',
     $  ' of ' ,i3,', but the dimensioning',/7x,'parameter for the',
     $  ' number of such switches (nopgpa) has a',/7x,'value of ',i3,
     $  '.')
        qstop = .true.
      endif
c
      if (nprxpa .ne. noprmx) then
        write (nttyo,1020) nprxpa,noprmx
 1020   format(/' * Error - (EQLIB/wr6pkd) The dimensioning parameter',
     $  ' for the',/7x,'number of iopr print option switches',
     $  ' with string definitions',/7x,'(nprxpa) has a value of ',
     $  i3,', but the dimensioning',/7x,'parameter for the number',
     $  ' of such switches (noprpa) has a',/7x,'value of ',i3,'.')
        qstop = .true.
      endif
c
      if (nptxpa .ne. noptmx) then
        write (nttyo,1030) nptxpa,noptmx
 1030   format(/' * Error - (EQLIB/wr6pkd) The dimensioning parameter',
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
c     Main title.
c
      write (newin,1050)
 1050 format('|',78('-'),'|')
      write (newin,1070)
 1070 format('| Main Title',13x,'| (utitl1(n))',41x,'|')
      write (newin,1050)
      do n = 1,ntitl1
        write (newin,1090) utitl1(n)
 1090   format('|',a78,'|')
      enddo
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Temperature.
c
      ux1(1) = ' '
      ux1(2) = ' '
      ux1(3) = ' '
      ux1(4) = ' '
      tx1 = 0.
      tx2 = 0.
      tx3 = 0.
      tx4 = 0.
      ttk12 = 0.
      ttk22 = 0.
      ttk14 = 0.
      ttk24 = 0.
      if (jtemp .eq. 0) then
        ux1(1) = 'x'
        tx1 = tempcb
      elseif (jtemp .eq. 1) then
        ux1(2) = 'x'
        tx2 = tempcb
        ttk12 = ttk(1)
      elseif (jtemp .eq. 2) then
        ux1(3) = 'x'
        tx3 = tempcb
        ttk22 = ttk(1)
      elseif (jtemp .eq. 3) then
        ux1(4) = 'x'
        tx4 = tempcb
        ttk14 = ttk(1)
        ttk24 = ttk(2)
      endif
      write (newin,1130) ux1(1),tx1,ux1(2),tx2,ttk12,ux1(3),tx3,ttk22,
     $ ux1(4),tx4,ttk24,ttk14
 1130 format('|Temperature option (jtemp):',51x,'|',
     $ /'|  [',a1,'] ( 0) Constant temperature:',46x,'|',
     $ /'|             Value (C)         |',1pe12.5,'| (tempcb)',24x,
     $ '|',/'|  [',a1,'] ( 1) Linear tracking in Xi:',t80,'|',
     $ /'|             Base Value (C)    |',e12.5,'| (tempcb)',24x,'|',
     $ /'|             Derivative        |',e12.5,'| (ttk(1))',24x,'|',
     $ /'|  [',a1,'] ( 2) Linear tracking in time:',t80,'|',
     $ /'|             Base Value (C)    |',e12.5,'| (tempcb)',24x,'|',
     $ /'|             Derivative        |',e12.5,'| (ttk(1))',24x,'|',
     $ /'|  [',a1,'] ( 3) Fluid mixing tracking (fluid 2 = special',
     $ ' reactant):',t80,'|',
     $ /'|             T of fluid 1 (C)  |',e12.5,'| (tempcb)',24x,'|',
     $ /'|             T of fluid 2 (C)  |',e12.5,'| (ttk(2))',24x,'|',
     $ /'|             Mass ratio factor |',e12.5,'| (ttk(1))',24x,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Pressure.
c
      ux1(1) = ' '
      ux1(2) = ' '
      ux1(3) = ' '
      ux1(4) = ' '
      ux1(5) = ' '
      px3 = 0.
      px4 = 0.
      px5 = 0.
      ptk14 = 0.
      ptk15 = 0.
      if (jpress .eq. 0) then
        ux1(1) = 'x'
      elseif (jpress .eq. 1) then
        ux1(2) = 'x'
      elseif (jpress .eq. 2) then
        ux1(3) = 'x'
        px3 = pressb
      elseif (jpress .eq. 3) then
        ux1(4) = 'x'
        px4 = pressb
        ptk14 = ptk(1)
      elseif (jpress .eq. 4) then
        ux1(5) = 'x'
        px5 = pressb
        ptk15 = ptk(1)
      endif
      write (newin,1150) ux1(1),ux1(2),ux1(3),px3,ux1(4),px4,ptk14,
     $ ux1(5),px5,ptk15
 1150 format('|Pressure option (jpress):',53x,'|',
     $ /'|  [',a1,'] ( 0) Follow the data file reference pressure',
     $ ' curve',t80,'|',
     $ /'|  [',a1,'] ( 1) Follow the 1.013-bar/steam-saturation curve',
     $ t80,'|',
     $ /'|  [',a1,'] ( 2) Constant pressure:',49x,'|',
     $ /'|             Value (bars)      |',1pe12.5,'| (pressb)',24x,
     $ '|',/'|  [',a1,'] ( 3) Linear tracking in Xi:',t80,'|',
     $ /'|             Base Value (bars) |',e12.5,'| (pressb)',24x,'|',
     $ /'|             Derivative        |',e12.5,'| (ptk(1))',24x,'|',
     $ /'|  [',a1,'] ( 4) Linear tracking in time:',t80,'|',
     $ /'|             Base Value (bars) |',e12.5,'| (pressb)',24x,'|',
     $ /'|             Derivative        |',e12.5,'| (ptk(1))',24x,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Reactants.
c
      write (newin,2010)
 2010 format('|Reactants (Irreversible Reactions) | (nrct)',t80,'|')
      write (newin,1050)
c
      nxr = 0
      nsr = 0
      ner = 0
c
      do nrc = 1,nrct
        write (newin,2030) ureac(nrc)
 2030   format('|Reactant        |',a24,'| (ureac(n))',t80,'|')
        write (newin,1050)
        write (newin,2050) urcjco(jcode(nrc))
 2050   format('|->|Type         |',a24,'| (urcjco(jcode(n)))',t80,'|')
        write (newin,1050)
        write (newin,2070) urcjre(jreac(nrc))
 2070   format('|->|Status       |',a24,'| (urcjre(jreac(n)))',t80,'|')
        write (newin,1050)
        write (newin,2090) morr(nrc)
 2090   format('|->|Amount remaining (moles) |',1pe12.5,
     $  '| (morr(n))',t80,'|')
        write (newin,1050)
        write (newin,2110) modr(nrc)
 2110   format('|->|Amount destroyed (moles) |',1pe12.5,
     $  '| (modr(n))',t80,'|')
        write (newin,1050)
c
        if (jcode(nrc) .eq. 1) then
c
c         Have a solid solution. Write its composition.
c
          nxr = nxr + 1
          write (newin,2210)
 2210     format('|->|Composition',t80,'|')
          write (newin,1050)
          write (newin,2230)
 2230     format('|--->|Component',15x,'|Mole frac.  | (this is a',
     $    ' table header)',10x,'|')
          write (newin,1050)
          do iki = 1,ixrti(nxr)
            write (newin,2250) ucxri(iki,nxr),rxbari(iki,nxr)
 2250       format('|--->|',a24,'|',1pe12.5,
     $      '| (ucxri(i,n), rxbari(i,n))',9x,'|')
          enddo
          write (newin,1050)
        endif
c
        if (jcode(nrc) .eq. 2) then
c
c         Have a special reactant. Write its molar volume, composition,
c         and reaction.
c
          nsr = nsr + 1
          write (newin,2310) vreac(nrc)
 2310     format('|->|Molar volume (cm3/mol)   |',1pe12.5,
     $    '| (vreac(n))',t80,'|')
          write (newin,1050)
c
          write (newin,2330)
 2330     format('|->|Composition',t80,'|')
          write (newin,1050)
          write (newin,2350)
 2350     format('|--->|Element |Stoich. Number',8x,'| (this is a',
     $    ' table header)',t80,'|')
          write (newin,1050)
          do nci = 1,iesrti(nsr)
            write (newin,2370) uesri(nci,nsr),cesri(nci,nsr)
 2370       format('|--->|',a8,'|',1pe22.15,
     $      '| (uesri(i,n), cesri(i,n))',16x,'|')
          enddo
          write (newin,1050)
c
          write (newin,2390)
 2390     format('|->|Reaction',t80,'|')
          write (newin,1050)
          write (newin,2410)
 2410     format('|--->|Species',17x,'|Reaction Coefficient',2x,
     $    '| (this is a table header)|')
          write (newin,1050)
          do nbi = 1,ibsrti(nsr)
            write (newin,2430) ubsri(nbi,nsr),cbsri(nbi,nsr)
 2430       format('|--->|',a24,'|',1pe22.15,
     $      '| (ubsri(i,n), cbsri(i,n))|')
          enddo
          write (newin,1050)
        endif
c
        if (jcode(nrc) .eq. 5) then
c
c         Have a generic ion exchanger. Write its exchange model
c         and composition.
c
          ner = ner + 1
c
          write (newin,2431) ugermo(ner)
 2431     format('|->|Exchange model  |',a24,'| (ugermo(n))',
     $    21x,'|')
          write (newin,1050)
c
          write (newin,2432)
 2432     format('|->|Composition',t80,'|')
          write (newin,1050)
          do jei = 1,jgerti(ner)
            write (newin,2440) ugerji(jei,ner)
 2440       format('|--->|Exchange site',3x,'|',a8,'| (ugerji(j,n))',
     $      33x,'|')
            write (newin,1050)
            write (newin,2450)
 2450       format('|----->|Component',15x,'|Eq. frac.   | (this is a',
     $      ' table header)',8x,'|')
            write (newin,1050)
            do iei = 1,igerti(jei,ner)
              write (newin,2460) ugersi(iei,jei,ner),egersi(iei,jei,ner)
 2460         format('|----->|',a24,'|',1pe12.5,
     $        '| (ugersi(i,j,n), egersi(i,j,n))',2x,'|')
            enddo
          enddo
          write (newin,1050)
        endif
c
        ux1(1) = ' '
        ux1(2) = ' '
        ux1(3) = ' '
        sx1 = 0.
        sx2 = 0.
        sx3 = 0.
        if (nsk(nrc) .le. 0) then
          ux1(1) = 'x'
          sx1 = sfcar(nrc)
        elseif (nsk(nrc) .eq. 1) then
          ux1(2) = 'x'
          sx2 = ssfcar(nrc)
        elseif (nsk(nrc) .eq. 2) then
          ux1(3) = 'x'
          sx3 = sfcar(nrc)
        endif
c
        write (newin,2610) ux1(1),sx1,ux1(2),sx2,ux1(3),sx3
 2610   format('|->|Surface area option (nsk(n)):',t80,'|',
     $  /'|->|  [',a1,'] ( 0) Constant surface area:',t80,'|',
     $  /'|->|             Value (cm2)       |',1pe12.5,
     $  '| (sfcar(n))',t80,'|',
     $  /'|->|  [',a1,'] ( 1) Constant specific surface area:',t80,
     $  '|',/'|->|             Value (cm2/g)     |',1pe12.5,
     $ '| (ssfcar(n))',t80,'|',
     $  /'|->|  [',a1,'] ( 2) n**2/3 growth law- current surface',
     $  ' area:',t80,'|',
     $  /'|->|             Value (cm2)       |',1pe12.5,
     $  '| (sfcar(n))',t80,'|')
        write (newin,1050)
c
        write (newin,2630) fkrc(nrc)
 2630   format('|->|Surface area factor      |',1pe12.5,
     $  '| (fkrc(n))',t80,'|')
        write (newin,1050)
c
c       Loop over the forward and backward directions.
c
        do jd = 1,2
          write (newin,2640) urcrld(jd),urcnrk(nrk(jd,nrc),jd),jd
 2640     format('|->|',a24,2x,'|',a24,
     $    '| (urcnrk(nrk(',i1,',n)))',t80,'|')
          write (newin,1050)
c
          if (nrk(jd,nrc) .eq. 1) then
c
c           Relative rate law.
c
            do i = 1,3
              write (newin,2730) urcrel(i,jd),rkb(i,jd,nrc),i,jd
 2730         format('|--->|',a24,2x,'|',1pe12.5,
     $        '| (rkb(',i1,',',i1,',n))',t80,'|')
              write (newin,1050)
            enddo
          endif
c
          if (nrk(jd,nrc) .eq. 2) then
c
c           TST rate law.
c
            do i = 1,imech(jd,nrc)
              ux8 = ' '
              write (ux8,'(i2)') i
              j2 = ilnobl(ux8)
              write (newin,2830) ux8(1:j2)
 2830         format('|--->|Mechanism ',a,t80,'|')
              write (newin,1050)
              write (newin,2850) urcrls(jd),csigma(i,jd,nrc),jd
 2850         format('|----->|sigma(i,',a1,',n)            |',1pe12.5,
     $        '| csigma(i,',i1,',n)',t80,'|')
              write (newin,1050)
              write (newin,2870) urcrls(jd),rkb(i,jd,nrc),jd
 2870         format('|----->|k(i,',a1,',n) (mol/cm2/sec)  |',1pe12.5,
     $        '| rkb(i,',i1,',n)',t80,'|')
              write (newin,1050)
              write (newin,2890) trkb(i,jd,nrc),jd
 2890         format('|----->|Ref. Temperature (C)    |',1pe12.5,
     $        '| trkb(i,',i1,',n)',t80,'|')
              write (newin,1050)
c
              ux1(1) = ' '
              ux1(2) = ' '
              ux1(3) = ' '
              sx1 = 0.
              sx2 = 0.
              sx3 = 0.
              if (iact(i,jd,nrc) .le. 0) then
                ux1(1) = 'x'
              elseif (iact(i,jd,nrc) .eq. 1) then
                ux1(2) = 'x'
                sx2 = eact(i,jd,nrc)
              else
                ux1(3) = 'x'
                sx3 = hact(i,jd,nrc)
              endif
              write (newin,2910) jd,ux1(1),ux1(2),sx2,jd,ux1(3),sx3,jd
 2910         format('|----->|Temperature dependence option (iact(i,',
     $        i1,',n)):',27x,'|',/'|----->|  [',a1,'] ( 0) No',
     $        ' temperature dependence',t80,'|',/'|----->|  [',a1,
     $        '] ( 1) Constant activation energy:',t80,'|',
     $        /'|----->|             Value (kcal/mol)  |',1pe12.5,
     $        '| (eact(i,',i1,',n))',t80,'|',/'|----->|  [',a1,
     $        '] ( 2) Constant activation enthalpy:',t80,'|',
     $        /'|----->|             Value (kcal/mol)  |',1pe12.5,
     $        '| (hact(i,',i1,',n))',t80,'|')
              write (newin,1050)
c
              write (newin,2930) jd
 2930         format('|----->|Kinetic activity product species',
     $        ' (ndact(i,',i1,',n))',t80,'|')
              write (newin,1050)
c
              write (newin,2950) urcrls(jd),jd,jd
 2950         format('|------->|Species',41x,'|-N(j,i,',a1,',n)',9x,'|',
     $        /'|------->| (udac(j,i,',i1,',n))',32x,
     $        '| (cdac(j,i,',i1,',n))',t80,'|')
              write (newin,1050)
c
              do j = 1,ndact(i,jd,nrc)
                write (newin,2970) udac(j,i,jd,nrc),cdac(j,i,jd,nrc)
 2970           format('|------->|',a48,'|',1pe12.5,t80,'|')
                write (newin,1050)
              enddo
c
              if (ndact(i,jd,nrc) .le. 0) then
                ux48 = 'None'
                xx = 0.
                write (newin,2970) ux48,xx
                write (newin,1050)
              endif
            enddo
          endif
c
          if (nrk(jd,nrc) .eq. 3) then
c
c           Linear rate law.
c
            i = 1
            ux8 = ' '
            write (ux8,'(i2)') i
            j2 = ilnobl(ux8)
            write (newin,2830) ux8(1:j2)
            write (newin,1050)
            write (newin,2870) urcrls(jd),rkb(i,jd,nrc),jd
            write (newin,1050)
            write (newin,2890) trkb(i,jd,nrc),jd
            write (newin,1050)
c
            ux1(1) = ' '
            ux1(2) = ' '
            ux1(3) = ' '
            sx1 = 0.
            sx2 = 0.
            sx3 = 0.
            if (iact(i,jd,nrc) .le. 0) then
              ux1(1) = 'x'
            elseif (iact(i,jd,nrc) .eq. 1) then
              ux1(2) = 'x'
              sx2 = eact(i,jd,nrc)
            else
              ux1(3) = 'x'
              sx3 = hact(i,jd,nrc)
            endif
            write (newin,2910) jd,ux1(1),ux1(2),sx2,jd,ux1(3),sx3,jd
            write (newin,1050)
          endif
        enddo
      enddo
c
      if (nrct .le. 0) then
        ux24 = 'None'
        xx = 0.
        write (newin,2030) ux24
        write (newin,1050)
        write (newin,2050) ux24
        write (newin,1050)
        write (newin,2070) ux24
        write (newin,1050)
        write (newin,2090) xx
        write (newin,1050)
        write (newin,2110) xx
        write (newin,1050)
      endif
c
      write (newin,3110)
 3110 format('* Valid reactant type strings (urcjco(jcode(n))) are:',
     $ t80,'*')
      ux80(1:1) = '*'
      ux80(2:79) = ' '
      ux80(80:80) = '*'
      j3 = 6
      k = 0
      do jcox = 0,5
        jj = j3 + 27
        ux80(j3:jj) = urcjco(jcox)
        j3 = jj + 1
        k = k + 1
        if (k .ge. 2) then
          write (newin,3130) ux80(1:80)
 3130     format(a)
          ux80(2:79) = ' '
          j3 = 6
          k = 0
        endif
      enddo
      if (k .gt. 0) write (newin,3130) ux80(1:80)
c
      write (newin,1052)
 1052 format('*',78('-'),'*')
c
      write (newin,3150)
 3150 format('* Valid reactant status strings (urcjre(jreac(n))) are:',
     $ t80,'*')
      ux80(1:1) = '*'
      ux80(2:79) = ' '
      ux80(80:80) = '*'
      j3 = 6
      k = 0
      do jrex = -1,2
        jj = j3 + 27
        ux80(j3:jj) = urcjre(jrex)
        j3 = jj + 1
        k = k + 1
        if (k .ge. 2) then
          write (newin,3170) ux80(1:80)
 3170     format(a)
          ux80(2:79) = ' '
          j3 = 6
          k = 0
        endif
      enddo
      if (k .gt. 0) write (newin,3130) ux80(1:80)
c
      write (newin,1052)
c
      do jd = 1,2
        ux24 = urcrld(jd)
        call locase(ux24)
        j2 = ilnobl(ux24)
        write (newin,3190) ux24(1:j2),jd
 3190   format('* Valid ',a,' strings (urcnrk(nrk(',i1,',n))) are:',
     $  t80,'*')
        ux80(1:1) = '*'
        ux80(2:79) = ' '
        ux80(80:80) = '*'
        j3 = 6
        k = 0
        do jnrk = -1,3
          ux24 = urcnrk(jnrk,jd)
          j4 = ilnobl(ux24)
          if (ux24(1:j4) .ne. 'Illegal value') then
            jj = j3 + 27
            ux80(j3:jj) = ux24
            j3 = jj + 1
            k = k + 1
            if (k .ge. 2) then
              write (newin,3210) ux80(1:80)
 3210         format(a)
              ux80(2:79) = ' '
              j3 = 6
              k = 0
            endif
          endif
        enddo
        if (k .gt. 0) write (newin,3130) ux80(1:80)
c
        write (newin,1052)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Starting, minimum, and maximum values of key run parameters.
c
      write (newin,1200)
 1200 format('|Starting, minimum, and maximum values of key run',
     $ ' parameters.',t80,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Starting value of Xi.
c
      write (newin,1210) xistti
 1210 format('|Starting Xi value        |',1pe12.5,'| (xistti)',t80,'|')
      write (newin,1050)
c
c     Maximum value of Xi.
c
      write (newin,1215) ximaxi
 1215 format('|Maximum Xi value         |',1pe12.5,'| (ximaxi)',t80,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Starting value of time.
c
      write (newin,1220) tistti
 1220 format('|Starting time (seconds)  |',1pe12.5,'| (tistti)',t80,'|')
      write (newin,1050)
c
c     Maximum value of time.
c
      write (newin,1225) timmxi
 1225 format('|Maximum time (seconds)   |',1pe12.5,'| (timmxi)',t80,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Minimum value of pH.
c
      write (newin,1230) phmini
 1230 format('|Minimum value of pH      |',1pe12.5,'| (phmini)',t80,'|')
      write (newin,1050)
c
c     Maximum value of pH.
c
      write (newin,1235) phmaxi
 1235 format('|Maximum value of pH      |',1pe12.5,'| (phmaxi)',t80,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Minimum value of Eh (v).
c
      write (newin,1240) ehmini
 1240 format('|Minimum value of Eh (v)  |',1pe12.5,'| (ehmini)',t80,'|')
      write (newin,1050)
c
c     Maximum value of Eh (v).
c
      write (newin,1245) phmaxi
 1245 format('|Maximum value of Eh (v)  |',1pe12.5,'| (ehmaxi)',t80,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Minimum value of log fO2.
c
      write (newin,1250) o2mini
 1250 format('|Minimum value of log fO2 |',1pe12.5,'| (o2mini)',t80,'|')
      write (newin,1050)
c
c     Maximum value of log fO2.
c
      write (newin,1255) o2maxi
 1255 format('|Maximum value of log fO2 |',1pe12.5,'| (o2maxi)',t80,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Minimum value of the activity of water.
c
      write (newin,1260) awmini
 1260 format('|Minimum value of aw      |',1pe12.5,'| (awmini)',t80,'|')
      write (newin,1050)
c
c     Maximum value of the activity of water.
c
      write (newin,1265) awmaxi
 1265 format('|Maximum value of aw      |',1pe12.5,'| (awmaxi)',t80,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum number of steps.
c
      write (newin,1270) kstpmx
 1270 format('|Maximum number of steps  |',i12,'| (kstpmx)',t80,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print interval parameters.
c
      write (newin,1300)
 1300 format('|Print interval parameters.',t80,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Xi print interval.
c
      write (newin,1310) dlxprn
 1310 format('|Xi print interval        |',1pe12.5,'| (dlxprn)',t80,'|')
      write (newin,1050)
c
c     Log Xi print interval.
c
      write (newin,1315) dlxprl
 1315 format('|Log Xi print interval    |',1pe12.5,'| (dlxprl)',t80,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Time print interval.
c
      write (newin,1320) dltprn
 1320 format('|Time print interval      |',1pe12.5,'| (dltprn)',t80,'|')
      write (newin,1050)
c
c     Log time print interval.
c
      write (newin,1325) dltprl
 1325 format('|Log time print interval  |',1pe12.5,'| (dltprl)',t80,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     pH print interval.
c
      write (newin,1330) dlhprn
 1330 format('|pH print interval        |',1pe12.5,'| (dlhprn)',t80,'|')
      write (newin,1050)
c
c     Eh (v) print interval.
c
      write (newin,1335) dleprn
 1335 format('|Eh (v) print interval    |',1pe12.5,'| (dleprn)',t80,'|')
      write (newin,1050)
c
c     Log fO2 print interval.
c
      write (newin,1340) dloprn
 1340 format('|Log fO2 print interval   |',1pe12.5,'| (dloprn)',t80,'|')
      write (newin,1050)
c
c     Activity of water print interval.
c
      write (newin,1345) dlaprn
 1345 format('|aw print interval        |',1pe12.5,'| (dlaprn)',t80,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Steps print interval.
c
      write (newin,1380) ksppmx
 1380 format('|Steps print interval     |',i12,'| (ksppmx)',t80,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Plot interval parameters.
c
      write (newin,1400)
 1400 format('|Plot interval parameters.',t80,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Xi plot interval.
c
      write (newin,1410) dlxplo
 1410 format('|Xi plot interval         |',1pe12.5,'| (dlxplo)',t80,'|')
      write (newin,1050)
c
c     Log Xi plot interval.
c
      write (newin,1415) dlxpll
 1415 format('|Log Xi plot interval     |',1pe12.5,'| (dlxpll)',t80,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Time plot interval.
c
      write (newin,1420) dltplo
 1420 format('|Time plot interval       |',1pe12.5,'| (dltplo)',t80,
     $ '|')
      write (newin,1050)
c
c     Log time plot interval.
c
      write (newin,1425) dltpll
 1425 format('|Log time plot interval   |',1pe12.5,'| (dltpll)',t80,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     pH plot interval.
c
      write (newin,1430) dlhplo
 1430 format('|pH plot interval         |',1pe12.5,'| (dlhplo)',t80,'|')
      write (newin,1050)
c
c     Eh (v) plot interval.
c
      write (newin,1435) dleplo
 1435 format('|Eh (v) plot interval     |',1pe12.5,'| (dleplo)',t80,'|')
      write (newin,1050)
c
c     Log fO2 plot interval.
c
      write (newin,1440) dloplo
 1440 format('|Log fO2 plot interval    |',1pe12.5,'| (dloplo)',t80,'|')
      write (newin,1050)
c
c     Activity of water plot interval.
c
      write (newin,1445) dlaplo
 1445 format('|aw plot interval         |',1pe12.5,'| (dlaplo)',t80,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Steps plot interval.
c
      write (newin,1450) ksplmx
 1450 format('|Steps plot interval      |',i12,'| (ksplmx)',t80,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iopt model option switches.
c     Note: iopt(1) = iopt1, etc.
c
      write (newin,1460)
 1460 format('|Iopt Model Option Switches',
     $ ' ("( 0)" marks default choices)',t80,'|')
      write (newin,1050)
c
      uxo = 'iopt('
      do n = 1,noptmx
        j6 = index(uoptcx(n),'EQ6')
        if (j6 .gt. 0) then
          uxn = ' '
          write (uxn,'(i2)') n
          call lejust(uxn)
          j5 = ilnobl(uxn)
          uxo(6:8) = uxn(1:j5)
          j4 = ilnobl(uxo) + 1
          uxo(j4:j4) = ')'
          j3 = ilnobl(uopttx(n))
          write (newin,1465) uxo(1:j4),uopttx(n)(1:j3)
 1465     format('|',a,' - ',a,':',t80,'|')
          i = iopt(n)
          do j = 1,jptxpa
            ij = ioptox(j,n)
            uxs = ' '
            if (ij .eq. i) uxs = 'x'
            j2 = ilnobl(uoptox(j,n))
            if (j2 .gt. 0) then
              write (newin,1470) uxs,ij,uoptox(j,n)(1:j2)
 1470         format('|  [',a1,'] (',i2,') ',a,t80,'|')
            endif
          enddo
          write (newin,1050)
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iopr print option switches.
c     Note: iopr(1) = iopr1, etc.
c
      write (newin,1490)
 1490 format('|Iopr Print Option Switches',
     $ ' ("( 0)" marks default choices)',t80,'|')
      write (newin,1050)
c
      uxo = 'iopr('
      do n = 1,noprmx
        j6 = index(uoprcx(n),'EQ6')
        if (j6 .gt. 0) then
          uxn = ' '
          write (uxn,'(i2)') n
          call lejust(uxn)
          j5 = ilnobl(uxn)
          uxo(6:8) = uxn(1:j5)
          j4 = ilnobl(uxo) + 1
          uxo(j4:j4) = ')'
          j3 = ilnobl(uoprtx(n))
          write (newin,1465) uxo(1:j4),uoprtx(n)(1:j3)
          i = iopr(n)
          do j = 1,jprxpa
            ij = ioprox(j,n)
            uxs = ' '
            if (ij .eq. i) uxs = 'x'
            j2 = ilnobl(uoprox(j,n))
            if (j2 .gt. 0) then
              write (newin,1470) uxs,ij,uoprox(j,n)(1:j2)
            endif
          enddo
          write (newin,1050)
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iodb debugging print option switches.
c     Note: iodb(1) = iodb1, etc.
c
      write (newin,1495)
 1495 format('|Iodb Debugging Print Option Switches',
     $ ' ("( 0)" marks default choices)',t80,'|')
      write (newin,1050)
c
      uxo = 'iodb('
      do n = 1,nodbmx
        j6 = index(uodbcx(n),'EQ6')
        if (j6 .gt. 0) then
          uxn = ' '
          write (uxn,'(i2)') n
          call lejust(uxn)
          j5 = ilnobl(uxn)
          uxo(6:8) = uxn(1:j5)
          j4 = ilnobl(uxo) + 1
          uxo(j4:j4) = ')'
          j3 = ilnobl(uodbtx(n))
          write (newin,1465) uxo(1:j4),uodbtx(n)(1:j3)
          i = iodb(n)
          do j = 1,jdbxpa
            ij = iodbox(j,n)
            uxs = ' '
            if (ij .eq. i) uxs = 'x'
            j2 = ilnobl(uodbox(j,n))
            if (j2 .gt. 0) then
              write (newin,1470) uxs,ij,uodbox(j,n)(1:j2)
            endif
          enddo
          write (newin,1050)
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Nxopt options.
c
      write (newin,2510)
 2510 format('|Mineral Sub-Set Selection Suppression Options |',
     $ ' (nxopt)',t80,'|')
      write (newin,1050)
c
        write (newin,2530)
 2530 format('|Option  |Sub-Set Defining Species|',
     $ ' (this is a table header)',t80,'|')
      write (newin,1050)
c
      if (nxopt .gt. 0) then
        do n = 1,nxopt
          write (newin,2550) uxopt(n),uxcat(n)
 2550     format('|',a8,'|',a24,'| (uxopt(n), uxcat(n))',t80,'|')
        enddo
      else
        ux8 = 'None'
        ux24 = 'None'
        write (newin,2550) ux8,ux24
      endif
      write (newin,1050)
c
      write (newin,2552)
 2552 format('* Valid mineral sub-set selection suppression option',
     $ ' strings (uxopt(n)) are:',t80,'*')
      ux80(1:1) = '*'
      ux80(2:79) = ' '
      ux80(80:80) = '*'
      j3 = 6
      k = 0
      do n = 1,4
        jj = j3 + 11
        ux80(j3:jj) = uxopti(n)
        j3 = jj + 1
        k = k + 1
        if (k .ge. 6) then
          write (newin,2554) ux80(1:80)
 2554       format(a)
          ux80(2:79) = ' '
          j3 = 6
          k = 0
        endif
      enddo
      if (k .gt. 0) write (newin,2554) ux80(1:80)
c
      write (newin,1052)
c
      write (newin,2570)
 2570 format('|Exceptions to the Mineral Sub-Set Selection Suppression',
     $ ' Options | (nxopex)',t80,'|')
      write (newin,1050)
c
        write (newin,2580)
 2580 format('|Mineral',17x,'| (this is a table header)',t80,'|')
      write (newin,1050)
c
      if (nxopex .gt. 0) then
        do n = 1,nxopex
          write (newin,2590) uxopex(n)
 2590     format('|',a24,'| (uxopex(n))',t80,'|')
        enddo
      else
        ux24 = 'None'
        write (newin,2590) ux24
      endif
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Nffg options.
c
      write (newin,4110)
 4110 format('|Fixed Fugacity Options | (nffg)',t80,'|')
      write (newin,1050)
c
        write (newin,4130)
 4130 format('|Gas',21x,'|Moles to Add |Log Fugacity | --',t80,'|',
     $ /'| (uffg(n))',14x,'| (moffg(n))  | (xlkffg(n)) | --',t80,'|')
      write (newin,1050)
c
      if (nffg .gt. 0) then
        do n = 1,nffg
          write (newin,4150) uffg(n),moffg(n),xlkffg(n)
 4150     format('|',a24,'| ',1pe12.5,'| ',e12.5,'| --',t80,'|')
        enddo
      else
        ux24 = 'None'
        xx = 0.
        write (newin,4150) ux24,xx,xx
      endif
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Numerical Parameters.
c
      write (newin,1510)
 1510 format('|Numerical Parameters',t80,'|')
      write (newin,1050)
c
      write (newin,1520) nordmx
 1520 format('|Max. finite-difference order',t45,'|',1x,i3,8x,
     $ '| (nordmx)',t80,'|')
c
      write (newin,1530) tolbt
 1530 format('|Beta convergence tolerance',t45,'|',1pe12.5,
     $ '| (tolbt)',t80,'|')
c
      write (newin,1550) toldl
 1550 format('|Del convergence tolerance',t45,'|',1pe12.5,
     $ '| (toldl)',t80,'|')
c
      write (newin,1570) itermx
 1570 format('|Max. No. of N-R iterations',t45,'|',1x,i3,8x,
     $ '| (itermx)',t80,'|')
c
      write (newin,1590) tolxsf
 1590 format('|Search/find convergence tolerance',t45,'|',1pe12.5,
     $ '| (tolxsf)',t80,'|')
c
      write (newin,1720) tolsat
 1720 format('|Saturation tolerance',t45,'|',1pe12.5,
     $ '| (tolsat)',t80,'|')
c
      write (newin,1730) ntrymx
 1730 format('|Max. No. of Phase Assemblage Tries',t45,'|',
     $ 1x,i3,8x,'| (ntrymx)',t80,'|')
c
      write (newin,1750) dlxmx0
 1750 format('|Zero order step size (in Xi)',t45,'|',1pe12.5,
     $ '| (dlxmx0)',t80,'|')
c
      write (newin,1770) dlxdmp
 1770 format('|Max. interval in Xi between PRS transfers',t45,'|',
     $ 1pe12.5,'| (dlxdmp)',t80,'|')
c
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Mark the start of the bottom half of the file.
c
      write (newin,7170)
 7170 format('* Start of the bottom half of the input file',t80,'*')
      write (newin,1052)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Secondary title.
c
      write (newin,7190)
 7190 format('| Secondary Title',8x,'| (utitl2(n))',41x,'|')
      write (newin,1050)
      do n = 1,ntitl2
        write (newin,1090) utitl2(n)
      enddo
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Special basis switches.
c
      write (newin,7210)
 7210 format('|Special Basis Switches (for model definition only)',
     $ 7x,'| (nsbswt)',11x,'|')
      write (newin,1050)
c
      do n = 1,nsbswt
        write (newin,7230) usbsw(1,n),usbsw(2,n)
 7230   format('|Replace |',a48,'| (usbsw(1,n))',7x,'|',
     $  /'|   with |',a48,'| (usbsw(2,n))',7x,'|')
        write (newin,1050)
      enddo
c
      if (nsbswt .le. 0) then
        ux48a = 'None'
        ux48b = 'None'
        write (newin,7230) ux48a,ux48b
        write (newin,1050)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Original temperature.
c
      write (newin,7310) tempci
 7310 format('|Original temperature (C) |',1pe12.5,'| (tempci)',t80,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Original pressure.
c
      write (newin,7330) pressi
 7330 format('|Original pressure (bars) |',1pe12.5,'| (pressi)',t80,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Ion exchanger creation.
c
      write (newin,1600)
 1600 format('|Create Ion Exchangers',2x,'| (net)',48x,'|')
      write (newin,1050)
c
c     Set the advisory for exchanger creation blocks to follow or not.
c
      qgexbf = qgexsh .or. net.gt.0
      if (qgexbf) then
        write (newin,1603)
 1603   format('|Advisory: at least one exchanger creation block'
     $ ' follows on this file.',t80,'|')
      else
        write (newin,1605)
 1605   format('|Advisory: no exchanger creation blocks follow on',
     $  ' this file.',
     $  t80,'|')
      endif
c
c     Flag: on further processing, show at least one such block.
c
      ux1(1) = ' '
      if (qgexsh) ux1(1) = 'x'
      write (newin,1607) ux1(1)
 1607 format('|Option: on further processing (writing a pickup file',
     $ ' or running XCON6 on the',t80,'|',/'|present file), force the',
     $ ' inclusion of at least one such block (qgexsh):',t80,'|',
     $ /'|  [',a1,'] (.true.)',t80,'|')
      write (newin,1050)
c
      do ne = 1,net
        write (newin,1610) ugexp(ne)
 1610   format('|Exchanger phase |',a24,'| (ugexp(n))',
     $  25x,'|')
        write (newin,1050)
        write (newin,1620) mwtges(ne)
 1620   format('|->|Mol. wt. (Z)    |',1pg12.5,'| (mwtges(n))',
     $  33x,'|')
        write (newin,1050)
        write (newin,1630) ugexmo(ne)
 1630   format('|->|Exchange model  |',a24,'| (ugexmo(n))',
     $  21x,'|')
        write (newin,1050)
        write (newin,1640) tgexp(ne)
 1640   format('|->|Ref. temp. (C)  |',1pg12.5,'| (tgexp(n))',
     $  34x,'|')
        write (newin,1050)
        do je = 1,jgext(ne)
          write (newin,1650) ugexj(je,ne)
 1650     format('|->|Exchange site   |',a8,'| (ugexj(j,n))',
     $    36x,'|')
          write (newin,1050)
          write (newin,1660) cgexj(je,ne)
 1660     format('|--->|Stoich. number  |',1pg12.5,'| (cgexj(j,n))',
     $    30x,'|')
          write (newin,1050)
          write (newin,1670) zgexj(je,ne)
 1670     format('|--->|Electr. charge  |',1pg12.5,'| (zgexj(j,n))',
     $    30x,'|')
          write (newin,1050)
          do n = 1,ngexrt(je,ne)
            write (newin,1680) ugexr(n,je,ne)
 1680       format('|--->|Reaction        |',a24,'| (ugexr(i,j,n)',
     $      17x,'|')
            write (newin,1050)
            write (newin,1682)
 1682       format('|----->|Parameter |Value',7x,'|Units   | (this',
     $      ' is a table header)',13x,'|')
            write (newin,1050)
            write (newin,1684) xlkgex(n,je,ne),uxkgex(n,je,ne)
 1684       format('|----->|K func.   |',1pe12.5,'|',a8,
     $      '| (xlkgex(i,j,n), uxkgex(i,j,n))',7x,'|')
            write (newin,1686) xhfgex(n,je,ne),uhfgex(n,je,ne)
 1686       format('|----->|DelH0r    |',1pe12.5,'|',a8,
     $      '| (xhfgex(i,j,n), uhfgex(i,j,n))',7x,'|')
            write (newin,1688) xvfgex(n,je,ne),uvfgex(n,je,ne)
 1688       format('|----->|DelV0r    |',1pe12.5,'|',a8,
     $      '| (xvfgex(i,j,n), uvfgex(i,j,n))',7x,'|')
            write (newin,1050)
          enddo
        enddo
      enddo
c
      if (net.le.0 .and. qgexsh) then
        ux24 = 'None'
        ux8 = 'None'
        xx = 0.
        write (newin,1610) ux24
        write (newin,1050)
        write (newin,1620) xx
        write (newin,1050)
        write (newin,1630) ux24
        write (newin,1050)
        write (newin,1640) xx
        write (newin,1050)
        write (newin,1650) ux8
        write (newin,1050)
        write (newin,1660) xx
        write (newin,1050)
        write (newin,1670) xx
        write (newin,1050)
        write (newin,1680) ux24
        write (newin,1050)
        write (newin,1682)
        write (newin,1050)
        write (newin,1684) xx,ux8
        write (newin,1686) xx,ux8
        write (newin,1688) xx,ux8
        write (newin,1050)
      endif
c
      if (net.gt.0 .or. qgexsh) then
        write (newin,1694)
 1694   format('* Valid exchange model strings (ugexmo(n)) are:',t80,
     $  '*')
        ux80(1:1) = '*'
        ux80(2:79) = ' '
        ux80(80:80) = '*'
        j3 = 6
        k = 0
        do n = 1,nvet_par
          if (ugexmv(n)(1:6) .ne. 'ERROR ') then
            jj = j3 + 27
            ux80(j3:jj) = ugexmv(n)
            j3 = jj + 1
            k = k + 1
            if (k .ge. 2) then
              write (newin,1696) ux80(1:80)
 1696         format(a)
              ux80(2:79) = ' '
              j3 = 6
              k = 0
            endif
          endif
        enddo
        if (k .gt. 0) write (newin,1696) ux80(1:80)
c
        write (newin,1052)
c
        write (newin,1690)
 1690   format('* Valid units strings (uxkgex(i,j,n)/uhfgex(i,j,n)/',
     $  'uvfgex(i,j,n)) are:',t80,'*')
        ux80(1:1) = '*'
        ux80(2:79) = ' '
        ux80(80:80) = '*'
        j3 = 6
        k = 0
        do n = 1,4
          jj = j3 + 19
          ux80(j3:jj) = uxfuni(n)
          j3 = jj + 1
          k = k + 1
          if (k .ge. 3) then
            write (newin,1692) ux80(1:80)
 1692         format(a)
            ux80(2:79) = ' '
            j3 = 6
            k = 0
          endif
        enddo
        if (k .gt. 0) write (newin,1692) ux80(1:80)
c
        write (newin,1052)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Nxmod options.
c
      write (newin,5080)
 5080 format('|Alter/Suppress Options  | (nxmod)',45x,'|')
      write (newin,1050)
      write (newin,5082)
 5082 format('|Species',41x,'|Option',10x,'|Alter value |',
     $ /'| (uxmod(n))',37x,'|(ukxm(kxmod(n)))| (xlkmod(n))|')
      write (newin,1050)
c
      do n = 1,nxmod
        write (newin,5090) uxmod(n),ukxm(kxmod(n)),xlkmod(n)
 5090   format('|',a48,'|',a16,'|',1pe12.5,'|')
      enddo
c
      if (nxmod .le. 0) then
        ux48 = 'None'
        ux16 = 'None'
        xx = 0.
        write (newin,5090) ux48,ux16,xx
      endif
c
      write (newin,1050)
c
      write (newin,5092)
 5092 format('* Valid alter/suppress strings (ukxm(kxmod(n))) are:',
     $ t80,'*')
      ux80(1:1) = '*'
      ux80(2:79) = ' '
      ux80(80:80) = '*'
      j3 = 6
      k = 0
      do kxmd = -1,2
        jj = j3 + 19
        ux80(j3:jj) = ukxm(kxmd)
        j3 = jj + 1
        k = k + 1
        if (k .ge. 3) then
          write (newin,5094) ux80(1:80)
 5094       format(a)
          ux80(2:79) = ' '
          j3 = 6
          k = 0
        endif
      enddo
      if (k .gt. 0) write (newin,5094) ux80(1:80)
c
      write (newin,1052)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iopg activity coefficient option switches.
c     Note: iopg(1) = iopg1, etc.
c
      write (newin,7510)
 7510 format('|Iopg Activity Coefficient Option Switches',
     $ ' ("( 0)" marks default choices)',t80,'|')
      write (newin,1050)
c
      uxo = 'iopg('
      do n = 1,nopgmx
        j6 = index(uopgcx(n),'EQ6')
        if (j6 .gt. 0) then
          uxn = ' '
          write (uxn,'(i2)') n
          call lejust(uxn)
          j5 = ilnobl(uxn)
          uxo(6:8) = uxn(1:j5)
          j4 = ilnobl(uxo) + 1
          uxo(j4:j4) = ')'
          j3 = ilnobl(uopgtx(n))
          write (newin,7530) uxo(1:j4),uopgtx(n)(1:j3)
 7530     format('|',a,' - ',a,':',t80,'|')
          i = iopg(n)
          do j = 1,jpgxpa
            ij = iopgox(j,n)
            uxs = ' '
            if (ij .eq. i) uxs = 'x'
            j2 = ilnobl(uopgox(j,n))
            if (j2 .gt. 0) then
              write (newin,7550) uxs,ij,uopgox(j,n)(1:j2)
 7550         format('|  [',a1,'] (',i2,') ',a,t80,'|')
            endif
          enddo
          write (newin,1050)
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Index limits.
c
      write (newin,7610)
 7610 format('|Matrix Index Limits',t80,'|')
      write (newin,1050)
      write (newin,7630) kct
 7630 format('|No. of chem. elements   |',i5,'| (kct)',t80,'|')
      write (newin,7650) kbt
 7650 format('|No. of basis species    |',i5,'| (kbt)',t80,'|')
      write (newin,7670) kmt
 7670 format('|Index of last pure min. |',i5,'| (kmt)',t80,'|')
      write (newin,7690) kxt
 7690 format('|Index of last sol-sol.  |',i5,'| (kxt)',t80,'|')
      write (newin,7710) kdim
 7710 format('|Matrix size             |',i5,'| (kdim)',t80,'|')
      write (newin,7730) kprs
 7730 format('|PRS data flag           |',i5,'| (kprs)',t80,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Species for which mass balances are defined.
c
      write (newin,7750)
 7750 format('|Mass Balance Species (Matrix Row Variables)',5x,
     $ '|Units/Constraint| --',9x,'|',
     $ /'| (ubmtbi(n))',36x,'|(ujf6(jflgi(n)))| --',9x,'|')
      write (newin,1050)
c
      qstop = .false.
      do nbi = 1,nbti
        jfl = jflgi(nbi)
        jfli = 1
        if (jfl .eq. 30) jfli = 2
c
        write (newin,7790) ubmtbi(nbi),ujf6(jfli)
 7790   format('|',a48,'|',a16,'| --',9x,'|')
      enddo
      if (qstop) stop
c
      if (nbti .le. 0) then
        ux48 = 'None'
        ux16 = 'None'
        xx = 0.
        write (newin,7790) ux48,xx,ux16
      endif
c
      write (newin,1050)
c
      write (newin,7810)
 7810 format('* Valid jflag strings (ujf6(jflgi(n))) are:',t80,'*')
      ux80(1:1) = '*'
      ux80(2:79) = ' '
      ux80(80:80) = '*'
      j3 = 6
      k = 0
      do jfli = 1,2
        jj = j3 + 19
        ux80(j3:jj) = ujf6(jfli)
        j3 = jj + 1
        k = k + 1
        if (k .ge. 3) then
          write (newin,7830) ux80(1:80)
 7830     format(a)
          ux80(2:79) = ' '
          j3 = 6
          k = 0
        endif
      enddo
      if (k .gt. 0) write (newin,7830) ux80(1:80)
c
      write (newin,1052)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Mass balance totals.
c
      write (newin,7850)
 7850 format('|Mass Balance Totals (moles)',51x,'|')
      write (newin,1050)
c
      write (newin,7870)
 7870 format('|Basis species (info. only)',6x,'|Equilibrium System',4x,
     $ '|Aqueous Solution',6x,'|',
     $ /'| (ubmtbi(n))',20x,'| (mtbi(n))',12x,'| (mtbaqi(n))',10x,'|')
      write (newin,1050)
c
      do nbi = 1,nbti
        write (newin,7890) ubmtbi(nbi),mtbi(nbi),mtbaqi(nbi)
 7890   format('|',a32,'|',1pe22.15,'|',e22.15,'|')
      enddo
c
      ux48 = 'Electrical imbalance'
      write (newin,7890) ux48,electr,electr
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Ordinary basis switches.
c
      write (newin,7970)
 7970 format('|Ordinary Basis Switches (for numerical purposes only)',
     $ 4x,'| (nobswt)',11x,'|')
      write (newin,1050)
c
      do n = 1,nobswt
        write (newin,7990) uobsw(1,n),uobsw(2,n)
 7990   format('|Replace |',a48,'| (uobsw(1,n))',7x,'|',
     $  /'|   with |',a48,'| (uobsw(2,n))',7x,'|')
        write (newin,1050)
      enddo
c
      if (nobswt .le. 0) then
        ux48a = 'None'
        ux48b = 'None'
        write (newin,7990) ux48a,ux48b
        write (newin,1050)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Matrix column variables and values.
c
      write (newin,7910)
 7910 format('|Matrix Column Variables and Values',44x,'|')
      write (newin,1050)
c
      write (newin,7930)
 7930 format('|Basis species (uzveci(n))',23x,
     $ '|Log moles (zvclgi(n)) | --   |')
      write (newin,1050)
c
      do krow = 1,kdim
        write (newin,7950) uzveci(krow),zvclgi(krow)
 7950   format('|',a48,'|',1pe22.15,'| --   |')
      enddo
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Phases and species in the PRS.
c
      write (newin,8110)
 8110 format('|Phases and Species in the PRS',49x,'|')
      write (newin,1050)
c
      qx = kprs.gt.0 .and. nprpti.gt.0 .and. nprsti.gt.0
      if (qx) then
        nr1 = 1
        do npi = 1,nprpti
          write (newin,8130) uprphi(npi)
 8130     format('|Phase',11x,'|',a24,'| (uprphi(n))',t80,'|')
          write (newin,1050)
          write (newin,8150) mprphi(npi)
 8150     format('|->|No. of Moles',4x,'|',1pe22.15,'| (mprphi(n))',
     $    t80,'|')
          write (newin,1050)
          write (newin,8170)
 8170     format('|--->|Species',17x,'|No. of Moles',10x,
     $    '| --',t80,'|',/'|--->| (uprspi(i,n))',10x,
     $    '| (mprspi(i,n))',8x,'| --',t80,'|')
          write (newin,1050)
c
          do nsi = nr1,nprsti
            if (uprspi(nsi)(25:48) .ne. uprphi(npi)(1:24)) then
              nr1 = nsi
              go to 700
            endif
            write (newin,8190) uprspi(nsi),mprspi(nsi)
 8190       format('|--->|',a24,'|',1pe22.15,'| --',t80,'|')
          enddo
  700     write (newin,1050)
        enddo
      endif
c
      if (.not.qx) then
        ux24 = 'None'
        xx = 0.
        write (newin,8130) ux24
        write (newin,1050)
        write (newin,8150) xx
        write (newin,1050)
        write (newin,8170)
        write (newin,1050)
        write (newin,8190) ux24,xx
        write (newin,1050)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     End of problem.
c
      write (newin,8310)
 8310 format('|End of problem',64x,'|')
      write (newin,1050)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
