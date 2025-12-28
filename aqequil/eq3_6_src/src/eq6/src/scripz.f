      subroutine scripz(abar,acflg,acfw,acfwlg,actlg,actw,actwlg,
     $ affpd,affsd,afrc1,aft1,ah,ahmes,ahnbs,ahrc,alki,alk1,
     $ alk2,awmax,awmin,a3bar,cbsr,cdrsd,cegexs,cesr,conc,conclg,csts,
     $ ctb,cteaq,dvoso,dwoso,egers,egexjc,egexjf,egexpa,egexpc,egexs,
     $ egexw,eh,ehmax,ehmes,ehmin,ehnbs,ehrc,elecsr,electr,fje,fjest,
     $ fo2,fo2lg,fo2lrc,fugac,fugalg,fxi,fxist,iaqsln,iemop,iemop0,
     $ iemos,iemos0,iern1,iern2,iexr,iexrt,ifrn1,ifrn2,ilrn1,ilrn2,
     $ imech,imrn1,imrn2,iopg,iopr,iopt,ipndx1,ixrn1,ixrn2,jcode,jcsort,
     $ jern1,jern2,jexr,jexrt,jflag,jflagd,jflgi,jgext,jgsort,jpflag,
     $ jreac,jsca,jscat,jscr,jscrt,jsflag,jsol,jssort,kbt,kern1,kern2,
     $ kgexsa,km1,kmt,kx1,kxt,kstep,kstpmx,loph,losp,mlmrra,modr,
     $ moph,mophg,mophj,mopht,morr,mosp,mospg,mospj,mospt,mprph,
     $ mprsp,mrgers,mrmlra,mwtrc,mwtsp,narn1,narn2,nat,nbasp,nbaspd,
     $ nbt,ncmpe,ncmpe0,ncmpr,nct,ndrsd,ndrsrd,nelect,nern1,nern2,
     $ nert,net,nfrn1,nfrn2,ngext,ngexsa,ngrn1,ngrn2,ngt,nhydr,nhydx,
     $ nlrn1,nlrn2,nlt,nmrn1,nmrn2,nmrt,nmt,noutpt,no2gaq,npet,npet0,
     $ npt,npts,nrct,nrdxsp,nrk,nrndex,nst,nsts,nstsr,ntf1t,ntf2t,
     $ nxridx,nxrn1,nxrn2,nxrt,nxt,osc,oscst,omega,o2max,o2min,pch,pe,
     $ pemes,penbs,perc,ph,phcl,phmax,phmes,phmin,phnbs,ppmwe,presg,
     $ press,qaft1,qftpr2,qmod,qphcl,qredox,qrho,qriinf,qstopx,qvhfxi,
     $ qvlsow,qzprnt,rho,rhoc,rhowc,rk,rreacn,rreac1,rrelr1,rxbar,
     $ sfcar,sidrph,sidrsp,sigmst,sigmam,ssfcar,tdays,tdsglw,tdspkc,
     $ tdsplc,tempc,thours,time1,timemx,tmins,tolsat,tolxsf,tolxst,
     $ tolxsu,tyears,uelem,ugermo,ugexj,ugexmo,uphase,ureac,uspec,
     $ uxtype,vodrt,voph,vophg,vophj,vopht,vosoct,vosol,vosp,vospg,
     $ vospj,vospt,vreac,wfh2o,wftds,wkgwi,woh2o,wodr,wodrt,woph,
     $ wophg,wophj,wopht,worr,worrt,wosoct,wosol,wosp,wospg,wospj,
     $ wospt,wotds,xbar,xbarlg,xbarw,xbrwlg,xgers,xgexw,xi1,xidump,
     $ ximax,xistsv,xirct,zchar)
c
c     This subroutine writes a detailed description on the output file
c     of the modeled system at the current value of reaction progress
c     (xi1).
c
c     NOTE: This subroutine is designed to carry out a largely pure
c     write function. Very few things to be printed by this subroutine
c     should be computed within it. Unless the data are already
c     available in EQ6/path.f, they should be calculated in
c     EQ6/cdappl.f, which is called by EQ6/path.f prior to printing or
c     plotting. This practice should be followed so that data needed for
c     printing are also available for plotting. It is acceptable to make
c     certain simple calculations in this subroutine, such as taking
c     logarithms or converting units. If need be, these can be repeated
c     elsewhere for plotting.
c
c     Note that qftpr2 is a logical flag denoting a print being made
c     only for a special fluid-centered flow-through open system just
c     after a discontinuity in the derivatives along the reaction
c     path.
c
c     This subroutine is called by:
c
c       EQ6/path.f
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
      include 'eqlib/eqlpar.h'
      include 'eqlib/eqldv.h'
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer noutpt
c
      integer jcode(nrctmx),jreac(nrctmx),nrndex(nrctmx),nxridx(nrctmx)
c
      integer iemop(npetmx),iemop0(npetmx),iemos(nsetmx),iemos0(nsetmx),
     $ iexr(nrctmx),imech(2,nrctmx),iopg(nopgmx),iopr(noprmx),
     $ iopt(noptmx),ipndx1(kmax),jcsort(nstmax),jern1(jetmax,netmax),
     $ jern2(jetmax,netmax),jexr(nrctmx),jflag(nstmax),jflagd(nstmax),
     $ jflgi(nbtmax),jgext(netmax),jgsort(ngtmax),jpflag(nptmax),
     $ jsca(nrctmx),jscr(nrctmx),jsflag(nstmax),jsol(nxtmax),
     $ jssort(nstmax),kern1(netmax),kern2(netmax),kgexsa(ketmax,netmax),
     $ nbasp(nbtmax),nbaspd(nbtmax),ncmpe(2,npetmx),ncmpe0(2,npetmx),
     $ ncmpr(2,nptmax),ndrsd(ndrsmx),ndrsrd(2,nstmax),
     $ ngexsa(ietmax,jetmax,netmax),ngext(jetmax,netmax),nrk(2,nrctmx),
     $ nsts(nstsmx),nstsr(2,nstmax)
c
      integer iern1,iern2,ifrn1,ifrn2,ilrn1,ilrn2,imrn1,imrn2,
     $ ixrn1,ixrn2
c
      integer nat,nbt,nct,net,ngt,nlt,nmt,npt,nst,nxt
c
      integer narn1,narn2,nern1,nern2,nfrn1,nfrn2,ngrn1,ngrn2,
     $ nlrn1,nlrn2,nmrn1,nmrn2,nxrn1,nxrn2
c
      integer iaqsln,iexrt,jexrt,jscat,jscrt,kbt,km1,kmt,kx1,kxt,kstep,
     $ kstpmx,nelect,nert,nhydr,nhydx,nmrt,no2gaq,npet,npet0,npts,nrct,
     $ nrdxsp,ntf1t,ntf2t,nxrt
c
      logical qaft1,qftpr2,qmod,qphcl,qredox,qrho,qriinf,qstopx,
     $ qvhfxi,qvlsow,qzprnt
c
      character*48 uspec(nstmax)
      character*32 uxtype(jsomax)
      character*24 ugermo(nertmx),ureac(nrctmx)
      character*24 ugexmo(netmax),uphase(nptmax)
      character*8 uelem(nctmax),ugexj(jetmax,netmax)
c
      real*8 cbsr(nbt1mx,nsrtmx),cesr(nctmax,nsrtmx),
     $ egers(ietmax,jetmax,nertmx),elecsr(nsrtmx),modr(nrctmx),
     $ mrgers(ietmax,jetmax,nertmx),morr(nrctmx),mwtrc(nrctmx),
     $ rreacn(nrctmx),rreac1(nrctmx),rrelr1(nrctmx),
     $ rxbar(iktmax,nxrtmx),sfcar(nrctmx),ssfcar(nrctmx),vreac(nrctmx),
     $ wodr(nrctmx),worr(nrctmx),xgers(ietmax,jetmax,nertmx)
c
      real*8 acflg(nstmax),actlg(nstmax),affpd(nptmax),affsd(nstmax),
     $ afrc1(nrctmx),ahrc(nbtmax),cdrsd(ndrsmx),
     $ cegexs(ietmax,jetmax,netmax),conc(nstmax),conclg(nstmax),
     $ csts(nstsmx),ctb(nbtmax),cteaq(nctmax),egexjc(jetmax,netmax),
     $ egexjf(jetmax,netmax),egexpa(netmax),egexpc(netmax),
     $ egexs(ietmax,jetmax,netmax),egexw(ketmax,netmax),ehrc(nbtmax),
     $ fo2lrc(nbtmax),fugac(ngtmax),fugalg(ngtmax),loph(nptmax),
     $ losp(nstmax),moph(nptmax),mophg(nptmax),mophj(nptmax),
     $ mopht(nptmax),mosp(nstmax),mospg(nstmax),mospj(nstmax),
     $ mospt(nstmax),mprph(nptmax),mprsp(nstmax),mwtsp(nstmax)
c
      real*8 perc(nbtmax),ppmwe(nctmax),rk(imchmx,2,nrctmx),
     $ sidrph(nptmax),sidrsp(nstmax),voph(nptmax),vophg(nptmax),
     $ vophj(nptmax),vopht(nptmax),vosp(nstmax),vospg(nstmax),
     $ vospj(nstmax),vospt(nstmax),woph(nptmax),wophg(nptmax),
     $ wophj(nptmax),wopht(nptmax),wosp(nstmax),wospg(nstmax),
     $ wospj(nstmax),wospt(nstmax),xbar(nstmax),xbarlg(nstmax),
     $ xgexw(ketmax,netmax),xirct(nrctmx),zchar(nstmax)
c
      real*8 abar,acfw,acfwlg,actw,actwlg,aft1,ah,ahmes,ahnbs,
     $ alki,alk1,alk2,awmax,awmin,a3bar,dvoso,dwoso,eh,ehmax,ehmin,
     $ ehmes,ehnbs,electr,fje,fjest,fo2,fo2lg,fxi,fxist,osc,oscst,omega,
     $ o2max,o2min,pch,pe,penbs,pemes,ph,phcl,phmax,phmes,phmin,phnbs,
     $ presg,press,rho,sigmam,sigmst,tdays,tempc,thours,time1,timemx,
     $ tmins,tolsat,tolxsf,tolxst,tolxsu,tyears,vodrt,vosoct,vosol,
     $ wfh2o,wftds,wkgwi,wodrt,woh2o,worrt,wosoct,wosol,wotds,xbarw,
     $ xbrwlg,xi1,xidump,ximax,xistsv
c
      real(8) mlmrra,mrmlra,rhoc,rhowc,tdsglw,tdspkc,tdsplc
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ie,j,je,j1,j2,j3,kcol,ncount,ne,np,npe,nplast,np1,ns,
     $ nsr,nrc,nr1,nr2,nse,nt
c
      integer ilnobl
c
      logical qrct,qtimmx
c
      character(len=16) ux16a,ux16b,ux16c
c
      real*8 betchb,chdrft,chdrfs,dp,elects,lxx,mxx,sigza,sigzc,sigzi,
     $ sigzia,sigzis,sigzm,sigzma,sigzms,tolspf,wconst,xilg
c
      real*8 tlg
c
c-----------------------------------------------------------------------
c
      wconst = omega/mosp(narn1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1000)
 1000 format(/' - - - - - - - - - - - - - - - - - - - - - - - - -',
     $ ' - - - - - - -',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      xilg = -99999.
      if (xi1 .gt. 0.) xilg = tlg(xi1)
c
      write (noutpt,1010) xi1,xilg
 1010 format(/20x,'Xi= ',1pe12.5,/16x,'Log Xi= ',0pf12.5,/)
c
      if (iopt(2) .ge. 1) then
        if (qriinf) then
          write (noutpt,1020)
 1020     format(/12x,'Time= infinity',/)
        else
          write (noutpt,1030) time1,tmins,thours,tdays,tyears
 1030     format(/12x,'Time= ',1pe10.3,' seconds',
     $    /16x,'= ',1pe10.3,' minutes',/16x,'= ',1pe10.3,' hours',
     $    /16x,'= ',1pe10.3,' days',/16x,'= ',1pe10.3,' years',/)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1050) tempc
 1050 format(/' Temperature= ',f6.2,' C')
c
      dp = press - presg
      if (abs(dp) .le. 1.e-4) then
        write (noutpt,1060) press
 1060   format(/' Pressure= ',1pg12.5,' bars',/)
      else
        write (noutpt,1070) press,presg,dp
 1070   format(/' Pressure= ',1pg12.5,' bars',
     $  /' Data file reference curve pressure= ',g12.5,' bars',
     $  /' Pressure difference= ',g12.5,' bars',/)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qftpr2) then
        write(noutpt,1090)
 1090   format(/' Fluid-centered flow-through open system print',
     $  /' following a discontinuity.')
c
c       If this path is taken, only a partial description of
c       the system is printed.
c
        go to 200
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qzprnt) write (noutpt,1100)
 1100 format(/' Required print point.')
c
      if (xi1.ge.xidump .and. xi1.gt.xistsv) write (noutpt,1110)
 1110 format(/' Required PRS shift.')
c
      if (qaft1) write (noutpt,1120)
 1120 format(/' The total affinity is repeatedly nearly zero.')
c
      if (iexrt .gt. 0) then
        write (noutpt,1130)
 1130   format(/' Exhaustion of the following reactant(s):')
        do i = 1,iexrt
          nrc = iexr(i)
          j2 = ilnobl(ureac(nrc))
          write (noutpt,1140) ureac(nrc)(1:j2)
 1140     format(5x,a)
        enddo
      endif
c
      if (jexrt .gt. 0) then
        write (noutpt,1150)
 1150   format(/' Step size is limited by the reactivation of',
     $  ' the following formerly',/1x,'exhausted reactant(s):',/)
        do j = 1,jexrt
          nrc = jexr(j)
          j2 = ilnobl(ureac(nrc))
          write (noutpt,1140) ureac(nrc)
        enddo
      endif
c
      if (jscat .gt. 0) then
        write (noutpt,1160)
 1160   format(/' Step size is limited by the sign change of',
     $  ' the affinity of the',/1x,'following reactant(s):',/)
        do j = 1,jscat
          nrc = jsca(j)
          j2 = ilnobl(ureac(nrc))
          write (noutpt,1140) ureac(nrc)
        enddo
      endif
c
      if (jscrt .gt. 0) then
        write (noutpt,1170)
 1170   format(/,' Step size is limited by the sign change of',
     $  ' the relative rate of the',/1x,'following reactant(s):',/)
        do j = 1,jscrt
          nrc = jscr(j)
          j2 = ilnobl(ureac(nrc))
          write (noutpt,1140) ureac(nrc)
        enddo
      endif
c
      if (xi1 .le. xistsv) write (noutpt,1200)
 1200 format(/' Start or restart of the run.')
c
      if (qmod) write (noutpt,1210)
 1210 format(/' Have a change in the ES phase assemblage.')
c
      if (xi1 .ge. ximax) write (noutpt,1220)
 1220 format(/' Maximum value of Xi.')
c
      qtimmx = iopt(2) .ge. 1 .and. time1.ge.( (1. - tolxst)*timemx )
      if (qtimmx) write (noutpt,1230)
 1230 format(/' Maximum value of time.')
c
      if (abs(ph - phmin) .le. tolxsu) then
        write (noutpt,1240)
 1240   format(/' Minimum value of pH.')
      elseif (ph. lt. phmin) then
        write (noutpt,1250)
 1250   format(/' pH is less than the minimum value.')
      endif
c
      if (abs(ph - phmax) .le. tolxsu) then
        write (noutpt,1260)
 1260   format(/' Maximum value of pH.')
      elseif (ph. gt. phmax) then
        write (noutpt,1270)
 1270   format(/' pH is greater than the maximum value.')
      endif
c
      if (qredox) then
        if (abs(eh - ehmin) .le. tolxsu) then
          write (noutpt,1280)
 1280     format(/' Minimum value of Eh (v).')
        elseif (eh. lt. ehmin) then
          write (noutpt,1290)
 1290     format(/' Eh (v) is less than the minimum value.')
        endif
c
        if (abs(eh - ehmax) .le. tolxsu) then
          write (noutpt,1300)
 1300     format(/' Maximum value of Eh (v).')
        elseif (eh. gt. ehmax) then
          write (noutpt,1310)
 1310     format(/' Eh (v) is greater than the maximum value.')
        endif
      endif
c
      if (qredox) then
        if (abs(fo2lg - o2min) .le. tolxsu) then
          write (noutpt,1320)
 1320     format(/' Minimum value of log fO2.')
        elseif (fo2lg. lt. o2min) then
          write (noutpt,1330)
 1330     format(/' Log fO2 is less than the minimum value.')
        endif
c
        if (abs(fo2lg - o2max) .le. tolxsu) then
          write (noutpt,1340)
 1340     format(/' Maximum value of log fO2.')
        elseif (fo2lg. gt. o2max) then
          write (noutpt,1350)
 1350     format(/' log fO2 is greater than the maximum value.')
        endif
      endif
c
      if (abs(actw - awmin) .le. tolxsu) then
        write (noutpt,1360)
 1360   format(/' Minimum value of the activity of water.')
      elseif (actw. lt. awmin) then
        write (noutpt,1370)
 1370   format(/' The activity of water is less than the minimum',
     $  ' value.')
      endif
c
      if (abs(actw - awmax) .le. tolxsu) then
        write (noutpt,1380)
 1380   format(/' Maximum value of the activity of water.')
      elseif (actw. gt. awmax) then
        write (noutpt,1390)
 1390   format(/' The activity of water is greater than the maximum',
     $  ' value.')
      endif
c
      if (kstep .ge. kstpmx) write (noutpt,1450)
 1450 format(/' Maximum number of steps.')
c
      if (qvlsow) write (noutpt,1460)
 1460 format(/' Have very nearly exhausted solvent water.')
c
      if (qvlsow) write (noutpt,1470)
 1470 format(/' Have extremely high ionic strength.')
c
      if (qstopx) write (noutpt,1480)
 1480 format(/' Early termination.')
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print tables of data for reactants and reaction rates.
c
      call prtrct(afrc1,aft1,imchmx,imech,iopt,modr,morr,
     $ noptmx,noutpt,nrct,nrctmx,nrk,rk,rreacn,rreac1,rrelr1,sfcar,
     $ ureac,wodr,wodrt,worr,worrt,xi1,xistsv)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print a table of the elemental composition of the aqueous
c     solution.
c
      call prteca(cteaq,mrmlra,nct,nctmax,noutpt,ppmwe,
     $ qrho,rho,uelem)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute and print a table of the numerical compostion of the
c     aqueous solution.
c
      call prtnca(ctb,jflag,jsflag,mrmlra,mwtsp,narn1,narn2,
     $ nbasp,nbaspd,nbt,nbtmax,noutpt,nstmax,qrho,rho,uspec,wfh2o)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute and print a table of the sensible compostion of the
c     aqueous solution.
c
      call prtsca(ctb,jflag,jsflag,mrmlra,mwtsp,narn1,narn2,
     $ nbasp,nbaspd,nbt,nbtmax,nelect,nhydr,nhydx,noutpt,no2gaq,
     $ nstmax,qrho,rho,uspec,wfh2o)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print pH, Eh, and pe-, all with reference to appropriate
c     pH scales. Also print the pHCl.
c
      call prpheh(ah,ahmes,ahnbs,eh,ehmes,ehnbs,iopg,
     $ nopgmx,noutpt,pch,pe,pemes,penbs,ph,phcl,phmes,phnbs,
     $ qphcl,qredox,qrho)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print various aqueous phase parameters.
c
      call prtvpa(abar,acfw,acfwlg,actw,actwlg,a3bar,fje,
     $ fjest,fo2,fo2lg,fxi,fxist,iopg,mlmrra,mrmlra,nopgmx,noutpt,
     $ osc,oscst,qrho,rhoc,rhowc,sigmam,sigmst,tdsglw,tdspkc,
     $ tdsplc,vosol,wfh2o,wftds,woh2o,wosol,wotds,xbarw,xbrwlg)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print more precise mass results for the aqueous phase.
c
      write (ux16a,'(1pg15.8)') woh2o
      call lejust(ux16a)
      j1 = ilnobl(ux16a)
      write (ux16b,'(1pg15.8)') wotds
      call lejust(ux16b)
      j2 = ilnobl(ux16b)
      write (ux16c,'(1pg15.8)') wosol
      call lejust(ux16c)
      j3 = ilnobl(ux16c)
      write (noutpt,1500) ux16a(1:j1),ux16b(1:j2),ux16c(1:j3)
 1500 format(//11x,'--- More Precise Aqueous Phase Masses ---',
     $ //23x,'Solvent mass= ',a,' g',
     $ /17x,'Solutes (TDS) mass= ',a,' g',
     $ /14x,'Aqueous solution mass= ',a,' g',//)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print a table of computed alkalinity parameters.
c
      call prtalk(alki,alk1,alk2,mrmlra,noutpt,ntf1t,ntf2t,
     $ qrho,rho,tempc,wfh2o)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute and print a table describing the aqueous phase charge
c     balance.
c
      call gszm(conc,jcsort,narn1,narn2,nstmax,sigza,sigzc,
     $ sigzi,sigzm,zchar)
c
      sigzia = sigzi/wconst
      sigzma = sigzm/wconst
c
      chdrft = sigzia - electr
      do nrc = 1,nrct
        if (jcode(nrc) .eq. 2) then
          nsr = nrndex(nrc)
          chdrft = chdrft - elecsr(nsr)*xirct(nrc)
        endif
      enddo
c
      betchb = chdrft/sigzma
      sigzis = wfh2o*sigzi
      elects = wfh2o*electr
      chdrfs = wfh2o*chdrft
      sigzms = wfh2o*sigzma
c
      write (noutpt,1510) sigzia,electr,chdrft,sigzma,sigzis,elects,
     $ chdrfs,sigzms,betchb
 1510 format(//7x,' --- Aqueous Solution Charge Balance ---',
     $ //6x,'    Actual Charge imbalance= ',1pe11.4,' eq',
     $  /6x,'  Expected Charge imbalance= ',e11.4,' eq',
     $  /6x,'         Charge discrepancy= ',e11.4,' eq',
     $  /6x,'        Sigma |equivalents|= ',e11.4,' eq',
     $ //6x,'    Actual Charge imbalance= ',e11.4,' eq/kg.solu',
     $  /6x,'  Expected Charge imbalance= ',e11.4,' eq/kg.solu',
     $  /6x,'         Charge discrepancy= ',e11.4,' eq/kg.solu',
     $  /6x,'        Sigma |equivalents|= ',e11.4,' eq/kg.solu',
     $ //6x,'Relative charge discrepancy= ',e11.4)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print detailed listing of aqueous species.
c
      call prtaqs(acflg,actlg,conc,conclg,iopr,jcsort,narn1,
     $ narn2,noprmx,noutpt,nstmax,uspec)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print the contributions of species to each aqueous mass total.
c
      call prtpct(conc,csts,ctb,iopr,jcsort,jflag,narn1,narn2,
     $ nbasp,nbaspd,nbt,nbtmax,nelect,nhydr,noprmx,no2gaq,noutpt,
     $ nstmax,nsts,nstsmx,nstsr,uspec)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopr(9) .gt. 0) then
c
c       Compute and print the mean ionic activities and activity
c       coefficients.
c
        call prtmip(acflg,actlg,conclg,ctb,nbaspd,nbt,nbtmax,
     $  nelect,nhydr,nhydx,noutpt,nstmax,uspec,zchar)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopr(5) .gt. 0) then
c
c       Compute and print activity ratios of aqueous species.
c
        call prtacr(actlg,iopr,jsflag,nbaspd,nbt,nbtmax,nelect,
     $  nhydr,noprmx,no2gaq,noutpt,nstmax,uspec,zchar)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qredox) then
c
c       Print data for redox reactions that are not constrained to be
c       at equilibrium.
c
        call prtrdx(ah,ahrc,cdrsd,eh,ehrc,fo2lg,fo2lrc,jflag,
     $  jsflag,narn1,nbasp,nbaspd,nbt,nbtmax,ndrsd,ndrsmx,ndrsrd,
     $  nelect,nhydr,no2gaq,noutpt,nstmax,pe,perc,uspec)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(1).eq.0 .or. iopt(1).eq.1) then
c
c       List the solid phases in the equilibrium system (ES).
c
        write (noutpt,1600)
 1600   format(//21x,'--- Summary of Solid Phases (ES) ---')
        write (noutpt,1610)
 1610   format(/'   Phase/End-member',12x,'Log moles',5x,'Moles',8x,
     $  'Grams',5x,'Volume, cm3',/)
c
        ncount = 0
c
        do kcol = km1,kmt
          ncount = ncount + 1
          np = ipndx1(kcol)
          write (noutpt,1620) uphase(np),loph(np),moph(np),woph(np),
     $    voph(np)
 1620     format(1x,a24,4x,f11.4,3(2x,1pe11.4))
        enddo
c
        nplast = 0
        do kcol = kx1,kxt
          np = ipndx1(kcol)
          if (np .ne. nplast) then
            if (ncount .gt. 0) write (noutpt,1630)
 1630       format(1x)
            ncount = ncount + 1
            write (noutpt,1620) uphase(np),loph(np),moph(np),woph(np),
     $      voph(np)
c
            nr1 = ncmpr(1,np)
            nr2 = ncmpr(2,np)
            nt = nr2 - nr1 + 1
            if (nt .ge. 2) then
              do ns = nr1,nr2
                if (jsflag(ns) .le. 0) then
                  write (noutpt,1640) uspec(ns),losp(ns),mosp(ns),
     $            wosp(ns),vosp(ns)
 1640             format(3x,a24,2x,f11.4,3(2x,1pe11.4))
                endif
              enddo
            endif
            nplast = np
          endif
        enddo
c
        nplast = 0
        do np = iern1,iern2
          if (moph(np) .gt. 0.) then
            if (np .ne. nplast) then
              if (ncount .gt. 0) write (noutpt,1630)
              ncount = ncount + 1
              write (noutpt,1620) uphase(np),loph(np),moph(np),
     $        woph(np),voph(np)
c
              ne = np - iern1 + 1
              do je = 1,jgext(ne)
                if (je .gt. 1) write (noutpt,1642)
 1642           format(1x)
                ns = jern1(je,ne) - 1
                do ie = 1,ngext(je,ne)
                  ns = ns + 1
                  if (jsflag(ns) .le. 0) then
                    write (noutpt,1640) uspec(ns),losp(ns),mosp(ns),
     $              wosp(ns),vosp(ns)
                  endif
                enddo
              enddo
              nplast = np
            endif
          endif
        enddo
c
        if (ncount .le. 0) write (noutpt,1650)
 1650   format(1x,'None')
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  200 continue
c
      if (iopt(1) .eq. 2) then
c
c       List the mineral assemblage which is forming instantaneously
c       with respect to Xi (on a per unit Xi basis). When this is
c       calculated at a phase boundary, the assemblage given is the one
c       to the left of the boundary. When this is calculated just
c       after the boundary, it is the one to the right of the boundary.
c
        if (npts .gt. 0) then
          write (noutpt,1700)
 1700     format(//3x,'--- Summary of Solid Phases (ES, Instantaneous',
     $    ' Basis) ---',/19x,'(defined by derivatives)')
          write (noutpt,1710)
 1710     format(//3x,'Phase/End-member',11x,'d mol/d Xi',2x,
     $    'd grams/d Xi',2x,'d V(cm3)/d Xi',/)
c
          ncount = 0
c
          do npe = 2,npet0
            np = iemop0(npe)
            ncount = ncount + 1
            if (np.ge.imrn1 .and. np.le.imrn2) then
c
c             Have a pure mineral.
c
              write (noutpt,1720) uphase(np),mophj(np),wophj(np),
     $        vophj(np)
 1720         format(1x,a24,2x,2(2x,1pe11.4),2x,1pe12.5,2x,1pe11.4)
            endif
          enddo
c
          do npe = 2,npet0
            np = iemop0(npe)
            ncount = ncount + 1
            if (np.ge.ixrn1 .and. np.le.ixrn2) then
c
c             Have a solid solution.
c
              write (noutpt,1730) uphase(np),mophj(np),wophj(np),
     $        vophj(np)
 1730         format(/1x,a24,2x,2(2x,1pe11.4),2x,1pe12.5,2x,1pe11.4)
              nr1 = ncmpe0(1,npe)
              nr2 = ncmpe0(2,npe)
              do nse = nr1,nr2
                ns = iemos0(nse)
                write (noutpt,1740) uspec(ns),mospj(ns),wospj(ns),
     $          vospj(ns)
 1740           format(3x,a24,2(2x,1pe11.4),2x,1pe12.5,2x,1pe11.4)
              enddo
            endif
          enddo
c
          do npe = 2,npet0
            np = iemop0(npe)
            ncount = ncount + 1
            if (np.ge.iern1 .and. np.le.iern2) then
c
c             Have a generic ion exchange phase.
c
              write (noutpt,1730) uphase(np),mophj(np),wophj(np),
     $        vophj(np)
              nr1 = ncmpe0(1,npe)
              nr2 = ncmpe0(2,npe)
              do nse = nr1,nr2
                ns = iemos0(nse)
                write (noutpt,1740) uspec(ns),mospj(ns),wospj(ns),
     $          vospj(ns)
              enddo
            endif
          enddo
c
          if (ncount .le. 0) write (noutpt,1650)
c
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(1) .eq. 2) then
c
c       Print a summary of the solid phases in the ES and the PRS.
c
        write (noutpt,1800)
 1800   format(//12x,'--- Summary of Solid Phases (ES + PRS) ---')
        write (noutpt,1610)
c
        ncount = 0
        do np = 1,npt
          if (np .ne. iaqsln) then
            if (jpflag(np) .lt. 2) then
              if (mopht(np) .gt. 0.) then
                nr1 = ncmpr(1,np)
                nr2 = ncmpr(2,np)
                nt = nr2 - nr1 + 1
                mxx = mopht(np)
                lxx = tlg(mxx)
c
                if (nt .lt. 2) then
                  ncount = ncount + 1
                  write (noutpt,1620) uphase(np),lxx,mxx,wopht(np),
     $            vopht(np)
                else
                  if (ncount .gt. 0) write (noutpt,1630)
                  ncount = ncount + 1
                  write (noutpt,1620) uphase(np),lxx,mxx,
     $            wopht(np),vopht(np)
                  do ns = nr1,nr2
                    if (mospt(ns) .gt. 0.) then
                      mxx = mospt(ns)
                      lxx = tlg(mxx)
                      write (noutpt,1640) uspec(ns),lxx,mxx,wospt(ns),
     $                vospt(ns)
                    endif
                  enddo
                endif
c
              endif
            endif
          endif
        enddo
c
        if (ncount .le. 0) write (noutpt,1650)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if ((nmrt + nxrt + nert) .gt. 0) then
c
        write (noutpt,1810)
 1810   format(//12x,'--- Grand Summary of Solid Phases (ES + PRS',
     $  ' + Reactants) ---')
        write (noutpt,1610)
c
        ncount = 0
        do np = imrn1,imrn2
          qrct = .false.
          if (nmrt .gt. 0) then
            do nrc = 1,nrct
              if (jcode(nrc) .eq. 0) then
                np1 = nrndex(nrc)
                if (np1 .eq. np) then
                  qrct = .true.
                  go to 300
                endif
              endif
            enddo
          endif
c
  300     if (mophg(np).gt.0. .or. qrct) then
            ncount = ncount + 1
            ns = ncmpr(1,np)
            mxx = mophg(np)
            lxx = tlg(abs(mxx))
            write (noutpt,1620) uphase(np),lxx,mxx,wophg(np),
     $      vophg(np)
          endif
        enddo
c
        do np = ifrn1,ifrn2
          qrct = .false.
          if (nmrt .gt. 0) then
            do nrc = 1,nrct
              if (jcode(nrc) .eq. 0) then
                np1 = nrndex(nrc)
                if (np1 .eq. np) then
                  qrct = .true.
                  go to 310
                endif
              endif
            enddo
          endif
c
  310     if (mophg(np).gt.0. .or. qrct) then
            ncount = ncount + 1
            ns = ncmpr(1,np)
            mxx = mophg(np)
            lxx = tlg(abs(mxx))
            write (noutpt,1620) uphase(np),lxx,mxx,wophg(np),
     $      vophg(np)
          endif
        enddo
c
        if (iopt(4) .gt. 0) then
          do np = ixrn1,ixrn2
            if (jpflag(np) .lt. 2) then
              qrct = .false.
c
              if (nxrt .gt. 0) then
                do nrc = 1,nrct
                  if (jcode(nrc) .eq. 1) then
                    np1 = nrndex(nrc)
                    if (np1 .eq. np) then
                      qrct = .true.
                      go to 320
                    endif
                  endif
                enddo
              endif
c
  320         if (mophg(np).gt.0. .or. qrct) then
                if (ncount .gt. 0) write (noutpt,1630)
                ncount = ncount + 1
                mxx = mophg(np)
                lxx = tlg(abs(mxx))
                write (noutpt,1620) uphase(np),lxx,mxx,wophg(np),
     $          vophg(np)
                nr1 = ncmpr(1,np)
                nr2 = ncmpr(2,np)
                do ns = nr1,nr2
                  mxx = mospg(ns)
                  lxx = tlg(abs(mxx))
                  write (noutpt,1640) uspec(ns),lxx,mxx,
     $            wospg(ns),vospg(ns)
                enddo
              endif
c
            endif
          enddo
        endif
c
        do np = iern1,iern2
          if (jpflag(np) .lt. 2) then
            qrct = .false.
c
            if (nert .gt. 0) then
              do nrc = 1,nrct
                if (jcode(nrc) .eq. 5) then
                  np1 = nrndex(nrc)
                  if (np1 .eq. np) then
                    qrct = .true.
                    go to 330
                  endif
                endif
              enddo
            endif
c
  330       if (mophg(np).gt.0. .or. qrct) then
              if (ncount .gt. 0) write (noutpt,1630)
              ncount = ncount + 1
              mxx = mophg(np)
              lxx = tlg(abs(mxx))
              write (noutpt,1620) uphase(np),lxx,mxx,wophg(np),
     $        vophg(np)
              nr1 = ncmpr(1,np)
              nr2 = ncmpr(2,np)
              do ns = nr1,nr2
                mxx = mospg(ns)
                if (mxx .gt. 0.) then
                  lxx = tlg(abs(mxx))
                  write (noutpt,1640) uspec(ns),lxx,mxx,
     $            wospg(ns),vospg(ns)
                endif
              enddo
            endif
c
          endif
        enddo
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1820) wosoct,vosoct,wodrt,vodrt,dwoso,dvoso
 1820 format(//26x,'Mass, grams',7x,'Volume, cm3',
     $ //11x,'Created  ',6x,1pe12.5,6x,e12.5,
     $ /11x,'Destroyed',6x,e12.5,6x,e12.5,
     $ /11x,'Net      ',6x,e12.5,6x,e12.5,
     $ //10x,'These volume totals may be incomplete because of missing',
     $ /10x,'partial molar volume data in the data base.',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qftpr2) then
        write (noutpt,1000)
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print saturation index and affinity tables for reactions in
c     the aqueous phase not constrained to be at equilibrium. The
c     data correspond to the reactions and thermodynamic data in
c     the 'd' set.
c
      call prtsia(affsd,jflagd,jflgi,jsflag,narn1,narn2,nbasp,
     $ nbaspd,nbt,nbtmax,ndrsd,ndrsmx,ndrsrd,nhydr,noutpt,nrdxsp,
     $ nstmax,sidrsp,uspec)
c
      if (iopr(7) .ge. 0) then
c
c       Print saturation index and affinity tables for the
c       various non-aqueous phases.
c
        call prtsat(affpd,iern1,iern2,ilrn1,ilrn2,imrn1,imrn2,
     $  iopr,iopt,ixrn1,ixrn2,jpflag,noutpt,noprmx,noptmx,nptmax,
     $  sidrph,tolsat,uphase)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(4) .gt. 0) then
c
c       Print tables of the compositions and affinities of solid
c       solution phases present in the equilbrium system (ES).
c
        write (noutpt,1900)
 1900   format(//11x,'--- Solid Solution Product Phases ---',/)
c
        ncount = 0
        do np = ixrn1,ixrn2
          if (jpflag(np) .eq. -1) then
            ncount = ncount + 1
            call prtsso(acflg,actlg,affpd,affsd,ixrn1,ixrn2,
     $      jsol,jsomax,ncmpr,noutpt,np,nptmax,nstmax,nxtmax,sidrph,
     $      sidrsp,tolsat,uspec,uphase,uxtype,xbar,xbarlg)
          endif
        enddo
        if (ncount .le. 0) write (noutpt,1910)
 1910   format(1x,'None')
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (net .gt. 0) then
c
c       Print tables of the compositions and affinities of generic
c       ion exchange phases present in the equilbrium system (ES).
c
        write (noutpt,1930)
 1930   format(//11x,'--- Generic Ion Exchange Phases ---',/)
c
        tolspf = tolsat
        ncount = 0
        do np = iern1,iern2
          ncount = ncount + 1
c
          call prtgex(acflg,actlg,affpd,affsd,cegexs,conc,
     $    egexjc,egexjf,egexpa,egexpc,egexs,egexw,iern1,iern2,ietmax,
     $    jern1,jern2,jetmax,jgext,kern1,kern2,ketmax,kgexsa,moph,
     $    mosp,netmax,ngexsa,ngext,noutpt,np,nptmax,nstmax,sidrph,
     $    sidrsp,tolspf,ugexj,ugexmo,uspec,uphase,xbar,xbarlg,xgexw,
     $    wkgwi)
c
        enddo
        if (ncount .le. 0) write (noutpt,1940)
 1940   format(1x,'None')
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopr(8) .ge. 0) then
c
c       Print table of equilibrium fugacities.
c
        call prtfug(jgsort,fugac,fugalg,jsflag,ngrn1,ngt,ngtmax,
     $  noutpt,nstmax,uspec)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,2000)
 2000 format(/' - - - - - - - - - - - - - - - - - - - - - - - - -',
     $ ' - - - - - - -',//)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
