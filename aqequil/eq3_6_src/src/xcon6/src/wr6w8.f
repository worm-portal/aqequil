      subroutine wr6w8(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,
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
c     This subroutine writes the EQ6 INPUT file in compact ("W") format
c     for version 8.0.
c
c     This subroutine is a near-clone of EQLIB/wr6pkw.f.
c
c     The calling sequence of this subroutine is identical to that of
c     XCON6/wr6d8.f, EQLIB/wr6pkw.f, and EQLIB/wr6pkd.f.
c
c     The calling sequence of this subroutine is identical to that of
c     EQ6/rd6inw.f, EQ6/rd6ind.f, XCON6/rd6w8.f, and XCON6/rd6d8.f,
c     except that newin is added and ninpts, nprob, noutpt, qend,
c     and qrderr are deleted.
c
c     This subroutine is called by:
c
c       XCON6/xcon6.f
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
      real*8 awmaxi,awmini,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,
     $ dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,
     $ dlxplo,dlxprl,dlxprn,ehmaxi,ehmini,electr,o2maxi,o2mini,phmaxi,
     $ phmini,pressb,pressi,tempcb,tempci,timmxi,tistti,tolbt,toldl,
     $ tolsat,tolxsf,ximaxi,xistti
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,iei,iki,je,jei,j2,j3,j4,kcol,krow,n,nbi,nci,ne,ner,npi,
     $ nrc,nsi,nsr,nxr
c
      integer ilnobl
c
      logical qx
c
      character*8 uendit
c
c-----------------------------------------------------------------------
c
      data uendit /'endit.  '/
c
c-----------------------------------------------------------------------
c
c     Main title.
c
      do n = 1,ntitl1
        j2 = ilnobl(utitl1(n))
        write (newin,1020) utitl1(n)(1:j2)
 1020   format(a)
      enddo
      j3 = ilnobl(uendit)
      if (ntitl1 .lt. ntitmx) write (newin,1020) uendit(1:j3)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Temperature parameters.
c     Note: ttk(1) = ttk1, etc.
c
      write (newin,1130) jtemp,tempcb,ttk(1),ttk(2)
 1130 format(5x,'jtemp= ',i2,/4x,'tempcb= ',1pe12.5,
     $ /6x,'ttk1= ',1pe12.5,6x,'ttk2= ',e12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Pressure parameters.
c     Note: ptk(1) = ptk1, etc.
c
      write (newin,1140) jpress,pressb,ptk(1),ptk(2)
 1140 format(4x,'jpress= ',i2,/4x,'pressb= ',1pe12.5,
     $ /6x,'ptk1= ',1pe12.5,6x,'ptk2= ',e12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Number of reactants.
c
      write (newin,1490) nrct
 1490 format(6x,'nrct= ',i2)
c
      write (newin,1500)
 1500 format('*-------------------------------------------------',
     $ '----------------------------')
c
c     Reactants.
c
      nxr = 0
      nsr = 0
      ner = 0
c
      do nrc = 1,nrct
c
c       Name, flags, and masses.
c
        j2 = ilnobl(ureac(nrc))
        write (newin,1540) ureac(nrc)(1:j2),jcode(nrc),jreac(nrc),
     $  morr(nrc),modr(nrc)
 1540   format(2x,'reactant= ',a,/5x,'jcode= ',i2,15x,'jreac= ',i2,
     $  /6x,'morr= ',1pe12.5,6x,'modr= ',e12.5)
c
        if (jcode(nrc) .eq. 1) then
c
c         Solid solution compositions.
c
          nxr = nxr + 1
          do iki = 1,ixrti(nxr)
            write (newin,1580) ucxri(iki,nxr),rxbari(iki,nxr)
 1580       format(3x,a24,3x,1pe12.5)
          enddo
          j3 = ilnobl(uendit)
          write (newin,1570) uendit(1:j3)
 1570     format(3x,a)
c
        elseif (jcode(nrc) .eq. 2) then
c
c         Special reactant compositions.
c
          nsr = nsr + 1
          write (newin,1620) vreac(nrc)
 1620     format(5x,'vreac= ',e12.5)
c
          write (newin,1640)
 1640     format ('* Elemental composition')
          do nci = 1,iesrti(nsr)
            write (newin,1650) uesri(nci,nsr),cesri(nci,nsr)
 1650       format(3x,a8,3x,1pe22.15)
          enddo
          j3 = ilnobl(uendit)
          write (newin,1570) uendit(1:j3)
c
          write (newin,1680)
 1680     format ('* Reaction')
          do nbi = 1,ibsrti(nsr)
            write (newin,1690) ubsri(nbi,nsr)(1:24),cbsri(nbi,nsr)
 1690       format(3x,a24,3x,1pe22.15)
          enddo
          write (newin,1570) uendit(1:j3)
c
        elseif (jcode(nrc) .eq. 5) then
c
c         Generic ion exchanger phase compositions.
c
          ner = ner + 1
c
          j3 = ilnobl(ugermo(ner))
          write (newin,1692) ugermo(ner)(1:j3)
 1692     format(4x,'ugermo= ',a)
c
          j4 = ilnobl(uendit)
          do jei = 1,jgerti(ner)
            j3 = ilnobl(ugerji(jei,ner))
            write (newin,1700) ugerji(jei,ner)(1:j3)
 1700       format(3x,a)
            do iei = 1,igerti(jei,ner)
              write (newin,1710) ugersi(iei,jei,ner),egersi(iei,jei,ner)
 1710         format(6x,a24,3x,1pe12.5)
            enddo
            write (newin,1720) uendit(1:j4)
 1720       format(6x,a)
          enddo
          write (newin,1570) uendit(1:j4)
c
        endif
c
c       Surface area parameters.
c
        write (newin,1750) nsk(nrc),sfcar(nrc),ssfcar(nrc),fkrc(nrc)
 1750   format(7x,'nsk= ',i2,15x,'sfcar= ',1pe12.5,4x,'ssfcar= ',e12.5,
     $  /6x,'fkrc= ',e12.5)
c
c       Kinetic rate laws.
c
        write (newin,1770) nrk(1,nrc),nrk(2,nrc)
 1770   format(6x,'nrk1= ',i2,16x,'nrk2= ',i2)
c
c       Rate law parameters, forward direction (destruction).
c
        if (nrk(1,nrc) .eq. 1) then
c
c         Arbitrary kinetics.
c
          write (newin,1800) (rkb(i,1,nrc), i = 1,3)
 1800     format(6x,'rkb1= ',1pe12.5,6x,'rkb2= ',e12.5,
     $    6x,'rkb3= ',e12.5)
c
        elseif (nrk(1,nrc) .eq. 2) then
c
c         Transition state theory.
c
          write (newin,1820) imech(1,nrc)
 1820     format(5x,'imech= ',i2)
          do i = 1,imech(1,nrc)
            write (newin,1840) rkb(i,1,nrc),trkb(i,1,nrc),
     $      iact(i,1,nrc)
 1840       format(7x,'rkb= ',1pe12.5,6x,'trkb= ',1pe12.5,6x,
     $      'iact= ',i2)
c
            write (newin,1860) eact(i,1,nrc),hact(i,1,nrc)
 1860       format(6x,'eact= ',1pe12.5,6x,'hact= ',e12.5)
c
            write (newin,1880) ndact(i,1,nrc),csigma(i,1,nrc)
 1880       format(5x,'ndact= ',i2,14x,'csigma= ',1pe12.5)
c
            do n = 1,ndact(i,1,nrc)
              write (newin,1920) udac(n,i,1,nrc),cdac(n,i,1,nrc)
 1920         format(6x,'udac= ',a24,6x,'cdac= ',1pe12.5)
            enddo
          enddo
c
        elseif (nrk(1,nrc) .eq. 3) then
c
c         Linear rate law.
c
          i = imech(1,nrc)
          write (newin,1840) rkb(i,1,nrc),trkb(i,1,nrc),
     $    iact(i,1,nrc)
          write (newin,1860) eact(i,1,nrc),hact(i,1,nrc)
c
        endif
c
c       Rate law parameters, backward direction (formation).
c
        if (nrk(2,nrc) .eq. 1) then
c
c         Arbitrary kinetics.
c
          write (newin,1800) (rkb(i,2,nrc), i = 1,3)
c
        elseif (nrk(2,nrc) .eq. 2) then
c
c         Transition state theory.
c
          write (newin,1820) imech(2,nrc)
          do i = 1,imech(2,nrc)
            write (newin,1840) rkb(i,2,nrc),trkb(i,2,nrc),
     $      iact(i,2,nrc)
c
            write (newin,1860) eact(i,2,nrc),hact(i,2,nrc)
c
            write (newin,1880) ndact(i,2,nrc),csigma(i,2,nrc)
c
            do n = 1,ndact(i,2,nrc)
               write (newin,1920) udac(n,i,2,nrc),cdac(n,i,2,nrc)
             enddo
          enddo
c
        elseif (nrk(2,nrc) .eq. 3) then
c
c         Linear rate law.
c
          i = imech(2,nrc)
          write (newin,1840) rkb(i,2,nrc),trkb(i,2,nrc),
     $    iact(i,1,nrc)
          write (newin,1860) eact(i,2,nrc),hact(i,2,nrc)
c
        endif
c
        write (newin,1500)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Starting, minimum, and maximum values of key run parameters.
c
      write (newin,1150) xistti,ximaxi,tistti,timmxi
 1150 format(4x,'xistti= ',1pe12.5,4x,'ximaxi= ',e12.5,
     $ /4x,'tistti= ',1pe12.5,4x,'timmxi= ',e12.5)
c
      write (newin,1160) phmini,phmaxi,ehmini,ehmaxi,o2mini,o2maxi,
     $ awmini,awmaxi,kstpmx
 1160 format(4x,'phmini= ',1pe12.5,4x,'phmaxi= ',e12.5,
     $ /4x,'ehmini= ',1pe12.5,4x,'ehmaxi= ',e12.5,
     $ /4x,'o2mini= ',1pe12.5,4x,'o2maxi= ',e12.5,
     $ /4x,'awmini= ',1pe12.5,4x,'awmaxi= ',e12.5,
     $ /4x,'kstpmx= ',i12)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print interval parameters.
c
      write (newin,1170) dlxprn,dlxprl,dltprn,dltprl
 1170 format(4x,'dlxprn= ',1pe12.5,4x,'dlxprl= ',e12.5,
     $ /4x,'dltprn= ',e12.5,4x,'dltprl= ',e12.5)
c
      write (newin,1180) dlhprn,dleprn,dloprn,dlaprn,ksppmx
 1180 format(4x,'dlhprn= ',1pe12.5,4x,'dleprn= ',e12.5,
     $ /4x,'dloprn= ',e12.5,4x,'dlaprn= ',e12.5,/4x,'ksppmx= ',i12)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Plot interval parameters.
c
      write (newin,1190) dlxplo,dlxpll,dltplo,dltpll
 1190 format(4x,'dlxplo= ',1pe12.5,4x,'dlxpll= ',e12.5,
     $ /4x,'dltplo= ',e12.5,4x,'dltpll= ',e12.5)
c
      write (newin,1195) dlhplo,dleplo,dloplo,dlaplo,ksplmx
 1195 format(4x,'dlhplo= ',1pe12.5,4x,'dleplo= ',e12.5,
     $ /4x,'dloplo= ',e12.5,4x,'dlaplo= ',e12.5,/4x,'ksplmx= ',i12)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Comment header for iopt, iopr, iodb, and iopg option switches.
c
      write (newin,1200)
 1200 format('*',15x,'1    2    3    4    5    6    7    8    9   10')
c
c     Iopt option switches.
c     Note: iopt(1) = iopt1, etc.
c
      write (newin,1220) (iopt(i), i = 1,20)
 1220 format(2x,'iopt1-10= ',10i5,/1x,'iopt11-20= ',10i5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iopr option switches.
c     Note: iopr(1) = iopr1, etc.
c
      write (newin,1240) (iopr(i), i = 1,20)
 1240 format(2x,'iopr1-10= ',10i5,/1x,'iopr11-20= ',10i5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iodb option switches.
c     Note: iodb(1) = iodb1, etc.
c
      write (newin,1260) (iodb(i), i = 1,20)
 1260 format(2x,'iodb1-10= ',10i5,/1x,'iodb11-20= ',10i5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Number of nxopt options.
c
      write (newin,1310) nxopt
 1310 format(5x,'nxopt= ',i2)
c
c     Nxopt options.
c
      if (nxopt .gt. 0) then
        do n = 1,nxopt
          j2 = ilnobl(uxcat(n))
          write (newin,1340) uxopt(n),uxcat(n)(1:j2)
 1340     format(4x,'option= ',a6,1x,a)
        enddo
        write (newin,1360) nxopex
 1360   format(4x,'nxopex= ',i2)
        if (nxopex .gt. 0) then
          do n = 1,nxopex
            j2 = ilnobl(uxopex(n))
            write (newin,1390) uxopex(n)(1:j2)
 1390       format(' exception= ',a)
          enddo
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Number of nffg options.
c
      write (newin,1400) nffg
 1400 format(6x,'nffg= ',i2)
c
c     Nffg options.
c
      if (nffg .gt. 0) then
        do n = 1,nffg
          j2 = ilnobl(uffg(n))
          write (newin,1440) uffg(n)(1:j2),moffg(n),xlkffg(n)
 1440     format (3x,'species= ',a,/5x,'moffg= ',1pe12.5,4x,
     $    'xlkffg= ',e12.5)
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum finite-difference order.
c
      write (newin,2020) nordmx
 2020 format(4x,'nordmx= ',i3)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Newton-Raphson convergence tolerances.
c
      write (newin,2030) tolbt,toldl
 2030 format(5x,'tolbt= ',1pe12.5,5x,'toldl= ',e12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum number of Newton-Raphson iterations.
c
      write (newin,2110) itermx
 2110 format(4x,'itermx= ',i3)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Search/find tolerance.
c
      write (newin,2032) tolxsf
 2032 format(4x,'tolxsf= ',1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Saturation tolerance.
c
      write (newin,2034) tolsat
 2034 format(4x,'tolsat= ',1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum number of phase assemblage tries.
c
      write (newin,2112) ntrymx
 2112 format(4x,'ntrymx= ',i3)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Zero-order step size (in Xi).
c
      write (newin,2090) dlxmx0
 2090 format(4x,'dlxmx0= ',1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     PRS transfer interval in Xi.
c
      write (newin,2010) dlxdmp
 2010 format(4x,'dlxdmp= ',1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (newin,1500)
c
c     Write the bottom half of the PICKUP file.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Old title.
c
      do n = 1,ntitl2
        j2 = ilnobl(utitl2(n))
        write (newin,1020) utitl2(n)(1:j2)
      enddo
      j3 = ilnobl(uendit)
      if (ntitl1 .lt. ntitmx) write (newin,1020) uendit(1:j3)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Special basis switches.
c
      write (newin,2610)
 2610 format('*   Special basis switches')
      write (newin,2630) nsbswt
 2630 format(4x,'nsbswt= ',i3)
c
      do n = 1,nsbswt
        j2 = ilnobl(usbsw(1,n))
        write (newin,2650) usbsw(1,n)(1:j2)
 2650   format('species= ',a)
        j2 = ilnobl(usbsw(2,n))
        write (newin,2670) usbsw(2,n)(1:j2)
 2670   format(2x,'switch with= ',a)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Original temperature.
c
      write (newin,2690) tempci
 2690 format(4x,'tempci= ',1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Original pressure.
c
      write (newin,2710) pressi
 2710 format(4x,'pressi= ',1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Ion exchanger creation.
c
      write (newin,2790)
 2790 format('* Ion exchanger creation')
c
      write (newin,2795) qgexsh
 2795 format(4x,'qgexsh= ',l8)
c
      write (newin,2810) net
 2810 format(7x,'net= ',i3)
c
      do ne = 1,net
        j2 = ilnobl(ugexp(ne))
        write (newin,2840) ugexp(ne)(1:j2)
 2840   format(5x,'ugexp= ',a)
c
        write (newin,2860) mwtges(ne)
 2860   format(4x,'mwtges= ',1pe12.5)
c
        j3 = ilnobl(ugexmo(ne))
        write (newin,2870) ugexmo(ne)(1:j3)
 2870   format(4x,'ugexmo= ',a)
c
        write (newin,2880) tgexp(ne)
 2880   format(5x,'tgexp= ',1pe12.5)
c
        write (newin,2900) jgext(ne)
 2900   format(5x,'jgext= ',i3)
c
        do je = 1,jgext(ne)
          j3 = ilnobl(ugexj(je,ne))
          write (newin,2940) ugexj(je,ne)(1:j3)
 2940     format(5x,'ugexj= ',a)
c
          write (newin,2960) cgexj(je,ne),zgexj(je,ne)
 2960     format(5x,'cgexj= ',1pe12.5,5x,'zgexj= ',e12.5)
c
          write (newin,3010) ngexrt(je,ne)
 3010     format(4x,'ngexrt= ',i3)
c
          do n = 1,ngexrt(je,ne)
            j2 = ilnobl(ugexr(n,je,ne))
            write (newin,3040) ugexr(n,je,ne)(1:j2)
 3040       format(5x,'ugexr= ',a)
            j4 = ilnobl(uxkgex(n,je,ne))
            write (newin,3060) xlkgex(n,je,ne),uxkgex(n,je,ne)(1:j4)
 3060       format(4x,'xlkgex= ',1pe12.5,5x,'units= ',a)
            j4 = ilnobl(uhfgex(n,je,ne))
            write (newin,3070) xhfgex(n,je,ne),uhfgex(n,je,ne)(1:j4)
 3070       format(4x,'xhfgex= ',1pe12.5,5x,'units= ',a)
            j4 = ilnobl(uvfgex(n,je,ne))
            write (newin,3080) xvfgex(n,je,ne),uvfgex(n,je,ne)(1:j4)
 3080       format(4x,'xvfgex= ',1pe12.5,5x,'units= ',a)
          enddo
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Number of nxmod options.
c
      write (newin,3130) nxmod
 3130 format(5x,'nxmod= ',i2)
c
c     Nxmod options.
c
      do n = 1,nxmod
        j2 = ilnobl(uxmod(n))
        write (newin,3160) uxmod(n)(1:j2),kxmod(n),xlkmod(n)
 3160   format(3x,'species= ',a,
     $  /4x,'option= ',i2,14x,'xlkmod= ',1pe12.5)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iopg options.
c     Note: iopg(1) = iopg1, etc.
c
      write (newin,1200)
      write (newin,3180) (iopg(i), i = 1,20)
 3180 format(2x,'iopg1-10= ',10i5,/1x,'iopg11-20= ',10i5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Index limits.
c
      write (newin,3230) kct,kbt,kmt,kxt,kdim,kprs
 3230 format(7x,'kct= ',i2,17x,'kbt= ',i2,17x,'kmt= ',i2,/
     $ 7x,'kxt= ',i2,16x,'kdim= ',i2,16x,'kprs= ',i2)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Species for which mass balances are defined.
c
      write (newin,3330)
 3330 format('* Data file basis species and jflgi values')
c
      do nbi = 1,nbti
        write (newin,3350) ubmtbi(nbi),jflgi(nbi)
 3350   format(3x,a48,3x,i2)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Mass balance totals.
c
      write (newin,3400)
 3400 format('*',7x,
     $ 'Mass balance totals     Aqueous mass balance totals')
c
      do nbi = 1,nbti
        write (newin,3420) mtbi(nbi),mtbaqi(nbi)
 3420   format(6x,1pe22.15,6x,e22.15)
      enddo
c
      write (newin,3440) electr
 3440 format(9x,'Electrical imbalance= ',3x,1pe22.15)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Ordinary basis switches.
c
      write (newin,3450)
 3450 format('*   Ordinary basis switches')
      write (newin,3460) nobswt
 3460 format(4x,'nobswt= ',i3)
c
      do n = 1,nobswt
        j2 = ilnobl(uobsw(1,n))
        write (newin,2650) uobsw(1,n)(1:j2)
        j2 = ilnobl(uobsw(2,n))
        write (newin,2670) uobsw(2,n)(1:j2)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Matrix column variables.
c
      write (newin,3470)
 3470 format('* Matrix species or entities')
c
      do krow = 1,kdim
        j2 = ilnobl(uzveci(krow))
        write (newin,3490) uzveci(krow)(1:j2)
 3490   format(3x,a)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Matrix variable values.
c
      write (newin,3510)
 3510 format('*   Values of matrix variables')
      do kcol = 1,kdim
        write (newin,3530) zvclgi(kcol)
 3530   format(3x,1pe22.15)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      qx = kprs.gt.0 .and. nprpti.gt.0 .and. nprsti.gt.0
      if (qx) then
        write (newin,3600)
 3600   format('*  Phases and species in the PRS')
        nsi = 0
        do npi = 1,nprpti
          write (newin,3630) uprphi(npi),mprphi(npi)
 3630     format(1x,a24,6x,1pe22.15)
          do iki = 1,iktmax
            if (uprspi(nsi + 1)(25:48) .eq. uprphi(npi)(1:24)) then
              nsi = nsi + 1
              write (newin,3670) uprspi(nsi)(1:24),mprspi(nsi)
 3670         format(3x,a24,6x,1pe22.15)
            else
              j3 = ilnobl(uendit)
              write (newin,3660) uendit(1:j3)
 3660         format(3x,a)
              go to 300
            endif
          enddo
  300     continue
        enddo
        j3 = ilnobl(uendit)
        write (newin,3620) uendit(1:j3)
 3620   format(1x,a)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
