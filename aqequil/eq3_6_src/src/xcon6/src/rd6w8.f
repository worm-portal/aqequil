      subroutine rd6w8(awmaxi,awmini,cbsri,cdac,cesri,cgexj,
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
c     This subroutine reads the EQ6 input file in compact ("W") format
c     for version 8.0.
c
c     This subroutine is a near-clone of EQ6/rd6inw.f. However, the
c     present subroutine embodies only a pure read function (it does
c     only minimal checking of what is read to ensure that what
c     follows is readable). EQ6/rd6inw.f differs in that it also
c     writes an instant echo of what is read to the EQ6 output file.
c
c     The calling sequence of this subroutine is identical to that of
c     EQ6/rd6inw.f, EQ6/rd6ind.f, and XCON6/rd6d8.f.
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
     $ kstpmx,kxt,nbti,nffg,nert,net,nobswt,nprob,nprpti,nprsti,nrct,
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
     $ egersi(ietmax,jetmax,nertmx),fkrc(nrctmx),
     $ hact(imchmx,2,nrctmx),modr(nrctmx),moffg(nffgmx),morr(nrctmx),
     $ mprphi(nprpmx),mprspi(nprsmx),mtbaqi(nbtmax),mtbi(nbtmax),
     $ mwtges(netmax),ptk(nptkmx),rkb(imchmx,2,nrctmx),
     $ rxbari(iktmax,nxrtmx),sfcar(nrctmx),ssfcar(nrctmx),tgexp(netmax),
     $ trkb(imchmx,2,nrctmx),ttk(nttkmx),vreac(nrctmx),
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
      integer i,iki,ikti,je,jei,jeti,jj,j2,j3,kcol,krow,n,nbi,nci,ne,
     $ nei,ner,neti,nn,npi,nrc
c
      integer ilnobl
c
      character*80 uline
      character*24 uxf,uxg,uxs
      character*8 uendit,uxe,ux8
c
      real*8 mx,xx
c
c-----------------------------------------------------------------------
c
      data uendit /'endit.  '/
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
c     qend   = .true if the end of the input file has been encountered
c     qrderr = .true if the current problem can't be read because of
c                a read format error or a dimensional overflow
c
      qend = .false.
      qrderr = .false.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Main title.
c
c     Note: if there are exactly ntitpa lines in the title, the title
c     need not be terminated by an 'endit.'. The 'endit.' if present
c     is here not considered to be part of the title.
c
      ntitl1 = 0
      read (ninpts,1000,end=100,err=990) uline
 1000 format(a80)
      go to 110
c
  100 qend = .true.
      go to 999
c
  110 continue
c
      if (uline(1:8) .ne. uendit(1:8)) then
        ntitl1 = 1
        utitl1(1) = uline
        do n = 2,ntitmx
          read (ninpts,1000,err=990) uline
          if (uline(1:8) .eq. uendit(1:8)) go to 120
          ntitl1 = n
          utitl1(n) = uline
        enddo
  120   continue
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Temperature parameters.
c     Note: ttk(1) = ttk1, etc.
c
      read (ninpts,1100,err=990) jtemp,tempcb,ttk(1),ttk(2)
 1100 format(12x,i2,/12x,e12.5,/2(12x,e12.5))
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Pressure parameters.
c     Note: ptk(1) = ptk1, etc.
c
      read (ninpts,1120,err=990) jpress,pressb,ptk(1),ptk(2)
 1120 format(12x,i2,/12x,e12.5,/2(12x,e12.5))
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Number of reactants.
c
      read (ninpts,1480,err=990) nrct
 1480 format(12x,i2)
c
      if (nrct .gt. nrctmx) then
        write (nttyo,1520) nrctmx
 1520   format(/' * Error - (XCON6/rd6w8) Have too many reactants',
     $  /7x,'The code is only dimensioned for ',i4,' reactants.',
     $  /7x,'Reduce the number of reactants or increase the',
     $  /7x,'dimensioning parameter nrctpa.')
        go to 990
      endif
c
c     Reactants.
c
      nsrt = 0
      nxrt = 0
      nert = 0
c
      do nrc = 1,nrct
c
c       Name, flags, and masses.
c
        read (ninpts,1530,err=990) ureac(nrc),jcode(nrc),jreac(nrc),
     $  morr(nrc),modr(nrc)
 1530   format(12x,a24,/12x,i2,22x,i2,/2(12x,e12.5))
        j2 = ilnobl(ureac(nrc))
c
        if (jcode(nrc) .eq. 1) then
c
c         Have a solid solution reactant.
c
          nxrt = nxrt + 1
c
          if (nxrt .gt. nxrtmx) then
            write (nttyo,1550) nxrtmx
 1550       format(/' * Error - (XCON6/rd6w8) Have too many solid',
     $      /7x,'solution reactants. The code is only dimensioned',
     $      /7x,'for ',i4,' such reactants. Reduce the number of',
     $      /7x,'such reactants or increase the dimensioning',
     $      ' parameter nxrtpa.')
            go to 990
          endif
c
          ikti = 0
          do i = 1,iktmax + 1
            read (ninpts,1560,err=990) uxs,xx
 1560       format(3x,a24,3x,e12.5)
c
            if (uxs(1:8) .eq. uendit(1:8)) go to 150
c
            ikti = ikti + 1
c
            if (ikti .gt. iktmax) then
              write (nttyo,1590) ureac(nrc)(1:j2),iktmax
 1590         format(/' * Error - (XCON6/rd6w8) Have too many',
     $        ' end-members',/7x,'in the solid solution reactant ',
     $        a,'.',/7x,'The code is only dimensioned for ',
     $        i4,' end-members per',/7x,'solid solution. Reduce',
     $        ' the number of end-members or',
     $        /7x,'increase the dimensioning parameter iktpar.')
              go to 990
            endif
c
            ucxri(ikti,nxrt) = uxs
            rxbari(ikti,nxrt) = xx
          enddo
  150     ixrti(nxrt) = ikti
        elseif (jcode(nrc) .eq. 2) then
c
c         Have a special reactant.
c
          nsrt = nsrt + 1
c
          if (nsrt .gt. nsrtmx) then
            write (nttyo,1600) nsrtmx
 1600       format(/' * Error - (XCON6/rd6w8) Have too many special',
     $      /7x,'reactants. The code is only dimensioned for ',i4,
     $      /7x,'such reactants. Reduce the number of such reactants',
     $      /7x,'or increase the dimensioning parameter nsrtpa.')
            go to 990
          endif
c
          read (ninpts,1610,err=990) vreac(nrc)
 1610     format(12x,e12.5)
c
          nci = 0
          do n = 1,nctmax + 1
            read (ninpts,1630,err=990) uxe,xx
 1630       format(3x,a8,3x,e22.15)
            if (uxe(1:8) .eq. uendit(1:8)) go to 160
c
            nci = nci + 1
c
            if (nci .gt. nctmax) then
              write (nttyo,1660) ureac(nrc)(1:j2),nctmax
 1660         format(/' * Error - (XCON6/rd6w8) Have too many',
     $        ' chemical',/7x,'elements in the special reactant ',a,
     $        '.',/7x,'The code is only dimensioned for ',i4,
     $        ' elements.',/7x,'Reduce the number of elements or',
     $        ' increase the',/7x,'dimensioning parameter nctpar.')
              go to 990
            endif
c
            uesri(nci,nsrt) = uxe
            cesri(nci,nsrt) = xx
          enddo
  160     iesrti(nsrt) = nci
c
          nbi = 0
          do n = 1,nbt1mx + 1
            read (ninpts,1670,err=990) uxs,xx
 1670       format(3x,a24,3x,e22.15)
            if (uxs(1:8) .eq. uendit(1:8)) go to 170
c
            nbi = nbi + 1
c
            if (nbi .gt. nbt1mx) then
              write (nttyo,1700) ureac(nrc)(1:j2),nbtmax
 1700         format(/' * Error - (XCON6/rd6w8) Have too many basis',
     $        ' basis species in the',/7x,'in the reaction for the',
     $        ' special reactant ',a,'.',/7x,'The code is only',
     $        ' dimensioned for ',i4,' basis species. Increase',
     $        /7x,'the dimensioning parameter nbtpar.')
              go to 990
            endif
c
            ubsri(nbi,nsrt) = uxs
            cbsri(nbi,nsrt) = xx
          enddo
  170     ibsrti(nsrt) = nbi
        elseif (jcode(nrc) .eq. 5) then
c
c         Have a generic ion exchanger reactant.
c
          nert = nert + 1
          ner = nert
c
          if (nert .gt. nertmx) then
            write (nttyo,1710) nertmx
 1710       format(/' * Error - (XCON6/rd6w8) Have too many generic',
     $      ' ion exchanger',/7x,'reactants. The code is only',
     $      ' dimensioned for ',i4,' such reactants.',/7x,'Reduce',
     $      ' the number of such reactants or increase the',
     $      /7x,'dimensioning parameter nertpa.')
            go to 990
          endif
c
          read (ninpts,1712,err=990) ugermo(ner)
 1712     format(12x,a24)
c
          jeti = 0
          do jj = 1,jetmax + 1
            read (ninpts,1716,err=990) uxs
 1716       format(3x,a24)
            if (uxs(1:8) .eq. uendit(1:8)) go to 154
c
            jeti = jeti + 1
            jei = jeti
c
            if (jeti .gt. jetmax) then
              write (nttyo,1718) ureac(nrc)(1:j2),jetmax
 1718         format(/' * Error - (XCON6/rd6w8) Have too many',
     $        ' exchange sites',/7x,'in the generic ion exchanger',
     $        ' reactant ',a,'.',/7x,'The code is only dimensioned',
     $        ' for ',i4,' exchange sites per',/7x,'generic ion',
     $        ' exchanger. Reduce the number of exchange sites or',
     $        /7x,'increase the dimensioning parameter jetpar.')
              go to 990
            endif
c
            ugerji(jei,ner) = uxs
c
            neti = 0
            do nn = 1,netmax + 1
              read (ninpts,1720,err=990) uxs,xx
 1720         format(6x,a24,3x,e12.5)
              if (uxs(1:8) .eq. uendit(1:8)) go to 152
c
              neti = neti + 1
              nei = neti
c
              if (neti .gt. netmax) then
                j3 = ilnobl(ugerji(jei,nert))
                write (nttyo,1724) ureac(nrc)(1:j2),
     $          ugerji(jei,ner)(1:j3),netmax
 1724           format(/' * Error - (XCON6/rd6w8) Have too many',
     $          ' species on exchange',/7x,'site ',a,' of the generic',
     $          ' ion exchanger reactant',/7x,a,'. The code is only',
     $          ' dimensioned for species per',/7x,'exchange site.',
     $          ' Reduce the number of end-members or increase',
     $          /7x,'the dimensioning parameter netpar.')
                go to 990
              endif
c
              ugersi(nei,jei,ner) = uxs
              egersi(nei,jei,ner) = xx
              xgersi(nei,jei,ner) = xx
            enddo
  152       igerti(jei,ner) = neti
          enddo
  154     jgerti(ner) = jeti
        endif
c
c       Surface area parameters.
c
        read (ninpts,1740,err=990) nsk(nrc),sfcar(nrc),ssfcar(nrc),
     $  fkrc(nrc)
 1740   format(12x,i2,22x,e12.5,12x,e12.5,/12x,e12.5)
c
c       Rate law codes: nrk(1,nrc) = forward rate law code,
c       nrk(2,nrc) = backward rate law code.
c
        read (ninpts,1760,err=990) nrk(1,nrc),nrk(2,nrc)
 1760   format(12x,i2,22x,i2)
c
c       Read forward (dissolution, dissociation) rate data.
c
        if (nrk(1,nrc) .eq. -1) then
c
c         Use the backward rate law.
c
          continue
        elseif (nrk(1,nrc) .eq. 0) then
c
c         Illegal value.
c
          write (nttyo,1775) ureac(nrc)(1:j2)
 1775     format(/' * Error - (XCON6/rd6w8) The forward rate law code',
     $    /7x,'has an illegal value of 0 for ',a,'.')
          go to 990
        elseif (nrk(1,nrc) .eq. 1) then
c
c         Arbitrary kinetics (relative rates, indifferent to time).
c
          imech(1,nrc) = 3
c
          if (imech(1,nrc) .gt. imchmx) then
            write (nttyo,1780) ureac(nrc)(1:j2),imchmx
 1780       format(/' * Error - (XCON6/rd6w8) Have too many rate',
     $      /7x,'constants or corresponding mechanisms in the forward',
     $      /7x,'rate law for reactant ',a,'. The code is only',
     $      /7x,'dimensioned for ',i2,' rate constants per rate law.',
     $      /7x,'Reduce the number of rate constants or increase the',
     $      /7x,'dimensioning parameter imchpa.')
            go to 990
          endif
c
          read (ninpts,1790,err=990) (rkb(i,1,nrc), i = 1,3)
 1790     format(3(12x,e12.5))
        elseif (nrk(1,nrc) .eq. 2) then
c
c         Transition state rate law. Up to imchmx parallel mechanisms
c         are allowed.
c
          read (ninpts,1810,err=990) imech(1,nrc)
 1810     format(12x,i2)
c
          if (imech(1,nrc) .gt. imchmx) then
            write (nttyo,1780) ureac(nrc)(1:j2),imchmx
            go to 990
          endif
c
          do i = 1,imech(1,nrc)
            read (ninpts,1830,err=990) rkb(i,1,nrc),trkb(i,1,nrc),
     $      iact(i,1,nrc)
 1830       format(12x,e12.5,12x,e12.5,12x,i2)
c
            read (ninpts,1850,err=990) eact(i,1,nrc),hact(i,1,nrc)
 1850       format(12x,e12.5,12x,e12.5)
c
            read (ninpts,1870,err=990) ndact(i,1,nrc),csigma(i,1,nrc)
 1870       format(12x,i2,22x,e12.5)
c
            if (ndact(i,1,nrc) .gt. ndctmx) then
              write (nttyo,1890) i,ureac(nrc)(1:j2),ndctmx
 1890         format(/' * Error - (XCON6/rd6w8) Have too many',
     $        /7x,'species in the activity product in term ',i2,
     $        /7x,'of the forward direction rate law for reactant',
     $        /7x,a,'. The code is only dimensioned for ',i3,
     $        /7x,'such species. Reduce the number of such species',
     $        /7x,'or increase the dimensioning parameter ndctpa.')
              go to 990
            endif
c
            do n = 1,ndact(i,1,nrc)
              read (ninpts,1910,err=990) udac(n,i,1,nrc),
     $        cdac(n,i,1,nrc)
 1910         format(12x,a24,12x,e12.5)
            enddo
          enddo
        elseif (nrk(1,nrc) .eq. 3) then
c
c         Linear rate law.
c
          imech(1,nrc) = 1
c
          if (imech(1,nrc) .gt. imchmx) then
            write (nttyo,1780) ureac(nrc)(1:j2),imchmx
            go to 990
          endif
c
          i = 1
          read (ninpts,1830,err=990) rkb(i,1,nrc),trkb(i,1,nrc),
     $    iact(i,1,nrc)
          read (ninpts,1850,err=990) eact(i,1,nrc),hact(i,1,nrc)
        else
c
          write (nttyo,1950) nrk(1,nrc),ureac(nrc)(1:j2)
 1950     format(/' * Error - (XCON6/rd6w8) The forward rate law code',
     $    /7x,'has an unrecognized value of ',i2,' for ',a,'.')
          go to 990
c
        endif
c
c       Read backward (precipitation, formation) rate data.
c
        if (nrk(2,nrc) .eq. -1) then
c
c         Use the forward rate law.
c
          continue
        elseif (nrk(2,nrc) .eq. 0) then
c
c         Use instantaneous partial equilibrium.
c
          continue
        elseif (nrk(2,nrc) .eq. 1) then
c
c         Arbitrary kinetics.
c
          imech(2,nrc) = 3
c
          if (imech(2,nrc) .gt. imchmx) then
            write (nttyo,1960) ureac(nrc)(1:j2),imchmx
 1960       format(/' * Error - (XCON6/rd6w8) Have too many rate',
     $      /7x,'constants or corresponding mechanisms in the backward',
     $      /7x,'rate law for reactant ',a,'. The code is only',
     $      /7x,'dimensioned for ',i2,' rate constants per rate law.',
     $      /7x,'Reduce the number of rate constants or increase the',
     $      /7x,'dimensioning parameter imchpa.')
            go to 990
          endif
c
          read (ninpts,1790,err=990) (rkb(i,2,nrc), i = 1,3)
        elseif (nrk(2,nrc) .eq. 2) then
c
c         Transition state rate law.
c
          read (ninpts,1810,err=990) imech(2,nrc)
c
          if (imech(2,nrc) .gt. imchmx) then
            write (nttyo,1960) ureac(nrc)(1:j2),imchmx
            go to 990
          endif
c
          do i = 1,imech(2,nrc)
            read (ninpts,1830,err=990) rkb(i,2,nrc),trkb(i,2,nrc),
     $      iact(i,2,nrc)
c
            read (ninpts,1850,err=990) eact(i,2,nrc),hact(i,2,nrc)
c
            read (ninpts,1870,err=990) ndact(i,2,nrc),csigma(i,2,nrc)
c
            if (ndact(i,2,nrc) .gt. ndctmx) then
              write (nttyo,1970) i,ureac(nrc)(1:j2),ndctmx
 1970         format(/' * Error - (XCON6/rd6w8) Have too many',
     $        /7x,'species in the activity product in term ',i2,
     $        /7x,'of the backward direction rate law for reactant',
     $        /7x,a,'. The code is only dimensioned for ',i3,
     $        /7x,'such species. Reduce the number of such species',
     $        /7x,'or increase the dimensioning parameter ndctpa.')
              go to 990
            endif
c
            do n = 1,ndact(i,2,nrc)
              read (ninpts,1910,err=990) udac(n,i,2,nrc),
     $        cdac(n,i,2,nrc)
            enddo
          enddo
        elseif (nrk(2,nrc) .eq. 3) then
c
c         Linear rate law.
c
          imech(2,nrc) = 1
c
          if (imech(2,nrc) .gt. imchmx) then
            write (nttyo,1960) ureac(nrc)(1:j2),imchmx
            go to 990
          endif
c
          i = 1
          read (ninpts,1830,err=990) rkb(i,2,nrc),trkb(i,2,nrc),
     $    iact(i,2,nrc)
          read (ninpts,1850,err=990) eact(i,2,nrc),hact(i,2,nrc)
        else
c
          write (nttyo,1980) nrk(2,nrc),ureac(nrc)(1:j2)
 1980     format(/' * Error - (XCON6/rd6w8) The backward rate law code',
     $    /7x,'has an unrecognized value of ',i2,' for ',a,'.')
          go to 990
c
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Starting, minimum, and maximum values of key run parameters.
c
      read (ninpts,1150,err=990) xistti,ximaxi,tistti,timmxi
 1150 format(2(12x,e12.5),/2(12x,e12.5))
c
      read (ninpts,1160,err=990) phmini,phmaxi,ehmini,ehmaxi,o2mini,
     $ o2maxi,awmini,awmaxi,kstpmx
 1160 format(2(12x,e12.5),/2(12x,e12.5),/2(12x,e12.5),/2(12x,e12.5),
     $ /12x,i12)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print interval parameters.
c
      read (ninpts,1170,err=990) dlxprn,dlxprl,dltprn,dltprl
 1170 format(2(12x,e12.5),12x,/2(12x,e12.5))
c
      read (ninpts,1180,err=990) dlhprn,dleprn,dloprn,dlaprn,ksppmx
 1180 format(2(12x,e12.5),12x,/2(12x,e12.5),/12x,i12)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Plot interval parameters.
c
      read (ninpts,1190,err=990) dlxplo,dlxpll,dltplo,dltpll
 1190 format(2(12x,e12.5),/2(12x,e12.5))
c
      read (ninpts,1195,err=990) dlhplo,dleplo,dloplo,dlaplo,ksplmx
 1195 format(2(12x,e12.5),/2(12x,e12.5),/12x,i12)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iopt model option switches.
c     Note: iopt(1) = iopt1, etc.
c
      read (ninpts,1210,err=990) (iopt(i), i = 1,20)
 1210 format(12x,10i5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iopr print option switches.
c     Note: iopr(1) = iopr1, etc.
c
      read (ninpts,1230,err=990) (iopr(i), i = 1,20)
 1230 format(12x,10i5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iodb debugging print option switches.
c     Note: iodb(1) = iodb1, etc.
c
      read (ninpts,1250,err=990) (iodb(i), i = 1,20)
 1250 format(12x,10i5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Number of nxopt options.
c
      read (ninpts,1300,err=990) nxopt
 1300 format(12x,i2)
c
      if (nxopt .gt. nxopmx) then
        write (nttyo,1320) nxopmx
 1320   format(/' * Error - (XCON6/rd6w8) Have too many mineral',
     $  /7x,'subset-selection suppression options. The code is',
     $  /7x,'only dimensioned for ',i3,' such options. Reduce the',
     $  /7x,'number of options or increase the dimensioning',
     $  /7x,'parameter nxoppa.')
        go to 990
      endif
c
c     Nxopt options.
c
      do n = 1,nxopt
        read (ninpts,1330,err=990) uxopt(n),uxcat(n)
 1330   format(12x,a6,1x,a24)
      enddo
c
      if (nxopt .gt. 0) then
        read (ninpts,1350,err=990) nxopex
 1350   format(12x,i2)
c
        if (nxopex .gt. nxpemx) then
          write (nttyo,1370) nxpemx
 1370     format(/' * Error - (XCON6/rd6w8) Have too many',
     $    /7x,'exceptions specified to the mineral subset-selection',
     $    /7x,'suppression options. The code is only dimensioned',
     $    /7x,'for ',i3,'exceptions. Reduce the number of exceptions',
     $    /7x,'or increase the dimensioning parameter nxpepa.')
          go to 990
        endif
c
        do n = 1,nxopex
          read (ninpts,1380,err=990) uxopex(n)
 1380     format(12x,a24)
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Number of nffg options.
c
      read (ninpts,1400,err=990) nffg
 1400 format(12x,i2)
c
      if (nffg .gt. nffgmx) then
        write (nttyo,1420) nffgmx
 1420   format(/' * Error - (XCON6/rd6w8) Have too many gases whose',
     $  /7x,'fugacities are to be fixed. The code is only dimensioned',
     $  /7x,'for ',i4,' such gases. Reduce the number of gases or',
     $  /7x,'increase the dimensioning parameter nffgpa.')
        go to 990
      endif
c
c     Nffg options.
c
      do n = 1,nffg
        read (ninpts,1430,err=990) uffg(n),moffg(n),xlkffg(n)
 1430   format(12x,a24,/12x,e12.5,12x,e12.5)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum finite-difference order.
c
      read (ninpts,2010,err=990) nordmx
 2010 format(12x,i3)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Newton-Raphson convergence tolerances.
c
      read (ninpts,2020,err=990) tolbt,toldl
 2020 format(2(12x,e12.5))
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum number of Newton-Raphson iterations.
c
      read (ninpts,2010,err=990) itermx
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Search/find tolerance.
c
      read (ninpts,2022,err=990) tolxsf
 2022 format(12x,e12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Saturation tolerance.
c
      read (ninpts,2024,err=990) tolsat
 2024 format(12x,e12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum number of phase assemblage tries.
c
      read (ninpts,2102,err=990) ntrymx
 2102 format(12x,i3)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Zero-order step size (in Xi).
c
      read (ninpts,2080,err=990) dlxmx0
 2080 format(12x,e12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     PRS transfer interval in Xi.
c
      read (ninpts,2000,err=990) dlxdmp
 2000 format(12x,e12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Process the bottom half of the current input file.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Secondary title.
c
      ntitl2 = 0
      read (ninpts,1000,err=990) uline
c
      if (uline(1:8) .ne. uendit(1:8)) then
        ntitl2 = 1
        utitl2(1) = uline
        do n = 2,ntitmx
          read (ninpts,1000,err=990) uline
          if (uline(1:8) .eq. uendit(1:8)) go to 220
          utitl2(n) = uline
          ntitl2 = n
        enddo
  220   continue
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Special basis switches.
c
      read (ninpts,2620,err=990) nsbswt
 2620 format(12x,i3)
c
      do n = 1,nsbswt
        read (ninpts,2640,err=990) usbsw(1,n)
 2640   format(9x,a48)
c
        read (ninpts,2660,err=990) usbsw(2,n)
 2660   format(15x,a48)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Original temperature.
c
      read (ninpts,2680,err=990) tempci
 2680 format(12x,e12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Original pressure.
c
      read (ninpts,2700,err=990) pressi
 2700 format(12x,e12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Ion exchanger creation. This section reads the directives for
c     creating ion exchange phases and species and their associated
c     intrinsic (including thermodynamic) properties. The data is
c     essentially that which might be read from a supporting data file.
c
      read (ninpts,2795) ux8
 2795 format(12x,a8)
      call lejust(ux8)
      call locase(ux8)
      qgexsh = ux8(1:2).eq.'.t' .or. ux8(1:1).eq. 't'
c
      read (ninpts,2800,err=990) net
 2800 format(12x,i3)
c
      if (net .gt. netmax) then
        write (nttyo,2820) netmax,net
 2820   format(/' * Error - (XCON6/rd6w8) Have exceeded the maximum',
     $  ' number of ',i3,/7x,'generic ion exchange phases while',
     $  ' reading the data to create',/7x,'such phases. Increase',
     $  ' the dimensioning parameter netpar',/7x,'to at least ',i3,'.')
        go to 990
      endif
c
      do ne = 1,net
        read (ninpts,2830,err=990) ugexp(ne)
 2830   format(12x,a24)
        j2 = ilnobl(ugexp(ne))
c
        read (ninpts,2850,err=990) mwtges(ne)
 2850   format(12x,e12.5)
c
        read (ninpts,2830,err=990) ugexmo(ne)
        j3 = ilnobl(ugexmo(ne))
c
        read (ninpts,2850,err=990) tgexp(ne)
c
        read (ninpts,2890,err=990) jgext(ne)
 2890   format(12x,i3)
c
        if (jgext(ne) .gt. jetmax) then
          write (nttyo,2910) jetmax,ugexp(ne)(1:j2),jgext(ne)
 2910     format(/' * Error - (XCON6/rd6w8) Have exceeded the maximum',
     $    ' number of ',i3,/7x,'exchange sites on a generic ion',
     $    ' exchange phase while reading',/7x,'the data to create',
     $    a,'. Increase the',/7x,'dimensioning parameter jetpar to',
     $    ' at least ',i3,'.')
          go to 990
        endif
c
        do je = 1,jgext(ne)
          read (ninpts,2930,err=990) ugexj(je,ne)
 2930     format(12x,a8)
          j3 = ilnobl(ugexj(je,ne))
c
          read (ninpts,2950,err=990) cgexj(je,ne),zgexj(je,ne)
 2950     format(12x,e12.5,12x,e12.5)
c
          read (ninpts,2890,err=990) ngexrt(je,ne)
c
          if (ngexrt(je,ne) .gt. ietmax) then
            write (nttyo,3020) netmax,ugexj(je,ne)(1:j3),
     $      ugexp(ne)(1:j2),ngexrt(je,ne)
 3020       format(/' * Error - (XCON6/rd6w8) Have exceeded the',
     $      ' maximum number of ',i3,/7x,'reactions for a site',
     $      ' belonging to a generic ion exchange',/7x,'phase while',
     $      ' reading the data for site ',a,' of exchange phase'
     $      /7x,a,'. Increase the dimensioning parameter',
     $      /7x,'ietpar to at least ',i3,'.')
            go to 990
          endif
c
          do n = 1,ngexrt(je,ne)
            read (ninpts,3030,err=990) ugexr(n,je,ne)
 3030       format(12x,a56)
            read (ninpts,3050,err=990) xlkgex(n,je,ne),uxkgex(n,je,ne)
 3050       format(12x,e12.5,12x,a8)
            read (ninpts,3050,err=990) xhfgex(n,je,ne),uhfgex(n,je,ne)
            read (ninpts,3050,err=990) xvfgex(n,je,ne),uvfgex(n,je,ne)
          enddo
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Number of nxmod options.
c
      read (ninpts,3120,err=990) nxmod
 3120 format(12x,i2)
c
      if (nxmod .gt. nxmdmx) then
        write (nttyo,3140) nxmdmx
 3140   format(/' * Error - (XCON6/rd6w8) Have too many nxmod',
     $  /7x,'alter/suppress options. The code is only dimensioned',
     $  /7x,'for ',i2,' such options. Reduce the number of such',
     $  /7x,'options or increase the dimensioning parameter nxmdpa.')
        go to 990
      endif
c
c     Nxmod options.
c
      do n = 1,nxmod
        read (ninpts,3150,err=990) uxmod(n),kxmod(n),xlkmod(n)
 3150   format(12x,a48,/12x,i2,22x,e12.5)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iopg options.
c     Note: iopg(1) = iopg1, etc.
c
      read (ninpts,3170,err=990) (iopg(i), i = 1,20)
 3170 format(12x,10i5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Index limits.
c
      read (ninpts,3220,err=990) kct,kbt,kmt,kxt,kdim,kprs
 3220 format(3(12x,i2,10x),/3(12x,i2,10x))
c
      nbti = kbt
c
      if (kct .gt. nctmax) then
        write (nttyo,3240) nctmax
 3240   format(/' * Error - (XCON6/rd6w8) Have too many chemical',
     $  /7x,'elements present. The code is only dimensioned',
     $  /7x,'for ',i3,' elements. Reduce the number of elements',
     $  /7x,'or increase the dimensioning parameter nctpar.')
        go to 990
      endif
c
      if (kbt .gt. nbtmax) then
        write (nttyo,3250) nbtmax
 3250   format(/' * Error - (XCON6/rd6w8) Have too many basis',
     $  /7x,'species present. The code is only dimensioned',
     $  /7x,'for ',i3,' basis species. Reduce the number of elements',
     $  /7x,'or increase the dimensioning parameter nctpar.')
        go to 990
      endif
c
      if (kdim .gt. kmax) then
        write (nttyo,3260) kmax
 3260   format(/' * Error - (XCON6/rd6w8) Have too many matrix',
     $  /7x,'variables. The code is only dimensioned for ',i3,
     $  /7x,'matrix variables. Reduce the number of such variables',
     $  /7x,'or increase the dimensioning parameter kpar.')
        go to 990
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Species for which mass balances are defined.
c
c       ubmtbi = names of the data file basis species for which mass
c                  balances are defined
c       jflgi  = jflag input for basis species
c          0 = Retain as an active basis species
c         30 = Convert to a dependent species; fold the mass balance
c                total for this species into the mass balance totals
c                of basis species which remain active
c
      do nbi = 1,nbti
        read (ninpts,3340,err=990) ubmtbi(nbi),jflgi(nbi)
 3340   format(3x,a48,3x,i2)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Mass balance totals.
c
c       mtbi   = total number of moles of basis species in the
c                  equilibrium system (ES)
c       mtbaqi = total number of moles of basis species in the
c                  aqueous solution
c
      do nbi = 1,nbti
        read (ninpts,3410,err=990) mtbi(nbi),mtbaqi(nbi)
 3410   format(6x,e22.15,6x,e22.15)
      enddo
      read (ninpts,3430,err=990) electr
 3430 format(34x,e22.15)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Ordinary basis switches.
c
      read (ninpts,2620,err=990) nobswt
c
      do n = 1,nobswt
        read (ninpts,2640,err=990) uobsw(1,n)
c
        read (ninpts,2660,err=990) uobsw(2,n)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Matrix column variables.
c
c       uzveci = names of species or properties associated with the
c                  matrix column variables
c
      do krow = 1,kdim
        read (ninpts,3480,err=990) uzveci(krow)
 3480   format(3x,a48)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Matrix variable values.
c
c       zvclg1 = values of matrix variables
c
      do kcol = 1,kdim
        read (ninpts,3520,err=990) zvclgi(kcol)
 3520   format(3x,e22.15)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize any non-zero values for the physically removed
c     system. Here mprphi and mprspi are arrays for the number of
c     moles of phases and species, respectively, in the PRS.
c
      nprpti = 0
      nprsti = 0
      if (kprs .gt. 0) then
        do npi = 1,nprpmx + 1
          read (ninpts,3610,err=990) uxf,mx
 3610     format(1x,a24,6x,e22.15)
          if (uxf(1:8) .eq. uendit(1:8)) go to 310
          nprpti = nprpti + 1
c
          if (nprpti .gt. nprpmx) then
            write (nttyo,3640) nprpmx
 3640       format(/' * Error - (XCON6/rd6w8) Have too many phases',
     $      /7x,'in the physically removed system (PRS) on the',
     $      /7x,'input file. The code is only dimensioned for ',i3,
     $      /7x,'such phases. Reduce the number of such phases',
     $      /7x,'or increase the dimensioning parameter nprppa.')
            go to 990
          endif
c
          uprphi(nprpti)(1:24) = uxf
          mprphi(nprpti) = mx
c
          do iki = 1,iktmax + 1
            read (ninpts,3650,err=990) uxg,mx
 3650       format(3x,a24,6x,e22.15)
            if (uxg(1:8) .eq. uendit(1:8)) go to 300
c
            if (iki .gt. iktmax) then
              j2 = ilnobl(uxf)
              write (nttyo,3680) uxf(1:j2),iktmax
 3680         format(/' * Error - (XCON6/rd6w8) Have too many',
     $        ' end-members',/7x,'in the PRS solid solution ',a,'.',
     $        /7x,'The code is only dimensioned for ',
     $        i4,' end-members per',/7x,'solid solution. Reduce',
     $        ' the number of end-members or',
     $        /7x,'increase the dimensioning parameter iktpar.')
              go to 990
            endif
c
            nprsti = nprsti + 1
c
            if (nprsti .gt. nprsmx) then
              write (nttyo,3690) nprsmx
 3690         format(/' * Error - (XCON6/rd6w8) Have too many species',
     $        /7x,'in the physically removed system (PRS) on the',
     $        /7x,'input file. The code is only dimensioned for ',i3,
     $        /7x,'such species. Reduce the number of such species',
     $        /7x,'or increase the dimensioning parameter nprspa.')
              go to 990
            endif
c
            uprspi(nprsti)(1:24) = uxg
            uprspi(nprsti)(25:48) = uxf
            mprspi(nprsti) = mx
          enddo
  300     continue
        enddo
c
  310   continue
      endif
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
