      subroutine rd6inw(awmaxi,awmini,cbsri,cdac,cesri,cgexj,
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
c     This subroutine reads the EQ6 input file in compact ("W") format.
c
c     This subroutine is a near-clone of XCON6/rd6w8.f. The present
c     subroutine differs from the latter, in that in addition to the
c     pure read function, it writes an instant echo of what is read to
c     the output file, and sandwiches this between prefatory and
c     ending messages. To maintain close consistency with XCON6/rd6w8.f,
c     this subroutine assigns no default values and only performs such
c     checking of what is read to ensure that what follows is readable.
c
c     The calling sequence of this subroutine is identical to that of
c     EQ6/rd6ind.f, XCON6/rd6w8.f, and XCON6/rd6d8.f.
c
c     This subroutine is called by:
c
c       EQ6/eq6.f
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
     $ ixrti(nxrtmx),jcode(nrctmx),jflgi(nbtmax),jgerti(nertmx),
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
      integer i,iki,ikti,je,jei,jeti,jj,j2,j3,j4,kcol,krow,n,nbi,nci,
     $ ne,nei,ner,neti,nn,npi,nrc
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
c       utitl1 = main title
c       ntitl1 = number of lines in the main title
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
  110 write (noutpt,1010) nprob
      write (nttyo,1010) nprob
 1010 format(//' Reading problem ',i3,' from the input file ...',/)
c
      j2 = ilnobl(uline)
      j2 = min(j2,79)
      write (noutpt,1020) uline
 1020 format(1x,a)
c
      if (uline(1:8) .ne. uendit(1:8)) then
        ntitl1 = 1
        utitl1(1) = uline
      endif
c
      if (ntitl1 .gt. 0) then
        do n = 2,ntitmx
          read (ninpts,1000,err=990) uline
          j2 = ilnobl(uline)
          j2 = min(j2,79)
          write (noutpt,1020) uline(1:j2)
          if (uline(1:8) .eq. uendit(1:8)) go to 120
          ntitl1 = n
          utitl1(n) = uline
        enddo
  120   continue
c
c       Write the first 5 lines of the input problem title to the
c       screen file.
c
        write (nttyo,1040)
 1040   format(3x,'The input problem title is (first 5 lines',
     $   ' maximum):',/)
c
        i = min(ntitl1,5)
        do n = 1,i
          j2 = ilnobl(utitl1(n))
          j2 = min(j2,74)
          write (nttyo,1050) utitl1(n)(1:j2)
 1050     format(5x,a)
        enddo
      endif
c
      write (nttyo,1060)
 1060 format(/3x,'Continuing to read the problem input ...')
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Temperature parameters.
c     Note: ttk(1) = ttk1, etc.
c
c       jtemp  = temperature tracking code:
c
c         0 = Constant temperature:
c               tempc = tempcb
c
c         1 = Linear tracking in Xi:
c               tempc = tempcb + ttk(1)*xi1
c
c         2 = Linear tracking in time (iopt(2) must be 1):
c               tempc = tempcb + ttk(1)*time1
c
c         3 = Fluid mixing tracking:
c               tempc= (tempcb*ttk(1) + xi*ttk(2))/(xi + ttk(1)), where:
c                 ttk(1) = ratio of the mass of the starting aqueous
c                   to that of the aqueous solution being treated as a
c                   reactant; usually the mass of each will be close to
c                   1 kilogram, hence the value of ttk(1) will then be
c                   about 1.0.
c                 ttk(2) = temperature of the fluid being treated as a
c                   reactant ("Fluid 2")
c
c       tempcb = the base temperature, C
c
c       ttk    = temperature tracking coefficients
c
      read (ninpts,1100,err=990) jtemp,tempcb,ttk(1),ttk(2)
 1100 format(12x,i2,/12x,e12.5,/2(12x,e12.5))
      write (noutpt,1110) jtemp,tempcb,ttk(1),ttk(2)
 1110 format(6x,'jtemp= ',i2,/5x,'tempcb= ',1pe12.5,
     $ /7x,'ttk1= ',e12.5,6x,'ttk2= ',e12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Pressure parameters.
c     Note: ptk(1) = ptk1, etc.
c
c       jpress = pressure tracking code:
c
c         0 = Follow the data file's reference pressure curve:
c               press = func(data file, T)
c
c         1 = Follow the 1.013-bar/steam-saturation curve:
c               press = built-in func(T)
c
c         2 = Constant pressure:
c               press = pressb
c
c         3 = Linear tracking in Xi:
c               press = pressb + ptk(1)*xi1
c
c         4 = Linear tracking in time (iopt(2) must be 1):
c               press = pressb + ptk(1)*time1
c
c       pressb = the base pressure, bars
c
c       ptk    = pressure tracking coefficients
c
      read (ninpts,1120,err=990) jpress,pressb,ptk(1),ptk(2)
 1120 format(12x,i2,/12x,e12.5,/2(12x,e12.5))
      write (noutpt,1130) jpress,pressb,ptk(1),ptk(2)
 1130 format(4x,'jpress= ',i2,/5x,'pressb= ',1pe12.5,
     $ /7x,'ptk1= ',e12.5,6x,'ptk2= ',e12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Number of reactants.
c
c       nrct   = number of reactants
c
      read (ninpts,1480,err=990) nrct
 1480 format(12x,i2)
      write (noutpt,1490) nrct
 1490 format(7x,'nrct= ',i2)
c
      write (noutpt,1500)
 1500 format(' *-------------------------------------------------',
     $ '----------------------------')
c
      if (nrct .gt. nrctmx) then
        write (noutpt,1520) nrctmx
        write (nttyo,1520) nrctmx
 1520   format(/' * Error - (EQ6/rd6inw) Have too many reactants',
     $  /7x,'The code is only dimensioned for ',i4,' reactants.',
     $  /7x,'Reduce the number of reactants or increase the',
     $  /7x,'dimensioning parameter nrctpa.')
        go to 990
      endif
c
c     Reactants.
c
c       nxrt  = number of solid solution reactants
c       nsrt  = number of special reactants
c       nert  = number of ion exchanger reactants
c       ureac = name of reactant
c       morr  = moles of reactant remaining when Xi = xistti
c                 (the starting number of moles for the process
c                 corresponds to Xi = 0)
c       modr  = moles of reactant irreversibly destroyed when
c                 Xi = xistti (usually zero at Xi = 0)
c       jcode = reactant type flag:
c                 0 = mineral
c                 1 = solid solution
c                 2 = special reactant
c                 3 = aqueous species
c                 4 = gas
c       jreac = reactant status flag:
c                  0 = set to react (this is the only value a user
c                        should set; the code will write other values
c                        on pickup files)
c                 -1 = saturated, but the remaining reactant mass
c                        continues to react irreversibly; there is
c                        usually also a secondary product mass of the
c                        same species, so that the net rate of reaction
c                        is zero and the solution does not saturate
c                  1 = exhausted
c                  2 = saturated; the status of any remaining reactant
c                        mass is changed to that of a product phase;
c                        this only happens if iopt(1) = 0
c
      nxrt = 0
      nsrt = 0
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
        write (noutpt,1540) ureac(nrc)(1:j2),jcode(nrc),jreac(nrc),
     $  morr(nrc),modr(nrc)
 1540   format(3x,'reactant= ',a,/6x,'jcode= ',i2,15x,'jreac= ',i2,
     $  /7x,'morr= ',1pe12.5,6x,'modr= ',e12.5)
c
        if (jcode(nrc) .eq. 1) then
c
c         Have a solid solution reactant.
c
c           ucxri  = name of end-member of a solid solution reactant
c           rxbari = mole fraction of end-member of a solid solution
c                      reactant
c
          nxrt = nxrt + 1
c
          if (nxrt .gt. nxrtmx) then
            write (noutpt,1550) nxrtmx
            write (nttyo,1550) nxrtmx
 1550       format(/' * Error - (EQ6/rd6inw) Have too many solid',
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
            if (uxs(1:8) .eq. uendit(1:8)) then
              j3 = ilnobl(uxs)
              write (noutpt,1570) uxs(1:j3)
 1570         format(4x,a)
              go to 150
            endif
            write (noutpt,1580) uxs,xx
 1580       format(4x,a24,3x,1pe12.5)
c
            ikti = ikti + 1
c
            if (ikti .gt. iktmax) then
              write (noutpt,1590) ureac(nrc)(1:j2),iktmax
              write (nttyo,1590) ureac(nrc)(1:j2),iktmax
 1590         format(/' * Error - (EQ6/rd6inw) Have too many',
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
c           vreac  = volume, cc/mol
c           uesri  = name of element
c           cesri  =  moles chemical element per mole reactant
c           ubsri  = name of basis species in the chemical reaction
c           cbsri =  corresponding reaction coefficient
c
          nsrt = nsrt + 1
c
          if (nsrt .gt. nsrtmx) then
            write (noutpt,1600) nsrtmx
            write (nttyo,1600) nsrtmx
 1600       format(/' * Error - (EQ6/rd6inw) Have too many special',
     $      /7x,'reactants. The code is only dimensioned for ',i4,
     $      /7x,'such reactants. Reduce the number of such reactants',
     $      /7x,'or increase the dimensioning parameter nsrtpa.')
            go to 990
          endif
c
          read (ninpts,1610,err=990) vreac(nrc)
 1610     format(12x,e12.5)
          write (noutpt,1620) vreac(nrc)
 1620     format(6x,'vreac= ',1pe12.5)
c
          write (noutpt,1622)
 1622     format (' * Elemental composition')
c
          nci = 0
          do n = 1,nctmax + 1
            read (ninpts,1630,err=990) uxe,xx
 1630       format(3x,a8,3x,e22.15)
            if (uxe(1:8) .eq. uendit(1:8)) then
              j3 = ilnobl(uxe)
              write (noutpt,1640) uxe(1:j3)
 1640         format(4x,a)
              go to 160
            endif
            write (noutpt,1650) uxe,xx
 1650       format(4x,a8,3x,1pe22.15)
c
            nci = nci + 1
c
            if (nci .gt. nctmax) then
              write (noutpt,1660) ureac(nrc)(1:j2),nctmax
              write (nttyo,1660) ureac(nrc)(1:j2),nctmax
 1660         format(/' * Error - (EQ6/rd6inw) Have too many',
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
          write (noutpt,1662)
 1662     format (' * Reaction')
c
c         Note: the special reactant itself should appear first in
c         the associated chemical reaction. Ordinarily, its reaction
c         coefficient should be -1.0.
c
          nbi = 0
          do n = 1,nbt1mx + 1
            read (ninpts,1670,err=990) uxs,xx
 1670       format(3x,a24,3x,e22.15)
            if (uxs(1:8) .eq. uendit(1:8)) then
              j3 = ilnobl(uxs)
              write (noutpt,1680) uxs(1:j3)
 1680         format(4x,a)
              go to 170
            endif
            write (noutpt,1690) uxs,xx
 1690       format(4x,a24,3x,1pe22.15)
c
            nbi = nbi + 1
c
            if (nbi .gt. nbt1mx) then
              write (noutpt,1700) ureac(nrc)(1:j2),nbtmax
              write (nttyo,1700) ureac(nrc)(1:j2),nbtmax
 1700         format(/' * Error - (EQ6/rd6inw) Have too many basis',
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
            write (noutpt,1710) nertmx
            write (nttyo,1710) nertmx
 1710       format(/' * Error - (EQ6/rd6inw) Have too many generic',
     $      ' ion exchanger',/7x,'reactants. The code is only',
     $      ' dimensioned for ',i4,' such reactants.',/7x,'Reduce',
     $      ' the number of such reactants or increase the',
     $      /7x,'dimensioning parameter nertpa.')
            go to 990
          endif
c
          read (ninpts,1712,err=990) ugermo(ner)
 1712     format(12x,a24)
          j3 = ilnobl(ugermo(ner))
          write (noutpt,1714) ugermo(ner)(1:j3)
 1714     format(5x,'ugermo= ',a)
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
              write (noutpt,1718) ureac(nrc)(1:j2),jetmax
              write (nttyo,1718) ureac(nrc)(1:j2),jetmax
 1718         format(/' * Error - (EQ6/rd6inw) Have too many',
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
                write (noutpt,1724) ureac(nrc)(1:j2),
     $          ugerji(jei,ner)(1:j3),netmax
                write (nttyo,1724) ureac(nrc)(1:j2),
     $          ugerji(jei,ner)(1:j3),netmax
 1724           format(/' * Error - (EQ6/rd6inw) Have too many',
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
c         nsk    = surface area flag
c           0 = fixed surface area
c           1 = fixed specific surface area
c           2 = n**2/3 growth law- current surface area
c         sfcar  = surface area, cm2
c         ssfcar = specific surface area, cm2/mol
c         fkrc   = ratio of effective surface area to total surface
c                    area of a reactant; this parameter can also be
c                    used as a generalized correction factor for the
c                    rate constant
c
        read (ninpts,1740,err=990) nsk(nrc),sfcar(nrc),ssfcar(nrc),
     $  fkrc(nrc)
 1740   format(12x,i2,22x,e12.5,12x,e12.5,/12x,e12.5)
        write (noutpt,1750) nsk(nrc),sfcar(nrc),fkrc(nrc)
 1750   format(8x,'nsk= ',i2,15x,'sfcar= ',1pe12.5,4x,'ssfcar= ',e12.5,
     $  /7x,'fkrc= ',e12.5)
c
c
c       Rate law codes: nrk(1,nrc) = forward rate law code,
c       nrk(2,nrc) = backward rate law code.
c
c         Forward rate law codes:
c
c           -1 = Use the backward rate law form (legal only if
c                  nrk(2,nrc) = 2)
c            0 = Illegal value
c            1 = Relative rate
c            2 = Transition state theory net rate
c            3 = Linear rate
c
c         For the case nrk(1,nrc) = 2, the reactant must be
c         either a pure mineral (jcode = 0) or a solid solution
c         (jcode = 1).
c
c         Backward rate law codes:
c
c           -1 = Use the forward rate law form (legal only if
c                  nrk(1,nrc) = 2)
c            0 = No rate law specified; the reaction may be controlled
c                  by partial equilibrium
c            1 = Relative rate
c            2 = Transition state theory net rate
c            3 = Linear rate
c
c         For the case nrk(2,nrc) = 2, the reactant must be
c         either a pure mineral (jcode = 0) or a solid solution
c         (jcode = 1)
c
        read (ninpts,1760,err=990) nrk(1,nrc),nrk(2,nrc)
 1760   format(12x,i2,22x,i2)
        write (noutpt,1770) nrk(1,nrc),nrk(2,nrc)
 1770   format(7x,'nrk1= ',i2,16x,'nrk2= ',i2)
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
          write (noutpt,1775) ureac(nrc)(1:j2)
          write (nttyo,1775) ureac(nrc)(1:j2)
 1775     format(/' * Error - (EQ6/rd6inw) The forward rate law code',
     $    /7x,'has an illegal value of 0 for ',a,'.')
          go to 990
        elseif (nrk(1,nrc) .eq. 1) then
c
c         Arbitrary kinetics (relative rates, indifferent to time).
c
c           relative rate = rkb(1,1,nrc) + rkb(2,1,nrc)*xi1
c                             + (rkb(3,1,nrc)/2.)*xi1**2
c
c         where
c
c           rkb    = relative rate constants
c           xi1    = the reaction progress variable
c
          imech(1,nrc) = 3
c
          if (imech(1,nrc) .gt. imchmx) then
            write (noutpt,1780) ureac(nrc)(1:j2),imchmx
            write (nttyo,1780) ureac(nrc)(1:j2),imchmx
 1780       format(/' * Error - (EQ6/rd6inw) Have too many rate',
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
          write (noutpt,1800) (rkb(i,1,nrc), i = 1,3)
 1800     format(7x,'rkb1= ',1pe12.5,6x,'rkb2= ',e12.5,
     $    6x,'rkb3= ',e12.5)
        elseif (nrk(1,nrc) .eq. 2) then
c
c         Transition state rate law. Up to imchmx parallel mechanisms
c         are allowed.
c
c           imech  = number of parallel transition state mechanisms
c           rkb    = rate constants, one for each rate law term
c           ndact  = number of aqueous species that appear in each
c                      rate law term
c           udac   = the name of such an aqueous species
c           cdac   = the power the activity of such a species is raised
c                      to in the rate law term
c           csigma = ratio of affinity of a macroscopic reaction to the
c                      affinity of the microcopic reaction for
c                      destruction of an activated complex
c
          read (ninpts,1810,err=990) imech(1,nrc)
 1810     format(12x,i2)
          write (noutpt,1820) imech(1,nrc)
 1820     format(6x,'imech= ',i2)
c
          if (imech(1,nrc) .gt. imchmx) then
            write (noutpt,1780) ureac(nrc)(1:j2),imchmx
            write (nttyo,1780) ureac(nrc)(1:j2),imchmx
            go to 990
          endif
c
          do i = 1,imech(1,nrc)
            read (ninpts,1830,err=990) rkb(i,1,nrc),trkb(i,1,nrc),
     $      iact(i,1,nrc)
 1830       format(12x,e12.5,12x,e12.5,12x,i2)
            write (noutpt,1840) rkb(i,1,nrc),trkb(i,1,nrc),
     $      iact(i,1,nrc)
 1840       format(8x,'rkb= ',1pe12.5,6x,'trkb= ',e12.5,
     $      6x,'iact= ',i2)
c
            read (ninpts,1850,err=990) eact(i,1,nrc),hact(i,1,nrc)
 1850       format(12x,e12.5,12x,e12.5)
            write (noutpt,1860) eact(i,1,nrc),hact(i,1,nrc)
 1860       format(7x,'eact= ',1pe12.5,6x,'hact= ',e12.5)
c
            read (ninpts,1870,err=990) ndact(i,1,nrc),csigma(i,1,nrc)
 1870       format(12x,i2,22x,e12.5)
            write (noutpt,1880) ndact(i,1,nrc),csigma(i,1,nrc)
 1880       format(6x,'ndact= ',i2,14x,'csigma= ',1pe12.5)
c
            if (ndact(i,1,nrc) .gt. ndctmx) then
              write (noutpt,1890) i,ureac(nrc)(1:j2),ndctmx
              write (nttyo,1890) i,ureac(nrc)(1:j2),ndctmx
 1890         format(/' * Error - (EQ6/rd6inw) Have too many',
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
              write (noutpt,1920) udac(n,i,1,nrc),cdac(n,i,1,nrc)
 1920         format(7x,'udac= ',a24,6x,'cdac= ',1pe12.5)
            enddo
          enddo
        elseif (nrk(1,nrc) .eq. 3) then
c
c         Linear rate law.
c
c           Rate = fkrc(nrc)*sfcar(nrc)*rkb(1,1,nrc)
c
          imech(1,nrc) = 1
c
          if (imech(1,nrc) .gt. imchmx) then
            write (noutpt,1780) ureac(nrc)(1:j2),imchmx
            write (nttyo,1780) ureac(nrc)(1:j2),imchmx
            go to 990
          endif
c
          i = 1
          read (ninpts,1830,err=990) rkb(i,1,nrc),trkb(i,1,nrc),
     $    iact(i,1,nrc)
          write (noutpt,1840) rkb(i,1,nrc),trkb(i,1,nrc),
     $    iact(i,1,nrc)
          read (ninpts,1850,err=990) eact(i,1,nrc),hact(i,1,nrc)
          write (noutpt,1860) eact(i,1,nrc),hact(i,1,nrc)
        else
c
          write (noutpt,1950) nrk(1,nrc),ureac(nrc)(1:j2)
          write (nttyo,1950) nrk(1,nrc),ureac(nrc)(1:j2)
 1950     format(/' * Error - (EQ6/rd6inw) The forward rate law code',
     $    /7x,'has an unrecognized value of ',i2,' for ',a,'.')
          go to 990
c
        endif
c
c       Read backward (precipitation, formation) rate data.
c
c         If nrk(2,nrc) = 0, then the backward may be governed by
c         instantaneous partial equilibrium.
c
c         If nrk(2,nrc) = -1, then the forward rate law will be used.
c         This is legal only for the cases of nrk(1,nrc) = 2
c         (transition state theory net rate).
c
c         If nrk(2,nrc) <= 0, don't enter any additional backward
c         rate law data.
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
            write (noutpt,1960) ureac(nrc)(1:j2),imchmx
            write (nttyo,1960) ureac(nrc)(1:j2),imchmx
 1960       format(/' * Error - (EQ6/rd6inw) Have too many rate',
     $      /7x,'constants or corresponding mechanisms in the backward',
     $      /7x,'rate law for reactant ',a,'. The code is only',
     $      /7x,'dimensioned for ',i2,' rate constants per rate law.',
     $      /7x,'Reduce the number of rate constants or increase the',
     $      /7x,'dimensioning parameter imchpa.')
            go to 990
          endif
c
          read (ninpts,1790,err=990) (rkb(i,2,nrc), i = 1,3)
          write (noutpt,1800) (rkb(i,2,nrc), i = 1,3)
        elseif (nrk(2,nrc) .eq. 2) then
c
c         Transition state rate law.
c
          read (ninpts,1810,err=990) imech(2,nrc)
          write (noutpt,1820) imech(2,nrc)
c
          if (imech(2,nrc) .gt. imchmx) then
            write (noutpt,1960) ureac(nrc)(1:j2),imchmx
            write (nttyo,1960) ureac(nrc)(1:j2),imchmx
            go to 990
          endif
c
          do i = 1,imech(2,nrc)
            read (ninpts,1830,err=990) rkb(i,2,nrc),trkb(i,2,nrc),
     $      iact(i,2,nrc)
            write (noutpt,1840) rkb(i,2,nrc),trkb(i,2,nrc),
     $      iact(i,2,nrc)
c
            read (ninpts,1850,err=990) eact(i,2,nrc),hact(i,2,nrc)
            write (noutpt,1860) eact(i,2,nrc),hact(i,2,nrc)
c
            read (ninpts,1870,err=990) ndact(i,2,nrc),csigma(i,2,nrc)
            write (noutpt,1880) ndact(i,2,nrc),csigma(i,2,nrc)
c
            if (ndact(i,2,nrc) .gt. ndctmx) then
              write (noutpt,1970) i,ureac(nrc)(1:j2),ndctmx
              write (nttyo,1970) i,ureac(nrc)(1:j2),ndctmx
 1970         format(/' * Error - (EQ6/rd6inw) Have too many',
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
              write (noutpt,1920) udac(n,i,2,nrc),cdac(n,i,2,nrc)
            enddo
          enddo
        elseif (nrk(2,nrc) .eq. 3) then
c
c         Linear rate law.
c
          imech(2,nrc) = 1
c
          if (imech(2,nrc) .gt. imchmx) then
            write (noutpt,1960) ureac(nrc)(1:j2),imchmx
            write (nttyo,1960) ureac(nrc)(1:j2),imchmx
            go to 990
          endif
c
          i = 1
          read (ninpts,1830,err=990) rkb(i,2,nrc),trkb(i,2,nrc),
     $    iact(i,2,nrc)
          write (noutpt,1840) rkb(i,2,nrc),trkb(i,2,nrc),
     $    iact(i,2,nrc)
          read (ninpts,1850,err=990) eact(i,2,nrc),hact(i,2,nrc)
          write (noutpt,1860) eact(i,2,nrc),hact(i,2,nrc)
        else
c
          write (noutpt,1980) nrk(2,nrc),ureac(nrc)(1:j2)
          write (nttyo,1980) nrk(2,nrc),ureac(nrc)(1:j2)
 1980     format(/' * Error - (EQ6/rd6inw) The backward rate law code',
     $    /7x,'has an unrecognized value of ',i2,' for ',a,'.')
          go to 990
c
        endif
c
        write (noutpt,1500)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Starting, minimum, and maximum values of key run parameters.
c
c       xistti  = starting value of reaction progress for this run
c       ximaxi  = maximum value of Xi for this run
c       tistti  = time at Xi = xistti (sec)
c       timmxi  = maximum time (sec) for this run
c
      read (ninpts,1150,err=990) xistti,ximaxi,tistti,timmxi
 1150 format(2(12x,e12.5),/2(12x,e12.5))
      write (noutpt,1155) xistti,ximaxi,tistti,timmxi
 1155 format(5x,'xistti= ',1pe12.5,4x,'ximaxi= ',e12.5,/
     $ 5x,'tistti= ',1pe12.5,4x,'timmxi= ',e12.5)
c
c       phmini  = minimum pH for this run
c       phmaxi  = maximum pH for this run
c       ehmini  = minimum Eh (v) for this run
c       ehmaxi  = maximum Eh (v) for this run
c       o2mini  = minimum log fO2 for this run
c       o2maxi  = maximum log fO2 for this run
c       awmini  = minimum aw for this run
c       awmaxi  = maximum aw for this run
c       kstpmx  = maximum number of steps for this run
c
      read (ninpts,1160,err=990) phmini,phmaxi,ehmini,ehmaxi,o2mini,
     $ o2maxi,awmini,awmaxi,kstpmx
 1160 format(2(12x,e12.5),/2(12x,e12.5),/2(12x,e12.5),/2(12x,e12.5),
     $ /12x,i12)
      write (noutpt,1165) phmini,phmaxi,ehmini,ehmaxi,o2mini,o2maxi,
     $ awmini,awmaxi,kstpmx
 1165 format(5x,'phmini= ',1pe12.5,4x,'phmaxi= ',e12.5,/
     $ 5x,'ehmini= ',1pe12.5,4x,'ehmaxi= ',e12.5,/
     $ 5x,'o2mini= ',1pe12.5,4x,'o2maxi= ',e12.5,/
     $ 5x,'awmini= ',1pe12.5,4x,'awmaxi= ',e12.5,/
     $ 5x,'kstpmx= ',i12)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print interval parameters.
c
c       dlxprn = print increment in terms of Xi
c       dlxprl = log Xi print interval
c       dltprn = print increment in terms of time (seconds)
c       dltprl = log time print interval
c
      read (ninpts,1170,err=990) dlxprn,dlxprl,dltprn,dltprl
 1170 format(2(12x,e12.5),12x,/2(12x,e12.5))
      write (noutpt,1175) dlxprn,dlxprl,dltprn,dltprl
 1175 format(5x,'dlxprn= ',1pe12.5,4x,'dlxprl= ',e12.5,
     $ /5x,'dltprn= ',e12.5,4x,'dltprl= ',e12.5)
c
c       dlhprn = print increment in terms of pH
c       dleprn = print increment in terms of Eh (v)
c       dloprn = print increment in terms of log fO2
c       dlaprn = print increment in terms of aw
c       ksppmx = limit on number of steps from the last print point
c                  at which a print will be made
c
      read (ninpts,1180,err=990) dlhprn,dleprn,dloprn,dlaprn,ksppmx
 1180 format(2(12x,e12.5),12x,/2(12x,e12.5),/12x,i12)
      write (noutpt,1185) dlhprn,dleprn,dloprn,dlaprn,ksppmx
 1185 format(5x,'dlhprn= ',1pe12.5,4x,'dleprn= ',e12.5,
     $ /5x,'dloprn= ',e12.5,4x,'dlaprn= ',e12.5,/5x,'ksppmx= ',i12)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Plot interval parameters.
c
c       dlxplo = plot increment in terms of Xi
c       dlxpll = log Xi plot interval
c       dltplo = plot increment in terms of time (seconds)
c       dltpll = log time plot interval
c
      read (ninpts,1190,err=990) dlxplo,dlxpll,dltplo,dltpll
 1190 format(2(12x,e12.5),/2(12x,e12.5))
      write (noutpt,1192) dlxplo,dlxpll,dltplo,dltpll
 1192 format(5x,'dlxplo= ',1pe12.5,4x,'dlxpll= ',e12.5,
     $ /5x,'dltplo= ',e12.5,4x,'dltpll= ',e12.5)
c
c       dlhplo = plot increment in terms of pH
c       dleplo = plot increment in terms of Eh (v)
c       dloplo = plot increment in terms of log fO2
c       dlaplo = plot increment in terms of aw
c       ksplmx = limit on number of steps from the last plot point
c                  at which a plot will be made
c
      read (ninpts,1195,err=990) dlhplo,dleplo,dloplo,dlaplo,ksplmx
 1195 format(2(12x,e12.5),/2(12x,e12.5),/12x,i12)
      write (noutpt,1197) dlhplo,dleplo,dloplo,dlaplo,ksplmx
 1197 format(5x,'dlhplo= ',1pe12.5,4x,'dleplo= ',e12.5,
     $ /5x,'dloplo= ',e12.5,4x,'dlaplo= ',e12.5,/5x,'ksplmx= ',i12)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iopt model option switches.
c     Note: iopt(1) = iopt1, etc.
c
c       iopt(1) = Switch to choose the physical model:
c          0 = Closed system
c          1 = Titration system
c          2 = Fluid-centered flow-through open system
c
c       iopt(2) = Kinetic mode:
c          0 = Reaction progress mode (arbitrary kinetics)
c          1 = Reaction progress/time mode (true kinetics)
c
c       iopt(3) = Phase boundary searches:
c          0 =  The step size is constrained by the locations of
c                 the predicted phase boundaries (iopt(3) is reset
c                 to 0 if other options require this)
c          1 =  The location of phase boundaries is estimated from
c                 Taylor's series and printed, but the step size is
c                 unconstrained
c          2  = The locations of phase boundaries are ignored
c
c       iopt(4) = Solid solution products:
c          0 = Solid solutions are ignored
c          1 = Solid solutions are permitted
c
c       iopt(5) = Clear the ES solids read from the input file:
c          0 = Don't do it
c          1 = Do it
c
c       iopt(6) = Clear the ES solids at the initial value of reaction
c                   progress; this is done after calculating the state
c                   of the system at this point:
c          0 = Don't do it
c          1 = Do it
c
c       iopt(7) = Clear the ES solids at the end of the run:
c          0 = Don't do it
c          1 = Do it, unless the run terminates early due to
c                calculational difficulties
c
c       iopt(8) = Not used
c
c       iopt(9) = Clear the PRS solids read from the input file:
c          0 = Don't do it
c          1 = Do it
c
c       iopt(10) = Clear the PRS solids at the end of the run:
c          0 = Don't do it
c          1 = Do it, unless the run terminates early due to
c                calculational difficulties
c
c       iopt(11) = Auto basis switching, in pre-N-R optimization:
c          0 = Off
c          1 = On
c
c       iopt(12) = Auto basis switching, after N-R iteration:
c          0 = Off
c          1 = On
c
c       iopt(13) = Switch to choose calculational mode:
c          0 = Normal path tracing
c          1 = Economy mode (if permissible)
c          2 = Super economy mode (if permissible)
c
c       iopt(14) = ODE integrator corrector mode:
c          0 = Allow stiff and simple correctors
c          1 = Allow only the simple corrector
c          2 = Allow only the stiff corrector
c          3 = Allow no correctors
c
c       iopt(15) = Force suppression of all redox reactions:
c          0 = Don't do it
c          1 = Do it (this is potentially dangerous)
c
c       iopt(16) = Backup file options:
c         -1 = Don't write backup files
c          0 = Write backup files
c          1 = Write a sequential backup file
c
c       iopt(17) = Pickup file options:
c         -1 = Don't write a pickup file
c          0 = Write a pickup file
c
c       iopt(18) = Tab file options:
c         -1 = Don't write the tab file
c          0 = Write the tab file
c          1 = The tab file output is appended to the tabx file from a
c                 previous run; if more than one problem is stacked on
c                 the input file, this option applies to only the first
c                 problem, for any subsequent any subsequent problems,
c                 iopt(18) is reset to 0
c
c       iopt(19) = Used only by EQ3NR
c
c       iopt(20) = Advanced EQ6 pickup file options:
c          0 = Write a normal EQ6 pickup file
c          1 = Write an EQ6 input file with "Fluid 2" set up as a
c                special reactant for mixing with another fluid or
c                reaction with a more generalized equilibrium system
c
      write (noutpt,1200)
 1200 format(' *',15x,'1    2    3    4    5    6    7    8    9   10')
c
      read (ninpts,1210,err=990) (iopt(i), i = 1,20)
 1210 format(12x,10i5)
      write (noutpt,1220) (iopt(i), i = 1,20)
 1220 format(3x,'iopt1-10= ',10i5,/2x,'iopt11-20= ',10i5)
c
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iopr print option switches.
c     Note: iopr(1) = iopr1, etc.
c
c       iopr(1) = List the names of all species read from the
c                   supporting data file:
c          0 = Don't list
c          1 = List
c
c       iopr(2) = List all reactions:
c          0 = Don't list
c          1 = List all reactions (this can be quite lengthy)
c          2 = Also print the log K values
c          3 = Also print the coefficients of the interpolating
c                polynomials
c
c       iopr(3) = List the hard core diameters of the aqueous species:
c          0 = Don't list
c          1 = List
c
c       iopr(4) = Print a table at each print point of the
c                   concentrations, activities, and activity
c                   coefficients of the aqueous species:
c         -3  = Omit species with molalities < 1.e-8
c         -2 =  Omit species with molalities < 1.e-12
c         -1 =  Omit species with molalities < 1.e-20
c          0 =  Omit species with molalities < 1.e-100
c          1 =  Include all species
c
c       iopr(5) = Print a table at each print point of the cation/H+
c                   activity ratios, anion-H+ activity products, and
c                   neutral species activities:
c          0 = Don't print
c          1 = Print cation/H+ activity ratios only
c          2 = Print cation/H+ activity ratios and anion-H+ activity
c                products only
c          3 = Print cation/H+ activity ratios, anion-H+ activity
c                products, and neutral species activities
c
c       iopr(6) = At each print point, print a table of the percentage
c                    contributions for each aqueous mass balance total:
c         -1 = Don't print any tables
c          0 = Print tables including 99% of all contributing species
c          1 = Print tables including all contributing species
c
c       iopr(7) = Print a table at each print point of the saturation
c                   indices and affinities of the various non-aqueous
c                   phases:
c         -1 = Don't print
c          0 = Print for those phases not undersaturated by
c                more than 10 kcal
c          1 = Print for all phases
c
c       iopr(8) = Print a table at each print point of the fugacities
c                   of the gas species:
c         -1 = Don't print
c          0 = Print
c          1 = Print
c
c       iopr(9) = Print a table at each print of the mean molal
c                   activity coefficients:
c         -1 = Don't print
c          0 = Don't print
c          1 = Print
c
c       iopr(10) = Print a tabulation at the start of running the
c                    current problem of the Pitzer interaction
c                    coefficients:
c          0 = Don't print
c          1 = Print a summary of the names of the species present and
c                the number of Pitzer interaction coefficients
c          2 = Print a summary of the names of the species present and
c                the number of Pitzer interaction coefficients
c
c       iopr(11) - iopr(16) = Not used
c
c       iopr(17) = Pickup file format:
c          0 = Use the same format ("D" or "W") as the input file
c          1 = Use "W" format
c          2 = Use "D" format
c
c       iopr(18) - iopr(20) = Not used
c
      read (ninpts,1230,err=990) (iopr(i), i = 1,20)
 1230 format(12x,10i5)
      write (noutpt,1240) (iopr(i), i = 1,20)
 1240 format(3x,'iopr1-10= ',10i5,/2x,'iopr11-20= ',10i5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iodb debugging print option switches.
c     Note: iodb(1) = iodb1, etc.
c
c       iodb(1) = General diagnostic messages:
c          0 = Don't print
c          1 = Print Level 1 diagnostic messages
c          2 = Print Level 1 and Level 2 diagnostic messages
c
c       iodb(2) = Kinetics-related diagnostic messages:
c          0 = Don't print
c          1 = Print Level 1 kinetics diagnostic messages
c          2 = Print Level 1 and Level 2 kinetics diagnostic messages
c
c       iodb(3) = Pre-Newton-Raphson optimization:
c          0 = Don't print
c          1 = Print summary information
c          2 = Print detailed information
c          3 = Print more detailed information
c          4 = Also print changes to activity coefficients
c
c       iodb(4) = Information describing Newton-Raphson iteration:
c          0 = Don't print
c          1 = Print summary information
c          2 = Print detailed information including the residual and
c                correction vectors
c          3 = Also print the Jacobian matrix
c          4 = Also print changes to activity coefficients
c
c       iodb(5) = Step size/order selection:
c          0 = Don't print
c          1 = Print the chosen scale factor
c          2 = Print the orders under consideration and their respective
c                step size scaling factors
c
c       iodb(6) = Hypothetical affinity calculations (for solid
c                   solutions, etc.):
c          0 = Don't print
c          1 = Print summary information
c          2 = Print detailed information
c
c       iodb(7) = Search iterations (for phase boundaries, etc.):
c          0 = Don't print
c          1 = Print summary information
c
c       iodb(8) = Information describing ODE corrector iteration:
c          0 = Don't print
c          1 = Print summary information
c          2 = Print detailed information
c
c       iodb(9) - iodb(20) = Not used
c
      read (ninpts,1250,err=990) (iodb(i), i = 1,20)
 1250 format(12x,10i5)
      write (noutpt,1260) (iodb(i), i = 1,20)
 1260 format(3x,'iodb1-10= ',10i5,/2x,'iodb11-20= ',10i5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Number of nxopt options.
c
c       nxopt = number of mineral subset-selection suppression options
c
      read (ninpts,1300,err=990) nxopt
 1300 format(12x,i2)
      write (noutpt,1310) nxopt
 1310 format(6x,'nxopt= ',i2)
c
      if (nxopt .gt. nxopmx) then
        write (noutpt,1320) nxopmx
        write (nttyo,1320) nxopmx
 1320   format(/' * Error - (EQ6/rd6inw) Have too many mineral',
     $  /7x,'subset-selection suppression options. The code is',
     $  /7x,'only dimensioned for ',i3,' such options. Reduce the',
     $  /7x,'number of options or increase the dimensioning',
     $  /7x,'parameter nxoppa.')
        go to 990
      endif
c
c     Nxopt options.
c
c       uxopt  = 'all' or 'alwith'
c       uxcat  = name of a chemical element
c       nxopex = number of exceptions
c       uxopex = name of the minerals that are the exceptions
c
      do n = 1,nxopt
        read (ninpts,1330,err=990) uxopt(n),uxcat(n)
 1330   format(12x,a6,1x,a24)
        j2 = ilnobl(uxcat(n))
        write (noutpt,1340) uxopt(n),uxcat(n)(1:j2)
 1340   format(5x,'option= ',a6,1x,a)
      enddo
c
      if (nxopt .gt. 0) then
        read (ninpts,1350,err=990) nxopex
 1350   format(12x,i2)
        write (noutpt,1360) nxopex
 1360   format(5x,'nxopex= ',i2)
c
        if (nxopex .gt. nxpemx) then
          write (noutpt,1370) nxpemx
          write (nttyo,1370) nxpemx
 1370     format(/' * Error - (EQ6/rd6inw) Have too many',
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
          j2 = ilnobl(uxopex(n))
          write (noutpt,1390) uxopex(n)(1:j2)
 1390     format(2x,'exception= ',a)
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Number of nffg options.
c
c       nffg   = number of gases with fixed fugacities
c
      read (ninpts,1400,err=990) nffg
 1400 format(12x,i2)
      write (noutpt,1410) nffg
 1410 format(7x,'nffg= ',i2)
c
      if (nffg .gt. nffgmx) then
        write (noutpt,1420) nffgmx
        write (nttyo,1420) nffgmx
 1420   format(/' * Error - (EQ6/rd6inw) Have too many gases whose',
     $  /7x,'fugacities are to be fixed. The code is only dimensioned',
     $  /7x,'for ',i4,' such gases. Reduce the number of gases or',
     $  /7x,'increase the dimensioning parameter nffgpa.')
        go to 990
      endif
c
c     Nffg options.
c
c       uffg   = name of species
c       moffg  = moles of gas species
c       xlkffg = log fugacity
c
      do n = 1,nffg
        read (ninpts,1430,err=990) uffg(n),moffg(n),xlkffg(n)
 1430   format(12x,a24,/12x,e12.5,12x,e12.5)
        write (noutpt,1440) uffg(n),moffg(n),xlkffg(n)
 1440   format(4x,'species= ',a24,/6x,'moffg= ',1pe12.5,4x,
     $  'xlkffg= ',e12.5)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Under normal circumstances, users should enter zeros and take
c     the code default values for dlxdmp, tolbt, toldl, tolxsf, tolsat,
c     dlxmx0, itermx, and ntrymx.
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
c       tolbt = convergence tolerance on residual magnitude
c                 (in Newton-Raphson iteration)
c       toldl = convergence tolerance on correction magnitude
c                 (in Newton-Raphson iteration)
c
      read (ninpts,2020,err=990) tolbt,toldl
 2020 format(2(12x,e12.5))
      write (noutpt,2030) tolbt,toldl
 2030 format(6x,'tolbt= ',1pe12.5,5x,'toldl= ',e12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum number of Newton-Raphson iterations.
c
c       itermx = limit on the number of Newton-Raphson iterations.
c                  Recommended values lie in the range 50 to 200.
c
      read (ninpts,2100,err=990) itermx
 2100 format(12x,i3)
      write (noutpt,2110) itermx
 2110 format(5x,'itermx= ',i3)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Search/find tolerance.
c
c       tolxsf = search/find tolerance
c
      read (ninpts,2022,err=990) tolxsf
 2022 format(12x,e12.5)
      write (noutpt,2032) tolxsf
 2032 format(5x,'tolxsf= ', 1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Saturation tolerance.
c
c       tolsat = supersaturation tolerance below which no attempt
c                  is made to precipitate a phase. A good value is
c                  0.0001 kcal. Don't make it too small.
c
      read (ninpts,2024,err=990) tolsat
 2024 format(12x,e12.5)
      write (noutpt,2034) tolsat
 2034 format(5x,'tolxsf= ', 1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Maximum number of phase assemblage tries.
c
c       ntrymx = limit on the number of attempted phase assemblages
c                  at a point of reaction progress. Recommended
c                  values lie in the range of 50 to 100.
c
      read (ninpts,2102,err=990) ntrymx
 2102 format(12x,i3)
      write (noutpt,2112) ntrymx
 2112 format(5x,'ntrymx= ',i3)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Zero-order step size (in Xi).
c
c       dlxmx0 = the normal step size for order zero; this is not the
c                  minimum step size (dlxmin). Recommended values lie
c                  in the range 1.e-10 to 1.e-6.
c
      read (ninpts,2080,err=990) dlxmx0
 2080 format(12x,e12.5)
      write (noutpt,2090) dlxmx0
 2090 format(5x,'dlxmx0= ',1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     PRS transfer interval in Xi.
c
c       dlxdmp = fixed reaction progress interval for the transfer of
c                  phases from the equilibrium system (ES) to the
c                  physically removed system (PRS); this only works
c                  in conjunction with the fluid-centered flow-through
c                  open system model (iopt(1) = 2)
c
      read (ninpts,2000,err=990) dlxdmp
 2000 format(12x,e12.5)
      write (noutpt,2005) dlxdmp
 2005 format(5x,'dlxdmp= ',1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1500)
c
c     Process the bottom half of the current input file.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Secondary title.
c
c     utitl2 = secondary title
c     ntitl2 = number of lines in the secondary title.
c
      ntitl2 = 0
      read (ninpts,1000,err=990) uline
c
      j2 = ilnobl(uline)
      j2 = min(j2,79)
      write (noutpt,1020) uline
c
      if (uline(1:8) .ne. uendit(1:8)) then
        ntitl2 = 1
        utitl2(1) = uline
        do n = 2,ntitmx
          read (ninpts,1000,err=990) uline
          j2 = ilnobl(uline)
          j2 = min(j2,79)
          write (noutpt,1020) uline(1:j2)
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
c       nsbswt = the number of special basis switches
c       usbsw(1,n) = the species to be switched from the strict
c                      basis to the auxiliary basis set in the
c                      in the n-th switch
c       usbsw(2,n) = the species to be switched from the auxiliary
c                      basis set to the strict basis set in the
c                      n-th switch
c
      write (noutpt,2610)
 2610 format(' *   Special basis switches')
c
      read (ninpts,2620,err=990) nsbswt
 2620 format(12x,i3)
      write (noutpt,2630) nsbswt
 2630 format(5x,'nsbswt= ',i3)
c
      do n = 1,nsbswt
        read (ninpts,2640,err=990) usbsw(1,n)
 2640   format(9x,a48)
        j2 = ilnobl(usbsw(1,n))
        write (noutpt,2650) usbsw(1,n)(1:j2)
 2650   format(1x,'species= ',a)
c
        read (ninpts,2660,err=990) usbsw(2,n)
 2660   format(15x,a48)
        j2 = ilnobl(usbsw(2,n))
        write (noutpt,2670) usbsw(2,n)(1:j2)
 2670   format(3x,'switch with= ',a)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Original temperature.
c
c       tempci = temperature, C, at the end of the previous
c                  EQ3NR or EQ6 run
c
      read (ninpts,2680,err=990) tempci
 2680 format(12x,e12.5)
      write (noutpt,2690) tempci
 2690 format(5x,'tempci= ',1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Original pressure.
c
c       pressi = pressure, bars, at the end of the previous
c                  EQ3NR or EQ6 run
c
      read (ninpts,2700,err=990) pressi
 2700 format(12x,e12.5)
      write (noutpt,2710) pressi
 2710 format(5x,'pressi= ',1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Ion exchanger creation. This section reads the directives for
c     creating ion exchange phases and species and their associated
c     intrinsic (including thermodynamic) properties. The data is
c     essentially that which might be read from a supporting data file.
c
c       qgexsh = logical flag. If .true., force the writing of at
c                  least one exchanger block on a "D" format input
c                  file (when running XCON6) or pickup file.
c       net    = the number of ion exchangers
c       ugexp  = array of ion exchange phase names
c       mwtges = array of molecular weights of substrates of ion
c                  exchange phases
c       ugexmo = array of strings identifying models for composing
c                  exact ion exchange species and corresponding
c                  reactions; examples include:
c                    'Vanselow' = Vanselow model
c                    'Gapon'    = Gapon model
c                    'Site-mixing' = the general site-mixing model
c       tgexp  = array of reference temperatures (C)  for the
c                  thermodynamic data for the reactions specified
c                  for the exchanger phases
c       jgext  = array of numbers of exchange sites in ion exchange
c                  phases
c       ugexj  = array of names specified for the sites on the various
c                 exchangers
c       cgexj  = array of numbers of moles of exchange sites per mole
c                  of substrate in ion exchange phases
c       zgexj  = array of electrical charges of exchange sites
c                  in ion exchange phases. The total charge number
c                  for a site is the product of this charge and
c                  and the number of moles of the site per mole
c                  of exchanger. The sign of the charge on a site
c                  is generally the opposite of that of the ions
c                  which exchange on it.
c       ngexrt = array of numbers of ion exchange reactions specified
c                  for a given site and exchanger
c       ugexr  = array of strings containing compact representations
c                  of the exchange reactions; e.g., 'Na+ = Ca++' for a
c                  reaction in which Na+ on the exchanger is replaced
c                  by Ca++. One may make specifications such as
c                  'Na+ = ' in which case the ion goes into solution
c                  leaving a bare substrate. All reactions are
c                  normalized to the exchange (or loss) of one
c                  equivalent. The exact form of the reaction is
c                  otherwise dependent on the mixing law specifed in
c                  the element of the ugexmo array for the current
c                  exchanger.
c       xlkgex = array of equilibrium constants or related data for
c                  the reaction specified in the above string
c       uxkgex = array of strings denoting the kind of data in the
c                  corresponding entry of the xlkgex array:
c                    ' '       = log K per equivalent
c                    'LogK/eq' = log K per equivalent
c                    'kcal/eq' = DeltaG0r, kcal, per equivalent
c                    'kJ/eq'   = DeltaG0r, kJ, per equivalent
c       xhfgex = array of enthalpy of reaction data for the reaction
c                  specified in the above string
c       uhfgex = array of strings denoting the kind of data in the
c                  corresponding entry of the xhfgex array:
c                    ' '       = DeltaG0r, kcal, per equivalent
c                    'kcal/eq' = DeltaG0r, kcal, per equivalent
c                    'kJ/eq'   = DeltaG0r, kJ, per equivalent
c       xvfgex = array of volume of reaction data for the reaction
c                  specified in the above string
c       uvfgex = array of strings denoting the kind of data in the
c                  corresponding entry of the xhfgex array:
c                    ' '      = cm3 per equivalent
c                    'cm3/eq' = cm3 per equivalent
c
      write (noutpt,2790)
 2790 format(' * Ion exchanger creation')
c
      read (ninpts,2795) ux8
 2795 format(12x,a8)
      call lejust(ux8)
      call locase(ux8)
      qgexsh = ux8(1:2).eq.'.t' .or. ux8(1:1).eq. 't'
      write (noutpt,2797) qgexsh
 2797 format(5x,'qgexsh= ',l8)
c
      read (ninpts,2800,err=990) net
 2800 format(12x,i3)
      write (noutpt,2810) net
 2810 format(8x,'net= ',i3)
c
      if (net .gt. netmax) then
        write (noutpt,2820) netmax,net
        write (nttyo,2820) netmax,net
 2820   format(/' * Error - (EQ6/rd6inw) Have exceeded the maximum',
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
        write (noutpt,2840) ugexp(ne)(1:j2)
 2840   format(6x,'ugexp= ',a)
c
        read (ninpts,2850,err=990) mwtges(ne)
 2850   format(12x,e12.5)
        write (noutpt,2860) mwtges(ne)
 2860   format(5x,'mwtges= ',1pe12.5)
c
        read (ninpts,2830,err=990) ugexmo(ne)
        j3 = ilnobl(ugexmo(ne))
        write (noutpt,2870) ugexmo(ne)(1:j3)
 2870   format(5x,'ugexmo= ',a)
c
        read (ninpts,2850,err=990) tgexp(ne)
        write (noutpt,2880) tgexp(ne)
 2880   format(6x,'tgexp= ',1pe12.5)
c
        read (ninpts,2890,err=990) jgext(ne)
 2890   format(12x,i3)
        write (noutpt,2900) jgext(ne)
 2900   format(6x,'jgext= ',i3)
c
        if (jgext(ne) .gt. jetmax) then
          write (noutpt,2910) jetmax,ugexp(ne)(1:j2),jgext(ne)
          write (nttyo,2910) jetmax,ugexp(ne)(1:j2),jgext(ne)
 2910     format(/' * Error - (EQ6/rd6inw) Have exceeded the maximum',
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
          write (noutpt,2940) ugexj(je,ne)(1:j3)
 2940     format(6x,'ugexj= ',a)
c
          read (ninpts,2950,err=990) cgexj(je,ne),zgexj(je,ne)
 2950     format(12x,e12.5,12x,e12.5)
          write (noutpt,2960) cgexj(je,ne),zgexj(je,ne)
 2960     format(6x,'cgexj= ',1pe12.5,5x,'zgexj= ',e12.5)
c
          read (ninpts,2890,err=990) ngexrt(je,ne)
          write (noutpt,3010) ngexrt(je,ne)
 3010     format(5x,'ngexrt= ',i3)
c
          if (ngexrt(je,ne) .gt. ietmax) then
            write (noutpt,3020) netmax,ugexj(je,ne)(1:j3),
     $      ugexp(ne)(1:j2),ngexrt(je,ne)
            write (nttyo,3020) netmax,ugexj(je,ne)(1:j3),
     $      ugexp(ne)(1:j2),ngexrt(je,ne)
 3020       format(/' * Error - (EQ6/rd6inw) Have exceeded the',
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
            j2 = ilnobl(ugexr(n,je,ne))
            write (noutpt,3040) ugexr(n,je,ne)(1:j2)
 3040       format(6x,'ugexr= ',a)
            read (ninpts,3050,err=990) xlkgex(n,je,ne),uxkgex(n,je,ne)
 3050       format(12x,e12.5,12x,a8)
            j4 = ilnobl(uxkgex(n,je,ne))
            write (noutpt,3060) xlkgex(n,je,ne),uxkgex(n,je,ne)(1:j4)
 3060       format(5x,'xlkgex= ',1pe12.5,5x,'units= ',a)
            read (ninpts,3050,err=990) xhfgex(n,je,ne),uhfgex(n,je,ne)
            j4 = ilnobl(uhfgex(n,je,ne))
            write (noutpt,3070) xhfgex(n,je,ne),uhfgex(n,je,ne)(1:j4)
 3070       format(5x,'xhfgex= ',1pe12.5,5x,'units= ',a)
            read (ninpts,3050,err=990) xvfgex(n,je,ne),uvfgex(n,je,ne)
            j4 = ilnobl(uvfgex(n,je,ne))
            write (noutpt,3080) xvfgex(n,je,ne),uvfgex(n,je,ne)(1:j4)
 3080       format(5x,'xvfgex= ',1pe12.5,5x,'units= ',a)
          enddo
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Number of nxmod options.
c
c       nxmod = the number of suppressed/altered species/reactions
c
      read (ninpts,3120,err=990) nxmod
 3120 format(12x,i2)
      write (noutpt,3130) nxmod
 3130 format(6x,'nxmod= ',i2)
c
      if (nxmod .gt. nxmdmx) then
        write (noutpt,3140) nxmdmx
        write (nttyo,3140) nxmdmx
 3140   format(/' * Error - (EQ6/rd6inw) Have too many nxmod',
     $  /7x,'alter/suppress options. The code is only dimensioned',
     $  /7x,'for ',i2,' such options. Reduce the number of such',
     $  /7x,'options or increase the dimensioning parameter nxmdpa.')
        go to 990
      endif
c
c     Nxmod options.
c
c       uxmod = the name (all 48 letters) of the species. If the
c                 phase part is not given, the option is applied to
c                 every species for which the species part of its name
c                 generates a match.
c       kxmod  = alter/suppress code
c         -1 = the species is suppressed
c          0 = the log K is replaced by xlkmod
c          1 = the log K is augmented by xlkmod
c          2 = same as kxmod=1, but xlkmod is input in units of kcal
c                per mole of the associated species
c       xlkmod = log K value alteration function as defined above
c
      do n = 1,nxmod
        read (ninpts,3150,err=990) uxmod(n),kxmod(n),xlkmod(n)
 3150   format(12x,a48,/12x,i2,22x,e12.5)
        j2 = ilnobl(uxmod(n))
        write (noutpt,3160) uxmod(n)(1:j2),kxmod(n),xlkmod(n)
 3160   format(4x,'species= ',a,
     $  /5x,'option= ',i2,14x,'xlkmod= ',1pe12.5)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iopg options.
c     Note: iopg(1) = iopg1, etc.
c
c       iopg(1) = Model for aqueous species activity coefficients:
c         -1 = The  Davies equation
c          0 = The B-dot equation
c          1 = Pitzer's  equations
c          2 = HC + DH equations
c
c       iopg(2) = Rescaling of aqueous ionic activity coefficients for
c                   consistency with a desired pH scale:
c         -1 = "Internal" pH scale; no scaling (e.g., the "Davies"
c                scale if iopg(1) = -1, the "B-dot" scale if
c                iopg(1) = 0, or the Pitzer scale if iopg(1) = 1)
c          0 = The NBS pH scale (log gamma(Cl-) is defined by the
c                Bates-Guggenheim equation)
c          1 = The Mesmer pH scale (log gamma(H+) = 0)
c
c       iopg(3) = iopg(10) = Not used
c
      read (ninpts,3170,err=990) (iopg(i), i = 1,20)
 3170 format(12x,10i5)
      write (noutpt,1200)
      write (noutpt,3180) (iopg(i), i = 1,20)
 3180 format(3x,'iopg1-10= ',10i5,/2x,'iopg11-20= ',10i5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Index limits.
c
c       kct    = the number of elements in the matrix
c       kbt    = the number of basis species in the matrix
c       kmt    = position of last pure mineral in the matrix
c       kxt    = position of last solid solution in the matrix
c       kdim   = size of the matrix
c       kprs   = flag to read input for initializing the physically
c                  removed system (PRS)
c
      read (ninpts,3220,err=990) kct,kbt,kmt,kxt,kdim,kprs
 3220 format(3(12x,i2,10x),/3(12x,i2,10x))
      write (noutpt,3230) kct,kbt,kmt,kxt,kdim,kprs
 3230 format(8x,'kct= ',i2,17x,'kbt= ',i2,17x,'kmt= ',i2,/
     $ 8x,'kxt= ',i2,16x,'kdim= ',i2,16x,'kprs= ',i2)
c
      nbti = kbt
c
      if (kct .gt. nctmax) then
        write (noutpt,3240) nctmax
        write (nttyo,3240) nctmax
 3240   format(/' * Error - (EQ6/rd6inw) Have too many chemical',
     $  /7x,'elements present. The code is only dimensioned',
     $  /7x,'for ',i3,' elements. Reduce the number of elements',
     $  /7x,'or increase the dimensioning parameter nctpar.')
        go to 990
      endif
c
      if (kbt .gt. nbtmax) then
        write (noutpt,3250) nbtmax
        write (nttyo,3250) nbtmax
 3250   format(/' * Error - (EQ6/rd6inw) Have too many basis',
     $  /7x,'species present. The code is only dimensioned',
     $  /7x,'for ',i3,' basis species. Reduce the number of elements',
     $  /7x,'or increase the dimensioning parameter nctpar.')
        go to 990
      endif
c
      if (kdim .gt. kmax) then
        write (noutpt,3260) kmax
        write (nttyo,3260) kmax
 3260   format(/' * Error - (EQ6/rd6inw) Have too many matrix',
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
      write (noutpt,3330)
 3330 format(' * Data file basis species and jflag values')
c
      do nbi = 1,nbti
        read (ninpts,3340,err=990) ubmtbi(nbi),jflgi(nbi)
 3340   format(3x,a48,3x,i2)
        write (noutpt,3350) ubmtbi(nbi),jflgi(nbi)
 3350   format(4x,a48,3x,i2)
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
c       electr = electrical imbalance
c
      write (noutpt,3400)
 3400 format(' *',7x,
     $ 'Mass balance totals     Aqueous mass balance totals')
c
      do nbi = 1,nbti
        read (ninpts,3410,err=990) mtbi(nbi),mtbaqi(nbi)
 3410   format(6x,e22.15,6x,e22.15)
        write (noutpt,3420) mtbi(nbi),mtbaqi(nbi)
 3420   format(7x,1pe22.15,6x,e22.15)
      enddo
      read (ninpts,3430,err=990) electr
 3430 format(34x,e22.15)
      write (noutpt,3440) electr
 3440 format(10x,'Electrical imbalance= ',3x,1pe22.15)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Ordinary basis switches.
c
c       nobswt = the number of ordinary basis switches
c       uobsw(1,n) = the species to be switched from the basis
c                      set in the n-th switch
c       uobsw(2,n) = the species to be switched into the basis
c                      set in the n-th switch
c
      write (noutpt,3450)
 3450 format(' *   Ordinary basis switches')
c
      read (ninpts,2620,err=990) nobswt
      write (noutpt,3460) nobswt
 3460 format(5x,'nobswt= ',i3)
c
      do n = 1,nobswt
        read (ninpts,2640,err=990) uobsw(1,n)
        j2 = ilnobl(uobsw(1,n))
        write (noutpt,2650) uobsw(1,n)(1:j2)
c
        read (ninpts,2660,err=990) uobsw(2,n)
        j2 = ilnobl(uobsw(2,n))
        write (noutpt,2670) uobsw(2,n)(1:j2)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Matrix column variables.
c
c       uzveci = names of species or properties associated with the
c                  matrix column variables
c
      write (noutpt,3470)
 3470 format(' * Matrix species or entities')
c
      do krow = 1,kdim
        read (ninpts,3480,err=990) uzveci(krow)
 3480   format(3x,a48)
        j2 = ilnobl(uzveci(krow))
        write (noutpt,3490) uzveci(krow)(1:j2)
 3490   format(4x,a)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Matrix variable values.
c
c       zvclg1 = values of matrix variables
c
      write (noutpt,3510)
 3510 format(' *   Values of matrix variables')
      do kcol = 1,kdim
        read (ninpts,3520,err=990) zvclgi(kcol)
 3520   format(3x,e22.15)
        write (noutpt,3530) zvclgi(kcol)
 3530   format(4x,1pe22.15)
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
        write (noutpt,3600)
 3600   format(' *  Phases and species in the PRS')
        j3 = ilnobl(uendit)
        do npi = 1,nprpmx + 1
          read (ninpts,3610,err=990) uxf,mx
 3610     format(1x,a24,6x,e22.15)
          if (uxf(1:8) .eq. uendit(1:8)) then
            write (noutpt,3620) uendit(1:j3)
 3620       format(2x,a)
            go to 310
          else
            write (noutpt,3630) uxf,mx
 3630       format(2x,a24,6x,1pe22.15)
          endif
          nprpti = nprpti + 1
c
          if (nprpti .gt. nprpmx) then
            write (noutpt,3640) nprpmx
            write (nttyo,3640) nprpmx
 3640       format(/' * Error - (EQ6/rd6inw) Have too many phases',
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
            if (uxg(1:8) .eq. uendit(1:8)) then
              write (noutpt,3660) uendit(1:j3)
 3660         format(4x,a)
              go to 300
            else
              write (noutpt,3670) uxg,mx
 3670         format(4x,a24,6x,1pe22.15)
            endif
c
            if (iki .gt. iktmax) then
              j2 = ilnobl(uxf)
              write (noutpt,3680) uxf(1:j2),iktmax
              write (nttyo,3680) uxf(1:j2),iktmax
 3680         format(/' * Error - (EQ6/rd6inw) Have too many',
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
              write (noutpt,3690) nprsmx
              write (nttyo,3690) nprsmx
 3690         format(/' * Error - (EQ6/rd6inw) Have too many species',
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
      write (noutpt,4000) nprob
      write (nttyo,4000) nprob
 4000 format(/'   Done reading problem ',i3,'.',/)
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
