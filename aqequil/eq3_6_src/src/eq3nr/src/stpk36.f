      subroutine stpk36(awmaxi,awmini,cbsri,cbsr1,cdac,cesri,cesr1,
     $ csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,
     $ dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,
     $ dlxplo,dlxprl,dlxprn,eact,ehmaxi,ehmini,fkrc,hact,iact,
     $ ibsrti,ibsrt1,iesrti,iesrt1,iktmax,imchmx,imech,iopt,itermx,
     $ ixrti,jcode,jpress,jreac,jtemp,ksplmx,ksppmx,kstpmx,modr,
     $ moffg,morr,mprphi,mprspi,nbt1mx,nctmax,ndact,ndctmx,nffg,
     $ nffgmx,noptmx,nordmx,noutpt,nprob,nprpmx,nprpti,nprsmx,
     $ nprsti,nptkmx,nrct,nrctmx,nrk,nsk,nsrt,nsrtmx,ntitl1,ntitmx,
     $ ntrymx,nttkmx,nttyo,nxopex,nxopmx,nxopt,nxpemx,nxrt,nxrtmx,
     $ o2maxi,o2mini,phmaxi,phmini,ptk,press,pressb,rkb,rxbari,
     $ sfcar,ssfcar,tempc,tempcb,tempc1,timmxi,tistti,tolsat,
     $ tolxsf,trkb,ttk,ubsri,ubsr1,ucxri,udac,uesri,uesr1,uffg,
     $ uprphi,uprspi,ureac,ureac1,utitl1,uxcat,uxopex,uxopt,vreac,
     $ ximaxi,xistti,xlkffg)
c
c     This subroutine sets up certain variables and arrays for writing
c     the top half of an EQ6 input file. There are currently three
c     options for the form and content of this part of this file.
c     (see below). EQ3NR/setpk3.f performs a somewhat similiar function
c     for data to be written on the normal EQ3NR pickup file, or the
c     bottom half of an EQ6 input file.
c
c     This subroutine should only be called if iopt(19) is greater than
c     zero.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
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
      integer iktmax,imchmx,nbt1mx,nctmax,ndctmx,nffgmx,noptmx,nordmx,
     $ nprpmx,nprsmx,nptkmx,nsrtmx,ntitmx,nttkmx,nrctmx,nxopmx,nxpemx,
     $ nxrtmx
c
      integer noutpt,nttyo
c
      integer iact(imchmx,2,nrctmx),ibsrti(nsrtmx),iesrti(nsrtmx),
     $ imech(2,nrctmx),iopt(noptmx),ixrti(nxrtmx),jcode(nrctmx),
     $ jreac(nrctmx),ndact(imchmx,2,nrctmx),nrk(2,nrctmx),nsk(nrctmx)
c
      integer ibsrt1,iesrt1,itermx,jpress,jtemp,ksplmx,ksppmx,
     $ kstpmx,nffg,nprob,nprpti,nprsti,nrct,nsrt,ntitl1,ntrymx,
     $ nxopex,nxopt,nxrt
c
      character(len=80) utitl1(ntitmx)
      character(len=48) uprspi(nprsmx)
      character(len=24) ubsri(nbt1mx,nsrtmx),ubsr1(nbt1mx),
     $ ucxri(iktmax,nxrtmx),udac(ndctmx,imchmx,2,nrctmx),uffg(nffgmx),
     $ uprphi(nprpmx),ureac(nrctmx),uxcat(nxopmx),uxopex(nxpemx)
      character(len=8) uesri(nctmax,nsrtmx),uesr1(nctmax),uxopt(nxopmx)
c
      character(len=24) ureac1
c
      real(8) cbsri(nbt1mx,nsrtmx),cbsr1(nbt1mx),
     $ cdac(ndctmx,imchmx,2,nrctmx),cesri(nctmax,nsrtmx),cesr1(nctmax),
     $ csigma(imchmx,2,nrctmx),eact(imchmx,2,nrctmx),fkrc(nrctmx),
     $ hact(imchmx,2,nrctmx),modr(nrctmx),moffg(nffgmx),morr(nrctmx),
     $ mprphi(nprpmx),mprspi(nprsmx),ptk(nptkmx),rkb(imchmx,2,nrctmx),
     $ rxbari(iktmax,nxrtmx),sfcar(nrctmx),ssfcar(nrctmx),
     $ trkb(imchmx,2,nrctmx),ttk(nttkmx),vreac(nrctmx),xlkffg(nffgmx)
c
      real(8) awmaxi,awmini,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,
     $ dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,
     $ dlxplo,dlxprl,dlxprn,ehmaxi,ehmini,o2maxi,o2mini,phmaxi,phmini,
     $ press,pressb,tempc,tempcb,tempc1,timmxi,tistti,tolsat,tolxsf,
     $ ximaxi,xistti
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,n,nmax
c
      integer ilnobl
c
      character(len=8) ux8
c
c-----------------------------------------------------------------------
c
      if (iopt(19) .le. 0) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set up the initial part of the new main title.
c
      call initcb(utitl1,ntitmx)
      ntitl1 = 25
      utitl1(1) = 'EQ6 input file name= sample.6i'
      utitl1(2) = 'Description= "Sample"'
      utitl1(3) = 'Version level= 8.0'
      utitl1(4) = 'Revised mm/dd/yy    Revisor= Username'
      utitl1(6) = '  This is a sample EQ6 input file written as an EQ3NR
     $ pickup file.'
      utitl1(7) = 'The EQ3NR input file used to generate this output is
     $identified in'
      utitl1(8) = 'the second title given below.'
      utitl1(10) = '  You are expected to modify this EQ6 input file to
     $meet your own needs.'
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set various EQ6 input file variables.
c
      ksplmx = 0
      ksppmx = 0
      kstpmx = 500
      ntrymx = 0
      itermx = 0
c
      xistti = 0.
      ximaxi = 1.0
c
      tistti = 0.
      timmxi = 1.e+38
c
      awmaxi = 1.e+38
      awmini = -1.e+38
      ehmaxi = 1.e+38
      ehmini = -1.e+38
      o2maxi = 1.e+38
      o2mini = -1.e+38
      phmaxi = 1.e+38
      phmini = -1.e+38
c
      dlxdmp = 0.
      dlxmx0 = 0.
c
      dlxprn = 0.
      dlxprl = 0.5
      dlxplo = 0.
      dlxpll = 0.
c
      dltprn = 1.e+38
      dltprl = 1.e+38
      dltplo = 1.e+38
      dltpll = 1.e+38
c
      dlaprn = 1.e+38
      dleprn = 1.e+38
      dlhprn = 1.e+38
      dloprn = 1.e+38
      dlaplo = 1.e+38
      dleplo = 1.e+38
      dlhplo = 1.e+38
      dloplo = 1.e+38
c
      tolsat = 0.
      tolxsf = 0.
      nordmx = 6
c
      jtemp = 0
      call initaz(ttk,nttkmx)
      tempcb = tempc
c
      if (iopt(19) .eq. 3) then
        if (abs(tempc - tempc1) .gt. 1.e-6) then
c
c         Set up for non-isothermal fluid-mixing.
c
          jtemp = 3
          ttk(2) = tempc1
          ttk(1) = 1.0
        endif
      endif
c
      jpress = 0
      call initaz(ptk,nptkmx)
      pressb = press
c
      nrct = 0
      call initiz(jcode,nrctmx)
      call initiz(jreac,nrctmx)
      call initiz(nsk,nrctmx)
c
      nsrt = 0
      call initiz(ibsrti,nsrtmx)
      call initiz(iesrti,nsrtmx)
c
      nmax = nbt1mx*nsrtmx
      call initaz(cbsri,nmax)
      call initcb(ubsri,nmax)
c
      nmax = nctmax*nsrtmx
      call initaz(cesri,nmax)
      call initcb(uesri,nmax)
c
      call initiz(ixrti,nxrtmx)
c
      nxrt = 0
      nmax = iktmax*nxrtmx
      call initaz(rxbari,nmax)
      call initcb(ucxri,nmax)
c
      call initcb(ureac,nrctmx)
      call initaz(fkrc,nrctmx)
      call initaz(sfcar,nrctmx)
      call initaz(ssfcar,nrctmx)
c
      nmax = 2*nrctmx
      call initiz(imech,nmax)
      call initiz(nrk,nmax)
c
      nmax = imchmx*2*nrctmx
      call initaz(csigma,nmax)
      call initaz(rkb,nmax)
      call initaz(trkb,nmax)
      call initaz(eact,nmax)
      call initaz(hact,nmax)
      call initiz(iact,nmax)
      call initiz(ndact,nmax)
c
      nmax = ndctmx*imchmx*2*nrctmx
      call initaz(cdac,nmax)
      call initcb(udac,nmax)
c
      nffg = 0
      call initcb(uffg,nffgmx)
      call initaz(moffg,nffgmx)
      call initaz(xlkffg,nffgmx)
c
      nxopt = 0
      call initcb(uxopt,nxopmx)
      call initcb(uxcat,nxopmx)
      nxopex = 0
      call initcb(uxopex,nxpemx)
c
      nprpti = 0
      call initaz(mprphi,nprpmx)
      call initcb(uprphi,nprpmx)
c
      nprsti = 0
      call initaz(mprspi,nprsmx)
      call initcb(uprspi,nprsmx)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(19) .eq. 1) then
c
c       Write an EQ6 input file with a single dissolving reactant,
c       Quartz. Its dissolution is controlled by a specified
c       relative rate of 1.0. The basic purpose of this option is
c       to provide a template for a file on which the reactants
c       react according to specified relative rates. It is intended
c       that the user then modify this template to describe the
c       problem actually desired.
c
        utitl1(12) = '  This sample file has Quartz dissolving at a rela
     $tive rate of 1.0.'
        utitl1(13) = 'It may or may not actually run. For example, Quart
     $z may not appear'
        utitl1(14) = 'on the supporting data file.'
c
        nrct = 1
        ureac(1) = 'Quartz'
        jcode(1) = 0
        jreac(1) = 0
        morr(1) = 1.0
        modr(1) = 0.
        nrk(1,1) = 1
        rkb(1,1,1) = 1.0
c
      elseif (iopt(19) .eq. 2) then
c
c       Write an EQ6 input file with a single dissolving reactant,
c       Albite. Its dissolution is controlled by a specified
c       TST rate law. The basic purpose of this option is to
c       provide a template for a file on which the reactants
c       react according to specified TST rate laws. It is intended
c       that the user then modify this template to describe the
c       problem actually desired.
c
        utitl1(12) = '  This sample file has Albite dissolving according
     $ to a TST rate law.'
        utitl1(13) = 'It may or may not actually run. For example, Albit
     $e may not appear'
        utitl1(14) = 'on the supporting data file.'
        utitl1(16) = '  The kinetic data shown here are taken from Knaus
     $s and Wolery (1986). You'
        utitl1(17) = 'may or may not wish to use these particular data.
     $In any case, you are'
        utitl1(18) = 'entirely responsible for the data that you do use.
     $'
        utitl1(20) = '                            References'
        utitl1(22) = 'Knauss, K.G., and Wolery, T.J., 1986, Dependence o
     $f albite dissolution'
        utitl1(23) = '  kinetics on pH and time at 25C and 70C, Geochimi
     $ca et Cosmochimica Acta,'
c
        utitl1(24) = '  v. 50, p. 2481-2497.'
        nrct = 1
        ureac(1) = 'Albite'
        jcode(1) = 0
        jreac(1) = 0
        morr(1) = 1.0
        modr(1) = 0.
        nsk(1) = 0
        sfcar(1) = 1.0e+4
        fkrc(1) = 1.0
c
        nrk(1,1) = 1
        nrk(2,1) = -1
        imech(1,1) = 3
c
c       Acid range.
c
        rkb(1,1,1) = 1.00e-15
        trkb(1,1,1) = 25.0
        iact(1,1,1) = 2
        hact(1,1,1) = 28.47
        csigma(1,1,1) = 1.0
        ndact(1,1,1) = 1
        udac(1,1,1,1) = 'H+'
        cdac(1,1,1,1) = 1.0
c
c       Neutral range.
c
        rkb(2,1,1) = 3.98e-17
        trkb(2,1,1) = 25.0
        iact(2,1,1) = 2
        hact(2,1,1) = 12.89
        csigma(2,1,1) = 1.0
        ndact(2,1,1) = 0
c
c       Alkaline range.
c
        rkb(3,1,1) = 5.01e-21
        trkb(3,1,1) = 25.0
        iact(3,1,1) = 2
        hact(3,1,1) = 7.68
        csigma(3,1,1) = 1.0
        ndact(3,1,1) = 1
        udac(1,3,1,1) = 'H+'
        cdac(1,3,1,1) = -0.5
c
      elseif (iopt(19) .eq. 3) then
c
c       Write an EQ6 input file for a fluid mixing calculation in
c       which "Fluid 2" is mixed into "Fluid 1". Note that "Fluid 2"
c       appears first on the EQ3NR input file, beflore "Fluid 1".
c       It appears as a special reactant ont he top half of the EQ6
c       input file. "Fluid 1" is the starting aqueous solution in the
c       equilibrium system, and appears in the bottom half. If this
c       option is chosen for the first and only problem on the EQ3NR
c       input file, the pickup file produced is set up for mixing of
c       this fluid with itself. This would normally be useful only
c       as a template.
c
c       Set iopt(17) for the first fluid ("Fluid 2") to -1, so no
c       pickup file is produced (iopt(19) can be set to 0). For the
c       second fluid ("Fluid 1"), set iopt(17) equal to 0 and
c       iopt(19) equal to 3.
c
        utitl1(12) = '  This sample file has "Fluid 2" mixing into "Flui
     $d 1". "Fluid 2" was the'
        utitl1(13) = 'first aqueous solution on the EQ3NR input file, "F
     $luid 1" was the second.'
        utitl1(14) = 'The option switch iopt(19) was set to 3 for the la
     $tter to produce the'
        utitl1(15) = 'present EQ6 input file.'
        if (jtemp .eq. 0) then
          utitl1(16) = ' '
          utitl1(17) = '  Mixing will be isothermal as specified below.'
        elseif (jtemp .eq. 3) then
          utitl1(16) = ' '
          utitl1(17) = '  Mixing will be non-isothermal as specified bel
     $ow.'
        endif
c
        if (nprob .eq. 1) then
          utitl1(18) = ' '
          utitl1(19) = '  Warning: This EQ6 input file is set up to mix
     $a fluid with itself.'
        endif
c
        nrct = 1
        nsrt = 1
        ureac(1) = ureac1
        jcode(1) = 2
        jreac(1) = 0
        morr(1) = 1.0
        modr(1) = 0.
        nrk(1,1) = 1
        rkb(1,1,1) = 1.0
c
        vreac(1) = 1.0
c
        iesrti(1) = iesrt1
        do n = 1,iesrt1
          uesri(n,1) = uesr1(n)
          cesri(n,1) = cesr1(n)
        enddo
c
        ibsrti(1) = ibsrt1
        do n = 1,ibsrt1
          ubsri(n,1) = ubsr1(n)
          cbsri(n,1) = cbsr1(n)
        enddo
c
      else
c
c       Unknown option switch value.
c
        write (ux8,"(i5)") iopt(19)
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1000) ux8(1:j2)
        write (nttyo,1000) ux8(1:j2)
 1000   format(/' * Error - (EQ3NR/stpk36) Programming error trap:',
     $  ' The "Advanced EQ3NR',/7x,'PICKUP File Options" switch',
     $  ' iopt(19) has an unknown value of ',a,'.')
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
