      subroutine echoz(axlks,awmaxi,awmini,azero,cbsr,cdac,cdrs,
     $ cesr,csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,
     $ dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmax,
     $ dlxmx0,dlxpll,dlxplo,dlxprl,dlxprn,eact,ehmaxi,ehmini,iktmax,
     $ imchmx,imech,iodb,iopg,iopr,iopt,itermx,jcode,jpress,jtemp,
     $ jsflag,ksplmx,ksppmx,kstpmx,kxmod,mwtrc,narn1,narn2,narxmx,
     $ nat,nata,natmax,nbaspd,nbt,nbta,nbtd,nbtmax,nbt1mx,ncmpr,nct,
     $ ncta,nctmax,ndact,ndctmx,ndrs,ndrsmx,ndrsr,nffgmx,ngt,ngta,
     $ ngtmax,nlt,nlta,nltmax,nmt,nmta,nmtmax,nodbmx,nopgmx,noprmx,
     $ noptmx,nordmx,noutpt,npslmx,npt,npta,nptkmx,nptmax,nrct,nrctmx,
     $ nrk,nrndex,nsk,nsrt,nsrtmx,nsscmx,nsslmx,nst,nsta,nstmax,ntprmx,
     $ ntrymx,nttkmx,nttyo,nxmdmx,nxmod,nxopmx,nxopex,nxopt,nxpemx,
     $ nxridx,nxrt,nxrtmx,nxt,nxta,nxtmax,o2maxi,o2mini,phmaxi,
     $ phmini,press,pressb,ptk,qredox,rkb,rxbar,sscrew,tempc,tempcb,
     $ tempk,timmxi,tistti,tolbt,toldl,tolsat,tolsst,tolxsf,tolxst,
     $ tolxsu,trkb,ttk,uactop,udac,uelem,uffg,ureac,uspec,uxcat,uxmod,
     $ uxopex,uxopt,vreac,xistti,ximaxi,xlkmod,xlks,zkfac,zklgmn,
     $ zklogl,zklogu)
c
c     This subroutine writes an echo of various parameters after any
c     default values or range corrections have been applied.
c
c     This subroutine is called by:
c
c       EQ6/eq6.f
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
      integer iktmax,imchmx,narxmx,natmax,nbtmax,nbt1mx,nctmax,ndctmx,
     $ ndrsmx,nffgmx,ngtmax,nltmax,nmtmax,nodbmx,nopgmx,noprmx,noptmx,
     $ nordmx,nptkmx,nptmax,nrctmx,nsrtmx,nsscmx,nstmax,ntprmx,nttkmx,
     $ nxmdmx,nxopmx,nxpemx,nxrtmx,nxtmax
c
      integer noutpt,nttyo
c
      integer imech(2,nrctmx),iodb(nodbmx),iopg(nopgmx),iopr(noprmx),
     $ iopt(noptmx),jcode(nrctmx),jsflag(nstmax),kxmod(nxmdmx),
     $ nbaspd(nbtmax),ncmpr(2,nptmax),ndact(imchmx,2,nrctmx),
     $ ndrs(ndrsmx),ndrsr(2,nstmax),nrk(2,nrctmx),nrndex(nrctmx),
     $ nsk(nrctmx),nxridx(nrctmx)
c
      integer itermx,jpress,jtemp,ksplmx,ksppmx,kstpmx,narn1,narn2,nat,
     $ nata,nbt,nbta,nbtd,nct,ncta,ngt,ngta,nlt,nlta,nmt,nmta,
     $ npslmx,npt,npta,nrct,nsrt,nsslmx,nst,nsta,ntrymx,nxmod,nxopex,
     $ nxopt,nxrt,nxt,nxta
c
      logical qredox
c
      character*48 uspec(nstmax),uxmod(nxmdmx)
      character*32 uactop
      character*24 udac(ndctmx,imchmx,2,nrctmx),uffg(nffgmx),
     $ ureac(nrctmx),uxcat(nxopmx),uxopex(nxpemx)
      character*8 uelem(nctmax),uxopt(nxopmx)
c
      real*8 axlks(narxmx,ntprmx,nstmax),azero(natmax),
     $ cbsr(nbt1mx,nsrtmx),cdac(ndctmx,imchmx,2,nrctmx),cdrs(ndrsmx),
     $ cesr(nctmax,nsrtmx),csigma(imchmx,2,nrctmx),
     $ eact(imchmx,2,nrctmx),mwtrc(nrctmx),ptk(nptkmx),
     $ rkb(imchmx,2,nrctmx),rxbar(iktmax,nxrtmx),sscrew(nsscmx),
     $ trkb(imchmx,2,nrctmx),ttk(nttkmx),vreac(nrctmx),xlkmod(nxmdmx),
     $ xlks(nstmax)
c
      real*8 awmaxi,awmini,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,
     $ dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmax,dlxmx0,
     $ dlxpll,dlxplo,dlxprl,dlxprn,ehmaxi,ehmini,o2maxi,o2mini,phmaxi,
     $ phmini,press,pressb,tempc,tempcb,tempk,timmxi,tistti,tolbt,
     $ toldl,tolsat,tolsst,tolxsf,tolxst,tolxsu,ximaxi,xistti,zkfac,
     $ zklgmn,zklogl,zklogu
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer nf
c
      integer i,ik,ilevel,j,j2,n,na,nb,nc,np,nrc,nr1,nr2,ns,nsr,nxr
c
      integer ilnobl
c
      real*8 cx,cxx
c
c-----------------------------------------------------------------------
c
      j2 = ilnobl(uactop)
      write (noutpt,1000) uactop(1:j2)
 1000 format(/' The activity coefficients of aqueous species will be',
     $ ' calculated using'/3x,a,'.'/)
c
      if (iopr(3) .gt. 1) then
        write (noutpt,1010)
 1010   format(/7x,'Species',17x,'HC Diameter',/)
        do ns = narn1 + 1,narn2
          na = ns - narn1 + 1
          write (noutpt,1020) uspec(ns),azero(na)
 1020     format(5x,a24,3x,f7.3)
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print a table showing for phases, species, and groups thereof,
c     the number of each of entity on the data base, the number the
c     software is dimensioned for, and the number appearing in the
c     current problem.
c
      call prtntt(nat,nata,natmax,nbt,nbta,nbtmax,nct,ncta,
     $ nctmax,ngt,ngta,ngtmax,nlt,nlta,nltmax,nmt,nmta,nmtmax,noutpt,
     $ npt,npta,nptmax,nst,nsta,nstmax,nxt,nxta,nxtmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (jtemp .eq. 0) then
        write (noutpt,1030) tempcb
 1030   format(/'   Temperature= ',f10.4,' C')
      elseif (jtemp .eq. 1) then
        write (noutpt,1040) tempcb,ttk(1)
 1040   format(/'   Temperature= ',f10.4,' + ',1pe12.5,' * Xi C')
      elseif (jtemp .eq. 2) then
        write (noutpt,1050) tempcb,ttk(1)
 1050   format(/'   Temperature= ',f10.4,' + ',1pe12.5,' * time C')
      elseif (jtemp .eq. 3) then
        write (noutpt,1060) tempcb,ttk(1),ttk(2),ttk(1)
 1060   format(/'   Temperature= (',f10.4,' * ',f10.4,' + Xi * ',
     $  f10.4,')',/20x,'/(Xi + ',f10.4,') C')
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (jpress .eq. 0) then
        write (noutpt,1100)
 1100   format(/'   Pressure= the data file reference curve value at',
     $  ' any temperature',/)
      elseif (jpress .eq. 1) then
        write (noutpt,1102)
 1102   format(/'   Pressure= the 1.013-bar/steam-saturation curve',
     $  ' value at any temperature',/)
      elseif (jpress .eq. 2) then
        write (noutpt,1110) pressb
 1110   format(/'   Pressure= ',1pg12.5,' bars',/)
      elseif (jpress .eq. 3) then
        write (noutpt,1120) pressb,ptk(1)
 1120   format(/'   Pressure= ',1pg12.5,' + (',g12.5,' * Xi) bars',/)
      elseif (jpress .eq. 4) then
        write (noutpt,1130) pressb,ptk(1)
 1130   format(/'   Pressure= ',1pg12.5,' + (',g12.5,' * time) bars',/)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopr(2) .ge. 1) then
        ilevel = 2
        nf = noutpt
        call echolk(axlks,cdrs,ilevel,jsflag,narxmx,ndrs,ndrsmx,
     $  ndrsr,nf,nst,ntprmx,nstmax,press,tempc,uspec,xlks)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1200) xistti,ximaxi,tistti,timmxi,phmini,phmaxi,
     $ ehmini,ehmaxi,o2mini,o2maxi,awmini,awmaxi,kstpmx
 1200 format(/' xistti= ',1pe12.5,' (Initial value of Xi)',
     $ /' ximaxi= ',e12.5,' (Maximum value of Xi)',
     $ /' tistti= ',e12.5,' (Initial value of time, sec)',
     $ /' timmxi= ',e12.5,' (Maximum value of time, sec)',
     $ /' phmini= ',e12.5,' (Minimum value of pH)',
     $ /' phmaxi= ',e12.5,' (Maximum value of pH)',
     $ /' ehmini= ',e12.5,' (Minimum value of Eh, v)',
     $ /' ehmaxi= ',e12.5,' (Maximum value of Eh, v)',
     $ /' o2mini= ',e12.5,' (Minimum value of log fO2)',
     $ /' o2maxi= ',e12.5,' (Maximum value of log fO2)',
     $ /' awmini= ',e12.5,' (Minimum value of aw)',
     $ /' awmaxi= ',e12.5,' (Maximum value of aw)',
     $ /' kstpmx= ',4x,i8,' (Maximum number of steps this run)',//)
      write (noutpt,1210) dlxprn,dlxprl,dltprn,dltprl,dlhprn,dleprn,
     $ dloprn,dlaprn,ksppmx
 1210 format('   dlxprn= ',1pe12.5,' (Print interval in Xi)',
     $ /'   dlxprl= ',e12.5,' (Print interval in log Xi)',
     $ /'   dltprn= ',e12.5,' (Print interval in time, sec)',
     $ /'   dltprl= ',e12.5,' (Print interval in log time)',
     $ /'   dlhprn= ',e12.5,' (Print interval in pH units)',
     $ /'   dleprn= ',e12.5,' (Print interval in Eh, v)',
     $ /'   dloprn= ',e12.5,' (Print interval in log fO2)',
     $ /'   dlaprn= ',e12.5,' (Print interval in aw)',
     $ /'   ksppmx= ',4x,i8,' (Print interval in steps)',//)
      write (noutpt,1220) dlxplo,dlxpll,dltplo,dltpll,dlhplo,dleplo,
     $ dloplo,dlaplo,ksplmx
 1220 format('     dlxplo= ',1pe12.5,' (Plot interval in Xi)',
     $ /'     dlxpll= ',e12.5,' (Plot interval in log Xi)',
     $ /'     dltplo= ',e12.5,' (Plot interval in time)',
     $ /'     dltpll= ',e12.5,' (Plot interval in log time)',
     $ /'     dlhplo= ',e12.5,' (Plot interval in pH units)',
     $ /'     dleplo= ',e12.5,' (Plot interval in Eh, v)',
     $ /'     dloplo= ',e12.5,' (Plot interval in log fO2)',
     $ /'     dlaplo= ',e12.5,' (Plot interval in aw)',
     $ /'     ksplmx= ',4x,i8,' (Plot interval in steps)',//)
      write (noutpt,1230) dlxdmp
 1230 format('   dlxdmp= ',1pe12.5,' (PRS transfer interval in Xi)',/)
      write (noutpt,1240) dlxmx0,dlxmax
 1240 format(/'     dlxmx0= ',1pe12.5,' (Zero-order step size in Xi)',
     $ /'     dlxmax= ',e12.5,' (Maximum step size)',/)
      write (noutpt,1250) nordmx
 1250 format(/' nordmx= ',i2,' (Maximum dimensioned order)',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1300) (iopt(i), i = 1,10)
 1300 format(/1x,'iopt(1)=  ',i2,' (Physical system model)',
     $ /1x,'iopt(2)=  ',i2,' (Kinetic mode)',
     $ /1x,'iopt(3)=  ',i2,' (Suppress phase boundary searches)',
     $ /1x,'iopt(4)=  ',i2,' (Solid solutions)',
     $ /1x,'iopt(5)=  ',i2,' (Clear ES solids read from the input',
     $ ' file)',
     $ /1x,'iopt(6)=  ',i2,' (Clear ES solids at the starting point)',
     $ /1x,'iopt(7)=  ',i2,' (Clear ES solids at the end of the run)',
     $ /1x,'iopt(8)=  ',i2,' (Not used)',
     $ /1x,'iopt(9)=  ',i2,' (Clear PRS solids read from the input',
     $ ' file)',
     $ /1x,'iopt(10)= ',i2,' (Clear PRS solids at the end of the run)')
      write (noutpt,1310) (iopt(i), i = 11,18)
 1310 format(1x,'iopt(11)= ',i2,' (Auto basis switching, in',
     $ ' pre-Newton-Raphson optimization)',
     $ /1x,'iopt(12)= ',i2,' (Auto basis switching, after',
     $ ' Newton-Raphson iteration)',
     $ /1x,'iopt(13)= ',i2,' (Calculational mode)',
     $ /1x,'iopt(14)= ',i2,' (ODE integrator corrector mode)',
     $ /1x,'iopt(15)= ',i2,' (Force global redox suppression)',
     $ /1x,'iopt(16)= ',i2,' (Backup file options)',
     $ /1x,'iopt(17)= ',i2,' (Pickup file options)',
     $ /1x,'iopt(18)= ',i2,' (Tab file options)',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1320) (iopg(i), i = 1,2)
 1320 format(/3x,'iopg(1)=  ',i2,' (Aqueous species activity',
     $ ' coefficient model)',
     $ /3x,'iopg(2)=  ',i2,' (pH scale)',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1330) (iopr(i), i = 1,11)
 1330 format(/1x,'iopr(1)=  ',i2,' (List all species)',
     $ /1x,'iopr(2)=  ',i2,' (List all reactions)',
     $ /1x,'iopr(3)=  ',i2,' (List HC diamaters)',
     $ /1x,'iopr 4)=  ',i2,' (Aqueous species concentration print',
     $ ' cut-off)',
     $ /1x,'iopr(5)=  ',i2,' (Ion/H+ activity ratios)',
     $ /1x,'iopr(6)=  ',i2,' (Mass balance percentages)',
     $ /1x,'iopr(7)=  ',i2,' (Affinity print cut-off)',
     $ /1x,'iopr(8)=  ',i2,' (Fugacities)',
     $ /1x,'iopr(9)=  ',i2,' (Mean molal activity coefficient)',
     $ /1x,'iopr(10)= ',i2,' (Pitzer coefficients tabulation)',
     $ /1x,'iopr(17)= ',i2,' (Pickup file format)',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1340) (iodb(i), i = 1,8)
 1340 format(/3x,'iodb(1)=  ',i2,' (General diagnostics)',
     $ /3x,'iodb(2)=  ',i2,' (Kinetics diagnostics)',
     $ /3x,'iodb(3)=  ',i2,' (Pre-Newton-Raphson optimization)',
     $ /3x,'iodb(4)=  ',i2,' (Newton-Raphson iterations)',
     $ /3x,'iodb(5)=  ',i2,' (Order/scaling calculations)',
     $ /3x,'iodb(6)=  ',i2,' (Hypothetical affinity iterations)',
     $ /3x,'iodb(7)=  ',i2,' (Search iterations)',
     $ /3x,'iodb(8)=  ',i2,' (ODE corrector iterations)',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1370) tolbt,toldl,tolxsf,tolxst,tolxsu,
     $ tolsat,tolsst
 1370 format(/' tolbt = ',1pe12.5,' (Residual function convergence',
     $ ' tolerance)',
     $ /' toldl = ',e12.5,' (Correction term convergence tolerance)',
     $ /' tolxsf= ',e12.5,' (Search/find tolerance (general,',
     $ ' relative))',
     $ /' tolxst= ',e12.5,' (Search/find tolerance on time (relative))',
     $ /' tolxsu= ',e12.5,' (Search/find tolerance on pH, Eh, etc.',
     $ ' (absolute))',
     $ /' tolsat= ',e12.5,' (Saturation tolerance)',
     $ /' tolsst= ',e12.5,' (Supersaturation tolerance)',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1380) (sscrew(i), i = 1,6)
 1380 format(/'   sscrew(1)= ',1pe10.3,
     $ ' (Matrix variable step size parameter)',
     $ /'   sscrew(2)= ',0pf10.5,' (Not used)',
     $ /'   sscrew(3)= ',1pe10.3,
     $ ' (Rate function step size parameter)',
     $ /'   sscrew(4)= ',1pe10.3,
     $ ' (Rate function corrector parameter)',
     $ /'   sscrew(5)= ',0pf10.5,
     $ ' (Under-relaxation parameter (Newton-Raphson))',
     $ /'   sscrew(6)= ',0pf10.5,' (economy mode step size)',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1400) zklogu,zklogl,zkfac,zklgmn
 1400 format(/' zklogu= ',f10.3,' (threshhold log mass for solids)',
     $ /' zklogl= ',f10.3,' (Log mass decrement for PRS shift)',
     $ /' zkfac = ',f10.3,' (Shift adjustment factor)',
     $ /' zklgmn= ',f10.3,' (Minimum log mass after a shift)',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1410) itermx,ntrymx,npslmx,nsslmx
 1410 format(/'   itermx= ',i3,' (Newton-Raphson iteration limit)',
     $ /'   ntrymx= ',i3,' (Phase assemblage try limit)',
     $ /'   npslmx= ',i3,' (Critical phase instability slide limit)',
     $ /'   nsslmx= ',i3,' (Critical redox instability slide limit)',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (.not.qredox) write (noutpt,1420)
 1420 format(/' There is no redox balance constraint; hence, the fO2,',
     $ /' Eh, pe, and Ah are all undefined.',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nsrt .gt. 0) then
        write (noutpt,1500)
 1500   format(//14x,'--- Special Reactants ---',/)
        do nrc = 1,nrct
          if (jcode(nrc) .eq. 2) then
            j2 = ilnobl(ureac(nrc))
            nsr = nrndex(nrc)
            write (noutpt,1510) ureac(nrc)(1:j2),mwtrc(nrc),vreac(nrc)
 1510       format(/3x,a,//5x,'Mol. wt. = ',f12.3,'g/mol',
     $      /5x,'Molar volume = ',f12.3,' cc',/)
c
            write (noutpt,1520)
 1520       format(/5x,'Elemental Composition:',
     $      //9x,'Element',5x,'mol/mol',/)
            do nc = 1,nct
              cx = cesr(nc,nsr)
              if (cx .ne. 0.) write (noutpt,1530) uelem(nc),cx
 1530         format(7x,a8,3x,1pe12.5)
            enddo
c
            write (noutpt,1540)
 1540       format(//5x,'Numerical Composition:',
     $      //9x,'Species',26x,'mol/mol',/)
            cxx = -cbsr(nbt1mx,nsr)
            do nb = 1,nbtd
              ns = nbaspd(nb)
              cx = cbsr(nb,nsr)/cxx
              if (cx .ne. 0.) write (noutpt,1550) uspec(ns),cx
 1550         format(7x,a24,3x,1pe22.15)
            enddo
c
            write (noutpt,1580)
 1580       format(//5x,'Reaction:')
c
            nf = noutpt
            call prrsr(cbsr,jcode,nbaspd,nbtd,nbtmax,nbt1mx,
     $      nf,noutpt,nrc,nrctmx,nrndex,nsrtmx,nstmax,nttyo,
     $      ureac,uspec)
c
            write (noutpt,1590)
 1590       format(1x)
          endif
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nxrt .gt. 0) then
        write (noutpt,1600)
 1600   format(/11x,'--- Solid Solution Reactants ---',/)
        do nrc = 1,nrct
          if (jcode(nrc) .eq. 1) then
            np = nrndex(nrc)
            j2 = ilnobl(ureac(nrc))
            write (noutpt,1610) ureac(nrc)(1:j2),mwtrc(nrc),vreac(nrc)
 1610       format(/3x,a,
     $      //5x,'Molecular weight = ',f12.3,' g/mol',
     $      /5x,'    Molar volume = ',f12.3,' cc',
     $      //9x,'Composition',16x,'Mole Fraction',/)
            nxr = nxridx(nrc)
            nr1 = ncmpr(1,np)
            nr2 = ncmpr(2,np)
            do ns = nr1,nr2
              ik = ns - nr1 + 1
              if (rxbar(ik,nxr) .ne. 0.)
     $        write (noutpt,1620) uspec(ns),rxbar(ik,nxr)
 1620         format(7x,a24,3x,1pe12.5)
            enddo
          endif
        enddo
        write (noutpt,1630)
 1630   format(1x)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c                 Rate laws and rate coefficients
c
      if (nrct .le. 0) go to 200
      write (noutpt,1700)
 1700 format(/21x,'--- Reactants/Rate Laws ---',/)
c
c     Forward direction.
c
      write (noutpt,1710)
 1710 format(/3x,'Forward direction')
c
      do nrc = 1,nrct
        if (nrk(1,nrc) .eq. -1) then
c
c         Use the backward net rate law form.
c
          write(noutpt,1720) ureac(nrc)
 1720     format(/1x,a24,3x,'Use backward net rate law ')
        elseif (nrk(1,nrc).eq.1) then
c
c         Specified relative rate.
c
          write (noutpt,1730) ureac(nrc),(rkb(j,1,nrc), j = 1,3)
 1730     format(/1x,a24,3x,'Specified relative rate',
     $    //23x,'rkb1= ',1pe12.5,/23x,'rkb2= ',e12.5,
     $    /23x,'rkb3= ',e12.5)
          write (noutpt,1740)
 1740   format(1x)
        elseif (nrk(1,nrc).eq.2) then
c
c         Transition state theory (TST) form.
c
          write (noutpt,1750) ureac(nrc)
 1750     format(/1x,a24,3x,'Transition state theory')
          do i = 1,imech(1,nrc)
            write (noutpt,1760) i,rkb(i,1,nrc),ndact(i,1,nrc),
     $      csigma(i,1,nrc)
 1760       format(/21x,'term= ',i2,
     $      //23x,'rkb= ',1pe12.5,3x,'ndact= ',i2,
     $      /23x,'csigma= ',1pe12.5)
            if (ndact(i,1,nrc).gt.0) then
              do n = 1,ndact(i,1,nrc)
                write (noutpt,1770) udac(n,i,1,nrc),cdac(n,i,1,nrc)
 1770           format(25x,'udac= ',a24,3x,'cdac= ',1pe12.5)
              enddo
            endif
          enddo
          write (noutpt,1740)
        elseif (nrk(1,nrc).eq.3) then
c
c         Linear form.
c
          write (noutpt,1780) ureac(nrc),rkb(1,1,nrc)
 1780     format(/1x,a24,3x,'Linear rate',//23x,'rkb= ',1pe12.5)
          write (noutpt,1740)
        endif
      enddo
c
c     Backward direction.
c
      write (noutpt,1810)
 1810 format(//3x,'Backward direction')
      do nrc = 1,nrct
        if (nrk(2,nrc).eq.-1) then
          write (noutpt,1820) ureac(nrc)
 1820     format(/1x,a24,3x,'Use forward net rate law ')
        elseif (nrk(2,nrc).eq.0) then
          write(noutpt,1830) ureac(nrc)
 1830     format(/1x,a24,3x,'Instantaneous equilibrium')
        elseif (nrk(2,nrc).eq.1) then
c
c         Specified relative rate.
c
          write (noutpt,1730) ureac(nrc),(rkb(j,2,nrc), j = 1,3)
        elseif (nrk(2,nrc).eq.2) then
c
c         Transition state theory (TST) form.
c
          write (noutpt,1750) ureac(nrc)
          do i = 1,imech(2,nrc)
            write (noutpt,1760) i,rkb(i,2,nrc),ndact(i,2,nrc),
     $      csigma(i,2,nrc)
            if (ndact(i,2,nrc).gt.0) then
              do n = 1,ndact(i,2,nrc)
                write (noutpt,1770) udac(n,i,2,nrc),cdac(n,i,2,nrc)
              enddo
            endif
          enddo
          write (noutpt,1740)
        elseif (nrk(2,nrc).eq.3) then
c
c         Linear form.
c
          write (noutpt,1780) ureac(nrc),rkb(1,2,nrc)
          write (noutpt,1740)
        endif
      enddo
c
  200 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1900)
 1900 format(/' - - - - - - - - - - - - - - - - - - - - - - - - ',
     $ '- - - - - - ',/)
      end
