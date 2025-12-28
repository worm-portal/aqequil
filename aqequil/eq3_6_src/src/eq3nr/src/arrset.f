      subroutine arrset(aamatr,abar,acflg,acflgo,act,actlg,adh,
     $ adhh,adhv,adhfs,adhfsx,advfs,advfsx,afcnst,alpha,al10,
     $ amtb,aphi,avcnst,azero,a3bar,a3bars,axhfs,axhfsx,axlks,
     $ axlksx,axvfs,axvfsx,bacfmx,bbig,bdh,bdhh,bdhv,bdot,bdoth,
     $ bdotv,beta,betamx,bfac,bgamx,bneg,bpx,bsigmm,bfje,bfxi,
     $ cco2,cdrs,cdrsx,cdrtw,cdrw,cegexs,cgexj,cjbasp,cnufac,
     $ conc,conclg,coval,cpgexs,csts,delam,delvec,dgpit,dhfs,
     $ dlogxw,dpelm,dpslm,dselm,dvfs,efac,egexjc,egexjf,egexs,
     $ eh,ehfac,elam,eps100,fje,fjeo,fo2,fo2lg,fsort,fugac,fugalg,
     $ fxi,fxio,gmmatr,gpit,ibetmx,ibpxt,ibswx,iction,iebal,ielam,
     $ iern1,iern2,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,
     $ ifrn1,ifrn2,ifzeta,igas,igstak,iindx1,ilcphi1,ilcphi2,
     $ ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,ilzeta,imrn1,imrn2,
     $ insgf,iopg,iodb,iopt,ipch,ipcv,ipivot,ipndx1,irdxc3,istack,
     $ ixbasp,ixrn1,ixrn2,izmax,jcsort,jern1,jern2,jflag,jgext,
     $ jgsort,jgstak,jjndex,jjsort,jpflag,jsflag,jsitex,jsol,
     $ jssort,jstack,ka1,kat,kbt,kct,kction,kdim,kebal,kelect,
     $ ker,ke1,ket,khydr,kkndex,km1,ko2gaq,kwater,kx1,kxt,loph,
     $ losp,lsort,mgext,moph,mosp,mrgexs,mtb,nalpha,napt,narn1,
     $ narn2,narxt,nbasp,nbaspd,nbaspx,nbt,nbtd,nbti,nbw,nchlor,
     $ ncmpr,ncosp,nct,ndecsp,ndrs,ndrsx,ndrsr,ndrsrd,ndrsrx,
     $ nelect,nern1,nern2,net,nfac,ngexsa,ngext,ngrn1,ngrn2,ngt,
     $ nhydr,nhydx,nmut,nmux,nmxi,nmxx,noutpt,no2gaq,nphasx,npt,
     $ nredox,nslt,nslx,nst,nsts,nstsr,nsxi,nsxx,ntfx,ntfxt,ntpr,
     $ nttyo,omega,omeglg,palpha,pe,pelm,pmu,presg,press,pslamn,
     $ pslm,qbassw,qchlor,qhawep,qpit75,qredox,q6mode,rhsvec,selm,
     $ sigmam,sigmmo,smp100,tempc,tempk,tfx,ubacmx,ubbig,ubgamx,
     $ ubneg,ubetmx,ucospi,ugexj,ugexmo,ujflls,uphase,uspec,
     $ uzvec1,weight,wfac,xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,
     $ xhfs,xlke,xlks,xvfs,zchar,zchsq2,zchcu6,zgexj,zvclg1,zvec1)
c
c     This subroutine builds the iindx1 array. It sets up the matrix
c     structure for hybrid Newton-Raphson iteration and computes
c     initial values for the iteration variables. It optimizes these
c     in preparation for hybrid Newton-Raphson iteration. The variable
c     ker is returned as 0 if all went well, as 1 if the input
c     constraints look suspiciously poor, and as 2 if they look
c     really bad.
c
c     This subroutine is somewhat analogous to EQ6/optmzr.f, which also
c     performs optimization of iteration variables prior to hybrid
c     Newton-Raphson iteration.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       narn1  = start of the range of aqueous species; this is
c                  also the index of the solvent, water
c       narn2  = end of the range of aqueous species
c
c     Principal output:
c
c       conc   = array of species concentrations
c       iindx1 = array of indices of components appearing as
c                  matrix variables
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      include 'eqlib/eqldv.h'
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer noutpt,nttyo
c
      integer ibpxt(nxtmax),ibswx(nbtmax),iction(nbtmax),
     $ igstak(ngtmax),iindx1(kmax),ipndx1(kmax),insgf(natmax),
     $ iodb(nodbmx),iopg(nopgmx),iopt(noptmx),ipivot(kmax),
     $ istack(nstmax),ixbasp(nbtmax),jcsort(nstmax),
     $ jern1(jetmax,netmax),jern2(jetmax,netmax),jflag(nstmax),
     $ jgext(netmax),jgsort(ngtmax),jgstak(ngtmax),
     $ jjndex(nbtmax),jjsort(nstmax),jsflag(nstmax),jsitex(nstmax),
     $ jsol(nxtmax),jpflag(nptmax),jssort(nstmax),jstack(nstmax),
     $ kction(nbtmax),kkndex(nbtmax)
c
      integer narxt(ntprmx),nbasp(nbtmax),nbaspd(nbtmax),nbaspx(nbtmax),
     $ ncmpr(2,nptmax),ncosp(nbtmax),ndecsp(nbtmax),ndrs(ndrsmx),
     $ ndrsx(ndrsmx),ndrsr(2,nstmax),ndrsrd(2,nstmax),ndrsrx(2,nstmax),
     $ nfac(nbtmax),ngexsa(ietmax,jetmax,netmax),ngext(jetmax,netmax),
     $ nphasx(nstmax),nsts(nstsmx),nstsr(2,nstmax),ntfx(ntfxmx)
c
      integer nalpha(nsltmx),nmux(3,nmutmx),nmxi(2,natmax),
     $ nmxx(3,nmxmax),nslx(2,nsltmx),nsxi(2,natmax),nsxx(2,nsxmax)
c
      integer napt,nmut,nslt
c
      integer ibetmx,iebal,ielam,iern1,iern2,ifrn1,ifrn2,igas,ilrn1,
     $ ilrn2,imrn1,imrn2,ipch,ipcv,irdxc3,ixrn1,ixrn2,izmax,ka1,kat,kbt,
     $ kct,kdim,kebal,kelect,ker,ke1,ket,khydr,km1,ko2gaq,kwater,kx1,
     $ kxt,narn1,narn2,nbt,nbtd,nbti,nbw,nchlor,nct,nelect,nern1,
     $ nern2,net,ngrn1,ngrn2,ngt,nhydr,nhydx,no2gaq,npt,nredox,nst,
     $ ntfxt,ntpr
c
      integer ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,
     $ ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta
c
      logical qbassw,qchlor,qhawep,qpit75,qredox,q6mode
c
      character(len=48) ucospi(nbtmax),uspec(nstmax),uzvec1(kmax)
      character(len=48) ubacmx,ubbig,ubetmx,ubgamx,ubneg
      character(len=32) ujflls(0:njfmax)
      character(len=24) ugexmo(netmax),uphase(nptmax)
      character(len=8) ugexj(jetmax,netmax)
c
      real(8) dgpit(2,ipbtmx,napmax),dpslm(2,nsltmx),
     $ gpit(ipbtmx,napmax),palpha(ipbtmx,napmax),pmu(nmutmx),
     $ pslamn(0:ipbtmx,nsltmx),pslm(nsltmx)
c
      real(8) delam(2,nazpmx,nazpmx),dpelm(2,nazpmx,nazpmx),
     $ dselm(2,nazmmx:nazpmx),elam(nazpmx,nazpmx),pelm(nazpmx,nazpmx),
     $ selm(nazmmx:nazpmx)
c
      real(8) xhfs(nstmax),xlks(nstmax),xvfs(nstmax)
c
      real(8) dhfs(ipchmx,nstmax),dvfs(ipcvmx,nstmax)
c
      real(8) axhfs(narxmx,ntprmx,nstmax),axhfsx(narxmx,ntprmx,nstmax),
     $ axlks(narxmx,ntprmx,nstmax),axlksx(narxmx,ntprmx,nstmax),
     $ axvfs(narxmx,ntprmx,nstmax),axvfsx(narxmx,ntprmx,nstmax)
c
      real(8) adhfs(narxmx,ntprmx,ipchmx,nstmax),
     $ adhfsx(narxmx,ntprmx,ipchmx,nstmax),
     $ advfs(narxmx,ntprmx,ipcvmx,nstmax),
     $ advfsx(narxmx,ntprmx,ipcvmx,nstmax)
c
      real(8) aamatr(kmax,kmax),acflg(nstmax),acflgo(nstmax),
     $ act(nstmax),actlg(nstmax),alpha(kmax),amtb(nbtmax),
     $ azero(natmax),a3bars(natmax),beta(kmax),bfac(nbtmax),
     $ bpx(ibpxmx,nxtmax),cco2(5),cdrs(ndrsmx),cdrsx(ndrsmx),
     $ cdrtw(nstmax),cdrw(nstmax),cegexs(ietmax,jetmax,netmax),
     $ cgexj(jetmax,netmax),cjbasp(nbtmax),cnufac(nstmax),conc(nstmax),
     $ conclg(nstmax),coval(nbtmax),cpgexs(ietmax,jetmax,netmax),
     $ csts(nstsmx),delvec(kmax),dlogxw(nbtmax),efac(nbtmax),
     $ egexjc(jetmax,netmax),egexjf(jetmax,netmax),
     $ egexs(ietmax,jetmax,netmax),fsort(ngtmax),fugac(ngtmax),
     $ fugalg(ngtmax),gmmatr(kmax,kmax),loph(nptmax),losp(nstmax),
     $ lsort(nstmax)
c
      real(8) mgext(jetmax,netmax),moph(nptmax),mosp(nstmax),
     $ mrgexs(ietmax,jetmax,netmax),mtb(nbtmax),rhsvec(kmax),
     $ tfx(ntfxmx),weight(nstmax),wfac(iktmax,nxtmax),xbar(nstmax),
     $ xbarlg(nstmax),zchar(nstmax),zchsq2(nstmax),zchcu6(nstmax),
     $ zgexj(jetmax,netmax),zvclg1(kmax),zvec1(kmax)
c
      real(8) adh,adhh,adhv,aphi,bdh,bdhh,bdhv,bdot,bdoth,bdotv
c
      real(8) abar,actwlc,afcnst,al10,avcnst,a3bar,bacfmx,bbig,betamx,
     $ bfje,bfxi,bgamx,bneg,bsigmm,btmxoe,btmxoo,eh,ehfac,eps100,fje,
     $ fjeo,fo2,fo2lg,fxi,fxio,omega,omeglg,pe,press,sigmam,sigmmo,
     $ smp100,tempc,tempk,xbarw,xbarwc,xbrwlc,xbrwlg
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer iter,jfl,jlen,jlen1,jlen2,j2,kb,kc,kcol,kount,krow,nb,
     $ nb2,ncycle,ncylim,negbfc,nloop,nlopmx,npass,nplim,nr1,ns,nse,
     $ nswtch,ns1,ns2
c
      integer ilnobl
c
      logical qabsw,qawfix,qbswx,qcfxi,qcgam,qcsigm,qloop,qpracf,
     $ qtestc,qtestp,qxbarw
c
      character(len=56) uspn56,usp156,usp256
      character(len=32) ujtp
      character(len=24) ux24
c
      real(8) presg,xlke
c
      real(8) av,azdel,betfnc,bxecor,bxp1,cecorr,chfacf,chfsgm,cx,cxl,
     $ cxn,cxo,dx,fjec,fxic,lx,rlxgam,sigmmc,sigza,sigzc,sigzi,sigzm,
     $ stx,tfxc,tolbig,tolbtf,tolgpt,tolneg,tolxpt,tolzpt,xecorr,
     $ zdel,zx1,zx2
c
      real(8) coefdr,texp,tlg
c
c-----------------------------------------------------------------------
c
c     The following are iteration limits:
c
c       nlopmx = the maximum number of auto basis switching loops
c       nplim  = the maximum number of passes
c       ncylim = the maximum number of cycles
c
c     Passes refine estimates of the ionic strength, etc., the
c     activity of water, and activity coefficients of aqueous species.
c     Cycles are embedded in passes. They refine estimates of species
c     concentrations before new estimates of ionic strength, etc.,
c     are made.
c
      data nlopmx /12/,nplim  /7/,ncylim /15/
c
c     The following are tolerance parameters for the optimization:
c
      data tolbtf /0.1/
      data tolzpt /0.25/
      data tolbig /0.5/,tolneg /-0.1/
      data tolxpt /0.5/,tolgpt /0.1/
c
c     The following is needed for ncmpex.f, but is only relevant to EQ6.
c
      data qxbarw/.false./
c
c-----------------------------------------------------------------------
c
      write (noutpt,1000)
      write (nttyo,1000)
 1000 format(/' Starting Pre-Newton-Raphson Optimization.',/)
c
      qbswx = .false.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do nb = 1,nbtmax
        kkndex(nb) = 0
        kction(nb) = 0
      enddo
c
      bbig = 0.
      bneg = 0.
      ubbig = 'None'
      ubneg = 'None'
c
      bgamx = 0.
      ubgamx = 'None'
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set up the structure of the iteration matrix. Its contents may be
c     altered subsequently by automatic basis switching.
c
c     Build the iindx1 array.
c
      kebal = 0
      kwater = 0
      khydr = 0
      kelect = 0
      ko2gaq = 0
c
      ka1 = 0
      kat = 0
      ke1 = 0
      ket = 0
c
      kb = 0
      kc = 0
      do nb = 1,nbt
        ns = nbasp(nb)
        jfl = jflag(ns)
        if (jfl.ne.-1 .and. jfl.ne.30) then
          kb = kb + 1
          iindx1(kb) = nb
          uzvec1(kb) = uspec(ns)
          if (nb .eq. iebal) kebal = kb
          if (ns .eq. narn1) kwater = kb
          if (ns .eq. nhydr) khydr = kb
          if (ns .eq. nelect) kelect = kb
          if (ns .eq. no2gaq) ko2gaq = kb
        endif
      enddo
      kbt = kb
      kdim = kbt
cXXX
      if (qchlor) then
        kct = nct -1
      else
        kct = nct
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Fill the kction array. This array is used to mark the columns
c     belonging to basis species which are involved in constraints
c     placed on other basis species.
c
      do nb = 1,nbt
        ns1 = ncosp(nb)
        do kcol = 1,kbt
          nb2 = iindx1(kcol)
          ns2 = nbasp(nb2)
          if (ns2 .eq. ns1) then
            kction(nb) = kcol
            go to 100
          endif
        enddo
  100   continue
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(3) .ge. 2) then
c
c       Print the active data file basis set.
c
        write (noutpt,1010)
 1010   format(/16x,'--- Active Data File Basis Set ---',
     $  //2x,'krow   Name',30x,'Constraint',/)
c
        do krow = 1,kdim
          nb = iindx1(krow)
          ns = nbaspd(nb)
          ujtp = 'Defining equation'
          if (nb .eq. iebal) then
              ujtp = 'Electrical balance'
            elseif (ns .eq. no2gaq) then
              if ((irdxc3 .eq. -1) .or. (irdxc3 .eq. -2)) then
                ujtp = 'Eh'
              elseif (irdxc3 .eq. 1) then
                ujtp = 'Aqueous redox reaction'
              endif
            else
              jfl = jflag(ns)
              ujtp = ujflls(jfl)
          endif
          j2 = ilnobl(ujtp)
c
c         Calling sequence substitutions:
c           uspec(ns) for unam48
c
          call fmspnx(jlen,uzvec1(krow),uspn56)
          write (noutpt,1020) krow,uspn56,ujtp(1:j2)
 1020     format(1x,i4,2x,a32,2x,a)
        enddo
c
        write (noutpt,1030)
 1030   format(/1x)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The label below is a return point after an automatic basis switch.
c     Here nloop is the loop counter for auto basis switching.
c
      nloop = -1
c
      qloop = .true.
  200 nloop = nloop + 1
      if (iodb(3) .ge. 1) write (noutpt,1040) nloop
 1040 format(6x,'nloop= ',i2)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(3) .ge. 2) then
        if (qbassw) then
c
c       Print the computational basis set.
c
          write (noutpt,1050)
 1050     format(16x,'--- Basis Set Changes ---',
     $    //2x,'krow   Data File',25x,'Current',/)
c
          do krow = 1,kbt
            nb = iindx1(krow)
            ns1 = nbaspd(nb)
            ns2 = nbasp(nb)
            if (ns1 .ne. ns2) then
c
c             Calling sequence substitutions:
c               jlen1 for jlen
c               uspec(ns1) for unam48
c               usp156 for uspn56
c
              call fmspnx(jlen1,uspec(ns1),usp156)
c
c             Calling sequence substitutions:
c               jlen2 for jlen
c               uspec(ns2) for unam48
c               usp256 for uspn56
c
              call fmspnx(jlen2,uspec(ns2),usp256)
              jlen2 = min(jlen2,32)
              write (noutpt,1060) krow,usp156,usp256(1:jlen2)
 1060         format(1x,i4,2x,a32,2x,a)
            endif
          enddo
          write (noutpt,1030)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize concentrations and masses of the active basis species.
c
      do kcol = 1,kbt
        nb = iindx1(kcol)
        nse = nbasp(nb)
        jfl = jflag(nse)
        if (jfl .eq. -1) then
          conc(nse) = 0.
        elseif (jfl.ge.0 .and. jfl.le.3) then
c
c         Concentrations.
c
          ns = nbaspd(nb)
          if (nse .eq. ns) then
            conc(nse) = coval(nb)
          else
            cx = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns,nstmax)
            nr1 = ndrsr(1,ns)
            stx = -cx/cdrs(nr1)
            conc(nse) = stx*coval(nb)
          endif
        elseif (jfl.ge.7 .and. jfl.le.11) then
c
c         Alkalinity.
c
          ux24 = uspec(nse)(1:24)
          tfxc = 1.
          if (ux24(1:6) .eq. 'HCO3- ') tfxc = 1.
          if (ux24(1:6) .eq. 'CO3-- ') tfxc = 2.
          conc(nse) = coval(nb)/tfxc
        elseif (jfl .eq. 16) then
c
c         Activity.
c
          cx = coval(nb)
          conc(nse) = texp(cx)
        elseif (jfl.eq.19 .or. jfl.eq.20) then
c
c         pX and pH.
c
          cx = -coval(nb)
          conc(nse) = texp(cx)
        elseif (jfl.eq.22 .or. jfl.eq.23) then
c
c         pmX and pmH.
c
          cx = -coval(nb)
          conc(nse) = texp(cx)
        else
c
c         All other cases.
c
          conc(nse) = 1.e-7
        endif
      enddo
c
      conc(narn1) = 0.
      if (nelect .gt. 0) conc(nelect) = 0.
      if (no2gaq .gt. 0) conc(no2gaq) = 0.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the charge imbalance.
c
      zdel = 0.
      do nb = 1,nbt
        ns = nbasp(nb)
        if (ns.ge.narn1 .and. ns.le.narn2) then
          zdel = zdel + zchar(ns)*conc(ns)
        endif
      enddo
      azdel = abs(zdel)
c
c     Calculate a starting value for the SUM(i) m(i) function
c     (sigmam). Treat the calculated charge imbalance among the
c     basis species as the equivalent of a monovalent ion.
c
      call csigm(conc,jcsort,narn1,narn2,nstmax,sigmmc)
      sigmam = sigmmc + azdel
c
c     Calculate a starting value for the ionic strength (fxi).
c     Treat the calculated charge imbalance among the basis species
c     as the equivalent of a monovalent ion.
c
      call cfxi(conc,fxic,jcsort,narn1,narn2,nstmax,zchsq2)
      fxi = fxic + 0.5*azdel
c
c     Calculate a starting value for the J electrostatic moment
c     function (fje). Treat the calculated charge imbalance among
c     the basis species as the equivalent of a monovalent ion.
c
      call cfje(conc,fjec,jcsort,narn1,narn2,nstmax,zchcu6)
      fje = fjec + (-zdel/6.)
c
c     Calculate the activity coefficients of aqueous species.
c     Note that this also gets the starting value for the mole
c     fraction of water.
c
c     Calling sequence substitutions:
c       acflg for acflgc
c
      call gcoeff(abar,acflg,actwlc,adh,adhh,adhv,al10,
     $ aphi,azero,a3bar,a3bars,bdh,bdhh,bdhv,bdot,bdoth,bdotv,
     $ cco2,conc,delam,dgpit,dpelm,dpslm,dselm,elam,fje,fxi,gpit,
     $ ielam,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,
     $ ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta,insgf,
     $ iopg,ipbtmx,izmax,jcsort,nalpha,napmax,napt,narn1,narn2,
     $ natmax,nazmmx,nazpmx,nchlor,nhydr,nmut,nmutmx,nmux,nmxi,
     $ nmxmax,nmxx,nopgmx,noutpt,nslt,nsltmx,nslx,nstmax,nsxi,
     $ nsxmax,nsxx,nttyo,omega,palpha,pelm,pmu,press,pslamn,
     $ pslm,qhawep,qpit75,selm,sigmam,tempk,uspec,xbarwc,xbrwlc,
     $ zchar,zchsq2,zchcu6)
c
c     Calculate the activity coefficients of exchanger species.
c
c     Calling sequence substitutions:
c       acflg for acflgc
c
      call lamgex(acflg,cgexj,jern1,jern2,jetmax,jgext,net,
     $ netmax,nstmax,xbarlg)
c
c     Initialize the mole fraction of water.
c
      xbrwlg = xbrwlc
      xbarw = xbarwc
c
c     Copy the ionic strength, etc., and the activity coefficients.
c
      sigmmo = sigmam
      fxio = fxi
      fjeo = fje
      do ns = narn1,narn2
        acflgo(ns) = acflg(ns)
      enddo
c
c     Load the entries of the conclg array corresponding
c     to active basis species.
c
      do kcol = 1,kbt
        nb = iindx1(kcol)
        ns = nbasp(nb)
        cx = conc(ns)
        conclg(ns) = tlg(cx)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Determine whether the constraints fix the activity of water.
c     If so, the variable qawfix is set to .true.
c
      call dawfix(aamatr,cdrs,eps100,gmmatr,iindx1,iodb,
     $ irdxc3,jflag,jjndex,kbt,kkndex,kmax,narn1,nbasp,nbtmax,ncosp,
     $ ndrs,ndrsmx,ndrsr,nelect,nhydr,nodbmx,no2gaq,noutpt,nstmax,
     $ qawfix,uspec)
c
cXX   Need new coding to deal with phases assemblages that fix a(w).
cXX   Such assemblages are now merely trapped.
c
c     Coding to deal with phase assemablages that fix the activity of
c     water has not yet been implemented. Stop if this condition has
c     been detected. The needed new coding is in EQ3NR/arrset.f and
c     EQ3NR/ arrsim.f. If the following trap were not in place, the
c     matrix constructed below in this subroutine would be singular.
c     There should be no problem, however, in the coding for the hybrid
c     Newton-Raphson method.
c
      if (qawfix) then
        write (noutpt,1100)
        write (nttyo,1100)
 1100   format(/' * Error - (EQ3NR/arrset) The phase assemblage',
     $  /7x,'corresponding to the specified solubility constraints',
     $  /7x,'fixes the activity of water. This code is presently',
     $  /7x,'unable to solve problems of this type.')
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Here npass is the pass counter.
c
      npass = 0
c
c     The label below is a return point for subsequent passes. A pass
c     is an adjustment for the ionic strength, etc., the activity of
c     water, and the activity coefficients of the solute species.
c
  210 npass = npass + 1
c
c     Note:
c
c       betfnc = convergence function
c       negbfc = the number of successive iterations that the
c                convergence function betfnc has been zero or negative
c
      betfnc = 0.
      negbfc = 0
      btmxoe = 0.
      btmxoo = 0.
c
      if (iodb(3) .ge. 1) then
        write (noutpt,1110) npass
 1110   format(/11x,'npass= ',i2)
        write (noutpt,1120) sigmam,fxi,fje,xbrwlc,xbarwc
 1120   format(/13x,'sigmam= ',1pe12.5,/13x,'fxi= ',1pe12.5,
     $  /13x,'fje= ',1pe12.5,
     $  //13x,'xbrwlc= ',0pf9.5,/13x,'xbarwc= ',1pe12.5,/)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Here ncycle is the cycle counter.
c
      ncycle = 0
c
c     The label below is a return point for beginning a new cycle.
c     A cycle is an structure within a pass in which the concentrations
c     of the basis species are adjusted, while the ionic strength, etc.,
c     and the activity coefficients are held constant.
c
  220 ncycle = ncycle + 1
      if (iodb(3) .ge. 1) write (noutpt,1200) ncycle
 1200 format(16x,'ncycle= ',i2)
      ker = 0
      qabsw = iopt(11).ge.1 .and. npass.eq.1 .and. ncycle.eq.1
c
c     Set up the zvclg1 and actlg array entries for the basis species
c     whose concentrations do not have to be estimated simultaneously.
c
      do kcol = 1,kbt
        nb = iindx1(kcol)
        ns = nbasp(nb)
        if (ns.ge.narn1 .and. ns.le.narn2) then
          jfl = jflag(ns)
          if (kcol .eq. kwater) then
            zvclg1(kcol) = xbrwlg
            actlg(narn1) = xbrwlg + acflg(narn1)
          elseif (jfl .le. 15) then
            zvclg1(kcol) = conclg(ns)
            actlg(ns) = conclg(ns) + acflg(ns)
          elseif (jfl .eq. 16) then
            zvclg1(kcol) = coval(ns) - acflg(ns)
            conclg(ns) = zvclg1(kcol)
            actlg(ns) = conclg(ns) + acflg(ns)
          elseif (jfl.eq.19 .or. jfl.eq.20) then
            zvclg1(kcol) = -coval(ns) - acflg(ns)
            conclg(ns) = zvclg1(kcol)
            actlg(ns) = conclg(ns) + acflg(ns)
          elseif (jfl.eq.22 .or. jfl.eq.23) then
            zvclg1(kcol) = -coval(ns)
            conclg(ns) = zvclg1(kcol)
            actlg(ns) = conclg(ns) + acflg(ns)
          endif
        elseif (ns.ge.nern1 .and. ns.le.nern2) then
          zvclg1(kcol) = conclg(ns)
          actlg(ns) = cjbasp(nb)*(xbarlg(ns) + acflg(ns))
        else
          zvclg1(kcol) = tlg(coval(nb)) + xbarlg(ns)
          actlg(ns) = xbarlg(ns) + acflg(ns)
        endif
      enddo
c
      if (irdxc3 .eq. 0) then
        if (ko2gaq .gt. 0) then
          zvclg1(ko2gaq) = fo2lg
          actlg(no2gaq) = fo2lg
        endif
        if (kelect .gt. 0) then
          zvclg1(kelect) = -pe
          actlg(nelect) = -pe
        endif
      endif
c
c     Make starting estimates that must be evaluated simultaneously.
c     These include all cases of equilibrium constraints and compount
c     activity constraints (e.g. pHCl) and any case in which log fO2
c     log fO2 is constrained by Eh, pe-, or a redox couple.
c
      call arrsim(aamatr,acflg,actlg,bbig,cdrs,cjbasp,cnufac,
     $ conc,conclg,coval,delvec,dlogxw,eh,ehfac,eps100,gmmatr,iction,
     $ iindx1,iodb,ipivot,irdxc3,ixbasp,jcsort,jflag,jjndex,kbt,ker,
     $ khydr,kkndex,kmax,kwater,narn1,narn2,nbasp,nbt,nbti,nbtmax,
     $ nbw,ncosp,ndecsp,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,nhydr,
     $ nodbmx,no2gaq,noutpt,npass,nredox,nstmax,nttyo,omega,qawfix,
     $ rhsvec,ucospi,uspec,xbar,xbarlg,xbarw,xbrwlg,xlke,xlks,
     $ zchar,zvclg1)
c
c     Recalculate the concentrations, etc., of dependent species.
c
      call ncmpex(acflg,act,actlg,cdrs,cegexs,cgexj,conc,
     $ conclg,cpgexs,egexjc,egexjf,egexs,eps100,fo2,fo2lg,fsort,
     $ fugac,fugalg,iern1,iern2,ietmax,ifrn1,ifrn2,igas,igstak,
     $ iindx1,ilrn1,ilrn2,imrn1,imrn2,istack,ixrn1,ixrn2,jcsort,
     $ jern1,jern2,jetmax,jflag,jgext,jgsort,jgstak,jjsort,jpflag,
     $ jsflag,jsitex,jssort,jstack,kbt,kdim,kelect,kmax,km1,ko2gaq,
     $ kwater,kxt,loph,losp,lsort,mgext,mrgexs,mtb,moph,mosp,narn1,
     $ narn2,nbasp,nbt,nbtmax,ncmpr,ndrs,ndrsmx,ndrsr,nelect,nern1,
     $ nern2,netmax,ngexsa,ngext,ngrn1,ngrn2,ngt,ngtmax,noutpt,
     $ no2gaq,nphasx,npt,nptmax,nst,nstmax,nttyo,omega,omeglg,
     $ press,qxbarw,q6mode,ugexj,ugexmo,uphase,uspec,xbar,xbarlg,
     $ xbarw,xbarwc,xbrwlc,xbrwlg,xlks,zchar,zgexj,zvclg1,zvec1)
c
      xbarw = xbar(narn1)
      xbrwlg = xbarlg(narn1)
c
c     Compute the residuals.
c
      call betas(acflg,actlg,afcnst,alpha,amtb,bbig,beta,
     $ betamx,bneg,cdrs,conc,conclg,coval,csts,eh,ehfac,fo2lg,
     $ ibetmx,iebal,iindx1,irdxc3,jcsort,jflag,jsflag,jssort,kbt,
     $ kdim,kelect,khydr,kmax,km1,ko2gaq,kwater,kxt,mtb,mosp,
     $ narn1,narn2,nbasp,nbtmax,ncosp,ndrs,ndrsmx,ndrsr,nelect,
     $ nern1,nern2,nhydr,noutpt,no2gaq,nredox,nst,nstmax,nsts,
     $ nstsmx,nstsr,ntfx,ntfxmx,ntfxt,nttyo,omega,qredox,q6mode,
     $ tfx,ubbig,ubneg,ubetmx,uspec,uzvec1,weight,xbrwlg,xlke,
     $ xlks,zchar)
c
c     Calculate the beta convergence function.
c
      betfnc = 0.
      if (mod(ncycle,2) .eq. 0) then
        if (btmxoe .ge. smp100) betfnc = (btmxoe - betamx)/btmxoe
        btmxoe = betamx
      else
        if (btmxoo .ge. smp100) betfnc = (btmxoo - betamx)/btmxoo
        btmxoo = betamx
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print values of master iteration variables.
c
      if (iodb(3) .ge. 2) then
        write (noutpt,1230)
 1230   format(//10x,'--- Pre-Newton-Raphson Optimization Summary ---',
     $  //2x,'kcol   Name',32x,'zvclg1      zvec1',/)
        do kcol = 1,kdim
          nb = iindx1(kcol)
          ns = nbasp(nb)
          zx1 = zvclg1(kcol)
          zx2 = texp(zx1)
c
c         Calling sequence substitutions:
c           uspec(ns) for unam48
c
          call fmspnx(jlen,uspec(ns),uspn56)
          write (noutpt,1240) kcol,uspn56,zx1,zx2
 1240     format(1x,i4,2x,a32,2x,f10.4,2x,1pe12.5)
        enddo
c
        write (noutpt,1250)
 1250   format(/2x,'krow   Name',32x,'Beta',/)
        do krow = 1,kdim
          nb = iindx1(krow)
          ns = nbaspd(nb)
c
c         Calling sequence substitutions:
c           uspec(ns) for unam48
c
          call fmspnx(jlen,uspec(ns),uspn56)
          write (noutpt,1260) krow,uspn56,beta(krow)
 1260     format(1x,i4,2x,a32,2x,1pe12.5)
        enddo
        write (noutpt,1030)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Identify the dominant species in each mass balance and
c     compute the corresponding exponent for a continued
c     fraction correction.
c
      call cfracf(cdrs,csts,efac,jcsort,jflag,jssort,kmax,mosp,
     $ narn1,narn2,nbasp,nbaspd,nbt,nbtmax,ndrs,ndrsmx,ndrsr,nern1,
     $ nern2,nfac,nst,nstmax,nsts,nstsmx,nstsr,q6mode,weight)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Filter the data obtained from EQLIB/cfracf. The following
c     filters are specific to EQ3NR.
c
      if (iebal .gt. 0) then
c
c       For this species, a charge balance constraint is used instead
c       of a mass balance constraint.
c
        nfac(iebal) = 0
        efac(iebal) = 1.0
      endif
c
      if (nbw .gt. 0) then
c
c       For water, a mole fraction equation is used instead of a mass
c       balance constraint.
c
        nfac(nbw) = 0
        efac(nbw) = 1.0
      endif
c
      do nb = 1,nbt
        ns = nbasp(nb)
        if (ns.eq.nhydr .or. ns.eq.nhydx) then
c
c         For H+ or OH-, mass balance constraints are rarely used.
c
          if (jflag(ns) .gt. 15) then
c
c           In the present case, a mass balance constraint is not
c           being used.
c
            nfac(nb) = 0
            efac(nb) = 1.0
          endif
        elseif (ns.eq.no2gaq .or. ns.eq.nelect) then
c
c         For O2(g,aq) or e-, mass balance contraints are not used.
c
          nfac(nb) = 0
          efac(nb) = 1.0
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(3) .ge. 3) then
c
c       Write a table containing the preliminary results.
c
        kount = 0
        do kcol = 1,kbt
          nb = iindx1(kcol)
          ns = nbasp(nb)
          ns2 = nfac(nb)
          if (ns2.ne.0 .and. ns2.ne.ns) kount = kount + 1
        enddo
c
        if (kount .gt. 0) then
          write (noutpt,1300)
 1300     format(16x,'--- Mass Balance Dominants ---',/)
          write (noutpt,1310)
 1310     format(4x,'Master Species',21x,'Dominant Species',/)
          do kcol = 1,kbt
            nb = iindx1(kcol)
            ns = nbasp(nb)
            ns2 = nfac(nb)
            if (ns2.ne.0 .and. ns2.ne.ns) then
c
c             Calling sequence substitutions:
c               jlen1 for jlen
c               uspec(ns) for unam48
c               usp156 for uspn56
c
              call fmspnx(jlen1,uspec(ns),usp156)
c
c             Calling sequence substitutions:
c               jlen2 for jlen
c               uspec(ns2) for unam48
c               usp256 for uspn56
c
              call fmspnx(jlen2,uspec(ns2),usp256)
              jlen2 = min(jlen2,32)
              write (noutpt,1320) usp156,usp256(1:jlen2)
 1320         format(2x,a32,3x,a)
            endif
          enddo
          write (noutpt,1030)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set up the bfac correction factor array. In the continued fraction
c     method, m(new) = m(old)/bfac, where bfac = (beta + 1)**efac.
c     If the same species dominates more than one mass balance,
c     then this algorithm can be applied to only one of the associated
c     basis species. Otherwise, oscillatory behavior will occur. In each
c     set of mass balances with a common dominating species, find the
c     mass balance with the greatest bfac factor. This is usually nearly
c     equivalent to finding the mass balance with the greater beta
c     residual, as efac often has a value of unity.
c
      call gbfac(beta,bfac,efac,iindx1,kbt,kmax,nbt,nbtmax,nfac)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(3).ge.3 .and. .not.qbswx) then
c
c       Write a table containing the modified results.
c
        write (noutpt,1330)
 1330   format(/16x,'--- Factors for Continued Fraction',
     $  ' Corrections ---',/)
        write (noutpt,1340)
 1340   format(4x,'Master Species',22x,'bfac',10x,'efac',/)
        do nb = 1,nbt
          if (bfac(nb) .gt. 0.) then
            ns = nbasp(nb)
c
c           Calling sequence substitutions:
c             uspec(ns) for unam48
c
            call fmspnx(jlen,uspec(ns),uspn56)
            write (noutpt,1350) uspn56,bfac(nb),efac(nb)
 1350       format(2x,a32,3x,1pe12.5,3x,1pe12.5)
          endif
        enddo
        write (noutpt,1030)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qabsw .and. qloop .and. nloop.lt.nlopmx) then
c
c       In automatic basis switching mode (iopt(11) .ge. 1), try to
c       first reduce the magntiude of large positive mass balance
c       residuals by making one or more basis switches.
c
        call absswa(adhfs,adhfsx,advfs,advfsx,avcnst,axhfs,
     $  axhfsx,axlks,axlksx,axvfs,axvfsx,beta,cdrs,cdrsx,cdrtw,cdrw,
     $  csts,dhfs,dvfs,efac,eps100,ibswx,iebal,iindx1,iodb,ipch,ipchmx,
     $  ipcv,ipcvmx,jcsort,jflag,jsflag,jssort,kbt,kmax,mosp,narn1,
     $  narn2,narxmx,narxt,nbasp,nbaspd,nbaspx,nbt,nbtmax,nbw,ncosp,
     $  ndrs,ndrsmx,ndrsr,ndrsrx,ndrsx,nelect,nhydr,nodbmx,no2gaq,
     $  noutpt,nst,nstmax,nsts,nstsmx,nstsr,nswtch,ntpr,ntprmx,nttyo,
     $  presg,press,qbassw,qbswx,q6mode,tempc,uspec,uzvec1,weight,
     $  xvfs,xlks,xhfs)
c
        if (nswtch .le. 0) then
c
c         No switches were made.
c
          qloop = .false.
          go to 250
        endif
c
c       Reset the ixbasp and cjbasp arrays. The former is a flag
c       array, each member of which denotes whether the
c       thermodynamic activity of the corresponding basis species
c       is defined in terms of molality (= 0) or mole fraction (= 1).
c       The cjbasp array contains any site stoichiometric factors
c       associated with the operational basis species.
c
        call gibasp(cgexj,cjbasp,iern1,ixbasp,jern1,jern2,
     $  jetmax,jgext,narn1,narn2,nbasp,nbt,nbtmax,nern1,nern2,
     $  netmax,nphasx,nstmax)
c
c       Null some arrays.
c
        do ns = 1,nstmax
          conc(ns) = 0.
        enddo
c
        av = -99999.
        call initav(conclg,nstmax,av)
c
        write (noutpt,1370) nloop,nswtch
        write (nttyo,1370) nloop,nswtch
 1370   format(8x,'Completed loop ',i3,' after making ',i3,
     $  ' basis switches.')
c
c       Go back for another loop.
c
        go to 200
      endif
  250 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the electrical balance residual.
c
      bxecor = 0.
      xecorr = 0.
      if (kebal .gt. 0) then
        ns = nbaspd(iebal)
        cecorr = -alpha(kebal)/zchar(ns)
        xecorr = zchsq2(ns)*cecorr
        call gszm(conc,jcsort,narn1,narn2,nstmax,sigza,sigzc,
     $  sigzi,sigzm,zchar)
        bxecor = abs(xecorr)/sigzm
      endif
c
      if (iodb(3) .ge. 1) then
c
c       Calling sequence substitutions:
c         jlen1 for jlen
c         ubbig for unam48
c         usp156 for uspn56
c
        call fmspnx(jlen1,ubbig,usp156)
c
c       Calling sequence substitutions:
c         jlen2 for jlen
c         ubneg for unam48
c         usp256 for uspn56
c
        call fmspnx(jlen2,ubneg,usp256)
c
        write (noutpt,1500) betamx,betfnc,bbig,usp156(1:jlen1),
     $  bneg,usp256(1:jlen2)
 1500   format(18x,'betamx= ',1pe12.5,', betfnc= ',1pe12.5,
     $  /18x,'  bbig= ',1pe12.5,', ubbig= ',a,
     $  /18x,'  bneg= ',1pe12.5,', ubneg= ',a,/)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Test the balance residuals for mass and alkalinity to see if
c     another cycle should be made before attempting to make an improved
c     estimate of the ionic strength.
c
      qtestc = bbig.le.tolbig .and. bneg.ge.tolneg .and.
     $  betamx.le.tolxpt .and. bxecor.le.tolzpt
c
c     Quit doing cycles if:
c
c       1. The cycle convergence criteria are met.
c       2. The maximum number of cycles have been done.
c       3. The convergence function betfnc indicates that
c          the cycles are not converging.
c
      if (qtestc) go to 300
      if (ncycle .ge. ncylim) go to 300
      if (ncycle.gt.2 .and. betfnc.le.tolbtf) then
        negbfc = negbfc + 1
        if (negbfc .ge. 3) go to 300
      else
        negbfc = 0
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make improvements in concentration estimates and go back for
c     another cycle. The algorithm employed here is the modified
c     continued fraction method. Note that some constraints are placed
c     on the use of the bfac correction factors. One is to avoid a
c     blow-out in taking the base ten logarithm of such a factor, which
c     would occur if the factor had a zero or negative value. Another
c     limits the magnitude of the change in one step.
c
      do kcol = 1,kbt
        nb = iindx1(kcol)
        ns = nbasp(nb)
        if (nb.ne.iebal .and. ns.ne.narn1) then
          jfl = jflag(ns)
          bxp1 = beta(kcol) + 1.
          if (jfl.ge.0 .and. jfl.le.3) then
            dx = bfac(nb)
            lx = tlg(dx)
            if (lx .gt. 20.) lx = 20.
            conclg(ns) = conclg(ns) - lx
          elseif (jfl.ge.7 .and. jfl.le.11) then
            if (bxp1 .le. 0.) bxp1 = 1.e-20
            lx = tlg(bxp1)
            if (lx .gt. 20.) lx = 20.
            conclg(ns) = conclg(ns) - lx
          endif
        endif
      enddo
c
      if (iebal .gt. 0) then
c
c       Electrical balance correction. Skip if mass balance
c       residuals are way off.
c
        if (bbig.le.0.25 .and. bneg.ge.-0.25) then
          ns = nbasp(iebal)
          cxo = conclg(ns)
          cxn = conc(ns) + cecorr
          if (cecorr .ge. 0.) then
            cx = tlg(cxn)
            cxl = cxo + 2.0
            conclg(ns) = min(cx,cxl)
          else
            if (cxn .lt. 0.) cxn = 0.
            cx = tlg(cxn)
            cxl = cxo - 2.0
            conclg(ns) = max(cx,cxl)
          endif
        endif
      endif
c
c     Go back for another cycle.
c
      go to 220
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The cycles for the current pass have been completed. Test to
c     see if another pass should be made.
c
  300 write (noutpt,1620) npass,ncycle
      write (nttyo,1620) npass,ncycle
 1620 format(13x,'Completed pass ',i3,' in ',i3,' cycles.')
c
c     Save the current values of the ionic strength, etc., and of the
c     the ativity coefficients.
c
      sigmmo = sigmam
      fxio = fxi
      fjeo = fje
      do ns = narn1,narn2
        acflgo(ns)=acflg(ns)
      enddo
c
c     Determine the maximum allowed change factors for "Sigma m",
c     etc. (chfsgm) and the log activity coefficients for aqueous
c     species (chfacf).
c
      chfsgm = 1.3
      if (sigmam .le. 2.e-1) chfsgm = 5.
      if (sigmam .le. 1.e-2) chfsgm = 10.
      if (sigmam .le. 1.e-3) chfsgm = 100.
      if (qtestc .and. chfsgm.lt.100.) chfsgm = 100.
      chfacf = 0.05
c
c     Note: setting iter = 0 would fix the activity coefficients
c     in concentrated solutions. Setting it to a high value (here
c     99999) insures that the activity coefficients are recalculated
c     as desired.
c
      iter = 99999
      rlxgam = 1.0
      qpracf = iodb(3) .ge. 4
c
      call ngcadv(abar,acflg,acflgo,actwlc,adh,adhh,adhv,
     $ afcnst,al10,aphi,azero,a3bar,a3bars,bacfmx,bdh,bdhh,bdhv,
     $ bdot,bdoth,bdotv,bgamx,bpx,bsigmm,bfje,bfxi,cco2,cgexj,
     $ chfacf,chfsgm,conc,delam,dgpit,dpelm,dpslm,dselm,elam,
     $ eps100,fje,fjeo,fxi,fxio,gpit,ibpxt,ielam,ifcphi1,ifcphi2,
     $ ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,ilcphi1,ilcphi2,ilnnn,
     $ iln2n,ilpsi1,ilpsi2,ilzeta,insgf,iopg,iter,ipndx1,ixrn1,
     $ ixrn2,izmax,jcsort,jern1,jern2,jgext,jsol,kx1,kxt,nalpha,
     $ napt,narn1,narn2,nchlor,ncmpr,net,nhydr,nmut,nmux,nmxi,
     $ nmxx,noutpt,nslt,nslx,nst,nsxi,nsxx,nttyo,omega,palpha,
     $ pelm,pmu,press,pslamn,pslm,qhawep,qpit75,qpracf,q6mode,
     $ rlxgam,selm,sigmam,sigmmo,tempk,ubacmx,ubgamx,uphase,
     $ uspec,wfac,xbar,xbarlg,xbarwc,xbrwlc,zchar,zchcu6,zchsq2)
c
c     Recalculate the concentrations, etc., of dependent species.
c
      call ncmpex(acflg,act,actlg,cdrs,cegexs,cgexj,conc,
     $ conclg,cpgexs,egexjc,egexjf,egexs,eps100,fo2,fo2lg,fsort,
     $ fugac,fugalg,iern1,iern2,ietmax,ifrn1,ifrn2,igas,igstak,
     $ iindx1,ilrn1,ilrn2,imrn1,imrn2,istack,ixrn1,ixrn2,jcsort,
     $ jern1,jern2,jetmax,jflag,jgext,jgsort,jgstak,jjsort,jpflag,
     $ jsflag,jsitex,jssort,jstack,kbt,kdim,kelect,kmax,km1,ko2gaq,
     $ kwater,kxt,loph,losp,lsort,mgext,mrgexs,mtb,moph,mosp,narn1,
     $ narn2,nbasp,nbt,nbtmax,ncmpr,ndrs,ndrsmx,ndrsr,nelect,nern1,
     $ nern2,netmax,ngexsa,ngext,ngrn1,ngrn2,ngt,ngtmax,noutpt,
     $ no2gaq,nphasx,npt,nptmax,nst,nstmax,nttyo,omega,omeglg,
     $ press,qxbarw,q6mode,ugexj,ugexmo,uphase,uspec,xbar,xbarlg,
     $ xbarw,xbarwc,xbrwlc,xbrwlg,xlks,zchar,zgexj,zvclg1,zvec1)
c
      xbarw = xbar(narn1)
      xbrwlg = xbarlg(narn1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(3) .ge. 1) then
        write (noutpt,1510) bsigmm
 1510   format(/13x,'bsigmm= ',1pe12.5)
        write (noutpt,1520) bfxi
 1520   format(13x,'bfxi= ',1pe12.5)
        write (noutpt,1522) bfje
 1522   format(13x,'bfje= ',1pe12.5)
        j2 = ilnobl(ubgamx(1:24))
        write (noutpt,1530) bgamx,ubgamx(1:j2)
 1530   format(13x,'bgamx= ',1pe12.5,', ubgamx= ',a)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      qcsigm =abs(bsigmm) .le. tolgpt
      qcfxi = abs(bfxi) .le. tolgpt
      qcgam = bgamx .le. tolgpt
c
      qtestp = qcsigm .and. qcfxi .and. qcgam
c
c     Are pass criteria satisfied?
c
      if (.not.qtestp) then
c
c       Pass criteria are not satisfied. Test for maximum number
c       of passes.
c
        if (npass .ge. nplim) then
c
c         Quit. Optimization ended outside requested limits
c         because the pass requirements were not satisfied.
c
          if (ker .lt. 2) then
            write (noutpt,1600)
            write (nttyo,1600)
 1600       format(/'   Done. Optimization ended outside requested',
     $      ' limits.',/)
          else
            write (noutpt,1610)
            write (nttyo,1610)
 1610       format(/'   Done. Optimization ended outside allowable',
     $      ' limits.',/)
          endif
          go to 999
        endif
c
c       Do another pass.
c
        go to 210
      endif
c
c     Are cycle criteria satisfied?
c
      if (qtestc) then
c
c       Yes, optimization succeeded.
c
        write (noutpt,1630)
        write (nttyo,1630)
 1630   format(/'   Done. Optimization ended within requested',
     $  ' limits.',/)
        go to 999
c
      elseif (npass .le. 2) then
c
c       The pass convergence criteria are satisfied, but the cycle
c       convergence criteria are not. Try another pass.
c
        go to 210
      else
c
c       Quit. Optimization ended outside requested limits
c       because cycle requirements were not met.
c
        if (ker .lt. 2) then
          write (noutpt,1600)
          write (nttyo,1600)
        else
          write (noutpt,1610)
          write (nttyo,1610)
        endif
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
