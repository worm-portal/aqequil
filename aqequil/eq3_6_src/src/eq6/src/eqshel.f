      subroutine eqshel(aadh,aadhh,aadhv,aamatr,aaphi,abar,abdh,
     $ abdhh,abdhv,abdot,abdoth,abdotv,acflg,acflgo,act,actlg,
     $ adadhh,adadhv,adbdhh,adbdhv,adbdth,adbdtv,adh,adhh,adhv,
     $ afcnst,affp,affs,alpha,al10,amtb,aphi,apx,avcnst,azero,
     $ a3bar,a3bars,bacfmx,bbig,bdh,bdhh,bdhv,bdot,bdoth,bdotv,
     $ beta,betamx,betao,bgamx,bneg,bpx,cbsr,cco2,cdac,cegexs,
     $ cesr,cess,cdrs,cdrsd,cdrsx,cdrtw,cdrw,cjbasp,cnufac,conc,
     $ conclg,cpgexs,cscale,csigma,csts,dadhh,dadhv,dbdhh,dbdhv,
     $ dbdth,dbdtv,deltim,delvco,delvec,delxi,dlogxw,dlxmin,
     $ drer0,drir0,dzvc0,d1zvc1,eact,egers,egexjc,egexjf,egexs,
     $ eh,ehfac,elecsr,eps100,farad,fje,fjeo,fkrc,fo2,fo2lg,
     $ fsort,fugac,fugalg,fxi,fxio,gmmatr,hact,iact,iapxt,ibpxt,
     $ ibswx,ielam,ier,iern1,iern2,ifcphi1,ifcphi2,ifnnn,ifn2n,
     $ ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igstak,iindx0,iindx1,
     $ ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,
     $ ilzeta,imech,imrn1,imrn2,insgf,iodb,iopg,iopt,ipch,
     $ ipivot,ipndx1,ipcv,istack,iter,itermx,ixbasp,ixrn1,ixrn2,
     $ izmax,jcode,jcsort,jflag,jgsort,jgstak,jjsort,jpflag,
     $ jpress,jptffl,jreac,jsflag,jsitex,jsol,jssort,jstack,
     $ jtemp,kbt,kction,kdim,kelect,khydr,khydx,km10,km1,kmt,
     $ kmt0,ko2gaq,kpsat,kpsst,krdxsp,kwater,kx1,kx10,kxt,kxt0,
     $ loph,losp,lsort,modr,modr0,moph,morr,morr0,mosp,mrgers,
     $ mrgexs,mtb,mtbaq,mtb0,mte,mteaq,mwtrc,narn1,narn2,narxt,
     $ nat,nbasp,nbaspd,nbaspx,nbt,nbtd,nbw,nchlor,ncmpr,ncorr,
     $ nct,ndac,ndact,ndrs,ndrsd,ndrsx,ndrsr,ndrsrd,ndrsrx,nelect,
     $ nern1,nern2,ness,nessr,net,nfrn1,nfrn2,ngrn1,ngrn2,ngt,
     $ nhydr,nhydx,nlrn1,nlrn2,nlt,nmrn1,nmrn2,nmt,nord,noutpt,
     $ no2gaq,npchk,nphasx,npslmx,npt,nrct,nrdxsp,nrk,nrndex,nsk,
     $ nsslmx,nst,nsts,nstsr,ntpr,ntrymx,ntprt,nttyo,nxridx,nxrn1,
     $ nxrn2,nxt,nweope,nwndpc,omega,omeglg,prcinf,press,pressb,
     $ pressd,ptk,qbassw,qbseqc,qbye,qcnpre,qcntmp,qhawep,qmod,
     $ qoptmz,qpit75,qredox,qriinf,qscon,qshoot,qstart,qtrch,qtvchk,
     $ qxknph,q6mode,rcnstv,rconst,rkb,rhsvec,rirec0,rk,rreacn,
     $ rreac1,rrelr0,rrelr1,rtcnst,rxbar,screwd,sidrph,sidrsp,
     $ sigmam,sigmmo,smp100,tdays,tempc,tempcb,tempcd,tempcu,
     $ tempc0,tempk,timemx,time0,time1,tiplol,tiplot,tiprnl,tiprnt,
     $ tistsv,tolbt,toldl,tolsat,tolsst,tolxst,trkb,ttk,ubacmx,
     $ ubgamx,udac,ufixf,ugermo,ulbeta,uldel,uphase,ureac,uspec,
     $ uzvec0,uzvec1,vreac,weight,wfac,wodr,worr,xbar,xbarlg,
     $ xbarw,xbarwc,xbrwlc,xbrwlg,xgers,xirct,xirct0,xistsv,xi0,
     $ xi1,zchar,zchcu6,zchsq2,zklogu,zvclg0,zvclg1,zvec0,zvec1)
c
c     This subroutine provides a shell around EQ6/eqcalc.f that allows
c     fine adjustments to the value of reaction progress in order to
c     step over fine singularities, such as may exist in the middle of
c     a major jump in the oxygen fugacity.
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
c       qbye comes in .true. if EQ6/path.f just changed the phase
c         assemblage by deleting one or more minerals. This
c         instructs this subroutine to print the index (iindx1)
c         structure of the system of equations as it does on the
c         starting call and whenever this subroutine itself adds or
c         deletes a phase.
c
c       iter is the number of Newton-Raphson iterations that
c         EQLIB/newton.f performed.
c
c       ier is an error parameter which is returned as 0 if the
c         calculations in this subroutine converged successfully and
c         there were no violations of possible open system constraints.
c
c       ier    = error flag:
c
c                  Values returned from EQ6/eqphas.f:
c                    =    0  Okay
c                    =   10  Go back and take a smaller step size to
c                              avoid exceeding the supersaturation
c                              tolerance (tolsst). This is done only
c                              when an appropriate set of conditions
c                              is satisified. It isn't necessary for
c                              EQ6/path.f to analyze the situation
c                              when this value is returned to it by
c                              EQ6/eqshel.f.
c                    =   20  Hit the maximum number of tries to find
c                              the correct phase assemblage
c                    =   30  Caught in a region of computational
c                              instability about the phase boundary for
c                              an appearing phase. Iteration fails when
c                              the phase is added to the phase
c                              assemblage.
c                    =   40  Caught in a region of computational
c                              instability about the phase boundary for
c                              a disappearing phase. The system is
c                              too supersaturated with respect to a
c                              phase if that phase is deleted from the
c                              phase assemblage.
c                    =   50  Caught in a region of computational
c                              instability associated with solvent
c                              water. The amount of water in the
c                              system is probably very low.
c                    =   60  Caught in a region of computational
c                              instability associated with the
c                              master redox variable. This is probably
c                              associated with a redox jump.
c                    =   90  Caught in a region of computational
c                              instability associated with the
c                              aqueous activity coefficient model.
c                    =  100  Detected out-of-range values for variables
c                              associated with the basis species before
c                              starting iteration.
c                    =  110  Detected out-of-range values for variables
c                              associated with the basis species after
c                              an iteration crash.
c                    =  150  Calculation failed, no diagnostics were
c                              generated.
c
c                  Values returned by the present subroutine:
c                    =    0  Okay
c                    =   10  Go back and take a smaller step size to
c                              avoid exceeding the supersaturation
c                              tolerance (tolsst)
c                    =  170  Too much of a phase was destroyed under
c                              the flow-through open system model;
c                              go back and first move part of the
c                              mass of protected phases in the ES to
c                              the PRS
c                    =  180  One of a number of problems occurred which
c                              may be resolvable, at least partially,
c                              by going back and cutting the step size
c                    =  190  Need to slide over a region of
c                              computational instability, but sliding
c                              is inhibited; go back, but terminate
c                              work on the current problem
c
c-----------------------------------------------------------------------
c
c     Modules.
c
c     The module mod6pt contains data required to evaluate Pitzer's
c     equations.
c
      use mod6pt
c
c     The module mod6xf contains most of the standard-state
c     thermodynamic data.
c
      use mod6xf
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      include 'eqlib/eqlpar.h'
      include 'eqlib/eqldv.h'
      include 'eqlib/eqlge.h'
      include 'eqlib/eql1s.h'
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer noutpt,nttyo
c
      integer iact(imchmx,2,nrctmx),imech(2,nrctmx),ixbasp(nbtmax),
     $ ndac(ndctmx,imchmx,2,nrctmx),ndact(imchmx,2,nrctmx),
     $ nrk(2,nrctmx),nsk(nrctmx)
c
      integer jcode(nrctmx),jreac(nrctmx),nrndex(nrctmx),nxridx(nrctmx)
c
      integer iapxt(nxtmax),ibpxt(nxtmax),ibswx(nbtmax),igstak(ngtmax),
     $ iindx0(kmax),iindx1(kmax),ipivot(kmax),ipndx1(kmax),
     $ insgf(natmax),iodb(nodbmx),iopg(nopgmx),iopt(noptmx),
     $ istack(nstmax),jcsort(nstmax),jflag(nstmax),jgsort(ngtmax),
     $ jgstak(ngtmax),jjsort(nstmax),jpflag(nptmax),jsflag(nstmax),
     $ jsitex(nstmax),jsol(nxtmax),jssort(nstmax),jstack(nstmax),
     $ kction(nbtmax)
c
      integer narxt(ntprmx),nbasp(nbtmax),nbaspd(nbtmax),nbaspx(nbtmax),
     $ ncmpr(2,nptmax),ness(nessmx),nessr(2,nstmax),ndrs(ndrsmx),
     $ ndrsd(ndrsmx),ndrsx(ndrsmx),ndrsr(2,nstmax),ndrsrd(2,nstmax),
     $ ndrsrx(2,nstmax),npchk(nptmax),nphasx(nstmax),nsts(nstsmx),
     $ nstsr(2,nstmax)
c
      integer iern1,iern2,ifrn1,ifrn2,ilrn1,ilrn2,imrn1,imrn2,
     $ ixrn1,ixrn2
c
      integer ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,
     $ ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta
c
      integer nat,nbt,nct,net,ngt,nlt,nmt,npt,nst,nxt
c
      integer narn1,narn2,nern1,nern2,nfrn1,nfrn2,ngrn1,ngrn2,
     $ nlrn1,nlrn2,nmrn1,nmrn2,nxrn1,nxrn2
c
      integer ielam,ier,igas,ipch,ipcv,iter,itermx,izmax,jpress,jptffl,
     $ jtemp,kbt,kdim,kelect,khydr,khydx,km1,km10,kmt,kmt0,ko2gaq,
     $ kpsat,kpsst,krdxsp,kwater,kx1,kx10,kxt,kxt0,nbtd,nbw,nchlor,
     $ ncorr,nelect,nhydr,nhydx,nord,no2gaq,npslmx,nrct,nrdxsp,nsslmx,
     $ ntpr,ntprt,ntrymx,nweope,nwndpc
c
      logical qxknph(nptmax)
c
      logical qbassw,qbseqc,qbye,qcnpre,qcntmp,qhawep,qmod,qoptmz,
     $ qpit75,qredox,qriinf,qshoot,qscon,qstart,qtrch,qtvchk,q6mode
c
      character*48 uspec(nstmax),uzvec0(kmax),uzvec1(kmax)
      character*48 ubacmx,ubgamx
      character*24 ugermo(nertmx),ureac(nrctmx)
      character*24 udac(ndctmx,imchmx,2,nrctmx)
      character*24 uphase(nptmax)
      character*8 ulbeta(kmax),uldel(kmax)
      character*8 ufixf
c
      real*8 aadh(narxmx,ntprmx),aadhh(narxmx,ntprmx),
     $ aadhv(narxmx,ntprmx),aaphi(narxmx,ntprmx),
     $ abdh(narxmx,ntprmx),abdhh(narxmx,ntprmx),
     $ abdhv(narxmx,ntprmx),abdot(narxmx,ntprmx),
     $ abdoth(narxmx,ntprmx),abdotv(narxmx,ntprmx)
c
      real*8 dadhh(ipchmx),dadhv(ipcvmx),dbdhh(ipchmx),dbdhv(ipcvmx),
     $ dbdth(ipchmx),dbdtv(ipcvmx)
c
      real*8 adadhh(narxmx,ntprmx,ipchmx),
     $ adadhv(narxmx,ntprmx,ipcvmx),adbdhh(narxmx,ntprmx,ipchmx),
     $ adbdhv(narxmx,ntprmx,ipcvmx),adbdth(narxmx,ntprmx,ipchmx),
     $ adbdtv(narxmx,ntprmx,ipcvmx)
c
       real*8 cdac(ndctmx,imchmx,2,nrctmx),csigma(imchmx,2,nrctmx),
     $ eact(imchmx,2,nrctmx),fkrc(nrctmx),hact(imchmx,2,nrctmx),
     $ rkb(imchmx,2,nrctmx),trkb(imchmx,2,nrctmx)
c
      real*8 cbsr(nbt1mx,nsrtmx),cesr(nctmax,nsrtmx),
     $ egers(ietmax,jetmax,nertmx),elecsr(nsrtmx),modr(nrctmx),
     $ mrgers(ietmax,jetmax,nertmx),morr(nrctmx),mwtrc(nrctmx),
     $ rreacn(nrctmx),rreac1(nrctmx),rrelr1(nrctmx),
     $ rxbar(iktmax,nxrtmx),vreac(nrctmx),wodr(nrctmx),
     $ worr(nrctmx),xgers(ietmax,jetmax,nertmx)
c
      real*8 aamatr(kmax,kmax),acflg(nstmax),acflgo(nstmax),
     $ act(nstmax),actlg(nstmax),affp(nptmax),affs(nstmax),
     $ alpha(kmax),amtb(nbtmax),apx(iapxmx,nxtmax),azero(natmax),
     $ a3bars(natmax),beta(kmax),betao(kmax),bpx(ibpxmx,nxtmax),
     $ cco2(5),cegexs(ietmax,jetmax,netmax),cess(nessmx),
     $ cdrs(ndrsmx),cdrsd(ndrsmx),cdrsx(ndrsmx),cdrtw(nstmax),
     $ cdrw(nstmax),cjbasp(nbtmax),cnufac(nstmax),conc(nstmax),
     $ conclg(nstmax),cpgexs(ietmax,jetmax,netmax),cscale(nstmax),
     $ csts(nstsmx),drer0(nrd1mx,nrctmx),delvco(kmax),delvec(kmax),
     $ dlogxw(nbtmax),drir0(nrd1mx),dzvc0(nrd1mx,kmax),d1zvc1(kmax),
     $ egexjc(jetmax,netmax),egexjf(jetmax,netmax),
     $ egexs(ietmax,jetmax,netmax)
c
      real*8 fsort(ngtmax),fugac(ngtmax),fugalg(ngtmax),
     $ gmmatr(kmax,kmax),loph(nptmax),losp(nstmax),lsort(nstmax),
     $ modr0(nrctmx),moph(nptmax),morr0(nrctmx),mosp(nstmax),
     $ mrgexs(ietmax,jetmax,netmax),mtb(nbtmax),mtbaq(nbtmax),
     $ mtb0(nbtmax),mte(nctmax),mteaq(nctmax),ptk(nttkmx),
     $ rhsvec(kmax),rk(imchmx,2,nrctmx),rrelr0(nrctmx),sidrph(nptmax),
     $ sidrsp(nstmax),tempcu(ntprmx),ttk(nttkmx),weight(nstmax),
     $ wfac(iktmax,nxtmax),xbar(nstmax),xbarlg(nstmax),xirct(nrctmx),
     $ xirct0(nrctmx),zchar(nstmax),zchcu6(nstmax),zchsq2(nstmax),
     $ zvclg0(kmax),zvclg1(kmax),zvec0(kmax),zvec1(kmax)
c
      real*8 adh,adhh,adhv,aphi,bdh,bdhh,bdhv,bdot,bdoth,bdotv
c
      real*8 abar,afcnst,al10,avcnst,a3bar,bacfmx,bbig,betamx,bgamx,
     $ bneg,deltim,delxi,dlxmin,eh,ehfac,eps100,farad,fje,fjeo,fo2,
     $ fo2lg,fxi,fxio,omega,omeglg,prcinf,press,pressb,pressd,rcnstv,
     $ rconst,rirec0,rtcnst,screwd,sigmam,sigmmo,smp100,tdays,tempc,
     $ tempcb,tempcd,tempk,tempc0,timemx,time0,time1,tiplol,tiplot,
     $ tiprnl,tiprnt,tistsv,tolbt,toldl,tolsat,tolsst,tolxst,xbarw,
     $ xbarwc,xbrwlc,xbrwlg,xistsv,xi0,xi1,zklogu
c
c-----------------------------------------------------------------------
c
c     Local variable declarations with global dimensioning.
c
c     Variables for restoring the entering configuration, including
c     the entering phase assemblage.
c
      integer isv_kmax,isv_nbtmax,isv_nstmax
c
      SAVE isv_kmax,isv_nbtmax,isv_nstmax
c
      integer, dimension(:), allocatable :: iindxs,ipndxs,nbasps
c
      SAVE iindxs,ipndxs,nbasps
c
      real(8), dimension(:), allocatable :: acflgs,zvclgs
c
      SAVE acflgs,zvclgs
c
c     The following do not need to be SAVEd.
c
      integer kdims,km1s,kmts,kx1s,kxts
c
      real(8) xbarws,xbrwls
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jlen,j2,k,kcol,krow,n,nb,nords,np,nphasl,nplast,np1,
     $ nrdxsl,ns,ns2,ntpr0
c
      integer ilnobl
c
      logical qbswok,qslmod,qsspgb,qztayl
c
      character*56 uspn56
c
      real*8 dlmoph,dltimd,lold,lxx,mold,mxx0
c
      real*8 texp,tlg
c
c-----------------------------------------------------------------------
c
c     Allocate or reallocate local work arrays as needed.
c
      if (.not.ALLOCATED(iindxs)) then
c
c       Local work arrays are not allocated. Zero the saved
c       array size variables. Note that only one array is tested
c       to see if it is allocated. It is assumed that all local
c       work arrays are either allocated or not.
c
        isv_kmax = 0
        isv_nbtmax = 0
        isv_nstmax = 0
      else
c
c       Local work arrays are allocated. Check to see if any of the
c       array size variables have changed. If so, deallocate
c       the corresponding local work arrays and zero the corresponding
c       saved size variables.
c
        if (kmax .ne. isv_kmax) then
          DEALLOCATE(iindxs,ipndxs)
          DEALLOCATE(zvclgs)
          isv_kmax = 0
        endif
c
        if (nbtmax .ne. isv_nbtmax) then
          DEALLOCATE(nbasps)
          isv_nbtmax = 0
        endif
c
        if (nstmax .ne. isv_nstmax) then
          DEALLOCATE(acflgs)
          isv_nstmax = 0
        endif
      endif
c
c     At this point, the saved array size values are zero if the
c     corresponding arrays need to be allocated.
c
      if (isv_kmax .eq. 0) then
        ALLOCATE(iindxs(kmax),ipndxs(kmax))
        ALLOCATE(zvclgs(kmax))
        isv_kmax = kmax
      endif
c
      if (isv_nbtmax .eq. 0) then
        ALLOCATE(nbasps(nbtmax))
        isv_nbtmax = nbtmax
      endif
c
      if (isv_nstmax .eq. 0) then
        ALLOCATE(acflgs(nstmax))
        isv_nstmax = nstmax
      endif
c
c     Zero the contents of the local work arrays.
c
      do k = 1,kmax
        iindxs(k) = 0
        ipndxs(k) = 0
      enddo
c
      do k = 1,kmax
        zvclgs(k) = -99999.
      enddo
c
      do n = 1,nbtmax
        nbasps(n) = 0
      enddo
c
      do n = 1,nstmax
        acflgs(n) = 0.
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize some parameters.
c
      qslmod = .false.
      qmod = .false.
      qbseqc = .false.
c
      ier = 0
      nrdxsl = 0
      nphasl = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Save the "kernel" description for the current equilibrium system.
c     This will be used to recover if the equilibrium calculation
c     overseen by EQ6/eqphas.f fails for any reason. This kernel
c     contains the minimum of information required to calculate the
c     fully expanded description using EQLIB/ncmpex.f. Note that
c     EQ6/eqphas.f (which the present subroutine calls) and EQ6/path.f
c     (which calls the present subroutine) have their own distinct
c     backup kernels.
c
      km1s = km1
      kmts = kmt
      kx1s = kx1
      kxts = kxt
      kdims = kdim
c
      call copyia(nbasp,nbasps,nbt)
      call copyia(iindx1,iindxs,kdim)
      call copyia(ipndx1,ipndxs,kdim)
c
      call copyaa(zvclg1,zvclgs,kdim)
      call copyaa(acflg,acflgs,nst)
c
      xbarws = xbarwc
      xbrwls = xbrwlc
c
c     Save the order of the finite differences.
c
      nords = nord
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  100 continue
c
c     This is a return point if the step size is changed within
c     the present subroutine. There are two cases:
c
c        1. Advance the step size through a region of critical
c           redox instability (redox jump) from the last nonsingular
c           point of reaction progress. Move the simulation over
c           the region of singularity by using only predictor functions
c           expanded from that last good point.
c
c        2. Advance the step size across a region of singularity
c           that is associated with the presence or absence of a phase.
c           There are two cases. first, the calculations may fail to
c           converge upon adding a newly supersaturated phase because
c           the number of moles of the phase that would be produced is
c           so small that the Jacobian matrix becomes effectively
c           singular. Second, when a phase is disappearing, its number
c           of moles may become so small that the above singlularity
c           occurs. This is a normal mode for dropping a phase.
c           Sometimes, the affinity of the phase exceeds the
c           supersaturation tolerance when the system is solved for
c           the case of the phase not present. In each case, the step
c           size must be advanced across the region of apparent
c           singularity.
c
      if (ncorr .le. 0) then
        write (noutpt,1000) xi1,delxi,nord
 1000   format(/' Stepping to Xi= ',1pe11.4,', delxi= ',1pe11.4,
     $  ', nord= ',i1)
      endif
c
      if (iopt(2) .gt. 0) then
        if (xi1 .gt. xistsv) then
          if (.not.qshoot) then
            call timeca(deltim,delxi,drir0,iodb,nodbmx,nord,noutpt,
     $      nrd1mx,nttyo,prcinf,qriinf,rirec0,time0,time1)
c
            if (delxi .le. dlxmin) then
c
c             Make sure that the calculated time does not exceed any
c             any specified limits such as the maximum time just because
c             delxi is at the minimum value.
c
              call tivchk(deltim,delxi,qtvchk,time1,time0,timemx,
     $        tiplol,tiplot,tiprnl,tiprnt,tolxst)
            endif
          endif
        else
          time1 = tistsv
          deltim = 0.
        endif
c
        tdays = time1/86400.
        dltimd = deltim/86400.
c
        write (noutpt,1010) ncorr,tdays,dltimd
 1010   format(3x,'ncorr= ',i2,', time= ',1pe11.4,' d, deltim= ',
     $  e11.4,' d')
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nord .gt. 0) then
c
c       Make a Taylor's series expansion of the dz/d(xi) vector.
c       This information is used to track how the reacting system
c       is changing.
c
        call d1ztay(delxi,dzvc0,d1zvc1,kdim,kmax,nord,nrd1mx)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nrct .ge. 1) then
c
c       Increment the irreversible reactions (those associated with
c       the "reactants"); update the ES mass balance totals accordingly.
c
        call reacts(cbsr,csts,delxi,drer0,iern1,ietmax,iktmax,
     $  iodb,jcode,jetmax,jgext,jreac,modr,modr0,morr,morr0,mrgers,
     $  mtb,mtb0,nbaspd,nbt,nbtmax,nbt1mx,ncmpr,nern1,nern2,nertmx,
     $  netmax,ngext,nodbmx,nord,noutpt,nptmax,nrct,nrctmx,nrd1mx,
     $  nrndex,nsrtmx,nstmax,nsts,nstsmx,nstsr,nttyo,nxridx,nxrtmx,
     $  rrelr0,rxbar,ureac,xirct,xirct0)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set a flag for a go back followed by a step size reduction when
c     the equilibrium calculation indicates a supersaturation that
c     exceeds the tolerance value used in locating the phase boundary
c     for a newly appearing phase.
c
      qsspgb = iopt(3).le.1 .and. .not.qslmod .and.(.not.qscon)
     $ .and. delxi.gt.dlxmin .and. .not.qstart
c
c     Make the equilibrium calculation. Basis switching may occur.
c     Also, the phase assemblage may be changed.
c
      call eqphas(aamatr,abar,acflg,acflgo,act,actlg,
     $ adh,adhh,adhv,afcnst,affp,affs,alpha,al10,amtb,aphi,
     $ apx,avcnst,azero,a3bar,a3bars,bacfmx,bbig,bdh,bdhh,bdhv,
     $ bdot,bdoth,bdotv,beta,betamx,betao,bgamx,bneg,bpx,cco2,
     $ cegexs,cess,cdrs,cdrsd,cdrsx,cdrtw,cdrw,cjbasp,cnufac,
     $ conc,conclg,cpgexs,cscale,csts,delvco,delvec,d1zvc1,dlogxw,
     $ egexjc,egexjf,egexs,eh,ehfac,eps100,farad,fje,fjeo,fo2,
     $ fo2lg,fsort,fugac,fugalg,fxi,fxio,gmmatr,iapxt,ibpxt,ibswx,
     $ ielam,ier,iern1,iern2,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,
     $ ifpsi2,ifrn1,ifrn2,ifzeta,igas,igstak,iindx0,iindx1,ilcphi1,
     $ ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,ilzeta,imrn1,
     $ imrn2,insgf,iodb,iopg,iopt,ipch,ipivot,ipndx1,ipcv,istack,
     $ iter,itermx,ixbasp,ixrn1,ixrn2,izmax,jcsort,jflag,jgsort,
     $ jgstak,jjsort,jpflag,jsflag,jsitex,jsol,jssort,jstack,kbt,
     $ kction,kdim,kelect,khydr,khydx,km1,km10,kmt,kmt0,ko2gaq,
     $ kpsat,kpsst,krdxsp,kwater,kx1,kx10,kxt0,kxt,loph,losp,lsort,
     $ moph,mosp,mrgexs,mtb,mtbaq,narn1,narn2,narxt,nat,nbasp,nbaspd,
     $ nbaspx,nbt,nbtd,nbw,nchlor,ncmpr,nct,ndrs,ndrsd,ndrsx,ndrsr,
     $ ndrsrd,ndrsrx,nelect,nern1,nern2,ness,nessr,net,nfrn1,nfrn2,
     $ ngrn1,ngrn2,ngt,nhydr,nhydx,nlrn1,nlrn2,nlt,nmrn1,nmrn2,
     $ nmt,nord,no2gaq,noutpt,npchk,nphasx,npt,nrdxsp,nst,nsts,
     $ nstsr,ntpr,ntrymx,nttyo,nxrn1,nxrn2,nxt,omega,omeglg,prcinf,
     $ press,qbassw,qbseqc,qbye,qcnpre,qcntmp,qhawep,qmod,qoptmz,
     $ qpit75,qredox,qsspgb,qstart,qxknph,q6mode,rcnstv,rconst,
     $ rhsvec,rtcnst,screwd,sidrph,sidrsp,sigmam,sigmmo,smp100,tempc,
     $ tempk,tolbt,toldl,tolsat,tolsst,ubacmx,ubgamx,ulbeta,uldel,
     $ uphase,uspec,uzvec1,weight,wfac,xbar,xbarlg,xbarw,xbarwc,
     $ xbrwlc,xbrwlg,zchar,zchcu6,zchsq2,zvclg1,zvec1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ier .eq. 8) then
c
c       Go back and reduce the step size to avoid exceeding the
c       dimensioned limit on the iindx1 array (this also corresponds to
c       the dimensions of the Jacobian matrix). EQ6/path.f will go back
c       and cut the step size, if this is possible.
c
        ier = 180
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ier .eq. 10) then
c
c       Go back and reduce the step size to avoid exceeding the
c       supersaturation tolerance (tolsst). This is only done under
c       an appropriate set of conditions, so no further analysis
c       is required. EQ6/path.f will go back and cut the step size.
c
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ier .eq. 20) then
c
c       EQ6/eqphas.f hit the maximum number of tries to find the correct
c       phase without succeeding. Go back and let EQ6/path.f deal with
c       the situation. It may be able to resolve the problem by reducing
c       the step size.
c
        ier = 180
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ier.eq.30 .or. ier.eq.40) then
c
c       The calculation appears to be in a region of critical
c       instability in the ES phase assemblage (can't solve to within
c       specified tolerances with a certain phase or phases in or out
c       of the assemblage). Try to slide over it.
c
        if (nphasl .le. 0) then
c
c         Am not currently engaged in such a slide. Need to start
c         the process.
c
          if (npslmx .le. 0) then
c
c           This kind of sliding is inhibited.
c
            write (noutpt,1100)
            write (nttyo,1100)
 1100       format(/' * Warning - (EQ6/eqshel) Need to slide Xi',
     $      ' forward to get over a region of',/7x,'critical',
     $      ' instability in the ES phase assemblage. However, the',
     $      ' maximum',/7x,'number of steps for this process (npslmx)',
     $      ' is zero. Sliding of this kind',/7x,'is therefore',
     $      ' inhibited.')
            ier = 190
            go to 999
          else
c
c           Set up to start sliding.
c
            write (noutpt,1110)
 1110       format(/' Setting up to slide Xi forward to get over a',
     $      ' region of critical'/3x,'instability in the ES phase',
     $      ' assemblage.')
            qslmod = .true.
          endif
        endif
c
        if (nphasl .ge. npslmx) then
c
c         Check for the maximum number of tries sliding forward.
c
          write (noutpt,1120)
          write (nttyo,1120)
 1120     format(/' * Warning - (EQ6/eqshel) Have done the maximum',
     $    ' number of tries to',/7x,'slide over a region of critical',
     $    ' instability in the ES phase assemblage.')
          ier = 190
          go to 999
        endif
c
c       Slide Xi forward one increment. The increment gets progressively
c       larger with the number of tries.
c
        nphasl = nphasl + 1
        delxi = delxi + 3.**(nphasl - 1)*dlxmin
        write (noutpt,1130)
 1130   format(/' Trying to slide over a region of critical',
     $  ' instability',/7x,'in the ES phase assemblage:')
        write (noutpt,1140) nphasl,delxi
 1140   format(/3x,'Try ',i2,': delxi= ',1pe12.5,/)
        go to 300
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ier .eq. 50) then
c
c       The calculation appears to be in a region of critical redox
c       instability. Try to slide over it.
c
        if (nrdxsl .le. 0) then
c
c         Am not currently engaged in such a slide. Need to start
c         the process.
c
          if (nsslmx .le. 0) then
c
c           This kind of sliding is inhibited.
c
            write (noutpt,1200)
            write (nttyo,1200)
 1200       format(/' * Warning - (EQ6/eqshel) Need to slide Xi',
     $      ' forward to get over a region of',/7x,'critical redox',
     $      ' instability. However, the maximum number of steps',
     $      /7x,'for this process (nsslmx) is zero. Sliding of this',
     $      ' kind is',/7x,'therefore inhibited.')
            ier = 190
            go to 999
          else
c
c           Set up to start sliding.
c
            write (noutpt,1210)
 1210       format(/' Setting up to slide Xi forward to get over a',
     $      ' region of critical'/3x,'redox instability.')
            qslmod = .true.
          endif
        endif
c
        if (nrdxsl .ge. nsslmx) then
c
c         Check for the maximum number of tries sliding forward.
c
          write (noutpt,1220)
          write (nttyo,1220)
 1220     format(/' * Warning - (EQ6/eqshel) Have done the maximum',
     $    ' number of tries',/7x,'to slide over a region of critical',
     $    ' redox instability.')
          ier = 190
          go to 999
        endif
c
c       Slide Xi forward one increment. The increment gets progressively
c       larger with the number of tries.
c
        nrdxsl = nrdxsl + 1
        delxi = delxi + 3.**(nrdxsl - 1)*dlxmin
        write (noutpt,1230)
 1230   format(/' Trying to slide over a region of critical redox',
     $  ' instability:')
        write (noutpt,1240) nrdxsl,delxi
 1240   format(/3x,'Try ',i2,': delxi= ',1pe12.5,/)
        go to 300
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ier.eq.60 .or. ier.eq.70) then
c
c       The calculation appears to be in a region of critical
c       instability associated with solvent water. Go back and let
c       EQ6/path.f deal with the situation. It may be able to resolve
c       the problem by reducing the step size.
c
        ier = 180
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ier .eq. 80) then
c
c       EQ6/eqphas.f detected critical instability associated with the
c       aqueous activity coefficient model. Go back and let EQ6/path.f
c       deal with the situation. It may be able to at least partially
c       resolve the problem by reducing the step size. However, it is
c       likely that the reaction path is causing the activity
c       coefficient model to be extrapolated outside its range of
c       validity.
c
        ier = 180
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ier .eq. 100) then
c
c       The calculation was about to start with out-of-range values
c       for the iteration variables associated with the aqueous basis
c       species. Go back and let EQ6/path.f deal with the situation.
c       It may be able to resolve the problem by reducing the step size.
c       Alternatively, there may be some kind of programming error.
c
        ier = 180
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ier .eq. 110) then
c
c       The iteration ended with out-of-range values for the iteration
c       variables associated with the aqueous basis species. Go back and
c       let EQ6/path.f deal with the situation. It may be able to
c       resolve the problem by reducing the step size. Alternatively,
c       there may be some kind of programming error.
c
        ier = 180
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ier .eq. 150) then
c
c       The iteration ended without any useful diagnostics being
c       produced. Go back and let EQ6/path.f deal with the situation.
c       It may be able to resolve the problem by reducing the step size.
c       Alternatively, there may be some kind of programming error.
c
        ier = 180
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ier .eq. 0) then
        if (iopt(1) .eq. 2) then
c
c         Check to see that no significant numbers of moles of solids
c         were destroyed unexpectedly if computing a fluid-centered
c         flow-through open system model.
c
          do kcol = km1s,kmts
            if ( uzvec0(kcol)(1:5) .ne. ufixf(1:5) ) then
              if (zvclg0(kcol) .le. zklogu) then
                ns = iindxs(kcol)
                lxx = losp(ns)
                dlmoph = zvec0(kcol) - texp(lxx)
                if (dlmoph .gt. 0.) then
                  if (tlg(dlmoph) .gt. zklogu) then
                    ier = 170
                    mxx0 = zvec0(kcol)
c
c                   Calling sequence substitutions:
c                     uspec(ns) for unam48
c
                    call fmspnm(jlen,uspec(ns),uspn56)
                    write (noutpt,1300) uspn56(1:jlen),mxx0,mosp(ns),
     $              dlmoph
 1300               format(' Some of ',a,' was unexpectedly destroyed.',
     $              /5x,'Previous number of moles was ',1pe12.5,
     $              /5x,'Current number of moles is   ',e12.5,
     $              /5x,'Amount destroyed was ',e12.5,/21x,
     $              /3x,"That's too much. Will go back and first",
     $              ' transfer some of the current',/3x,'amount to the',
     $              ' physically removed system (PRS).')
                  endif
                endif
              endif
            endif
          enddo
c
          nplast = 0
c
          do kcol = kx1s,kxts
            ns = iindxs(kcol)
            np = ipndxs(kcol)
            if (np .ne. nplast) then
              nplast = np
              mold = 0.
c
              do krow = kcol,kxts
                np1 = ipndxs(krow)
                if (np1 .ne. np) go to 200
                mxx0 = zvec0(krow)
                mold = mold + mxx0
              enddo
  200         continue
c
              lold = tlg(mold)
              if (lold .gt. zklogu) then
                dlmoph = mold - moph(np)
                if (dlmoph .gt. 0.) then
                  if (tlg(dlmoph) .gt. zklogu) then
                    ier = 170
                    j2 = ilnobl(uphase(np))
                    write (noutpt,1300) uphase(np)(1:j2),mold,moph(np),
     $              dlmoph
                  endif
                endif
              endif
            endif
          enddo
          go to 999
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The calculation finished with no problems.
c
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  300 continue
c
c     Restore the kernel for the equilibrium sytem from the backup.
c     Undo any basis switches that may have been made by EQ6/optmzr.f.
c
c     Note: some aspects of the backup kernel are not needed here.
c     For example, the zvclg1 array will be recomputed further below
c     by the call to EQ6/ztaylr.f, overwriting the values restored from
c     the zvclgs array. However, a full restoration is done here for the
c     sake of completeness.
c
      if (qbseqc) then
        do kcol = 1,kbt
          nb = iindx1(kcol)
          ns = nbasp(nb)
          ns2 = nbasps(nb)
          if (ns2 .ne. ns) then
            call switch(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,
     $      axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,eps100,ipch,ipchmx,
     $      ipcv,ipcvmx,jflag,jsflag,narn1,narxmx,nbasp,nbaspd,nbaspx,
     $      nb,nbt,nbtmax,nbw,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,noutpt,
     $      ns2,nst,nstmax,ntprmx,nttyo,qbassw,qbswok,uspec)
          endif
        enddo
        qbseqc = .false.
        nord = nords
      endif
c
      km1 = km1s
      kmt = kmts
      kx1 = kx1s
      kxt = kxts
      kdim = kdims
c
      call copyia(nbasps,nbasp,nbt)
      call copyia(iindxs,iindx1,kdim)
      call copyia(ipndxs,ipndx1,kdim)
c
      call copyaa(zvclgs,zvclg1,kdim)
      call copyaa(acflgs,acflg,nst)
c
      xbarwc = xbarws
      xbrwlc = xbrwls
c
c     Reset the ixbasp and cjbasp arrays. The former is a flag
c     array, each member of which denotes whether the
c     thermodynamic activity of the corresponding basis species
c     is defined in terms of molality (= 0) or mole fraction (= 1).
c     The cjbasp array contains any site stoichiometric factors
c     associated with the operational basis species.
c
      call gibasp(cgexj,cjbasp,iern1,ixbasp,jern1,jern2,
     $ jetmax,jgext,narn1,narn2,nbasp,nbt,nbtmax,nern1,nern2,
     $ netmax,nphasx,nstmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Slightly increase the value of reaction progress to try to slide
c     over a region of critical instability.
c
      xi1 = xi0 + delxi
c
c     Make a Taylor's series expansion of the z vector, applying
c     change limits. This provides a protected set of values for
c     the AE solver; i.e., a set not containing any values that
c     will cause the solver to fail in an unrecoverable fashion.
c
      qztayl = .true.
      call ztaylr(delxi,dzvc0,kdim,kmax,km1,kxt,nord,nrd1mx,
     $ qztayl,zklogu,zvclg0,zvclg1,zvec0,zvec1)
c
c     Save the new z vector expansion.
c
      call copyaa(zvclg1,zvclgs,kdim)
c
      if (.not.qcntmp .or. .not.qcnpre) then
c
c       Recompute the temperature and pressure. Then recompute the
c       thermodynamic and kinetic quantities which depend these
c       variables.
c
        ntpr0 = ntpr
        call tpadv(aadh,aadhh,aadhv,aaphi,abdh,abdhh,abdhv,
     $  abdoth,abdot,abdotv,adadhh,adadhv,adbdhh,adbdhv,adbdth,adbdtv,
     $  adh,adhfe,adhh,adhv,adhfs,adhfsd,advfe,advfs,advfsd,afcnst,
     $  al10,amu,aslm,aphi,aprehw,apresg,apresh,apx,avcnst,axhfe,axhfs,
     $  axhfsd,axlke,axlks,axlksd,axvfe,axvfs,axvfsd,bdh,bdhh,bdhv,
     $  bdot,bdoth,bdotv,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,dhfe,dhfs,
     $  dhfsd,dvfe,dvfs,dvfsd,eact,ehfac,farad,hact,iact,iapxmx,iktmax,
     $  imchmx,imech,iopg,iopt,ipbtmx,ipch,ipchmx,ipcv,ipcvmx,ixrn1,
     $  ixrn2,jpfcmx,jpress,jptffl,jsol,jtemp,narxmx,narxt,narxth,nbasp,
     $  nbaspd,nbt,nbtd,nbtmax,ncmpr,ndrsr,ndrsrd,nmut,nmutmx,nopgmx,
     $  noptmx,noutpt,nptkmx,nptmax,nrct,nrctmx,nrk,nslt,nsltmx,nst,
     $  nstmax,ntpr,ntprmx,ntprt,nttkmx,nttyo,nweope,nwndpc,nxt,nxtmax,
     $  pmu,presg,presh,press,pressb,pressd,pslamn,ptk,rcnstv,rconst,rk,
     $  rkb,rtcnst,tempc,tempcb,tempcd,tempcu,tempk,time1,trkb,ttk,
     $  uphase,uspec,wfac,xhfe,xhfs,xhfsd,xi1,xlke,xlks,xlksd,xvfe,
     $  xvfs,xvfsd)
c
        qtrch = ntpr .ne. ntpr0
      endif
      go to 100
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
