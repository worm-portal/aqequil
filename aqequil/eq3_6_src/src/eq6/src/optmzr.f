      subroutine optmzr(aamatr,abar,acflg,acflgo,act,actlg,
     $ adh,adhh,adhv,afcnst,al10,alpha,amtb,aphi,avcnst,azero,
     $ a3bar,a3bars,bacfmx,bbig,bdh,bdhh,bdhv,bdot,bdoth,bdotv,
     $ beta,betamx,bgamx,bneg,bpx,cco2,cdrs,cdrsx,cdrtw,cdrw,
     $ cegexs,cgexj,cjbasp,cnufac,conc,conclg,cpgexs,cscale,
     $ csts,coval,delvec,dlogxw,egexjc,egexjf,egexs,ehfac,
     $ eps100,fje,fjeo,fo2,fo2lg,fsort,fugac,fugalg,fxi,fxio,
     $ gmmatr,ibetmx,ibpxt,ibswx,ielam,iern1,iern2,ifcphi1,
     $ ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,
     $ igas,igstak,iindx1,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,
     $ ilpsi2,ilrn1,ilrn2,ilzeta,imrn1,imrn2,insgf,iodb,iopg,
     $ iopt,ipch,ipcv,ipivot,ipndx1,irdxc3,istack,ixbasp,ixrn1,
     $ ixrn2,izmax,jcsort,jern1,jern2,jflag,jgext,jgsort,jgstak,
     $ jjsort,jpflag,jsflag,jsitex,jsol,jssort,jstack,kbt,kction,
     $ kdim,kelect,khydr,km1,kmt,ko2gaq,kwater,kx1,kxt,loph,losp,
     $ lsort,mgext,moph,mosp,mrgexs,mtb,narn1,narn2,narxt,nat,
     $ nbasp,nbaspd,nbaspx,nbw,nbt,nbtd,nchlor,ncmpr,ncosp,ndrs,
     $ ndrsx,ndrsr,ndrsrd,ndrsrx,nelect,nern1,nern2,net,ngexsa,
     $ nfac,ngext,nhydr,nhydx,ngrn1,ngrn2,ngt,noutpt,no2gaq,
     $ nphasx,npt,nst,nsts,nstsr,ntfx,ntfxt,ntpr,nttyo,omega,
     $ omeglg,press,qbassw,qblamx,qhawep,qpit75,qredox,q6mode,
     $ rhsvec,sigmam,sigmmo,smp100,tempc,tempk,tfx,tolbig,tolneg,
     $ tolxpt,ubacmx,ubbig,ubetmx,ubgamx,ubneg,ugexj,ugexmo,
     $ uphase,uspec,uzvec1,weight,wfac,xbar,xbarlg,xbarw,xbarwc,
     $ xbrwlc,xbrwlg,zchar,zchcu6,zchsq2,zgexj,zvclg1,zvec1)
c
c     This subroutine optimizes the starting values prior to hybrid
c     Newton-Raphson iteration. The structure is somewhat similar to
c     that in EQ3NR/arrset.f. Automatic basis switching is done in
c     "loops." Inside loops are "passes," in which the activity
c     coefficients are readjusted. Inside passes are "cycles," in which
c     the chief iteration variables are adjusted.
c
c     In the present subroutine, the cycle algorithm is based on
c     minimizing the function:
c
c       aleph = Sum(i) zeta(i)**2
c
c     where i spans the set of primary iteration variables and zeta
c     is a residual vector derived from the beta residual vectors
c     normally used in EQ3/6. In EQ3NR/arrset.f, the cycle algorithm
c     is the continued-fraction algorithm. That algorithm is not
c     well suited to dealing with nonpphysical mass balances,
c     which must be dealt with for basis species like H2O(l), H+,
c     O2(g,aq), and e-. Nor is it well-suited for dealing with
c     non-aqueous species (e.g., pure minerals) when they are not
c     switched into the active basis set but are instead treated as
c     "extended" basis species. Hence the change of algorithm.
c
c     This subroutine does not change the presumed phase assemblage.
c
c     This subroutine is called by:
c
c       EQ6/eqcalc.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c
c     Principal output:
c
c       nhydr  = species index of the hydrogen ion
c       nhydx  = species index of the hydroxide ion
c       nchlor = species index of the chloride ion
c       qblamx = .true. if solid solutions are in the matrix
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
      include 'eqlib/eqldv.h'
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer noutpt,nttyo
c
      integer ibpxt(nxtmax),ibswx(nbtmax),igstak(ngtmax),iindx1(kmax),
     $ insgf(natmax),iodb(nodbmx),iopg(nopgmx),iopt(noptmx),
     $ ipivot(kmax),ipndx1(kmax),istack(nstmax),ixbasp(nbtmax),
     $ jcsort(nstmax),jern1(jetmax,netmax),jern2(jetmax,netmax),
     $ jflag(nstmax),jgext(netmax),jgsort(ngtmax),jgstak(ngtmax),
     $ jjsort(nstmax),jpflag(nptmax),jsflag(nstmax),jsitex(nstmax),
     $ jsol(nxtmax),jssort(nstmax),jstack(nstmax),kction(nbtmax),
     $ narxt(ntprmx),nbasp(nbtmax),nbaspd(nbtmax),nbaspx(nbtmax),
     $ ncmpr(2,nptmax),ncosp(nbtmax),ndrs(ndrsmx),ndrsx(ndrsmx),
     $ ndrsr(2,nstmax),ndrsrd(2,nstmax),ndrsrx(2,nstmax),
     $ ngexsa(ietmax,jetmax,netmax),ngext(jetmax,netmax),nphasx(nstmax),
     $ nsts(nstsmx),nstsr(2,nstmax),ntfx(ntfxmx)
c
      integer nfac(nbtmax)
c
      integer ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,
     $ ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta
c
      integer ibetmx,ielam,iern1,iern2,ifrn1,ifrn2,igas,ilrn1,ilrn2,
     $ imrn1,imrn2,ipch,ipcv,irdxc3,ixrn1,ixrn2,izmax,kbt,kdim,kelect,
     $ khydr,km1,kmt,ko2gaq,kwater,kx1,kxt,narn1,narn2,nat,nbw,nbt,
     $ nbtd,nchlor,nelect,nern1,nern2,net,nhydr,nhydx,ngrn1,ngrn2,ngt,
     $ no2gaq,npt,nst,ntpr,ntfxt
c
      logical qabsw,qbassw,qblamx,qhawep,qpit75,qredox,q6mode
c
      character*48 uspec(nstmax),uzvec1(kmax)
      character*48 ubacmx,ubbig,ubetmx,ubgamx,ubneg
      character*24 ugexmo(netmax),uphase(nptmax)
      character*8 ugexj(jetmax,netmax)
c
      real*8 aamatr(kmax,kmax),acflg(nstmax),acflgo(nstmax),act(nstmax),
     $ actlg(nstmax),alpha(kmax),amtb(nbtmax),azero(natmax),
     $ a3bars(natmax),beta(kmax),bpx(ibpxmx,nxtmax),cco2(5),
     $ cdrs(ndrsmx),cdrsx(ndrsmx),cdrtw(nstmax),cdrw(nstmax),
     $ cegexs(ietmax,jetmax,netmax),cgexj(jetmax,netmax),
     $ cjbasp(nbtmax),cnufac(nstmax),conc(nstmax),conclg(nstmax),
     $ coval(nbtmax),cpgexs(ietmax,jetmax,netmax),cscale(nstmax),
     $ csts(nstsmx)
c
      real*8 delvec(kmax),dlogxw(nbtmax),egexjc(jetmax,netmax),
     $ egexjf(jetmax,netmax),egexs(ietmax,jetmax,netmax),
     $ fsort(ngtmax),fugac(ngtmax),fugalg(ngtmax),gmmatr(kmax,kmax),
     $ mgext(jetmax,netmax),mrgexs(ietmax,jetmax,netmax),
     $ mtb(nbtmax),loph(nptmax),losp(nstmax),lsort(nstmax),
     $ moph(nptmax),mosp(nstmax),rhsvec(kmax),tfx(ntfxmx),
     $ weight(nstmax),wfac(iktmax,nxtmax),xbar(nstmax),xbarlg(nstmax),
     $ zchar(nstmax),zchcu6(nstmax),zchsq2(nstmax),
     $ zgexj(jetmax,netmax),zvclg1(kmax),zvec1(kmax)
c
      real*8 adh,adhh,adhv,aphi,bdh,bdhh,bdhv,bdot,bdoth,bdotv
c
      real*8 abar,afcnst,al10,avcnst,a3bar,bacfmx,bbig,betamx,bgamx,
     $ bneg,ehfac,eps100,fje,fjeo,fo2,fo2lg,fxi,fxio,omega,omeglg,press,
     $ sigmam,sigmmo,smp100,tempc,tempk,tolbig,tolneg,tolxpt,xbarw,
     $ xbarwc,xbrwlc,xbrwlg,zetsqm,zetsqo
c
c-----------------------------------------------------------------------
c
c     Local variable declarations with global dimensioning.
c
      integer, dimension(:), allocatable :: iastak,jasort,jastak
      integer, dimension(:), allocatable :: kmvar,knflag
c
      real(8), dimension(:), allocatable :: daleph,delprc,dijmaj
      real(8), dimension(:), allocatable :: zesort,zeta,zetasq
c
      real(8), dimension(:), allocatable :: efac
c
      integer isv_nbtmax,isv_kmax
c
      SAVE isv_nbtmax,isv_kmax
c
      SAVE daleph,delprc,dijmaj,efac,iastak,jasort,jastak,kmvar,knflag,
     $ zesort,zeta,zetasq
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer icorr,kdirec,iebal,iter,j,jdim,jlen,jlen1,jlen2,j2,k,
     $ kbig,kcol,kcorr,kcscan,ker,kk,kkdim,kount,krow,krscan,n,nb,
     $ nbigger,nchange,ncycle,ncylim,negafc,nloop,nlopmx,npass,nplim,
     $ nredox,ns,nscan,nsd,nswtch,ns1,ns2
c
      integer ilnobl
c
      logical qbswx,qcgam,qcfxi,qcsigm,qcycnc,qpracf,qscan2,qsclim,
     $ qscosc,qstops,qtestc,qtestp,qxbarw
c
      character*56 uspn56,usp156,usp256
      character*48 uzebig,uzeneg
      character*24 ujtp
c
      real*8 abig,actwlc,adij,alefnc,alepoe,alepoo,aleph,alephm,av,ax,
     $ azex,betfnc,bfje,bfxi,bsigmm,btmxoe,btmxoo,bx,chfacf,chfsgm,cx,
     $ dij,dscan,dscmin,dx,dzx,dzxclm,eh,eps1hi,rlxgam,tolatf,tolgpt,
     $ zebig,zeneg,zetamx,zetsum,zex,zvclim,zvclgi,zvclgo,zvcnew,zx,
     $ zx1,zx2
c
      real*8 texp
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
c
      data nlopmx/12/,nplim  /7/,ncylim /15/
c
c     The following are local tolerance parameters for the
c     optimization:
c
      data tolatf /0.005/,tolgpt /0.1/
c
c     The following is the minimum scan increment.
c
      data dscmin /0.05/
c
c     The following are not used by EQ6, but are needed to satisfy
c     some subroutine calls.
c
      data iebal  /0/,nredox /0/
      data eh     /0./
c
c     The following is set to .false. to turn off iterative improvement
c     of xbarw in ncmpex.f.
c
      data qxbarw/.false./
c
c-----------------------------------------------------------------------
c
c     Allocate or reallocate local work arrays as needed.
c
      if (.not.ALLOCATED(iastak)) then
c
c       Local work arrays are not allocated. Zero the saved
c       array size variables. Note that only one array is tested
c       to see if it is allocated. It is assumed that all local
c       work arrays are either allocated or not.
c
        isv_kmax = 0
        isv_nbtmax = 0
      else
c
c       Local work arrays are allocated. Check to see if any of the
c       array size variables have changed. If so, deallocate
c       the corresponding local work arrays and zero the corresponding
c       saved size variables.
c
        if (kmax .ne. isv_kmax) then
          DEALLOCATE(iastak,jasort,jastak)
          DEALLOCATE(kmvar,knflag)
          DEALLOCATE(daleph,delprc,dijmaj)
          DEALLOCATE(zesort,zeta,zetasq)
          isv_kmax = 0
        endif
c
        if (nbtmax .ne. isv_nbtmax) then
          DEALLOCATE(efac)
          isv_nbtmax = 0
        endif
      endif
c
c     At this point, the saved array size values are zero if the
c     corresponding arrays need to be allocated.
c
      if (isv_kmax .eq. 0) then
        ALLOCATE(iastak(kmax),jasort(kmax),jastak(kmax))
        ALLOCATE(kmvar(kmax),knflag(kmax))
        ALLOCATE(daleph(kmax),delprc(kmax),dijmaj(kmax))
        ALLOCATE(zesort(kmax),zeta(kmax),zetasq(kmax))
        isv_kmax = kmax
      endif
c
      if (isv_nbtmax .eq. 0) then
        ALLOCATE(efac(nbtmax))
        isv_nbtmax = nbtmax
      endif
c
c     Zero the contents of the local work arrays.
c
      do k = 1,kmax
        iastak(k) = 0
        jasort(k) = 0
        jastak(k) = 0
        kmvar(k) = 0
        knflag(k) = 0
      enddo
c
      do k = 1,kmax
        daleph(k) = 0.
        delprc(k) = 0.
        dijmaj(k) = 0.
        zesort(k) = 0.
        zeta(k) = 0.
        zetasq(k) = 0.
      enddo
c
      do n = 1,nbtmax
        efac(n) = 0.
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(3) .ge. 1) then
        write (noutpt,1000)
 1000   format(/' Starting Pre-Newton-Raphson Optimization.',/)
      endif
c
      qbswx = .false.
      eps1hi = 1.0/eps100
c
      bsigmm = 0.
      bfxi = 0.
      bgamx = 0.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(3) .ge. 3) then
c
c       Print the active data file basis set.
c
        write (noutpt,1010)
 1010   format(/16x,'--- Active Data File Basis Set ---',
     $  //2x,'krow   Name',30x,'Constraint',/)
c
        do krow = 1,kdim
          if (krow .le. kbt) then
            ujtp = 'Mass balance, moles'
          elseif (krow.ge.km1 .and. krow.le.kxt) then
            ujtp = 'Mass action'
          else
            ujtp = 'Defining equation'
          endif
          j2 = ilnobl(ujtp)
c
          if (krow .le. kbt) then
            nb = iindx1(krow)
            ns = nbaspd(nb)
          else
            ns = iindx1(krow)
          endif
c
c         Calling sequence substitutions:
c           uspec(ns) for unam48
c
          call fmspnx(jlen,uspec(ns),uspn56)
c
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
  200 nloop = nloop + 1
      if (iodb(3) .ge. 1) write (noutpt,1040) nloop
 1040 format(6x,'nloop= ',i2)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(3) .ge. 3) then
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
c
          if (kxt .gt. kbt) write (noutpt,1070)
 1070     format(16x,'--- Extended Basis Species ---',
     $    //2x,'krow   Species',/)
c
          do krow = km1,kxt
            ns = iindx1(krow)
            call fmspnx(jlen,uspec(ns),uspn56)
            jlen = min(jlen,32)
            write (noutpt,1080) krow,uspn56(1:jlen)
 1080       format(1x,i4,2x,a)
          enddo
c
          write (noutpt,1030)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Note: At this point, do not recalculate values of the following:
c
c       The SUM(i) m(i) function (sigmam)
c       The ionic strength (fxi)
c       The J electrostatic moment function (fje)
c       The activity coefficients of aqueous species
c       The mole fraction of water
c       The activity coefficients of exchanger species
c       The activity coefficients of solid solution components
c         present in the equilibrium system.
c
c     This subroutine is always called with initial values for all of
c     these parameters. At the start of an EQ6 run, all these parameters
c     are initialized by EQ6/exivar.f. After that, new values are
c     subsequently generated by normal stepping procedures in
c     EQ6/path.f.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Here npass is the pass counter.
c
      npass = -1
c
c     The label below is a return point for subsequent passes. A pass
c     is an adjustment for the ionic strength, etc., the activity of
c     water, and the activity coefficients of the solute species.
c
  210 npass = npass + 1
c
c     Note:
c
c       alefnc = aleph convergence function
c       negafc = the number of successive iterations that the
c                convergence function alefnc
c       betfnc = beta convergence function
c
      alefnc = 0.
      alepoe = 0.
      alepoo = 0.
      negafc = 0
c
      betfnc = 0.
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
      ncycle = -1
c
      qcsigm =abs(bsigmm) .le. tolgpt
      qcfxi = abs(bfxi) .le. tolgpt
      qcgam = bgamx .le. tolgpt
c
      qtestp = qcsigm .and. qcfxi .and. qcgam
c
      do krow = 1,kdim
        delprc(krow) = 0.
      enddo
c
      if (iodb(3) .ge. 1) write (noutpt,1190)
 1190 format(/14x,'Beginning cycle corrections',/)
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
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Expand the description of the equilibrium system from the
c     extended basis variables.
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
c     Compute the usual residuals (the alpha vector and beta vectors).
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
      do krow = 1,kdim
        bx = beta(krow)
        if (krow .le. kbt) then
c
c         Mass balance element.
c
          if (bx .ge. 0.0) then
c
c           Have a non-negative beta(i) value. Set zeta(i) = beta(i).
c
            zeta(krow) = bx
          else
c
c           Have a negative beta(i) value.
c
            dx = 1.0 + bx
            if (dx .gt. eps100) then
c
c             Have (-1.0 + eps100) < beta(i) < 0.0.
c             Set zeta(i) = 1.0 - 1.0/(1.0 + beta(i)).
c
              zeta(krow) = 1.0 - (1.0/dx)
            else
c
c             Have beta(i) < (-1.0 + eps100).
c             Set zeta(i) = 1.0 - (1.0/eps100).
c
              zeta(krow) = 1.0 - eps1hi
            endif
          endif
        else
c
c         Mass action element.
c
c         zeta(i) = 0.001*beta(i)/cscale(i).
c
          ns = iindx1(krow)
          zeta(krow) = 0.001*bx/cscale(ns)
        endif
      enddo
c
c     Calculate the aleph and zeta-squared functions.
c
      aleph = 0.
      do krow = 1,kdim
        zetasq(krow) = zeta(krow)*zeta(krow)
        aleph = aleph + zetasq(krow)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the aleph convergence function.
c
      alefnc = 0.
      if (mod(ncycle,2) .eq. 0) then
        if (alepoe .ge. smp100) alefnc = (alepoe - aleph)/alepoe
        alepoe = aleph
      else
        if (alepoo .ge. smp100) alefnc = (alepoo - aleph)/alepoo
        alepoo = aleph
      endif
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
      if (iodb(3) .ge. 3) then
        write (noutpt,1230)
 1230   format(//10x,'--- Pre-Newton-Raphson Optimization Summary ---',
     $  //2x,'kcol   Name',32x,'zvclg1      zvec1',/)
c
        do kcol = 1,kbt
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
        do kcol = km1,kxt
          ns = iindx1(kcol)
          zx1 = zvclg1(kcol)
          zx2 = texp(zx1)
c
c         Calling sequence substitutions:
c           uspec(ns) for unam48
c
          call fmspnx(jlen,uspec(ns),uspn56)
          write (noutpt,1240) kcol,uspn56,zx1,zx2
        enddo
c
        write (noutpt,1250)
 1250   format(/1x,'krow   Name',28x,'Alpha',9x,'Beta',9x,'Zeta',/)
c
        do krow = 1,kbt
          nb = iindx1(krow)
          ns = nbaspd(nb)
c
c         Calling sequence substitutions:
c           uspec(ns) for unam48
c
          call fmspnx(jlen,uspec(ns),uspn56)
          write (noutpt,1260) krow,uspn56,alpha(krow),beta(krow),
     $    zeta(krow)
 1260     format(1x,i3,2x,a28,2x,1pe12.5,2x,e12.5,2x,e12.5)
        enddo
c
        do krow = km1,kxt
          ns = iindx1(krow)
c
c         Calling sequence substitutions:
c           uspec(ns) for unam48
c
          call fmspnx(jlen,uspec(ns),uspn56)
          write (noutpt,1260) krow,uspn56,alpha(krow),beta(krow),
     $    zeta(krow)
        enddo
        write (noutpt,1030)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(3) .ge. 1) then
c
c       Calculate some data on the state of the zeta residual
c       functions.
c
        zetamx = 0.
        zebig = 0
        uzebig = 'None'
        zeneg = 0
        uzeneg = 'None'
        do kcol = 1,kdim
          zex = zeta(kcol)
          if (zex .gt. zebig) then
            zebig = zex
            uzebig = uzvec1(kcol)
          endif
          if (zex .lt. zeneg) then
            zeneg = zex
            uzeneg = uzvec1(kcol)
          endif
          azex = abs(zex)
          if (azex .gt. zetamx) then
            zetamx = azex
          endif
        enddo
c
c       Write some data on the state of the aleph and beta residual
c       functions.
c
        write (noutpt,1270) aleph,alefnc
 1270   format(/19x,'aleph= ',1pe12.5,', alefnc= ',e12.5)
c
c       Write some data on the state of the beta residual functions.
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
c         usp246 for uspn56
c
        call fmspnx(jlen2,ubneg,usp256)
c
        write (noutpt,1280) betamx,betfnc,bbig,usp156(1:jlen1),
     $  bneg,usp256(1:jlen2)
 1280   format(/18x,'betamx= ',1pe12.5,', betfnc= ',e12.5,
     $  /18x,'  bbig= ',1pe12.5,', ubbig= ',a,
     $  /18x,'  bneg= ',1pe12.5,', ubneg= ',a)
c
c       Write some data on the state of the zeta residual functions.
c
c       Calling sequence substitutions:
c         jlen1 for jlen
c         uzebig for unam48
c         usp156 for uspn56
c
        call fmspnx(jlen1,uzebig,usp156)
c
c       Calling sequence substitutions:
c         jlen2 for jlen
c         uzeneg for unam48
c         usp256 for uspn56
c
        call fmspnx(jlen2,uzeneg,usp256)
c
        write (noutpt,1290) zetamx,zebig,usp156(1:jlen1),
     $  zeneg,usp256(1:jlen2)
 1290   format(/18x,'zetamx= ',1pe12.5,
     $  /18x,' zebig= ',1pe12.5,', uzebig= ',a,
     $  /18x,' zeneg= ',1pe12.5,', uzeneg= ',a,/)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Test the mass balance residuals see if another cycle should be
c     made before attempting to make an improved estimate of the
c     ionic strength.
c
      qtestc = betamx.le.tolxpt .and. bbig.le.tolbig .and.
     $  bneg.ge.tolneg
c
      if (ncycle.gt.2 .and. alefnc.le.tolatf) then
        negafc = negafc + 1
      else
        negafc = 0
      endif
c
      qcycnc = negafc .ge. 4
c
c     Quit doing cycles if:
c
c       1. The cycle convergence criteria are met.
c       2. The maximum number of cycles have been done.
c       3. The convergence function alefnc indicates that the cycles
c            are not converging.
c
      if (qtestc) go to 400
      if (ncycle .ge. ncylim) go to 400
      if (qcycnc) then
        if (iodb(3) .ge. 1) then
          write (noutpt,1292)
 1292     format(11x,'The cycles are not converging.',/)
        endif
        go to 400
       endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Identify the dominant species in each mass balance and
c     compute the corresponding exponent for a continued
c     fraction correction.
c
c     Presently the species determined as the dominant can only be an
c     aqueous species or an exchanger species. Minerals and solid
c     solution component species are ignored.
c
      call cfracf(cdrs,csts,efac,jcsort,jflag,jssort,kmax,mosp,
     $ narn1,narn2,nbasp,nbaspd,nbt,nbtmax,ndrs,ndrsmx,ndrsr,nern1,
     $ nern2,nfac,nst,nstmax,nsts,nstsmx,nstsr,q6mode,weight)
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
          nsd = nfac(nb)
          if (nsd.ne.0 .and. nsd.ne.ns) kount = kount + 1
        enddo
c
        if (kount .gt. 0) then
          write (noutpt,1300)
 1300     format(16x,'--- Mass Balance Dominants ---',/)
          write (noutpt,1310)
 1310     format(4x,'Master Species',21x,'Dominant Species',/)
c
          do kcol = 1,kbt
            nb = iindx1(kcol)
            ns = nbasp(nb)
            nsd = nfac(nb)
            if (nsd.ne.0 .and. nsd.ne.ns) then
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
c               uspec(nsd) for unam48
c               usp256 for uspn56
c
              call fmspnx(jlen2,uspec(nsd),usp256)
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
      if (qabsw .and. nloop.lt.nlopmx) then
c
c       In automatic basis switching mode (iopt(11) .ge. 1), try to
c       first reduce the magntiude of large positive mass balance
c       residuals by making one or more basis switches.
c
c       Note that this block utilizes the efac array calculated above.
c       This array is primarily thought of as being connected with
c       the continued fraction (cycle) algorithm, but it has its use
c       in automatic basis switching mode.
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
        if (nswtch .gt. 0) then
c
          do ns = 1,nstmax
            mosp(ns) = 0.
          enddo
c
          av = -99999.
          call initav(losp,nstmax,av)
c
c         Reset the ixbasp and cjbasp arrays. The former is a flag
c         array, each member of which denotes whether the
c         thermodynamic activity of the corresponding basis species
c         is defined in terms of molality (= 0) or mole fraction (= 1).
c         The cjbasp array contains any site stoichiometric factors
c         associated with the operational basis species.
c
          call gibasp(cgexj,cjbasp,iern1,ixbasp,jern1,jern2,
     $    jetmax,jgext,narn1,narn2,nbasp,nbt,nbtmax,nern1,nern2,
     $    netmax,nphasx,nstmax)
c
c         Go back for another loop.
c
          go to 200
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Sort the zeta-squared terms.
c
c     Calling sequence substitutions:
c       zesort for asort
c       zetasq for aval
c       jasort for jsort
c       iastak for istack
c       jastak for jstack
c       kmax for nmax
c       kdim for nval
c
c     Insure that each sort made by calling EQLIBU/qsortw.f is made
c     starting from scratch. That is, the jasort array left over
c     from a previous call to EQLIBU/qsortw is not used as a starting
c     point for the current sort.
c
      do k = 1,kdim
        jasort(k) = k
      enddo
c
      call qsortw(zesort,zetasq,iastak,jasort,jastak,kmax,noutpt,
     $ nttyo,kdim)
c
c     Reverse the sorting order, so that the values proceed from
c     greatest to smallest.
c
      do k = 1,kdim
        jastak(k) = jasort(kdim + 1 - k)
      enddo
c
      do k = 1,kdim
        krow = jastak(k)
        jasort(k) = krow
        zesort(k) = zetasq(krow)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Find the set of zeta-squared terms to work on. Consider only
c     the kkdim largest ones, those that account for 99% of aleph.
c     If there are more than three, truncate to three.
c
      ax = 0.
      kkdim = 0
      zetsum = 0.
      do k = 1,kdim
        kkdim = kkdim + 1
        zetsum = zetsum + zesort(k)
        ax = zetsum/aleph
        if (ax .ge. 0.99) go to 270
      enddo
  270 continue
      if (kkdim .gt. 3) kkdim = 3
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the vector of partial derivatives d aleph/d log z(j),
c     where aleph is the sum of squares of the zeta residuals and
c     log z(j) is j-th primary iteration variable. At present, z(j)
c     is always the number of moles of the j-th member of the
c     "extended" basis set. Note that:
c
c       d aleph/d log z(j) = 2 * Sum(i) zeta(i)*(d zeta(i)/d log z(j))
c
c     where:
c
c       d zeta(i)/d log z(j) = d alpha(i)/d log z(j)
c
c            for zeta(i) = alpha(i);
c
c
c       d zeta(i)/d log z(j) = d beta(i)/d log z(j)
c
c         = [1./amtb(i)] * d alpha(i)/d log z(j)
c
c            for zeta(i) = beta(i)
c
c       where amtb(i) is treated as a constant, which it really isn't
c       when i refers to a species like H2O, H+, OH-, O2(g,aq), or e-;
c       and
c
c
c       d zeta(i)/d log z(j) = [0.001/cscale(i)] * d alpha(i)/d log z(j)
c
c            for zeta(i) = [0.001/cscale(i)] * alpha(i)
c
c
c     Note that  d alpha(i)/d log z(j) = Jacob(i,j), where [Jacob] is
c     the Jacobian matrix used in Newton-Raphson iteration.
c
c     Get the Jacobian.
c
      call matrix(aamatr,al10,bpx,cdrs,cdrtw,cdrw,cjbasp,
     $ cnufac,conc,csts,dlogxw,eps100,ibpxmx,iebal,iern1,ietmax,
     $ iindx1,ipndx1,irdxc3,ixbasp,ixrn1,jcsort,jern1,jern2,jetmax,
     $ jflag,jjsort,jsitex,kbt,kction,kdim,kelect,khydr,kmax,kmt,
     $ km1,ko2gaq,kwater,kxt,kx1,mosp,narn1,narn2,nbasp,nbt,nbtmax,
     $ nbw,ncosp,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,netmax,
     $ noutpt,no2gaq,nphasx,nredox,nst,nstmax,nsts,nstsmx,nstsr,
     $ ntfx,ntfxmx,ntfxt,nttyo,nxtmax,omega,qredox,q6mode,tfx,
     $ ugexmo,uspec,weight,xbar,xbarw,xbarwc,zchar)
c
c     Compute the d aleph/d log z(j) vector.
c
      do kcol = 1,kdim
        daleph(kcol) = 0.
      enddo
c
      do krow = 1,kbt
        nb = iindx1(krow)
        bx = beta(krow)
        if (bx .ge. 0.0) then
c
c         Have a non-negative beta(i) value: zeta(i) = beta(i).
c
          if (amtb(nb) .gt. eps100) then
            cx = 2.*zeta(krow)/amtb(nb)
          else
            cx = 2.*zeta(krow)*eps1hi
          endif
        else
c
c         Have a negative beta(i) value.
c
          dx = 1.0 + bx
          if (dx .gt. eps100) then
c
c           Have (-1.0 + eps100) < beta(i) < 0.0:
c           zeta(i) = 1.0 - 1.0/(1.0 + beta(i)).
c
            cx = 2.*zeta(krow)**3/amtb(nb)
          else
c
c           Have beta(i) < (-1.0 + eps100):
c           zeta(i) = 1.0 - (1.0/eps100).
c
            if (amtb(nb) .gt. eps100) then
              cx = 2.*eps1hi**3/amtb(nb)
            else
              cx = 2.*eps1hi**4
            endif
          endif
        endif
        do kcol = 1,kdim
          daleph(kcol) = daleph(kcol) + cx*aamatr(krow,kcol)
        enddo
      enddo
c
      do krow = km1,kxt
        ns = iindx1(krow)
        cx = 0.001/cscale(ns)
        do kcol = 1,kdim
          daleph(kcol) = daleph(kcol) + cx*aamatr(krow,kcol)
        enddo
      enddo
c
      if (iodb(3) .ge. 4) then
        write (noutpt,1330)
 1330   format(/7x,'--- Partial Derivatives of Aleph ---',
     $  //'  krow    Variable',25x,'daleph(krow)',/)
        do krow = 1,kdim
c
c         Calling sequence substitutions:
c           uzvec1(krow) for unam48
c
          call fmspnx(jlen,uzvec1(krow),uspn56)
c
          write (noutpt,1340) krow,uspn56,daleph(krow)
 1340     format(2x,i3,3x,a32,3x,1pe12.5)
        enddo
        write (noutpt,1030)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Now find the primary iteration variable that is the principal
c     determining variable for each of the most significant zeta-
c     squared variables.
c
      do k = 1,kkdim
        kmvar(k) = 0
      enddo
c
      do k = 1,kkdim
        krow = jasort(k)
        if (krow .le. kbt) then
          nb = iindx1(krow)
          bx = beta(krow)
          if (bx .ge. 0.0) then
c
c           Have a non-negative beta(i) value: zeta(i) = beta(i).
c
            if (amtb(nb) .gt. eps100) then
              cx = 2.*zeta(krow)/amtb(nb)
            else
              cx = 2.*zeta(krow)*eps1hi
            endif
          else
c
c           Have a negative beta(i) value.
c
            dx = 1.0 + bx
            if (dx .gt. eps100) then
c
c             Have (-1.0 + eps100) < beta(i) < 0.0:
c             zeta(i) = 1.0 - 1.0/(1.0 + beta(i)).
c
              if (amtb(nb) .gt. eps100) then
                cx = 2.*zeta(krow)**3/amtb(nb)
              else
                cx = 2.*zeta(krow)**4
              endif
            else
c
c             Have beta(i) < (-1.0 + eps100):
c             zeta(i) = 1.0 - (1.0/eps1000.
c
              if (amtb(nb) .gt. eps100) then
                cx = 2.*eps1hi**3/amtb(nb)
              else
                cx = 2.*eps1hi**4
              endif
            endif
          endif
        else
          ns = iindx1(krow)
          cx = 0.002*zeta(krow)/cscale(ns)
        endif
        kbig = 0
        abig = 0.
        do kcol = 1,kdim
          dij = cx*aamatr(krow,kcol)
          adij = abs(dij)
          if (adij .gt. abig) then
            if (krow.gt.kbt .or. kcol.ne.kwater .or. krow.eq.kwater)
     $        then
c
c             For a mass balance row other than that for H2O, disallow
c             H2O as the principal determining species.
c
              kbig = kcol
              abig = adij
              dijmaj(krow) = dij
            endif
          endif
        enddo
        kmvar(k) = kbig
      enddo
c
c     Special restrictions apply to H2O, H+, O2(g,aq), and e-.
c
      do k = 1,kkdim
        krow = jasort(k)
        kcol = kmvar(k)
c
        if (krow.eq.kwater .or. krow.eq.khydr .or.
     $    krow.eq.ko2gaq .or. krow.eq.kelect) then
c
c         Only H2O may be adjusted for the H2O row, only H+ for the H+
c         row, only O2(g,aq) for the O2 row, and e- for the e- row.
c
          if (kcol .ne. krow) kmvar(k) = krow
          go to 290
        endif
c
        if (krow .le. kbt) then
          kcol = kmvar(k)
          if (kcol.eq.kwater .or. kcol.eq.khydr .or.
     $      kcol.eq.ko2gaq .or. kcol.eq.kelect) then
c
c           H2O may be adjusted only for the H2O row, among mass
c           balance rows, H+ for the H+ row, O2(g,aq) for the O2 row,
c           and e- for the e- row.
c
            if (kcol .ne. krow) kmvar(k) = krow
          endif
        endif
c
  290   continue
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The following is a return point if no cycle corrections were
c     generated and an alternate set of variables to correct is to be
c     tried.
c
  350 continue
c
      if (iodb(3) .ge. 2) then
        write (noutpt,1400)
 1400   format(/11x,'--- Principal Dependencies ---',
     $  //4x,'Row Variable (i)',19x,'Column Variable (j)',/)
c
        do k = 1, kkdim
          krow = jasort(k)
          if (krow .le. kbt) then
            nb = iindx1(krow)
            ns = nbaspd(nb)
          else
            ns = iindx1(krow)
          endif
          kcol = kmvar(k)
c
c         Calling sequence substitutions:
c           jlen1 for jlen
c           uspec(ns) for unam48
c           usp156 for uspn56
c
          call fmspnx(jlen1,uspec(ns),usp156)
c
          if (kcol .gt. 0) then
c
c           Calling sequence substitutions:
c             jlen2 for jlen
c             uzvec1(kcol) for unam48
c             usp256 for uspn56
c
            call fmspnx(jlen2,uzvec1(kcol),usp256)
          else
            usp256 = 'None'
            jlen2 = 4
          endif
c
          write (noutpt,1410) usp156,usp256(1:jlen2)
 1410     format(2x,a32,3x,a)
        enddo
      endif
c
      if (iodb(3) .ge. 2) write (noutpt,1030)
c
      do k = 1,kkdim
        if (kmvar(k) .eq. 0) then
          krow = jasort(k)
c
          if (krow .le. kbt) then
            nb = iindx1(krow)
            ns = nbaspd(nb)
          elseif (krow.ge.km1 .and. krow.le.kxt) then
            ns = iindx1(krow)
          endif
          kcol = krow
          kmvar(k) = kcol
c
c         Calling sequence substitutions:
c           jlen1 for jlen
c           uspec(ns) for unam48
c           usp156 for uspn56
c
          call fmspnx(jlen1,uspec(ns),usp156)
c
c         Calling sequence substitutions:
c           jlen2 for jlen
c           uzvec1(kcol) for unam48
c           usp256 for uspn56
c
          call fmspnx(jlen2,uzvec1(kcol),usp256)
c
          write (noutpt,1440) usp156(1:jlen1),usp256(1:jlen2)
          write (nttyo,1440) usp156(1:jlen1),usp256(1:jlen2)
 1440     format(/" * Warning- (EQ6/optmzr) Couldn't find a variable",
     $    ' to adjust to reduce',/7x,'the aleph residual for ',
     $    a,'. Will use the variable associated',/7x,'with ',a,'.')
        endif
      enddo
c
c     Find and flag all instances in which a principal determining
c     variable duplicates a prior entry in the list. Set knflag(k) = 0
c     for such duplications, otherwise set knflag(k) = 1.
c
      do k = 1,kkdim
        knflag(k) = 1
      enddo
c
      do k = 1,kkdim - 1
        do kk = k + 1,kkdim
          if (kmvar(kk) .eq. kmvar(k)) knflag(kk) = 0
        enddo
      enddo
c
c     Find the number of variables to be corrected in the current
c     cycle.
c
      jdim = 0
      do k = 1,kkdim
        if (knflag(k) .eq. 1) jdim = jdim + 1
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(3) .ge. 3) then
c
c       Write a table containing data on the structure of the
c       corrections to be made for the current step.
c
        write (noutpt,1530)
 1530   format(/6x,'--- Structure for Current Cycle Corrections',
     $  ' ---',/)
c
        write (noutpt,1540)
 1540   format(4x,'Row Variable (i)',19x,'zeta(i)**2',/)
c
        do k = 1,kkdim
          krow = jasort(k)
          if (krow .le. kbt) then
            nb = iindx1(krow)
            ns = nbaspd(nb)
          else
            ns = iindx1(krow)
          endif
c
c         Calling sequence substitutions:
c           jlen1 for jlen
c           uspec(ns) for unam48
c           usp156 for uspn56
c
          call fmspnx(jlen1,uspec(ns),usp156)
c
          write (noutpt,1550) usp156,zesort(k)
 1550     format(2x,a32,3x,1pe12.5)
        enddo
c
        write (noutpt,1560)
 1560   format(//4x,'Column Variable (j)',9x,'d zeta(i)**2/d log z(j)',
     $  /)
c
        do k = 1,kkdim
c
c         Note: write all duplications, so don't test to see if
c         knflag(k) = 1.
c
          krow = jasort(k)
          kcol = kmvar(k)
c
c         Calling sequence substitutions:
c           jlen2 for jlen
c           uzvec1(kcol) for unam48
c           usp256 for uspn56
c
          call fmspnx(jlen,uzvec1(kcol),usp256)
c
          write (noutpt,1550) usp256,dijmaj(krow)
        enddo
        write (noutpt,1030)
      endif
c
      if (jdim .gt. 1) then
c
c       If the magnitude of aleph is quite large, correct only one
c       variable per cycle.
c
        if (aleph .gt. 1.e+10) then
          if (iodb(3) .ge. 3) then
            write (noutpt,1570)
 1570       format(' Because of the large magnitude of aleph, only',
     $      ' one variable will be',/' corrected in this cycle.',/)
          endif
          jdim = 1
          do k = 2,kkdim
            knflag(k) = 0
          enddo
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Execute the cycle algorithm.
c
c     Testing has shown that the alpha residual for the O2(g,aq) mass
c     balance taken as a function of the log fO2 variable can have
c     a shape like the following:
c
c
c                   -----------------------------------
c                  |    x                              |
c                  |    x                              |
c                  |    x                              |
c                  |    x                              |
c                  |    x                              |
c            ^     |    x                              |
c            |     |     x                             |
c                  |      x                            |
c          alpha   |       xxxxxxxxxxxxxx              |
c                0 |-----------------------x-----------|
c                  |                          x        |
c                  |                           x       |
c                  |                            x      |
c                  |                            x      |
c                  |                            x      |
c                  |                            x      |
c                  |                            x      |
c                  |                            x      |
c                  |                            x      |
c                   -----------------------------------
c
c                               log fO2 ->
c
c     For this variable, the zeta residual is the alpha residual.
c     The function is nearly flat over about an 8-10 unit range of
c     log fO2, and very steep on the sides. It looks like a step
c     function rotated 90 degrees. The corresponding plot of aleph
c     versus log fO2 looks like a square well.
c
c     On the limbs, the slopes are so steep that attempting to
c     extrapolate to a zero residual using first-order partial
c     derivatives produces corrections that are practically
c     infinitesimal. Thus, any method that generates the values of
c     corrections from first-order partial derivatives is pretty much
c     guaranteed to fail. This rules out using whole classes of
c     methods.
c
c     Within the "flat" region, the Newton-Raphson method is
c     pretty much guaranteed to succeed. In the present subroutine
c     what we need is something that will get us in or near the flat
c     region when we are literally "out on a limb."
c
c     One possibility is to use a scanning procedure. Another is to
c     use a non-scanning method that does not use derivatives. An
c     example of the latter would be the cycle algorithm that
c     is employed for pre-Newton-Raphson optimiztion by EQ3NR. However,
c     that algorithm does not extend well to EQ6 problems. In EQ6, mass
c     balances must be considered for basis species such as H2O(l), H+,
c     and O2(g,aq) (or alternatives such as OH- and e-). These mass
c     balances are distinct in that some species can make negative
c     contributions. Also, EQ6 problems require consideration of some
c     heterogeneous mass action equations when the active basis set is
c     "extended" to include non-aqueous species (e.g., minerals).
c
c     In the present version of this subroutine, we will use a simple
c     scannning approach. Although the first-order partial derivatives
c     appear to be not useful in estimating the actual magnitude of
c     needed corrections, they do appear to identify the primary
c     iteration variables that are most in need of correction. They
c     also appear to define the direction in which one should be
c     scanning.
c
      kcorr = 0
      do krow = 1,kdim
        delvec(krow) = 0.
      enddo
c
      j = 0
      do k = 1,kkdim
        if (knflag(k) .eq. 1) then
          j = j + 1
          krscan = jasort(k)
          kcscan = kmvar(k)
c
          if (iodb(3) .ge. 2) then
c
c           Calling sequence substitutions:
c             uzvec1(kcscan) for unam48
c
            call fmspnx(jlen,uzvec1(kcscan),uspn56)
c
            write (noutpt,1610) uspn56(1:jlen)
 1610       format(/'   --- Scanning on ',a,' ---',
     $      //'  zvclg1(kcscan)  zetasq(krscan)',3x,'aleph',10x,
     $     'dscan',/)
            dscan = 0.
            write (noutpt,1620) zvclg1(kcscan),zetasq(krscan),aleph,
     $      dscan
 1620       format(5x,f9.4,3x,1pe12.5,3x,e12.5,3x,e12.5)
          endif
c
c         Get the direction in which to start scanning. Look at the
c         sign of the appropriate partial derivative to get this.
c         The other direction will be scanned if no improvement can
c         be made in the initial direction.
c
c         Note: could try dijmaj(krscan) in place of dalpeh(kcscan)
c         below.
c
          if (daleph(kcscan) .gt. 0.) then
            kdirec = -1
          else
            kdirec = 1
          endif
c
          if (kcscan .eq. krscan) then
            if (zeta(krscan) .gt. 0.) then
              kdirec = -1
            else
              kdirec = 1
            endif
          endif
c
          qscan2 = .false.
c
c         The label below is a return point if the scan is being
c         continued in the opposite direction.
c
  300     continue
c
c         Set the initial value of the scan increment.
c
          if (kcscan .eq. kwater) then
            dscan = 0.0625
            if (abs(zeta(krscan)) .gt. 10.) dscan = 0.125
            if (abs(zeta(krscan)) .gt. 1.e+2) dscan = 0.25
            if (abs(zeta(krscan)) .gt. 1.e+4) dscan = 0.5
            if (abs(zeta(krscan)) .gt. 1.e+8) dscan = 1.0
          elseif (kcscan .eq. khydr) then
            dscan = 0.0625
            if (abs(zeta(krscan)) .gt. 10.) dscan = 0.125
            if (abs(zeta(krscan)) .gt. 1.e+2) dscan = 0.25
            if (abs(zeta(krscan)) .gt. 1.e+4) dscan = 0.5
            if (abs(zeta(krscan)) .gt. 1.e+8) dscan = 1.0
            if (abs(zeta(krscan)) .gt. 1.e+16) dscan = 2.0
          elseif (kcscan .eq. ko2gaq) then
            dscan = 0.25
            if (abs(zeta(krscan)) .gt. 10.) dscan = 0.5
            if (abs(zeta(krscan)) .gt. 1.e+2) dscan = 1.0
            if (abs(zeta(krscan)) .gt. 1.e+4) dscan = 2.0
            if (abs(zeta(krscan)) .gt. 1.e+8) dscan = 4.0
            if (abs(zeta(krscan)) .gt. 1.e+12) dscan = 6.0
            if (abs(zeta(krscan)) .gt. 1.e+16) dscan = 8.0
            if (abs(zeta(krscan)) .gt. 1.e+24) dscan = 12.0
            if (abs(zeta(krscan)) .gt. 1.e+32) dscan = 16.0
          elseif (kcscan .eq. kelect) then
            dscan = 0.125
            if (abs(zeta(krscan)) .gt. 10.) dscan = 0.5
            if (abs(zeta(krscan)) .gt. 1.e+2) dscan = 1.0
            if (abs(zeta(krscan)) .gt. 1.e+4) dscan = 2.0
            if (abs(zeta(krscan)) .gt. 1.e+8) dscan = 4.0
            if (abs(zeta(krscan)) .gt. 1.e+12) dscan = 6.0
            if (abs(zeta(krscan)) .gt. 1.e+16) dscan = 8.0
            if (abs(zeta(krscan)) .gt. 1.e+24) dscan = 12.0
            if (abs(zeta(krscan)) .gt. 1.e+32) dscan = 16.0
          else
            dscan = 0.125
            if (abs(zeta(krscan)) .gt. 10.) dscan = 0.5
            if (abs(zeta(krscan)) .gt. 1.e+2) dscan = 1.0
            if (abs(zeta(krscan)) .gt. 1.e+4) dscan = 2.0
            if (abs(zeta(krscan)) .gt. 1.e+8) dscan = 4.0
            if (abs(zeta(krscan)) .gt. 1.e+12) dscan = 6.0
            if (abs(zeta(krscan)) .gt. 1.e+16) dscan = 8.0
            if (abs(zeta(krscan)) .gt. 1.e+24) dscan = 12.0
            if (abs(zeta(krscan)) .gt. 1.e+32) dscan = 16.0
          endif
c
          icorr = 0
          zvclgi = zvclg1(kcscan)
c
          if (kdirec .gt. 0) then
            zvclim = 99999.
          else
            zvclim = -99999.
          endif
c
          nscan = 0
          qstops = .false.
c
          zvclgo = zvclg1(kcscan)
          zetsqo = zetasq(krscan)
          zetsqm = zetsqo
          alephm = aleph
          nbigger = 0
c
c         The following is a return point for continuing the scan.
c
  410     nscan = nscan + 1
c
c         Increment the variable being scanned.
c
          zvclg1(kcscan) = zvclg1(kcscan) + kdirec*dscan
c
c         Apply a total change limit for a given cycle.
c
          dzxclm = 40.0
          if (kcscan .eq. kwater) dzxclm = 2.0
          if (kcscan .eq. khydr) dzxclm = 4.0
          if (kcscan.ge.km1 .and.kcscan.le.kxt) then
            dzxclm = 0.25
          endif
c
          qsclim = .false.
          dzx = zvclg1(kcscan) - zvclgi
          if (abs(dzx) .gt. dzxclm) then
            qsclim = .true.
            if (kdirec .gt. 0) then
              zvclg1(kcscan) = zvclgi + dzxclm
            else
              zvclg1(kcscan) = zvclgi - dzxclm
            endif
          endif
c
c         Expand the description of the equilibrium system using the
c         incremented variable.
c
          call ncmpex(acflg,act,actlg,cdrs,cegexs,cgexj,conc,
     $    conclg,cpgexs,egexjc,egexjf,egexs,eps100,fo2,fo2lg,fsort,
     $    fugac,fugalg,iern1,iern2,ietmax,ifrn1,ifrn2,igas,igstak,
     $    iindx1,ilrn1,ilrn2,imrn1,imrn2,istack,ixrn1,ixrn2,jcsort,
     $    jern1,jern2,jetmax,jflag,jgext,jgsort,jgstak,jjsort,jpflag,
     $    jsflag,jsitex,jssort,jstack,kbt,kdim,kelect,kmax,km1,ko2gaq,
     $    kwater,kxt,loph,losp,lsort,mgext,mrgexs,mtb,moph,mosp,narn1,
     $    narn2,nbasp,nbt,nbtmax,ncmpr,ndrs,ndrsmx,ndrsr,nelect,nern1,
     $    nern2,netmax,ngexsa,ngext,ngrn1,ngrn2,ngt,ngtmax,noutpt,
     $    no2gaq,nphasx,npt,nptmax,nst,nstmax,nttyo,omega,omeglg,
     $    press,qxbarw,q6mode,ugexj,ugexmo,uphase,uspec,xbar,xbarlg,
     $    xbarw,xbarwc,xbrwlc,xbrwlg,xlks,zchar,zgexj,zvclg1,zvec1)
c
          xbarw = xbar(narn1)
          xbrwlg = xbarlg(narn1)
c
c         Recompute the residuals.
c
          call betas(acflg,actlg,afcnst,alpha,amtb,bbig,beta,
     $    betamx,bneg,cdrs,conc,conclg,coval,csts,eh,ehfac,fo2lg,
     $    ibetmx,iebal,iindx1,irdxc3,jcsort,jflag,jsflag,jssort,kbt,
     $    kdim,kelect,khydr,kmax,km1,ko2gaq,kwater,kxt,mtb,mosp,
     $    narn1,narn2,nbasp,nbtmax,ncosp,ndrs,ndrsmx,ndrsr,nelect,
     $    nern1,nern2,nhydr,noutpt,no2gaq,nredox,nst,nstmax,nsts,
     $    nstsmx,nstsr,ntfx,ntfxmx,ntfxt,nttyo,omega,qredox,q6mode,
     $    tfx,ubbig,ubneg,ubetmx,uspec,uzvec1,weight,xbrwlg,xlke,
     $    xlks,zchar)
c
c         Recompute the zeta residuals and the aleph residual.
c
          do krow = 1,kdim
            bx = beta(krow)
            if (krow .le. kbt) then
c
c             Mass balance element.
c
              if (bx .ge. 0.0) then
c
c               Have a non-negative beta(i) value. Set zeta(i) = beta(i).
c
                zeta(krow) = bx
              else
c
c               Have a negative beta(i) value.
c
                dx = 1.0 + bx
                if (dx .gt. eps100) then
c
c                 Have (-1.0 + eps100) < beta(i) < 0.0.
c                 Set zeta(i) = 1.0 - 1.0/(1.0 + beta(i)).
c
                  zeta(krow) = 1.0 - (1.0/dx)
                else
c
c                 Have beta(i) < (-1.0 + eps100).
c                 Set zeta(i) = 1.0 - (1.0/eps100).
c
                  zeta(krow) = 1.0 - eps1hi
                endif
              endif
            else
c
c             Mass action element.
c
c             zeta(i) = 0.001*beta(i)/cscale(i).
c
              ns = iindx1(krow)
              zeta(krow) = 0.001*bx/cscale(ns)
            endif
          enddo
c
c         Calculate the aleph and zeta-squared functions.
c
          aleph = 0.
          do krow = 1,kdim
            zetasq(krow) = zeta(krow)*zeta(krow)
            aleph = aleph + zetasq(krow)
          enddo
c
          if (iodb(3) .ge. 2) then
            write (noutpt,1620) zvclg1(kcscan),zetasq(krscan),aleph,
     $      dscan
            if (nscan .le. 0) then
c
c             The best point for this scan within the limits of the
c             minimum non-zero scan increment was the initial point.
c
              write (noutpt,1630)
 1630         format(/5x,'No correction was made.')
            endif
          endif
c
c         If the scan has stepped back to a previous point because
c         it gave a smaller residual, exit the scan now that all
c         residual functions have been recalculated.
c
          if (qstops) go to 420
c
c         Now check the results at the current point of the scan.
c         Is the new point better or not? If either the current
c         zeta-squared function or the aleph function has increased,
c         then the scan has over-corrected.
c
          if (aleph .gt. alephm) nbigger = nbigger + 1
          if (zetasq(krscan).gt.zetsqm .or.
     $      (aleph.gt.alephm .and. nbigger.ge.2)) then
c
c           Have over-corrected. Back up.
c
            nscan = nscan - 1
            nbigger = nbigger - 1
            zvclim = zvclg1(kcscan)
            zvclg1(kcscan) = zvclgo
c
c           Note: not everything has been reset here. It is still
c           necessary to go back and re-expand the description of the
c           equilibrium system.
c
            if (dscan .gt. dscmin) then
c
c             Reduce the scan increment.
c
              dscan = 0.5*dscan
c
c             Reduce the scan increment again if the new value of the
c             iteration variable would match the previously tried
c             extreme value.
c
              zvcnew = zvclg1(kcscan) + kdirec*dscan
              if (kdirec .gt. 0) then
                if (zvcnew .gt. (zvclim - eps100)) dscan = 0.5*dscan
              else
                if (zvcnew .lt. (zvclim + eps100)) dscan = 0.5*dscan
              endif
              if (dscan .lt. dscmin) dscan = dscmin
c
c             Go back and try the new, reduced scan increment.
c
            else
c
c             The previous value for this scan was better than any
c             subsequent adjusted value within the limitations imposed
c             by the minimum scan increment.
c
              dscan = 0.
              nscan = nscan  -1
c
c             Set up to stop the current scan.
c
              qstops = .true.
            endif
c
c           Go back and re-expand the system and recalculate all
c           residual functions. If qstops is .true., the scan increment
c           is zero and the only purpose in going back is to ensure
c           that all the residual functions are consistent with the
c           last point of the scan before continuing with the next
c           optimization step (e.g., next scan, cycle, pass), if any.
c
            qsclim = .false.
            go to 410
          endif
c
c         At this point, a non-zero correction has resulted in an
c         improvement, at least in the current zeta-squared function.
c
          icorr = icorr + 1
          zvclgo = zvclg1(kcscan)
          zetsqo = zetasq(krscan)
          zetsqm = min(zetsqm,zetsqo)
          if (aleph.lt. alephm) then
c
c           The aleph function was also improved.
c
            alephm = aleph
            nbigger = 0
          endif
c
          if (.not.qsclim .and. nscan .lt. 10) then
c
            zx = zvclg1(kcscan) - zvclgi
            if (abs(zx) .lt. dzxclm) then
c
c             Continue the current scan, using the same scan increment.
c
              go to 410
            endif
          endif
c
  420     if (iodb(3) .ge. 2) write (noutpt,1030)
c
c         At this point, a scan in one direction is complete.
c
          if (icorr .eq. 0) then
c
c           No correction was made by the scan in the current
c           direction.
c
            if (.not.qscan2) then
c
c             Go back and scan in the other direction.
c
              qscan2 = .true.
              kdirec = -kdirec
              if (iodb(3) .ge. 2) then
                write (noutpt,1640)
 1640           format(5x,'Scanning in the opposite direction',/)
              endif
              go to 300
            endif
          endif
c
c         End corrections for the current variable on this cycle.
c         Calculate the total correction for this scan.
c
          delvec(kcscan) = zvclg1(kcscan) - zvclgi
c
c         Count the number of variables which were actually corrected.
c
          if (icorr .gt. 0) kcorr = kcorr + 1
        endif
      enddo
c
      if (kcorr .le. 0) then
        if (iodb(3) .ge. 2) then
          write (noutpt,1650)
 1650     format(11x,'No corrections were generated.')
        endif
c
c       See if an alternate set of variables to correct should be
c       tried. For example, if H+ was to be corrected to reduce the
c       residual for Al+++, try correcting Al+++.
c
        nchange = 0
        do k = 1,kkdim
          krow = jasort(k)
          kcol = kmvar(k)
          if (kcol .ne. krow) then
            kmvar(k) = krow
            nchange = nchange + 1
          endif
        enddo
        if (nchange .gt. 0) then
          if (iodb(3) .ge. 2) then
            write (noutpt,1660)
 1660       format(11x,'Changing the set of variables to be',
     $      ' corrected',/11x,'and trying again.',/)
          endif
          go to 350
        endif
        write (noutpt,1662)
 1662   format(1x)
      endif
c
c     Check for an oscillating scan.
c
      qscosc = .false.
c
      do kcol = 1,kdim
        if (abs(delvec(kcol) + delprc(kcol)) .gt. eps100) go to 370
      enddo
c
c     The scan is oscillating. See if an alternate set of variables
c     to correct should be tried. For example, if H+ was to be
c     corrected to reduce the residual for Al+++, try correcting Al+++.
c
      if (iodb(3) .ge. 2) then
        write (noutpt,1670)
 1670   format(/11x,'The corrections are oscillating.')
      endif
c
      nchange = 0
      do k = 1,kkdim
        krow = jasort(k)
        kcol = kmvar(k)
        if (kcol .ne. krow) then
          kmvar(k) = krow
          nchange = nchange + 1
        endif
      enddo
      if (nchange .gt. 0) then
        if (iodb(3) .ge. 2) then
          write (noutpt,1660)
        endif
        go to 350
      endif
      write (noutpt,1662)
c
      qscosc = .true.
c
  370 continue
c
      do kcol = 1,kdim
        delprc(kcol) = delvec(kcol)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(3).ge.1 .and. kcorr.gt.0) then
c
c       Summarize the scan corrections for the current cycle.
c
        write (noutpt,1680)
 1680   format(/16x,'--- Cycle Corrections ---',
     $  //2x,'krow   Name',32x,'zvclgi',6x,'zvclg1',6x,'delvec',/)
c
        do k = 1,kkdim
          if (knflag(k) .eq. 1) then
            krow = jasort(k)
            kcol = kmvar(k)
            if (delvec(kcol) .ne. 0.) then
c
c             Calling sequence substitutions:
c               uzvec1(kcol) for unam48
c
              call fmspnx(jlen,uzvec1(kcol),uspn56)
              zvclgi = zvclg1(kcol) - delvec(kcol)
              write (noutpt,1690) krow,uspn56,zvclgi,zvclg1(kcol),
     $        delvec(kcol)
 1690         format(1x,i4,2x,a32,2x,f10.4,2x,f10.4,2x,f10.4)
            endif
          endif
        enddo
        write (noutpt,1030)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (kcorr.le.0 .or. qscosc .or. qsclim) go to 400
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
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
  400 if (iodb(3) .ge. 1) then
        write (noutpt,1700) npass,ncycle
 1700   format(13x,'Completed pass ',i3,' in ',i3,' cycles.')
      endif
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
        write (noutpt,1720) bsigmm
 1720   format(/13x,'bsigmm= ',1pe12.5)
        write (noutpt,1730) bfxi
 1730   format(13x,'bfxi= ',1pe12.5)
        write (noutpt,1740) bfje
 1740   format(13x,'bfje= ',1pe12.5)
        j2 = ilnobl(ubgamx(1:24))
        write (noutpt,1750) bgamx,ubgamx(1:j2)
 1750   format(13x,'bgamx= ',1pe12.5,', ubgamx= ',a)
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
          if (iodb(3) .ge. 1) then
            if (ker .lt. 2) then
              write (noutpt,1800)
 1800         format(/'   Done. Optimization ended outside requested',
     $        ' limits.',/)
            else
              write (noutpt,1810)
 1810         format(/'   Done. Optimization ended outside allowable',
     $        ' limits.',/)
            endif
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
        if (iodb(3) .ge. 1) then
          write (noutpt,1830)
 1830     format(/'   Done. Optimization ended within requested',
     $   ' limits.',/)
        endif
        go to 999
c
      elseif (ncycle.eq.0 .and. kcorr .le. 0) then
c
c       No cycle corrections were generated after starting a new
c       pass. Quit.
c
        if (iodb(3) .ge. 1) then
          if (ker .lt. 2) then
            write (noutpt,1800)
          else
            write (noutpt,1810)
          endif
        endif
        go to 999
c
      elseif (npass .le. 2) then
c
c       The pass convergence criteria are satisfied, but the cycle
c       convergence criteria are not.
c
c       Try another pass.
c
        go to 210
c
      else
c
c       Quit. Optimization ended outside requested limits
c       because cycle requirements were not met.
c
        if (iodb(3) .ge. 1) then
          if (ker .lt. 2) then
            write (noutpt,1800)
          else
            write (noutpt,1810)
          endif
        endif
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
