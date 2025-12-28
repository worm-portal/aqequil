      subroutine eqcalc(aamatr,abar,acflg,acflgo,act,actlg,adh,
     $ adhh,adhv,afcnst,alpha,al10,amtb,aphi,apx,avcnst,azero,
     $ a3bar,a3bars,bacfmx,bbig,bdh,bdhh,bdhv,bdot,bdoth,bdotv,
     $ beta,betamx,betao,bfje,bfxi,bgamx,bneg,bpx,bsigmm,cco2,
     $ cegexs,cess,cdrs,cdrsd,cdrsx,cdrtw,cdrw,cjbasp,cnufac,conc,
     $ conclg,cpgexs,cscale,csts,delmax,delvco,delvec,dlogxw,
     $ egexjc,egexjf,egexs,eh,ehfac,eps100,farad,fje,fjeo,fo2,
     $ fo2lg,fsort,fugac,fugalg,fxi,fxio,gmmatr,iapxt,ibetmx,ibpxt,
     $ ibswx,idelmx,ielam,ier,iern1,iern2,ifcphi1,ifcphi2,ifnnn,
     $ ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igstak,iindx1,
     $ ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,
     $ ilzeta,imrn1,imrn2,insgf,iodb,iopg,iopt,ipch,ipcv,ipivot,
     $ ipndx1,istack,iter,itermx,ixbasp,ixrn1,ixrn2,izmax,jcsort,
     $ jflag,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jsol,
     $ jssort,jstack,kbt,kction,kdim,kelect,khydr,khydx,km1,kmt,
     $ ko2gaq,krdxsp,kwater,kx1,kxt,loph,losp,lsort,moph,mosp,
     $ mrgexs,mtb,narn1,narn2,narxt,nat,nbasp,nbaspd,nbaspx,nbt,
     $ nbtd,nbw,nchlor,ncmpr,nct,ndrs,ndrsd,ndrsx,ndrsr,ndrsrd,
     $ ndrsrx,nelect,nern1,nern2,ness,nessr,net,nfac,nfrn1,nfrn2,
     $ ngrn1,ngrn2,ngt,nhydr,nhydx,nlrn1,nlrn2,nlt,nmrn1,nmrn2,
     $ nmt,noutpt,no2gaq,nphasx,npt,nrdxsp,nst,nsts,nstsr,ntpr,
     $ nttyo,nxrn1,nxrn2,nxt,omega,omeglg,press,qbassw,qhawep,
     $ qoptmz,qpit75,qredox,q6mode,rhsvec,screwd,sigmam,sigmmo,
     $ smp100,tempc,tempk,tolbt,toldl,ubacmx,ubgamx,ulbeta,uldel,
     $ uphase,uspec,uzvec1,weight,wfac,xbar,xbarlg,xbarw,xbarwc,
     $ xbrwlc,xbrwlg,zchar,zchcu6,zchsq2,zvclg1,zvec1)
c
c     This subroutine attempts to calculate the chemical equilibrium
c     state of a system for a pre-defined phase assemblage. Unlike in
c     previous versions of EQ6, this subroutine no longer has the
c     function of determining the actual or most stable phase
c     assemblage. Adjustments to the phase assemblage are now carried
c     out by EQ6/eqphas.f, which as necessary makes repeated calls to
c     the present subroutine. The present subroutine drives
c     EQ6/optmzr.f and EQLIB/newton.f, which respectively carry out
c     pre-Newton-Raphson optimization of starting estimates and hybrid
c     Newton-Raphson iteration.
c
c     Note that there are no attempts in this subroutine to diagnose
c     the cause of any failure to converge. Such a failure might occur,
c     for example, because a phase needs to be deleted from the phase
c     assemblage. It is intended that this subroutine normally be
c     called by EQ6/eqphas.f, which does perform such diagnostics.
c
c     This subroutine is called by:
c
c       EQ6/eqphas.f
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
c       iter   = the number of hydrbid Newton-Raphson iterations
c                  done by EQLIB/newton.f
c
c       ier    = error flag (returned from EQLIB/newton.f):
c                  =  0  Okay
c                  =  1  Encountered a zero matrix
c                  =  2  Encountered a non-zero, computationally
c                          singular matrix
c                  =  3  Iteration was diverging
c                  =  4  Hit the maximum number of iterations (itermx)
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
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer noutpt,nttyo
c
      integer nfac(nbtmax)
c
      integer iapxt(nxtmax),ibpxt(nxtmax),ibswx(nbtmax),igstak(ngtmax),
     $ iindx1(kmax),ipndx1(kmax),insgf(natmax),iodb(nodbmx),
     $ iopg(nopgmx),iopt(noptmx),ipivot(kmax),istack(nstmax),
     $ ixbasp(nbtmax),jcsort(nstmax),jflag(nstmax),jgsort(ngtmax),
     $ jgstak(ngtmax),jjsort(nstmax),jpflag(nptmax),jsflag(nstmax),
     $ jsitex(nstmax),jsol(nxtmax),jssort(nstmax),jstack(nstmax),
     $ kction(nbtmax)
c
      integer narxt(ntprmx),nbasp(nbtmax),nbaspd(nbtmax),nbaspx(nbtmax),
     $ ncmpr(2,nptmax),ndrs(ndrsmx),ndrsd(ndrsmx),ndrsx(ndrsmx),
     $ ndrsr(2,nstmax),ndrsrd(2,nstmax),ndrsrx(2,nstmax),ness(nessmx),
     $ nessr(2,nstmax),nphasx(nstmax),nsts(nstsmx),nstsr(2,nstmax)
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
      integer ibetmx,idelmx,ielam,ier,igas,ipch,ipcv,iter,itermx,izmax,
     $ kbt,kdim,kelect,khydr,khydx,km1,kmt,ko2gaq,krdxsp,kwater,kx1,
     $ kxt,nbtd,nbw,nchlor,nelect,nhydr,nhydx,no2gaq,nrdxsp,ntpr
c
      logical qbassw,qhawep,qoptmz,qpit75,qredox,q6mode
c
      character*48 uspec(nstmax),uzvec1(kmax)
      character*48 ubacmx,ubgamx
      character*24 uphase(nptmax)
      character*8 ulbeta(kmax),uldel(kmax)
c
      real*8 aamatr(kmax,kmax),acflg(nstmax),acflgo(nstmax),
     $ act(nstmax),actlg(nstmax),alpha(kmax),amtb(nbtmax),
     $ apx(iapxmx,nxtmax),azero(natmax),a3bars(natmax),
     $ beta(kmax),betao(kmax),bpx(ibpxmx,nxtmax),cco2(5),
     $ cegexs(ietmax,jetmax,netmax),cess(nessmx),cdrs(ndrsmx),
     $ cdrsd(ndrsmx),cdrsx(ndrsmx),cdrtw(nstmax),cdrw(nstmax),
     $ cjbasp(nbtmax),cnufac(nstmax),conc(nstmax),conclg(nstmax),
     $ cpgexs(ietmax,jetmax,netmax),cscale(nstmax),
     $ csts(nstsmx),delvco(kmax),delvec(kmax),dlogxw(nbtmax),
     $ egexjc(jetmax,netmax),egexjf(jetmax,netmax),
     $ egexs(ietmax,jetmax,netmax)
c
      real*8 fsort(ngtmax),fugac(ngtmax),fugalg(ngtmax),
     $ gmmatr(kmax,kmax),loph(nptmax),losp(nstmax),lsort(nstmax),
     $ moph(nptmax),mosp(nstmax),mrgexs(ietmax,jetmax,netmax),
     $ mtb(nbtmax),rhsvec(kmax),weight(nstmax),wfac(iktmax,nxtmax),
     $ xbar(nstmax),xbarlg(nstmax),zchar(nstmax),zchcu6(nstmax),
     $ zchsq2(nstmax),zvclg1(kmax),zvec1(kmax)
c
      real*8 adh,adhh,adhv,aphi,bdh,bdhh,bdhv,bdot,bdoth,bdotv
c
      real*8 abar,afcnst,al10,avcnst,a3bar,bacfmx,bbig,betamx,bfje,
     $ bfxi,bgamx,bneg,bsigmm,delmax,eh,ehfac,eps100,farad,fje,fjeo,
     $ fo2,fo2lg,fxi,fxio,omega,omeglg,press,screwd,sigmam,sigmmo,
     $ smp100,tempc,tempk,tolbt,toldl,xbarw,xbarwc,xbrwlc,xbrwlg
c
c-----------------------------------------------------------------------
c
c     Local variable declarations with global dimensioning.
c
c     Note: the following variables are not used by EQ6 except to
c     satisfy calls to the following subroutines shared with EQ3NR:
c
c       EQLIB/betas.f
c       EQLIB/newton.f
c
c     These subroutines are called directly in this subroutine, and
c     indirectly through EQ6/optmzr.f. Some of these variables are
c     not arrays.
c
      integer isv_nbtmax,isv_ntfxmx
c
      SAVE isv_nbtmax,isv_ntfxmx
c
      integer, dimension(:), allocatable :: ncosp
      integer, dimension(:), allocatable :: ntfx
c
      SAVE ncosp,ntfx
c
      real(8), dimension(:), allocatable :: coval
      real(8), dimension(:), allocatable :: tfx
c
      SAVE coval,tfx
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,j3,kcol,n
c
      integer iebal,irdxc3,nredox,ntfxt
c
      integer ilnobl
c
      logical qblamx,qxbarw
c
      character*48 ubbig,ubetmx,ubneg
      character*8 uaff,umoles,ux8,ux8a,ux8b
c
      real*8 actwlc,screwn,tolxpt,tolbig,tolneg
c
c-----------------------------------------------------------------------
c
      data iebal  /0/
      data irdxc3 /0/
      data nredox /0/
      data ntfxt  /0/
c
      data qxbarw/.false./
c
c-----------------------------------------------------------------------
c
c     Set the screwn parameter.
c
      data screwn /0.10/
c
c-----------------------------------------------------------------------
c
c     The following are tolerance parameters for the optimization
c     carried out by EQ6/optmzr.f. These are convergence tolerances
c     within that subroutine; however, they are used in the present
c     subroutine to determine if pre-Newton-Raphson iteration is
c     necessary.
c
      data tolbig /10.0/,tolneg /-0.9/
      data tolxpt /10.0/
c
c-----------------------------------------------------------------------
c
      data uaff /'aff     '/,umoles /'moles   '/
c
c-----------------------------------------------------------------------
c
c     Allocate or reallocate local work arrays as needed.
c
      if (.not.ALLOCATED(ncosp)) then
c
c       Local work arrays are not allocated. Zero the saved
c       array size variables. Note that only one array is tested
c       to see if it is allocated. It is assumed that all local
c       work arrays are either allocated or not.
c
        isv_nbtmax = 0
        isv_ntfxmx = 0
      else
c
c       Local work arrays are allocated. Check to see if any of the
c       array size variables have changed. If so, deallocate
c       the corresponding local work arrays and zero the corresponding
c       saved size variables.
c
        if (nbtmax .ne. isv_nbtmax) then
          DEALLOCATE(ncosp)
          DEALLOCATE(coval)
          isv_nbtmax = 0
        endif
c
        if (ntfxmx .ne. isv_ntfxmx) then
          DEALLOCATE(ntfx)
          DEALLOCATE(tfx)
          isv_ntfxmx = 0
        endif
      endif
c
c     At this point, the saved array size values are zero if the
c     corresponding arrays need to be allocated.
c
      if (isv_nbtmax .eq. 0) then
        ALLOCATE(ncosp(nbtmax))
        ALLOCATE(coval(nbtmax))
        isv_nbtmax = nbtmax
      endif
c
      if (isv_ntfxmx .eq. 0) then
        ALLOCATE(ntfx(ntfxmx))
        ALLOCATE(tfx(ntfxmx))
        isv_ntfxmx = ntfxmx
      endif
c
c     Zero the contents of the local work arrays.
c
      do n = 1,nbtmax
        ncosp(n) = 0
      enddo
c
      do n = 1,nbtmax
        dlogxw(n) = 0.
        coval(n) = 0.
      enddo
c
      do n = 1,ntfxmx
        ntfx(n) = 0
      enddo
c
      do n = 1,ntfxmx
        tfx(n) = 0.
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set up the ulbeta and uldel arrays.
c
      call initcv(ulbeta,kbt,umoles)
      call initcv(uldel,kbt,umoles)
c
      do kcol = km1,kxt
        ulbeta(kcol) = uaff
        uldel(kcol) = umoles
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set the qblamx flag, which indicates the presence of solid
c     solutions in the currently specified phase assemblage.
c
      qblamx = kxt .ge. kx1
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Expand the system description.
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
c     Compute the residual functions prior to iteration.
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
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qoptmz) then
        if (betamx.gt.tolxpt .or. bbig.gt.tolbig .or. bneg.lt.tolneg)
     $    then
cXXX
c         write (noutpt,3000) betamx,tolxpt
c         write (nttyo,3000) betamx,tolxpt
c3000     format(//' ***** betamx= ',1pe11.4,' and tolxpt= ',1pe11.4,
c    $    ' ****')
c         j2 = ilnobl(ubetmx)
c         write (noutpt,3005) ubetmx(1:j2)
c         write (nttyo,3005) ubetmx(1:j2)
c3005     format('   ***** ubetmx= ',a,' *****',/)
c         write (noutpt,3010) bbig,tolbig
c         write (nttyo,3010) bbig,tolbig
c3010     format(/' ***** bbig= ',1pe11.4,' and tolbig= ',1pe11.4,
c    $    ' ****')
c         j2 = ilnobl(ubbig)
c         write (noutpt,3015) ubbig(1:j2)
c         write (nttyo,3015) ubbig(1:j2)
c3015     format('   ***** ubbig= ',a,' *****',/)
c         write (noutpt,3020) bneg,tolneg
c         write (nttyo,3020) bneg,tolneg
c3020     format(/' ***** bneg= ',1pe11.4,' and tolneg= ',1pe11.4,
c    $    ' ****')
c         j2 = ilnobl(ubneg)
c         write (noutpt,3025) ubneg(1:j2)
c         write (nttyo,3025) ubneg(1:j2)
c3025     format('   ***** ubbneg= ',a,' *****',/)
cXXX
c
c         Optimize the iteration variables before starting hybrid
c         Newton-Raphson iteration. The phase assemblage is fixed
c         in this process.
c
          call optmzr(aamatr,abar,acflg,acflgo,act,actlg,
     $    adh,adhh,adhv,afcnst,al10,alpha,amtb,aphi,avcnst,azero,
     $    a3bar,a3bars,bacfmx,bbig,bdh,bdhh,bdhv,bdot,bdoth,bdotv,
     $    beta,betamx,bgamx,bneg,bpx,cco2,cdrs,cdrsx,cdrtw,cdrw,
     $    cegexs,cgexj,cjbasp,cnufac,conc,conclg,cpgexs,cscale,
     $    csts,coval,delvec,dlogxw,egexjc,egexjf,egexs,ehfac,
     $    eps100,fje,fjeo,fo2,fo2lg,fsort,fugac,fugalg,fxi,fxio,
     $    gmmatr,ibetmx,ibpxt,ibswx,ielam,iern1,iern2,ifcphi1,
     $    ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,
     $    igas,igstak,iindx1,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,
     $    ilpsi2,ilrn1,ilrn2,ilzeta,imrn1,imrn2,insgf,iodb,iopg,
     $    iopt,ipch,ipcv,ipivot,ipndx1,irdxc3,istack,ixbasp,ixrn1,
     $    ixrn2,izmax,jcsort,jern1,jern2,jflag,jgext,jgsort,jgstak,
     $    jjsort,jpflag,jsflag,jsitex,jsol,jssort,jstack,kbt,kction,
     $    kdim,kelect,khydr,km1,kmt,ko2gaq,kwater,kx1,kxt,loph,losp,
     $    lsort,mgext,moph,mosp,mrgexs,mtb,narn1,narn2,narxt,nat,
     $    nbasp,nbaspd,nbaspx,nbw,nbt,nbtd,nchlor,ncmpr,ncosp,ndrs,
     $    ndrsx,ndrsr,ndrsrd,ndrsrx,nelect,nern1,nern2,net,ngexsa,
     $    nfac,ngext,nhydr,nhydx,ngrn1,ngrn2,ngt,noutpt,no2gaq,
     $    nphasx,npt,nst,nsts,nstsr,ntfx,ntfxt,ntpr,nttyo,omega,
     $    omeglg,press,qbassw,qblamx,qhawep,qpit75,qredox,q6mode,
     $    rhsvec,sigmam,sigmmo,smp100,tempc,tempk,tfx,tolbig,tolneg,
     $    tolxpt,ubacmx,ubbig,ubetmx,ubgamx,ubneg,ugexj,ugexmo,
     $    uphase,uspec,uzvec1,weight,wfac,xbar,xbarlg,xbarw,xbarwc,
     $    xbrwlc,xbrwlg,zchar,zchcu6,zchsq2,zgexj,zvclg1,zvec1)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Do hybrid Newton-Raphson iteration for a given phase assemblage.
c
      if (iodb(4) .ge. 1) then
        write (noutpt,1600)
 1600   format(/,' Starting hybrid Newton-Raphson iteration.',/)
      endif
c
      call newton(aamatr,abar,acflg,acflgo,act,actlg,actwlc,
     $ adh,adhh,adhv,afcnst,alpha,al10,amtb,aphi,azero,a3bar,a3bars,
     $ bacfmx,bbig,beta,betamx,betao,bdh,bdhh,bdhv,bdot,bdoth,bdotv,
     $ bfje,bfxi,bgamx,bneg,bpx,bsigmm,cco2,cdrs,cdrtw,cdrw,cegexs,
     $ cgexj,cjbasp,cnufac,conc,conclg,coval,cpgexs,csts,delam,
     $ delmax,delvco,delvec,dgpit,dlogxw,dpelm,dpslm,dselm,egexjc,
     $ egexjf,egexs,eh,ehfac,elam,eps100,fje,fjeo,fo2,fo2lg,fsort,
     $ fugac,fugalg,fxi,fxio,gmmatr,gpit,ibpxt,idelmx,iebal,ielam,
     $ ier,iern1,iern2,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,
     $ ifrn1,ifrn2,ifzeta,igas,igstak,iindx1,ilcphi1,ilcphi2,ilnnn,
     $ iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,ilzeta,imrn1,imrn2,insgf,
     $ iodb,iopg,ipivot,ipndx1,irdxc3,istack,iter,itermx,ixbasp,
     $ ixrn1,ixrn2,izmax,jcsort,jern1,jern2,jflag,jgext,jgsort,
     $ jgstak,jjsort,jpflag,jsflag,jsitex,jsol,jssort,jstack,kbt,
     $ kction,kdim,kelect,khydr,km1,kmt,ko2gaq,kwater,kx1,kxt,loph,
     $ losp,lsort,mgext,moph,mosp,mrgexs,mtb,nalpha,napt,narn1,
     $ narn2,nbasp,nbt,nbw,nchlor,ncmpr,ncosp,ndrs,ndrsr,nelect,
     $ nern1,nern2,net,ngexsa,ngext,ngrn1,ngrn2,ngt,nhydr,nmut,nmux,
     $ nmxi,nmxx,noutpt,no2gaq,nphasx,npt,nredox,nslt,nslx,nst,nsts,
     $ nstsr,nsxi,nsxx,ntfx,ntfxt,nttyo,omega,omeglg,palpha,pelm,
     $ pmu,press,pslamn,pslm,qhawep,qpit75,qredox,q6mode,rhsvec,
     $ screwd,screwn,selm,sigmam,sigmmo,tempk,tfx,tolbt,toldl,ubacmx,
     $ ubbig,ubetmx,ubgamx,ubneg,ugexj,ugexmo,ulbeta,uldel,uphase,
     $ uspec,uzvec1,weight,wfac,xbar,xbarlg,xbarw,xbarwc,xbrwlc,
     $ xbrwlg,xlke,xlks,zchar,zchsq2,zchcu6,zgexj,zvclg1,zvec1)
c
      if (ier .gt. 0) then
        if (iodb(4) .gt. 0) then
          write (noutpt,1610)
 1610     format(/' * Note - (EQ6/eqcalc) Hybrid Newton-Raphson',
     $    ' iteration failed.')
          ux8a = ' '
          write (ux8a,'(i7)') iter
          call lejust(ux8a)
          j2 = ilnobl(ux8a)
          if (ier. eq. 1) then
            write (noutpt,1630) ux8a(1:j2)
 1630       format(7x,'after ',a,' iterations because a zero matrix',
     $      ' was encountered. This is',/7x,'probably due to a',
     $      ' programming error.')
          elseif (ier .eq. 2) then
            write (noutpt,1640) ux8a(1:j2)
 1640       format(7x,'after ',a,' iterations because a non-zero,',
     $      ' computationally singular',/7x,'matrix was encountered.')
          elseif (ier .eq. 3) then
            write (noutpt,1650) ux8a(1:j2)
 1650       format(7x,'after ',a,' iterations because the code',
     $      ' detected that',/7x,'iteration was diverging.')
          elseif (ier .eq. 4) then
            write (noutpt,1660) ux8a(1:j2)
 1660       format(7x,'after ',a,' iterations because the maximum',
     $      ' number of iterations',/7x,'was done.')
          else
            ux8b = ' '
            write (ux8b,'(i7)') ier
            call lejust(ux8b)
            j3 = ilnobl(ux8b)
            write (noutpt,1670) ux8a(1:j2),ux8b(1:j3)
 1670       format(7x,'after ',a,' iterations because an unknown event',
     $      ' occurred. The ier',/7x,'error code has the unknown value',
     $      ' ',a,'. This condition is a',/7x,'programming error.')
          endif
        endif
        go to 999
      endif
c
      if (iodb(4) .ge. 1) then
        ux8 = ' '
        write (ux8,'(i7)') iter
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1690) ux8(1:j2)
 1690   format('   Done. Hybrid Newton-Raphson iteration converged in ',
     $  a,' iterations.',/)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
