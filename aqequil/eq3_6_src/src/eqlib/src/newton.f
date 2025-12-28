      subroutine newton(aamatr,abar,acflg,acflgo,act,actlg,actwlc,
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
c     This subroutine performs hybrid Newton-Raphson iteration to solve
c     for the equilibrium state of a chemical system.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ6/eqcalc.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ielam  = flag to not use (-1) or use (0) "higher order"
c                  (higher than 2nd order) electrostatic terms in
c                  those activity coefficient models which contain
c                  provision for such terms
c       screwd = under-relaxation control parameter. It is used to
c                  reduce the magnitude of the del vector, if necessary
c                  so that the magnitude of the largest element of that
c                  vector does not exceed screwd.
c       screwn = under-relaxation control parameter
c       tolbt  = convergence tolerance on betamx
c       toldl  = convergence tolerance on delmax
c       ulbeta = label for the quantity type for beta vector;
c                  e.g., 'conc' or 'moles'
c       uldel  = label for the quantity type for del corrections;
c                 e.g., 'conc' or 'moles'
c       uspec  = aqueous species name array
c       uzvec1 = name array corresponding to zvclg1
c       iodb   = array of debugging print options
c       itermx = maximum number of iterations
c       kdim   = dimension of aamatr
c       narn1  = start of species range for aqueous solution
c       narn2  = end of species range for aqueous solution
c       nbasp  = array of basis species
c       nbt    = number of basis species
c       ncmpr  = species range array for phases
c       nern1  = start of species range for generic ion exchangers
c       nern2  = end of species range for generic ion exchangers
c       nst    = number of aqeuous species
c       q6mode = flag denoting usage for EQ3NR or EQ6:
c                  .false. = EQ3NR
c                  .true.  = EQ6NR
c       wfac   = array of solid solution non-ideality parameters
c
c     Principal input/output:
c
c       zvclg1 = the 'log z' array, the array corrected by
c                  Newton-Raphson iteration
c       conc   = concentration array
c       acflg  = activity coefficient array
c       fje    = the ionic asymmetry (the 3rd-order electrostatic
c                  moment function J)
c       fxi    = the ionic strength (the 2nd-order electrostatic
c                  moment function I)
c
c     Principal work space/output:
c
c       aamatr = Jacobian matrix
c       abar   = average ion size
c       a3bar  = average cube of the distance of closest approach
c                  for all solute species
c       a3bars = characteristic average cube of the distance of closest
c                  approach for a solute species
c       gmmatr = copy of aamatr
c       delvec = correction array
c       rhsvec = right hand side vector
c       beta   = normalized residual function array
c       alpha  = residual function array
c       acflgo = copy of activity coefficient array
c       betao  = old beta array
c       delvco = old del array
c       betamx = norm (largest magnitude) of the beta array
c       bbig   = largest magnitude positive mass balance residual
c       bneg   = largest magnitude negative mass balance residual
c       bacfmx = norm (largest magnitude) activity coefficient residual
c       ubbig  = name of species corresponding to bbig
c       ubneg  = name of species corresponding to bneg
c       ubacmx = name of species corresponding to bacfmx
c       ipivot = the pivot vector
c       iter   = Newton-Raphson iteration number
c       idelmx = kcol index corresponding to delmax
c       ier    = error flag (mostly returned from EQLIB/nrstep.f):
c                  =  0  Okay
c                  =  1  Encountered a zero matrix
c                  =  2  Encountered a non-zero, computationally
c                          singular matrix
c                  =  3  Iteration was diverging
c                  =  4  Hit the maximum number of iterations (itermx)
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
      integer ibpxt(nxtmax),igstak(ngtmax),iindx1(kmax),insgf(natmax),
     $ iodb(nodbmx),iopg(nopgmx),ipivot(kmax),ipndx1(kmax),
     $ istack(nstmax),ixbasp(nbtmax),jcsort(nstmax),
     $ jern1(jetmax,netmax),jern2(jetmax,netmax),jflag(nstmax),
     $ jgext(netmax),jgsort(ngtmax),jgstak(ngtmax),jjsort(nstmax),
     $ jpflag(nptmax),jsflag(nstmax),jsitex(nstmax),jsol(nxtmax),
     $ jssort(nstmax),jstack(nstmax),kction(nbtmax),nbasp(nbtmax),
     $ ncmpr(2,nptmax),ncosp(nbtmax),ndrs(ndrsmx),ndrsr(2,nstmax),
     $ ngexsa(ietmax,jetmax,netmax),ngext(jetmax,netmax),nphasx(nstmax),
     $ nsts(nstsmx),nstsr(2,nstmax),ntfx(ntfxmx)
c
      integer nalpha(nsltmx),nmux(3,nmutmx),nmxi(2,natmax),
     $ nmxx(3,nmxmax),nslx(2,nsltmx),nsxi(2,natmax),nsxx(2,nsxmax)
c
      integer napt,nmut,nslt
c
      integer ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,
     $ ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta
c
      integer idelmx,iebal,ielam,ier,iern1,iern2,ifrn1,ifrn2,igas,ilrn1,
     $ ilrn2,imrn1,imrn2,iter,itermx,irdxc3,ixrn1,ixrn2,izmax,kbt,kdim,
     $ kelect,khydr,km1,kmt,ko2gaq,kwater,kx1,kxt,narn1,narn2,nbt,nbw,
     $ nchlor,nelect,nern1,nern2,net,ngrn1,ngrn2,ngt,nhydr,no2gaq,npt,
     $ nredox,nst,ntfxt
c
      logical qhawep,qpit75,qredox,q6mode
c
      character*48 uspec(nstmax),uzvec1(kmax)
      character*48 ubbig,ubetmx,ubgamx,ubneg,ubacmx
      character*24 ugexmo(netmax),uphase(nptmax)
      character*8 ugexj(jetmax,netmax)
      character*8 ulbeta(kmax),uldel(kmax)
c
      real*8 dgpit(2,ipbtmx,napmax),dpslm(2,nsltmx),
     $ gpit(ipbtmx,napmax),palpha(ipbtmx,napmax),pmu(nmutmx),
     $ pslamn(0:ipbtmx,nsltmx),pslm(nsltmx)
c
      real*8 delam(2,nazpmx,nazpmx),dpelm(2,nazpmx,nazpmx),
     $ dselm(2,nazmmx:nazpmx),elam(nazpmx,nazpmx),pelm(nazpmx,nazpmx),
     $ selm(nazmmx:nazpmx)
c
      real*8 aamatr(kmax,kmax),acflg(nstmax),acflgo(nstmax),
     $ act(nstmax),actlg(nstmax),alpha(kmax),amtb(nbtmax),azero(natmax),
     $ a3bars(natmax),beta(kmax),betao(kmax),bpx(ibpxmx,nxtmax),
     $ cco2(5),cdrs(ndrsmx),cdrtw(nstmax),cdrw(nstmax),
     $ cegexs(ietmax,jetmax,netmax),cgexj(jetmax,netmax),cjbasp(nbtmax),
     $ cnufac(nstmax),conc(nstmax),conclg(nstmax),coval(nbtmax),
     $ cpgexs(ietmax,jetmax,netmax),csts(nstsmx),delvco(kmax),
     $ delvec(kmax),dlogxw(nbtmax),egexjc(jetmax,netmax),
     $ egexs(ietmax,jetmax,netmax),egexjf(jetmax,netmax)
c
      real*8 fsort(ngtmax),fugac(ngtmax),fugalg(ngtmax),
     $ gmmatr(kmax,kmax),loph(nptmax),losp(nstmax),lsort(nstmax),
     $ mgext(jetmax,netmax),moph(nptmax),mosp(nstmax),
     $ mrgexs(ietmax,jetmax,netmax),mtb(nbtmax),rhsvec(kmax),
     $ tfx(ntfxmx),weight(nstmax),wfac(iktmax,nxtmax),
     $ xbar(nstmax),xbarlg(nstmax),xlks(nstmax),zchar(nstmax),
     $ zchsq2(nstmax),zchcu6(nstmax),zgexj(jetmax,netmax),
     $ zvclg1(kmax),zvec1(kmax)
c
      real*8 adh,adhh,adhv,aphi,bdh,bdhh,bdhv,bdot,bdoth,bdotv
c
      real*8 abar,actwlc,afcnst,al10,a3bar,bacfmx,bbig,betamx,bfje,
     $ bfxi,bgamx,bneg,bsigmm,delmax,eh,ehfac,eps100,fje,fjeo,fo2,
     $ fo2lg,fxi,fxio,omega,omeglg,press,screwd,screwn,sigmam,sigmmo,
     $ tempk,tolbt,toldl,xbarw,xbarwc,xbrwlc,xbrwlg,xlke
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ibetmx,jlen,j2,kcol,negbfc,negdfc,negxfc,noibfc,npconv,
     $ npobfc,npodfc,ns
c
      integer ilnobl
c
      logical qcacf,qcbeta,qcdel,qcgam,qconv,qcsigm,qcfxi,qpracf,qxbarw
c
      character*56 uspn56
      character*8 udelmx,utb,utd
c
      real*8 betfnc,betmxo,btfcnr,bx,dx,chfacf,chfsgm,delfnc,rdx,
     $ rlxfac,rlxgam
c
c-----------------------------------------------------------------------
c
      iter = 0
      ier = 0
      ibetmx = 0
      idelmx = 0
      npconv = 0
c
      delmax = 0.
      betmxo = 0.
      betfnc = 0.
      delfnc = 0.
      btfcnr = 0.
      rlxfac = 1.
      rlxgam = 1.
      negdfc = 0
      negbfc = 0
      negxfc = 0
      noibfc = 0
      npodfc = 0
      npobfc = 0
c
      do kcol = 1,kdim
        delvec(kcol) = 0.
      enddo
c
      qpracf = iodb(4) .ge. 4
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(4) .ge. 2) then
c
c       Debugging prints- state of the system prior to hybrid
c       Newton-Raphson iteration.
c
        write (noutpt,1000)
 1000   format(/8x,'Name',36x,'zvclg1',/)
        do kcol = 1,kdim
c
c         Calling sequence substitutions:
c           uzvec1(kcol) for unam48
c
          call fmspnx(jlen,uzvec1(kcol),uspn56)
          jlen = min(jlen,38)
          write (noutpt,1010) kcol,uspn56(1:jlen),zvclg1(kcol)
 1010     format(1x,i3,2x,a,t46,1pe12.5)
        enddo
        write (noutpt,1015)
 1015   format(1x)
      endif
c
      if (qpracf) then
        write (noutpt,1020) sigmam,fxi,fje,xbrwlc,xbarwc
 1020   format(/3x,'sigmam= ',1pe12.5,/6x,'fxi= ',e12.5,
     $  /6x,'fje= ',e12.5,
     $  //6x,'xbrwlc= ',0pf10.5,/6x,'xbarwc= ',1pe12.5)
c
        write (noutpt,1030)
 1030   format(/9x,'Species',21x,'gamma',/)
        do ns = narn1,narn2
          write (noutpt,1040) ns,uspec(ns),acflg(ns)
 1040     format(1x,i4,2x,a24,2x,1pe12.5)
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the Newton-Raphson residual functions.
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
c     The label below is the return point for iter >= 1.
c
  100 sigmmo = sigmam
      fxio = fxi
      fjeo = fje
c
c     Save the activity coefficients.
c
      call copyaa(acflg,acflgo,nst)
c
      bx = 0.
      dx = 0.
      if (ibetmx .gt. 0) bx = beta(ibetmx)
      if (idelmx .gt. 0) dx = delvec(idelmx)
c
c     Calculate the beta convergence function.
c
      betfnc = 0.
      if (betmxo .gt. 0.) betfnc = (betmxo - betamx)/betmxo
      betmxo = betamx
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(4) .ge. 2) then
        write (noutpt,1100)
 1100   format(/8x,'Name',37x,'Beta',/)
        do kcol = 1,kdim
c
c         Calling sequence substitutions:
c           uzvec1(kcol) for unam48
c
          call fmspnx(jlen,uzvec1(kcol),uspn56)
          jlen = min(jlen,38)
          write (noutpt,1110) kcol,uspn56(1:jlen),beta(kcol)
 1110     format(1x,i3,2x,a,t46,1pe12.5)
        enddo
        write (noutpt,1120)
 1120   format(/1x)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(4) .ge. 1) then
        ubetmx = ' '
        utb = ' '
        udelmx = ' '
        utd = ' '
c
        if (ibetmx .gt. 0) then
          utb = ulbeta(ibetmx)
          ubetmx = uzvec1(ibetmx)
        endif
c
        if (idelmx .gt. 0) then
          utd = uldel(idelmx)
          udelmx = uzvec1(idelmx)
        endif
c
        write (noutpt,1200) iter
 1200   format(' iter= ',i3)
c
        if (abs(rlxgam) .le. eps100) then
          write (noutpt,1210)
 1210     format(7x,'Gammas are fixed')
        elseif (abs(rlxgam - 1.0) .gt. eps100) then
          write (noutpt,1220) rlxgam
 1220     format(7x,'Gamma relaxation factor= ',g12.5)
        endif
c
        write (noutpt,1230) utd,udelmx,dx,delfnc
 1230   format(5x,'delvec(',2a8,')= ',1pe12.5,', delfnc= ',e12.5)
c
        if (abs(rlxfac - 1.0) .gt. eps100) then
          rdx = rlxfac*dx
          write (noutpt,1240) rlxfac,utd,udelmx,rdx
 1240     format(7x,'Relaxation factor= ',g12.5,
     $    /5x,'Relaxed delvec(',2a8,')= ',1pe12.5)
        endif
c
        write (noutpt,1250) utb,ubetmx,bx,betfnc
 1250   format(7x,'beta(',2a8,')= ',1pe12.5,', betfnc= ',e12.5)
c
c       Calling sequence substitutions:
c         ubbig for unam48
c
        call fmspnx(jlen,ubbig,uspn56)
        write (noutpt,1260) bbig,uspn56(1:jlen)
 1260   format(9x,'bbig= ',1pe12.5,',   ubbig= ',a)
c
c       Calling sequence substitutions:
c         ubneg for unam48
c
        call fmspnx(jlen,ubneg,uspn56)
        write (noutpt,1270) bneg,uspn56(1:jlen)
 1270   format(9x,'bneg= ',1pe12.5,',   ubneg= ',a)
c
        j2 = ilnobl(ubgamx(1:24))
        write (noutpt,1280) bgamx,ubgamx(1:j2)
 1280   format(8x,'bgamx= ',1pe12.5,',  ubgamx= ',a)
        write (noutpt,1290) bsigmm,bfxi,bfje
 1290   format(7x,'bsigmm= ',1pe12.5,/9x,'bfxi= ',e12.5,
     $  /9x,'bfje= ',1pe12.5)
        write (noutpt,1300) btfcnr
 1300   format(7x,'btfcnr= ',1pe12.5,/)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check to see if the iteration satisfies the convergence criteria.
c     Both residual functions and correction terms are tested.
c     Iteration may terminate acceptably without satisfying the
c     constraint on the correction terms (see below).
c
      qcbeta = betamx .le. tolbt
      qcdel = delmax .le. toldl
      qconv = qcbeta .and. qcdel
      qcsigm = abs(bsigmm) .le. tolbt
      qcfxi = abs(bfxi) .le. tolbt
      qcgam = bgamx .le. tolbt
      qcacf = qcsigm .and. qcfxi .and. qcgam
c
cXX   Redefine qcacf when the activity coefficients of species in
cXX   non-aqueous phases are updated numerically the same as those
cXX   of aqueous species.
c     qcacf = bacfmx .le. tolbt
c
c     Force to run at least one iteration
c
      if (iter .ge. 1) then
c
c       Test for convergence
c
        if (qconv .and. qcacf) go to 999
      endif
c
      qxbarw = .false.
cXX   qxbarw = abs(bsigmm).le.0.05 .and. bgamx.le.0.05 .and.
cXX  $ abs(rlxgam - 1.0).le.eps100
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Do a Newton-Raphson step.
c
      call nrstep(aamatr,acflg,act,actlg,afcnst,alpha,al10,amtb,
     $ bbig,beta,betamx,betao,betfnc,betmxo,bneg,bpx,btfcnr,cdrs,cdrtw,
     $ cdrw,cegexs,cgexj,cjbasp,cnufac,conc,conclg,cpgexs,csts,coval,
     $ delfnc,delmax,delvco,delvec,dlogxw,egexjc,egexjf,egexs,eh,ehfac,
     $ eps100,fo2,fo2lg,fsort,fugac,fugalg,gmmatr,ibetmx,ibpxmx,idelmx,
     $ iebal,ier,iern1,iern2,ietmax,ifrn1,ifrn2,igas,igstak,iindx1,
     $ ilrn1,ilrn2,imrn1,imrn2,iodb,ipivot,ipndx1,irdxc3,istack,iter,
     $ itermx,ixbasp,ixrn1,ixrn2,jcsort,jern1,jern2,jetmax,jflag,jgext,
     $ jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jssort,jstack,kbt,
     $ kction,kdim,kelect,khydr,kmax,kmt,km1,ko2gaq,kwater,kxt,kx1,
     $ loph,losp,lsort,mgext,moph,mosp,mrgexs,mtb,narn1,narn2,nbasp,
     $ nbt,nbtmax,nbw,ncmpr,ncosp,ndrs,ndrsmx,ndrsr,negbfc,negdfc,
     $ negxfc,nelect,nern1,nern2,netmax,ngexsa,ngext,ngrn1,ngrn2,ngt,
     $ ngtmax,nhydr,nodbmx,noibfc,noutpt,no2gaq,npconv,nphasx,npobfc,
     $ npodfc,npt,nptmax,nredox,nst,nstmax,nsts,nstsmx,nstsr,ntfx,
     $ ntfxmx,ntfxt,nttyo,nxtmax,omega,omeglg,press,qcacf,qcbeta,
     $ qredox,qxbarw,q6mode,rhsvec,rlxfac,screwd,screwn,sigmam,sigmmo,
     $ tfx,ubetmx,ubbig,ubneg,ugexj,ugexmo,uphase,uspec,uzvec1,weight,
     $ xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xlke,xlks,zchar,zgexj,
     $ zvclg1,zvec1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ier .gt. 0) then
c
c       Iteration has been stopped.
c
        if (qcbeta .and. qcgam) then
c
c         Have pseudo-convergence. Reset ier to 0.
c
          ier = 0
          if (iodb(4) .ge. 1) then
            write (noutpt,1310)
 1310       format('  Hybrid Newton-Raphson iteration has',
     $      ' pseudo-converged.')
            do kcol = 1,kdim
              if (abs(delvec(kcol)) .gt. toldl) then
c
c               Calling sequence substitutions:
c                 uzvec1(kcol) for unam48
c
                call fmspnx(jlen,uzvec1(kcol),uspn56)
                write (noutpt,1320) uspn56(1:jlen),delvec(kcol)
 1320           format(11x,'delvec(',a,') = ',g12.5)
              endif
            enddo
          endif
          go to 999
        else
c
c         Iteration has failed.
c
          write (noutpt,1330) iter
 1330     format('  Hybrid Newton-Raphson iteration has gone sour',
     $    ' (iter= ',i3,').')
          go to 999
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Recompute sigmam, fxi, and fje. Recompute the activity
c     coefficients. Then recompute the concentrations of dependent
c     species. After that, recalculate the Newton-Raphson residual
c     functions. Finally, go back and see whether or not to do another
c     iteration.
c
      if (iodb(4) .ge. 2) write (noutpt,1340)
 1340 format(/3x,'--- Post-Newton-Raphson update of activity',
     $ ' coefficients ---',/)
c
c     Determine the maximum allowed change factors for "Sigma m",
c     etc. (chfsgm) and the log activity coefficients for aqueous
c     species (chfacf).
c
      chfsgm = 1.5
      if (sigmam .le. 2.e-1) chfsgm = 5.
      if (sigmam .le. 1.e-2) chfsgm = 10.
      chfacf = 0.50
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
c     Recalculate the Newton-Raphson residual functions.
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
      go to 100
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
