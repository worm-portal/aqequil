      subroutine path(aadh,aadhh,aadhv,aaphi,abdh,abdhh,abdhv,
     $ abdot,abdoth,abdotv,adadhh,adadhv,adbdhh,adbdhv,adbdth,
     $ adbdtv,adh,adhh,adhv,afcnst,aftarg,al10,aphi,apx,atwt,
     $ avcnst,awmax,awmaxi,awmin,awmini,azero,bdh,bdhh,bdhv,bdot,
     $ bdoth,bdotv,bpx,cbsr,cbsri,cco2,cdac,cdrs,cdrsd,cdrsx,
     $ cegexs,cesr,cesri,cess,cpgexs,cscale,csigma,csts,dadhh,
     $ dadhv,dbdhh,dbdhv,dbdth,dbdtv,dlaplo,dlaprn,dleplo,dleprn,
     $ dlhplo,dlhprn,dloplo,dloprn,dltplo,dltpll,dltprl,dltprn,
     $ dlxdmp,dlxmax,dlxmin,dlxmx0,dlxplo,dlxpll,dlxprl,dlxprn,
     $ eact,egers,egersi,egexjf,ehfac,ehmax,ehmaxi,ehmin,ehmini,
     $ elecsr,electr,eps100,farad,fkrc,hact,iact,iapxt,iaqsln,
     $ ibpxt,ibsrti,ielam,iern1,iern2,iesrti,ifcphi1,ifcphi2,
     $ ifnnn,ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igerti,
     $ iindx1,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,
     $ ilrn2,ilzeta,imech,imrn1,imrn2,insgf,iodb,iopg,iopr,iopt,
     $ ipch,ipndx1,ipcv,irang,itermx,ixrn1,ixrn2,ixrti,izmax,
     $ jcode,jffg,jflag,jflagd,jflgi,jgerti,jpflag,jpress,jptffl,
     $ jreac,jsflag,jsitex,jsol,jtemp,kbt,kct,kdim,kelect,khydr,
     $ khydx,km1,kmt,ko2gaq,kprs,krdxsp,ksplmx,ksppmx,kstpmx,kwater,
     $ kxmod,kx1,kxt,loph,losp,modr,moffg,moph,morr,mosp,mprph,
     $ mprphi,mprsp,mprspi,mrgers,mrgexs,mtb,mtbi,mtbaq,mtbaqi,
     $ mte,mteaq,mwtrc,mwtsp,narn1,narn2,narxt,nat,nbasp,nbaspd,
     $ nbaspi,nbaspx,nbkupa,nbkupb,nbt,nbtd,nbti,nbw,nchlor,
     $ ncmpr,nct,ndac,ndact,nelect,nern1,nern2,ness,nessr,net,
     $ ndrs,ndrsd,ndrsx,ndrsr,ndrsrd,ndrsrx,nert,newin,nffg,
     $ nfrn1,nfrn2,ngrn1,ngrn2,ngt,nhydr,nhydx,nllnmx,nlrn1,
     $ nlrn2,nlt,nmrn1,nmrn2,nmrt,nmt,nobswt,noutpt,no2gaq,npchk,
     $ nphasx,nprob,nprpti,nprsti,npslmx,npt,nrct,nrdxsp,nrk,
     $ nrndex,nsbswt,nsk,nsrt,nsslmx,nst,nsts,nstsr,ntabx,ntf1,
     $ ntf1t,ntf2,ntf2t,ntitl1,ntitl2,ntitld,ntpr,ntprt,ntrymx,
     $ nttyo,nxmod,nxopex,nxopt,nxridx,nxrn1,nxrn2,nxrt,nxt,o2max,
     $ o2maxi,o2min,o2mini,phmax,phmaxi,phmin,phmini,prcinf,
     $ press,pressb,pressd,pressi,ptk,qcnpre,qcntmp,qdwipp,
     $ qecon,qgexsh,qhawep,qoptmz,qpit75,qredox,qscon,qtatxt,
     $ rconst,rcnstv,rk,rkb,rtcnst,rxbar,rxbari,sfcar,smp100,
     $ sscrew,ssfcar,tempc,tempcb,tempcd,tempci,tempcu,tempk,
     $ tf1,tf2,timemx,time1,timmxi,tistrt,tistti,trkb,ttk,tolaft,
     $ tolbt,toldl,tolsat,tolsst,tolxsf,tolxst,tolxsu,uaqsln,
     $ ubmtbi,ubsri,ucxri,udac,uelem,uesri,uffg,ufixf,ugerji,
     $ ugermo,ugersi,uinfor,ulinex,uobsw,uphase,uplatm,uprphi,
     $ uprspi,ureac,usbsw,uspec,usteq6,utitl1,utitl2,utitld,
     $ uveeq6,uxcat,uxmod,uxopex,uxopt,uzvec1,uzveci,vosp0,vreac,
     $ wfac,xgers,xgersi,ximax,ximaxi,xistrt,xistti,xi1,xlkffg,
     $ xlkmod,zchar,zchcu6,zchsq2,zklgmn,zklogl,zklogu,zvclgi,
     $ zvclg1,zvec1)
c
c     This subroutine is the control subroutine for tracing a reaction
c     path. It determines a variable step size, which is in terms of an
c     overall reaction progress variable. The corresponding time step,
c     if defined, is computed as a dependent variable. This subroutine
c     also detects and responds to various events, such as saturation
c     of irreversible reactions and exhaustion of reactants.
c
c     The mathematical treatment of the algebraic relations is
c     analogous to the predictor-corrector method of integrating
c     differential equations. Finite difference expressions of order
c     up to nordmx are used as predictor functions. These are
c     manipulated into the form of truncated Taylor's series.
c     EQ6/eqcalc.f is used as the "corrector." This subroutine does
c     equilbrium calculations using a hybrid Newton-Raphson algorithm.
c     Note that it corrects to satisfy algebraic equations, not
c     differential equations.
c
c     In kinetic mode (iopt(2) > 0), rate laws are integrated using an
c     actual predictor-corrector algorithm. Algebraic equations continue
c     to be dealt with as described above. Rate law integration is
c     done using finite-difference expressions of the same order
c     as those used to represent algebraic variables.
c
c     Apart from no time/time, there are three distinct calculational
c     modes:
c
c       (1) Normal mode- The step size is constrained to keep the
c         predictor functions fairly accurate. There is less burden
c         on the Newton-Raphson thermodynamic calculations.
c         Normal mode requires longer run times and is more
c         expensive than economy mode. Normal mode must be used
c         for some kinds of problems. Examples include the cases
c         iopt(2) > 0 (kinetic mode) and iopt(1) = 2, (flow-through
c         open system mode). Setting iopt(13) = 0 insures that
c         normal mode will be selected. It will automatically be
c         selected if economy or super-economy modes are not valid
c         options.
c
c       (2) Economy mode- The step size is allowed to become large
c         more quickly than in normal mode. The predictor functions
c         are limited to second order. There is no attempt to
c         constrain the step size in order to keep these functions
c         accurate. The burden on the equilibrium calculations is
c         heavier than in normal mode, because the initial values of
c         the unknowns will generally be farther off the mark.
c         Economy mode is intended to be faster and cheaper than
c         normal mode. The useful information density along the
c         reaction path with economy mode is essentially equivalent to
c         that obtained with normal mode. Use iopt(13) = 1 to select
c         economy mode.
c
c       (3) Super economy mode- This is a special from of economy
c         mode. The step size is typically large. It defaults to
c         dlxprn, the linear print interval. Phase boundaries are
c         ignored. The order of the finite differences is restricted
c         to zero. Super economy mode provides much less information
c         density along the reaction path than does economy mode or
c         normal mode. Use iopt(13) = 2 to select super economy mode.
c
c     The step size is controlled by various constraints. One of these
c     is the estimated accuracy of the predictor functions. If the
c     system is rapidly changing, delxi will be small and the run time
c     per unit reaction progress will be long. This is generally the
c     case at the start of a reaction path simulation. When the system
c     is changing slowly with respect to reaction progress, delxi may
c     become large. Other controls on the step size include:
c
c        (1) The initial step size (for order zero)
c        (2) An arbitrary upper limit on step size
c        (3) The print interval requirements
c        (4) The plot interval requirements
c        (5) An arbitrary dump interval (flow-through model only)
c        (6) The exhaustion of any reactant
c        (7) The appearance or disappearance of any secondary
c            phase (unless iopt(3) .ge. 1)
c        (8) A maximum in the mass of any secondary mineral
c            (fluid-centered flow-through model only)
c        (9) Crossing 100 degrees C, if temperature is changing with
c            reaction progress
c       (10) The terminal value of reaction progress
c       (11) The terminal value of time, or infinite time
c
c     In the flow-through open system model (iopt(1) = 2), masses of
c     secondary minerals are shifted to the physically removed
c     subsystem after any of (5) through (11) above.
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
      include 'eqlib/eqlwd.h'
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nbkupa,nbkupb,newin,noutpt,ntabx,nttyo
c
      integer jcode(nrctmx),jreac(nrctmx),nrndex(nrctmx),nxridx(nrctmx)
c
      integer iact(imchmx,2,nrctmx),imech(2,nrctmx),
     $ ndac(ndctmx,imchmx,2,nrctmx),ndact(imchmx,2,nrctmx),
     $ nrk(2,nrctmx),nsk(nrctmx)
c
      integer ibsrti(nsrtmx),iesrti(nsrtmx),igerti(jetmax,nertmx),
     $ ixrti(nxrtmx),jflgi(nbtmax),jgerti(nertmx)
c
      integer iapxt(nxtmax),ibpxt(nxtmax),iindx1(kmax),
     $ ipndx1(kmax),insgf(natmax),iodb(nodbmx),iopg(nopgmx),
     $ iopr(noprmx),iopt(noptmx),jffg(nffgmx),jflag(nstmax),
     $ jflagd(nstmax),jpflag(nptmax),jsflag(nstmax),jsitex(nstmax),
     $ jsol(nxtmax),kxmod(nxmdmx)
c
      integer narxt(ntprmx),nbasp(nbtmax),nbaspd(nbtmax),
     $ nbaspi(nbtmax),nbaspx(nbtmax),ncmpr(2,nptmax),ness(nessmx),
     $ nessr(2,nstmax),ndrs(ndrsmx),ndrsd(ndrsmx),ndrsx(ndrsmx),
     $ ndrsr(2,nstmax),ndrsrd(2,nstmax),ndrsrx(2,nstmax),npchk(nptmax),
     $ nphasx(nstmax),nsts(nstsmx),nstsr(2,nstmax),ntf1(ntf1mx),
     $ ntf2(ntf2mx)
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
      integer kprs,nbti,nprpti,nprsti
c
      integer iaqsln,ielam,igas,ipch,ipcv,itermx,irang,izmax,jpress,
     $ jptffl,jtemp,kbt,kct,kdim,kelect,khydr,khydx,km1,kmt,ko2gaq,
     $ krdxsp,ksplmx,ksppmx,kstpmx,kwater,kx1,kxt,nbtd,nbw,nchlor,
     $ nelect,nert,nffg,nhydr,nhydx,nllnmx,nmrt,nobswt,no2gaq,nprob,
     $ npslmx,nrct,nrdxsp,nsbswt,nsrt,nsslmx,ntf1t,ntf2t,ntitl1,
     $ ntitl2,ntitld,ntpr,ntprt,ntrymx,nxmod,nxopex,nxopt,nxrt
c
      logical qcnpre,qcntmp,qdwipp,qecon,qgexsh,qhawep,qoptmz,qpit75,
     $ qredox,qscon,qstart,qtatxt,qvhfxi,qvlsow
c
      character(len=nllnmx) ulinex
      character(len=80) utitl1(ntitmx),utitl2(ntitmx),utitld(ntidmx)
      character(len=48) ubmtbi(nbtmax),uprspi(nprsmx),uzveci(kmax)
      character(len=48) uobsw(2,nbtmax),usbsw(2,nbtmax),uspec(nstmax),
     $ uxmod(nxmdmx),uzvec1(kmax)
      character (len=24) ubsri(nbt1mx,nsrtmx),ucxri(iktmax,nxrtmx),
     $ ugersi(ietmax,jetmax,nertmx),uprphi(nprpmx)
      character(len=24) ugermo(nertmx),ureac(nrctmx)
      character(len=24) uffg(nffgmx),uphase(nptmax),uxcat(nxopmx),
     $ uxopex(nxpemx)
      character(len=24) udac(ndctmx,imchmx,2,nrctmx)
      character(len=24) uaqsln
      character(len=8) uesri(nctmax,nsrtmx),ugerji(jetmax,nertmx)
      character(len=8) uelem(nctmax),uxopt(nxopmx)
      character(len=8) ufixf,uinfor,uplatm,usteq6,uveeq6
c
      real(8) aadh(narxmx,ntprmx),aadhh(narxmx,ntprmx),
     $ aadhv(narxmx,ntprmx),aaphi(narxmx,ntprmx),
     $ abdh(narxmx,ntprmx),abdhh(narxmx,ntprmx),
     $ abdhv(narxmx,ntprmx),abdot(narxmx,ntprmx),
     $ abdoth(narxmx,ntprmx),abdotv(narxmx,ntprmx)
c
      real(8) dadhh(ipchmx),dadhv(ipcvmx),dbdhh(ipchmx),dbdhv(ipcvmx),
     $ dbdth(ipchmx),dbdtv(ipcvmx)
c
      real(8) cbsri(nbt1mx,nsrtmx),cesri(nctmax,nsrtmx),
     $ egersi(ietmax,jetmax,nertmx),mprphi(nprpmx),mprspi(nprsmx),
     $ mtbi(nbtmax),mtbaqi(nbtmax),rxbari(iktmax,nxrtmx),
     $ xgersi(ietmax,jetmax,nertmx),zvclgi(kmax)
c
      real(8) adadhh(narxmx,ntprmx,ipchmx),
     $ adadhv(narxmx,ntprmx,ipcvmx),adbdhh(narxmx,ntprmx,ipchmx),
     $ adbdhv(narxmx,ntprmx,ipcvmx),adbdth(narxmx,ntprmx,ipchmx),
     $ adbdtv(narxmx,ntprmx,ipcvmx)
c
       real(8) cdac(ndctmx,imchmx,2,nrctmx),csigma(imchmx,2,nrctmx),
     $ eact(imchmx,2,nrctmx),fkrc(nrctmx),hact(imchmx,2,nrctmx),
     $ rkb(imchmx,2,nrctmx),trkb(imchmx,2,nrctmx)
c
      real(8) cbsr(nbt1mx,nsrtmx),cesr(nctmax,nsrtmx),
     $ egers(ietmax,jetmax,nertmx),elecsr(nsrtmx),modr(nrctmx),
     $ mrgers(ietmax,jetmax,nertmx),morr(nrctmx),mwtrc(nrctmx),
     $ rxbar(iktmax,nxrtmx),sfcar(nrctmx),ssfcar(nrctmx),vreac(nrctmx),
     $ xgers(ietmax,jetmax,nertmx)
c
      real(8) apx(iapxmx,nxtmax),atwt(nctmax),
     $ azero(natmax),bpx(ibpxmx,nxtmax),cco2(5),
     $ cegexs(ietmax,jetmax,netmax),cpgexs(ietmax,jetmax,netmax),
     $ cess(nessmx),cdrs(ndrsmx),cdrsd(ndrsmx),cdrsx(ndrsmx),
     $ cscale(nstmax),csts(nstsmx),
     $ egexjf(jetmax,netmax),loph(nptmax),losp(nstmax),
     $ moph(nptmax),mosp(nstmax),moffg(nffgmx),mprph(nptmax),
     $ mprsp(nstmax),mrgexs(ietmax,jetmax,netmax),mtb(nbtmax),
     $ mtbaq(nbtmax),mte(nctmax),mteaq(nctmax),mwtsp(nstmax),
     $ ptk(nptkmx),rk(imchmx,2,nrctmx),
     $ sscrew(nsscmx),tempcu(ntprmx),tf1(ntf1mx),tf2(ntf2mx),
     $ ttk(nttkmx),vosp0(nstmax),wfac(iktmax,nxtmax),xlkffg(nffgmx),
     $ xlkmod(nxmdmx),zchar(nstmax),zchcu6(nstmax),zchsq2(nstmax),
     $ zvclg1(kmax),zvec1(kmax)
c
      real(8) adh,adhh,adhv,aphi,bdh,bdhh,bdhv,bdot,bdoth,bdotv
c
      real(8) afcnst,aftarg,al10,avcnst,awmax,awmaxi,awmin,awmini,
     $ dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,
     $ dltplo,dltprl,dltprn,dlxdmp,dlxmax,dlxmin,dlxmx0,dlxpll,dlxplo,
     $ dlxprl,dlxprn,ehfac,ehmax,ehmaxi,ehmin,ehmini,electr,eps100,
     $ farad,o2max,o2maxi,o2min,o2mini,phmax,phmaxi,phmin,phmini,
     $ prcinf,press,pressb,pressd,pressi,rcnstv,rconst,rtcnst,
     $ smp100,tempc,tempcb,tempcd,tempci,tempk,timemx,time1,timmxi,
     $ tistrt,tistti,tolaft,tolbt,toldl,tolsat,tolsst,tolxsf,tolxst,
     $ tolxsu, ximax,ximaxi,xistrt,xistti,xi1,zklgmn,zklogl,zklogu
c
c-----------------------------------------------------------------------
c
c     Local variable declarations with static global dimensioning.
c
      character(len=32) uxtype(jso_par)
c
      real(8) egexjc(jet_par,net_par),egexpa(net_par),egexpc(net_par),
     $ egexs(iet_par,jet_par,net_par),egexw(ket_par,net_par),
     $ xgexw(ket_par,net_par)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations with variable global dimensioning.
c
      integer, dimension(:), allocatable :: ibswx,idirec
c
      integer, dimension(:), allocatable :: igstak,ipivot,ipivtr,
     $ istack,ixbasp,jcsort,jgstak,jgsort,jjsort,jssort,jstack,
     $ kction
c
      integer, dimension(:), allocatable :: iemop,iemop0,iemos,iemos0,
     $ iexr,iindx0,ipndx0,jexr,jreac0,jsca,jscr
c
      integer, dimension(:,:), allocatable :: ncmpe,ncmpe0
c
      integer, dimension(:,:,:), allocatable :: ndactb
c
      integer, dimension(:,:,:,:), allocatable :: ndacb
c
      logical, dimension(:), allocatable :: qxknph
c
      character(len=48), dimension(:), allocatable :: uzvec0
      character(len=8), dimension(:), allocatable :: ulbeta,uldel
c
      real(8), dimension(:,:), allocatable :: akmat0,akmat1,daffp0,
     $ dafrc0,demop0,demos0,drer0,drer0s,dzvc0,dzvc0s,fdafm1,fdaf0,
     $ fdarm1,fdar0,fdpem1,fdpe0,fdrem1,fdre0,fdre1,fdrrm1,fdrr0,
     $ fdrr1,fdsem1,fdse0,fdzvm1,fdzv0
c
      real(8), dimension(:), allocatable :: affp,affp0,affs,afrc0,
     $ afrc1,afrcp,cjbasp,cnufac,ctb,daw0,deh0,do20,dph0,drir0,drir0s,
     $ dxsm00,dxsm10,dxsm11,dxval0,d1zvc1,d1emp1,d2emp1,ehrc,emop,
     $ emop0,emos,emos0,fdaw0,fdawm1,fdeh0,fdehm1,fdo20,fdo2m1,
     $ fdph0,fdphm1,fdrim1,fdri0,fdri1,fo2lrc
c
      real(8), dimension(:), allocatable :: acflg,acflgo,acflg0,act,
     $ actlg,ahrc,affpd,affsd,alpha,amtb,a3bars,beta,betao,cdrtw,cdrw,
     $ conc,conclg,cteaq,delvco,delvec,dlogxw,fsort,fugac,fugalg,
     $ lsort,ppmwb,ppmwe,rhsvec,rreacn,rreac0,rreac1,rrelrp,rrelr0,
     $ rrelr1,sfcar0,sidrph,sidrsp,wodr,worr,xbar,xbarlg
c
      real(8), dimension(:), allocatable :: modr0,mophg,mophj,mopht,
     $ moph0,morr0,mospg,mospj,mospt,mosp0,mtb0,perc,voph,vophg,
     $ vophj,vopht,vosp,vospg,vospj,vospt,weight,woph,wophg,wophj,
     $ wopht,wosp,wospg,wospj,wospt,xirct,xirct0,zvclg0,zvec0
c
      real(8), dimension(:), allocatable :: alphar,betar,delvcr,rhsvcr,
     $ dvjdte
c
      real(8), dimension(:), allocatable :: hhcvec,xhcvec
c
      real(8), dimension(:,:), allocatable :: aamatr,gmmatr
      real(8), dimension(:,:), allocatable :: armatr,grmatr
      real(8), dimension(:,:), allocatable :: aimatr,mmmatr,sgmatr,
     $ xxmatr,xymatr,rrxfi1
c
      real(8), dimension(:,:,:,:), allocatable :: cdacb
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer nbkupn
c
      integer nrct1
c
      integer i,iavkdm,ibtrmx,idlrmx,ier,iexrt,inmax,iter,iwdh2o,
     $ iwnffg,j,jcut,jexrt,jexrtx,jordlm,jsawth,jscat,jscatx,jscrt,
     $ jscrtx,j2,j3,kaft1,kcol,kdim0,kly,km10,kmt0,kord,kordlm,kordp1,
     $ kordsv,kpsat,kpsst,kstep,kstpab,kstppl,kstppr,kstpmn,kstpze,kx10,
     $ kxt0,kzmax,n,naft1,nb,nc,ncorr,ncut,ndelay,ne,ng,nlwffg,nmax,
     $ nord,nordr,nordrs,nordsv,nordz,nordzs,np,npe,npet,npet0,npts,
     $ nrc,nr1,nr2,ns,nsawth,nse,nset,nset0,nswtch,ntpr0,nwdh2o,
     $ nweope,nwndpc,nwnffg
c
      integer iarmxn,ilnobl
c
      logical qabswx,qadjdx,qaflip,qaft1,qbassw,qbseqc,qbswx,qbye,
     $ qconst,qodeok,qftpr2,qdmpr1,qdmpr2,qdump,qhcon,qmin,qmod,qmod1,
     $ qmod2,qphcl,qplaw0,qplaw1,qpleh0,qpleh1,qplolt,qplolx,qplott,
     $ qplotx,qplo20,qplo21,qplph0,qplph1,qpr,qpraw0,qpraw1,qpreh0,
     $ qpreh1,qprnlt,qprnlx,qprntt,qprntx,qpro20,qpro21,qprph0,qprph1,
     $ qrapch,qreax,qreq,qrho,qriinf,qrpcfl,qsawth,qshoot,qskip,
     $ qstabl,qstabr,qstabz,qstop,qstopx,qtplo,qtprn,qtrch,qtvchk,
     $ qwhcfa,qx,qxbarw,qzdump,qzplot,qzprnt,qztayl,q6mode
c
      character(len=48) ubacmx,ubgamx
      character(len=24) ustr24,ux24
      character(len=16) ux16,ux16a,ux16b
      character(len=8) upkfor
c
      real(8) abar,acfw,acfwlg,actw,actwlc,actwlg,adel,adlzlg,adx,aft1,
     $ aft0,aftm1,ah,ahmes,ahnbs,alk,alki,alk1,alk2,av,avkdim,avdlxi,
     $ awstrt,aw0,aw0plo,aw0prn,aw1,aw1plo,aw1prn,a3bar,bacfmx,bbig,
     $ betamx,bgamx,bneg,deltim,delxi,delxia,dlxilm,dlximx,dlxipl,
     $ dlxipr,dlxisv,dlxis2,dlxlim,dlxode,dlxtmx,dlxtpr,dlxtpl,dvoso,
     $ dwoso,dx,dxdmp,dxsave,dx1,dx2,eh,ehmes,ehnbs,ehstrt,eh0,eh0plo,
     $ eh0prn,eh1,eh1plo,eh1prn,fdlim,fdx,fje,fjeo,fje0,fjest,fo2,
     $ fo2lg,fo2lg0,fo2lg1,fxi,fxio,fxi0,fxist,fxprpl,morrw1,mx,lprcin,
     $ lx,osc,oscst,omega,omeglg,o2strt,o20plo,o20prn,o21plo,o21prn
c
      real(8) pch,pe,pemes,penbs,ph,phcl,phmes,phnbs,phstrt,ph0,ph0plo,
     $ ph0prn,ph1,ph1plo,ph1prn,prminf,rho,rirecp,rirec0,rirec1,rx0,
     $ rx1,rxx,scale,scalim,scfcr,scfcrs,scfcz,scfczs,scnsti,scnstd,
     $ screwd,sigmam,sigmmo,sigmm0,sigmst,tempc0,thours,time0,tiplol,
     $ tiplot,tiplxx,tiprnl,tiprnt,tiprxx,tistrd,tistry,tistsv,tmins,
     $ tolsar,tolsrr,tyears,vodrt,tdays,vosoct,vosol,wfh2o,wftds,
     $ whcfac,wkgh2o,wwstrt,wkgsol,wkgwi,wodrt,woh2o,worrt,wosoct,
     $ wosol,wotds,wx,xbarw,xbarwc,xbrwlc,xbrwlg,xf,xidump,xilim,
     $ xim1,xiplol,xiplot,xiprnl,xiprnt,xistsv,xi0,xlf,xval0,xx
c
      real(8) dxe0mx,dxe1mx,dxe0pl,dxe1pl,dxe0pr,dxe1pr,dxh0mx,dxh1mx,
     $ dxh0pl,dxh1pl,dxh0pr,dxh1pr,dxo0mx,dxo1mx,dxo0pl,dxo1pl,
     $ dxo0pr,dxo1pr,dxw0mx,dxw1mx,dxw0pl,dxw1pl,dxw0pr,dxw1pr
c
      real(8) btrfnc,btrmax,btrmxo,dlrfnc,dlrmax,dlrmxo
c
      real(8) mlmrra,mrmlra,rhoc,rhowc,tdsgks,tdsglw,tdspkc,tdsplc
c
      real(8) fctrl,texp,tlg
c
c-----------------------------------------------------------------------
c
c     Variable declarations: Local data needed to assist in writing
c     an EQ6 pickup file. These variables do not appear on that file
c     itself.
c
      integer ibsrt1,iesrt1
c
      character(len=24), dimension(:), allocatable :: ubsr1
      character(len=8), dimension(:), allocatable :: uesr1
c
      character(len=24) ureac1
c
      real(8), dimension(:), allocatable :: cbsr1,cesr1
c
c-----------------------------------------------------------------------
c
c     Names of solid solution models.
c
      data uxtype(1) /'Ideal solution                  '/
      data uxtype(2) /'Binary, third-order Maclaurin   '/
      data uxtype(3) /'Binary, parabolic Maclaurin     '/
      data uxtype(4) /'Binary, cubic Maclaurin (P,T)   '/
      data uxtype(5) /'Binary, Guggenheim (T)          '/
      data uxtype(6) /'Ternary, regular                '/
      data uxtype(7) /'Plagioclase (Newton et al. 1980)'/
c
c     Limit on the number of points of reaction progress beyond the
c     starting point at which warnings are issued regarding fugacities
c     not being at specified values.
c
      data nlwffg / 3 /
c
c     The variable "alk" is not used in this subroutine except to
c     satisfy certain EQLIB subroutine calls.
c
      data alk /0./
c
      data fxprpl /0.2/
c
      data qxbarw/.false./
c
c-----------------------------------------------------------------------
c
c     Allocate additional arrays needed to compute the reaction path.
c
      ALLOCATE(ctb(nbtmax))
      ALLOCATE(ppmwb(nbtmax))
      ALLOCATE(ehrc(nbtmax))
      ALLOCATE(cjbasp(nbtmax))
      ALLOCATE(ibswx(nbtmax))
      ALLOCATE(ixbasp(nbtmax))
      ALLOCATE(dlogxw(nbtmax))
c
      ALLOCATE(aamatr(kmax,kmax))
      ALLOCATE(gmmatr(kmax,kmax))
c
      ALLOCATE(iindx0(kmax))
      ALLOCATE(ipivot(kmax))
      ALLOCATE(ipndx0(kmax))
c
      ALLOCATE(alpha(kmax))
      ALLOCATE(beta(kmax))
      ALLOCATE(betao(kmax))
      ALLOCATE(delvco(kmax))
      ALLOCATE(delvec(kmax))
      ALLOCATE(d1zvc1(kmax))
      ALLOCATE(rhsvec(kmax))
      ALLOCATE(zvclg0(kmax))
      ALLOCATE(zvec0(kmax))
c
      ALLOCATE(uzvec0(kmax))
      ALLOCATE(uldel(kmax))
      ALLOCATE(ulbeta(kmax))
c
      ALLOCATE(a3bars(natmax))
c
      ALLOCATE(kction(nbtmax))
c
      ALLOCATE(ahrc(nbtmax))
      ALLOCATE(amtb(nbtmax))
      ALLOCATE(fo2lrc(nbtmax))
      ALLOCATE(mtb0(nbtmax))
      ALLOCATE(perc(nbtmax))
c
      ALLOCATE(igstak(ngtmax))
      ALLOCATE(jgstak(ngtmax))
      ALLOCATE(jgsort(ngtmax))
c
      ALLOCATE(cteaq(nctmax))
      ALLOCATE(ppmwe(nctmax))
c
      ALLOCATE(fsort(ngtmax))
      ALLOCATE(fugac(ngtmax))
      ALLOCATE(fugalg(ngtmax))
c
      ALLOCATE(ncmpe(2,npetmx))
      ALLOCATE(ncmpe0(2,npetmx))
c
      ALLOCATE(dzvc0(nrd1mx,kmax))
      ALLOCATE(dzvc0s(nrd1mx,kmax))
      ALLOCATE(fdzvm1(nrd1mx,kmax))
      ALLOCATE(fdzv0(nrd1mx,kmax))
c
      ALLOCATE(akmat0(nrd1mx,nrd1mx))
      ALLOCATE(akmat1(nrd1mx,nrd1mx))
c
      ALLOCATE(daffp0(nordmx,nptmax))
      ALLOCATE(fdafm1(nordmx,nptmax))
      ALLOCATE(fdaf0(nordmx,nptmax))
c
      ALLOCATE(dafrc0(nordmx,nrctmx))
      ALLOCATE(drer0(nrd1mx,nrctmx))
      ALLOCATE(drer0s(nrd1mx,nrctmx))
      ALLOCATE(fdarm1(nordmx,nrctmx))
      ALLOCATE(fdar0(nordmx,nrctmx))
      ALLOCATE(fdrem1(nordmx,nrctmx))
      ALLOCATE(fdre0(nordmx,nrctmx))
      ALLOCATE(fdre1(nordmx,nrctmx))
      ALLOCATE(fdrrm1(nrd1mx,nrctmx))
      ALLOCATE(fdrr0(nrd1mx,nrctmx))
      ALLOCATE(fdrr1(nrd1mx,nrctmx))
c
      ALLOCATE(demop0(nordmx,npetmx))
      ALLOCATE(fdpem1(nordmx,npetmx))
      ALLOCATE(fdpe0(nordmx,npetmx))
c
      ALLOCATE(demos0(nordmx,nsetmx))
      ALLOCATE(fdsem1(nordmx,nsetmx))
      ALLOCATE(fdse0(nordmx,nsetmx))
c
      ALLOCATE(daw0(nordmx))
      ALLOCATE(fdawm1(nordmx))
      ALLOCATE(fdaw0(nordmx))
      ALLOCATE(deh0(nordmx))
      ALLOCATE(fdehm1(nordmx))
      ALLOCATE(fdeh0(nordmx))
      ALLOCATE(do20(nordmx))
      ALLOCATE(fdo2m1(nordmx))
      ALLOCATE(fdo20(nordmx))
      ALLOCATE(dph0(nordmx))
      ALLOCATE(fdphm1(nordmx))
      ALLOCATE(fdph0(nordmx))
c
      ALLOCATE(dxsm00(nrd1mx))
      ALLOCATE(dxsm10(nrd1mx))
      ALLOCATE(dxsm11(nrd1mx))
c
      ALLOCATE(drir0(nrd1mx))
      ALLOCATE(drir0s(nrd1mx))
c
      ALLOCATE(fdrim1(nrd1mx))
      ALLOCATE(fdri0(nrd1mx))
      ALLOCATE(fdri1(nrd1mx))
c
      ALLOCATE(dxval0(nrd1mx))
c
      ALLOCATE(iemop(npetmx))
      ALLOCATE(iemop0(npetmx))
c
      ALLOCATE(d1emp1(npetmx))
      ALLOCATE(d2emp1(npetmx))
      ALLOCATE(emop(npetmx))
      ALLOCATE(emop0(npetmx))
c
      ALLOCATE(qxknph(nptmax))
c
      ALLOCATE(affp(nptmax))
      ALLOCATE(affp0(nptmax))
      ALLOCATE(affpd(nptmax))
      ALLOCATE(mophg(nptmax))
      ALLOCATE(mophj(nptmax))
      ALLOCATE(mopht(nptmax))
      ALLOCATE(moph0(nptmax))
      ALLOCATE(sidrph(nptmax))
      ALLOCATE(voph(nptmax))
      ALLOCATE(vophj(nptmax))
      ALLOCATE(vophg(nptmax))
      ALLOCATE(vopht(nptmax))
      ALLOCATE(woph(nptmax))
      ALLOCATE(wophg(nptmax))
      ALLOCATE(wophj(nptmax))
      ALLOCATE(wopht(nptmax))
c
      ALLOCATE(iexr(nrctmx))
      ALLOCATE(jexr(nrctmx))
      ALLOCATE(jreac0(nrctmx))
      ALLOCATE(jsca(nrctmx))
      ALLOCATE(jscr(nrctmx))
c
      ALLOCATE(afrc0(nrctmx))
      ALLOCATE(afrc1(nrctmx))
      ALLOCATE(afrcp(nrctmx))
      ALLOCATE(idirec(nrctmx))
      ALLOCATE(sfcar0(nrctmx))
      ALLOCATE(modr0(nrctmx))
      ALLOCATE(morr0(nrctmx))
      ALLOCATE(rreacn(nrctmx))
      ALLOCATE(rreac0(nrctmx))
      ALLOCATE(rreac1(nrctmx))
      ALLOCATE(rrelrp(nrctmx))
      ALLOCATE(rrelr0(nrctmx))
      ALLOCATE(rrelr1(nrctmx))
      ALLOCATE(wodr(nrctmx))
      ALLOCATE(worr(nrctmx))
      ALLOCATE(xirct(nrctmx))
      ALLOCATE(xirct0(nrctmx))
c
      ALLOCATE(rrxfi1(imchmx,nrctmx))
c
      ALLOCATE(iemos(nsetmx))
      ALLOCATE(iemos0(nsetmx))
c
      ALLOCATE(emos(nsetmx))
      ALLOCATE(emos0(nsetmx))
c
      ALLOCATE(istack(nstmax))
      ALLOCATE(jstack(nstmax))
      ALLOCATE(jcsort(nstmax))
      ALLOCATE(jjsort(nstmax))
      ALLOCATE(jssort(nstmax))
c
      ALLOCATE(acflg(nstmax))
      ALLOCATE(acflgo(nstmax))
      ALLOCATE(acflg0(nstmax))
      ALLOCATE(act(nstmax))
      ALLOCATE(actlg(nstmax))
      ALLOCATE(affs(nstmax))
      ALLOCATE(affsd(nstmax))
      ALLOCATE(cdrtw(nstmax))
      ALLOCATE(cdrw(nstmax))
      ALLOCATE(cnufac(nstmax))
      ALLOCATE(conc(nstmax))
      ALLOCATE(conclg(nstmax))
      ALLOCATE(lsort(nstmax))
      ALLOCATE(mospg(nstmax))
      ALLOCATE(mospj(nstmax))
      ALLOCATE(mospt(nstmax))
      ALLOCATE(mosp0(nstmax))
      ALLOCATE(sidrsp(nstmax))
      ALLOCATE(vosp(nstmax))
      ALLOCATE(vospg(nstmax))
      ALLOCATE(vospj(nstmax))
      ALLOCATE(vospt(nstmax))
      ALLOCATE(weight(nstmax))
      ALLOCATE(wosp(nstmax))
      ALLOCATE(wospg(nstmax))
      ALLOCATE(wospj(nstmax))
      ALLOCATE(wospt(nstmax))
      ALLOCATE(xbar(nstmax))
      ALLOCATE(xbarlg(nstmax))
c
c-----------------------------------------------------------------------
c
c     Allocate some additional arrays associated with the higher-order
c     (stiff) ODE integrator.
c
      nrct1 = nrct + 1
c
      ALLOCATE(ipivtr(nrct1))
c
      ALLOCATE(alphar(nrct1))
      ALLOCATE(betar(nrct1))
      ALLOCATE(delvcr(nrct1))
      ALLOCATE(rhsvcr(nrct1))
      ALLOCATE(dvjdte(nrct))
c
      ALLOCATE(armatr(nrct1,nrct1))
      ALLOCATE(grmatr(nrct1,nrct1))
c
      ALLOCATE(aimatr(kmax,kmax))
      ALLOCATE(mmmatr(nrct,kmax))
      ALLOCATE(sgmatr(nrct,kmax))
      ALLOCATE(xxmatr(kmax,nrct))
      ALLOCATE(xymatr(nrct,nrct1))
c
      ALLOCATE(cdacb(nbt,imchmx,2,nrct))
      ALLOCATE(ndacb(nbt,imchmx,2,nrct))
      ALLOCATE(ndactb(imchmx,2,nrct))
c
      ALLOCATE(hhcvec(nordmx))
      ALLOCATE(xhcvec(nordmx))
c
c-----------------------------------------------------------------------
c
c     Allocate some arrays needed to assist in writing an EQ6 pickup
c     file. These arrays do not appear on that file itself.
c
      ALLOCATE(uesr1(nctmax))
      ALLOCATE(cesr1(nctmax))
      ALLOCATE(ubsr1(nbt1mx))
      ALLOCATE(cbsr1(nbt1mx))
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set the format ("W" or "D") for the pickup and backup files.
c
      if (iopr(17) .le. 0) then
        upkfor = uinfor
      elseif (iopr(17) .eq. 1) then
        upkfor = 'W'
      elseif (iopr(17) .ge. 2) then
        upkfor = 'D'
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set the EQ6 calculational mode flag (.true. in EQ6).
c
      q6mode = .true.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set counters for warnings of problems associated with pressure
c     corrections.
c
c       nwndpc = number of warnings of no data to support pressure
c                  corrections
c       nweope = number of warnings of excursions outside the
c                  recommended pressure envelope
c
      nwndpc = 0
      nweope = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set counter for the number of points of reaction progress at
c     which warnings are issued regarding fugacities not being fixed at
c     specified values.
c
      nwnffg = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set counters which control the writing of warnings that little
c     solvent water remains in the equilibrium system.
c
      iwdh2o = 0
      nwdh2o = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set some constants that depend on the molecular weight of
c     water.
c
      omega = 1000./mwtsp(narn1)
      omeglg = log10(omega)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize a value for the density of the aqueous solution.
c
      qrho = .false.
      rho = 0.0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Reset the jflgi array so it contains the jflag values for
c     the basis species.
c
      do nb = 1,nbt
        ns = nbasp(nb)
        jflgi(nb) = jflag(ns)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set up the cdrw array. This provides a fast way to get the
c     reaction coefficient of H2O (Aqueous solution) in any reaction.
c
      call gcdrw(cdrs,cdrw,narn1,ndrs,ndrsmx,ndrsr,nst,nstmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set up the cdrtw array. This contains a stoichiometric sum of the
c     number of times H2O (Aqueous solution) is implied as a solvent
c     in any reaction.
c
      call gcdrtw(cdrs,cdrtw,narn1,narn2,ndrs,ndrsmx,ndrsr,
     $ nelect,no2gaq,nst,nstmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize concentration, activity coefficient, activity, etc.,
c     arrays.
c
      do ns = 1,nstmax
        acflg(ns) = 0.
        acflgo(ns) = 0.
        act(ns) = 0.
        conc(ns) = 0.
        mosp(ns) = 0.
        xbar(ns) = 0.
        affs(ns) = 0.
      enddo
c
      av = -99999.
      call initav(actlg,nstmax,av)
      call initav(conclg,nstmax,av)
      call initav(xbarlg,nstmax,av)
      call initav(loph,nptmax,av)
c
      nmax = nptmax
      call initaz(moph,nmax)
      call initaz(affp,nmax)
c
      nmax = nbtmax
      call initaz(amtb,nmax)
c
      nmax = 2*npetmx
      call initiz(ncmpe,nmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize (null) some other arrays.
c
      do ne = 1,netmax
        egexpa(ne) = 0.
        egexpc(ne) = 0.
      enddo
c
      nmax = jetmax*netmax
      call initaz(egexjc,nmax)
c
      nmax = ietmax*jetmax*netmax
      call initaz(egexs,nmax)
c
      nmax = ketmax*netmax
      call initaz(egexw,nmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize the mole fractions and activities of pure solids and
c     liquids.
c
      do ns = nlrn1,nlrn2
        act(ns) = 1.
        actlg(ns) = 0.
        xbar(ns) = 1.
        xbarlg(ns) = 0.
      enddo
c
      do ns = nmrn1,nmrn2
        act(ns) = 1.
        actlg(ns) = 0.
        xbar(ns) = 1.
        xbarlg(ns) = 0.
      enddo
c
      do ns = nfrn1,nfrn2
        act(ns) = 1.
        actlg(ns) = 0.
        xbar(ns) = 1.
        xbarlg(ns) = 0.
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize the flag array which marks if compositions are known
c     which maximize the affinities of the phases. Here assume that
c     this is the case for all phases except solid solutions.
c
      do np = 1,npt
        qxknph(np) = .true.
      enddo
      do np = ixrn1,ixrn2
        qxknph(np) = .false.
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (net .gt. 0) then
c
c       Initialize the numbers of moles of exchanger phases in the
c       moph and loph arrays.
c
        do nb = 1,nbt
          ns = nbaspd(nb)
          if (ns.ge.nern1 .and. ns.le.nern2) then
            np = nphasx(ns)
            mx = mtb(nb)
            lx = tlg(mx)
            moph(np) = mx
            loph(np) = lx
          endif
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize the jjsort, jssort, jcsort, and jgsort arrays by
c     setting each element equal to its index.
c
      call initii(jcsort,nst)
      call initii(jjsort,nst)
      call initii(jssort,nst)
      call initii(jgsort,ngt)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Expand the system description from the data read from the input
c     file. This includes estimating the numbers of moles of all phases
c     and species present, the concentrations, activity coefficients,
c     and activities of all the species, the ionic strength, and the
c     sum of the molalities of all aqueous aqueous solute species.
c
      call exivar(abar,acflg,acflgo,act,actlg,actwlc,adh,
     $ adhh,adhv,al10,aphi,azero,a3bar,a3bars,bdh,bdhh,bdhv,bdot,
     $ bdoth,bdotv,cco2,cdrs,cegexs,cgexj,conc,conclg,cpgexs,
     $ egexjc,egexjf,egexs,eps100,fje,fjeo,fo2,fo2lg,fsort,fugac,
     $ fugalg,fxi,fxio,ielam,iern1,iern2,ietmax,ifcphi1,ifcphi2,
     $ ifnnn,ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igstak,
     $ iindx1,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,
     $ ilrn2,ilzeta,imrn1,imrn2,insgf,iodb,iopg,ipbtmx,istack,
     $ ixrn1,ixrn2,izmax,jcsort,jern1,jern2,jetmax,jflag,jgext,
     $ jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jssort,jstack,
     $ kbt,kdim,kelect,kmax,km1,ko2gaq,kwater,kxt,loph,losp,
     $ lsort,mgext,moph,mosp,mrgexs,mtb,napmax,narn1,narn2,natmax,
     $ nazmmx,nazpmx,nbasp,nbaspd,nbt,nbtmax,nchlor,ncmpr,ndrs,
     $ ndrsmx,ndrsr,nelect,nern1,nern2,net,netmax,ngexsa,ngext,
     $ ngrn1,ngrn2,ngt,ngtmax,nhydr,nmutmx,nmxmax,nodbmx,nopgmx,
     $ noutpt,no2gaq,nphasx,npt,nptmax,nsltmx,nst,nstmax,nsxmax,
     $ nttyo,omega,omeglg,press,qhawep,qpit75,q6mode,sigmam,sigmmo,
     $ tempk,ugexj,ugexmo,uphase,uspec,xbar,xbarlg,xbarw,xbarwc,
     $ xbrwlc,xbrwlg,xlks,zchar,zchcu6,zchsq2,zgexj,zvclg1,zvec1)
c
c     Set the composition known flag for solid solutions present
c     in the ES.
c
      do np = ixrn1,ixrn2
        if (moph(np) .gt. 0.) qxknph(np) = .true.
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set up to make the initial equilibrium calculation.
c
c       kord   = the largest possible order on a given step
c       kstep  = step counter
c       nord   = the actual order used
c       npts   = the number of points available; restricts the
c                  value of kord
c
      kstep = 0
      kord = 0
      nord = 0
      npts = 1
      ncorr = 0
c
      km10 = km1
      kmt0 = kmt
      kx10 = kx1
      kxt0 = kxt
      kdim0 = kdim
c
      call copyia(iindx1,iindx0,kmax)
      call copyia(ipndx1,ipndx0,kmax)
      call copyca(uzvec1,uzvec0,kmax)
      call copyaa(zvclg1,zvclg0,kmax)
      call copyaa(zvec1,zvec0,kmax)
c
      qstart = .true.
      nordsv = 0
c
      xim1 = xistrt
      xi0 = xistrt
      xi1 = xistrt
      xistsv = xistrt
c
      xidump = 0.
c
      qriinf = .false.
c
      time0 = tistrt
      time1 = tistrt
      tistsv = tistrt
c
      delxi = 0.
      deltim = 0.
c
      rirec0 = 0.
      rirec1 = 0.
c
      do n = 1,nrctmx
        afrc0(n) = 0.
        afrc1(n) = 0.
        rreac0(n) = 0.
        rreac1(n) = 0.
        rrelr0(n) = 0.
        rrelr1(n) = 0.
      enddo
c
      call copyia(jreac,jreac0,nrctmx)
      call copyaa(morr,morr0,nrctmx)
      call copyaa(modr,modr0,nrctmx)
      call copyaa(mtb,mtb0,nbtmax)
      call copyaa(sfcar,sfcar0,nrctmx)
c
      call initaz(xirct,nrctmx)
      call initaz(xirct0,nrctmx)
c
      iexrt = 0
      jexrt = 0
      jscat = 0
      jscrt = 0
c
      betamx = 0.
      ubacmx = 'None'
      ubgamx = 'None'
c
      avkdim = real(kdim)
      prminf = 1./prcinf
      screwd = sscrew(5)
c
      tolsar = 1.e-8
      tolsrr = 1.e-8
c
      tiprnt = prcinf
      tiprnl = prcinf
      tiplot = prcinf
      tiplol = prcinf
c
      qbye = .false.
      qconst = .false.
      qdump = .false.
      qftpr2 = .false.
      qmod = .false.
      qrapch = .false.
      qriinf = .false.
      qtvchk = .false.
      qtrch = .false.
      qzprnt = .false.
      qstopx = .false.
      qshoot = .false.
      qvlsow = .false.
      qvhfxi = .false.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Zero various arrays pertaining to finite differences and
c     equivalent derivatives.
c
      do n = 1,nrd1mx
        dxsm00(n) = 0.
        dxsm10(n) = 0.
        dxsm11(n) = 0.
        fdrim1(n) = 0.
        fdri0(n) = 0.
        fdri1(n) = 0.
        drir0(n) = 0.
        drir0s(n) = 0.
      enddo
c
      nmax = nrd1mx*kmax
      call initaz(fdzv0,nmax)
      call initaz(fdzvm1,nmax)
      call initaz(dzvc0,nmax)
      call initaz(dzvc0s,nmax)
c
      do n = 1,nordmx
        fdawm1(n) = 0.
        fdaw0(n) = 0.
        daw0(n) = 0.
        fdehm1(n) = 0.
        fdeh0(n) = 0.
        deh0(n) = 0.
        fdo2m1(n) = 0.
        fdo20(n) = 0.
        do20(n) = 0.
        fdphm1(n) = 0.
        fdph0(n) = 0.
        dph0(n) = 0.
      enddo
c
      nmax = nordmx*nptmax
      call initaz(fdaf0,nmax)
      call initaz(fdafm1,nmax)
      call initaz(daffp0,nmax)
c
      nmax = nordmx*npetmx
      call initaz(fdpe0,nmax)
      call initaz(fdpem1,nmax)
      call initaz(demop0,nmax)
c
      nmax = 2*npetmx
      call initiz(ncmpe0,nmax)
c
      nmax = nordmx*nsetmx
      call initaz(fdse0,nmax)
      call initaz(fdsem1,nmax)
      call initaz(demos0,nmax)
c
      nmax = nordmx*nrctmx
      call initaz(fdarm1,nmax)
      call initaz(fdar0,nmax)
      call initaz(dafrc0,nmax)
      call initaz(fdrem1,nmax)
      call initaz(fdre0,nmax)
      call initaz(fdre1,nmax)
c
      nmax = nrd1mx*nrctmx
      call initaz(fdrrm1,nmax)
      call initaz(fdrr0,nmax)
      call initaz(fdrr1,nmax)
      call initaz(drer0,nmax)
      call initaz(drer0s,nmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set up the ixbasp and cjbasp arrays. The former is a flag
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
c     Make the equilibrium calculation at the initial point of reaction
c     progress.
c
      call eqshel(aadh,aadhh,aadhv,aamatr,aaphi,abar,abdh,
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
c     Note on ier codes returned by EQ6/eqshel.f:
c
c       =    0  Okay
c       =   10  Go back and take a smaller step size to avoid exceeding
c                 the supersaturation tolerance (tolsst)
c       =  170  Too much of a phase was destroyed under the flow-through
c                 open system model; go back and first move part of the
c                 mass of protected phases in the ES to the PRS
c       =  180  One of a number of problems occurred which may be
c                 resolvable, at least partially, by going back and
c                 cutting the step size
c       =  190  Need to slide over a region of computational
c                 instability, but sliding is inhibited; go back, but
c                 terminate work on the current problem
c
c     Note: EQ6/eqshel.f can't return ier = 10 or 170 when called at
c     the initial point of reaction progress.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ier .eq. 180) then
c
c       Normally this error code results in a reduction in the step
c       size. This can't be done at the initial point of reaction
c       progress.
c
        write (noutpt,1010)
        write (nttyo,1010)
 1010   format(/' * Error- (EQ6/path) The equilibrium calculation',
     $  ' failed at the initial value',/7x,'of reaction progress.',
     $   " Can't cut the step size to try to recover.",/7x,'See',
     $   ' previous notes and warnings for more information.')
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ier .gt. 0) then
c
c       Here ier should have a value of 190.
c
        write (noutpt,1020)
        write (nttyo,1020)
 1020   format(/' * Error - (EQ6/path) The equilibrium calculation',
     $  ' failed at the initial value',/7x,'of reaction progress.')
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the activity of water.
c
      actwlg = actlg(narn1)
      actw = texp(actwlg)
      awstrt = actw
      aw0 = actw
      aw1 = actw
c
c     Get the weights (masses) of solvent, total dissolved solutes,
c     and aqueous solution, and get the aqeuous solution density.
c
      call gwdenp(adwipp,bdwipp,jcsort,mlmrra,mosp,mrmlra,
     $ mwtsp,narn1,narn2,nstmax,qdwipp,rhoc,rhowc,tdsgks,tdsglw,
     $ tdspkc,tdsplc,tempc,vosol,wfh2o,wftds,wkgwi,woh2o,
     $ wosol,wotds)
c
      if (qdwipp) then
        qrho = .true.
        rho = rhoc
      else
        qrho = .false.
        rho = 0.0
      endif
c
c     Compute pH, Eh, and pe-, all with reference to appropriate
c     pH scales. Also compute the pHCl.
c
      call gpheh(acflg,actlg,actwlg,adh,ah,ahmes,ahnbs,conc,
     $ eh,ehfac,ehmes,ehnbs,farad,fo2lg,fxi,iopg,mrmlra,nchlor,nhydr,
     $ nopgmx,noutpt,nstmax,nttyo,pch,pe,pemes,penbs,ph,phcl,phmes,
     $ phnbs,qphcl,qredox,qrho,xlke)
c
      wkgh2o = 1.e-3*woh2o
      wkgsol = 1.e-3*wosol
      wwstrt = wkgh2o
      phstrt = ph
      ph0 = ph
      ph1 = ph
      ehstrt = eh
      eh0 = eh
      eh1 = eh
      o2strt = fo2lg
      fo2lg0 = fo2lg
      fo2lg1 = fo2lg
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Do automatic basis switching. Optimize the basis set so that
c     the basis species for each mass balance tends to be the species
c     which dominates that mass balance.
c
      kstpab = 0
      if (iopt(12) .gt. 0) then
        call absswb(adhfs,adhfsx,advfs,advfsx,avcnst,axhfs,
     $  axhfsx,axlks,axlksx,axvfs,axvfsx,cdrs,cdrtw,cdrsx,cdrw,
     $  csts,dhfs,dvfs,eps100,ibswx,iindx1,iodb,ipch,ipchmx,ipcv,
     $  ipcvmx,jcsort,jflag,jsflag,kbt,kmax,mosp,mtb,narn1,narn2,
     $  narxmx,narxt,nbasp,nbaspd,nbaspx,nbt,nbtmax,nbw,ndrs,ndrsmx,
     $  ndrsr,ndrsrx,ndrsx,nelect,nhydr,nodbmx,no2gaq,noutpt,nst,
     $  nstmax,nsts,nstsmx,nstsr,nswtch,ntpr,ntprmx,nttyo,presg,
     $  press,qbassw,qbswx,tempc,uspec,uzvec1,weight,xhfs,xvfs,xlks)
c
        if (nswtch .gt. 0) then
c
c         It is not necessary here to set the qabswx flag, which
c         indicate that automatic basis switching has just been done.
c
          write (noutpt,1030) nswtch
 1030     format(/' ',i2,' basis switches were executed automatically',
     $    ' after solving',/'at the initial value of reaction',
     $    ' progress.')
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
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1100) kstep,iter
 1100 format(' Steps completed= ',i5,', iter= ',i3)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Save some data at the initial point, as it is after the
c     calculations at this point are complete. Note that the set of
c     matrix variables may have been changed due to mineral
c     precipitation, etc.
c
      call copyia(iindx1,iindx0,kdim)
      call copyia(ipndx1,ipndx0,kdim)
c
      km10 = km1
      kmt0 = kmt
      kx10 = kx1
      kxt0 = kxt
      kdim0 = kdim
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check gases for which the fugacity is supposed to be fixed.
c
      do n = 1,nffg
        ng = jffg(n)
        xlf = xlkffg(n)
        xf = texp(xlf)
        dx = fugalg(ng) - xlkffg(n)
        if (abs(dx) .gt. toldl) then
          j2 = ilnobl(uffg(n))
          write (noutpt,1150) uffg(n)(1:j2),fugac(ng),xf
          write (nttyo,1150) uffg(n)(1:j2),fugac(ng),xf
 1150     format(/' * Warning - (EQ6/eq6) The fugacity of ',a,' is',
     $    /7x,1pg12.5,' bars at the start of the run. It is supposed',
     $    /7x,'to be fixed at ',e12.5,' bars. There is an insufficient',
     $    /7x,'mass of the gas component present in the system to',
     $    /7x,'"saturate" it at the desired fugacity. You may wish to',
     $    /7x,'restart this run, adding such mass using the moffg',
     $    /7x,'parameter on the input file. Add only about 0.5 to 1.0',
     $    /7x,'mole at a time.')
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nrct .gt. 0) then
c
c       Calculate the initial affinities and rates of the irreversible
c       reactions.
c
        call raff(acflg,actlg,afcnst,affp,afrc1,bpx,cdrs,cgexj,
     $  ibpxmx,ibpxt,iern1,ietmax,iktmax,ixrn1,ixrn2,jcode,jern1,jern2,
     $  jetmax,jflag,jgext,jpflag,jsflag,jsol,ncmpr,ndrs,ndrsmx,ndrsr,
     $  nertmx,net,netmax,ngext,noutpt,nptmax,nrct,nrctmx,nrndex,
     $  nstmax,nttyo,nxridx,nxrtmx,nxtmax,rxbar,uphase,uspec,wfac,
     $  xbar,xbarlg,xgers,xlks)
c
c       For each reactant, check the status flag against the affinity
c       and the number of moles present to ensure that any reactant
c       tagged as available to react, having an affinity favoring
c       dissolution, and having no moles remaining is retagged as
c       exhausted.
c
        do nrc = 1,nrct
          if (afrc1(nrc) .gt. -tolsar) then
            if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
              if (morr(nrc) .le. 0.) then
c
c               The following exception is intended to facilitate the
c               simulation of evaporation by declaring H2O as a
c               reactant with a negative relative rate and no initial
c               moles remaining. However, any reactant may be treated
c               analogously if desired.
c
                if (.not.(nrk(1,nrc).eq.1 .and. rk(1,1,nrc).lt.0.)) then
                  jreac(nrc) = 1
                endif
              endif
            endif
          endif
        enddo
c
c       Check the reactants for saturation.
c
        call rsatch(csts,egers,egexs,iern1,ietmax,iindx1,iktmax,
     $  iopt,ipndx1,jcode,jern1,jern2,jetmax,jgext,jpflag,jreac,kmax,
     $  km1,kmt,kx1,kxt,loph,losp,moph,morr,mosp,mrgers,mtb,mtb0,
     $  nbaspd,nbtmax,ncmpr,nern1,nern2,nert,nertmx,netmax,ngext,
     $  noptmx,noutpt,nptmax,nrct,nrctmx,nrk,nrndex,nstmax,nsts,nstsmx,
     $  nstsr,nttyo,nxridx,nxrt,nxrtmx,qreq,rxbar,tolxsf,uphase,ureac,
     $  uspec,xbar,xbarlg,zvclg1,zvec1)
c
c       Calculate rates at the initial point by evaluating the
c       rate laws.
c
        call rtcalc(act,afrc1,cdac,csigma,eps100,fkrc,idirec,
     $  imchmx,imech,iodb,iopt,jcode,jreac,morr,morr0,mwtrc,ndac,
     $  ndact,ndctmx,nodbmx,noptmx,nord,noutpt,nrk,nrct,nrctmx,nsk,
     $  nstmax,nttyo,prcinf,prminf,qriinf,rirec1,rk,rreac1,rrelr1,
     $  rtcnst,rrxfi1,sfcar,sfcar0,ssfcar,udac,ureac)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the total affinity (aft1).
c
      call caft1(afrc1,aft1,nrct,nrctmx,rrelr1)
      aft0 = aft1
      aftm1 = aft1
      qaft1 = .false.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute apparent "whole-phase" equivalent fractions and mole
c     fractions of the exchange ions present in generic ion exchanger
c     phases. Cations and anions are treated separately in these
c     calculations.
c
      call gegexw(cegexs,egexpc,egexpa,egexw,iern1,iern2,
     $ ietmax,jern1,jetmax,jgext,kern1,kern2,ketmax,kgexsa,moph,
     $ mosp,netmax,ngexsa,ngext,noutpt,nptmax,nstmax,nttyo,
     $ xgexw,zchar)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize arrays associated with finite-difference description
c     of the number of mole of phases and species in the Equilibrium
c     System (ES). Here the iemop and emop arrays respectively contain
c     the indices and numbers of moles of phases in the ES. The iemos
c     and emos array are the analogs for the species of these phases
c     (but for the aqueous solution phase, only the species H2O(l)
c     is tracked by this mechanism). The ncmpe array is a species range
c     pointer array for the phases, analogous to ncmpr. The fdpe0 and
c     fdse0 arrays contain the finite differences for the phases and
c     species, respectively. The demop and demos arrays contain the
c     corresponding derivatives.
c
c     Initialize the index arrays.
c
      call iiemop(iemop,iemos,iindx1,ipndx1,jsflag,kdim,kmax,
     $ ncmpe,ncmpr,noutpt,npet,npetmx,npt,nptmax,nset,nsetmx,nstmax,
     $ nttyo,uaqsln,uspec,uphase)
c
c     Load the corresponding numbers of moles.
c
      do npe = 1,npet
        np = iemop(npe)
        emop(npe) = moph(np)
        nr1 = ncmpe(1,npe)
        nr2 = ncmpe(2,npe)
        do nse = nr1,nr2
          ns = iemos(nse)
          emos(nse) = mosp(ns)
        enddo
      enddo
c
      npet0 = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate various secondary data at the initial point.
c
      call cdappl(acflg,acfw,acfwlg,actlg,actw,actwlg,adwipp,
     $ afcnst,affpd,affsd,ah,ahrc,alk,alk1,alk2,alki,atwt,bdwipp,
     $ cdrsd,cess,conc,csts,ctb,cteaq,dvoso,dwoso,eh,ehfac,ehrc,
     $ eps100,farad,fdpe0,fdse0,fjest,fo2lg,fo2lrc,fxist,iaqsln,
     $ iemop0,iemos0,iern1,iern2,ifrn1,ifrn2,ilrn1,ilrn2,imrn1,
     $ imrn2,iopt,ixrn1,ixrn2,jcode,jcsort,jern1,jflag,jflagd,
     $ jgext,jpflag,jsflag,jssort,modr,moph,mophg,mophj,mopht,morr,
     $ mosp,mospg,mospj,mospt,mprph,mprsp,mrgers,mrmlra,mtb,
     $ mtbaq,mte,mteaq,mwtges,mwtrc,mwtsp,narn1,narn2,nat,nbasp,
     $ nbaspd,nbt,nchlor,ncmpe0,ncmpr,nct,ndrsd,ndrsrd,nelect,
     $ nern1,nern2,nert,ness,nessr,net,nfrn1,nfrn2,ngext,ngrn1,
     $ ngrn2,ngt,nhydr,nlrn1,nlrn2,nlt,nmrn1,nmrn2,nmrt,nmt,no2gaq,
     $ npchk,npet,npet0,npt,npts,nrct,nrndex,nst,nsts,nstsr,ntf1,
     $ ntf1t,ntf2,ntf2t,nxridx,nxrn1,nxrn2,nxrt,nxt,osc,oscst,
     $ omega,pe,perc,ph,phmes,ppmwb,ppmwe,qriinf,qxknph,rreacn,
     $ rreac1,rxbar,sfcar,sidrsp,sidrph,sigmam,sigmst,tdays,tempc,
     $ tf1,tf2,thours,time1,tmins,tyears,uphase,uspec,vodrt,voph,
     $ vophg,vophj,vopht,vosoct,vosp,vospg,vospj,vospt,vosp0,
     $ vreac,wfh2o,wkgwi,wodr,wodrt,woph,wophg,wophj,wopht,
     $ worr,worrt,wosoct,wosp,wospg,wospj,wospt,xbar,xlke,
     $ xlksd,zchcu6,zchsq2)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print results at the initial point.
c
      call scripz(abar,acflg,acfw,acfwlg,actlg,actw,actwlg,
     $ affpd,affsd,afrc1,aft1,ah,ahmes,ahnbs,ahrc,alki,alk1,
     $ alk2,awmax,awmin,a3bar,cbsr,cdrsd,cegexs,cesr,conc,conclg,csts,
     $ ctb,cteaq,dvoso,dwoso,egers,egexjc,egexjf,egexpa,egexpc,egexs,
     $ egexw,eh,ehmax,ehmes,ehmin,ehnbs,ehrc,elecsr,electr,fje,fjest,
     $ fo2,fo2lg,fo2lrc,fugac,fugalg,fxi,fxist,iaqsln,iemop,iemop0,
     $ iemos,iemos0,iern1,iern2,iexr,iexrt,ifrn1,ifrn2,ilrn1,ilrn2,
     $ imech,imrn1,imrn2,iopg,iopr,iopt,ipndx1,ixrn1,ixrn2,jcode,jcsort,
     $ jern1,jern2,jexr,jexrt,jflag,jflagd,jflgi,jgext,jgsort,jpflag,
     $ jreac,jsca,jscat,jscr,jscrt,jsflag,jsol,jssort,kbt,kern1,kern2,
     $ kgexsa,km1,kmt,kx1,kxt,kstep,kstpmx,loph,losp,mlmrra,modr,
     $ moph,mophg,mophj,mopht,morr,mosp,mospg,mospj,mospt,mprph,
     $ mprsp,mrgers,mrmlra,mwtrc,mwtsp,narn1,narn2,nat,nbasp,nbaspd,
     $ nbt,ncmpe,ncmpe0,ncmpr,nct,ndrsd,ndrsrd,nelect,nern1,nern2,
     $ nert,net,nfrn1,nfrn2,ngext,ngexsa,ngrn1,ngrn2,ngt,nhydr,nhydx,
     $ nlrn1,nlrn2,nlt,nmrn1,nmrn2,nmrt,nmt,noutpt,no2gaq,npet,npet0,
     $ npt,npts,nrct,nrdxsp,nrk,nrndex,nst,nsts,nstsr,ntf1t,ntf2t,
     $ nxridx,nxrn1,nxrn2,nxrt,nxt,osc,oscst,omega,o2max,o2min,pch,pe,
     $ pemes,penbs,perc,ph,phcl,phmax,phmes,phmin,phnbs,ppmwe,presg,
     $ press,qaft1,qftpr2,qmod,qphcl,qredox,qrho,qriinf,qstopx,qvhfxi,
     $ qvlsow,qzprnt,rho,rhoc,rhowc,rk,rreacn,rreac1,rrelr1,rxbar,
     $ sfcar,sidrph,sidrsp,sigmst,sigmam,ssfcar,tdays,tdsglw,tdspkc,
     $ tdsplc,tempc,thours,time1,timemx,tmins,tolsat,tolxsf,tolxst,
     $ tolxsu,tyears,uelem,ugermo,ugexj,ugexmo,uphase,ureac,uspec,
     $ uxtype,vodrt,voph,vophg,vophj,vopht,vosoct,vosol,vosp,vospg,
     $ vospj,vospt,vreac,wfh2o,wftds,wkgwi,woh2o,wodr,wodrt,woph,
     $ wophg,wophj,wopht,worr,worrt,wosoct,wosol,wosp,wospg,wospj,
     $ wospt,wotds,xbar,xbarlg,xbarw,xbrwlg,xgers,xgexw,xi1,xidump,
     $ ximax,xistsv,xirct,zchar)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write results at the initial point on the tabx file. This file
c     is used to create the plot file tab.
c
      if (iopt(18) .ge. 0) then
        if (qtatxt) then
c
c         The TAB file is an ordinary text file.
c
          call wrtabx(actlg,afrc1,aft1,alk,cteaq,dvoso,dwoso,
     $    eh,fo2lg,iindx1,iktmax,iopt,ipndx1,kmax,km1,kmt,kstep,kx1,
     $    kxt,loph,ncmpr,modr,mopht,narn1,mosp,nct,nctmax,noptmx,
     $    nptmax,nrct,nrctmx,nstmax,ntabx,ntidmx,ntitl2,ntitld,ntitmx,
     $    nxtmax,pe,ph,ppmwe,prcinf,press,prminf,qbye,qmod,qriinf,
     $    tempc,time1,uelem,uphase,uplatm,ureac,uspec,usteq6,utitl2,
     $    utitld,uveeq6,vodrt,vosoct,wodrt,woh2o,wosoct,xbar,xi1)
        else
c
c         The TAB file is a .csv file.
c
          call wrtabc(acflg,actlg,actw,afrc1,aft1,alk,conclg,
     $    cteaq,ctb,dvoso,dwoso,eh,fje,fo2lg,fugac,fxi,iktmax,iopt,
     $    jflag,jsflag,kmax,kstep,kx1,kxt,mrmlra,modr,mosp,mospt,
     $    moph,mopht,narn1,narn2,nbasp,nbaspd,nbt,nbtmax,ncmpr,nct,
     $    nctmax,nelect,ngrn1,ngrn2,ngtmax,nhydr,nhydx,nllnmx,
     $    no2gaq,noptmx,noutpt,npt,nptmax,nrct,nrctmx,nstmax,ntabx,
     $    ntidmx,ntitl2,ntitld,ntitmx,nttyo,nxrn1,nxrn2,nxtmax,pe,ph,
     $    phmes,ppmwb,ppmwe,prcinf,press,prminf,qrho,qriinf,rho,rhowc,
     $    sidrph,sigmam,tdsgks,tdsglw,tempc,time1,uelem,ulinex,
     $    uphase,uplatm,ureac,uspec,usteq6,utitl2,utitld,uveeq6,
     $    vodrt,vosoct,wkgh2o,wodrt,wosoct,xbar,xbarlg,xi1)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (nttyo,1200)
 1200 format(1x)
c
c     Write entertainment for the user, showing the progress
c     of the current run.
c
      call wrentu(actw,eh,fo2lg,iopg,iopt,kstep,nopgmx,noptmx,
     $ nttyo,ph,qredox,time1,xi1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(6) .gt. 0) then
c
c       Clear ES solids at the starting point. Fictive fugacity-fixing
c       minerals are not cleared.
c
        write (noutpt,1230)
        write (nttyo,1230)
 1230   format(/' * Note (EQ6/path) Clearing solid phases from the',
     $  ' equilibrium system (ES)',/7x,'at the starting point of',
     $  ' reaction progress.')
c
        call clress(csts,iindx1,ipndx1,jpflag,jsflag,kdim,kmax,
     $  km1,kmt,kx1,kxt,loph,losp,moph,mosp,mtb,mtbaq,nbt,nbtmax,
     $  nptmax,nstmax,nsts,nstsmx,nstsr,ufixf,uzvec1,zvec1,zvclg1)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make sure that there is something to define a reaction path
c     (reactants, changing temperature, changing pressure).
c
      if (nrct.le.0 .and. qcntmp .and. qcnpre) then
        write (noutpt,1250)
        write (nttyo,1250)
 1250   format(/' * Note - (EQ6/path) No reaction path has been',
     $  ' defined by on the input file.',/7x,'There are no specified',
     $  ' reactants (irreversible reactions), and',/7x,'the',
     $  ' temperature and pressure are fixed. No steps will be taken.')
        go to 990
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check stop conditions.
c
      call chkstc(actw,awmax,awmin,eh,ehmax,ehmin,fo2lg,iopt,
     $ jreac,kstep,kstpmx,noptmx,noutpt,nrct,nrctmx,nttyo,o2max,o2min,
     $ ph,phmax,phmin,prcinf,qaft1,qcnpre,qcntmp,qconst,qredox,qstop,
     $ qvhfxi,qvlsow,timemx,time1,tolxst,tolxsu,ximax,xi1)
c
      if (qstop) go to 990
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      nbkupn = nbkupa
c
      if (iopt(16) .ge. 0) then
c
c       Prepare to write results at the initial point to the
c       backup file.
c
        call setpk6(actwlg,awmax,awmaxi,awmin,awmini,eh,ehmax,
     $  ehmaxi,ehmin,ehmini,fo2lg,iindx1,jflag,jflgi,kbt,kdim,kmax,
     $  kprs,mprph,mprphi,mprsp,mprspi,mtb,mtbi,mtbaq,mtbaqi,nbasp,
     $  nbaspd,nbaspi,nbti,nbtmax,ncmpr,nobswt,noutpt,nprpmx,nprpti,
     $  nprsmx,nprsti,npt,nptmax,nttyo,nstmax,o2max,o2maxi,o2min,
     $  o2mini,ph,phmax,phmaxi,phmin,phmini,prcinf,press,pressi,
     $  tempc,tempci,time1,timemx,timmxi,tistti,ubmtbi,uobsw,uphase,
     $  uprphi,uprspi,uspec,uzveci,uzvec1,xi1,ximax,ximaxi,xistti,
     $  zvclgi,zvclg1)
c
c       Rewind the current backup file (BAKUPA or BAKUPB).
c
        rewind (nbkupn)
c
c       Write results at the initial point to the backup file.
c
        if (upkfor(1:1) .eq. 'W') then
c
c         Compact (W) format.
c
c         Calling sequence substitutions:
c           nbkupn for newin
c
          call wr6pkw(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,
     $    dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,
     $    dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,
     $    dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,iact,
     $    ibsrti,iesrti,ietmax,igerti,iktmax,imchmx,imech,iodb,iopg,
     $    iopr,iopt,itermx,ixrti,jcode,jetmax,jflgi,jgerti,jgext,
     $    jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,
     $    kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,
     $    mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,
     $    nertmx,net,netmax,nbkupn,nffg,nffgmx,ngexrt,nobswt,nodbmx,
     $    nopgmx,noprmx,noptmx,nordmx,nprpmx,nprpti,nprsmx,nprsti,
     $    nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,
     $    ntitmx,ntrymx,nttkmx,nxmdmx,nxmod,nxopex,nxopmx,nxopt,
     $    nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,pressb,
     $    pressi,ptk,qgexsh,rkb,rxbari,sfcar,ssfcar,tempcb,tempci,
     $    tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,trkb,ttk,
     $    ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,ugersi,
     $    ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,ureac,
     $    usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,uxopt,
     $    uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,
     $    xlkmod,xvfgex,zgexj,zvclgi)
        else
c
c         Menu-style (D) format.
c
c         Calling sequence substitutions:
c           nbkupn for newin
c
          call wr6pkd(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,
     $    dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,
     $    dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,
     $    dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,iact,
     $    ibsrti,iesrti,ietmax,igerti,iktmax,imchmx,imech,iodb,iopg,
     $    iopr,iopt,itermx,ixrti,jcode,jetmax,jflgi,jgerti,jgext,
     $    jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,
     $    kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,
     $    mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,
     $    nertmx,net,netmax,nbkupn,nffg,nffgmx,ngexrt,nobswt,nodbmx,
     $    nopgmx,noprmx,noptmx,nordmx,nprpmx,nprpti,nprsmx,nprsti,
     $    nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,
     $    ntitmx,ntrymx,nttkmx,nxmdmx,nxmod,nxopex,nxopmx,nxopt,
     $    nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,pressb,
     $    pressi,ptk,qgexsh,rkb,rxbari,sfcar,ssfcar,tempcb,tempci,
     $    tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,trkb,ttk,
     $    ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,ugersi,
     $    ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,ureac,
     $    usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,uxopt,
     $    uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,
     $    xlkmod,xvfgex,zgexj,zvclgi)
        endif
c
        if (iopt(16) .eq. 0) then
c
c         Switch to write on the other backup file the next time.
c
          nbkupn = nbkupb
        endif
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize logical flags used in computing the reaction path.
c
c       qzdump = .true. if the last value of Xi was a point
c                  corresponding to the PRS transfer interval.
c       qdump  = .true. if a PRS transfer has just taken place.
c
      qzdump = .false.
      qzplot = .false.
      qskip = .false.
      qreq = .false.
      qreax = .false.
      qrpcfl = .false.
      qmod1 = .false.
      qmod2 = .false.
      qdmpr1 = .false.
      qdmpr2 = .false.
c
      qprntx = .false.
      qprntt = .false.
      qprnlx = .false.
      qprnlt = .false.
c
      qplotx = .false.
      qplott = .false.
      qplolx = .false.
      qplolt = .false.
c
      qprph0 = .false.
      qprph1 = .false.
      qpreh0 = .false.
      qpreh1 = .false.
      qpro20 = .false.
      qpro21 = .false.
      qpraw0 = .false.
      qpraw1 = .false.
c
      qplph0 = .false.
      qplph1 = .false.
      qpleh0 = .false.
      qpleh1 = .false.
      qplo20 = .false.
      qplo21 = .false.
      qplaw0 = .false.
      qplaw1 = .false.
c
      qhcon = .false.
      if (iopt(14) .eq. 2) qhcon = .true.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set or initialize other variables used in computing the
c     reaction path.
c
      kstpze = 0
      kstpmn = 0
      lprcin = tlg(prcinf)
c
      dlxtmx = prcinf
      xilim = prcinf
c
      scalim = 1.0
      scnsti = texp(sscrew(6)) - 1.
      scnstd = texp(-sscrew(6)) - 1.
c
      i = irang - 8
      fdx = real(i)
      xx = 100.
      fdx = min(fdx,xx)
      fdlim = texp(fdx)
      if (iodb(1) .ge. 1) then
        write (noutpt,1300) fdlim
        write (nttyo,1300) fdlim
 1300   format(/' * Note - (EQ6/path) The magnitude of the finite',
     $  ' difference functions',/7x,'is bounded by ',1pe10.3,'.',/)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Re-set the arrays associated with finite-difference description
c     of the number of mole of phases and species in the Equilibrium
c     System (ES).
c
c     Re-set the index arrays.
c
      call iiemop(iemop,iemos,iindx1,ipndx1,jsflag,kdim,kmax,
     $ ncmpe,ncmpr,noutpt,npet,npetmx,npt,nptmax,nset,nsetmx,nstmax,
     $ nttyo,uaqsln,uspec,uphase)
c
c     Re-load the corresponding numbers of moles.
c
      do npe = 1,npet
        np = iemop(npe)
        emop(npe) = moph(np)
        nr1 = ncmpe(1,npe)
        nr2 = ncmpe(2,npe)
        do nse = nr1,nr2
          ns = iemos(nse)
          emos(nse) = mosp(ns)
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize the akmat0 (F * B) and akmat1 (F * C) matrices. Each
c     relates finite differences to the equivalent derivatives in a
c     truncated Taylor's series. The akmat0 matrix is part of the
c     predictor function. the akmat1 matrix is part of the corrector.
c
      nmax = nrd1mx*nrd1mx
      call initaz(akmat0,nmax)
      call initaz(akmat1,nmax)
c
      do i = 1,nrd1mx
        akmat0(i,i) = fctrl(i)
        akmat1(i,i) = fctrl(i)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set remaining step counters and order parameters.
c
c       kordlm = the largest allowed value of kord (largest order
c                  considered for a given step), hence the largest
c                  allowed value of nord (the actual order used
c                  for a given step)
c       ndelay = the number of points that must be done before
c                  allowing the order to build up from zero
c       kly    = a counter used in restricting the growth of
c                   the step size; if it is >= 0, the scale limit
c                   (scalim) is set to a smaller value than it would be
c                   otherwise
c       kstppr = counter of steps since the last print point
c       kstppl = counter of steps since the last plot point
c
      kordlm = nordmx
      if (qecon) kordlm = min(2,nordmx)
      if (qscon) kordlm = 0
c
      ndelay = 1
c
      kstppr = 0
      kstppl = 0
c
      kly = 0
      delxi = dlxmx0
c
      naft1 = 0
      kaft1 = 0
c
      nsawth = 0
      jsawth = 0
      qsawth = .false.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the first print, plot and dump points.
c
      call sippdp(actw,aw0plo,aw0prn,aw1plo,aw1prn,dlaplo,
     $ dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,
     $ dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,dlxprl,
     $ dlxprn,eh,eh0plo,eh0prn,eh1plo,eh1prn,eps100,fo2lg,lprcin,
     $ o20plo,o20prn,o21plo,o21prn,ph,ph0plo,ph0prn,ph1plo,ph1prn,
     $ prcinf,qredox,tiplol,tiplot,tiprnl,tiprnt,tistsv,xidump,
     $ xiplol,xiplot,xiprnl,xiprnt,xistsv)
c
       qtplo = dltplo.lt.prcinf .or. dltpll.lt.prcinf
       qtprn = dltprn.lt.prcinf .or. dltprl.lt.prcinf
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make a new step in reaction progress.
c
  100 kstep = kstep + 1
      kstppr = kstppr + 1
      kstppl = kstppl + 1
      kstpab = kstpab + 1
      avkdim = (avkdim*kstep + real(kdim))/real(kstep + 1)
c
      kly = kly - 1
      jsawth = jsawth - 1
      qmod2 = qmod1
      qmod1 = qmod .or. qreax
      qrapch = .false.
c
      aftm1 = aft0
      aft0 = aft1
      ncorr = 0
      dlxode = prcinf
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for instability (sawtoothing) in the ODE integration.
c
      if (iopt(2) .gt. 0) then
        rx0 = fdri0(1)
        rx1 = rirec1 - rirec0
        rxx = rx0*rx1
        if (rxx .lt. 0.) then
          nsawth = nsawth + 1
          if (nsawth .ge. 8) then
c
c           Have sawtoothing. Take special action the next jsawth steps.
c
            jsawth = 20
            write (noutpt,'(" Sawtoothing detected")')
            write (nttyo,'(" Sawtoothing detected")')
          endif
        else
          nsawth = 0
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Pick the maximum allowed order for the current step (jordlm).
c
      if (jsawth .le. 0) then
        jordlm = kordlm
      else
        jordlm = min(jordlm,2)
      endif
      jordlm = min(kord,jordlm)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Save information at the current point of reaction progress
c     and update the finite differences for variables that are
c     tracked using such.
c
      call stepfd(acflg,acflg0,affp0,affp,afrc0,afrc1,aw0,aw1,
     $ delxi,dxsm00,eh0,eh1,emop,emop0,emos,emos0,fdafm1,fdaf0,
     $ fdarm1,fdar0,fdawm1,fdaw0,fdehm1,fdeh0,fdlim,fdo2m1,fdo20,
     $ fdpem1,fdpe0,fdphm1,fdph0,fdrem1,fdre0,fdrim1,fdri0,fdrrm1,
     $ fdrr0,fdsem1,fdse0,fdzvm1,fdzv0,fje,fje0,fo2lg0,fo2lg1,fxi,
     $ fxi0,iemop,iemop0,iemos,iemos0,iindx0,iindx1,iodb,iopt,ipndx0,
     $ ipndx1,jcode,jreac,jreac0,jpflag,kdim,kdim0,kmax,km1,km10,kmt,
     $ kmt0,kord,kx1,kx10,kxt,kxt0,modr,modr0,moph,moph0,morr,morr0,
     $ mosp,mosp0,mtb,mtb0,nbt,nbtmax,ncmpe,ncmpe0,nodbmx,noptmx,
     $ nordmx,noutpt,npet,npetmx,npet0,npt,nptmax,npts,nrct,nrctmx,
     $ nrd1mx,nrndex,nset,nsetmx,nset0,nstmax,ph0,ph1,qredox,rirec0,
     $ rirec1,rreac0,rreac1,rrelr0,rrelr1,sfcar,sfcar0,sigmam,sigmm0,
     $ tempc,tempc0,time0,time1,uzvec0,uzvec1,xi0,xi1,xim1,xirct,
     $ xirct0,zvclg0,zvclg1,zvec0,zvec1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(1) .eq. 2) then
c
c       Doing the fluid-centered flow-through open system model.
c
        if (qdump) then
c
c         Make a general partial shift of mass of the minerals in
c         the ES to the PRS. This is part of the fluid-centered
c         flow-through open system model. A partial shift of a phase
c         leaves some of its mass in the ES.
c
          call pshfta(csts,emop,emop0,emos,emos0,fdpe0,fdpem1,
     $    fdse0,fdsem1,iemop,iemos,iern1,iern2,ietmax,iindx1,imrn1,
     $    imrn2,iodb,ipndx1,ixrn1,ixrn2,jcsort,jern1,jetmax,jgext,
     $    jpflag,jsflag,kbt,km1,kmax,kmt,kx1,kxt,loph,losp,moph,mosp,
     $    mprph,mprsp,mrgexs,mtb,mtb0,nbasp,nbaspd,nbt,nbtmax,ncmpe,
     $    ncmpr,netmax,ngext,nodbmx,nordmx,noutpt,npet,npetmx,npt,
     $    nptmax,nsetmx,nstmax,nsts,nstsmx,nstsr,nttyo,uaqsln,ufixf,
     $    uphase,uspec,xbar,xbarlg,zklgmn,zklogl,zvclg0,zvclg1,
     $    zvec0,zvec1)
c
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(1) .eq. 2) then
c
c       Doing the fluid-centered flow-through open system model.
c
        qbye = .false.
        if (qmod2 .or. kstep.eq.2) then
c
c         On the previous step, a reactant became saturated or
c         exhausted, a new phase was added to the ES, or the previous
c         step was the first step of the run.
c
          if (kmt.ge.km1 .or. kxt.ge.kx1 .or. net.gt.0) then
c
c           Make a total PRS shift of any mineral which is disappearing
c           and purge it from the ES. The flag variable qbye is set
c           to true if any such shift is made here.
c
            call dumpdp(csts,demop0,demos0,emop,emop0,emos,
     $      emos0,fdpe0,fdpem1,fdse0,fdsem1,iemop,iemos,iern1,iern2,
     $      ietmax,iindx0,iindx1,imrn1,imrn2,iodb,ipndx0,ipndx1,ixrn1,
     $      ixrn2,jcsort,jern1,jetmax,jgext,jpflag,jsflag,kbt,kdim,
     $      kdim0,kmax,km1,km10,kmt,kmt0,kord,kstep,kx1,kx10,kxt,
     $      kxt0,loph,losp,moph,mosp,mprph,mprsp,mrgexs,mtb,mtb0,
     $      nbasp,nbaspd,nbt,nbtmax,ncmpe,ncmpr,ndelay,netmax,
     $      ngext,nodbmx,nordmx,noutpt,npet,npetmx,npet0,npt,nptmax,
     $      npts,nset,nsetmx,nset0,nstmax,nsts,nstsmx,nstsr,nttyo,
     $      qbye,uaqsln,ufixf,uspec,uphase,uzvec0,uzvec1,xbar,
     $      xbarlg,zklgmn,zklogl,zvclg0,zvclg1,zvec0,zvec1)
c
c           Note: if such a shift has been made here, the maximum
c           order of the finite-differences (kord) has been set
c           to zero.
c
            qmod2 = .false.
c
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the step size and the order for the next increment. The
c     order may be dropped if it appears that the high-order derivative
c     estimates are not numerically significant. In economy mode, this
c     algorithm is used to choose the order, but another algorithm (see
c     below) is used to choose the step size.
c
  110 if (qzdump) delxi = max(dlxis2,delxi)
c
      qzprnt = .false.
      qprntx = .false.
      qprntt = .false.
      qprnlx= .false.
      qprnlt= .false.
c
      qzplot = .false.
      qplotx = .false.
      qplott = .false.
      qplolx = .false.
      qplolt = .false.
c
      qprph0 = .false.
      qprph1 = .false.
      qpreh0 = .false.
      qpreh1 = .false.
      qpro20 = .false.
      qpro21 = .false.
      qpraw0 = .false.
      qpraw1 = .false.
c
      qplph0 = .false.
      qplph1 = .false.
      qpleh0 = .false.
      qpleh1 = .false.
      qplo20 = .false.
      qplo21 = .false.
      qplaw0 = .false.
      qplaw1 = .false.
c
      qzdump = .false.
      qstabl = .false.
      qsawth = .false.
      if (qskip) go to 120
      if (qshoot) go to 125
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qecon) then
c
c       Choose the step size and order when operating in economy
c       mode.
c
        call chdxec(delxi,dlxmx0,dzvc0,iodb,kdim,kmax,km1,kxt,
     $  nodbmx,nord,noutpt,nrd1mx,qmin,qscon,scale,scalim,scnstd,
     $  scnsti,uzvec1,zvec0)
c
      else
c
c       Choose the step size and order when operating in normal mode.
c       Here control is based on the Gear approach. An order greater
c       than zero is evaluated on the basis of the contributions
c       of terms of one higher order. A hidden order is carried
c       in the finite differences and derivatives of the z and r
c       vectors to facilitate this approach.
c
c       Get choices based on finite differences for the z vector.
c
        call chdxgz(delxi,dlxmx0,fdlim,fdzv0,iodb,iopt,jordlm,
     $  kdim,kmax,km1,kord,kxt,nodbmx,noptmx,nordmx,nordz,noutpt,
     $  nrd1mx,nsscmx,scalim,scfcz,sscrew,qmin,smp100,uzvec1,zklogu,
     $  zvec0,zvclg0)
c
c       Get choices based on finite differences for the r vector.
c
        call chdxgr(delxi,dlxmx0,fdri0,fdrr0,iodb,jordlm,jreac,
     $  kord,nodbmx,nordmx,nordr,noutpt,nrct,nrctmx,nrd1mx,nsscmx,
     $  scalim,scfcr,sscrew,qriinf,rirec0,rrelr0,ureac)
c
        qstabz = .false.
        qstabr = .false.
c
        if (iopt(2).gt.0 .and. qsawth) then
          if (kord .ge. 2) then
c
c           Compute derivatives for the z and r vectors from the
c           corresponding finite differences.
c
            kordp1 = kord + 1
c
c           First compute the matrix (akmat0 = F * B) used to
c           calculate derivatives for truncated Taylor's series,
c           given predictor-based finite differences.
c
c           Calling sequence substitutions:
c             kordp1 for nord
c
            call gakmat(akmat0,dxsm00,kordp1,nrd1mx)
c
c           Now compute the corresponding derivatives.
c
c           Calling sequence substitutions:
c             kordp1 for nord
c
            call zderiv(akmat0,dzvc0,fdzv0,kdim,kmax,kordp1,nrd1mx)
c
c           Calling sequence substitutions:
c             kordp1 for nord
c
            call rderiv(akmat0,drer0,drir0,fdri0,fdrr0,jreac,kordp1,
     $      nrct,nrctmx,nrd1mx)
c
c           Now compute the corresponding average derivatives.
c
            delxia = dxsm00(2)
c
c           Calling sequence substitutions:
c             kordp1 for nord
c
            call gsmdez(delxia,dzvc0,dzvc0s,kdim,kmax,kordp1,nrd1mx)
c
c           Calling sequence substitutions:
c             kordp1 for nord
c
            call gsmder(delxia,drer0,drer0s,drir0,drir0s,jreac,kordp1,
     $      nrct,nrctmx,nrd1mx)
c
c           Choose step sizes and orders.
c
            call chdxtz(delxi,dlxmx0,dzvc0,iodb,iopt,jordlm,kdim,
     $      kmax,km1,kord,kxt,nodbmx,noptmx,nordmx,nordzs,noutpt,
     $      nrd1mx,nsscmx,scalim,scfczs,sscrew,qmin,smp100,uzvec1,
     $      zklogu,zvec0,zvclg0)
c
            call chdxtr(delxi,dlxmx0,drer0,drir0,iodb,jordlm,jreac,
     $      kord,nodbmx,nordmx,nordrs,noutpt,nrct,nrctmx,nrd1mx,
     $      nsscmx,scalim,scfcrs,sscrew,qriinf,rirec0,rrelr0,ureac)
c
c           Resolve regular vs. "stable" results.
c
            if (scfczs .gt. scfcz) then
              scfcz = scfczs
              nordz = nordzs
              qstabz = .true.
            endif
c
            if (scfcrs .gt. scfcr) then
              scfcr = scfcrs
              nordr = nordrs
              qstabr = .true.
            endif
c
          endif
        endif
c
c       Resolve z vector vs. r vector results.
c
        scale = scfcz
        nord = nordz
        qstabl = qstabz
        if (scfcr .lt. scfcz) then
          scale = scfcr
          nord = nordr
          qstabl = qstabr
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate derivatives from finite differences. These are
c     consistent with the actual order to be used (nord).
c
      if (nord .gt. 0) then
c
c       First compute the matrix (akmat0 = F * B) used to calculate
c       derivatives for truncated Taylor's series, given predictor-
c       based finite differences.
c
        call gakmat(akmat0,dxsm00,nord,nrd1mx)
c
c       Now compute the derivatives.
c
        call zderiv(akmat0,dzvc0,fdzv0,kdim,kmax,nord,nrd1mx)
        if (iopt(2) .gt. 0) then
          call rderiv(akmat0,drer0,drir0,fdri0,fdrr0,jreac,nord,
     $    nrct,nrctmx,nrd1mx)
        endif
c
        call aderiv(akmat0,daffp0,fdaf0,nord,nordmx,npt,
     $  nptmax,nrd1mx)
        call bderiv(akmat0,dafrc0,fdar0,jreac,nord,nordmx,nrct,
     $  nrctmx,nrd1mx)
        call pderiv(akmat0,demop0,fdpe0,nord,nordmx,npet,
     $  npetmx,nrd1mx)
        call sderiv(akmat0,demos0,fdse0,nord,nordmx,nrd1mx,
     $  nset,nsetmx)
c
c       Calling sequence substitutions:
c         dph0 for dxx0
c         fdph0 for fdxx0
c
        call xderiv(akmat0,dph0,fdph0,nord,nordmx,nrd1mx)
c
        if (qredox) then
c
c         Calling sequence substitutions:
c           deh0 for dxx0
c           fdeh0 for fdxx0
c
          call xderiv(akmat0,deh0,fdeh0,nord,nordmx,nrd1mx)
c
c         Calling sequence substitutions:
c           do20 for dxx0
c           fdo20 for fdxx0
c
          call xderiv(akmat0,do20,fdo20,nord,nordmx,nrd1mx)
        endif
c
c       Calling sequence substitutions:
c         daw0 for dxx0
c         fdaw0 for fdxx0
c
        call xderiv(akmat0,daw0,fdaw0,nord,nordmx,nrd1mx)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        if (qstabl) then
          delxia = dxsm00(nord)
c
c         Compute average derivatives.
c
          call gsmdez(delxia,dzvc0,dzvc0s,kdim,kmax,nord,nrd1mx)
          call gsmder(delxia,drer0,drer0s,drir0,drir0s,jreac,nord,
     $    nrct,nrctmx,nrd1mx)
c
c         Replace the actual with the average derivatives.
c
          do kcol = 1,kdim
            do n = 1,nord
              dzvc0(n,kcol) = dzvc0s(n,kcol)
            enddo
          enddo
          do n = 1,nord
            drir0(n) = drir0s(n)
          enddo
          do nrc = 1,nrct
            do n = 1,nord
              drer0(n,nrc) = drer0s(n,nrc)
            enddo
          enddo
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the new value of delxi.
c
      delxi = delxi*scale
c
c     If the step size was previously cut because a phase boundary
c     was overstepped by an unacceptable amount, limit delxi so that
c     the new value of Xi will be less than the value at the last
c     known point of overstep.
c
      dlxlim = 0.75*(xilim - xi0)
      if (dlxlim .lt. dlxmin) then
c
c       The calculation has has now gotten very close to the previously
c       determined overstep location. After the present step is
c       completed, the calculation will be within a minimum step size
c       distance of that location. Failure to detect either the
c       phase boundary or the overstep at this point means that the
c       overstep location previously calculated is not actually valid.
c       What this means is that as the calculation gets closer to the
c       phase boundary, its calculated position is moving to a greater
c       value of reaction progress. This can occur when the mass
c       transfer associated with a step must be calculated from an
c       approximate, not an exact integration. Such is the case for
c       kinetic mode calculations. When the calculation goes back
c       and cuts the step size, the resulting calculation of the
c       position of the phase boundary then becomes more accurate.
c
        dlxlim = dlxmin
        xilim = prcinf
      endif
      delxi = min(delxi,dlxlim)
c
c     Make sure that delxi is not less than the minimum normally
c     allowed value.
c
      delxi = max(delxi,dlxmin)
c
      if (iodb(5) .ge. 1) write (noutpt,1320) nord,delxi
 1320 format(/' --- Selection: nord= ',i2,', delxi= ',1pe11.4,
     $ ' ---')
c
  120 qskip = .false.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set other simple limits on the step size.
c
      if (iopt(1) .eq. 2) then
        dlxis2 = delxi
        dxdmp = xidump - xi0
        if (delxi .gt. dxdmp) then
          dxdmp = max(dxdmp,dlxmin)
          delxi = dxdmp
          if (iodb(5) .ge. 2) write (noutpt,1330) delxi
 1330     format(' --- delxi cut to ',1pe11.4,' to match the Xi PRS',
     $    ' transfer interval ---')
        endif
      endif
c
      dlxisv = delxi
      dx1 = xiprnt - xi0
      dx2 = xiprnl - xi0
      dlxipr = min(dx1,dx2)
      if (delxi .gt. dlxipr) then
        delxi = dlxipr
        if (iodb(5) .ge. 2) write (noutpt,1340) delxi
 1340   format(' --- delxi cut to ',1pe11.4,' to match the Xi print',
     $  ' interval ---')
      endif
c
      dx1 = xiplot - xi0
      dx2 = xiplol - xi0
      dlxipl = min(dx1,dx2)
      if (delxi .gt. dlxipl) then
        delxi = dlxipl
        if (iodb(5) .ge. 2) write (noutpt,1350) delxi
 1350   format(' --- delxi cut to ',1pe11.4,' to match the Xi plot',
     $  ' interval ---')
      endif
c
      dlximx = ximax - xi0
      if (dlximx .lt. delxi) then
        delxi = dlximx
        if (iodb(5) .ge. 2) write (noutpt,1360) delxi
 1360   format(' --- delxi cut to  ',1pe11.4,' to match the maximum',
     $  ' value of Xi  ---')
      endif
c
      dlxilm = dlxmax
      if (nord .le. 0) dlxilm = dlxmx0
      if (dlxilm .lt. delxi) then
        delxi = dlxilm
        if (iodb(5) .ge. 2) write (noutpt,1370) delxi
 1370   format(' --- delxi cut to ',1pe11.4,' to match the upper limit',
     $  ' for any step ---')
      endif
c
  125 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nord.gt.0 .and. .not.qscon) then
c
c       Limit delxi when the variables corresponding to certain basis
c       species are changing rapidly. Such variables may include the pH
c       or the log fO2. This limit allows some minimum of information
c       density to be obtained in the neighborhood of the rapid
c       change.
c
        call ldlxrc(al10,delxi,dlxmin,dzvc0,iodb,iindx1,kbt,kdim,
     $  kelect,khydr,khydx,km1,kmax,ko2gaq,krdxsp,kwater,kxt,nbasp,
     $  nbtmax,nodbmx,nord,noutpt,nrd1mx,nstmax,nttyo,qrapch,uspec,
     $  zklogu,zvclg0,zvclg1,zvec0,zvec1)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nord.gt.0 .and. .not.qscon .and. iopt(3).lt.1
     $  .and. kxt.ge.km1 .and. iopt(1).le.1) then
c
c       Limit delxi when a phase is rapidly disapparing from the
c       equilibrium system. This is mechanism is only designed to
c       help other mechanisms efficiently locate a phase disappearance
c       boundary. It is not designed to increase information density
c       near the boundary.
c
        call ldlxrd(delxi,dlxmin,dzvc0,fdzv0,iodb,ipndx1,kdim,
     $  km1,kmax,kxt,loph,nodbmx,nord,noutpt,nptmax,nrd1mx,nttyo,
     $  uphase,zklogu,zvec0)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      qdump = .false.
      if (iopt(1).eq.2 .and. nord.ge.1 .and. npet.gt.0) then
c
c       Limit delxi by approximate position of significant maxima
c       of mineral masses. Because the cost of restricting delxi to
c       keep high accuracy in the estimates of first derivatives is
c       prohibitive, control can not be based strictly on derivatives,
c       but must be based mainly on how much of a mineral mass would
c       be destroyed on a given step. This permits avoiding maxima
c       which are false (due to Taylor's series inaccuracy) or trivial
c       (if the mass of a mineral is very small). If the Taylor's
c       series fail to detect a phase whose mass is going over a maximum
c       and iopt(1) = 2, then go back instructions in EQ6/eqcalc.f will
c       prevent a skip over the maximum.
c
        call fpbflo(al10,delxi,demop0,dlxmin,dxval0,d1emp1,
     $  d2emp1,emop,emop0,eps100,fdpe0,iemop,ier,iodb,nodbmx,nord,
     $  nordmx,noutpt,npet,npetmx,nptmax,nrd1mx,nttyo,qdump,toldl,
     $  uaqsln,ufixf,uphase,xim1,xi0,xi1,xval0,zklogu)
c
        if (ier .gt. 0) then
c
c         A search for a maximum in subroutine fpbflo failed. Go back
c         to the base point, shift mass to the PRS, and drop back
c         to order zero.
c
          qdump = .true.
          kord = 0
          nord = 0
          npts = 1
          delxi = 0.
          xi1 = xi0
          do npe = 1,npet
            emop(npe) = emop0(npe)
          enddo
          do nse = 1,nset
            emos(nse) = emos0(nse)
          enddo
          do kcol = 1,kxt
            zvclg1(kcol) = zvclg0(kcol)
            zvec1(kcol) = zvec0(kcol)
          enddo
c
          call pshfta(csts,emop,emop0,emos,emos0,fdpe0,fdpem1,
     $    fdse0,fdsem1,iemop,iemos,iern1,iern2,ietmax,iindx1,imrn1,
     $    imrn2,iodb,ipndx1,ixrn1,ixrn2,jcsort,jern1,jetmax,jgext,
     $    jpflag,jsflag,kbt,km1,kmax,kmt,kx1,kxt,loph,losp,moph,mosp,
     $    mprph,mprsp,mrgexs,mtb,mtb0,nbasp,nbaspd,nbt,nbtmax,ncmpe,
     $    ncmpr,netmax,ngext,nodbmx,nordmx,noutpt,npet,npetmx,npt,
     $    nptmax,nsetmx,nstmax,nsts,nstsmx,nstsr,nttyo,uaqsln,ufixf,
     $    uphase,uspec,xbar,xbarlg,zklgmn,zklogl,zvclg0,zvclg1,
     $    zvec0,zvec1)
c
          delxi = dlxmx0
          go to 110
        endif
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nord.gt.0 .and. .not.qscon) then
c
c       Find the phase boundary for a newly appearing phase. This is
c       here defined as a point at which the affinity to precipitate
c       has a target value (aftarg) that is slightly greater than zero.
c       If the Taylor's series fail to detect the phase boundary at
c       which a new phase appears, this condition will be subsequently
c       detected by a test in EQLIB/newton.f. This subroutine will then
c       issue a go back instruction, which will prevent the boundary
c       from being skipped over.
c
        call fpbnpp(affp,affp0,aftarg,daffp0,delxi,dlxmin,dxval0,
     $  eps100,iodb,iopt,jpflag,nodbmx,noptmx,nord,nordmx,noutpt,npchk,
     $  npt,nptmax,nrd1mx,nttyo,tolaft,tolsat,uphase,xi0,xi1,xval0)
c
        if (iopt(1).ne.2 .and. npet.gt.0) then
c
c         Find the phase boundary for a disappearing phase in the ES.
c
          call fpbdpp(delxi,demop0,dlxmin,dxval0,emop,emop0,
     $    eps100,iemop,iodb,iopt,nodbmx,noptmx,nord,nordmx,noutpt,
     $    npet,nrd1mx,npetmx,nptmax,nttyo,uphase,xi0,xi1,xval0)
        endif
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nrct.gt.0 .and. nord.gt.0) then
        if (aft0 .gt. 0.01) then
c
c         Check the signs of the predicted reactant affinities. If
c         necessary, reduce delxi so that none of the affinities
c         change sign. For the purposes of determining a crossover,
c         a reactant affinity is taken to be zero if its magnitude is
c         less than or equal to tolsar.
c
c         This check is made to help deal with rate expressions that
c         do not approach zero as the corresponding equilibria are
c         approached. For example, one might specify a constant rate
c         of dissolution that applies to any value of the corres-
c         ponding affinity as long as it implies undersaturation. It
c         is important not to extrapolate such a rate to any other
c         condition (in this example, saturation or supersaturation).
c
c         If close to overall equilibrium, this check is best avoided.
c         Here the rates of "badly behaving" reactions may tend to
c         oscillate about their equilibrium points. As long as the
c         magnitudes of the oscillations are small, it is better to
c         not try to pin down the cross-overs.
c
          call chksar(afrc0,afrcp,dafrc0,delxi,dlxmin,dxval0,
     $    eps100,iodb,jreac,nodbmx,noutpt,nord,nordmx,nrct,nrctmx,
     $    nrd1mx,nrk,nttyo,tolsar,ureac,xi0,xi1,xval0)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(2) .gt. 0) then
        if (nrct.gt.0 .and. nord.gt.0) then
          if (aft0 .gt. 0.01) then
c
c           Check the sign of the predicted inverse rate. If necessary,
c           reduce delxi so that the inverse rate does not become
c           less than eps100. Physically, it makes no sense that
c           the inverse rate should approach zero or become negative,
c           but errors in finite differences that could cause a
c           predicted inverse rate to do this. If so, the step size
c           needs to be reduced.
c
            call chksir(delxi,dlxmin,drir0,dxval0,eps100,iodb,
     $      nodbmx,nord,noutpt,nrd1mx,nttyo,rirec0,rirecp,xi0,xi1,xval0)
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nrct.gt.0 .and. nord.gt.0) then
        if (aft0 .gt. 0.01) then
c
c         Check the signs of the predicted relative rates. If,
c         necessary reduce delxi so that none of the relative rates
c         change sign. For the purposes of determining a crossover,
c         a relative rate is taken to be zero if its magnitude is
c         less than tolsrr.
c
c         This check is best avoided near overall equilibrium. See
c         comments made above in regard to a similar check on sign
c         changes of reactant affinities.
c
          call chksrr(delxi,dlxmin,drer0,dxval0,eps100,iodb,jreac,
     $    nodbmx,noutpt,nord,nrct,nrctmx,nrd1mx,nttyo,rrelr0,rrelrp,
     $    tolsrr,ureac,xi0,xi1,xval0)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(2) .gt. 0) then
        if (nrct .gt. 0) then
c
c         Make absolutely sure that the time increment is positive.
c         Errors in numerical integration have the potential to
c         produce a negative or zero result.
c
          call chksti(akmat0,drer0,drir0,deltim,delxi,dlxmin,
     $    fdri0,fdrr0,iodb,jreac,kly,kmax,nodbmx,kord,nord,noutpt,
     $    npts,nrct,nrctmx,nrd1mx,nttyo,prcinf,qriinf,rirec0,smp100,
     $    time0,time1,xi0,xi1)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nrct .gt. 0) then
c
c       Limit delxi by exhaustion of reactants.
c
        call fpexrc(delxi,dlxmin,drer0,dxval0,eps100,iodb,
     $  jreac,morr,morr0,nodbmx,nord,noutpt,nrct,nrctmx,nrd1mx,
     $  nttyo,qdump,rrelr0,ureac,xi0,xi1,xval0)
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(2) .gt. 0) then
c
c       Limit delxi so that the next time-based print point is not
c       exceeded.
c
        if (qtprn) then
          call chktpr(delxi,dlxmin,dlxtpr,drir0,dxval0,eps100,
     $    iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,prcinf,qdump,qriinf,
     $    rirec0,tiprnl,tiprnt,time0,time1,tolxst,xi0,xi1,xval0)
        endif
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Limit delxi so that the next pH-based print point is not exceeded.
c
      if (dlhprn .lt. prcinf) then
        call ckphpr(delxi,dlxmin,dph0,dxh0pr,dxh1pr,dxval0,
     $  eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,ph0,ph1,
     $  ph0prn,ph1prn,prcinf,qdump,tolxsu,xi0,xi1,xval0)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qredox) then
c
c       Limit delxi so that the next Eh-based print point is not
c       exceeded.
c
        if (dlhprn .lt. prcinf) then
          call ckehpr(delxi,dlxmin,deh0,dxe0pr,dxe1pr,dxval0,
     $    eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,eh0,eh1,
     $    eh0prn,eh1prn,prcinf,qdump,tolxsu,xi0,xi1,xval0)
        endif
c
c       Limit delxi so that the next log fO2-based print point is not
c       exceeded.
c
        if (dloprn .lt. prcinf) then
          call cko2pr(delxi,dlxmin,do20,dxo0pr,dxo1pr,dxval0,
     $    eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,fo2lg0,
     $    fo2lg1,o20prn,o21prn,prcinf,qdump,tolxsu,xi0,xi1,xval0)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Limit delxi so that the next aw-based print point is not exceeded.
c
      if (dlaprn .lt. prcinf) then
        call ckawpr(delxi,dlxmin,daw0,dxw0pr,dxw1pr,dxval0,
     $  eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,aw0,aw1,
     $  aw0prn,aw1prn,prcinf,qdump,tolxsu,xi0,xi1,xval0)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(2) .gt. 0) then
c
c       Limit delxi so that the next time-based plot point is not
c       exceeded.
c
        if (qtplo) then
          call chktpl(delxi,dlxmin,dlxtpl,drir0,dxval0,eps100,
     $    iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,prcinf,qdump,qriinf,
     $    rirec0,tiplol,tiplot,time0,time1,tolxst,xi0,xi1,xval0)
        endif
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Limit delxi so that the next pH-based plot point is not exceeded.
c
      if (dlhplo .lt. prcinf) then
        call ckphpl(delxi,dlxmin,dph0,dxh0pl,dxh1pl,dxval0,
     $  eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,ph0,ph1,
     $  ph0plo,ph1plo,prcinf,qdump,tolxsu,xi0,xi1,xval0)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qredox) then
c
c       Limit delxi so that the next Eh-based plot point is not
c       exceeded.
c
        if (dlhplo .lt. prcinf) then
          call ckehpl(delxi,dlxmin,deh0,dxe0pl,dxe1pl,dxval0,
     $    eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,eh0,eh1,
     $    eh0plo,eh1plo,prcinf,qdump,tolxsu,xi0,xi1,xval0)
        endif
c
c       Limit delxi so that the next log fO2-based plot point is not
c       exceeded.
c
        if (dloplo .lt. prcinf) then
          call cko2pl(delxi,dlxmin,do20,dxo0pl,dxo1pl,dxval0,
     $    eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,fo2lg0,
     $    fo2lg1,o20plo,o21plo,prcinf,qdump,tolxsu,xi0,xi1,xval0)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Limit delxi so that the next aw-based plot point is not exceeded.
c
      if (dlaplo .lt. prcinf) then
        call ckawpl(delxi,dlxmin,daw0,dxw0pl,dxw1pl,dxval0,
     $  eps100,iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,aw0,aw1,
     $  aw0plo,aw1plo,prcinf,qdump,tolxsu,xi0,xi1,xval0)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(2) .gt. 0) then
c
c       Limit delxi so that the requested maximum value of time
c       is not exceeded.
c
        if (timemx .lt. prcinf) then
          call chktmx(delxi,dlxmin,dlxtmx,drir0,dxval0,eps100,
     $    iodb,nodbmx,nord,noutpt,nrd1mx,nttyo,prcinf,qdump,qriinf,
     $    rirec0,timemx,time0,time1,tolxst,xi0,xi1,xval0)
        endif
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Limit delxi so that the requested minimum value of pH
c     is not exceeded.
c
      if (phmin .gt. -prcinf) then
        call ckphmn(delxi,dlxmin,dph0,dxh0mx,dxval0,eps100,iodb,
     $  nodbmx,nord,noutpt,nrd1mx,nttyo,ph0,ph1,phmin,prcinf,qdump,
     $  tolxsu,xi0,xi1,xval0)
      endif
c
c     Limit delxi so that the requested maximum value of pH
c     is not exceeded.
c
      if (phmax .lt. prcinf) then
        call ckphmx(delxi,dlxmin,dph0,dxh1mx,dxval0,eps100,iodb,
     $  nodbmx,nord,noutpt,nrd1mx,nttyo,ph0,ph1,phmax,prcinf,qdump,
     $  tolxsu,xi0,xi1,xval0)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qredox) then
c
c       Limit delxi so that the requested minimum value of Eh
c       is not exceeded.
c
        if (ehmin .gt. -prcinf) then
          call ckehmn(delxi,dlxmin,deh0,dxe0mx,dxval0,eps100,iodb,
     $    nodbmx,nord,noutpt,nrd1mx,nttyo,eh0,eh1,ehmin,prcinf,qdump,
     $    tolxsu,xi0,xi1,xval0)
        endif
c
c       Limit delxi so that the requested maximum value of Eh
c       is not exceeded.
c
        if (ehmax .lt. prcinf) then
          call ckehmx(delxi,dlxmin,deh0,dxe1mx,dxval0,eps100,iodb,
     $    nodbmx,nord,noutpt,nrd1mx,nttyo,eh0,eh1,ehmax,prcinf,qdump,
     $    tolxsu,xi0,xi1,xval0)
        endif
c
c       Limit delxi so that the requested minimum value of log fO2
c       is not exceeded.
c
        if (o2min .gt. -prcinf) then
          call cko2mn(delxi,dlxmin,do20,dxo0mx,dxval0,eps100,iodb,
     $    nodbmx,nord,noutpt,nrd1mx,nttyo,fo2lg0,fo2lg1,o2min,prcinf,
     $    qdump,tolxsu,xi0,xi1,xval0)
        endif
c
c       Limit delxi so that the requested maximum value of log fO2
c       is not exceeded.
c
        if (o2max .lt. prcinf) then
          call cko2mx(delxi,dlxmin,do20,dxo1mx,dxval0,eps100,iodb,
     $    nodbmx,nord,noutpt,nrd1mx,nttyo,fo2lg0,fo2lg1,o2max,prcinf,
     $    qdump,tolxsu,xi0,xi1,xval0)
        endif
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Limit delxi so that the requested minimum value of the activity
c     of water is not exceeded.
c
      if (awmin .gt. -prcinf) then
        call ckawmn(delxi,dlxmin,daw0,dxw0mx,dxval0,eps100,iodb,
     $  nodbmx,nord,noutpt,nrd1mx,nttyo,aw0,aw1,awmin,prcinf,qdump,
     $  tolxsu,xi0,xi1,xval0)
      endif
c
c     Limit delxi so that the requested maximum value of the activity
c     of water is not exceeded.
c
      if (awmax .lt. prcinf) then
        call ckawmx(delxi,dlxmin,daw0,dxw1mx,dxval0,eps100,iodb,
     $  nodbmx,nord,noutpt,nrd1mx,nttyo,aw0,aw1,awmax,prcinf,qdump,
     $  tolxsu,xi0,xi1,xval0)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make sure that delxi is not less than dlxmin unless one of the
c     following special conditions is satisfied.
c
      if (abs(delxi - dlxipr).gt.0. .and. abs(delxi - dlxipl).gt.0.
     $ .and. abs(delxi - dlximx).gt.0. .and. abs(delxi - dlxtmx).gt.0.)
     $ delxi = max(delxi,dlxmin)
c
      ncut = 0
      jcut = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  130 xi1 = xi0 + delxi
c
      if (.not.qshoot) then
c
c       The following trap keeps the code from doing lengthy cycling
c       with delxi equal to zero.
c
        if (delxi .le. 0.) then
          kstpze = kstpze + 1
          qx = kstpze .ge. 2
          if (qx) then
            write (noutpt,1380) kstpze
            write (nttyo,1380) kstpze
 1380       format(/' * Error - (EQ6/path) The step size has been zero',
     $      ' now ',i2,' times in a row.')
            stop
          endif
        else
          kstpze = 0
        endif
c
c       The following trap keeps the code from doing lengthy cycling
c       with  delxi equal to dlxmin.
c
        if (delxi .le. dlxmin) then
          kstpmn = kstpmn + 1
          qx = kstpmn .ge. 100
          if (qx) then
            write (noutpt,1390) kstpmn
            write (nttyo,1390) kstpmn
 1390       format(/' * Error - (EQ6/path) The step size has not',
     $      ' exceeded the minimum',/7x,'value now ',i3,
     $      ' times in a row.')
            stop
          endif
        else
          kstpmn = 0
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The following is a return point for doing an ODE corrector
c     iteration.
c
  140 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nrct .gt. 0) then
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
      if (iopt(2) .gt. 0) then
c
        if (qshoot) then
          time1 = prcinf
          deltim = prcinf
        else
c
c         Increment the time and related functions.
c
          call timeca(deltim,delxi,drir0,iodb,nodbmx,nord,noutpt,
     $    nrd1mx,nttyo,prcinf,qriinf,rirec0,time0,time1)
        endif
c
        if (delxi .le. dlxmin) then
c
c         Make sure that the calculated time does not exceed any
c         any specified limits such as the maximum time just because
c         delxi is at the minimum value.
c
          call tivchk(deltim,delxi,qtvchk,time1,time0,timemx,
     $    tiplol,tiplot,tiprnl,tiprnt,tolxst)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
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
        qtrch = ntpr0 .ne. ntpr
        if (qtrch) nord = 0
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (delxi .le. dlxmin) then
        nord = 0
        kord = 0
      endif
c
c     Make a Taylor's series expansion of the z vector, applying
c     change limits. It is important to make a protected expansion
c     here to clear out unrecoverable values in the zvclg1 and
c     zvec1 vectors that may have been generated in previous
c     order/step size adjustments for the current step. This is
c     true even if the current order is zero, as the order may
c     have been higher earlier in the current step.
c
      qztayl = .true.
      call ztaylr(delxi,dzvc0,kdim,kmax,km1,kxt,nord,nrd1mx,
     $ qztayl,zklogu,zvclg0,zvclg1,zvec0,zvec1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      qstart = xi1 .le. xistsv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make the equilibrium calculation at the current point of reaction
c     progress. Subroutine eqshel forms a shell around EQ6/eqphas.f.
c     It may make small advances in the value of the reaction progress
c     variable in order to step over a small region (commonly on the
c     order of 1 x 10-8 mol) in which the equilibrium state can not be
c     calculated satisfying the normal tolerances. One category of such
c     regions is associated with phase boundaries (where a phase such
c     as a mineral either appears or disappears). The problem is that
c     a calculation assuming the phase is present crashes, indicating
c     that the phase should be deleted. One the other hand, a
c     calculation assuming it is not present converges, but gives a
c     supersaturation in excess of the supersaturation tolerance (a
c     condition which normally instructs the code to add the phase to
c     the phase assemblage and try again). Another category of such
c     regions is associated with a rapid change in the value of a
c     variable associated with an aqueous basis species. Typically
c     this is the oxygen fugacity or equivalent redox variable. Such
c     a change is called a redox jump. A similar jump might occur for
c     some non-redox variable, such as pH. However, this is not likely
c     for any non-redox variable because of the general presence of
c     significant buffering effects of one kind or another (true
c     buffering, mass buffering).
c
      call eqshel(aadh,aadhh,aadhv,aamatr,aaphi,abar,abdh,
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
c     Note on ier codes returned by EQ6/eqshel.f:
c
c       =    0  Okay
c       =   10  Go back and take a smaller step size to avoid exceeding
c                 the supersaturation tolerance (tolsst)
c       =  170  Too much of a phase was destroyed under the flow-through
c                 open system model; go back and first move part of the
c                 mass of protected phases in the ES to the PRS
c       =  180  One of a number of problems occurred which may be
c                 resolvable, at least partially, by going back and
c                 cutting the step size
c       =  190  Need to slide over a region of computational
c                 instability, but sliding is inhibited; go back, but
c                 terminate work on the current problem
c
      if (ier .le. 0) go to 200
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ier .eq. 10) then
c
c       An incoming phase boundary was stepped over by too much.
c       Go back and step forward again, using a smaller step size.
c       Note: this error code is not returned if the step size is less
c       than or equal to the minimum value.
c
        xilim = xi1
        kly = 6
        go to 150
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set up to go back to the previous point. Try a smaller step
c     size or zero order Taylor's expansions, or write a pickup file
c     for the last good point and stop.
c
  150 continue
c
      call goback(acflg,acflg0,emop,emop0,emos,emos0,fje,
     $ fje0,fxi,fxi0,iemop,iemop0,iemos,iemos0,iindx0,iindx1,ipndx0,
     $ ipndx1,jpflag,jsflag,jreac,jreac0,kdim,kdim0,kmax,km1,km10,
     $ kmt,kmt0,kx1,kx10,kxt,kxt0,loph,losp,moph,moph0,mosp,mosp0,
     $ ncmpe,ncmpe0,npet,npetmx,npet0,npt,nptmax,nrct,nrctmx,nset,
     $ nsetmx,nset0,nst,nstmax,qreq,qriinf,sigmam,sigmm0,uzvec0,
     $ uzvec1,xi0,xi1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ier.eq.180 .and. delxi.le.dlxmin) then
c
c       Normally this error code results in a reduction in the step
c       size. This can't be done if the step size is already at the
c       minimum value. Write a pickup file at the last good point
c       and terminate work on the current problem.
c
        write (ux16,'(g12.5)') xi1
        call lejust(ux16)
        j2 = ilnobl(ux16)
        write (noutpt,1400) ux16(1:j2)
        write (nttyo,1400) ux16(1:j2)
 1400   format(/' * Error- (EQ6/path) The equilibrium calculation',
     $  ' failed at Xi= ',a,'.',/7x,"Can't cut the step size to try",
     $  ' to recover. See previous notes',/7x,'and warnings for',
     $  ' suggestions.')
        go to 180
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     If a slide forward in EQ6/eqshel.f failed to get over a region of
c     computational instability, write a pickup file at the last good
c     point and terminate work on the current problem.
c
      if (ier .eq. 190) go to 180
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ier .eq. 170) then
c
c       The equilibrium calculation was completed successfully, but too
c       much of a phase in the equilibrium system (ES) was destroyed to
c       satisfy the requirements of the fluid-centered flow-through open
c       system model. Go back and first transfer some of the phase's
c       mass to the physically removed system (PRS).
c
        km1 = km10
        kmt = kmt0
        kx1 = kx10
        kxt = kxt0
        kdim = kdim0
c
        call copyaa(zvclg0,zvclg1,kdim)
        call copyaa(zvec0,zvec1,kdim)
c
        do kcol = 1,kdim
          uzvec1(kcol) = uzvec0(kcol)
          iindx1(kcol) = iindx0(kcol)
          ipndx1(kcol) = ipndx0(kcol)
        enddo
c
        dxsave = 0.25*delxi
        kly = 6
        delxi = 0.
        deltim = 0.
        time1 = time0
c
c       Recompute the temperature and pressure. Then recompute the
c       thermodynamic and kinetic quantities which depend these
c       variables.
c
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
        if (nrct .ge. 1) then
c
c         Increment the irreversible reactions (those associated with
c         the "reactants"); update the ES mass balance totals
c         accordingly.
c
          call reacts(cbsr,csts,delxi,drer0,iern1,ietmax,iktmax,
     $    iodb,jcode,jetmax,jgext,jreac,modr,modr0,morr,morr0,mrgers,
     $    mtb,mtb0,nbaspd,nbt,nbtmax,nbt1mx,ncmpr,nern1,nern2,nertmx,
     $    netmax,ngext,nodbmx,nord,noutpt,nptmax,nrct,nrctmx,nrd1mx,
     $    nrndex,nsrtmx,nstmax,nsts,nstsmx,nstsr,nttyo,nxridx,nxrtmx,
     $    rrelr0,rxbar,ureac,xirct,xirct0)
        endif
c
c       Expand the system description.
c
        call ncmpex(acflg,act,actlg,cdrs,cegexs,cgexj,conc,
     $  conclg,cpgexs,egexjc,egexjf,egexs,eps100,fo2,fo2lg,fsort,
     $  fugac,fugalg,iern1,iern2,ietmax,ifrn1,ifrn2,igas,igstak,
     $  iindx1,ilrn1,ilrn2,imrn1,imrn2,istack,ixrn1,ixrn2,jcsort,
     $  jern1,jern2,jetmax,jflag,jgext,jgsort,jgstak,jjsort,jpflag,
     $  jsflag,jsitex,jssort,jstack,kbt,kdim,kelect,kmax,km1,ko2gaq,
     $  kwater,kxt,loph,losp,lsort,mgext,mrgexs,mtb,moph,mosp,narn1,
     $  narn2,nbasp,nbt,nbtmax,ncmpr,ndrs,ndrsmx,ndrsr,nelect,nern1,
     $  nern2,netmax,ngexsa,ngext,ngrn1,ngrn2,ngt,ngtmax,noutpt,
     $  no2gaq,nphasx,npt,nptmax,nst,nstmax,nttyo,omega,omeglg,
     $  press,qxbarw,q6mode,ugexj,ugexmo,uphase,uspec,xbar,xbarlg,
     $  xbarw,xbarwc,xbrwlc,xbrwlg,xlks,zchar,zgexj,zvclg1,zvec1)
c
        call pshfta(csts,emop,emop0,emos,emos0,fdpe0,fdpem1,
     $  fdse0,fdsem1,iemop,iemos,iern1,iern2,ietmax,iindx1,imrn1,
     $  imrn2,iodb,ipndx1,ixrn1,ixrn2,jcsort,jern1,jetmax,jgext,
     $  jpflag,jsflag,kbt,km1,kmax,kmt,kx1,kxt,loph,losp,moph,mosp,
     $  mprph,mprsp,mrgexs,mtb,mtb0,nbasp,nbaspd,nbt,nbtmax,ncmpe,
     $  ncmpr,netmax,ngext,nodbmx,nordmx,noutpt,npet,npetmx,npt,
     $  nptmax,nsetmx,nstmax,nsts,nstsmx,nstsr,nttyo,uaqsln,ufixf,
     $  uphase,uspec,xbar,xbarlg,zklgmn,zklogl,zvclg0,zvclg1,
     $  zvec0,zvec1)
c
        delxi = max(dxsave,dlxmin)
        nord = min(nord,kord)
        jcut = jcut + 1
        if (jcut .le. 10) go to 130
        go to 180
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Cut the step size or order and try again or terminate.
c
      iexrt = 0
      jexrt = 0
      jscat = 0
      jscrt = 0
      qdump = .false.
      kly = 6
c
      if (delxi .gt. dlxmin) then
        write (noutpt,1440)
 1440   format(' --- Cutting the step size and trying again ---')
        ncut = ncut + 1
        delxi = 0.25*delxi
        delxi = max(delxi,dlxmin)
        go to 130
      endif
c
      if (nord .gt. 0) then
        nord = 0
        write (noutpt,1430)
 1430   format(' --- Cutting to order zero and trying again ---')
        go to 130
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Reaction path tracing has failed. Write a pickup file for the last
c     good point and stop.
c
  180 if (.not.qrpcfl .and. iopt(17).ge.0) then
        write (ux16,'(g12.5)') xi1
        call lejust(ux16)
        j2 = ilnobl(ux16)
        write (ux16a,'(g12.5)') xi0
        call lejust(ux16a)
        j3 = ilnobl(ux16a)
        write (noutpt,1450) ux16(1:j2),ux16a(1:j3)
        write (nttyo,1450)  ux16(1:j2),ux16a(1:j3)
 1450   format(/' * Error - (EQ6/path) Reaction path tracing has',
     $  ' failed',/7x,'at Xi= ',a,'. Will try to go back to',
     $  ' Xi= ',a,',',/7x,'the last point of reaction progress',
     $  ' at which the calculations',/7x,'were successful. If this',
     $  ' succeeds, will then stop and write a',/7x,'pickup file.',
     $  ' Try using this to restart the calculation.')
        qrpcfl = .true.
        iexrt = 0
        jexrt = 0
        jscat = 0
        jscrt = 0
        delxi = 0.
        qstopx = .true.
        go to 130
      else
        write (noutpt,1460)
        write (nttyo,1460)
 1460   format(/" * Error - (EQ6/path) Can't recover by going back",
     $  ' to the last point',/7x,'of reaction progress at which',
     $  ' calculations were successful.',/7x,"Can't write a pickup",
     $  ' file describing the system at that point.',/7x,'Try',
     $  ' re-running the problem, setting the maximum number',
     $  ' of steps',/7x,'(kstpmx) to the value which corresponds to',
     $  ' that point. Then try',/7x,'restarting the calculation',
     $  ' from that point.')
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The equilibrium calculation at the current point succeeded.
c
  200 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the activity of water.
c
      actwlg = actlg(narn1)
      actw = texp(actwlg)
      aw1 = actw
c
c     Get the weights (masses) of solvent, total dissolved solutes,
c     and aqueous solution, and get the aqeuous solution density.
c
      call gwdenp(adwipp,bdwipp,jcsort,mlmrra,mosp,mrmlra,
     $ mwtsp,narn1,narn2,nstmax,qdwipp,rhoc,rhowc,tdsgks,tdsglw,
     $ tdspkc,tdsplc,tempc,vosol,wfh2o,wftds,wkgwi,woh2o,
     $ wosol,wotds)
c
      if (qdwipp) then
        qrho = .true.
        rho = rhoc
      else
        qrho = .false.
        rho = 0.0
      endif
c
c     Compute pH, Eh, and pe-, all with reference to appropriate
c     pH scales. Also compute the pHCl.
c
      call gpheh(acflg,actlg,actwlg,adh,ah,ahmes,ahnbs,conc,
     $ eh,ehfac,ehmes,ehnbs,farad,fo2lg,fxi,iopg,mrmlra,nchlor,nhydr,
     $ nopgmx,noutpt,nstmax,nttyo,pch,pe,pemes,penbs,ph,phcl,phmes,
     $ phnbs,qphcl,qredox,qrho,xlke)
c
      wkgh2o = 1.e-3*woh2o
      wkgsol = 1.e-3*wosol
      ph1 = ph
      eh1 = eh
      fo2lg1 = fo2lg
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the new affinities and rates of the irreversible
c     reactions. Determine if the kinetic mode requires a corrector
c     step or a cut in the step size.
c
      if (nrct .le. 0) go to 300
c
c     Update the affinities of the irreversible reactions (afrc1).
c
      call raff(acflg,actlg,afcnst,affp,afrc1,bpx,cdrs,cgexj,
     $ ibpxmx,ibpxt,iern1,ietmax,iktmax,ixrn1,ixrn2,jcode,jern1,jern2,
     $ jetmax,jflag,jgext,jpflag,jsflag,jsol,ncmpr,ndrs,ndrsmx,ndrsr,
     $ nertmx,net,netmax,ngext,noutpt,nptmax,nrct,nrctmx,nrndex,
     $ nstmax,nttyo,nxridx,nxrtmx,nxtmax,rxbar,uphase,uspec,wfac,
     $ xbar,xbarlg,xgers,xlks)
c
c     Check the irreversible reactions for saturation.
c
      call rsatch(csts,egers,egexs,iern1,ietmax,iindx1,iktmax,
     $ iopt,ipndx1,jcode,jern1,jern2,jetmax,jgext,jpflag,jreac,kmax,
     $ km1,kmt,kx1,kxt,loph,losp,moph,morr,mosp,mrgers,mtb,mtb0,
     $ nbaspd,nbtmax,ncmpr,nern1,nern2,nert,nertmx,netmax,ngext,
     $ noptmx,noutpt,nptmax,nrct,nrctmx,nrk,nrndex,nstmax,nsts,nstsmx,
     $ nstsr,nttyo,nxridx,nxrt,nxrtmx,qreq,rxbar,tolxsf,uphase,ureac,
     $ uspec,xbar,xbarlg,zvclg1,zvec1)
c
c     Calculate rates of irreversible reactions at the new point by
c     evaluating the corresponding rate laws.
c
      call rtcalc(act,afrc1,cdac,csigma,eps100,fkrc,idirec,
     $ imchmx,imech,iodb,iopt,jcode,jreac,morr,morr0,mwtrc,ndac,
     $ ndact,ndctmx,nodbmx,noptmx,nord,noutpt,nrk,nrct,nrctmx,nsk,
     $ nstmax,nttyo,prcinf,prminf,qriinf,rirec1,rk,rreac1,rrelr1,
     $ rtcnst,rrxfi1,sfcar,sfcar0,ssfcar,udac,ureac)
c
      if (qshoot) go to 220
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(2) .le. 0) go to 230
c
      if (qriinf) then
c
c       Trap an infinite time step.
c
        if (delxi .le. dlxmin) then
c
c         Provisionally, take an infinite time step.
c
          deltim = prcinf
          time1 = prcinf
c
c         Make sure that the calculated time does not exceed any
c         any specified limits such as the maximum time just because
c         delxi is at the minimum value.
c
          call tivchk(deltim,delxi,qtvchk,time1,time0,timemx,
     $    tiplol,tiplot,tiprnl,tiprnt,tolxst)
c
          if (qtvchk) then
            qriinf = .false.
          else
            write (noutpt,1490)
            write (nttyo,1490)
 1490       format(/' --- Infinite time step ---',/)
          endif
c
          go to 230
        else
c
c         Go back and cut the step size before taking an infinite
c         time step.
c
          ier = 0
          go to 150
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ncorr .le. 0) then
c
c       If using the predictor function, get predicted values for the
c       inverse rate and the relative rates to compare with values
c       obtained from rate law expressions. If using a corrector
c       function, these variables are already available.
c
        call rtaylr(delxi,drer0,drir0,jreac,nord,nrct,nrctmx,
     $  nrd1mx,rirec0,rirecp,rrelr0,rrelrp)
c
c       Zero the "old" values of btrmax and dlrmax, which are used
c       to estimate the convergence rate for the ODE integrator.
c
        btrmxo = 0.
        dlrmxo = 0.
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(14).le.0 .and. .not.qhcon) then
c
c       If allowing both predictors, turn on the higher-order (stiff)
c       corrector if the simple predictor is currently on but not
c       rapidly converging.
c
        if (ncorr .ge. 4) then
          qhcon = .true.
          if (iodb(2) .gt. 1) then
            write (noutpt,1492)
            write (nttyo,1492)
 1492       format(' Switching to the higher-order (stiff) corrector')
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute residual functions for the ODE integrator.
c
      call betars(alphar,betar,btrfnc,btrmax,btrmxo,ibtrmx,
     $ nrct,nrct1,nrctmx,rirec1,rirecp,rrelr1,rrelrp)
c
      if (iodb(8) .ge. 2) then
        write (noutpt,1500)
 1500   format(/5x,'ODE residual functions:',
     $  //9x,'Reactant',22x,'betar',11x,'alphar',/)
        do nrc = 1,nrct
          write (noutpt,1510) ureac(nrc),betar(nrc),alphar(nrc)
 1510     format(7x,a24,5x,1pe11.4,5x,e11.4)
        enddo
        ux24 = 'Inverse rate'
        write (noutpt,1510) ux24,betar(nrct1),alphar(nrct1)
        write (noutpt,1512)
 1512   format(1x)
      endif
c
c     Test the accuracy of rate law integration. If the accuracy
c     is sufficient, qodeok will be returned with a value of .true.
c
      call tstari(afrc1,alphar,betar,delxi,iodb,modr,nodbmx,
     $ noutpt,nrct,nrctmx,nrct1,nsscmx,qodeok,rirec1,rirecp,rrelr1,
     $ rrelrp,sscrew,time1,tistrt,ureac)
c
      if (iodb(8) .ge. 1) then
        write(noutpt,1520) btrmax,btrfnc
 1520   format(5x,'btrmax= ',1pe10.3,5x,'btrfnc= ',e10.3)
      endif
c
c     Satisfying any of the following tests causes the current
c     results to be considered acceptable. Corrector iteration,
c     if begun, is terminated.
c
c     Force at a couple of corrector iterations.
c
      if (ncorr .lt. 2) then
        if (btrmax .gt. eps100) then
          qodeok = .false.
        endif
      endif
c
      if (qodeok) go to 220
      if (qmod .or. qbye) go to 220
      if (nord .le. 0) go to 220
c
c     Check for affinity sign flips.
c
      qaflip = .false.
      do j = 1,nrct
        if (afrc1(j).ge.0. .and. afrc0(j) .lt.0 .) then
          qaflip = .true.
          if (iodb(2) .gt. 1) then
            j2 = ilnobl(ureac(j))
            write (noutpt,1522) ureac(j)(1:j2)
            write (nttyo,1522) ureac(j)(1:j2)
 1522       format(' Affinity sign flip (',a,': - to 0/+)')
          endif
        elseif (afrc1(j).lt.0. .and. afrc0(j) .ge.0 .) then
          qaflip = .true.
          if (iodb(2) .gt. 1) then
            j2 = ilnobl(ureac(j))
            write (noutpt,1524) ureac(j)(1:j2)
            write (nttyo,1524) ureac(j)(1:j2)
 1524       format(' Affinity sign flip (',a,': 0/+ to -)')
          endif
        endif
      enddo
c
c     Satisfying any of the following tests causes ODE correction
c     to be terminated or foregone in favor of cutting the step
c     size.
c
      if (qaflip) go to 210
      if (iopt(14) .ge. 3) go to 210
      if (ncorr.ge.6 .and. delxi.gt.dlxmin) go to 210
      if (ncorr .ge. 8) go to 210
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set up to make an ODE corrector iteration step.
c
      ncorr = ncorr + 1
c
      if (ncorr .eq. 1) then
        qwhcfa = .true.
      elseif (delxi .le. 0.) then
        qwhcfa = .true.
      else
        adx = abs(dlxode - delxi)
        qwhcfa = (adx/delxi) .gt. eps100
      endif
c
      if (qwhcfa) then
c
c       This if block is normally entered only when ncorr = 1. If
c       the step size is cut during ODE correction, the calculations
c       in this block must be redone. When ncorr = 1, dlxode should
c       be set to prcinf (practical infinity). This value is set at
c       the start of a new step.
c
c       Save the current step size.
c
        dlxode = delxi
c
c       Calculate the dxsm10 vector from the step size (delxi) and the
c       dxsm00 vector. Note that dxsm10(1) = delXi(0,-1) = -delxi, but
c       dxsm10(2) = delXi(0,1), dxsm10(3) = delXi(0,2), etc. More
c       simply put, dxsm10 is the vector of cumulative step sizes to
c       the new point, but centered at point 0. In constrast, dxsm00
c       is the vector of cumulative step sizes to point 0, centered
c       at point 0.
c
        do j = 1,kordlm
          dxsm10(j + 1) = dxsm00(j)
        enddo
        dxsm10(1) = -delxi
c
c       Now calclate dxsm11, the vector of cumulative step sizes to
c       the new point, centered at the new point.
c
        do n = 1,kordlm
          j = kordlm - n + 1
          dxsm11(j + 1) = dxsm00(j) + delxi
        enddo
        dxsm11(1) = delxi
c
c       Compute the matrix (akmat1 = FC) used to calculate derivatives
c       for truncated Taylor's series, given corrector-based finite
c       differences.
c
c       Calling sequence substitutions:
c         akmat1 for akmat0
c         dxsm10 for dxsm00
c
        call gakmat(akmat1,dxsm10,nord,nrd1mx)
c
        if (qhcon) then
c
c         Calculate the w factor (whcfac) needed for higher-order
c         (stiff) ODE corrections. This depends only on the recent
c         step size history, including the current step size.
c
          call gwhcfa(akmat1,delxi,dxsm11,hhcvec,nord,nordmx,
     $    nrd1mx,whcfac,xhcvec)
        endif
      endif
c
      if (.not. qhcon) then
c
c       Use the simple corrector. The rates (rrelr1 and rirec1)
c       calculated from the governing rate equations  will be used
c       below as the rates at the new point when constructing new
c       finite differences for the corrector. Here just record the
c       correction in the rates at the new point.
c
        do j = 1,nrct1
          delvcr(j) = alphar(j)
        enddo
      endif
c
      if (qhcon) then
c
c       Use the higher-order (stiff) ODE integrator. Adjust rrelr1
c       and rirec1 using this scheme. These new values will be used
c       to compute a new finite-difference-based Taylor's-series-
c       formatted corrector function. Basically, the current adjustment
c       represents an application of the Newton-Raphson method.
c
c       Set up the right-hand-size vector.
c
        do j = 1,nrct1
          rhsvcr(j) = -alphar(j)
        enddo
c
c       Recompute the cdacb, ndacb, and ndactb arrays prior to
c       calculating the Jacobian matrix J[r].
c
cXX - Note: this can be changed by basis switching. Am probably
cXX   recalculating these arrays many times unnecessarily.
c
        call gndacb(cdac,cdacb,cdrs,eps100,imech,imchmx,
     $  jflag,nbasp,nbt,nbtmax,ndac,ndacb,ndact,ndactb,ndctmx,
     $  ndrs,ndrsmx,ndrsr,nrct,nrctmx,nstmax)
c
c       Compute the Jacobian matrix J[r] (armatr).
c
        call garmat(act,afrc1,aimatr,al10,armatr,cdac,cdacb,
     $  cdrs,csigma,csts,delvec,dlogxw,dvjdte,eact,eps100,fkrc,gmmatr,
     $  hact,iact,idirec,iindx1,iktmax,imchmx,imech,ipivot,jcode,jreac,
     $  jtemp,kbt,kdim,kmax,mmmatr,morr,mwtrc,nbasp,nbt,nbtmax,ncmpr,
     $  ndac,ndacb,ndact,ndctmx,ndrs,ndrsmx,ndrsr,noutpt,nptmax,nrct,
     $  nrctmx,nrct1,nrk,nrndex,nsk,nstmax,nsts,nstsmx,nstsr,nttkmx,
     $  nttyo,nxridx,nxrtmx,rirec1,rk,rkb,rreac1,rrelr1,rrxfi1,rtcnst,
     $  rxbar,sfcar,sgmatr,ssfcar,tempc,tempcb,tempk,ttk,ureac,whcfac,
     $  xi1,xlks,xxmatr,xymatr)
c
c       Solve for the correction vector (delvcr). If EQLIBU/msolvr.f
c       can't solve the matrix, it is because the matrix is either
c       zero (ier = 1) or non-zero, but computationally singular
c       (ier = 2).
c
c       Calling sequence substitutions:
c         armatr for aamatr
c         delvcr for delvec
c         grmatr for gmmatr
c         ipivtr for ipivot
c         nrct1  for kdim
c         nrct1  for kmax
c         rhsvcr for rhsvec
c
        qpr = .false.
        call msolvr(armatr,delvcr,grmatr,ier,ipivtr,nrct1,nrct1,
     $  noutpt,nttyo,qpr,rhsvcr)
c
        if (ier .le. 0) then
c
c         Correct the predicted values. Note that the results are
c         stored in rrelr1 and rirec1, not rrelrp and rirecp.
c
          do j = 1,nrct
            rrelr1(j) = rrelrp(j) + delvcr(j)
          enddo
          rirec1 = rirecp + delvcr(nrct1)
        else
c
c         Higher-order (stiff) ODE correction failed on this corrector
c         iteration. Drop back to the simple corrector.
c
          do j = 1,nrct1
            delvcr(j) = alphar(j)
          enddo
        endif
      endif
c
      if (iodb(8) .ge. 2) then
        ux24 = 'Simple'
        if (qhcon) ux24 = 'Higher-order (stiff)'
        j2 = ilnobl(ux24)
        write (noutpt,1530) ux24(1:j2)
 1530   format(/5x,a,' ODE corrector:',
     $  //9x,'Reactant',22x,'delvcr',7x,'Corrected r',/)
        do nrc = 1,nrct
          write (noutpt,1540) ureac(nrc),delvcr(nrc),rrelr1(nrc)
 1540     format(7x,a24,5x,1pe11.4,5x,e11.4)
        enddo
        ux24 = 'Inverse rate'
        write (noutpt,1540) ux24,delvcr(nrct1),rirec1
        write (noutpt,1512)
      endif
c
c     Find the max norm of the correction vector (delvcr). In the limit
c     of the solution, the (unrelaxed) correction term bounds the error
c     in the variables being calculated.
c
      idlrmx = iarmxn(delvcr,nrct1)
      dlrmax = 0.
      if (idlrmx .gt. 0) dlrmax = abs(delvcr(idlrmx))
c
c     Calculate the delvcr improvement function (dlrfnc).
c
      dlrfnc = 0.
      if (dlrmxo .gt. 0.) dlrfnc = (dlrmxo -dlrmax)/dlrmxo
c
      dlrmxo = dlrmax
c
      if (iodb(8) .ge. 2) then
        write(noutpt,1550) dlrmax,dlrfnc
 1550   format(9x,'dlrmax= ',1pe10.3,5x,'dlrfnc= ',e10.3)
      endif
c
c     Compute the new set of finite differences for a corrector
c     step.
c
      call corrfd(delxi,dxsm11,fdlim,fdre0,fdre1,fdri0,
     $ fdri1,fdrr0,fdrr1,iodb,iopt,jreac,nodbmx,noptmx,nord,
     $ nordmx,noutpt,npts,nrct,nrctmx,nrd1mx,rirec0,rirec1,
     $ rreac0,rreac1,rrelr0,rrelr1)
c
c     Now recompute the corresponding set of derivatives for the
c     original base point.
c
c     Calling sequence substitutions:
c       akmat1 for akmat0
c       fdri1 for fdri0
c       fdrr1 for fdrr0
c
      call rderiv(akmat1,drer0,drir0,fdri1,fdrr1,jreac,nord,
     $ nrct,nrctmx,nrd1mx)
c
c     Save the current corrected rate values. They will be the
c     predicted rates for the next corrector iteration. This
c     saves the effort of recalculating them from finite
c     differences.
c
      do j = 1,nrct
        rrelrp(j) = rrelr1(j)
      enddo
      rirecp = rirec1
c
c     Go back to the base point and try again with the current
c     corrector function.
c
      go to 140
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  210 continue
c
c     The corrector algorithm has failed to satsify the specified
c     criteria at the current order and step size.
c
      if (delxi .gt. dlxmin) then
c
c       Go back and cut the step size to satisfy the ODE corrector
c       tolerance. When the step size is cut to the minimum value
c       (dlxmin), the order will be reduced to zero and the
c       ODE integration will be considered sufficiently accurate.
c
        ncorr = 0
        ier = 0
        go to 150
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  220 continue
c
c     The AE/ODE algorithms have succeeded for the current step.
c     However, it may be necessary to go back and cut the step size
c     if an event has not been sufficiently accurately located.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make sure that the currently targeted time-based print point
c     has not been exceeded.
c
      tiprxx = min(tiprnl,tiprnt)
      if (((time1 - tiprxx)/tiprxx) .gt. tolxst) then
        if (delxi .gt. dlxmin) then
c
c         Go back and cut the step size.
c
          ncorr = 0
          ier = 0
          go to 150
        else
c
c         Set the current time equal to the targeted time-based print
c         point
c
          time1 = tiprxx
          deltim = tiprxx - time0
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make sure that the currently targeted time-based plot point
c     has not been exceeded.
c
      tiplxx = min(tiplol,tiplot)
      if (((time1 - tiplxx)/tiplxx) .gt. tolxst) then
        if (delxi .gt. dlxmin) then
c
c         Go back and cut the step size.
c
          ncorr = 0
          ier = 0
          go to 150
        else
c
c         Set the current time equal to the targeted time-based plot
c         point
c
          time1 = tiplxx
          deltim = tiplxx - time0
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make sure that the maximum time has not been exceeded.
c
      if (((time1 - timemx)/timemx) .gt. tolxst) then
        if (delxi .gt. dlxmin) then
c
c         Go back and cut the step size.
c
          ncorr = 0
          ier = 0
          go to 150
        else
c
c         Set the current time equal to the maximum time.
c
          time1 = timemx
          deltim = timemx - time0
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(14).le.0 .and. qhcon) then
c
c       If allowing both predictors, turn on the simple corrector
c       if the higher-order (stiff)  predictor is currently on
c       and rapidly converging.
c
        if (ncorr .le. 2) then
          qhcon = .false.
          if (iodb(2) .gt. 1) then
            write (noutpt,1560)
            write (nttyo,1560)
 1560       format(' Switching to the low-order (simple) corrector')
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  230 continue
c
c     Check for reactants that are newly exhausted and for formerly
c     exhausted reactants that are starting to precipitate. Note
c     that iexrt is the number of newly exhausted reactants, not
c     the total number of currently exhausted reactants. Here
c     jexrt is the number of formerly exhausted reactants that
c     have been reactivated.
c
      iexrt = 0
      jexrt = 0
      do nrc = 1,nrct
        if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
c
c         The reactant is active. See if there are any moles
c         remaining. If not, the reactant is exhausted if the
c         affinity is positive, the relative rate is positive,
c         or a finite amount was present at the previous point
c         of reaction progress (all three conditions should
c         agree).
c
          if (morr(nrc) .le. 0.) then
            if (afrc1(nrc).gt.eps100 .or. rrelr1(nrc).gt.eps100
     $        .or. morr0(nrc).gt.0.) then
              jreac(nrc) = 1
              iexrt = iexrt + 1
              iexr(iexrt) = nrc
            endif
          endif
        elseif (jreac(nrc) .eq. 1) then
c
c         The reactant is exhausted (has no moles remaining). See if
c         the affinity favors formation. If so, activate the reactant
c         so that this may proceed. Note that reactivation is not
c         appropriate for the case in which nrk(2,nrc) = 0. In that
c         case, precipitation is governed by partial equilibrium,
c         not a rate law.
c
          if (nrk(2,nrc) .ne. 0) then
            if (afrc0(nrc).gt.eps100 .and. afrc1(nrc).le.eps100) then
              jreac(nrc) = 0
              jexrt = jexrt + 1
              jexr(jexrt) = nrc
            endif
          endif
        endif
      enddo
c
c     Check to see that the point of reactivation of a formerly
c     exhausted reactant has been located sufficiently accurately.
c     A finite-difference based search does not guarantee this.
c
      jexrtx = 0
      do j = 1,jexrt
        nrc = jexr(j)
        if (abs(afrc1(nrc)) .gt. eps100) then
          jexrtx = jexrtx + 1
          if (iodb(1) .gt. 0) then
            j2 = ilnobl(ureac(nrc))
            write (noutpt,1570) ureac(nrc)(1:j2),afrc1(nrc)
 1570       format(/3x,'The affinity of formerly exhausted reactant',
     $      ' ',a,/5x,'is ',1pe11.4,' kcal. This is outside the',
     $      ' initially permitted',/5x,'range for reactivation to',
     $      ' allow precipitation according to',/5x,'a rate law.',/)
          endif
        endif
      enddo
c
      if (jexrtx .gt. 0) then
c
c       Determine whether or not to reduce delxi and go back and try
c       again to better locate an event in question.
c
        call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)
        if (qadjdx) then
          ier = 0
          go to 150
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for reactants that have affinities or reaction rates which
c     are about to change sign. Here jscat is the number of reactants
c     with affinities about to undergo a sign change, and jscrt is
c     the number with rates about to undergo a sign change. Note that
c     a reactant with an affinity undergoing a sign change may not have
c     reaction rate undergoing a sign change, as the rate may be
c     truncated to zero for a given sign of the affinity. Also, even
c     if theoretically both the affinity and the rate should undergo
c     a sign change simultaneously, the testing tolerances may be such
c     that a sign change for only one or the other condition is caught
c     and reported.
c
      jscat = 0
      do nrc = 1,nrct
        if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
          if (nrk(2,nrc) .ne. 0) then
            if (afrc0(nrc).gt.tolsar .and. afrc1(nrc).le.tolsar) then
              jscat = jscat + 1
              jsca(jscat) = nrc
            elseif (afrc0(nrc).lt.-tolsar .and. afrc1(nrc).ge.-tolsar)
     $        then
              jscat = jscat + 1
              jsca(jscat) = nrc
            endif
          endif
        endif
      enddo
c
      jscrt = 0
      do nrc = 1,nrct
        if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
          if (nrk(2,nrc) .ne. 0) then
            if (rrelr0(nrc).gt.tolsrr .and. rrelr1(nrc).le.tolsrr) then
              jscrt = jscrt + 1
              jscr(jscrt) = nrc
            elseif (rrelr0(nrc).lt.-tolsrr .and. rrelr1(nrc).ge.-tolsrr)
     $        then
              jscrt = jscrt + 1
              jscr(jscrt) = nrc
            endif
          endif
        endif
      enddo
c
c     Check for oversteps. The finite-difference based searches do
c     not guarantee sufficient accuracy.
c
      jscatx = 0
      do j = 1,jscat
        nrc = jsca(j)
        if (abs(afrc1(nrc)) .gt. tolsar) then
          jscatx = jscatx + 1
          if (iodb(1) .gt. 0) then
            j2 = ilnobl(ureac(nrc))
            write (noutpt,1580) ureac(nrc)(1:j2),afrc1(nrc)
 1580       format(/3x,'The affinity of reactant ',a,' is ',1pe11.4,
     $      ' kcal. This is',/5x,'too much of an overstep past the',
     $      ' point at which the affinity of this',/5x,'reactant',
     $      ' changes sign.',/)
          endif
        endif
      enddo
c
      jscrtx = 0
      do j = 1,jscrt
        nrc = jscr(j)
        if (abs(rrelr1(nrc)) .gt. tolsrr) then
          jscrtx = jscrtx + 1
          if (iodb(1) .gt. 0) then
            j2 = ilnobl(ureac(nrc))
            write (noutpt,1590) ureac(nrc)(1:j2),rrelr1(nrc)
 1590       format(/3x,'The relative rate of reactant ',a,' is ',
     $      1pe11.4,'. This is',/5x,'too much of an overstep past the',
     $      ' point at which the relative rate of',/5x,'this',
     $      ' reactant changes sign.',/)
          endif
        endif
      enddo
c
      if (jscatx.gt.0 .or. jscrtx.gt.0) then
c
c       Determine whether or not to reduce delxi and go back and try
c       again to better locate an event in question.
c
        call dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)
        if (qadjdx) then
          ier = 0
          go to 150
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  300 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for oversteps of the pH, Eh, log fO2, or aw with respect to
c     currently defined print point values.
c
      call cophpr(actw,aw0prn,aw1prn,delxi,dlxmin,eh,eh0prn,
     $ eh1prn,fo2lg,iodb,nodbmx,noutpt,o20prn,o21prn,ph,ph0prn,
     $ ph1prn,qadjdx,qredox,tolxsu)
      if (qadjdx) then
        ier = 0
        go to 150
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for oversteps of the pH, Eh, log fO2, or aw with respect to
c     currently defined plot point values.
c
      call cophpl(actw,aw0plo,aw1plo,delxi,dlxmin,eh,eh0plo,
     $ eh1plo,fo2lg,iodb,nodbmx,noutpt,o20plo,o21plo,ph,ph0plo,
     $ ph1plo,qadjdx,qredox,tolxsu)
      if (qadjdx) then
        ier = 0
        go to 150
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for oversteps of the pH, Eh, log fO2, or aw with respect to
c     requested minimum and maximum values.
c
      call cophlm(actw,awmax,awmin,delxi,dlxmin,eh,ehmax,
     $ ehmin,fo2lg,iodb,nodbmx,noutpt,o2max,o2min,ph,phmax,phmin,
     $ qadjdx,qredox,tolxsu)
      if (qadjdx) then
        ier = 0
        go to 150
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Reset the iemop, emop, iemos, emos, etc., arrays.
c
      if (qmod) then
c
c       The ES phase assemblage has changed. Re-set the index arrays
c       associated with finite-difference description of the numbers of
c       mole of phases and species in the Equilibrium System (ES).
c
        call iiemop(iemop,iemos,iindx1,ipndx1,jsflag,kdim,kmax,
     $  ncmpe,ncmpr,noutpt,npet,npetmx,npt,nptmax,nset,nsetmx,nstmax,
     $  nttyo,uaqsln,uspec,uphase)
      endif
c
c     Update the numbers of moles variables for phases present in
c     the ES.
c
      do npe = 1,npet
        np = iemop(npe)
        emop(npe) = moph(np)
        nr1 = ncmpe(1,npe)
        nr2 = ncmpe(2,npe)
        do nse = nr1,nr2
          ns = iemos(nse)
          emos(nse) = mosp(ns)
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Recalculate the total affinity (aft1).
c
      call caft1(afrc1,aft1,nrct,nrctmx,rrelr1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write entertainment for the user, showing the progress
c     of the current run.
c
      call wrentu(actw,eh,fo2lg,iopg,iopt,kstep,nopgmx,noptmx,
     $ nttyo,ph,qredox,time1,xi1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      qprntx = (xi1 + eps100) .ge. xiprnt
      qprntt = (time1 + eps100) .ge. tiprnt
      qprnlx = (xi1 + eps100) .ge. xiprnl
      qprnlt = (time1 + eps100) .ge. tiprnl
c
      qprph0 = abs(ph - ph0prn) .le. tolxsu
      qprph1 = abs(ph - ph1prn) .le. tolxsu
      qpreh0 = abs(eh - eh0prn) .le. tolxsu
      qpreh1 = abs(eh - eh1prn) .le. tolxsu
      qpro20 = abs(fo2lg - o20prn) .le. tolxsu
      qpro21 = abs(fo2lg - o21prn) .le. tolxsu
      qpraw0 = abs(actw - aw0prn) .le. tolxsu
      qpraw1 = abs(actw - aw1prn) .le. tolxsu
c
      qzprnt = qprntx .or. qprntt .or. qprnlx .or. qprnlt .or.
     $ qprph0 .or. qprph1 .or. qpreh0 .or. qpreh1 .or.
     $ qpro20 .or. qpro21 .or. qpraw0 .or. qpraw1 .or.
     $ kstppr.ge.ksppmx
c
      qplotx = (delxi + eps100) .ge. (xiplot - xi0)
      qplott = (deltim + eps100) .ge. (tiplot - time0)
      qplolx = (delxi + eps100) .ge. (xiplol - xi0)
      qplolt = (deltim + eps100) .ge. (tiplol - time0)
c
      qplph0 = abs(ph - ph0plo) .le. tolxsu
      qplph1 = abs(ph - ph1plo) .le. tolxsu
      qpleh0 = abs(eh - eh0plo) .le. tolxsu
      qpleh1 = abs(eh - eh1plo) .le. tolxsu
      qplo20 = abs(fo2lg - o20plo) .le. tolxsu
      qplo21 = abs(fo2lg - o21plo) .le. tolxsu
      qplaw0 = abs(actw - aw0plo) .le. tolxsu
      qplaw1 = abs(actw - aw1plo) .le. tolxsu
c
      qzplot = qplotx .or. qplott .or. qplolx .or. qplolt .or.
     $ qplph0 .or. qplph1 .or. qpleh0 .or. qpleh1 .or.
     $ qplo20 .or. qplo21 .or. qplaw0 .or. qplaw1 .or.
     $ kstppl.ge.ksplmx
c
      qzdump = (xi1 + eps100) .ge. xidump
c
      if (qzdump) qdump = .true.
      if (qrapch .and. .not.qmod .and. .not.qbye) kly = max(kly,6)
c
c     Check to see which variable is changing most rapidly.
c
      if (qmod .or. qbye) then
        write (noutpt,1900) kstep,iter,ncorr
 1900   format(' Steps completed= ',i5,', iter= ',i3,', ncorr= ',i1)
      else
        inmax = 0
        kzmax = 0
        adlzlg = 0.
        do kcol = 1,kdim
          adel = abs(zvclg1(kcol) - zvclg0(kcol))
          if (adel .gt. adlzlg) then
            inmax = iindx1(kcol)
            kzmax = kcol
            adlzlg = adel
          endif
        enddo
        adlzlg = 0.
        if (kzmax .le. 0) then
          write (noutpt,1900) kstep,iter,ncorr
        else
          adlzlg = zvclg1(kzmax)
          ustr24 = 'Error'
          if (kzmax .le. kbt) then
            ns = nbasp(inmax)
            ustr24 = uspec(ns)
          elseif (kzmax .le. kxt) then
            ustr24 = uspec(inmax)
          endif
          j2 = ilnobl(ustr24)
          write (noutpt,1910) kstep,iter,ncorr,ustr24(1:j2),adlzlg
 1910     format(' Steps completed= ',i5,', iter= ',i3,', ncorr= ',i1,
     $     /' Most rapidly changing is zvclg1(',a,')= ',f11.4)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check gases for which the fugacity is supposed to be fixed.
c
      iwnffg = 0
      do n = 1,nffg
        ng = jffg(n)
        xlf = xlkffg(n)
        xf = texp(xlf)
        dx = fugalg(ng) - xlkffg(n)
        if (abs(dx) .gt. toldl) then
          if (iwnffg .le. 0) nwnffg = nwnffg + 1
          iwnffg = iwnffg + 1
          if (nwnffg .le. nlwffg) then
            j2 = ilnobl(uffg(n))
            write (noutpt,2000) uffg(n)(1:j2),fugac(ng),xf,xi1
            write (nttyo,2000) uffg(n)(1:j2),fugac(ng),xf,xi1
 2000       format(/' * Warning - (EQ6/path) The fugacity of ',a,
     $      /7x,'is now ',1pg12.5,' bars. It is supposed to be fixed',
     $      /7x,'at ',e12.5,' bars. At the present value of reaction',
     $      /7x,'progress (Xi= ',e12.5,'), there is no longer a',
     $      /7x,'sufficient mass of the gas component in the system',
     $      /7x,'to "saturate" it at the desired fugacity value. If',
     $      /7x,'you restart the run, you can add such mass using the',
     $      /7x,'moffg parameter on the input file. Add only about',
     $      /7x,'0.5 to 1.0 mole at a time. You can also increase',
     $      /7x,'the amount of gas component present as the run',
     $      /7x,'progresses by specifying the gas as a reactant.')
          endif
        endif
      enddo
      if (nwnffg .eq. nlwffg) then
        write (noutpt,2010)
        write (nttyo,2010)
 2010   format(/' * Warning - (EQ6/path) No more warnings will be',
     $  /7x,'issued regarding gas fugacities not being at desired',
     $  /7x,'fixed values.')
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Test for very little remaining solvent water.
c
      if (wkgh2o .lt. 1.e-10) then
        if (zvec1(1) .lt. zvec0(1)) then
          write (noutpt,2020) wkgh2o,wkgsol
          write (nttyo,2020) wkgh2o,wkgsol
 2020     format(/' * Note - (EQ6/path) The amount of remaining',
     $    ' solvent water',/7x,'is only ',1pe12.5,' kg and is',
     $    ' decreasing with reaction progress.',/7x,'The amount of',
     $    ' remaining aqueous solution is ',e12.5,' kg.',
     $    /7x,'This code is not designed to deal with fully dry',
     $    ' systems. This',/7x,'run will now be terminated.')
          qvlsow = .true.
        endif
      endif
c
      if (wkgh2o .lt. 1.e-8) then
        if (iwdh2o .ge. 0) then
          nwdh2o = nwdh2o + 1
          if (nwdh2o .le. 8) then
            if (zvec1(1) .lt. zvec0(1)) then
              write (noutpt,2030) wkgh2o,wkgsol
              write (nttyo,2030) wkgh2o,wkgsol
 2030         format(/' * Warning - (EQ6/path) The amount of remaining',
     $        ' solvent water',/7x,'is only ',1pe12.5,' kg and is',
     $        ' decreasing with reaction progress.',/7x,'The amount of',
     $        ' remaining aqueous solution is ',e12.5,' kg.',
     $        /7x,'This code is not designed to deal with fully dry',
     $        ' systems. This',/7x,"run can't continue much farther.")
            else
              write (noutpt,2040) wkgh2o,wkgsol
              write (nttyo,2040) wkgh2o,wkgsol
 2040         format(/' * Warning - (EQ6/path) The amount of remaining',
     $        ' solvent water',/7x,'is only ',1pe12.5,' kg. The amount',
     $        ' of remaining aqueous solution',/7x,'is only ',e12.5,
     $        ' kg. This code is not designed to deal with',/7x,
     $        ' fully dry systems. This run is near the limit of the',
     $        " code's capability.")
            endif
            iwdh2o = -10
          endif
        else
          iwdh2o = iwdh2o + 1
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Test for very high ionic strength.
c
      if (fxi .gt. 100.) then
        write (ux16,'(f10.3)') fxi
        call lejust(ux16)
        j2 = ilnobl(ux16)
        write (noutpt,2050) ux16(1:j2)
        write (nttyo,2050) ux16(1:j2)
 2050   format(/' * Note - (EQ6/path) The ionic strength is ',a,
     $  ' molal.',/7x,'This run will now be terminated.')
        qvhfxi = .true.
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Do automatic basis switching. Optimize the basis set so that
c     the basis species for each mass balance tends to be the species
c     which dominates that mass balance.
c
      qabswx = .false.
      if (iopt(12).gt.0 .and. kstpab.ge.20) then
c
c       Automatic basis switching after solving at the current point
c       of reaction progress is turned on, and this is the 20th step
c       since the last time the need for this was checked.
c
        call absswb(adhfs,adhfsx,advfs,advfsx,avcnst,axhfs,
     $  axhfsx,axlks,axlksx,axvfs,axvfsx,cdrs,cdrtw,cdrsx,cdrw,
     $  csts,dhfs,dvfs,eps100,ibswx,iindx1,iodb,ipch,ipchmx,ipcv,
     $  ipcvmx,jcsort,jflag,jsflag,kbt,kmax,mosp,mtb,narn1,narn2,
     $  narxmx,narxt,nbasp,nbaspd,nbaspx,nbt,nbtmax,nbw,ndrs,ndrsmx,
     $  ndrsr,ndrsrx,ndrsx,nelect,nhydr,nodbmx,no2gaq,noutpt,nst,
     $  nstmax,nsts,nstsmx,nstsr,nswtch,ntpr,ntprmx,nttyo,presg,
     $  press,qbassw,qbswx,tempc,uspec,uzvec1,weight,xhfs,xvfs,xlks)
c
        if (nswtch .gt. 0) then
c
c         Set a flag noting that automatic basis switching has just
c         been done.
c
          qabswx = .true.
c
          write (noutpt,2100) nswtch
 2100     format(/' ',i2,' basis switches were executed automatically',
     $    ' after solving',/'at the current value of reaction',
     $    ' progress.')
        endif
c
c       Reset the step counter which contains the number of steps since
c       the last time a check was made of the need for automatic basis
c       switching of this kind.
c
        kstpab = 0
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Save the order parameters.
c
      nordsv = nord
      kordsv = kord
c
c     Check to see if it is necessary to drop the order of the
c     Taylor's series to zero because a point has been reached
c     at which the derivatives of the variables being described
c     by finite differences are not continuous.
c
      if (qmod .or. qbye .or. qabswx .or. iexrt.gt.0 .or.
     $  jexrt.gt.0 .or. jscat.gt.0 .or. jscrt.gt.0 .or.
     $  qtrch .or. qreq) then
c
c       Drop the order of the Taylor's series to zero. Zero all
c       finite-difference and derivative data.
c
        if (iopt(1) .eq. 2) then
c
c         Set flags for dump and to control special print for
c         the case of the fluid-centered flow-through open system.
c
          qdump = .true.
          qdmpr1 = .true.
          qdmpr2 = .true.
        endif
c
        xilim = prcinf
        delxi = dlxmx0
        npts = 0
c
        call initaz(dxsm00,nrd1mx)
        call initaz(dxsm10,nrd1mx)
        call initaz(dxsm11,nrd1mx)
c
        nmax = nrd1mx*kmax
        call initaz(dzvc0,nmax)
        call initaz(dzvc0s,nmax)
        call initaz(fdzv0,nmax)
        call initaz(fdzvm1,nmax)
c
        call initaz(fdri0,nrd1mx)
        call initaz(fdrim1,nrd1mx)
        call initaz(drir0,nrd1mx)
        call initaz(drir0s,nrd1mx)
c
        nmax = nordmx*nrctmx
        call initaz(fdar0,nmax)
        call initaz(fdarm1,nmax)
        call initaz(dafrc0,nmax)
c
        nmax = nrd1mx*nrctmx
        call initaz(fdrr0,nmax)
        call initaz(fdrrm1,nmax)
        call initaz(drer0,nmax)
        call initaz(drer0s,nmax)
      endif
c
      npts = npts + 1
      kord = npts - 1 - ndelay
      kord = min(kord,kordlm)
      kord = max(0,kord)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check to see that the reaction path is not constant.
c
      qconst = .false.
      if (qcntmp .and. qcnpre) then
        if (nrct .gt. 0) then
          do nrc = 1,nrct
            if (rrelr1(nrc) .ne. 0.) go to 330
          enddo
        endif
        qconst = .true.
      endif
  330 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(2) .le. 0) then
c
c       If the temperature and pressure are constant, check to see if
c       the total affinity is repeatedly nearly zero.
c
        if (qcntmp .and. qcnpre) then
          if (aft1 .le. tolsar) then
            naft1 = naft1 + 1
            qaft1 = naft1 .ge. 20
          else
            naft1 = 0
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(2) .gt. 0) then
        if (qshoot) then
          jscat = 0
          jscrt = 0
        else
c
c         Check to see if the code is just making small oscillations
c         about the final equilibrium point.
c
          if (nord .le. 0) then
            if (aft1.gt.aft0 .and. aftm1.gt.aft0) then
              kaft1 = kaft1 + 1
              if (aft1 .le. 0.001) then
                jscat = 0
                jscrt = 0
              endif
            elseif (aft1.lt.aft0 .and. aftm1.lt.aft0) then
              kaft1 = kaft1 + 1
              if (aft1 .le. 0.001) then
                jscat = 0
                jscrt = 0
              endif
            else
              kaft1 = 0
            endif
          else
            kaft1 = 0
          endif
          if (kaft1 .ge. 6) then
            if (aft1 .le. 0.001) then
c
c             The total affinity is oscillating about some finite
c             value close to zero.
c
              time1 = min(timemx,tiprnt,tiprnl,tiplot,tiplol)
              deltim = time1 - time0
              qshoot = .true.
              qprntt = (time1 + eps100) .ge. tiprnt
              qprnlt = (time1 + eps100) .ge. tiprnl
              qzprnt = qzprnt .or. qprntt .or. qprnlt
              qplott = (time1 + eps100) .ge. tiplot
              qplolt = (time1 + eps100) .ge. tiplol
              qzplot = qzplot .or. qplott .or. qplolt
            endif
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     If doing a fluid-centered flow-through system mode run, check
c     to see if a "short" print description is required after a transfer
c     to the PRS.
c
      if (iopt(1).eq.2 .and. nord.gt.0)
     $ qftpr2 = .not.qdmpr1 .and. qdmpr2
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Determine if a print and/or plot of the system should
c     be made at the current point of reaction progress.
c
      if (xi1 .ge. ximax) go to 340
      if (iopt(2) .gt. 0) then
        if (time1 .ge. ((1. - tolxst)*timemx)) go to 340
      endif
      if (ph .le. (phmin + tolxsu)) go to 340
      if (ph .ge. (phmax - tolxsu)) go to 340
      if (qredox) then
        if (eh .le. (ehmin + tolxsu)) go to 340
        if (eh .ge. (ehmax - tolxsu)) go to 340
        if (fo2lg .le. (o2min + tolxsu)) go to 340
        if (fo2lg .ge. (o2max - tolxsu)) go to 340
      endif
      if (actw .le. (awmin + tolxsu)) go to 340
      if (actw .ge. (awmax - tolxsu)) go to 340
      if (kstep .ge. kstpmx) go to 340
c
      if (qstopx) go to 340
      if (iexrt.gt.0 .or. jexrt.gt.0) go to 340
      if (jscat.gt.0 .or. jscrt.gt.0) go to 340
      if (qconst .or. qmod .or. qbye .or. qreq .or. qtrch) go to 340
      if (qaft1) go to 340
      if (qvlsow .or. qvhfxi) go to 340
c
      if (qzprnt .or. qftpr2 .or. qdump) go to 350
c
      if (qzplot) go to 360
c
c     Set upper limits on the scale factor for the next step. Here
c     kly is a delay factor that comes into play when the step size
c     has been cut recently in connection with a search.
c
      scalim = 10.
      if (kly .gt. 0) scalim = 2.0
      go to 100
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  340 kly = 0
      if (qftpr2) qdmpr2 = .false.
      qftpr2 = .false.
c
  350 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute apparent "whole-phase" equivalent fractions and mole
c     fractions of the exchange ions present in generic ion exchanger
c     phases. Cations and anions are treated separately in these
c     calculations.
c
      call gegexw(cegexs,egexpc,egexpa,egexw,iern1,iern2,
     $ ietmax,jern1,jetmax,jgext,kern1,kern2,ketmax,kgexsa,moph,
     $ mosp,netmax,ngexsa,ngext,noutpt,nptmax,nstmax,nttyo,
     $ xgexw,zchar)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute various secondary parameters at the current point of
c     reaction progress.
c
      call cdappl(acflg,acfw,acfwlg,actlg,actw,actwlg,adwipp,
     $ afcnst,affpd,affsd,ah,ahrc,alk,alk1,alk2,alki,atwt,bdwipp,
     $ cdrsd,cess,conc,csts,ctb,cteaq,dvoso,dwoso,eh,ehfac,ehrc,
     $ eps100,farad,fdpe0,fdse0,fjest,fo2lg,fo2lrc,fxist,iaqsln,
     $ iemop0,iemos0,iern1,iern2,ifrn1,ifrn2,ilrn1,ilrn2,imrn1,
     $ imrn2,iopt,ixrn1,ixrn2,jcode,jcsort,jern1,jflag,jflagd,
     $ jgext,jpflag,jsflag,jssort,modr,moph,mophg,mophj,mopht,morr,
     $ mosp,mospg,mospj,mospt,mprph,mprsp,mrgers,mrmlra,mtb,
     $ mtbaq,mte,mteaq,mwtges,mwtrc,mwtsp,narn1,narn2,nat,nbasp,
     $ nbaspd,nbt,nchlor,ncmpe0,ncmpr,nct,ndrsd,ndrsrd,nelect,
     $ nern1,nern2,nert,ness,nessr,net,nfrn1,nfrn2,ngext,ngrn1,
     $ ngrn2,ngt,nhydr,nlrn1,nlrn2,nlt,nmrn1,nmrn2,nmrt,nmt,no2gaq,
     $ npchk,npet,npet0,npt,npts,nrct,nrndex,nst,nsts,nstsr,ntf1,
     $ ntf1t,ntf2,ntf2t,nxridx,nxrn1,nxrn2,nxrt,nxt,osc,oscst,
     $ omega,pe,perc,ph,phmes,ppmwb,ppmwe,qriinf,qxknph,rreacn,
     $ rreac1,rxbar,sfcar,sidrsp,sidrph,sigmam,sigmst,tdays,tempc,
     $ tf1,tf2,thours,time1,tmins,tyears,uphase,uspec,vodrt,voph,
     $ vophg,vophj,vopht,vosoct,vosp,vospg,vospj,vospt,vosp0,
     $ vreac,wfh2o,wkgwi,wodr,wodrt,woph,wophg,wophj,wopht,
     $ worr,worrt,wosoct,wosp,wospg,wospj,wospt,xbar,xlke,
     $ xlksd,zchcu6,zchsq2)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print a description of the system at the current point of
c     reaction progress.
c
      call scripz(abar,acflg,acfw,acfwlg,actlg,actw,actwlg,
     $ affpd,affsd,afrc1,aft1,ah,ahmes,ahnbs,ahrc,alki,alk1,
     $ alk2,awmax,awmin,a3bar,cbsr,cdrsd,cegexs,cesr,conc,conclg,csts,
     $ ctb,cteaq,dvoso,dwoso,egers,egexjc,egexjf,egexpa,egexpc,egexs,
     $ egexw,eh,ehmax,ehmes,ehmin,ehnbs,ehrc,elecsr,electr,fje,fjest,
     $ fo2,fo2lg,fo2lrc,fugac,fugalg,fxi,fxist,iaqsln,iemop,iemop0,
     $ iemos,iemos0,iern1,iern2,iexr,iexrt,ifrn1,ifrn2,ilrn1,ilrn2,
     $ imech,imrn1,imrn2,iopg,iopr,iopt,ipndx1,ixrn1,ixrn2,jcode,jcsort,
     $ jern1,jern2,jexr,jexrt,jflag,jflagd,jflgi,jgext,jgsort,jpflag,
     $ jreac,jsca,jscat,jscr,jscrt,jsflag,jsol,jssort,kbt,kern1,kern2,
     $ kgexsa,km1,kmt,kx1,kxt,kstep,kstpmx,loph,losp,mlmrra,modr,
     $ moph,mophg,mophj,mopht,morr,mosp,mospg,mospj,mospt,mprph,
     $ mprsp,mrgers,mrmlra,mwtrc,mwtsp,narn1,narn2,nat,nbasp,nbaspd,
     $ nbt,ncmpe,ncmpe0,ncmpr,nct,ndrsd,ndrsrd,nelect,nern1,nern2,
     $ nert,net,nfrn1,nfrn2,ngext,ngexsa,ngrn1,ngrn2,ngt,nhydr,nhydx,
     $ nlrn1,nlrn2,nlt,nmrn1,nmrn2,nmrt,nmt,noutpt,no2gaq,npet,npet0,
     $ npt,npts,nrct,nrdxsp,nrk,nrndex,nst,nsts,nstsr,ntf1t,ntf2t,
     $ nxridx,nxrn1,nxrn2,nxrt,nxt,osc,oscst,omega,o2max,o2min,pch,pe,
     $ pemes,penbs,perc,ph,phcl,phmax,phmes,phmin,phnbs,ppmwe,presg,
     $ press,qaft1,qftpr2,qmod,qphcl,qredox,qrho,qriinf,qstopx,qvhfxi,
     $ qvlsow,qzprnt,rho,rhoc,rhowc,rk,rreacn,rreac1,rrelr1,rxbar,
     $ sfcar,sidrph,sidrsp,sigmst,sigmam,ssfcar,tdays,tdsglw,tdspkc,
     $ tdsplc,tempc,thours,time1,timemx,tmins,tolsat,tolxsf,tolxst,
     $ tolxsu,tyears,uelem,ugermo,ugexj,ugexmo,uphase,ureac,uspec,
     $ uxtype,vodrt,voph,vophg,vophj,vopht,vosoct,vosol,vosp,vospg,
     $ vospj,vospt,vreac,wfh2o,wftds,wkgwi,woh2o,wodr,wodrt,woph,
     $ wophg,wophj,wopht,worr,worrt,wosoct,wosol,wosp,wospg,wospj,
     $ wospt,wotds,xbar,xbarlg,xbarw,xbrwlg,xgers,xgexw,xi1,xidump,
     $ ximax,xistsv,xirct,zchar)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write results at the current point of reaction progress on the
c     tabx file. This file is used to create the plot file tab.
c
 355  continue
      if (ximax.gt.xistsv .and. iopt(18).ge.0) then
        if (qtatxt) then
c
c         The TAB file is an ordinary text file.
c
          call wrtabx(actlg,afrc1,aft1,alk,cteaq,dvoso,dwoso,
     $    eh,fo2lg,iindx1,iktmax,iopt,ipndx1,kmax,km1,kmt,kstep,kx1,
     $    kxt,loph,ncmpr,modr,mopht,narn1,mosp,nct,nctmax,noptmx,
     $    nptmax,nrct,nrctmx,nstmax,ntabx,ntidmx,ntitl2,ntitld,ntitmx,
     $    nxtmax,pe,ph,ppmwe,prcinf,press,prminf,qbye,qmod,qriinf,
     $    tempc,time1,uelem,uphase,uplatm,ureac,uspec,usteq6,utitl2,
     $    utitld,uveeq6,vodrt,vosoct,wodrt,woh2o,wosoct,xbar,xi1)
        else
c
c         The TAB file is a .csv file.
c
          call wrtabc(acflg,actlg,actw,afrc1,aft1,alk,conclg,
     $    cteaq,ctb,dvoso,dwoso,eh,fje,fo2lg,fugac,fxi,iktmax,iopt,
     $    jflag,jsflag,kmax,kstep,kx1,kxt,mrmlra,modr,mosp,mospt,
     $    moph,mopht,narn1,narn2,nbasp,nbaspd,nbt,nbtmax,ncmpr,nct,
     $    nctmax,nelect,ngrn1,ngrn2,ngtmax,nhydr,nhydx,nllnmx,
     $    no2gaq,noptmx,noutpt,npt,nptmax,nrct,nrctmx,nstmax,ntabx,
     $    ntidmx,ntitl2,ntitld,ntitmx,nttyo,nxrn1,nxrn2,nxtmax,pe,ph,
     $    phmes,ppmwb,ppmwe,prcinf,press,prminf,qrho,qriinf,rho,rhowc,
     $    sidrph,sigmam,tdsgks,tdsglw,tempc,time1,uelem,ulinex,
     $    uphase,uplatm,ureac,uspec,usteq6,utitl2,utitld,uveeq6,
     $    vodrt,vosoct,wkgh2o,wodrt,wosoct,xbar,xbarlg,xi1)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(16) .ge. 0) then
c
c       Prepare to write results at the current point to the
c       backup file.
c
        call setpk6(actwlg,awmax,awmaxi,awmin,awmini,eh,ehmax,
     $  ehmaxi,ehmin,ehmini,fo2lg,iindx1,jflag,jflgi,kbt,kdim,kmax,
     $  kprs,mprph,mprphi,mprsp,mprspi,mtb,mtbi,mtbaq,mtbaqi,nbasp,
     $  nbaspd,nbaspi,nbti,nbtmax,ncmpr,nobswt,noutpt,nprpmx,nprpti,
     $  nprsmx,nprsti,npt,nptmax,nttyo,nstmax,o2max,o2maxi,o2min,
     $  o2mini,ph,phmax,phmaxi,phmin,phmini,prcinf,press,pressi,
     $  tempc,tempci,time1,timemx,timmxi,tistti,ubmtbi,uobsw,uphase,
     $  uprphi,uprspi,uspec,uzveci,uzvec1,xi1,ximax,ximaxi,xistti,
     $  zvclgi,zvclg1)
c
c       Rewind the current backup file (BAKUPA or BAKUPB).
c
        if (iopt(16) .eq. 0) rewind (nbkupn)
c
c       Write results at the current point to the backup file.
c
        if (upkfor(1:1) .eq. 'W') then
c
c         Compact (W) format.
c
c         Calling sequence substitutions:
c           nbkupn for newin
c
          call wr6pkw(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,
     $    dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,
     $    dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,
     $    dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,iact,
     $    ibsrti,iesrti,ietmax,igerti,iktmax,imchmx,imech,iodb,iopg,
     $    iopr,iopt,itermx,ixrti,jcode,jetmax,jflgi,jgerti,jgext,
     $    jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,
     $    kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,
     $    mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,
     $    nertmx,net,netmax,nbkupn,nffg,nffgmx,ngexrt,nobswt,nodbmx,
     $    nopgmx,noprmx,noptmx,nordmx,nprpmx,nprpti,nprsmx,nprsti,
     $    nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,
     $    ntitmx,ntrymx,nttkmx,nxmdmx,nxmod,nxopex,nxopmx,nxopt,
     $    nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,pressb,
     $    pressi,ptk,qgexsh,rkb,rxbari,sfcar,ssfcar,tempcb,tempci,
     $    tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,trkb,ttk,
     $    ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,ugersi,
     $    ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,ureac,
     $    usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,uxopt,
     $    uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,
     $    xlkmod,xvfgex,zgexj,zvclgi)
        else
c
c         Menu-style (D) format.
c
c         Calling sequence substitutions:
c           nbkupn for newin
c
          call wr6pkd(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,
     $    dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,
     $    dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,
     $    dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,iact,
     $    ibsrti,iesrti,ietmax,igerti,iktmax,imchmx,imech,iodb,iopg,
     $    iopr,iopt,itermx,ixrti,jcode,jetmax,jflgi,jgerti,jgext,
     $    jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,
     $    kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,
     $    mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,
     $    nertmx,net,netmax,nbkupn,nffg,nffgmx,ngexrt,nobswt,nodbmx,
     $    nopgmx,noprmx,noptmx,nordmx,nprpmx,nprpti,nprsmx,nprsti,
     $    nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,
     $    ntitmx,ntrymx,nttkmx,nxmdmx,nxmod,nxopex,nxopmx,nxopt,
     $    nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,pressb,
     $    pressi,ptk,qgexsh,rkb,rxbari,sfcar,ssfcar,tempcb,tempci,
     $    tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,trkb,ttk,
     $    ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,ugersi,
     $    ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,ureac,
     $    usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,uxopt,
     $    uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,
     $    xlkmod,xvfgex,zgexj,zvclgi)
        endif
c
        if (iopt(16) .eq. 0) then
c
c         Switch to write on the other backup file the next time.
c
          if (nbkupn .eq. nbkupa) then
            nbkupn = nbkupb
          elseif (nbkupn .eq. nbkupb) then
            nbkupn = nbkupa
          endif
        endif
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      qdmpr1 = .false.
      if (qftpr2) qdmpr2 = .false.
      kstppr = 0
c
c     Plot a description of the system at the current point of
c     reaction progress.
c
  360 continue
      kstppl = 0
c 360 if (iplot .ge. 1) then
c
c       Note: grafz has been removed from EQ6. Plotting is now done
c       via the tab file
c
c       call grafz
c       kstppl = 0
c     endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check stop conditions.
c
      call chkstc(actw,awmax,awmin,eh,ehmax,ehmin,fo2lg,iopt,
     $ jreac,kstep,kstpmx,noptmx,noutpt,nrct,nrctmx,nttyo,o2max,o2min,
     $ ph,phmax,phmin,prcinf,qaft1,qcnpre,qcntmp,qconst,qredox,qstop,
     $ qvhfxi,qvlsow,timemx,time1,tolxst,tolxsu,ximax,xi1)
c
      if (qstop .or. qstopx) go to 990
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Reset print, plot, and PRS transfer points as necessary.
c
      call adprpl(actw,aw0plo,aw0prn,aw1plo,aw1prn,dlaplo,
     $ dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,
     $ dltplo,dltprl,dltprn,dlxpll,dlxplo,dlxprl,dlxprn,dlxdmp,
     $ eh,eh0plo,eh0prn,eh1plo,eh1prn,fo2lg,fxprpl,o20plo,o20prn,
     $ o21plo,o21prn,ph,ph0plo,ph0prn,ph1plo,ph1prn,qplaw0,qplaw1,
     $ qpleh0,qpleh1,qplolt,qplolx,qplott,qplotx,qplo20,qplo21,
     $ qplph0,qplph1,qpraw0,qpraw1,qpreh0,qpreh1,qprnlx,qprnlt,
     $ qprntx,qprntt,qpro20,qpro21,qprph0,qprph1,qredox,time1,
     $ tiplol,tiplot,tiprnl,tiprnt,xidump,xiplol,xiplot,xiprnl,
     $ xiprnt,xi1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     If shooting to the end, go do the next step.
c
      if (qshoot) then
        nord = 0
        delxi = dlxmin
        go to 100
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Continue the reaction path simulation. If the current point
c     was hit only because of a print or plot requirement, and delxi
c     is very small relative to the preceding step size, the current
c     point is dropped as the base point (the point about which the
c     Taylor's series expansion is made) in favor of the preceding
c     point. This is done to protect the integrity of the finite
c     difference functions, by avoiding the case in which two points
c     are nearly indistinguishable. Here qskip is a logical flag to
c     skip recalculating the finite differences and derivatives.
c
      qreax = iexrt.gt.0 .or. jexrt.gt.0 .or. jscat.gt.0 .or.
     $ jscrt.gt.0 .or. qreq
      qtrch = .false.
      iexrt = 0
      jexrt = 0
      jscat = 0
      jscrt = 0
      qskip = .false.
      if (qzprnt. or. qzplot) then
        if (qdump) go to 100
        if (kord .le. 0) go to 100
        if (delxi .ge. 0.01*(xi0 - xim1)) go to 100
      else
        go to 100
      endif
c
      qskip = .true.
      nord = nordsv
      kord = kordsv
      delxi = max(delxi,dlxisv)
      xi1 = xi0
      time1 = time0
c
      km1 = km10
      kmt = kmt0
      kx1 = kx10
      kxt = kxt0
      kdim = kdim0
c
      call copyaa(zvclg0,zvclg1,kdim)
      call copyaa(zvec0,zvec1,kdim)
c
      do kcol = 1,kdim
        uzvec1(kcol) = uzvec0(kcol)
        iindx1(kcol) = iindx0(kcol)
        ipndx1(kcol) = ipndx0(kcol)
      enddo
c
      if (iopt(2) .gt. 0) then
        rirec1 = rirec0
c
        call copyaa(rreac0,rreac1,nrct)
        call copyaa(rrelr0,rrelr1,nrct)
        call copyaa(sfcar0,sfcar,nrct)
      endif
c
      call copyaa(afrc0,afrc1,nrct)
c
      go to 110
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  990 if (kstpmx .le. 0) then
c
c       Normal termination at the initial point.
c
        write (noutpt,2200)
        write (nttyo,2200)
 2200   format(/,' ---  The reaction path terminated normally at',
     $  ' the initial point ---',//)
c
        write (noutpt,2210) xistrt
        write (nttyo,2210) xistrt
 2210   format(7x,'Xi is ',1pe12.5)
c
        write (noutpt,2220) kdim
        write (nttyo,2220) kdim
 2220   format(7x,'The matrix dimension is ',i4,/)
c
        if (iopt(2) .gt. 0) then
          tistrd = tistrt/86400.
          tistry = tistrd/365.25
          write (noutpt,2230) tistrt,tistrd,tistry
          write (nttyo,2230) tistrt,tistrd,tistry
 2230     format(7x,'The time is:',
     $      /14x,1pe12.5,' seconds',
     $      /14x,1pe12.5,' days',
     $      /14x,1pe12.5,' years',/)
          if (qriinf) then
            write (noutpt,2240)
            write (nttyo,2240)
 2240       format(7x,'The system is now indistinguishable from',
     $      /7x,'what it would be at time equals infinity.')
          endif
        endif
      else
c
        if (.not.qstopx) then
c
c         Normal termination of a reaction path calculation.
c
          write (noutpt,2300)
          write (nttyo,2300)
 2300     format(/,' ---  The reaction path has terminated',
     $    ' normally ---',//)
        else
c
c         Early termination of a reaction path calculation.
c
          write (noutpt,2310)
          write (nttyo,2310)
 2310     format(/,' ---  The reaction path has terminated',
     $    ' early ---',//)
        endif
c
        write (noutpt,2320) kstep,xistsv,xi1
        write (nttyo,2320) kstep,xistsv,xi1
 2320   format(7x,i5,' steps were taken',/7x,'Xi increased from ',
     $  /14x,1pe12.5,' to ',1pe12.5)
c
        avdlxi = 0.
        if (kstep .gt. 0) avdlxi = (xi1 - xistsv)/kstep
        iavkdm = nint(avkdim)
        write (noutpt,2330) avdlxi,iavkdm
        write (nttyo,2330) avdlxi,iavkdm
 2330   format(7x,'The average value of delxi was ',1pe12.5,
     $  /7x,'The average matrix dimension was ',i4,/)
c
        if (iopt(2) .gt. 0) then
          tistrd = tistrt/86400.
          tistry = tistrd/365.25
          if (qriinf) then
            write (noutpt,2340) tistrt,tistrd,tistry
            write (nttyo,2340) tistrt,tistrd,tistry
 2340       format(7x,'The time increased from ',
     $      /14x,1pe12.5,' seconds to infinity',
     $      /14x,1pe12.5,' days to infinity',
     $      /14x,1pe12.5,' years to infinity',/)
          else
            tdays = time1/86400.
            tyears = tdays/365.25
            write (noutpt,2350) tistrt,time1,tistrd,tdays,tistry,tyears
            write (nttyo,2350) tistrt,time1,tistrd,tdays,tistry,tyears
 2350       format(7x,'The time increased from ',
     $      /14x,1pe12.5,' to ',1pe12.5,' seconds',
     $      /14x,1pe12.5,' to ',1pe12.5,' days',
     $      /14x,1pe12.5,' to ',1pe12.5,' years',/)
          endif
        endif
      endif
c
      write (ux16a,'(f9.4)') phstrt
      call lejust(ux16a)
      j2 = ilnobl(ux16a)
      write (ux16b,'(f9.4)') ph
      call lejust(ux16b)
      j3 = ilnobl(ux16b)
      if (ph .gt. phstrt) then
        write (noutpt,2360) ux16a(1:j2),ux16b(1:j3)
        write (nttyo,2360) ux16a(1:j2),ux16b(1:j3)
 2360   format(7x,'The pH increased from ',a,' to ',a)
      elseif (ph .lt. phstrt) then
        write (noutpt,2362) ux16a(1:j2),ux16b(1:j3)
        write (nttyo,2362) ux16a(1:j2),ux16b(1:j3)
 2362   format(7x,'The pH decreased from ',a,' to ',a)
      endif
c
      if (qredox) then
        write (ux16a,'(f9.4)') ehstrt
        call lejust(ux16a)
        j2 = ilnobl(ux16a)
        write (ux16b,'(f9.4)') eh
        call lejust(ux16b)
        j3 = ilnobl(ux16b)
        if (eh .gt. ehstrt) then
          write (noutpt,2370) ux16a(1:j2),ux16b(1:j3)
          write (nttyo,2370) ux16a(1:j2),ux16b(1:j3)
 2370     format(7x,'The Eh increased from ',a,' to ',a,' v')
        elseif (eh .lt. ehstrt) then
          write (noutpt,2372) ux16a(1:j2),ux16b(1:j3)
          write (nttyo,2372) ux16a(1:j2),ux16b(1:j3)
 2372     format(7x,'The Eh decreased from ',a,' to ',a,' v')
        endif
c
        write (ux16a,'(f9.4)') o2strt
        call lejust(ux16a)
        j2 = ilnobl(ux16a)
        write (ux16b,'(f9.4)') fo2lg
        call lejust(ux16b)
        j3 = ilnobl(ux16b)
        if (fo2lg .gt. o2strt) then
          write (noutpt,2380) ux16a(1:j2),ux16b(1:j3)
          write (nttyo,2380) ux16a(1:j2),ux16b(1:j3)
 2380     format(7x,'The log fO2 increased from ',a,' to ',a)
        elseif (eh .lt. ehstrt) then
          write (noutpt,2382) ux16a(1:j2),ux16b(1:j3)
          write (nttyo,2382) ux16a(1:j2),ux16b(1:j3)
 2382     format(7x,'The log fO2 decreased from ',a,' to ',a)
        endif
      endif
c
      write (ux16a,'(f9.4)') awstrt
      call lejust(ux16a)
      j2 = ilnobl(ux16a)
      write (ux16b,'(f9.4)') actw
      call lejust(ux16b)
      j3 = ilnobl(ux16b)
      if (actw .gt. awstrt) then
        write (noutpt,2390) ux16a(1:j2),ux16b(1:j3)
        write (nttyo,2390) ux16a(1:j2),ux16b(1:j3)
 2390   format(7x,'The aw increased from ',a,' to ',a)
      elseif (actw .lt. awstrt) then
        write (noutpt,2392) ux16a(1:j2),ux16b(1:j3)
        write (nttyo,2392) ux16a(1:j2),ux16b(1:j3)
 2392   format(7x,'The aw decreased from ',a,' to ',a)
      endif
c
      write (ux16a,'(1pg12.5)') wwstrt
      call lejust(ux16a)
      j2 = ilnobl(ux16a)
      write (ux16b,'(1pg12.5)') wkgh2o
      call lejust(ux16b)
      j3 = ilnobl(ux16b)
      if (wkgh2o .gt. wwstrt) then
        write (noutpt,2394) ux16a(1:j2),ux16b(1:j3)
        write (nttyo,2394) ux16a(1:j2),ux16b(1:j3)
 2394   format(7x,'The mass of solvent water increased from ',a,
     $  ' to ',a,' kg')
      elseif (wkgh2o .lt. wwstrt) then
        write (noutpt,2396) ux16a(1:j2),ux16b(1:j3)
        write (nttyo,2396) ux16a(1:j2),ux16b(1:j3)
 2396   format(7x,'The mass of solvent water decreased from ',a,
     $  ' to ',a,' kg')
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (.not.qstopx) then
        if (iopt(7).gt.0 .and. iopt(1).ne.2) then
c
c         Clear ES solids at the end of the run, unless the run
c         terminated early due to calculational problems.
c         Fictive fugacity-fixing minerals are not cleared.
c
          write (noutpt,2400)
          write (nttyo,2400)
 2400     format(/' * Note (EQ6/path) Clearing equilibrium system',
     $    ' (ES) solids',/7x,'at the end of the run.')
c
          call clress(csts,iindx1,ipndx1,jpflag,jsflag,kdim,kmax,
     $    km1,kmt,kx1,kxt,loph,losp,moph,mosp,mtb,mtbaq,nbt,nbtmax,
     $    nptmax,nstmax,nsts,nstsmx,nstsr,ufixf,uzvec1,zvec1,zvclg1)
        endif
c
        if (iopt(10) .gt. 0) then
c
c         Clear PRS solids at the end of the run, unless the run
c         terminated early due to calculational problems.
c
          write (noutpt,2410)
          write (nttyo,2410)
 2410     format(/' * Note (EQ6/path) Clearing physically removed',
     $    ' system (PRS) solids',/7x,'at the end of the run.')
          call initaz(mprsp,nstmax)
          call initaz(mprph,nptmax)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nprob .le. 1) then
c
c       Calculate data to describe the aqueous solution in the first
c       problem on the input file as a special reactant. This defines
c       the so-called "Fluid 2" to be described on the EQ6 pickup
c       file under the iopt(20) = 1 option.
c
        ureac1 = 'Fluid 2'
c
c       Scale the composition and reaction so that one "mole" of the
c       fluid contains 1 kg of solvent water.
c
        morrw1 = mosp(narn1)/omega
        wx = 1./morrw1
c
        do nc = 1,nct
          n = nc
          uesr1(n) = uelem(nc)
          cesr1(n) = wx*mteaq(nc)
        enddo
        iesrt1 = nct
c
        n = 1
        ubsr1(n) = ureac1
        cbsr1(n) = -1.0
        do nb = 1,nbt
          ns = nbaspd(nb)
          if (jflag(ns) .lt. 30) then
            n = n + 1
            ubsr1(n) = uspec(ns)
            cbsr1(n) = wx*mtbaq(nb)
          endif
        enddo
        ibsrt1 = n
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(17) .ge. 0) then
c
c       Prepare to write results at the final point to the
c       pickup file.
c
        call setpk6(actwlg,awmax,awmaxi,awmin,awmini,eh,ehmax,
     $  ehmaxi,ehmin,ehmini,fo2lg,iindx1,jflag,jflgi,kbt,kdim,kmax,
     $  kprs,mprph,mprphi,mprsp,mprspi,mtb,mtbi,mtbaq,mtbaqi,nbasp,
     $  nbaspd,nbaspi,nbti,nbtmax,ncmpr,nobswt,noutpt,nprpmx,nprpti,
     $  nprsmx,nprsti,npt,nptmax,nttyo,nstmax,o2max,o2maxi,o2min,
     $  o2mini,ph,phmax,phmaxi,phmin,phmini,prcinf,press,pressi,
     $  tempc,tempci,time1,timemx,timmxi,tistti,ubmtbi,uobsw,uphase,
     $  uprphi,uprspi,uspec,uzveci,uzvec1,xi1,ximax,ximaxi,xistti,
     $  zvclgi,zvclg1)
c
        if (iopt(20) .gt. 0) then
c
c         Prepare to write an advanced EQ6 pickup file. There is
c         currently only one option. It is determined by the iopt(20)
c         option switch. See comments in EQ6/stpkmd.f.
c
          call stpkmd(cbsri,cbsr1,cdac,cesri,cesr1,csigma,eact,
     $    fkrc,hact,iact,ibsrti,ibsrt1,iesrti,iesrt1,iktmax,imchmx,
     $    imech,iopt,ixrti,jcode,jreac,modr,morr,morrw1,nbt1mx,nctmax,
     $    ndact,ndctmx,noptmx,noutpt,nprob,nrct,nrctmx,nrk,nsk,nsrt,
     $    nsrtmx,ntitl1,ntitmx,nttyo,nxrt,nxrtmx,rkb,rxbari,sfcar,
     $    ssfcar,trkb,ubsri,ubsr1,ucxri,udac,uesri,uesr1,ureac,ureac1,
     $    utitl1,vreac)
        endif
c
c       Write results at the final point to the pickup file.
c
        if (upkfor(1:1) .eq. 'W') then
c
          call wr6pkw(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,
     $    dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,
     $    dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,
     $    dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,iact,
     $    ibsrti,iesrti,ietmax,igerti,iktmax,imchmx,imech,iodb,iopg,
     $    iopr,iopt,itermx,ixrti,jcode,jetmax,jflgi,jgerti,jgext,
     $    jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,
     $    kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,
     $    mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,
     $    nertmx,net,netmax,newin,nffg,nffgmx,ngexrt,nobswt,nodbmx,
     $    nopgmx,noprmx,noptmx,nordmx,nprpmx,nprpti,nprsmx,nprsti,
     $    nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,
     $    ntitmx,ntrymx,nttkmx,nxmdmx,nxmod,nxopex,nxopmx,nxopt,
     $    nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,pressb,
     $    pressi,ptk,qgexsh,rkb,rxbari,sfcar,ssfcar,tempcb,tempci,
     $    tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,trkb,ttk,
     $    ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,ugersi,
     $    ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,ureac,
     $    usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,uxopt,
     $    uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,
     $    xlkmod,xvfgex,zgexj,zvclgi)
        else
c
c         Menu-style (D) format.
c
          call wr6pkd(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,
     $    dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,
     $    dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,
     $    dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,iact,
     $    ibsrti,iesrti,ietmax,igerti,iktmax,imchmx,imech,iodb,iopg,
     $    iopr,iopt,itermx,ixrti,jcode,jetmax,jflgi,jgerti,jgext,
     $    jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,
     $    kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,
     $    mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,
     $    nertmx,net,netmax,newin,nffg,nffgmx,ngexrt,nobswt,nodbmx,
     $    nopgmx,noprmx,noptmx,nordmx,nprpmx,nprpti,nprsmx,nprsti,
     $    nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,
     $    ntitmx,ntrymx,nttkmx,nxmdmx,nxmod,nxopex,nxopmx,nxopt,
     $    nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,pressb,
     $    pressi,ptk,qgexsh,rkb,rxbari,sfcar,ssfcar,tempcb,tempci,
     $    tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,trkb,ttk,
     $    ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,ugersi,
     $    ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,ureac,
     $    usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,uxopt,
     $    uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,
     $    xlkmod,xvfgex,zgexj,zvclgi)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
