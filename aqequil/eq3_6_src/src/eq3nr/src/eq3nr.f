      program eq3nr
c
c     This is the main program of the EQ3NR code. Configuration
c     identification, the copyright statement, legal disclaimers, and
c     similar statements are contained in EQ3NR/aaaeq3.f, the lead-off
c     subroutine in the EQ3NR source code. A short description of this
c     program is also contained in that subroutine.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      include 'eqlib/eqldef.h'
      include 'eqlib/eqlpar.h'
      include 'eqlib/eqldv.h'
      include 'eqlib/eqlge.h'
      include 'eqlib/eql1s.h'
      include 'eqlib/eqlwd.h'
c
      include 'eqlib/eqlj8.h'
      include 'eqlib/eqlk8.h'
      include 'eqlib/eqlo8.h'

c-----------------------------------------------------------------------
c     File path parameters
c-----------------------------------------------------------------------
      integer :: numargs
      character(1024) :: temppath
      character(:), allocatable :: data1path
      character(:), allocatable :: threeipath
      integer :: pathindices(2)
      character(:), allocatable :: basename
      character(:), allocatable :: ofile
      character(:), allocatable :: ifile
      character(:), allocatable :: pfile
c
c-----------------------------------------------------------------
c
c
c     Array allocation size variables used in EQ3NR.
c
c       General size variables:
c
c         ntid_asv  = the number of lines in the data file title
c         ntit_asv  = the number of lines in any title appearing
c                       on an input or pickup file
c
c
c       Size variables associated with the description of thermodynamic
c       data:
c
c         ipch_asv  = the order for pressure corrections to enthalpy
c                       functions
c         ipcv_asv  = the order for pressure corrections to volume
c                       functions; the maximum order for pressure
c                       corrections to log K and other Gibbs-energy-
c                       based functions is one greater than this
c         ipbt_asv  = the maximum number of Pitzer alpha parameters
c                       for any species pair
c         jpfc_asv  = the number of coefficients in the Pitzer parameter
c                       temperature function
c         narx_asv  = maximum number of coefficients per temperature
c                       range
c         ntpr_asv  = number of temperature ranges
c
c
c       Size variables required to read other data from the
c       supporting data file:
c
c         ikta_asv  = the maximum number of end-member components
c                       of any solid solution on the data file
c         nata_asv  = the number of aqueous species on the data file
c         nbta_asv  = the number of basis species on the data file
c         ncta_asv  = the number of chemical elements on the data file
c         ngta_asv  = the number of gas species on the data file
c         nlta_asv  = the number of pure liquids on the data file
c         nmta_asv  = the number of pure minerals on the data file
c
c         napa_asv  = number of distinct sets of Pitzer alpha
c                       parameters on the data file
c         nmuta_asv = the number of triplets of ions on the data file
c                       for which distinct Pitzer mu coefficients
c                       are defined
c         nslta_asv = the number of pairs of ions on the data file
c                       for which distinct Pitzer lambda coefficients
c                       are defined
c
c
c       Size variables required to hold the compressed set of such
c       data:
c
c         ikt_asv   = the maximum number of end-member components
c                       of any solid solution in the compressed set
c         nat_asv   = the number of aqueous species in the compressed
c                       set
c         nbt_asv   = the number of basis species in the compressed set
c         nct_asv   = the number of chemical elements in the compressed
c                       set
c         ngt_asv   = the number of gas species in the compressed set
c         nlt_asv   = the number of pure liquids in the compressed set
c         nmt_asv   = the number of pure minerals in the compressed set
c
c         nap_asv   = the number of distinct sets of Pitzer alpha
c                       parameters in the compressed set
c         nmut_asv  = the number of triplets of ions in the compressed
c                       set for which distinct Pitzer mu coefficients
c                       are defined
c         nslt_asv  = the number of pairs of ions in the compressed set
c                       for which distinct Pitzer lambda coefficients
c                       are defined
c
c
      integer ntid_asv
c
      integer ipch_asv,ipcv_asv,ipbt_asv,jpfc_asv,narx_asv,ntpr_asv
c
      integer nazm_asv,nazp_asv
c
      integer ikta_asv,napa_asv,nata_asv,nbta_asv,ncta_asv,ndrsa_asv,
     $ nessa_asv,ngta_asv,nlta_asv,nmta_asv,nmuta_asv,nslta_asv,
     $ npta_asv,nsta_asv,nxta_asv
c
      integer iapxa_asv,ibpxa_asv
c
      integer nbta1_asv
c
      integer ntf1a_asv,ntf2a_asv
c
      integer ikt_asv,nap_asv,nat_asv,nbt_asv,nct_asv,ndrs_asv,
     $ ness_asv,ngt_asv,nlt_asv,nmt_asv,nmut_asv,nslt_asv,npt_asv,
     $ nst_asv,nxt_asv
c
      integer iapx_asv,ibpx_asv
c
      integer ntit_asv,nxmd_asv
c
      integer nbt1_asv,nmx_asv,nsts_asv,nsx_asv,nxic_asv
c
      integer ntf1_asv,ntf2_asv,ntfx_asv
c
      integer k_asv
c
      integer nxti_asv
c
      integer imch_asv,ndct_asv,nert_asv,nffg_asv,nprp_asv,nprs_asv,
     $ nptk_asv,nrct_asv,nsrt_asv,nttk_asv,nxop_asv,nxpe_asv,nxrt_asv
c
c-----------------------------------------------------------------------
c
      integer newin,ninpt,ninpts,noutpt,nttyo
c
      integer iodb(nodb_par),iopg(nopg_par),iopr(nopr_par),
     $ iopt(nopt_par)
c
      integer jgexti(net_par),ngexpi(net_par),ngexti(jet_par,net_par)
c
      integer, dimension(:), allocatable :: iapxt,ibpxt,ibswx,iction,
     $ igstak,iindx1,insgf,ipivot,ipndx1,istack,ixbasp,jcsort,jflag,
     $ jflagd,jgsort,jgstak,jjndex,jjsort,jpflag,jsflag,jsitex,jsol,
     $ jssort,jstack,kction,kkndex,kxmod
c
      integer, dimension(:), allocatable :: jflgi
c
      integer, dimension(:), allocatable :: narxt,nbasp,nbaspd,nbaspx,
     $ nbmap,ncmap,ncosp,ndecsp,ndrs,ndrsd,ndrsx,ness,nfac,nphasx,
     $ npnxp,nsmap,nsts,ntfx,ntf1,ntf2
c
      integer, dimension(:), allocatable :: nbaspi
c
      integer, dimension(:,:), allocatable :: ncmpr,ndrsr,ndrsrd,
     $ ndrsrx,nessr,nstsr
c
      integer, dimension(:,:), allocatable :: ncmpri
c
      integer, dimension(:), allocatable :: nalpha
c
      integer, dimension(:,:), allocatable :: nmux,nmxi,nmxx,nslx,
     $ nsxi,nsxx
c
      integer ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,
     $ ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta
c
      integer i,iaqsla,iaqsln,ibetmx,icarb,idelmx,ie,iebal,iebal3,ielam,
     $ ier,iern1,iern2,iexec0,ifrn1,ifrn2,igas,ilevel,ilrn1,ilrn2,imrn1,
     $ imrn2,io2gaq,ipch,ipcv,irang,irdxc3,itdsf3,iter,itermx,itmxsv,
     $ ixrn1,ixrn2,ixrn1a,ixrn2a,izmax,j,je,jexec0,jfl,jfleba,jlen,
     $ jpdblo,jpres3,jptffl,j1,j2,j3,j4,k,ka1,kat,kbt,kcarb,kcol,kct,
     $ kdim,kebal,kelect,ker,ke1,ket,khydr,km1,kmt,kobswt,ko2gaq,kprs,
     $ kwater,kx1,kxt
c
      integer n,nad1,napt,napta,narn1,narn2,narn1a,narn2a,nat,nata,
     $ nb,nbi,nbt,nbta,nbtafd,nbtd,nbti,nbw,nbwa,nb1,nb2,nc,ncarb,
     $ nchloa,nchlor,nct,ncta,ne,neleca,nelect,nern1,nern2,nerr,net,
     $ neti,ngrn1,ngrn1a,ngrn2a,ngrn2,ngt,ngta,nhydr,nhydra,nhydx,
     $ nhydxa,nlrn1,nlrn2,nlrn1a,nlrn2a,nlt,nlta,nmax,nmrn1,nmrn2,
     $ nmrn1a,nmrn2a,nmt,nmta,nmut,nmuta,nobswt,no2gaa,no2gai,no2gaq,
     $ np,npt,npta,nprob,nrecl,nrdxsa,nrdxsp,nredox,nrr1,nrr2,nr1,nr2,
     $ ns,nsbsw,nsbswt,nse,nslt,nslta,nss,nst,nsta,ns1,ns2,nt,ntfxt,
     $ ntf1t,ntf1ta,ntf2t,ntf2ta,ntitl,ntitl2,ntitld,ntpr,ntprt,nxi,
     $ nxic,nxmod,nxrn1,nxrn2,nxrn1a,nxrn2a,nxt,nxta,nxti
c
      integer ilnobl,nbasis
c
      logical, dimension(:), allocatable :: qxknph
c
      logical qbassw,qchlor,qbswok,qcwrpj,qdwipp,qelim,qend,qhawep,qop,
     $ qpit75,qredox,qrderr,qrho,q6mode
c
      character(len=80), dimension(:), allocatable :: utitl,utitl2
      character(len=48), dimension(:), allocatable :: ubasp,ubmtbi,
     $ ucospi,uspec,uspeci,uxmod,uzveci,uzvec1
      character(len=48), dimension(:,:), allocatable :: uobsw,usbsw
      character(len=24), dimension(:), allocatable :: umemi,usoli
      character(len=24), dimension(:), allocatable :: uphase,uptype
      character(len=8), dimension(:), allocatable :: uelem,uldel,ulbeta
c
      character(len=32) ujflls(0:njf_par),uxtype(jso_par)
      character(len=24) ugexpi(net_par),ugexsi(iet_par,jet_par,net_par)
      character(len=8) ugexji(jet_par,net_par)
c
      character(len=56) uspn56
      character(len=48) ubacmx,ubbig,ubetmx,ubgamx,ubneg
      character(len=32) uactop
      character(len=24) uaqsln,ublk24,uebal,uredox,ux24
      character(len=16) ux16a,ux16b
      character(len=11) utime0,utime1
      character(len=9) udate0,udate1
      character(len=8) udatfi,udakey,uinfor,upkfor,uplatc,uplatm,ustelg,
     $ ustelu,usteql,usteq3,uveelg,uveelu,uveeql,uveeq3
c
      character(len=8) uv,ux8,ux8a,ux8b
      character(len=80) ux80
c
      real(8), dimension(:), allocatable :: covali
c
      real(8), dimension(:), allocatable :: acflg,acflgo,act,actlg,
     $ affpd,affsd,ahrc,alpha,amtb,atwt,azero,a3bars,beta,betao,bfac,
     $ cdrs,cdrsd,cdrsx,cdrtw,cdrw,cess,cjbasp,cnufac,conc,conclg,coval,
     $ csts,ctb,cteaq,delvco,delvec,dlogxw,efac,ehrc,fo2lrc,fsort,fugac,
     $ fugalg,loph,losp,lsort,moph,mosp,mtb,mtbaq,mtbaqi,mtbi,mte,
     $ mteaq,mwtsp,perc,ppmwe,rhsvec,sidrph,sidrsp,tempcu,tfx,tf1,tf2,
     $ vosp0,weight,xbar,xbari,xbarlg,xlkmod,zchar,zchcu6,zchsq2,
     $ zvclgi,zvclg1,zvec1
c
      real(8), dimension(:,:), allocatable :: amu
      real(8), dimension(:,:,:), allocatable :: aslm
c
      real(8), dimension(:), allocatable :: pmu,pslm
      real(8), dimension(:,:), allocatable :: dpslm,gpit,palpha,pslamn
      real(8), dimension(:,:,:), allocatable :: dgpit
c
      real(8), dimension(:,:), allocatable :: elam,pelm
      real(8), dimension(:,:,:), allocatable :: delam,dpelm
c
      real(8), dimension(:), allocatable :: selm
      real(8), dimension(:,:), allocatable :: dselm
c
      real(8), dimension(:,:), allocatable :: aprehw,apresg
c
      real(8), dimension(:), allocatable :: dadhh,dadhv,dbdhh,dbdhv,
     $ dbdth,dbdtv
c
      real(8), dimension(:), allocatable :: xhfs,xhfsd,xlks,xlksd,
     $ xvfs,xvfsd
      real(8), dimension(:), allocatable :: dhfe,dvfe
      real(8), dimension(:,:), allocatable :: axhfe,axlke,axvfe
      real(8), dimension(:,:,:), allocatable :: adhfe,advfe
c
      real(8), dimension(:,:), allocatable :: dhfs,dhfsd,dvfs,dvfsd
      real(8), dimension(:,:,:), allocatable :: axhfs,axhfsd,axhfsx,
     $ axlks,axlksd,axlksx,axvfs,axvfsd,axvfsx
      real(8), dimension(:,:,:,:), allocatable :: adhfs,adhfsd,adhfsx,
     $ advfs,advfsd,advfsx
c
      real(8), dimension(:,:), allocatable :: apx,bpx,wfac
      real(8), dimension(:,:), allocatable :: aamatr,gmmatr
c
      real(8), dimension(:,:), allocatable :: aadh,aadhh,aadhv,aaphi,
     $ abdh,abdhh,abdhv,abdot,abdoth,abdotv
c
      real(8), dimension(:,:,:), allocatable :: adadhh,adadhv,adbdhh,
     $ adbdhv,adbdth,adbdtv
c
      real(8) cco2(5)
      real(8) cegexs(iet_par,jet_par,net_par),cgexpi(net_par),
     $ cgexp(net_par),cpgexs(iet_par,jet_par,net_par),
     $ egexjc(jet_par,net_par),egexjf(jet_par,net_par),
     $ egexpa(net_par),egexpc(net_par),egexs(iet_par,jet_par,net_par),
     $ egexsi(iet_par,jet_par,net_par),egexw(ket_par,net_par),
     $ mrgexs(iet_par,jet_par,net_par),
     $ xgexsi(iet_par,jet_par,net_par),xgexw(ket_par,net_par)
c
      real(8) adh,adhh,adhv,aphi,bdh,bdhh,bdhv,bdot,bdoth,bdotv
c
      real(8) prehw,presg,xhfe,xlke,xvfe
c
      real(8) abar,actwlc,afcnst,alki,al10,avcnst,a3bar,azch,azchmx,
     $ bacfmx,bbig,betamx,bfje,bgamx,bneg,bsigmm,bfxi,delmax,eh,
     $ ehfac,ehi,electr,eps100,farad,fje,fjeo,fo2,fo2lg,fo2lgi,fxi,
     $ fxio,omega,omeglg,pe,pei,presmx,press,pressi,rconst,rcnstv,
     $ rho,rtcnst,scamas,screwd,screwn,sigmam,sigmmo,sigzi,smp100,
     $ tcpu,tdamax,tdamin,tdspkg,tdspl,tempc,tempci,tempk,texec0,
     $ tolbt,toldl,tolspf,trun,tuser,vosol,wfh2o,wftds,wkgwi,woh2o,
     $ wosol,wotds,xbarw,xbarwc,xbrwlc,xbrwlg,x10,zx
cXX   real(8) vpgstp
c
      real(8) mlmrra,mrmlra,rhoc,rhowc,tdsgks,tdsglw,tdspkc,tdsplc
c
      real(8) axx,av,dp,pxl,pxu,xx
c
      real(8) texp,tlg
c
c-----------------------------------------------------------------------
c
c     Variable declarations: Data file original contents. These
c     variables and arrays remain unchanged after being filled by
c     reading the data file. These data comprise the 'a' set.
c
      integer, dimension(:), allocatable :: iapxta,ibpxta,insgfa,jsola,
     $ nbaspa,ndrsa,nessa
c
      integer, dimension(:,:), allocatable :: ncmpra,ndrsra,nessra
c
      real(8), dimension(:,:,:), allocatable :: axhfsa,axlksa,axvfsa
c
      real(8), dimension(:,:,:,:), allocatable :: adhfsa,advfsa
c
      real(8), dimension(:), allocatable :: atwta,azeroa,cdrsa,cessa,
     $ mwtspa,vosp0a,zchara
c
      real(8), dimension(:,:), allocatable :: apxa,bpxa
c
      character(len=80), dimension(:), allocatable :: utitld
      character(len=48), dimension(:), allocatable :: uspeca
      character(len=24), dimension(:), allocatable :: uphasa,uptypa
      character(len=8), dimension(:), allocatable :: uelema
c
c     Pitzer's equations subset.
c
      integer, dimension(:), allocatable :: nalpaa
      integer, dimension(:,:), allocatable :: nmuxa,nslxa
c
      real(8) , dimension(:,:), allocatable :: amua,palpaa
      real(8) , dimension(:,:,:), allocatable :: aslma
c
c-----------------------------------------------------------------------
c
c     Variable declarations: Variables closely associated with the
c     'a' set.
c
      logical, dimension(:), allocatable :: qclnsa
c
c-----------------------------------------------------------------------
c
c     Variable declarations: Pseudo-data file original contents.
c
      integer, dimension(:), allocatable :: ntf1a,ntf2a
c
      real(8), dimension(:), allocatable :: tf1a,tf2a
c
c-----------------------------------------------------------------------
c
c     Variable declarations: Data needed to write a complete EQ6 input
c     file as the EQ3NR pickup file.
c
      integer igerti(jet_par,nert_par),jgerti(nert_par)
c
      integer, dimension(:), allocatable :: ibsrti,iesrti,ixrti,
     $ jcode,jreac,nsk
c
      integer, dimension(:,:), allocatable :: imech,nrk
c
      integer, dimension(:,:,:), allocatable :: iact,ndact
c
      integer jpress,jtemp,ksplmx,ksppmx,kstpmx,nert,nffg,nprpti,nprsti,
     $ nrct,nsrt,ntitl1,ntrymx,nxopex,nxopt,nxrt
c
      logical qgexsh
c
      character(len=24) ugersi(iet_par,jet_par,nert_par),
     $ ugermo(nert_par)
      character(len=8) ugerji(jet_par,nert_par)
c
      character(len=80), dimension(:), allocatable :: utitl1
      character(len=48), dimension(:), allocatable :: uprspi
      character(len=24), dimension(:), allocatable :: uffg,uprphi,ureac,
     $ uxcat,uxopex
      character(len=8), dimension(:), allocatable :: uxopt
c
      character(len=24), dimension(:,:), allocatable :: ubsri,ucxri
      character(len=8), dimension(:,:), allocatable :: uesri
c
      character(len=24), dimension(:,:,:,:), allocatable :: udac
c
      real(8), dimension(:), allocatable ::  fkrc,modr,moffg,morr,
     $ mprphi,mprspi,ptk,sfcar,ssfcar,ttk,vreac,xlkffg
      real(8), dimension(:,:), allocatable :: cbsri,cesri,rxbari
      real(8), dimension(:,:,:), allocatable :: csigma,eact,hact,
     $ rkb,trkb
      real(8), dimension(:,:,:,:), allocatable :: cdac
c
      real(8) egersi(iet_par,jet_par,nert_par),
     $ xgersi(iet_par,jet_par,nert_par)
c
      real(8) awmaxi,awmini,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,
     $ dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,
     $ dlxplo,dlxprl,dlxprn,ehmaxi,ehmini,o2maxi,o2mini,phmaxi,phmini,
     $ pressb,tempcb,tempc1,timmxi,tistti,tolsat,tolxsf,ximaxi,xistti
c
c-----------------------------------------------------------------------
c
c     Variable declarations: Local data needed to assist in writing a
c     complete EQ6 input file as the EQ3NR pickup file. These variables
c     do not appear on that file itself.
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
c     BEGIN_MACHINE_DEPENDENT_CODE
c
c       On some systems, a BLOCK DATA subroutine must be declared in an
c       EXTERNAL statement to assure proper loading. On some other
c       systems, this is not necessary, but neither it is not harmful.
c       On yet some other systems, the EXTERNAL statement below may
c       cause a problem. If so, try commenting it out. If you still
c       have trouble, consult your local system documentation or
c       experiment to find out how to correctly handle a BLOCK DATA
c       subroutine on your system. The EXTERNAL statement below should
c       not cause a problem if you are using a compiler which is fully
c       compliant with the Fortran 90 standard. However, there is
c       no guarantee that it will be adequate to assure correct loading
c       of the BLOCK DATA subroutine.
c
        external bkdeq3
c
c     END_MACHINE_DEPENDENT_CODE
c
c-----------------------------------------------------------------------
c
c     Special data base of coefficients for the HCO3-CO3-OH total
c     alkalinity.
c
      include 'eqlib/eqltf1.h'
c
c     Special data base of coefficients for the extended total
c     alkalinity.
c
      include 'eqlib/eqltf2.h'
c
c-----------------------------------------------------------------------
c
c     Other data statements.
c
      data uaqsln /'Aqueous solution        '/
      data ublk24 /'                        '/
c
      data ujflls(0)  /'Total molality                  '/
      data ujflls(1)  /'Total molarity                  '/
      data ujflls(2)  /'Total mg/L                      '/
      data ujflls(3)  /'Total mg/kg.sol                 '/
      data ujflls(4)  /'ERROR                           '/
      data ujflls(5)  /'ERROR                           '/
      data ujflls(6)  /'ERROR                           '/
      data ujflls(7)  /'Total alkalinity, eq/kg.H2O     '/
      data ujflls(8)  /'Total alkalinity, eq/L          '/
      data ujflls(9)  /'Total alkalinity, eq/kg.sol     '/
      data ujflls(10) /'Total alkalinity, mg/L CaCO3    '/
      data ujflls(11) /'Total alkalinity, mg/L HCO3-    '/
      data ujflls(12) /'ERROR                           '/
      data ujflls(13) /'ERROR                           '/
      data ujflls(14) /'ERROR                           '/
      data ujflls(15) /'ERROR                           '/
      data ujflls(16) /'Log activity                    '/
      data ujflls(17) /'|zj| log ai +/- |zi| log aj     '/
      data ujflls(18) /'Log a(+/-,ij)                   '/
      data ujflls(19) /'pX                              '/
      data ujflls(20) /'pH                              '/
      data ujflls(21) /'pHCl                            '/
      data ujflls(22) /'pmH                             '/
      data ujflls(23) /'pmX                             '/
      data ujflls(24) /'ERROR                           '/
      data ujflls(25) /'Heterogenous equilibrium        '/
      data ujflls(26) /'ERROR                           '/
      data ujflls(27) /'Homogenous equilibrium          '/
      data ujflls(28) /'ERROR                           '/
      data ujflls(29) /'ERROR                           '/
      data ujflls(30) /'Make non-basis                  '/
c
      data uxtype(1)  /'Ideal solution                  '/
      data uxtype(2)  /'Binary, third-order Maclaurin   '/
      data uxtype(3)  /'Binary, parabolic Maclaurin     '/
      data uxtype(4)  /'Binary, cubic Maclaurin (P,T)   '/
      data uxtype(5)  /'Binary, Guggenheim  (T)         '/
      data uxtype(6)  /'Ternary, regular                '/
c
      data nrecl /0/
c
c     Under-relaxation control parameters:
c
c       screwd = bound on max norm of applied part of the delvec vector
c       screwn = factor bounding the increase in the max norm of the
c                beta vector
c
      data screwd/2.0/,screwn/0.50/
c
c-----------------------------------------------------------------------
c
c     Fixed fugacity phase range delimiters (not currently used
c     by EQ3NR):
c
      data ifrn1,ifrn2 /0,0/
c
c-----------------------------------------------------------------------
c
c     Maximum allowed pressure (bars):
c
      data presmx /90000./
c
c-----------------------------------------------------------------------
c
c     BEGIN_MACHINE_DEPENDENT_CODE
c
c       Define the console output device number. A value of 6 applies
c       in most cases.
c
        data nttyo  /6/
c
c     END_MACHINE_DEPENDENT_CODE
c
c-----------------------------------------------------------------------
c
c     Get time and date at the start of execution. This information
c     will be used for time and date stamping.
c
      call initim(iexec0,jexec0,texec0,noutpt,nttyo,
     $ udate0,utime0)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Disable underflow trapping, if any.
c
       call undflw
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set file unit numbers to zero.
c
      noutpt = 0
      nad1 = 0
      ninpt = 0
      ninpts = 0
      newin = 0
c
c     Open all files execpt pickup.
c
      numargs = COMMAND_ARGUMENT_COUNT()
      if (numargs.eq.2) then
          call GET_COMMAND_ARGUMENT(1,temppath)
          data1path = TRIM(temppath)
          temppath(:)='\0'
          call GET_COMMAND_ARGUMENT(2,temppath)
          threeipath = TRIM(temppath)
      else
          write (0, *) 'usage: eq3nr <data1> <3i>'
          stop
      end if

      call getbasename(threeipath, pathindices)
      basename = threeipath(pathindices(1):pathindices(2))

      ofile = basename // '.3o'
      ifile = basename // '.3ib'
      pfile = basename // '.3p'

      call openin(noutpt,nttyo,data1path,'unformatted',nad1)
      call openin(noutpt,nttyo,threeipath,'formatted',ninpt)

      call openou(noutpt,nttyo,ofile,'formatted',nrecl,noutpt)
      call openou(noutpt,nttyo,ifile,'formatted',nrecl,ninpts)
c
c     Make a copy of the input file, stripped of comments.
c
      call stripl(ninpt,ninpts)
      close (ninpt)
      ninpt = 0
      rewind ninpts
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get configuration identification data.
c
      call aaaeq3(usteq3,uveeq3)
      call aaaeql(usteql,uveeql)
      call aaaelg(ustelg,uveelg)
      call aaaelu(ustelu,uveelu)
      call platfd(uplatc,uplatm)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write configuration identification data, the copyright statement,
c     and any remaining statements or disclaimers.
c
      i = index(uveeq3,' ') - 1
      j = index(uveeq3,'.') - 1
      k = index(uplatc,' ') - 1
      if (i .le. 0) i = 8
      if (j .le. 0) j = 8
      if (k .le. 0) k = 8
      write (nttyo,1000) uveeq3(1:i),uveeq3(1:j),uveeq3(1:i),
     $ uplatc(1:k)
      write (noutpt,1000) uveeq3(1:i),uveeq3(1:j),uveeq3(1:i),
     $ uplatc(1:k)
 1000 format(/' EQ3/6, Version ',a,' (EQ3/6-V',a,'-REL-V',a,'-',a,')')
c
      i = j
      j = index(usteq3,' ') - 1
      k = index(uplatm,' ') - 1
      if (j .le. 0) j = 8
      if (k .le. 0) k = 8
      write (nttyo,1001) uveeq3(1:i),usteq3(1:j),uplatm(1:k)
      write (noutpt,1001) uveeq3(1:i),usteq3(1:j),uplatm(1:k)
 1001 format(' EQ3NR Speciation-Solubility Code (EQ3/6-V',a,
     $ '-EQ3NR-EXE-',a,'-',a,')')
c
      write (noutpt,1002)
      write (nttyo,1002)
 1002 format(' Supported by the following EQ3/6 libraries:')
c
      i = index(uveeql,'.') - 1
      j = index(usteql,' ') - 1
      if (i .le. 0) i = 8
      if (j .le. 0) j = 8
      write (nttyo,1003) uveeql(1:i),usteql(1:j),uplatm(1:k)
      write (noutpt,1003) uveeql(1:i),usteql(1:j),uplatm(1:k)
 1003 format('   EQLIB (EQ3/6-V',a,'-EQLIB-LIB-',a,'-',a,')')
c
      i = index(uveelg,'.') - 1
      j = index(ustelg,' ') - 1
      if (i .le. 0) i = 8
      if (j .le. 0) j = 8
      write (nttyo,1004) uveelg(1:i),ustelg(1:j),uplatm(1:k)
      write (noutpt,1004) uveelg(1:i),ustelg(1:j),uplatm(1:k)
 1004 format('   EQLIBG (EQ3/6-V',a,'-EQLIBG-LIB-',a,'-',a,')')
c
      i = index(uveelu,'.') - 1
      j = index(ustelu,' ') - 1
      if (i .le. 0) i = 8
      if (j .le. 0) j = 8
      write (nttyo,1005) uveelu(1:i),ustelu(1:j),uplatm(1:k)
      write (noutpt,1005) uveelu(1:i),ustelu(1:j),uplatm(1:k)
 1005 format('   EQLIBU (EQ3/6-V',a,'-EQLIBU-LIB-',a,'-',a,')',/)
c
      write (nttyo,1010)
      write (noutpt,1010)
 1010 format(' Copyright (c) 1987, 1990-1993, 1995, 1997, 2002 The',
     $ ' Regents of the',/' University of California, Lawrence',
     $ ' Livermore National Laboratory.',/' All rights reserved.',/)
c
c     Write additional statements and disclaimers.
c
      call prcndi(noutpt,nttyo)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the time and date on the output.
c
      j2 = ilnobl(udate0)
      write (noutpt,1080) utime0,udate0(1:j2)
      write (nttyo,1080) utime0,udate0(1:j2)
 1080 format(' Run',2x,a8,2x,a,/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the platform's real(8) floating-point parameters.
c
      call flpars(eps100,irang,noutpt,nttyo,smp100)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize array dimension variables corresponding to fixed
c     parameters.
c
c     The following are common with EQ6.
c
      ietmax = iet_par
      jetmax = jet_par
      ketmax = ket_par
      netmax = net_par
c
      jsomax = jso_par
      nodbmx = nodb_par
      nopgmx = nopg_par
      noprmx = nopr_par
      noptmx = nopt_par
      nvetmx = nvet_par
c
      nxmd_asv = nxmd_par
      nxmdmx = nxmd_asv
c
      ntit_asv = ntit_par
      ntitmx = ntit_asv
c
      iapxa_asv = iapxa_par
      ibpxa_asv = ibpxa_par
c
c     The following are common with EQ6, but are only needed to write
c     a full EQ6 input file as the EQ3NR pickup file.
c
      nrct_asv = nrct_par
      nert_asv = nert_par
      nsrt_asv = nsrt_par
      nxrt_asv = nxrt_par
      imch_asv = imch_par
      ndct_asv = ndct_par
      nffg_asv = nffg_par
      nprp_asv = nprp_par
      nprs_asv = nprs_par
      nptk_asv = nptk_par
      nttk_asv = nttk_par
      nxop_asv = nxop_par
      nxpe_asv = nxpe_par
c
c     Special over-rides for EQ3NR.
c
      nrct_asv = 2
      nsrt_asv = 1
      nxrt_asv = 1
      nprp_asv = 1
      nprs_asv = 1
      nxop_asv = 1
      nxpe_asv = 1
c
      nrctmx = nrct_asv
      nertmx = nert_asv
      nsrtmx = nsrt_asv
      nxrtmx = nxrt_asv
      imchmx = imch_asv
      ndctmx = ndct_asv
      nffgmx = nffg_asv
      nprpmx = nprp_asv
      nprsmx = nprs_asv
      nptkmx = nptk_asv
      nttkmx = nttk_asv
      nxopmx = nxop_asv
      nxpemx = nxpe_asv
c
c     The following are unique to EQ3NR.
c
      nxti_asv = nxti_par
c
      njfmax = njf_par
c
      nxtimx = nxti_asv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize the following floating-point constants:
c
c       rconst = the gas constant: 1.98726 cal/mol-K
c       rcnstv = the gas constant: 83.14510 bar-cm3/mol-K
c       vpgstp = the molar volume of a perfect gas at STP:
c                  22413.6 cm3/mol
c       farad  = the Faraday constant: 23062.3 cal/equiv-volt
c       al10   = ln 10; =~ 2.3026
c
      rconst = 1.98726
      rcnstv = 83.14510
cXX   vpgstp = 22413.6
      farad = 23062.3
      x10 = 10.
      al10 = log(x10)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set the EQ6 calculational mode flag (.false. in EQ3NR).
c
      q6mode = .false.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the header section of the data1 file. This section consists
c     of a record containing the string 'data1' (to ensure that the file
c     is indeed a data1 file), a record containing the keystring for the
c     activity coefficient model for aqueous species to which this data
c     file corresponds, and the array allocation size variables required
c     to allocate sufficient array space to store the rest of the data
c     on this data file.
c
      call indath(ikta_asv,ipbt_asv,ipch_asv,ipcv_asv,
     $ jpfc_asv,nad1,napa_asv,narx_asv,nata_asv,nbta_asv,
     $ ncta_asv,ngta_asv,nlta_asv,nmta_asv,npta_asv,nmuta_asv,
     $ noutpt,nslta_asv,nsta_asv,ntid_asv,ntpr_asv,nttyo,
     $ nxta_asv,udakey)
c
c     The value of nbta_asv at this point matches the number of
c     formally declared basis species on the data file. That is
c     the number of species in the "strict" and "auxiliary" basis
c     set sections of the aqueous species superblock. Save this
c     as nbtafd. Other species can be added to the basis set in
c     this version of EQ3/6. Thus, nbta_asv at this point is not
c     the "final answer."
c
      nbtafd = nbta_asv
c
c     Bump up the values of some allocation size variables to allow
c     for some species, phases, etc., to created by this software.
c     Examples: generic ion exchanger phases and species.
c
      nbta_asv = nbta_asv + 10
      npta_asv = npta_asv + net_par + 10
      nsta_asv = nsta_asv + iet_par*jet_par*net_par + 10
c
c     Set the values of some secondary allocation size variables.
c
      nbta1_asv = nbta_asv + 1
      ndrsa_asv = 7*nsta_asv
      nessa_asv = 5*nsta_asv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Allocate arrays to store the data read from the rest of the data1
c     file.
c
      ALLOCATE(iapxta(nxta_asv))
      ALLOCATE(ibpxta(nxta_asv))
      ALLOCATE(jsola(nxta_asv))
      ALLOCATE(narxt(ntpr_asv))
      ALLOCATE(nbaspa(nbta_asv))
      ALLOCATE(ndrsa(ndrsa_asv))
      ALLOCATE(nessa(nessa_asv))
c
      ALLOCATE(ncmpra(2,npta_asv))
      ALLOCATE(ndrsra(2,nsta_asv))
      ALLOCATE(nessra(2,nsta_asv))
c
      ALLOCATE(qclnsa(nsta_asv))
c
      ALLOCATE(utitld(ntid_asv))
      ALLOCATE(uspeca(nsta_asv))
      ALLOCATE(uphasa(npta_asv))
      ALLOCATE(uptypa(npta_asv))
      ALLOCATE(uelema(ncta_asv))
      ALLOCATE(ubasp(nbta_asv))
c
      ALLOCATE(tempcu(ntpr_asv))
c
      ALLOCATE(aprehw(narx_asv,ntpr_asv))
      ALLOCATE(apresg(narx_asv,ntpr_asv))
c
      ALLOCATE(axhfe(narx_asv,ntpr_asv))
      ALLOCATE(axlke(narx_asv,ntpr_asv))
      ALLOCATE(axvfe(narx_asv,ntpr_asv))
c
      ALLOCATE(aadh(narx_asv,ntpr_asv))
      ALLOCATE(aadhh(narx_asv,ntpr_asv))
      ALLOCATE(aadhv(narx_asv,ntpr_asv))
      ALLOCATE(aaphi(narx_asv,ntpr_asv))
      ALLOCATE(abdh(narx_asv,ntpr_asv))
      ALLOCATE(abdhh(narx_asv,ntpr_asv))
      ALLOCATE(abdhv(narx_asv,ntpr_asv))
      ALLOCATE(abdot(narx_asv,ntpr_asv))
      ALLOCATE(abdoth(narx_asv,ntpr_asv))
      ALLOCATE(abdotv(narx_asv,ntpr_asv))
c
      ALLOCATE(adadhh(narx_asv,ntpr_asv,ipch_asv))
      ALLOCATE(adadhv(narx_asv,ntpr_asv,ipcv_asv))
      ALLOCATE(adbdhh(narx_asv,ntpr_asv,ipch_asv))
      ALLOCATE(adbdhv(narx_asv,ntpr_asv,ipcv_asv))
      ALLOCATE(adbdth(narx_asv,ntpr_asv,ipch_asv))
      ALLOCATE(adbdtv(narx_asv,ntpr_asv,ipcv_asv))
      ALLOCATE(adhfe(narx_asv,ntpr_asv,ipch_asv))
      ALLOCATE(advfe(narx_asv,ntpr_asv,ipcv_asv))
c
      ALLOCATE(axhfsa(narx_asv,ntpr_asv,nsta_asv))
      ALLOCATE(axlksa(narx_asv,ntpr_asv,nsta_asv))
      ALLOCATE(axvfsa(narx_asv,ntpr_asv,nsta_asv))
c
      ALLOCATE(adhfsa(narx_asv,ntpr_asv,ipch_asv,nsta_asv))
      ALLOCATE(advfsa(narx_asv,ntpr_asv,ipcv_asv,nsta_asv))
c
      ALLOCATE(atwta(ncta_asv))
      ALLOCATE(cdrsa(ndrsa_asv))
      ALLOCATE(cessa(nessa_asv))
      ALLOCATE(mwtspa(nsta_asv))
      ALLOCATE(vosp0a(nsta_asv))
      ALLOCATE(zchara(nsta_asv))
c
      ALLOCATE(apxa(iapxa_asv,nxta_asv))
      ALLOCATE(bpxa(ibpxa_asv,nxta_asv))
c
c     SEDH data.
c
      ALLOCATE(azeroa(nata_asv))
      ALLOCATE(insgfa(nata_asv))
c
c     Pitzer data.
c
      ALLOCATE(nalpaa(nslta_asv))
      ALLOCATE(nmuxa(3,nmuta_asv))
      ALLOCATE(nslxa(2,nslta_asv))
c
      ALLOCATE(amua(jpfc_asv,nmuta_asv))
c
      ALLOCATE(aslma(jpfc_asv,0:ipbt_asv,nslta_asv))
      ALLOCATE(palpaa(ipbt_asv,napa_asv))
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the remainder of the data1 file. The image of the data file
c     is stored in arrays and variables that typically end in 'a.'
c     A subset of this is the reaction data. The variable nbta and
c     the arrays nbaspa, cdrsa, ndrsa, ndrsra, and axlksa comprise
c     what is called the 'a' set of reaction data. These are not
c     modified by any code manipulations.
c
      call indata(aadh,aadhh,aadhv,aaphi,abdh,abdhh,abdhv,
     $ abdot,abdoth,abdotv,adadhh,adadhv,adbdhh,adbdhv,adbdth,
     $ adbdtv,adhfe,adhfsa,advfe,advfsa,amua,aprehw,apresg,
     $ apxa,aslma,atwta,axhfe,axhfsa,axlke,axlksa,axvfe,axvfsa,
     $ azeroa,bpxa,cco2,cdrsa,cessa,eps100,iapxa_asv,iapxta,
     $ iaqsla,ibpxa_asv,ibpxta,ielam,igas,ikta_asv,insgfa,ipbt_asv,
     $ ipch,ipch_asv,ipcv,ipcv_asv,ixrn1a,ixrn2a,jpdblo,jpfc_asv,
     $ jptffl,jsola,mwtspa,nad1,nalpaa,napa_asv,napta,narn1a,narn2a,
     $ narx_asv,narxt,nata,nata_asv,nbaspa,nbta,nbta_asv,nbtafd,
     $ nbta1_asv,ncmpra,ncta,ncta_asv,ndrsa,ndrsa_asv,ndrsra,nessa,
     $ nessa_asv,nessra,ngrn1a,ngrn2a,ngta,ngta_asv,nlrn1a,nlrn2a,
     $ nlta,nlta_asv,nmrn1a,nmrn2a,nmta,nmta_asv,nmuta,nmuta_asv,
     $ nmuxa,noutpt,npta,npta_asv,nslta,nslta_asv,nslxa,nsta,
     $ nsta_asv,ntid_asv,ntitld,ntpr_asv,ntprt,nttyo,nxrn1a,nxrn2a,
     $ nxta,nxta_asv,qclnsa,palpaa,tdamax,tdamin,tempcu,ubasp,
     $ udakey,udatfi,uelema,uspeca,uphasa,uptypa,utitld,vosp0a,
     $ zchara)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Allocate arrays to store the data extracted from the pseudo-data
c     files.
c
      ntf1a_asv = ntf1a_par
      ALLOCATE(ntf1a(ntf1a_asv))
      ALLOCATE(tf1a(ntf1a_asv))
c
      ntf2a_asv = ntf2a_par
      ALLOCATE(ntf2a(ntf2a_asv))
      ALLOCATE(tf2a(ntf2a_asv))
c
      ntf1mx = ntf1a_asv
      ntf2mx = ntf2a_asv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the coefficients for calculating the HCO3-CO3-OH total
c     alkalinity. These are actually set in a pseudo-data file in the
c     EQLIB INCLUDE file eqltf1.h. The subroutine called here assigns
c     the coefficients to the species read from the data file.
c
c     Calling sequence substitutions:
c       ntf1a for ntfxa
c       ntf1mx for ntfxmx
c       ntf1ta for ntfxta
c       tf1a for tfxa
c       utf1xd for utfxxd
c
      call inttfx(narn1a,narn2a,noutpt,nsta_asv,ntf1a,ntf1mx,
     $ ntf1ta,nttyo,tf1a,uspeca,utf1xd)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the coefficients for calculating the extended total
c     alkalinity. These are actually set in a pseudo-data file in the
c     EQLIB INCLUDE file eqltf2.h. The subroutine called here assigns
c     the coefficients to the species read from the data file.
c
c     Calling sequence substitutions:
c       ntf2a for ntfxa
c       ntf2mx for ntfxmx
c       ntf2ta for ntfxta
c       tf2a for tfxa
c       utf2xd for utfxxd
c
      call inttfx(narn1a,narn2a,noutpt,nsta_asv,ntf2a,ntf2mx,
     $ ntf2ta,nttyo,tf2a,uspeca,utf2xd)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the indices of H+, OH-, Cl-, aqueous O2(g), and aqueous e-.
c     This will be repeated after data compression.
c
c     Calling sequence substitutions:
c       narn1a for narn1
c       narn2a for narn2
c       nchloa for nchlor
c       nhydra for nhydr
c       neleca for nelect
c       nhydxa for nhydx
c       no2gaa for no2gaq
c       nsta_asv for nstmax
c       uspeca for uspec
c
      call gspion(narn1a,narn2a,nchloa,neleca,nhydra,nhydxa,noutpt,
     $ no2gaa,nsta_asv,nttyo,uspeca)
c
c     Calling sequence substitutions:
c       nbaspa for nbasp
c       nbta_asv for nbtmax
c       nbta for nbt
c       ncta for nct
c       ndrsra for ndrsr
c       nrdxsa for nrdxsp
c       nsta_asv for nstmax
c       uspeca for uspec
c
      call grdxsp(nbaspa,nbta,nbta_asv,ncta,ndrsra,noutpt,
     $ nrdxsa,nsta_asv,nttyo,uspeca)
c
c     Get the basis index of water.
c
c     Calling sequence substitutions:
c       nbaspa for nbasp
c       nbta for nbt
c       nbta_asv for nbtmax
c       narn1a for ns
c
      nbwa = nbasis(nbaspa,nbta,nbta_asv,narn1a)
c
c     Get some constants that depend on the molecular weight of
c     water.
c
      omega = 1000./mwtspa(narn1a)
      omeglg = log10(omega)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Setup up the remaining array allocation variables required to
c     allocate variables needed to store the dta to be read from the
c     input file.
c
c     The following value of k_asv is special to EQ3NR.
c     A larger value is necessary for EQ6.
c
      k_asv = nbta_asv + 5
c
c     The following is unique to EQ3NR.
c
      nxic_asv = ikta_asv*nxti_asv
      nxicmx = nxic_asv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Allocate arrays to store the data to be read from the input
c     file.
c
      ALLOCATE(jflgi(nbta_asv))
      ALLOCATE(nbaspi(nbta_asv))
c
      ALLOCATE(npnxp(nxti_asv))
c
      ALLOCATE(ncmpri(2,nxti_asv))
c
      ALLOCATE(uspeci(nbta_asv))
      ALLOCATE(umemi(nxic_asv))
      ALLOCATE(usoli(nxti_asv))
c
      ALLOCATE(covali(nbta_asv))
c
      ALLOCATE(mtbaqi(nbta_asv))
      ALLOCATE(mtbi(nbta_asv))
      ALLOCATE(xbari(nxic_asv))
      ALLOCATE(zvclgi(k_asv))
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Setup up the remaining required array allocation variables.
c     At the present time, just set the regular variables which
c     correspond to specific variables required for the 'a' set
c     of data (e.g., nbt_asv to nbta_asv) to the values of the
c     corresponding 'a' set variables.
c
c     Primary allocation variables.
c
      nat_asv = nata_asv
      nbt_asv = nbta_asv
      nct_asv = ncta_asv
      ngt_asv = ngta_asv
      nlt_asv = nlta_asv
      nmt_asv = nmta_asv
      npt_asv = npta_asv
      nst_asv = nsta_asv
      nxt_asv = nxta_asv
c
      nap_asv = napa_asv
      nslt_asv = nslta_asv
      nmut_asv = nmuta_asv
c
      ikt_asv = ikta_asv
      iapx_asv = iapxa_asv
      ibpx_asv = ibpxa_asv
c
      ntf1_asv = ntf1a_asv
      ntf2_asv = ntf2a_asv
c
c     Secondary allocation variables.
c
      nbt1_asv = nbt_asv + 1
      ndrs_asv = 7*nst_asv
      ness_asv = 5*nst_asv
      ntfx_asv = max(ntf1_asv,ntf2_asv)
      nmx_asv = 3*nmut_asv
      nsx_asv = 2*nslt_asv
c
c     Note: nazp_asv and nazm_asv could be better calculated after
c     compression of the set of aqueous species.
c
      azchmx = 0.
      do ns = narn1a,narn2a
        azch = abs(zchara(ns))
        azchmx = max(azch,azchmx)
      enddo
      nazp_asv = nint(azchmx)
      nazm_asv = -nazp_asv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize array dimension variables corresponding to allocation
c     size variables associated with run-time dimensioning.
c
      iapxmx = iapx_asv
      ibpxmx = ibpx_asv
      iktmax = ikt_asv
      ipchmx = ipch_asv
      ipcvmx = ipcv_asv
      ipbtmx = ipbt_asv
      jpfcmx = jpfc_asv
c
      kmax = k_asv
c
      napmax = nap_asv
      narxmx = narx_asv
      natmax = nat_asv
      nazpmx = nazp_asv
      nazmmx = nazm_asv
      nbtmax = nbt_asv
      nbt1mx = nbt1_asv
      nctmax = nct_asv
      ndrsmx = ndrs_asv
      nessmx = ness_asv
      ngtmax = ngt_asv
      nltmax = nlt_asv
      nmtmax = nmt_asv
      nmutmx = nmut_asv
      nptmax = npt_asv
      nsltmx = nslt_asv
      nstmax = nst_asv
      ntidmx = ntid_asv
      ntprmx = ntpr_asv
      nxtmax = nxt_asv
c
      nmxmax = nmx_asv
      nsxmax = nsx_asv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Allocate the rest of the necessary arrays.
c
      ALLOCATE(iindx1(k_asv))
      ALLOCATE(ipndx1(k_asv))
      ALLOCATE(ipivot(k_asv))
c
      ALLOCATE(insgf(nat_asv))
c
      ALLOCATE(ibswx(nbt_asv))
      ALLOCATE(ixbasp(nbt_asv))
      ALLOCATE(nbasp(nbt_asv))
      ALLOCATE(nbaspd(nbt_asv))
      ALLOCATE(nbaspx(nbt_asv))
      ALLOCATE(nbmap(nbt_asv))
      ALLOCATE(ncosp(nbt_asv))
      ALLOCATE(ndecsp(nbt_asv))
      ALLOCATE(nfac(nbt_asv))
c
      ALLOCATE(ncmap(nct_asv))
c
      ALLOCATE(iction(nbt_asv))
      ALLOCATE(kction(nbt_asv))
      ALLOCATE(jjndex(nbt_asv))
      ALLOCATE(kkndex(nbt_asv))
c
      ALLOCATE(ndrs(ndrs_asv))
      ALLOCATE(ndrsd(ndrs_asv))
      ALLOCATE(ndrsx(ndrs_asv))
c
      ALLOCATE(ness(ness_asv))
c
      ALLOCATE(igstak(ngt_asv))
      ALLOCATE(jgstak(ngt_asv))
      ALLOCATE(jgsort(ngt_asv))
c
      ALLOCATE(jpflag(npt_asv))
      ALLOCATE(ncmpr(2,npt_asv))
c
      ALLOCATE(istack(nst_asv))
      ALLOCATE(jstack(nst_asv))
      ALLOCATE(jcsort(nst_asv))
      ALLOCATE(jjsort(nst_asv))
      ALLOCATE(jssort(nst_asv))
      ALLOCATE(jflag(nst_asv))
      ALLOCATE(jflagd(nst_asv))
      ALLOCATE(jsflag(nst_asv))
      ALLOCATE(jsitex(nst_asv))
      ALLOCATE(nphasx(nst_asv))
      ALLOCATE(nsmap(nst_asv))
c
      ALLOCATE(ntfx(ntfx_asv))
      ALLOCATE(ntf1(ntf1_asv))
      ALLOCATE(ntf2(ntf2_asv))
c
      ALLOCATE(kxmod(nxmd_asv))
c
      ALLOCATE(iapxt(nxt_asv))
      ALLOCATE(ibpxt(nxt_asv))
      ALLOCATE(jsol(nxt_asv))
c
      ALLOCATE(ndrsr(2,nst_asv))
      ALLOCATE(ndrsrd(2,nst_asv))
      ALLOCATE(ndrsrx(2,nst_asv))
      ALLOCATE(nessr(2,nst_asv))
c
      ALLOCATE(qxknph(npt_asv))
c
      ALLOCATE(utitl(ntit_asv))
      ALLOCATE(utitl2(ntit_asv))
c
      ALLOCATE(ubmtbi(nbt_asv))
      ALLOCATE(ucospi(nbt_asv))
c
      ALLOCATE(uobsw(2,nbt_asv))
      ALLOCATE(usbsw(2,nbt_asv))
c
      ALLOCATE(uspec(nst_asv))
      ALLOCATE(uxmod(nxmd_asv))
      ALLOCATE(uzveci(k_asv))
      ALLOCATE(uzvec1(k_asv))
c
      ALLOCATE(uphase(npt_asv))
      ALLOCATE(uptype(npt_asv))
c
      ALLOCATE(uelem(nct_asv))
c
      ALLOCATE(uldel(k_asv))
      ALLOCATE(ulbeta(k_asv))
c
      ALLOCATE(apx(iapx_asv,nxt_asv))
      ALLOCATE(bpx(ibpx_asv,nxt_asv))
      ALLOCATE(wfac(ikt_asv,nxt_asv))
c
      ALLOCATE(aamatr(k_asv,k_asv))
      ALLOCATE(gmmatr(k_asv,k_asv))
c
      ALLOCATE(acflg(nst_asv))
      ALLOCATE(acflgo(nst_asv))
      ALLOCATE(act(nst_asv))
      ALLOCATE(actlg(nst_asv))
      ALLOCATE(affsd(nst_asv))
c
      ALLOCATE(affpd(npt_asv))
      ALLOCATE(loph(npt_asv))
      ALLOCATE(moph(npt_asv))
      ALLOCATE(sidrph(npt_asv))
c
      ALLOCATE(azero(nat_asv))
      ALLOCATE(a3bars(nat_asv))
c
      ALLOCATE(ahrc(nbt_asv))
      ALLOCATE(amtb(nbt_asv))
      ALLOCATE(cjbasp(nbt_asv))
c
      ALLOCATE(atwt(nct_asv))
      ALLOCATE(cteaq(nct_asv))
      ALLOCATE(mte(nct_asv))
      ALLOCATE(mteaq(nct_asv))
      ALLOCATE(ppmwe(nct_asv))
c
      ALLOCATE(alpha(k_asv))
      ALLOCATE(beta(k_asv))
      ALLOCATE(betao(k_asv))
      ALLOCATE(delvco(k_asv))
      ALLOCATE(delvec(k_asv))
      ALLOCATE(rhsvec(k_asv))
      ALLOCATE(zvclg1(k_asv))
      ALLOCATE(zvec1(k_asv))
c
      ALLOCATE(bfac(nbt_asv))
      ALLOCATE(coval(nbt_asv))
      ALLOCATE(ctb(nbt_asv))
      ALLOCATE(dlogxw(nbt_asv))
      ALLOCATE(efac(nbt_asv))
      ALLOCATE(ehrc(nbt_asv))
      ALLOCATE(fo2lrc(nbt_asv))
      ALLOCATE(mtb(nbt_asv))
      ALLOCATE(mtbaq(nbt_asv))
      ALLOCATE(perc(nbt_asv))
c
      ALLOCATE(cdrs(ndrs_asv))
      ALLOCATE(cdrsd(ndrs_asv))
      ALLOCATE(cdrsx(ndrs_asv))
      ALLOCATE(cdrtw(nst_asv))
      ALLOCATE(cdrw(nst_asv))
      ALLOCATE(cnufac(nst_asv))
      ALLOCATE(conc(nst_asv))
      ALLOCATE(conclg(nst_asv))
      ALLOCATE(losp(nst_asv))
      ALLOCATE(lsort(nst_asv))
      ALLOCATE(mosp(nst_asv))
      ALLOCATE(mwtsp(nst_asv))
      ALLOCATE(sidrsp(nst_asv))
      ALLOCATE(vosp0(nst_asv))
      ALLOCATE(weight(nst_asv))
      ALLOCATE(xbar(nst_asv))
      ALLOCATE(xbarlg(nst_asv))
      ALLOCATE(zchar(nst_asv))
      ALLOCATE(zchcu6(nst_asv))
      ALLOCATE(zchsq2(nst_asv))
c
      ALLOCATE(cess(ness_asv))
      ALLOCATE(tfx(ntfx_asv))
      ALLOCATE(tf1(ntf1_asv))
      ALLOCATE(tf2(ntf2_asv))
      ALLOCATE(xlkmod(nxmd_asv))
c
      ALLOCATE(fsort(ngt_asv))
      ALLOCATE(fugac(ngt_asv))
      ALLOCATE(fugalg(ngt_asv))
c
      ALLOCATE(xhfs(nst_asv))
      ALLOCATE(xhfsd(nst_asv))
      ALLOCATE(xlks(nst_asv))
      ALLOCATE(xlksd(nst_asv))
      ALLOCATE(xvfs(nst_asv))
      ALLOCATE(xvfsd(nst_asv))
c
      ALLOCATE(dhfe(ipch_asv))
      ALLOCATE(dvfe(ipcv_asv))
c
      ALLOCATE(dadhh(ipch_asv))
      ALLOCATE(dadhv(ipcv_asv))
      ALLOCATE(dbdhh(ipch_asv))
      ALLOCATE(dbdhv(ipcv_asv))
      ALLOCATE(dbdth(ipch_asv))
      ALLOCATE(dbdtv(ipcv_asv))
c
      ALLOCATE(dhfs(ipch_asv,nst_asv))
      ALLOCATE(dhfsd(ipch_asv,nst_asv))
      ALLOCATE(dvfs(ipcv_asv,nst_asv))
      ALLOCATE(dvfsd(ipcv_asv,nst_asv))
c
      ALLOCATE(axhfs(narx_asv,ntpr_asv,nst_asv))
      ALLOCATE(axhfsd(narx_asv,ntpr_asv,nst_asv))
      ALLOCATE(axhfsx(narx_asv,ntpr_asv,nst_asv))
      ALLOCATE(axlks(narx_asv,ntpr_asv,nst_asv))
      ALLOCATE(axlksd(narx_asv,ntpr_asv,nst_asv))
      ALLOCATE(axlksx(narx_asv,ntpr_asv,nst_asv))
      ALLOCATE(axvfs(narx_asv,ntpr_asv,nst_asv))
      ALLOCATE(axvfsd(narx_asv,ntpr_asv,nst_asv))
      ALLOCATE(axvfsx(narx_asv,ntpr_asv,nst_asv))
c
      ALLOCATE(adhfs(narx_asv,ntpr_asv,ipch_asv,nst_asv))
      ALLOCATE(adhfsd(narx_asv,ntpr_asv,ipch_asv,nst_asv))
      ALLOCATE(adhfsx(narx_asv,ntpr_asv,ipch_asv,nst_asv))
      ALLOCATE(advfs(narx_asv,ntpr_asv,ipcv_asv,nst_asv))
      ALLOCATE(advfsd(narx_asv,ntpr_asv,ipcv_asv,nst_asv))
      ALLOCATE(advfsx(narx_asv,ntpr_asv,ipcv_asv,nst_asv))
c
c     Pitzer data.
c
      ALLOCATE(nalpha(nslt_asv))
c
      ALLOCATE(nmux(3,nmut_asv))
      ALLOCATE(nmxi(2,nat_asv))
      ALLOCATE(nmxx(3,nmx_asv))
      ALLOCATE(nslx(2,nslt_asv))
      ALLOCATE(nsxi(2,nat_asv))
      ALLOCATE(nsxx(2,nsx_asv))
c
      ALLOCATE(amu(jpfc_asv,nmut_asv))
      ALLOCATE(pmu(nmut_asv))
c
      ALLOCATE(aslm(jpfc_asv,0:ipbtmx,nslt_asv))
      ALLOCATE(pslm(nslt_asv))
      ALLOCATE(dpslm(2,nslt_asv))
      ALLOCATE(pslamn(0:ipbtmx,nslt_asv))
c
      ALLOCATE(gpit(ipbtmx,nap_asv))
      ALLOCATE(dgpit(2,ipbtmx,nap_asv))
      ALLOCATE(palpha(ipbtmx,nap_asv))
c
      ALLOCATE(elam(nazp_asv,nazp_asv))
      ALLOCATE(delam(2,nazp_asv,nazp_asv))
      ALLOCATE(pelm(nazp_asv,nazp_asv))
      ALLOCATE(dpelm(2,nazp_asv,nazp_asv))
      ALLOCATE(selm(nazm_asv:nazp_asv))
      ALLOCATE(dselm(2,nazm_asv:nazp_asv))
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Allocate arrays needed to write a complete EQ6 input file as
c     the EQ3NR pickup file.
c
      ALLOCATE(jcode(nrct_asv))
      ALLOCATE(jreac(nrct_asv))
      ALLOCATE(nsk(nrct_asv))
c
      ALLOCATE(imech(2,nrct_asv))
      ALLOCATE(nrk(2,nrct_asv))
c
      ALLOCATE(iact(imch_asv,2,nrct_asv))
      ALLOCATE(ndact(imch_asv,2,nrct_asv))
c
      ALLOCATE(ibsrti(nsrt_asv))
      ALLOCATE(iesrti(nsrt_asv))
c
      ALLOCATE(ixrti(nxrt_asv))
c
      ALLOCATE(utitl1(ntit_asv))
c
      ALLOCATE(uprspi(nprs_asv))
c
      ALLOCATE(uprphi(nprp_asv))
      ALLOCATE(ubsri(nbt1_asv,nsrt_asv))
      ALLOCATE(ucxri(ikt_asv,nxrt_asv))
      ALLOCATE(udac(ndct_asv,imch_asv,2,nrct_asv))
      ALLOCATE(uffg(nffg_asv))
      ALLOCATE(ureac(nrct_asv))
      ALLOCATE(uxcat(nxop_asv))
      ALLOCATE(uxopt(nxop_asv))
      ALLOCATE(uxopex(nxpe_asv))
      ALLOCATE(uesri(nct_asv,nsrt_asv))
c
      ALLOCATE(cbsri(nbt1_asv,nsrt_asv))
      ALLOCATE(cdac(ndct_asv,imch_asv,2,nrct_asv))
      ALLOCATE(cesri(nct_asv,nsrt_asv))
      ALLOCATE(csigma(imch_asv,2,nrct_asv))
      ALLOCATE(eact(imch_asv,2,nrct_asv))
      ALLOCATE(hact(imch_asv,2,nrct_asv))
      ALLOCATE(rkb(imch_asv,2,nrct_asv))
      ALLOCATE(trkb(imch_asv,2,nrct_asv))
      ALLOCATE(fkrc(nrct_asv))
      ALLOCATE(modr(nrct_asv))
      ALLOCATE(morr(nrct_asv))
      ALLOCATE(sfcar(nrct_asv))
      ALLOCATE(ssfcar(nrct_asv))
      ALLOCATE(vreac(nrct_asv))
      ALLOCATE(moffg(nffg_asv))
      ALLOCATE(xlkffg(nffg_asv))
      ALLOCATE(mprphi(nprp_asv))
      ALLOCATE(mprspi(nprs_asv))
      ALLOCATE(ptk(nptk_asv))
      ALLOCATE(ttk(nttk_asv))
      ALLOCATE(rxbari(ikt_asv,nxrt_asv))
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Allocate some arrays needed to assist in writing a complete
c     EQ6 input file as the EQ3NR pickup file. These arrays do not
c     appear on that file itself.
c
      ALLOCATE(uesr1(nct_asv))
      ALLOCATE(cesr1(nct_asv))
      ALLOCATE(cbsr1(nbt1_asv))
      ALLOCATE(ubsr1(nbt1_asv))
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Allocate arrays to handle mass balance coefficients.
c
      nsts_asv = 10*nst_asv
      nstsmx = nsts_asv
c
      ALLOCATE(nsts(nsts_asv))
      ALLOCATE(nstsr(2,nst_asv))
      ALLOCATE(csts(nsts_asv))
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize the problem counter.
c
      nprob = 0
c
c     The label below is a return point for processing an additional
c     problem specifed on the input file.
c
   20 nprob = nprob + 1
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set special species indices.
c
      nchlor = nchloa
      nhydr = nhydra
      nhydx = nhydxa
      no2gaq = no2gaa
      nelect = neleca
      nrdxsp = nrdxsa
      nbw = nbwa
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize certain arrays to zeros, blanks, etc.
c
      qbassw = .false.
c
      call initiz(iopt,noptmx)
      call initiz(iopg,nopgmx)
      call initiz(iopr,noprmx)
      call initiz(iodb,nodbmx)
c
      call initiz(jflgi,nbtmax)
      call initaz(covali,nbtmax)
c
      call initcb(uspeci,nbtmax)
      call initcb(ucospi,nbtmax)
c
      nmax = 2*nbtmax
      call initcb(uobsw,nmax)
      call initcb(usbsw,nmax)
c
      call initiz(ncosp,nbtmax)
      call initiz(ndecsp,nbtmax)
      call initaz(coval,nbtmax)
c
      call initiz(nbasp,nbtmax)
c
      call initiz(jflag,nstmax)
      call initiz(jflagd,nstmax)
      call initcb(uspec,nstmax)
      call initaz(xlks,nstmax)
      call initaz(zchar,nstmax)
      call initaz(zchsq2,nstmax)
      call initaz(zchcu6,nstmax)
c
cxxxxxxxxxxx
      call initaz(cdrs,ndrsmx)
      call initaz(cdrsx,ndrsmx)
      call initiz(ndrs,ndrsmx)
      call initiz(ndrsx,ndrsmx)
c
cxxxxxxxxxxx
      nmax = 2*nstmax
      call initiz(ndrsr,nmax)
      call initiz(ndrsrx,nmax)
c
      nmax = narxmx*ntprmx*nstmax
      call initaz(axlks,nmax)
      call initaz(axhfs,nmax)
      call initaz(axvfs,nmax)
cxxxxxxxxxxx
      call initaz(axlksx,nmax)
      call initaz(axhfsx,nmax)
      call initaz(axvfsx,nmax)
c
      nmax = narxmx*ntprmx*ipchmx*nstmax
      call initaz(adhfs,nmax)
      call initaz(adhfsx,nmax)
c
      nmax = narxmx*ntprmx*ipcvmx*nstmax
      call initaz(advfs,nmax)
      call initaz(advfsx,nmax)
cxxxxxxxxxxx
c
      call initcb(ugexp,netmax)
      call initaz(cgexpi,netmax)
      call initiz(jgext,netmax)
c
      nmax = jetmax*netmax
      call initiz(ngexrt,nmax)
c
      nmax = ietmax*jetmax*netmax
      call initaz(xgexsi,nmax)
c
      nmax = iapxmx*nxtmax
      call initaz(apx,nmax)
c
      nmax = ibpxmx*nxtmax
      call initaz(bpx,nmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Determine the problem input format.
c
      read (ninpts,1030,end=105,err=107) ux8
 1030 format(a8)
      backspace ninpts
      uinfor = 'W'
      if (ux8(1:8) .eq. '|-------') uinfor = 'D'
c
c     Read the problem input.
c
      if (uinfor(1:1) .eq. 'W') then
c
c       Compact (W) format.
c
        call rd3inw(cgexj,cgexpi,covali,ehi,egexsi,fo2lgi,iebal3,
     $  ietmax,iodb,iopg,iopr,iopt,irdxc3,itdsf3,itermx,jgext,jetmax,
     $  jflgi,jgexti,jpres3,kxmod,mwtges,nbti,nbtmax,ncmpri,net,neti,
     $  netmax,ngexti,ninpts,ngexrt,nobswt,nodbmx,nopgmx,noprmx,
     $  noptmx,noutpt,nprob,nsbswt,ntitl,ntitmx,nttyo,nxicmx,nxmdmx,
     $  nxmod,nxti,nxtimx,pei,press,qend,qgexsh,qrderr,rho,scamas,
     $  tdspkg,tdspl,tempc,tgexp,tolbt,toldl,tolspf,ucospi,uebal,
     $  ugexj,ugexji,ugexmo,ugexp,ugexpi,ugexr,ugexsi,umemi,uobsw,
     $  uredox,usbsw,usoli,uspeci,utitl,uhfgex,uvfgex,uxkgex,uxmod,
     $  xbari,xgexsi,xhfgex,xlkgex,xvfgex,xlkmod,zgexj)
      else
c
c       Menu-style (D) format.
c
        call rd3ind(cgexj,cgexpi,covali,ehi,egexsi,fo2lgi,iebal3,
     $  ietmax,iodb,iopg,iopr,iopt,irdxc3,itdsf3,itermx,jgext,jetmax,
     $  jflgi,jgexti,jpres3,kxmod,mwtges,nbti,nbtmax,ncmpri,net,neti,
     $  netmax,ngexti,ninpts,ngexrt,nobswt,nodbmx,nopgmx,noprmx,
     $  noptmx,noutpt,nprob,nsbswt,ntitl,ntitmx,nttyo,nxicmx,nxmdmx,
     $  nxmod,nxti,nxtimx,pei,press,qend,qgexsh,qrderr,rho,scamas,
     $  tdspkg,tdspl,tempc,tgexp,tolbt,toldl,tolspf,ucospi,uebal,
     $  ugexj,ugexji,ugexmo,ugexp,ugexpi,ugexr,ugexsi,umemi,uobsw,
     $  uredox,usbsw,usoli,uspeci,utitl,uhfgex,uvfgex,uxkgex,uxmod,
     $  xbari,xgexsi,xhfgex,xlkgex,xvfgex,xlkmod,zgexj)
      endif
c
      go to 109
  105 qend = .true.
      go to 109
  107 qrderr = .true.
  109 continue
c
      if (qrderr) then
        write (noutpt,1203)
        write (nttyo,1203)
 1203   format(/' * Error - (EQ3NR/eq3nr) Encountered an error while',
     $  ' reading the input file.',/7x,'The line associated with',
     $  ' the problem should be the last one',/7x,'echoed on the',
     $  ' output file or the first line immediately thereafter.')
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qend) then
c
        write (noutpt,1205)
        write (nttyo,1205)
 1205   format(/' No further input found.',/)
c
c       Get end time and date. Also get the run, user, and cpu times.
c
        call runtim(iexec0,jexec0,texec0,noutpt,nttyo,trun,
     $  tuser,tcpu,udate1,utime1)
c
        j2 = ilnobl(udate0)
        j3 = ilnobl(udate1)
        write (noutpt,1208) utime0,udate0(1:j2),utime1,udate1(1:j3)
        write (nttyo,1208) utime0,udate0(1:j2),utime1,udate1(1:j3)
 1208   format(/10x,'Start time = ',a8,2x,a,/12x,
     $  'End time = ',a8,2x,a)
c
c       Print the run, user, and cpu times.
c
        write (noutpt,1020) trun
        write (nttyo,1020)  trun
 1020   format(/10x,' Run time = ',g10.3,' seconds')
        if (tuser .gt. 0.) then
          write (noutpt,1022) tuser
          write (nttyo,1022)  tuser
 1022     format(10x,'User time = ',g10.3,' seconds')
        endif
        if (tcpu .gt. 0.) then
          write (noutpt,1024) tcpu
          write (nttyo,1024)  tcpu
 1024     format(10x,' Cpu time = ',g10.3,' seconds')
        endif
c
        write (noutpt,1090)
        write (nttyo,1090)
 1090   format(/' Normal exit')
c
c       Clear the IEEE flag for floating-point underflow, if such a
c       flag is present, to avoid getting an unnecessary system
c       warning message. Underflow is a normal condition in EQ3/6.
c       Make porting changes in the EQLIBU subroutine that is called
c       in this section. Do not make the porting changes here.
c
        call cliefu()
c
c       Close and delete the stripped input file.
c
        close (ninpts,status='delete')
c
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Scan the input file title for option strings.
c
c     Scan the input file title for the USEOLDPITZERMU string.
c     If found and if iopg(1) = 1 (use Pitzer's equations), set up
c     set to calculate activity coefficients using recalculated
c     observable third-order parameters (Cphi, psi, zeta) instead
c     of their conventional mu equivalents.
c
      qhawep = .true.
      if (iopg(1) .eq. 1) then
        do n = 1,ntitl
          ux80 = utitl(n)
          call locase(ux80)
          i = index(ux80,'useoldpitzermu')
          if (i .gt. 0) then
            qhawep = .false.
            go to 106
          endif
        enddo
  106   continue
      endif
c
      if (.not.qhawep) then
        write (noutpt,1190)
        write (nttyo,1190)
 1190   format(/' * Note - (EQ3NR/eq3nr) Found the USEOLDPITZERMU',
     $  ' option string',/7x,'in one of the input file titles.',
     $  ' Will evaluate the Pitzer',/7x,'equations using the',
     $  ' lambda-mu format. Implied psi coefficients',/7x,'not',
     $  ' on the supporting data file will not be effectively',
     $  ' treated',/7x,'as having zero values. Rather, the',
     $  ' corresponding mu coefficients',/7x,'will be treated as',
     $  ' having zero values. In effect, an implied',/7x,'psi',
     $  ' coefficient value will include contributions from',
     $  ' the',/7x,'two cognate Cphi coefficients that appear',
     $  ' in the psi-mu',/7x,'breakdown equation.')
      endif
c
c     Scan the input file title for the USEOLDPITZER75 string.
c     If found and if iopg(1) = 1 (use Pitzer's equations), use
c     the older approximation given by Pitzer (1975) for calculating
c     the higher-order electrostatic terms.
c
      qpit75 = .false.
      if (iopg(1) .eq. 1) then
        do n = 1,ntitl
          ux80 = utitl(n)
          call locase(ux80)
          i = index(ux80,'useoldpitzer75')
          if (i .gt. 0) then
            qpit75 = .true.
            go to 170
          endif
        enddo
  170   continue
      endif
c
      if (qpit75) then
        write (noutpt,1193)
        write (nttyo,1193)
 1193   format(/' * Warning - (EQ3NR/eq3nr) Found the USEOLDPITZER75',
     $  ' option string',/7x,'in one of the input file titles.',
     $  ' Will evaluate the Pitzer',/7x,'equations using the',
     $  ' old Pitzer (1975) approximation for',/7x,'higher order',
     $  ' electrostatic terms, not the later Harvie (1981)',
     $  /7x,'approximation that is now used nearly universally',
     $  ' to evaluate',/7x,'the Pitzer equations.')
      endif
c
c     Scan the input file title for the WRITEPITZERJTABLES string.
c     If found and if iopg(1) = 1 (use Pitzer's equations), calculate
c     and output tables of the J(x) and J'(x) functions. This is a
c     one-time only option. It is not carried forward on the PICKUP
c     file.
c
      qcwrpj = .false.
      if (iopg(1) .eq. 1) then
        do n = 1,ntitl
          ux80 = utitl(n)
          call locase(ux80)
          i = index(ux80,'writepitzerjtables')
          if (i .gt. 0) then
            qcwrpj = .true.
            go to 180
          endif
        enddo
  180   continue
      endif
c
      if (qcwrpj) then
        write (noutpt,1194)
        write (nttyo,1194)
 1194   format(/' * Note - (EQ3NR/eq3nr) Found the WRITEPITZERJTABLES',
     $  ' option string',/7x,'in the input file title.',
     $  ' Will calculate the Pitzer',/7x,"J(x) and J'(x) functions",
     $  ' for higher order electrostatic terms and',/7x,'write output',
     $  ' tables for both the Pitzer (1975) and Harvie (1981)',/7x,
     $  'approximations.')
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     If the phase part of certain species names read from the
c     input file is blank, make this part 'Aqueous solution'.
c
      do n = 1,nsbswt
        if (usbsw(1,n)(25:48) .eq. ublk24(1:24))
     $  usbsw(1,n)(25:48) = uaqsln(1:24)
        if (usbsw(2,n)(25:48) .eq. ublk24(1:24))
     $  usbsw(2,n)(25:48) = uaqsln(1:24)
      enddo
c
      do nbi = 1,nbti
        if (uspeci(nbi)(25:48) .eq. ublk24(1:24))
     $  uspeci(nbi)(25:48) = uaqsln(1:24)
      enddo
c
      do n = 1,nobswt
        if (uobsw(1,n)(25:48) .eq. ublk24(1:24))
     $  uobsw(1,n)(25:48) = uaqsln(1:24)
        if (uobsw(2,n)(25:48) .eq. ublk24(1:24))
     $  uobsw(2,n)(25:48) = uaqsln(1:24)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set the redox variables.
c
      qredox = .true.
c
      if (irdxc3 .ne. -3) then
        do nbi = 1,nbti
          if (uspeci(nbi)(1:6).eq.'O2(g) ' .and.
     $      uspeci(nbi)(25:48).eq.uaqsln(1:24)) then
            write (noutpt,1100)
            write (nttyo,1100)
 1100       format(/' * Note - (EQ3NR/eq3nr) The input line for O2(g)',
     $      /7x,'will be ignored because irdxc3 is not equal to -3.')
            jflgi(nbi) = 0
            covali(nbi) = 0.
            go to 102
          endif
        enddo
  102   continue
      endif
c
      fo2lg = -99999.
      fo2 = 0.
      if (irdxc3 .eq. 0) then
        fo2lg = fo2lgi
        fo2 = texp(fo2lg)
      elseif (irdxc3 .eq. -1) then
        eh = ehi
      elseif (irdxc3 .eq. -2) then
        pe = pei
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Normalize the mole fractions of any solid solution components.
c
      nerr = 0
      do nxi = 1,nxti
        nr1 = ncmpri(1,nxi)
        nr2 = ncmpri(2,nxi)
c
        xx = 0.
        do nxic = nr1,nr2
          xx = xx + xbari(nxic)
        enddo
c
        if (xx .le. 0.) then
          j2 = ilnobl(usoli(nxi))
          write (noutpt,1110) usoli(nxi)(1:j2)
          write (nttyo,1110) usoli(nxi)(1:j2)
 1110     format(/' * Error - (EQ3NR/eq3nr) The solid solution',
     $    a,/7x,'was given on the input file with a null composition.')
          nerr = nerr + 1
          go to 104
        endif
c
        do nxic = nr1,nr2
          xbari(nxic) = xbari(nxic)/xx
        enddo
  104   continue
      enddo
      if (nerr .gt. 0) go to 20
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check consistency between the activity coefficient option and the
c     data1 file.
c
      call cdakey(iopg,nopgmx,noutpt,nttyo,udakey,udatfi)
c
c     Get the name of the option for calculating the activity
c     coefficients of aqueous species.
c
      call nactop(iopg,nopgmx,noutpt,nttyo,uactop)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Open the pickup file if this is to be used.
c
      if (iopt(17) .ne. -1) then
        inquire(file=pfile,opened=qop)
        if (.not.qop) call openou(noutpt,nttyo,pfile,'formatted',
     $  nrecl,newin)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Copy the 'a' set of reaction data into the 'd' set. The 'a'
c     set remains a faithful image of the data file. The 'd' set
c     at this stage is subject to modification. The number of basis
c     species may be increased in response to input file directives.
c     Thus, nbtd may become greater than nbta, and the basis species
c     pointer array nbaspd lengthened with respect to nbaspa. The
c     'd' set at this stage may also be modified by directives on the
c     input file to execute special basis switching. This redefines
c     which of the basis species are the strict basis species, and
c     according causes reactions and attendant thermodynamic data to
c     be rewritten. Modifications of the 'd' set at this stage
c     support options affecting problem definition. After the data
c     have been compressed to eliminate unnecessary species, the 'd'
c     set will be redefined. This second stage form is equivalent to
c     a compressed version of the first stage form with additional
c     input file directed modifications.
c
c     Calling sequence substitutions:
c       adhfsa for adhfs
c       advfsa for advfs
c       axhfsa for axhfs
c       axlksa for axlks
c       axvfsa for axvfs
c       cdrsa for cdrs
c       nbaspa for nbasp
c       ndrsa for ndrs
c       ndrsra for ndrsr
c
      call cdrssd(adhfsa,adhfsd,advfsa,advfsd,axhfsa,axhfsd,
     $ axlksa,axlksd,axvfsa,axvfsd,cdrsa,cdrsd,ipch,ipchmx,ipcv,ipcvmx,
     $ narxmx,nbaspa,nbaspd,nbtmax,ndrsa,ndrsd,ndrsmx,ndrsra,ndrsrd,
     $ nstmax,ntprmx)
      nbtd = nbta
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Execute special basis switching to redefine which basis species
c     are the strict basis species in the 'd' set of reaction data
c     (nbtd/nbaspd/cdrsd/ndrsd/ndrsrd/axlksd). This type of basis
c     switching is executed to assist in problem definition. Ordinary
c     basis switching, executed later in this code, is done only
c     for numerical purposes.
c
      if (nsbswt .gt. 0) then
c
        do nsbsw = 1,nsbswt
c
c         Interpret the directive for the nsbsw-th switch.
c
          call intsbs(nb1,nb2,nbaspd,nbtd,nbtmax,noutpt,ns1,ns2,
     $    nsbsw,nstmax,nttyo,usbsw,uspeca)
c
c         Execute the switch.
c
          call swtchb(adhfsd,adhfsx,advfsd,advfsx,axhfsd,axhfsx,
     $    axlksd,axlksx,axvfsd,axvfsx,cdrsd,cdrsx,ipch,ipchmx,ipcv,
     $    ipcvmx,narxmx,nbaspd,nbtmax,nbw,nb1,nb2,ndrsd,ndrsmx,ndrsx,
     $    ndrsrd,ndrsrx,noutpt,ns1,ns2,nsta,nstmax,ntprmx,nttyo,uspeca)
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Force jflag to default to 27 instead of 30 for the species
c     O2(aq) and H2(aq).
c
      call jfloha(jflgi,nbtd,nbti,nbtmax,noutpt,nttyo,
     $ nxmdmx,nxmod,ubasp,uspeci,uxmod)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Interpret the basis species information listed on the input file.
c
      call intbs3(covali,ier,jflag,jflgi,narn1a,narn2a,nbaspd,
     $ nbtd,nbti,nbtmax,ndrsrd,ndecsp,noutpt,nrdxsp,nsta,nstmax,
     $ nttyo,uspeca,uspeci)
c
      if (ier .gt. 0) go to 20
c
c     Set jflag to -1 for species that can't appear in the system.
c
      call jflaux(jflag,nbaspd,nbtd,nbtmax,ndrsd,ndrsmx,
     $ ndrsrd,nstmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set up temperature variables.
c
      tempk = tempc + 273.15
      call gntpr(ntpr,ntprmx,ntprt,tempc,tempcu)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set species status flags.
c
      call flgset(axlksd,iopt,jflag,jpflag,jsflag,kxmod,narn1a,
     $ narn2a,narxmx,nbaspd,nbtd,nbtmax,ncmpra,ncta,ndrsd,ndrsmx,
     $ ndrsrd,noptmx,noutpt,npta,nptmax,nrdxsp,nsta,nstmax,ntpr,
     $ ntprmx,nttyo,nxmdmx,nxmod,uphasa,uptypa,uspeca,uxmod)
c
c     Check the jpflag array to make sure it is consistent with the
c     jsflag array.
c
      call flgchk(jpflag,jsflag,ncmpra,npta,nptmax,nstmax,qclnsa)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Look at each active auxiliary basis species. Write a warning if
c     any other species in the corresponding dissociation reaction is
c     not present in the model.
c
      call bspchk(jsflag,nbaspd,nbtd,nbtmax,ndrsd,ndrsmx,ndrsrd,
     $ noutpt,nrdxsp,nstmax,nttyo,uspeca)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Do data array compression. Write working data arrays that
c     do not include phases and species that are not necessary
c     for the current problem.
c
      call cmpdat(adhfs,adhfsd,advfs,advfsd,amu,amua,apx,
     $ apxa,aslm,aslma,atwt,atwta,axhfs,axhfsd,axlks,axlksd,axvfs,
     $ axvfsd,azero,azeroa,bpx,bpxa,cdrs,cdrsd,cess,cessa,iapxmx,
     $ iapxt,iapxta,iaqsla,iaqsln,ibpxmx,ibpxt,ibpxta,ilrn1,ilrn2,
     $ imrn1,imrn2,ipbtmx,ipch,ipchmx,ipcv,ipcvmx,insgf,insgfa,
     $ iopg,ixrn1,ixrn1a,ixrn2,ixrn2a,jflag,jpfcmx,jpflag,jsflag,
     $ jsitex,jsol,jsola,mwtsp,mwtspa,nalpaa,nalpha,napmax,napt,
     $ napta,narn1,narn1a,narn2,narn2a,narxmx,nat,natmax,nbasp,
     $ nbaspd,nbmap,nbt,nbtd,nbti,nbtmax,nchlor,ncmap,ncmpr,ncmpra,
     $ nct,ncta,nctmax,ndecsp,ndrs,ndrsd,ndrsmx,ndrsr,ndrsrd,ness,
     $ nessa,nessmx,nessr,nessra,ngrn1,ngrn1a,ngrn2,ngrn2a,ngt,
     $ nlrn1,nlrn1a,nlrn2,nlrn2a,nlt,nmrn1,nmrn1a,nmrn2,nmrn2a,nmt,
     $ nmut,nmuta,nmutmx,nmux,nmuxa,nopgmx,nslt,nsltmx,nslta,nslx,
     $ nslxa,nphasx,npt,npta,nptmax,nsmap,nsta,nst,nstmax,ntf1,ntf1a,
     $ ntf1mx,ntf1t,ntf1ta,ntf2,ntf2a,ntf2mx,ntf2t,ntf2ta,ntprmx,
     $ nxrn1,nxrn1a,nxrn2,nxrn2a,nxt,nxtmax,palpaa,palpha,qchlor,
     $ tf1,tf1a,tf2,tf2a,uelem,uelema,uphasa,uphase,uptypa,uptype,
     $ uspec,uspeca,vosp0,vosp0a,zchar,zchara)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do ne = 1,netmax
        cgexp(ne) = 0.
      enddo
c
      do ns = 1,nstmax
        acflg(ns) = 0.
        acflgo(ns) = 0.
        act(ns) = 0.
        conc(ns) = 0.
        mosp(ns) = 0.
        xbar(ns) = 0.
        zchsq2(ns) = 0.
        zchcu6(ns) = 0.
      enddo
c
      av = -99999.
      call initav(conclg,nstmax,av)
      call initav(losp,nstmax,av)
      call initav(actlg,nstmax,av)
      call initav(xbarlg,nstmax,av)
c
      do np = 1,nptmax
        moph(np) = 0.
      enddo
c
      av = -99999.
      call initav(loph,nptmax,av)
c
      do kcol = 1,kmax
        zvec1(kcol) = 0.
      enddo
c
      av = -99999.
      call initav(zvclg1,kmax,av)
c
      uv = 'conc'
      call initcv(ulbeta,kmax,uv)
      call initcv(uldel,kmax,uv)
c
      do ne = 1,netmax
        kern1(ne) = 0
        kern2(ne) = 0
      enddo
c
      do ne = 1,netmax
        egexpa(ne) = 0.
        egexpc(ne) = 0.
      enddo
c
      nmax = jetmax*netmax
      call initiz(jern1,nmax)
      call initiz(jern2,nmax)
      call initiz(ngext,nmax)
      call initaz(egexjc,nmax)
      call initaz(egexjf,nmax)
      call initaz(mgext,nmax)
c
      nmax = ietmax*jetmax*netmax
      call initaz(cegexs,nmax)
      call initaz(cpgexs,nmax)
      call initaz(mrgexs,nmax)
      call initiz(ngexro,nmax)
      call initiz(ngexso,nmax)
      call initiz(ngexsa,nmax)
      call initaz(egexs,nmax)
c
      nmax = ietmax*jetmax*nertmx
      call initaz(egersi,nmax)
      call initaz(xgersi,nmax)
c
      nmax = ketmax*netmax
      call initiz(kgexsa,nmax)
      call initaz(egexw,nmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the indices of special species after compression.
c
c     Get the indices of H+, OH-, Cl-, fictive aqueous O2(g), and
c     fictive aqueous e-.
c
      call gspion(narn1,narn2,nchlor,nelect,nhydr,nhydx,noutpt,
     $ no2gaq,nstmax,nttyo,uspec)
c
c     Get the index of the redox basis species.
c
c     Calling sequence substitutions:
c       nbaspd for nbasp
c       nbtd for nbt
c       ndrsrd for ndrsr
c
      call grdxsp(nbaspd,nbtd,nbtmax,nct,ndrsrd,noutpt,
     $ nrdxsp,nstmax,nttyo,uspec)
c
c     Get the basis index of water (nbw).
c
c     Calling sequence substitutions:
c       nbaspd for nbasp
c       nbtd for nbt
c       narn1 for ns
c
      nbw = nbasis(nbaspd,nbtd,nbtmax,narn1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Look at each active auxiliary basis species. Change the
c     corresponding log K polynomial coefficients so that log K is
c     fixed at a value of -9999999. if any other species in the
c     corresponding dissociation reaction is not in the model.
c
      call bsplkp(axlks,narxmx,nbasp,nbt,nbtmax,ndrs,ndrsmx,
     $ ndrsr,nstmax,ntprmx)
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
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Load coefficients for interpreting the input alkalinity, if any.
c
      ntfxt = ntf1t
      ntfxmx = ntf1mx
      do n = 1,ntf1t
        ntfx(n) = ntf1(n)
        tfx(n) = tf1(n)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Interpret input file directives to create generic ion-exchange
c     phases and species.
c
      call intexi(al10,axhfs,axlks,axvfs,cegexs,cess,cdrs,
     $ cgexj,cpgexs,egexjf,iern1,iern2,ietmax,jern1,jern2,jetmax,
     $ jflag,jgext,jpflag,jsflag,jsitex,kern1,kern2,ketmax,kgexsa,
     $ mwtges,mwtsp,narn1,narn2,narxmx,narxt,nbasp,nbt,nbtmax,ncmpr,
     $ ndrs,ndrsmx,ndrsr,ness,nessmx,nessr,netmax,net,nern1,nern2,
     $ ngexro,ngexrt,ngexsa,ngexso,ngext,noutpt,nphasx,npt,nptmax,
     $ nst,nstmax,ntprt,ntprmx,nttyo,nvetmx,rconst,tgexp,ugexj,
     $ ugexmo,ugexmv,ugexp,ugexr,ugexs,ugexsr,uhfgex,uphase,uspec,
     $ uvfgex,uxkgex,xhfgex,xlkgex,xvfgex,zchar,zgexj)
c
      if (ier .gt. 0) go to 20
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (net .gt. 0) then
c
c       Set jflag values for the exchanger species.
c
        do ne = 1,net
          do je = 1,jgext(ne)
            nrr1 = jern1(je,ne)
            nrr2 = jern2(je,ne)
            do ns = nrr1,nrr2
c
c             Calling sequence substitutions:
c               nbaspd for nbasp
c
              nb = nbasis(nbasp,nbt,nbtmax,ns)
              if (nb .eq. 0) jflag(ns) = 30
            enddo
          enddo
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Alter any log K values as directed by the input file.
c
      if (nxmod .gt. 0) then
        rtcnst = 0.001*rconst*tempk
        afcnst = al10*rtcnst
        call alters(afcnst,apresg,axlks,cdrs,kxmod,narxmx,narxt,
     $  ndrs,ndrsmx,ndrsr,noutpt,npt,nptmax,nst,nstmax,ntpr,ntprmx,
     $  nttyo,nxmdmx,nxmod,tempc,uphase,uspec,uxmod,xlkmod)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make up (z**2)/2 and (z**3)/6 for later use.
c
      do ns = 1,nst
        zx = 0.5*zchar(ns)*zchar(ns)
        zchsq2(ns) = zx
        zchcu6(ns) = (zx*zchar(ns))/3.
      enddo
c
c     Get the max norm of the electrical charges of the aqueous
c     species (izmax).
c
      call zsrt(izmax,narn1,narn2,nstmax,zchar)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Interpret input file data which defines the indices of species
c     required to evaluate such constraints as a mean activity or a
c     heterogeneous equilibrium placed on a basis species.
c
      call intnsp(coval,covali,ier,jflag,narn1,narn2,nbasp,nbt,
     $ nbti,nbtmax,nchlor,ncosp,ndecsp,nhydr,noutpt,nst,nstmax,nttyo,
     $ ucospi,uspec)
c
      if (ier .gt. 0) go to 20
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (neti .gt. 0) then
c
c       Interpret input file data for concentrations and compositions
c       of generic ion exchange phases.
c
        call intge3(cgexp,cgexpi,ier,iern1,iern2,ietmax,jern1,
     $  jern2,jetmax,jgext,jgexti,net,neti,netmax,ngexpi,ngexti,
     $  noutpt,nptmax,nstmax,nttyo,ugexj,ugexji,ugexp,ugexpi,ugexsi,
     $  uphase,uspec,xbar,xbarlg,xgexsi)
c
        if (ier .gt. 0) go to 20
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (net .gt. 0) then
c
c       Set default values for the concentrations (mol/kg.H2O) of
c       generic exchange phases. These values represent trace amounts.
c       Currently, the default concentration of an exchanger phase
c       is 1.e-12 molal.
c
        do ne = 1,net
          if (cgexp(ne) .le. 0.) cgexp(ne) = 1.e-12
        enddo
c
c       Map the concentrations of generic exchange phases into the
c       coval array. Set correpsonding values in the mtb, moph, and
c       loph arrays.
c
        do nb = 1,nbt
          ns = nbasp(nb)
          if (ns.ge.nern1 .and. ns.le.nern2) then
            np = nphasx(ns)
            ne = np - iern1 + 1
            coval(nb) = cgexp(ne)
            mtb(nb) = cgexp(ne)
            moph(np) = cgexp(ne)
            loph(np) = tlg(moph(np))
          endif
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Interpret input file data for compositions of solid solutions.
c     file.
c
      if (iopt(4).ge.1 .and. nxti.gt.0) then
        call intinx(ier,ixrn1,ixrn2,ncmpr,ncmpri,noutpt,npnxp,
     $  nptmax,nstmax,nttyo,nxicmx,nxti,nxtimx,umemi,uphase,usoli,
     $  uspec,xbar,xbari,xbarlg)
c
        if (ier .gt. 0) go to 20
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check to see that no species in the strict basis has a jflag
c     value of 30 at this point.
c
      nerr = 0
      do nb = 1,nbt
        ns = nbasp(nb)
        nt = ndrsr(2,ns) - ndrsr(1,ns) + 1
        if (nt .lt. 2) then
          if (jflag(ns) .eq. 30) then
c
c           Calling sequence substitutions:
c             uspec(ns) for unam48
c
            call fmspnm(jlen,uspec(ns),uspn56)
            write (noutpt,1240) uspn56(1:jlen)
            write (nttyo,1240) uspn56(1:jlen)
 1240       format(/' * Error - (EQ3NR/eq3nr) The strict basis species'
     $      /7x,a,' has a jflag value of 30. This implies',
     $      /7x,'that this species is to be treated as a dependent',
     $      ' species. However,',/7x,"a strict basis species can't",
     $      ' be treated in this manner.')
            nerr = nerr + 1
          endif
        endif
      enddo
c
      if (nerr .gt. 0) go to 20
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Redefine the 'd' set of reactions and attendant data to
c     match the ordinary set as it presently exists. The ordinary
c     set will be futher subjected to eliminations from the active
c     basis set and ordinary basis switching.
c
      call cdrssd(adhfs,adhfsd,advfs,advfsd,axhfs,axhfsd,
     $ axlks,axlksd,axvfs,axvfsd,cdrs,cdrsd,ipch,ipchmx,ipcv,ipcvmx,
     $ narxmx,nbasp,nbaspd,nbtmax,ndrs,ndrsd,ndrsmx,ndrsr,ndrsrd,
     $ nstmax,ntprmx)
      nbtd = nbt
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopg(1) .eq. 1) then
c
c       Build the S-lambda index arrays nsxi and nsxx.
c
        call bdslx(narn1,narn2,natmax,noutpt,nslt,nsltmx,
     $  nslx,nsxi,nsxx,nsxmax,nttyo)
c
c       Build the mu index arrays nmxi and nmxx.
c
        call bdmlx(narn1,narn2,natmax,nmut,nmutmx,nmux,nmxi,
     $  nmxmax,nmxx,noutpt,nttyo)
c
        if (iopr(10) .gt. 0) then
c
c         Write tables concerning Pitzer coefficients.
c
          call ptztab(iopr,narn1,narn2,natmax,nmutmx,nmux,nmxi,
     $    nmxmax,nmxx,noprmx,noutpt,nsltmx,nslx,nstmax,nsxi,nsxmax,
     $    nsxx,uspec)
        endif
c
c       Write warnings for species lacking Pitzer coefficients.
c
        call ptzchk(narn1,narn2,natmax,nmxi,noutpt,nstmax,nsxi,
     $  nttyo,uspec)
c
c       Transform conventional mu data to corresponding C, psi,
c       and zeta data (data originally defined in mu form is
c       not affected).
c
        if (qhawep) then
          call rc3ocf(amu,jpfcmx,ifcphi1,ifcphi2,ifnnn,ifn2n,
     $    ifpsi1,ifpsi2,ifzeta,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,
     $    ilpsi2,ilzeta,iodb,nmux,nmut,nmutmx,nodbmx,noutpt,nstmax,
     $    nttyo,uspec,zchar)
        endif
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the thermodynamic parameters that are functions of
c     temperature and pressure.
c
      call evdata(aadh,aadhh,aadhv,aaphi,abdh,abdhh,abdhv,
     $ abdot,abdoth,abdotv,adadhh,adadhv,adbdhh,adbdhv,adbdth,adbdtv,
     $ adh,adhfe,adhh,adhv,adhfs,adhfsd,advfe,advfs,advfsd,afcnst,
     $ al10,amu,aslm,aphi,aprehw,apresg,apresh,apx,avcnst,axhfe,
     $ axhfs,axhfsd,axlke,axlks,axlksd,axvfe,axvfs,axvfsd,bdh,bdhh,
     $ bdhv,bdot,bdoth,bdotv,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,
     $ dhfe,dhfs,dhfsd,dvfe,dvfs,dvfsd,ehfac,farad,iapxmx,iktmax,
     $ iopg,iopt,ipbtmx,ipch,ipchmx,ipcv,ipcvmx,ixrn1,ixrn2,jpfcmx,
     $ jptffl,jsol,narxmx,narxt,narxth,ncmpr,nmut,nmutmx,nopgmx,
     $ noptmx,noutpt,nptmax,nslt,nsltmx,nst,nstmax,ntpr,ntprmx,nttyo,
     $ nxt,nxtmax,pmu,prehw,presg,presh,press,pslamn,rconst,rcnstv,
     $ rtcnst,tempc,tempk,uphase,uspec,wfac,xhfe,xhfs,xhfsd,xlke,xlks,
     $ xlksd,xvfe,xvfs,xvfsd)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (jpres3 .eq. 0) then
c
c       The pressure corresponds to the data file reference presssure
c       curve.
c
        press = presg
      elseif (jpres3 .eq. 1) then
c
c       The pressure corresponds to the 1.013-bar/steam-saturation
c       presssure curve.
c
        press = presh
      elseif (jpres3 .eq. 2) then
c
c       Constant pressure.
c
        continue
      else
c
c       Error.
c
        write (noutpt,1242) jpres3
        write (nttyo,1242) jpres3
 1242   format(/' * Error - (EQ3NR/eq3nr) The pressure option flag',
     $  ' (jpres3)',/7x,'has an unknown value of ',i2,'.')
      endif
c
c     Stop if the pressure isn't greater than 0.
c
      if (press .le. 0) then
        write (noutpt,1244) press
        write (nttyo,1244) press
 1244   format(/' * Error - (EQ3NR/eq3nr) The pressure must be',
     $  ' greater than zero.',/7x,'The current pressure is ',1pg12.5,
     $  ' bars.')
        stop
      endif
c
c     Stop if press is greater than presmx.
c
      if (press .gt. presmx) then
        write (noutpt,1246) press,presmx
        write (nttyo,1246) press,presmx
 1246   format(/' * Error - (EQ3NR/eq3nr) The calculated pressure is',
     $  ' ',1pg12.5,' bars,',/7x,"greater than the code's built-in",
     $  ' maximum value of ',1pg12.5,' bars.',/7x,'Other limits',
     $  ' associated with the supporting data file may also apply.')
        stop
      endif
c
c     Note if press is less than the 1.013-bar/steam-saturation curve
c     value.
c
      if ((press - presh) .lt. -1.e-4) then
        write (noutpt,1248) press,presh
        write (nttyo,1248) press,presh
 1248   format(/' * Note - (EQ3NR/eq3nr) The calculated pressure is',
     $  ' ',1pg12.5,' bars,',/7x,'less than the 1.013-bar',
     $  ' steam-saturation curve value',/7x,'of ',g12.5,' bars.')
      endif
c
c     If the pressure is nearly identical to the data file reference
c     pressure curve value, set it equal to that value.
c
      if (abs(press - presg) .le. 1.e-4) press = presg
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for thermodynamic pressure corrections.
c
      dp = press - presg
      if (abs(dp) .ge. 1.e-4) then
        if (ipcv .lt. 0) then
c
c         There are no data to support needed pressure corrections.
c
          write (noutpt,1250) press,presg,dp
          write (nttyo,1250) press,presg,dp
 1250     format(/' * Warning - (EQ3NR/eq3nr) The supporting data file',
     $    /7x,'contains no data to support making thermodynamic',
     $    /7x,'pressure corrections. No such corrections will be made.',
     $    /7x,'The current pressure is ',1pg12.5,' bars, the standard',
     $    /7x,'grid pressure is ',g12.5,' bars, and the pressure',
     $    /7x,'difference is ',g12.5,' bars.')
        else
          if (press .le. 0.) then
c
c           The pressure is zero or negative.
c
            write (noutpt,1260) press
            write (nttyo,1260) press
 1260       format(/' * Error - (EQ3NR/eq3nr) The pressure must be',
     $      /7x,'greater than zero. The current pressure is ',1pg12.5,
     $      ' bars.')
            stop
          endif
          if (abs(dp) .gt. prehw) then
c
c           The pressure is outside the recommended envelope.
c
            pxu = presg + prehw
            pxl = presg - prehw
            if (pxl .le. 0.) pxl = 0.
            write (noutpt,1270) press,pxl,pxu,tempc
            write (nttyo,1270) press,pxl,pxu,tempc
 1270       format(/' * Warning - (EQ3NR/eq3nr) The current pressure',
     $      /7x,'of ',1pg12.5,' bars is outside the recommended',
     $      /7x,'pressure envelope of ',g12.5,' to ',g12.5,' bars',
     $      /7x,'at ',0pf6.2,' C.')
          endif
c
c         Make pressure corrections to the thermodynamic data.
c
          call pcorrm(adh,adhh,adhv,al10,aphi,avcnst,bdh,bdhh,bdhv,
     $    bdot,bdoth,bdotv,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,dhfe,
     $    dvfe,iopg,ipch,ipchmx,ipcv,ipcvmx,nopgmx,presg,press,rcnstv,
     $    tempk,xhfe,xlke,xvfe)
c
          call pcorrx(avcnst,dhfs,dvfs,ipch,ipchmx,ipcv,ipcvmx,
     $    nbasp,nbt,nbtmax,ndrsr,nst,nstmax,presg,press,xhfs,
     $    xlks,xvfs)
c
c         Calling sequence substitutions:
c           dhfsd for dhfs
c           dvfsd for dvfs
c           nbaspd for nbasp
c           nbtd for nbt
c           ndrsrd for ndrsr
c           xhfsd for xhfs
c           xlksd for xlks
c           xvfsd for xvfs
c
          call pcorrx(avcnst,dhfsd,dvfsd,ipch,ipchmx,ipcv,ipcvmx,
     $    nbaspd,nbtd,nbtmax,ndrsrd,nst,nstmax,presg,press,xhfsd,
     $    xlksd,xvfsd)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check default redox constraints.
c
      if (irdxc3 .eq. -3) then
        continue
      elseif (irdxc3 .eq. -2) then
        continue
      elseif (irdxc3 .eq. -1) then
        continue
      elseif (irdxc3 .eq. 0) then
        continue
      elseif (irdxc3 .eq. 1) then
        continue
      else
c
c       Error.
c
        write (noutpt,1280) irdxc3
        write (nttyo,1280) irdxc3
 1280   format(/' * Error - (EQ3NR/eq3nr) The default redox option',
     $  /7x,'flag (irdxc3) has an unknown value of ',i2,'.')
        go to 20
      endif
c
c     Make sure that aqueous O2(g) was included in the list of active
c     basis species on the data file.
c
      no2gai = 0
      do nbi = 1,nbti
        nb = ndecsp(nbi)
        if (nb .gt. 0) then
          ns = nbaspd(nb)
          if (ns .eq. no2gaq) then
            no2gai = nbi
            go to 110
          endif
        endif
      enddo
  110 continue
c
      if (irdxc3 .ne.-3 .and. jflag(no2gaq).ne.0) then
        write (noutpt,1300)
 1300   format(/' * Note - (EQ3NR/eq3nr) Have conflicting redox',
     $  ' options:')
        write (noutpt,1310) irdxc3,jflag(no2gaq)
 1310   format(9x,'irdxc3 = ',i2,' overrides jflag(O2(g)) = ',i2,/)
        jflag(no2gaq) = 0
c
c       Calling sequence substitutions:
c         no2gaq for ns
c
        io2gaq = nbasis(nbasp,nbt,nbtmax,no2gaq)
        coval(io2gaq) = -99999.
        ncosp(io2gaq) = 0
        if (no2gai .gt. 0) then
          jflgi(no2gai) = 0
          covali(no2gai) = -99999.
          ucospi(no2gai)(1:48) = ' '
        endif
      endif
c
      j2 = ilnobl(uredox)
      if (irdxc3.ne.1 .and. uredox(1:5).ne.'None ' .and. j2.gt.0) then
        write (noutpt,1300)
        write (noutpt,1320) irdxc3,uredox(1:j2)
 1320   format(9x,'irdxc3 = ',i2,' overrides uredox = ',a,/)
        uredox(1:24) = 'None'
      endif
c
c     If irdxc3 .ge. 1, find the index of the species whose name is
c     uredox.
c
      nredox = 0
      if (irdxc3 .eq. 1) then
        j2 = ilnobl(uredox)
        if (j2.le.0 .or. uspec(ns)(1:4).eq.'None') then
          write (noutpt,1330)
          write (nttyo,1330)
 1330     format(/' * Error - (EQ3NR/eq3nr) The default redox state is',
     $    /7x,'controlled by a couple, but the auxiliary basis species',
     $    /7x,"which defines the couple (uredox) isn't specified on",
     $    /7x,'the input file.')
          go to 20
        else
          do nb = 1,nbt
            ns = nbasp(nb)
            if (uspec(ns)(1:24) .eq. uredox(1:24)) then
              nredox  = ns
              go to 120
            endif
          enddo
c
          write (noutpt,1340) uredox(1:j2)
          write (nttyo,1340) uredox(1:j2)
 1340     format(/' * Error - (EQ3NR/eq3nr) The default redox state is',
     $    /7x,'controlled by a couple, but the species ',a,', whose',
     $    /7x,"associated reaction defines this couple, isn't",
     $    /7x,'in the basis set, as is required.')
          go to 20
        endif
      endif
c
  120 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check the total dissolved solutes option.
c
      if (itdsf3 .eq. 0) then
        continue
      elseif (itdsf3 .eq. 1) then
        continue
      else
c
c       Error.
c
        write (noutpt,1342) itdsf3
        write (nttyo,1342) itdsf3
 1342   format(/' * Error - (EQ3NR/eq3nr) The total dissolved solutes',
     $  /7x,'option flag (itdsf3) has an unknown value of ',i2,'.')
        go to 20
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check the electrical balancing option.
c
      if (iebal3 .eq. 0) then
        continue
      elseif (iebal3 .eq. 1) then
        continue
      else
c
c       Error.
c
        write (noutpt,1344) iebal3
        write (nttyo,1344) iebal3
 1344   format(/' * Error - (EQ3NR/eq3nr) The electrical balancing',
     $  /7x,'option flag (iebal3) has an unknown value of ',i2,'.')
        go to 20
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set default values.
c
      qrho = rho .gt. 0.
      call dfaltx(itermx,rho,scamas,tdspkg,tdspl,tolbt,
     $ toldl,tolspf)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (net .gt. 0) then
c
c       Echo a table for the generic ion exchangers, describing the
c       setup of species, reactions, and corresponding thermodynamic
c       data.
c
        call echgex(axlks,cdrs,cgexj,iern1,iern2,jern1,jern2,
     $  jetmax,jgext,jpflag,jsflag,narxmx,narxt,ndrs,ndrsmx,ndrsr,
     $  netmax,noutpt,nptmax,ntprmx,ntprt,nstmax,press,tempc,ugexj,
     $  ugexmo,uphase,uspec,xlks)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print an echo of the processed input.
c
      call echox(azero,cdrs,covali,eh,fo2lg,iebal3,iodb,iopg,
     $ iopr,iopt,irdxc3,itdsf3,itermx,ixrn1,ixrn2,jflag,jflgi,jpres3,
     $ narn1,narn2,nat,nata,natmax,nbt,nbta,nbti,nbtmax,nct,ncta,
     $ nctmax,ncmpr,ncosp,ndecsp,ndrs,ndrsmx,ndrsr,nelect,ngt,ngta,
     $ ngtmax,nhydr,nhydx,njfmax,nlt,nlta,nltmax,nmt,nmta,nmtmax,
     $ nodbmx,nopgmx,noprmx,noptmx,noutpt,no2gaq,npt,npta,nptmax,
     $ nredox,nst,nsta,nstmax,nxt,nxta,nxti,nxtmax,pe,presg,press,
     $ rho,scamas,tdspkg,tdspl,tempc,tolbt,toldl,tolspf,uactop,
     $ ucospi,uebal,ujflls,uphase,uredox,uspec,uspeci,xbar)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qcwrpj) then
c
c       Write tables of the Pitzer J(x) and J'(x) functions.
c
        call cwrpjt(noutpt)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Convert input concentration data not on the molal scale to
c     that scale. Convert pe, if input, to Eh.
c
      call setup(coval,eh,ehfac,ier,irdxc3,itdsf3,jflag,mwtsp,
     $ narn1,nbaspd,nbt,nbtmax,noutpt,nstmax,nttyo,pe,rho,tdspkg,
     $ tdspl,uspec)
c
      if (ier .gt. 0) go to 20
c
c     Echo input data after modification.
c
      write (noutpt,1400)
 1400 format(/21x,'--- Modified Input Constraints ---',
     $ //5x,'Species',20x,'coval   jflag   Type of Input',/)
      do nb = 1,nbt
        ns1 = nbasp(nb)
        if (jsflag(ns1) .le. 0) then
          jfl = jflag(ns1)
          j2 = ilnobl(ujflls(jfl))
c
          if (jfl.eq.17 .or. jfl.eq.18) then
            nse = ncosp(nb)
            write (noutpt,1410) uspec(ns1),coval(nb),jfl,
     $      ujflls(jfl)(1:j2)
 1410       format(2x,a24,2x,1pe12.5,2x,i2,2x,a)
            j3 = ilnobl(uspec(nse)(1:24))
            write (noutpt,1420) uspec(nse)(1:j3)
 1420       format(39x,'Counterion= ',a)
          elseif (jfl .eq. 25) then
            nse = ncosp(nb)
            write (noutpt,1422) uspec(ns1),jfl,ujflls(jfl)(1:j2)
 1422       format(2x,a24,16x,i2,2x,a)
            j3 = ilnobl(uspec(nse)(1:24))
            ux24 = uspec(nse)(25:48)
            j4 = ilnobl(ux24)
            write (noutpt,1425) uspec(nse)(1:j3),ux24(1:j4)
 1425       format(46x,'Species= ',a,/48x,'Phase= ',a)
c
c           Calling sequence substitutions:
c             noutpt for nf
c             nse for ns
c
            call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,nse,nstmax,uspec)
            write (noutpt,1426)
 1426       format(1x)
c
          elseif (jfl.eq.27 .or. jfl.eq.30) then
            write (noutpt,1422) uspec(ns1),jfl,ujflls(jfl)(1:j2)
          elseif (ns1 .eq. no2gaq) then
            write (noutpt,1427) uspec(ns1),coval(nb),jfl
 1427       format(2x,a24,2x,1pe12.5,2x,i2,2x,'Log fO2')
          else
            if (jfl .ge. 0) then
              write (noutpt,1410) uspec(ns1),coval(nb),jfl,
     $        ujflls(jfl)(1:j2)
            else
              write (noutpt,1430) uspec(ns1)
 1430         format(2x,a24,16x,'Not present in the model')
            endif
          endif
        endif
      enddo
c
      write (noutpt,1440)
 1440 format(1x)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Find the index of the species to be adjusted for electrical
c     balance, if any.
c
      call intieb(iebal,iebal3,ier,jsflag,nbasp,nbt,nbtmax,noutpt,
     $ nstmax,nttyo,uebal,uspec,zchar)
c
      if (ier .gt. 0) go to 20
c
      if (uebal(1:5) .ne. 'None ') then
        j2 = ilnobl(uebal)
        write (noutpt,1450) uebal(1:j2)
 1450   format(/' Electrical balance will be achieved by adjusting',
     $  /'   the concentration of ',a,'.',/)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write a list of inactive species.
c
      write (noutpt,1470)
 1470  format(/21x,'--- Inactive Species ---',/)
c
      k = 0
      do ns = 1,nst
        if (jsflag(ns) .eq. 1) then
          if (ns .ne. no2gaq) then
            k = k + 1
            j2 = ilnobl(uspec(ns)(1:24))
            write (noutpt,1480) uspec(ns)(1:j2)
 1480       format(4x,a)
          endif
        endif
      enddo
c
      if (k .le. 0) write (noutpt,1490)
 1490 format(4x,'None',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Call EQLIB/elim.f to rewrite the reaction equations (cdrs/ndrs/
c     ndrsr arrays) so that auxiliary basis variables with jflag = 30
c     are eliminated from the active basis set.
c
      qelim = .false.
      do nb = 1,nbt
        nse = nbasp(nb)
        nt = ndrsr(2,nse) - ndrsr(1,nse) + 1
        if (nt.ge.2 .and. jflag(nse).eq.30) then
          if (jsflag(nse) .gt. 0) then
            coval(nb) = 0.
          else
            qelim = .true.
            call elim(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,
     $      axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,eps100,ipch,ipchmx,
     $      ipcv,ipcvmx,jsflag,narxmx,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,
     $      nse,nst,nstmax,ntprmx,noutpt,nttyo,uspec)
          endif
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute mass balance coefficients using the active part of the
c     original basis set.
c
      call gcsts(cdrs,csts,jflag,nbaspd,nbt,nbtmax,ndrs,ndrsmx,
     $ ndrsr,noutpt,nsts,nstsmx,nstsr,nst,nstmax,nttyo,uspec)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Copy the existing nbasp array into the nbaspi array. The
c     main purpose of this is to recall what the basis set was
c     when mass balance relationships were defined. In particular,
c     this array holds the basis set as it was prior to the switching
c     out of bare exchangerspecies. This information will be used
c     in describing mass balance totals. In particular, it will be
c     needed to write the pickup file.
c
      do nb = 1,nbt
        nbaspi(nb) = nbasp(nb)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (net .gt. 0) then
c
c       Switch bare site species of generic ion exchange phases with
c       certain models such as Gapon and Vanselow out of the basis set.
c       Then suppress these bare site species so that they appear in the
c       model only in association with corresponding mass balances.
c       Changes are made here to the ion exchanger sections of both the
c       ordinary and the 'd' sets of reactions and associated data.
c
        call chsgex(adhfs,adhfsd,adhfsx,advfs,advfsd,advfsx,
     $  axhfs,axhfsd,axhfsx,axlks,axlksd,axlksx,axvfs,axvfsd,axvfsx,
     $  cdrs,cdrsd,cdrsx,eps100,iern1,ipch,ipchmx,ipcv,ipcvmx,jern1,
     $  jetmax,jflag,jgext,jsflag,narn1,narxmx,narxt,nbasp,nbaspd,
     $  nbaspx,nbt,nbtmax,nbw,ndrs,ndrsd,ndrsmx,ndrsr,ndrsrd,ndrsrx,
     $  ndrsx,nern1,nern2,net,netmax,ngext,noutpt,nphasx,nst,nstmax,
     $  ntprmx,ntprt,nttyo,qbassw,qbswok,ugexmo,uspec)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Save the present configuration of the jflag array. This will
c     be used to support calculation of saturation indices and
c     affinities for the 'd' set of reactions.
c
      call copyia(jflag,jflagd,nbtmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Execute any ordinary basis switching directives from the input
c     file.
c
      kobswt = 0
      if (nobswt .gt. 0) then
c
c       Copy the existing nbasp set into the nbaspx array.
c
        call copyia(nbasp,nbaspx,nbt)
c
c       Interpret the switches.
c
        call intbsw(nbasp,nbaspx,nbt,nbtmax,nobswt,noutpt,nst,
     $  nstmax,nttyo,uobsw,uspec)
c
c       Execute the switches.
c
        do nb = 1,nbt
          ns1 = nbaspx(nb)
          ns2 = nbasp(nb)
          if (ns1 .ne. ns2) then
            call switch(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,
     $      axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,eps100,ipch,ipchmx,
     $      ipcv,ipcvmx,jflag,jsflag,narn1,narxmx,nbasp,nbaspd,nbaspx,
     $      nb,nbt,nbtmax,nbw,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,noutpt,
     $      ns2,nst,nstmax,ntprmx,nttyo,qbassw,qbswok,uspec)
            kobswt = kobswt + 1
          endif
        enddo
      endif
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
c     Call lambda to compute activity coefficients for components of
c     solid solutions for which compositions were specified on the
c     input file.
c
      if (iopt(4).ge.1 .and. nxti.gt.0) then
        do np = ixrn1,ixrn2
c
c         Calling sequence substitutions:
c           acflg for acflgc
c
          call lambda(acflg,afcnst,bpx,ibpxmx,ibpxt,iktmax,ixrn1,
     $    ixrn2,jsol,ncmpr,noutpt,np,nptmax,nstmax,nttyo,nxtmax,wfac,
     $    xbar,xbarlg,uphase,uspec)
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Find which of HCO3- and CO3-- is in the basis set. If both are
c     in the basis set (potentially a bad situation), take the first.
c
      ncarb = 0
      do nb = 1,nbt
        ns = nbasp(nb)
        if (uspec(ns)(1:6).eq.'HCO3- ' .or.
     $    uspec(ns)(1:6).eq.'CO3-- ') then
          icarb = nb
          ncarb = ns
          go to 210
        endif
      enddo
  210 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check the input file for inconsistencies and errors.
c     Here ier accumulates the number of errors caught.
c
      call chkinx(cdrs,coval,ier,irdxc3,jflag,jsflag,narn1,
     $ narn2,nbasp,nbt,nbtmax,ncosp,ndrs,ndrsmx,ndrsr,nelect,nhydr,
     $ nhydx,noutpt,no2gaq,nredox,nstmax,nttyo,tempc,ucospi,uredox,
     $ uspec,zchar)
c
      if (ier .ge. 1) then
        write (noutpt,1510) ier
        write (nttyo,1510) ier
 1510   format(/' * Error - (EQ3NR/eq3nr) The input has failed ',i4,
     $  ' checks.')
        go to 20
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Save the input value of the alkalinity, if any.
c
      alki = 0.
      if (ncarb .gt. 0) then
        jfl = jflag(ncarb)
        if (jfl .eq. 7) alki = coval(icarb)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Recompute equilibrium constants, etc., for the ordinary set
c     if the reactions have been modified by elimination of
c     auxiliary basis species or by basis switching. The presence
c     of generic ion exchangers (net > 0) implies intrinsic or
c     hidden basis switching.
c
      if (qelim .or. kobswt.gt.0 .or. net.gt.0) then
c
        call evdatr(adhfs,advfs,axhfs,axlks,axvfs,dhfs,dvfs,
     $  ipch,ipchmx,ipcv,ipcvmx,narxmx,narxt,nst,nstmax,ntpr,ntprmx,
     $  tempc,xhfs,xlks,xvfs)
c
        if (ipcv .ge. 0) then
          call pcorrx(avcnst,dhfs,dvfs,ipch,ipchmx,ipcv,ipcvmx,
     $    nbasp,nbt,nbtmax,ndrsr,nst,nstmax,presg,press,xhfs,
     $    xlks,xvfs)
        endif
      endif
      kobswt = 0
c
c     Recompute equilibrium constants, etc., for the 'd' set
c     if generic ion exchangers are present. If so, the subset
c     of the 'd' set reactions dealing with such exchangers
c     was modified above by the call to EQLIB/chsgex.f.
c
      if (net .gt. 0) then
c
c       Calling sequence substitutions:
c         adhfsd for adhfs
c         advfsd for advfs
c         axhfsd for axhfs
c         axlksd for axlks
c         axvfsd for avhfs
c         dhfsd for dhfs
c         dvfsd for dvfs
c         xhfsd for xhfs
c         xlksd for xlks
c         xvfsd for xvfs
c
        call evdatr(adhfsd,advfsd,axhfsd,axlksd,axvfsd,dhfsd,dvfsd,
     $  ipch,ipchmx,ipcv,ipcvmx,narxmx,narxt,nst,nstmax,ntpr,ntprmx,
     $  tempc,xhfsd,xlksd,xvfsd)
c
c       Calling sequence substitutions:
c         dhfsd for dhfs
c         dvfsd for dvfs
c         nbaspd for nbasp
c         nbtd for nbt
c         ndrsrd for ndrsr
c         xhfsd for xhfs
c         xlksd for xlks
c         xvfsd for xvfs
c
        if (ipcv .ge. 0) then
          call pcorrx(avcnst,dhfsd,dvfsd,ipch,ipchmx,ipcv,ipcvmx,
     $    nbaspd,nbtd,nbtmax,ndrsrd,nst,nstmax,presg,press,xhfsd,
     $    xlksd,xvfsd)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopr(2) .ge. 1) then
        ilevel = iopr(2)
c
c     Calling sequence substitutions:
c       noutpt for nf
c
        call echolk(axlks,cdrs,ilevel,jsflag,narxmx,ndrs,ndrsmx,
     $  ndrsr,noutpt,nst,ntprmx,nstmax,presg,tempc,uspec,xlks)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize the jcsort, jjsort, jssort, and jgsort arrays by
c     setting each element equal to its index.
c
      call initii(jcsort,nst)
      call initii(jjsort,nst)
      call initii(jssort,nst)
      call initii(jgsort,ngt)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set up the cdrw array. This provides a fast way to get the
c     reaction coefficient of H2O (Aqueous solution) in any reaction.
c
      call gcdrw(cdrs,cdrw,narn1,ndrs,ndrsmx,ndrsr,nst,nstmax)
c
c     The cdrtw array is not actually used by EQ3NR, only EQ6.
c     It is included here for the sake of consistency.
c
      call gcdrtw(cdrs,cdrtw,narn1,narn2,ndrs,ndrsmx,ndrsr,
     $ nelect,no2gaq,nst,nstmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1520)
 1520 format(/' - - BEGIN ITERATIVE CALCULATIONS  - - - - - - - - - -',
     $ ' - - - - - - - - - - - -',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set up the matrix structure (key index array iindx1) for hybrid
c     Newton-Raphson iteration and compute starting values of the matrix
c     variables.
c
      call arrset(aamatr,abar,acflg,acflgo,act,actlg,adh,
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
      if (ker .ge. 2) go to 20
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Find the matrix position of HCO3- or CO3--, if any.
c
      kcarb = 0
      if (ncarb .gt. 0) then
        do kcol = 1,kbt
          nb = iindx1(kcol)
          ns = nbasp(nb)
          if (ns .eq. ncarb) then
            kcarb = kcol
            go to 220
          endif
        enddo
  220   continue
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the speciation-solubility model.
c
      write (noutpt,1600)
      write (nttyo,1600)
 1600 format(/,' Starting hybrid Newton-Raphson iteration.',/)
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
        write (noutpt,1610)
        write (nttyo,1610)
 1610   format(/' * Error - (EQ3NR/eq3nr) Hybrid Newton-Raphson',
     $  ' iteration failed')
        ux8a = ' '
        write (ux8a,1620) iter
 1620   format(i7)
        call lejust(ux8a)
        j2 = ilnobl(ux8a)
        if (ier. eq. 1) then
          write (noutpt,1630) ux8a(1:j2)
          write (nttyo,1630) ux8a(1:j2)
 1630     format(7x,'after ',a,' iterations because a zero matrix',
     $    ' was encountered. This is',/7x,'probably due to a',
     $    ' programming error.')
        elseif (ier .eq. 2) then
          write (noutpt,1640) ux8a(1:j2)
          write (nttyo,1640) ux8a(1:j2)
 1640     format(7x,'after ',a,' iterations because a non-zero,',
     $    ' computationally singular',/7x,'matrix was encountered.')
        elseif (ier .eq. 3) then
          write (noutpt,1650) ux8a(1:j2)
          write (nttyo,1650) ux8a(1:j2)
 1650     format(7x,'after ',a,' iterations because the code',
     $    ' detected that',/7x,'iteration was diverging.')
        elseif (ier .eq. 4) then
          write (noutpt,1660) ux8a(1:j2)
          write (nttyo,1660) ux8a(1:j2)
 1660     format(7x,'after ',a,' iterations because the maximum',
     $    ' number of iterations',/7x,'was done.')
        else
          ux8b = ' '
          write (ux8b,1620) ier
          call lejust(ux8b)
          j3 = ilnobl(ux8b)
          write (noutpt,1670) ux8a(1:j2),ux8b(1:j3)
          write (nttyo,1670) ux8a(1:j2),ux8b(1:j3)
 1670     format(7x,'after ',a,' iterations because an unknown event',
     $    ' occurred. The ier',/7x,'error code has the unknown value',
     $    ' ',a,'. This condition is a',/7x,'programming error.')
        endif
c
        call ndiagx(delmax,delvec,eps100,idelmx,iebal,iindx1,
     $  irdxc3,jflag,kcarb,kebal,khydr,kmax,ko2gaq,nbasp,nbtmax,nhydr,
     $  noutpt,nstmax,nttyo,screwd,uspec)
        go to 20
      endif
c
      write (noutpt,1690) iter
      write (nttyo,1690) iter
 1690 format('   Done. Hybrid Newton-Raphson iteration converged in ',
     $ i3,' iterations.',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Save the jflag value of the species being adjusted for electrical
c     balance.
c
      jfleba = 0
      if (iebal .gt. 0) then
        ns = nbaspd(iebal)
        jfleba = jflag(ns)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Reset jflag values which are not 0 or 30 to one of these values.
c     Species for which jflag will be set to 0 now have associated
c     mass balance totals, those for which jflag will be set to 30
c     do not.
c
      if (nredox .ge. 1) jflag(nredox) = 0
c
      do nb = 1,nbt
        ns = nbasp(nb)
        jflgi(nb) = jflag(ns)
        if (jflgi(nb) .eq. 27) jflgi(nb) = 30
        if (jflag(ns) .ne. 30) jflag(ns) = 0
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the weights (actually, the masses) of total dissolved
c     solutes, the solvent, and the solution, and also some related
c     quantities. Also get the solution density (rhowc) from the
c     WIPP brine density model if the temperature lies in the range
c     20-30 C. This model is a fit to NaCl solution data having
c     the form:
c
c       density (g/L) = a + b x TDS (g/L)
c
c     where a = 1000.96 and b = 0.639963. Also get the molarity/
c     molality conversion ratio (mrmlra, molarity = mrmlra x molality)
c     and its inverse.
c
      call gwdenp(adwipp,bdwipp,jcsort,mlmrra,mosp,mrmlra,
     $ mwtsp,narn1,narn2,nstmax,qdwipp,rhoc,rhowc,tdsgks,tdsglw,
     $ tdspkc,tdsplc,tempc,vosol,wfh2o,wftds,wkgwi,woh2o,
     $ wosol,wotds)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check potential difference between input and calculated
c     densities.
c
      if (rho.gt.0. .and. rhoc.gt.0.) then
        axx = abs((rhoc - rho)/rho)
        if (axx .gt. 0.01) then
          write (ux16a,'(1pg12.5)') rhoc
          call lejust(ux16a)
          j1 = ilnobl(ux16a)
          write (ux16b,'(1pg12.5)') rho
          call lejust(ux16b)
          j2 = ilnobl(ux16b)
          write (noutpt,1720) ux16a(1:j1),ux16b(1:j2)
          write (nttyo,1720) ux16a(1:j1),ux16b(1:j2)
 1720     format(/' * Warning - (EQ3NR/eq3nr) The calculated density',
     $    ' of ',a,' g/mL',/7x'differs from the input file/default',
     $    ' value of ',a,' g/mL',/7x,'by more than 1%.',
     $    ' The calculated value will be used',/7x,'in subsequent',
     $    ' calculations.')
        endif
      endif
c
c     Check potential difference between input and calculated
c     TDS values.
c
      if (itdsf3 .le. 0) then
        if (tdspkg.gt.0. .and. tdspkc.gt.0.) then
          axx = abs((tdspkg - tdspkc)/tdspkc)
          if (axx .gt. 0.01) then
            write (ux16a,'(1pg12.5)') tdspkc
            call lejust(ux16a)
            j1 = ilnobl(ux16a)
            write (ux16b,'(1pg12.5)') tdspkg
            call lejust(ux16b)
            j2 = ilnobl(ux16b)
            write (noutpt,1730) ux16a(1:j1),ux16b(1:j2)
            write (nttyo,1730) ux16a(1:j1),ux16b(1:j2)
 1730       format(/' * Warning - (EQ3NR/eq3nr) The calculated TDS',
     $      ' of ',a,' mg/kg.sol',/7x'differs from the input',
     $      ' file/default value of ',a,' mg/kg.sol',/7x,'by more',
     $      ' than 1%. The calculated value will be used in',
     $      /7x,'subsequent calculations.')
          endif
        elseif (tdspkg.le.0. .and. tdspkc.gt.0.) then
          write (ux16a,'(1pg12.5)') tdspkc
          call lejust(ux16a)
          j1 = ilnobl(ux16a)
          write (ux16b,'(1pg12.5)') tdspkg
          call lejust(ux16b)
          j2 = ilnobl(ux16b)
          write (noutpt,1732) ux16a(1:j1),ux16b(1:j2)
          write (nttyo,1732) ux16a(1:j1),ux16b(1:j2)
 1732     format(/' * Warning - (EQ3NR/eq3nr) The calculated TDS',
     $    ' of ',a,' mg/kg.sol',/7x'differs from the input',
     $    ' file/default value of ',a,' mg/kg.sol.',/7x,
     $    'The calculated value will be used in subsequent',
     $    ' calculations.')
        endif
      else
        if (tdspl.gt.0. .and. tdsplc.gt.0.) then
          axx = abs((tdspl - tdsplc)/tdsplc)
          if (axx .gt. 0.01) then
            write (ux16a,'(1pg12.5)') tdsplc
            call lejust(ux16a)
            j1 = ilnobl(ux16a)
            write (ux16b,'(1pg12.5)') tdspl
            call lejust(ux16b)
            j2 = ilnobl(ux16b)
            write (noutpt,1740) ux16a(1:j1),ux16b(1:j2)
            write (nttyo,1740) ux16a(1:j1),ux16b(1:j2)
 1740       format(/' * Warning - (EQ3NR/eq3nr) The calculated TDS',
     $      ' of ',a,' mg/L',/7x'differs from the input',
     $      ' file value of ',a,' mg/L',/7x,'by more',
     $      ' than 1%. The calculated value will be used in',
     $      /7x,'subsequent calculations.')
          endif
        elseif (tdspl.le.0. .and. tdsplc.gt.0.) then
          write (ux16a,'(1pg12.5)') tdsplc
          call lejust(ux16a)
          j1 = ilnobl(ux16a)
          write (ux16b,'(1pg12.5)') tdspl
          call lejust(ux16b)
          j2 = ilnobl(ux16b)
          write (noutpt,1742) ux16a(1:j1),ux16b(1:j2)
          write (nttyo,1742) ux16a(1:j1),ux16b(1:j2)
 1742     format(/' * Warning - (EQ3NR/eq3nr) The calculated TDS',
     $    ' of ',a,' mg/L',/7x'differs from the input',
     $    ' file/default value of ',a,' mg/L.',/7x,
     $    'The calculated value will be used in subsequent',
     $    ' calculations.')
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     If there is a density calculated from the WIPP brine density
c     model, use it as the density in subsequent calculations. If not
c     and a density value greater than zero was entered on the
c     input file, use the input value.
c
      if (qdwipp) then
        qrho = .true.
        rho = rhoc
        tdspkg = tdspkc
        tdspl = tdsplc
      elseif (qrho) then
        mrmlra = 0.001*wfh2o*rho
        mlmrra = 1./mrmlra
      endif
      if (.not.qrho) rho = 0.0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute mass balance totals for chemical elements.
c
      do nc = 1,nctmax
        mte(nc) = 0.
        mteaq(nc) = 0.
        cteaq(nc) = 0.
      enddo
c
      do nss = narn1,narn2
        ns = jcsort(nss)
        if (jsflag(ns) .le. 0) then
          nr1 = nessr(1,ns)
          nr2 = nessr(2,ns)
          do n = nr1,nr2
            nc = ness(n)
            mteaq(nc) = mteaq(nc) + cess(n)*mosp(ns)
          enddo
        endif
      enddo
c
      do nc = 1,nct
        mte(nc) = mteaq(nc)
        cteaq(nc) = wkgwi*mteaq(nc)
        ppmwe(nc) = 1000.*cteaq(nc)*atwt(nc)*wfh2o
      enddo
c
      do nss = nern1,nern2
        ns = jcsort(nss)
        if (jsflag(ns) .le. 0) then
          nr1 = nessr(1,ns)
          nr2 = nessr(2,ns)
          do n = nr1,nr2
            nc = ness(n)
            if (nc .gt. 0) mte(nc) = mte(nc) + cess(n)*mosp(ns)
          enddo
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute mass balance totals for active basis species.
c
      do nb = 1,nbtmax
        mtb(nb) = 0.
        mtbaq(nb) = 0.
        ctb(nb) = 0.
      enddo
c
      do nss = narn1,narn2
        ns = jcsort(nss)
        if (jsflag(ns) .le. 0) then
          nr1 = nstsr(1,ns)
          nr2 = nstsr(2,ns)
          do n = nr1,nr2
            nb = nsts(n)
            mtbaq(nb) = mtbaq(nb) + csts(n)*mosp(ns)
          enddo
        endif
      enddo
c
      do nb = 1,nbt
        mtb(nb) = mtbaq(nb)
        ctb(nb) = wkgwi*mtbaq(nb)
      enddo
c
      do nss = nern1,nern2
        ns = jcsort(nss)
        if (jsflag(ns) .le. 0) then
          nr1 = nstsr(1,ns)
          nr2 = nstsr(2,ns)
          do n = nr1,nr2
            nb = nsts(n)
            mtb(nb) = mtb(nb) + csts(n)*mosp(ns)
          enddo
        endif
      enddo
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
c     Initialize the flag array which marks if compositions are known
c     which maximize the affinities of the phases. Here assume that
c     this is the case for all phases except solid solutions. The
c     compositions for these phases will be determined by EQLIB/hpsat.f,
c     which will be called below by scripx.f. Note that input solid
c     solution compositions do not in general maximize the affinities.
c     Thus, they don't count here.
c
      do np = 1, npt
        qxknph(np) = .true.
      enddo
      do np = ixrn1,ixrn2
        qxknph(np) = .false.
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print the computed description of the aqueous solution.
c
      call scripx(abar,acflg,act,actlg,adh,afcnst,affpd,
     $ affsd,ahrc,alki,apx,atwt,a3bar,a3bars,bpx,cdrsd,cegexs,
     $ cess,conc,conclg,coval,csts,ctb,cteaq,egexjc,egexjf,egexpa,
     $ egexpc,egexs,egexw,ehfac,ehrc,eps100,eh,farad,fje,fo2,fo2lg,
     $ fo2lrc,fugac,fugalg,fxi,iapxmx,ibpxmx,iebal,iern1,iern2,ietmax,
     $ igas,iktmax,iopg,iopr,iopt,ilrn1,ilrn2,imrn1,imrn2,ixrn1,
     $ ixrn2,jcsort,jern1,jern2,jetmax,jflag,jflagd,jflgi,jfleba,jgext,
     $ jgsort,jpflag,jsflag,jsol,jsomax,kern1,kern2,ketmax,kgexsa,
     $ mlmrra,mrmlra,moph,mosp,mte,mteaq,mwtsp,narn1,narn2,natmax,
     $ nbasp,nbaspd,nbt,nbtmax,nchlor,ncmpr,nct,nctmax,ndrsd,ndrsmx,
     $ ndrsrd,nelect,ness,nessmx,nessr,net,neti,netmax,ngexpi,ngexsa,
     $ ngext,ngrn1,ngrn2,ngt,ngtmax,nhydr,nhydx,nopgmx,noprmx,noptmx,
     $ noutpt,no2gaq,npnxp,npt,nptmax,nrdxsp,nst,nstmax,nsts,nstsmx,
     $ nstsr,ntf1,ntf1mx,ntf1t,ntf2,ntf2mx,ntf2t,nttyo,nxrn1,nxrn2,
     $ nxt,nxti,nxtimx,nxtmax,omega,pe,perc,ppmwe,qrho,qxknph,rho,
     $ rhoc,rhowc,sidrph,sidrsp,sigmam,sigzi,tdsglw,tdspkc,tdspkg,
     $ tdspl,tdsplc,tempc,tf1,tf2,tolspf,uelem,ugexj,ugexmo,uphase,
     $ uspec,uxtype,vosol,wfac,wfh2o,wftds,wkgwi,woh2o,wosol,wotds,
     $ xbar,xbarlg,xbarw,xbrwlg,xgexw,xlke,xlksd,zchar,zchcu6,zchsq2)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nprob .le. 1) then
c
c       Calculate data to describe the aqueous solution in the first
c       problem on the input file as a special reactant. This defines
c       the so-called 'Fluid 2' to be described on the EQ3NR pickup
c       file under the iopt(19) = 3 option.
c
        ureac1 = 'Fluid 2'
        tempc1 = tempc
c
        do nc = 1,nct
          n = nc
          uesr1(n) = uelem(nc)
          cesr1(n) = mteaq(nc)
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
            cbsr1(n) = mtbaq(nb)
          endif
        enddo
        ibsrt1 = n
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(17) .ge. 0) then
c
c       Write the pickup file.
c
c       First set up some needed data.
c
        call setpk3(electr,iindx1,jflag,jflgi,kbt,kdim,kmax,
     $  kmt,kprs,kwater,kxt,mtb,mtbi,mtbaq,mtbaqi,narn1,narn2,nbasp,
     $  nbaspd,nbaspi,nbti,nbtmax,ndrsrd,nern1,nern2,nobswt,nstmax,
     $  ntitl,ntitl2,ntitmx,omeglg,press,pressi,scamas,sigzi,tempc,
     $  tempci,ubmtbi,uobsw,uspec,utitl,utitl2,uzveci,uzvec1,
     $  zvclgi,zvclg1)
c
        if (iopr(17) .le. 0) then
          upkfor = uinfor
        elseif (iopr(17) .eq. 1) then
          upkfor = 'W'
        elseif (iopr(17) .ge. 2) then
          upkfor = 'D'
        endif
c
        if (iopt(19) .le. 0) then
c
c         Write a normal EQ3NR pickup file. This corresponds to the
c         bottom half of a full EQ6 input file. To use this pickup
c         file in an EQ6 input file, you must concatenate it with the
c         top half. Usually this top half is extracted from an
c         existing EQ6 input file and edited as necessary.
c
          if (upkfor(1:1) .eq. 'W') then
c
c           Compact (W) format.
c
            call wr3pkw(electr,cgexj,ietmax,iopg,jetmax,jflgi,
     $      jgext,kbt,kct,kdim,kmax,kmt,kxmod,kxt,mtbi,mtbaqi,mwtges,
     $      nbti,nbtmax,net,netmax,newin,ngexrt,nobswt,nopgmx,nsbswt,
     $      ntitl2,ntitmx,nxmdmx,nxmod,qgexsh,pressi,tempci,tgexp,
     $      ubmtbi,ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,usbsw,utitl2,
     $      uvfgex,uxkgex,uxmod,uzveci,xhfgex,xlkgex,xlkmod,xvfgex,
     $      zgexj,zvclgi)
          else
c
c           Menu-style (D) format.
c
            call wr3pkd(electr,cgexj,ietmax,iopg,jetmax,jflgi,
     $      jgext,kbt,kct,kdim,kmax,kmt,kxmod,kxt,mtbi,mtbaqi,mwtges,
     $      nbti,nbtmax,net,netmax,newin,ngexrt,nobswt,nopgmx,nsbswt,
     $      ntitl2,ntitmx,nxmdmx,nxmod,qgexsh,pressi,tempci,tgexp,
     $      ubmtbi,ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,usbsw,utitl2,
     $      uvfgex,uxkgex,uxmod,uzveci,xhfgex,xlkgex,xlkmod,xvfgex,
     $      zgexj,zvclgi)
          endif
        else
c
c         Write an advanced EQ3NR pickup file. This corresponds to
c         a full EQ6 input file. There are currently three options
c         for the form and content of the top half of this file.
c         These are determined by the iopt(19) option switch. See
c         comments in EQ3NR/stpk36.f.
c
          itmxsv = itermx
c
          call stpk36(awmaxi,awmini,cbsri,cbsr1,cdac,cesri,cesr1,
     $    csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,
     $    dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,
     $    dlxplo,dlxprl,dlxprn,eact,ehmaxi,ehmini,fkrc,hact,iact,
     $    ibsrti,ibsrt1,iesrti,iesrt1,iktmax,imchmx,imech,iopt,itermx,
     $    ixrti,jcode,jpress,jreac,jtemp,ksplmx,ksppmx,kstpmx,modr,
     $    moffg,morr,mprphi,mprspi,nbt1mx,nctmax,ndact,ndctmx,nffg,
     $    nffgmx,noptmx,nordmx,noutpt,nprob,nprpmx,nprpti,nprsmx,
     $    nprsti,nptkmx,nrct,nrctmx,nrk,nsk,nsrt,nsrtmx,ntitl1,ntitmx,
     $    ntrymx,nttkmx,nttyo,nxopex,nxopmx,nxopt,nxpemx,nxrt,nxrtmx,
     $    o2maxi,o2mini,phmaxi,phmini,ptk,press,pressb,rkb,rxbari,
     $    sfcar,ssfcar,tempc,tempcb,tempc1,timmxi,tistti,tolsat,
     $    tolxsf,trkb,ttk,ubsri,ubsr1,ucxri,udac,uesri,uesr1,uffg,
     $    uprphi,uprspi,ureac,ureac1,utitl1,uxcat,uxopex,uxopt,vreac,
     $    ximaxi,xistti,xlkffg)
c
          itermx = itmxsv
c
          if (upkfor(1:1) .eq. 'W') then
c
c           Compact (W) format.
c
            call wr6pkw(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,
     $      dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,
     $      dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,
     $      dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,iact,
     $      ibsrti,iesrti,ietmax,igerti,iktmax,imchmx,imech,iodb,iopg,
     $      iopr,iopt,itermx,ixrti,jcode,jetmax,jflgi,jgerti,jgext,
     $      jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,
     $      kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,
     $      mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,
     $      nertmx,net,netmax,newin,nffg,nffgmx,ngexrt,nobswt,nodbmx,
     $      nopgmx,noprmx,noptmx,nordmx,nprpmx,nprpti,nprsmx,nprsti,
     $      nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,
     $      ntitmx,ntrymx,nttkmx,nxmdmx,nxmod,nxopex,nxopmx,nxopt,
     $      nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,pressb,
     $      pressi,ptk,qgexsh,rkb,rxbari,sfcar,ssfcar,tempcb,tempci,
     $      tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,trkb,ttk,
     $      ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,ugersi,
     $      ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,ureac,
     $      usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,uxopt,
     $      uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,
     $      xlkmod,xvfgex,zgexj,zvclgi)
          else
c
c           Menu-style (D) format.
c
            call wr6pkd(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,
     $      dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,
     $      dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,dlxplo,
     $      dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,iact,
     $      ibsrti,iesrti,ietmax,igerti,iktmax,imchmx,imech,iodb,iopg,
     $      iopr,iopt,itermx,ixrti,jcode,jetmax,jflgi,jgerti,jgext,
     $      jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,
     $      kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,
     $      mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,
     $      nertmx,net,netmax,newin,nffg,nffgmx,ngexrt,nobswt,nodbmx,
     $      nopgmx,noprmx,noptmx,nordmx,nprpmx,nprpti,nprsmx,nprsti,
     $      nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,
     $      ntitmx,ntrymx,nttkmx,nxmdmx,nxmod,nxopex,nxopmx,nxopt,
     $      nxpemx,nxrt,nxrtmx,o2maxi,o2mini,phmaxi,phmini,pressb,
     $      pressi,ptk,qgexsh,rkb,rxbari,sfcar,ssfcar,tempcb,tempci,
     $      tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,trkb,ttk,
     $      ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,ugersi,
     $      ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,ureac,
     $      usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,uxopt,
     $      uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,
     $      xlkmod,xvfgex,zgexj,zvclgi)
          endif
c
        endif
c
        write (noutpt,1800)
        write (nttyo,1800)
 1800   format(/' The pickup file has been written.')
      else
c
        write (noutpt,1810)
        write (nttyo,1810)
 1810   format(/' No pickup file was written.')
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Go look for another problem on the input file.
c
      go to 20
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
