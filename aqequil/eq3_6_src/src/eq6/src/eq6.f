      program eq6
c
c     This is the main program of the EQ6 code. Configuration
c     identification, the copyright statement, legal disclaimers,
c     and similar statements are contained in EQ6/aaaeq6.f, the
c     lead-off subroutine in the EQ3NR source code. A short description
c     of this program is also contained in that subroutine.
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
c
      include 'eqlib/eqldv.h'
      include 'eqlib/eqlge.h'
      include 'eqlib/eql1s.h'
c
      include 'eqlib/eqlo8.h'

c-----------------------------------------------------------------------
c     File path parameters
c-----------------------------------------------------------------------
      integer :: numargs
      character(1024) :: temppath
      character(:), allocatable :: data1path
      character(:), allocatable :: sixipath
      integer :: pathindices(2)
      character(:), allocatable :: basename
      character(:), allocatable :: ofile
      character(:), allocatable :: pfile
      character(:), allocatable :: bafile
      character(:), allocatable :: bbfile
      character(:), allocatable :: ifile
      character(:), allocatable :: tfile
      character(:), allocatable :: txfile
      character(:), allocatable :: tsfile
c
c-----------------------------------------------------------------------
c
      integer, parameter :: nllnpa = 7200
c
c-----------------------------------------------------------------------
c
c     Array allocation size variables used in EQ6.
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
      integer nbt1_asv,nmx_asv,nsts_asv,nsx_asv
c
      integer ntf1_asv,ntf2_asv
c
      integer k_asv,npet_asv,nset_asv
c
      integer imch_asv,ndct_asv,nert_asv,nffg_asv,nprp_asv,nprs_asv,
     $ nptk_asv,nrct_asv,nsrt_asv,nttk_asv,nxop_asv,nxpe_asv,nxrt_asv
c
c-----------------------------------------------------------------------
c
      integer nad1,nbkupa,nbkupb,newin,ninpt,ninpts,noutpt,ntab,ntabs,
     $ ntabx,nttyo
c
      integer iodb(nodb_par),iopg(nopg_par),iopr(nopr_par),
     $ iopt(nopt_par)
c
      integer igerti(jet_par,nert_par),jgerti(nert_par)
c
      integer, dimension(:), allocatable :: iapxt,ibpxt,ibsrti,iesrti,
     $ iffg,iindx1,insgf,ipndx1,ixrti,jcode,jffg,jflag,jflagd,jflgi,
     $ jpflag,jreac,jsflag,jsitex,jsol,kxmod
c
      integer, dimension(:), allocatable :: narxt,nbasp,nbaspd,nbaspi,
     $ nbaspx,nbmap,ncmap,ndecsp,ndrs,ndrsd,ndrsx,ness,npchk,nphasx,
     $ nrndex,nsk,nsmap,nsts,ntf1,ntf2,nxridx
c
      integer, dimension(:,:), allocatable :: imech,ncmpr,ndrsr,ndrsrd,
     $ ndrsrx,nessr,nrk,nstsr
c
      integer, dimension(:,:,:), allocatable :: iact,ndact
c
      integer, dimension(:,:,:,:), allocatable :: ndac
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
      integer i,ie,ier,iaqsla,iaqsln,ielam,iexec0,igas,ilevel,ipch,ipcv,
     $ irang,itermx,ixrn1a,ixrn2a,izmax,j,je,jexec0,jlen,jpdblo,jpress,
     $ jptffl,jtemp,j2,j3,k,kbt,kcol,kct,kdim,kelect,khydr,khydx,klim,
     $ km1,kmt,kobswt,ko2gaq,krdxsp,ksplmx,ksppmx,kstpmx,kwater,kx1,kxt
c
      integer n,napta,narn1a,narn2a,nart,nata,nb,nbta,nbtafd,nbtd,nbw,
     $ nbwa,nb1,nb2,nchloa,nchlor,ncta,ne,neleca,nelect,ner,nerr,nert,
     $ nffg,ngrn1a,ngrn2a,ngrt,ngta,nhydr,nhydra,nhydx,nhydxa,nllnmx,
     $ nlrn1a,nlrn2a,nlta,nmax,nmrn1a,nmrn2a,nmrt,nmta,nmuta,nobswt,
     $ no2gaa,no2gaq,np,npi,npslmx,npta,nprob,nrc,nrct,nrdxsa,nrdxsp,
     $ nrecl,nrr1,nrr2,nr1,nr2,nsbsw,nsbswt,ns,nse,nsi,nslta,nsrt,nss,
     $ nsslmx,nsta,ns1,ns2,nt,ntf1t,ntf1ta,ntf2t,ntf2ta,ntitld,ntitl1,
     $ ntitl2,ntpr,ntprh,ntprt,ntrymx,nxmod,nxopex,nxopt,nxrn1a,nxrn2a,
     $ nxrt,nxta
c
      integer ilnobl,nbasis
c
      logical qbassw,qbswok,qchlor,qcnpre,qcntmp,qcwrpj,qdwipp,qecon,
     $ qelim,qend,qex,qgexsh,qhawep,qloffg,qop,qoptmz,qpit75,qrderr,
     $ qredox,qscon,qtatxt
c
      character*24 ugersi(iet_par,jet_par,nert_par),ugermo(nert_par)
      character*8 ugerji(jet_par,nert_par)
c
      character(len=80), dimension(:), allocatable :: utitl1,utitl2
      character(len=48), dimension(:), allocatable :: ubasp,ubmtbi,
     $ uprspi,uspec,uxmod,uzveci,uzvec1
      character(len=24), dimension(:), allocatable :: uffg,uphase,
     $ uprphi,uptype,ureac,uxcat,uxopex
      character(len=8), dimension(:), allocatable :: uelem,uxopt
c
      character(len=48), dimension(:,:), allocatable :: uobsw,usbsw
      character(len=24), dimension(:,:), allocatable :: ubsri,ucxri
      character(len=8), dimension(:,:), allocatable :: uesri
c
      character(len=24), dimension(:,:,:,:), allocatable :: udac
c
      character(len=nllnpa) ulinex
      character(len=80) ux80
      character(len=56) uspn56
      character(len=32) uactop
      character(len=24) unamsp,unamph,ux
      character(len=24) uaqsln,ublk24
      character(len=11) utime0,utime1
      character(len=9) udate0,udate1
      character(len=8) udatfi,ufixf,udakey,uinfor,uplatc,uplatm,
     $ ustelg,ustelu,usteql,usteq6,uveelg,uveelu,uveeql,uveeq6,ux8
c
      real*8 sscrew(nssc_par)
      real*8 cco2(5)
c
      real*8 egers(iet_par,jet_par,nert_par),
     $ egersi(iet_par,jet_par,nert_par),
     $ mrgers(iet_par,jet_par,nert_par),
     $ xgers(iet_par,jet_par,nert_par),
     $ xgersi(iet_par,jet_par,nert_par)
c
      real*8 cegexs(iet_par,jet_par,net_par),
     $ cpgexs(iet_par,jet_par,net_par),egexjf(jet_par,net_par),
     $ mrgexs(iet_par,jet_par,net_par)
c
      real(8), dimension(:), allocatable :: atwt,azero,cdrs,
     $ cdrsd,cdrsx,cess,cscale,csts,elecsr,fkrc,loph,losp,
     $ modr,moffg,moph,morr,mosp,mprph,mprphi,mprsp,mprspi,mtb,
     $ mtbaq,mtbaqi,mtbi,mte,mteaq,mwtrc,mwtsp,ptk,sfcar,ssfcar,
     $ tempcu,tf1,tf2,ttk,xlkffg,xlkmod,vosp0,vreac,zchar,zchcu6,
     $ zchsq2,zvclgi,zvclg1,zvec1
c
      real(8), dimension(:,:), allocatable :: cbsr,cesr,rxbar
      real(8), dimension(:,:), allocatable :: cbsri,cesri,rxbari
      real(8), dimension(:,:,:), allocatable :: csigma,eact,hact,
     $ rk,rkb,trkb
      real(8), dimension(:,:,:,:), allocatable :: cdac
c
      real(8), dimension(:), allocatable :: dadhh,dadhv,dbdhh,dbdhv,
     $ dbdth,dbdtv
c
      real(8), dimension(:,:), allocatable :: apx,bpx,wfac
c
      real(8), dimension(:,:), allocatable :: aadh,aadhh,aadhv,aaphi,
     $ abdh,abdhh,abdhv,abdot,abdoth,abdotv
c
      real(8), dimension(:,:,:), allocatable :: adadhh,adadhv,adbdhh,
     $ adbdhv,adbdth,adbdtv
c
      real*8 adh,adhh,adhv,aphi,bdh,bdhh,bdhv,bdot,bdoth,bdotv
c
      real*8 afcnst,aftarg,al10,avcnst,awmax,awmaxi,awmin,awmini,
     $ dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,dloprn,dltpll,
     $ dltplo,dltprl,dltprn,dlxmax,dlxdmp,dlxmin,dlxmx0,dlxpll,dlxplo,
     $ dlxprl,dlxprn,ehfac,ehmax,ehmaxi,ehmin,ehmini,electr,eps100,
     $ farad,o2max,o2maxi,o2min,o2mini,phmax,phmaxi,phmin,phmini,
     $ prcinf,press,pressb,pressd,pressi,rconst,rcnstv,rtcnst,
     $ smp100,tcpu,tdamax,tdamin,tempc,tempcb,tempcd,tempci,tempk,
     $ texec0,timemx,timmxi,time1,tistrt,tistti,tolaft,tolbt,toldl,
     $ tolsat,tolsst,tolxsf,tolxst,tolxsu,trun,tuser,ximax,ximaxi,
     $ xistrt,xistti,xi1,x10,zkfac,zklgmn,zklogl,zklogu
c
      real*8 av,azch,azchmx,dp,dt,ex,gx,xx,zx
c
cXX   real*8 vpgstp
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
        external bkdeq6
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
      data ufixf  /'fix_f   '/
c
c     Set practical infinity.
c
      data prcinf/1.e+38/
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
      nbkupa = 0
      nbkupb = 0
      ninpt = 0
      ninpts = 0
      newin = 0
      ntab = 0
      ntabs = 0
      ntabx = 0
      nrecl = 0
c
c     Open all files execpt pickup, tab, tabs, and tabx.
c
      numargs = COMMAND_ARGUMENT_COUNT()
      if (numargs.eq.2) then
          call GET_COMMAND_ARGUMENT(1,temppath)
          data1path = TRIM(temppath)
          temppath(:)=''
          call GET_COMMAND_ARGUMENT(2,temppath)
          sixipath = TRIM(temppath)
      else
          write (0, *) 'usage: eq6 <data1> <6i>'
          stop
      end if

      call getbasename(sixipath, pathindices)
      basename = sixipath(pathindices(1):pathindices(2))

      ofile = basename // '.6o'
      pfile = basename // '.6p'
      bafile = basename // '.6ba'
      bbfile = basename // '.6bb'
      ifile = basename // '.6ib'
      tfile = basename // '.6t'
      txfile = basename // '.6tx'
      tsfile = basename // '.6ts'

      call openin(noutpt,nttyo,data1path,'unformatted',nad1)
      call openin(noutpt,nttyo,sixipath,'formatted',ninpt)

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
      call aaaeq6(usteq6,uveeq6)
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
      i = index(uveeq6,' ') - 1
      j = index(uveeq6,'.') - 1
      k = index(uplatc,' ') - 1
      if (i .le. 0) i = 8
      if (j .le. 0) j = 8
      if (k .le. 0) k = 8
      write (nttyo,1000) uveeq6(1:i),uveeq6(1:j),uveeq6(1:i),
     $ uplatc(1:k)
      write (noutpt,1000) uveeq6(1:i),uveeq6(1:j),uveeq6(1:i),
     $ uplatc(1:k)
 1000 format(/' EQ3/6, Version ',a,' (EQ3/6-V',a,'-REL-V',a,'-',a,')')
c
      i = j
      j = index(usteq6,' ') - 1
      k = index(uplatm,' ') - 1
      if (j .le. 0) j = 8
      if (k .le. 0) k = 8
      write (nttyo,1010) uveeq6(1:i),usteq6(1:j),uplatm(1:k)
      write (noutpt,1010) uveeq6(1:i),usteq6(1:j),uplatm(1:k)
 1010 format(' EQ6 Reaction-Path Code (EQ/36-V',a,'-EQ6-EXE-',
     $ a,'-',a,')')
c
      write (noutpt,1020)
      write (nttyo,1020)
 1020 format(' Supported by the following EQ3/6 libraries:')
c
      i = index(uveeql,'.') - 1
      j = index(usteql,' ') - 1
      if (i .le. 0) i = 8
      if (j .le. 0) j = 8
      write (nttyo,1030) uveeql(1:i),usteql(1:j),uplatm(1:k)
      write (noutpt,1030) uveeql(1:i),usteql(1:j),uplatm(1:k)
 1030 format('   EQLIB (EQ3/6-V',a,'-EQLIB-LIB-',a,'-',a,')')
c
      i = index(uveelg,'.') - 1
      j = index(ustelg,' ') - 1
      if (i .le. 0) i = 8
      if (j .le. 0) j = 8
      write (nttyo,1040) uveelg(1:i),ustelg(1:j),uplatm(1:k)
      write (noutpt,1040) uveelg(1:i),ustelg(1:j),uplatm(1:k)
 1040 format('   EQLIBG (EQ3/6-V',a,'-EQLIBG-LIB-',a,'-',a,')')
c
      i = index(uveelu,'.') - 1
      j = index(ustelu,' ') - 1
      if (i .le. 0) i = 8
      if (j .le. 0) j = 8
      write (nttyo,1050) uveelu(1:i),ustelu(1:j),uplatm(1:k)
      write (noutpt,1050) uveelu(1:i),ustelu(1:j),uplatm(1:k)
 1050 format('   EQLIBU (EQ3/6-V',a,'-EQLIBU-LIB-',a,'-',a,')',/)
c
      write (nttyo,1060)
      write (noutpt,1060)
 1060 format(' Copyright (c) 1987, 1990-1993, 1995, 1997, 2002 The',
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
      write(noutpt,1080) utime0,udate0(1:j2)
      write(nttyo,1080) utime0,udate0(1:j2)
 1080 format(' Run',2(2x,a8,2x,a),//)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the platform's real*8 floating-point parameters.
c
      call flpars(eps100,irang,noutpt,nttyo,smp100)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize array dimension variables.
c
c     The following are common with EQ3NR.
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
c     The following are common with EQ3NR, but are only needed by that
c     code to write a full EQ6 input file as its pickup file.
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
c     The following are unique to EQ6.
c
      nsscmx = nssc_par
      nllnmx = nllnpa
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize the following floating-point constants:
c
c       rconst = the gas constant: 1.98726 cal/mole-K
c       rcnstv = the gas constant: 83.14510 bar-cm3/mol-K
c       vpgstp  = the volume of a perfect gas at STP: 22413.6 cm3
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
c     Read the header section of the DATA1 file. This section consists
c     of a record containing the string 'data1' (to ensure that the file
c     is indeed a DATA1 file), a record containing the keystring for the
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
c     Increase the values of some allocation size variables to allow
c     for some species, phases, etc., to created by this software.
c     Examples: generic ion exchanger phases and species.
c
      nbta_asv = nbta_asv + jet_par*net_par + 10
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
c     Allocate arrays to store the data read from the rest of the DATA1
c     file.
c
      ALLOCATE(iapxta(nxta_asv))
      ALLOCATE(ibpxta(nxta_asv))
      ALLOCATE(insgfa(nata_asv))
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
      ALLOCATE(axhfe(narx_asv,ntpr_asv))
      ALLOCATE(axlke(narx_asv,ntpr_asv))
      ALLOCATE(axvfe(narx_asv,ntpr_asv))
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
c     Read the remainder of the DATA1 file. The image of the data file
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
      ntfxmx = ntf1mx
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
c     Get the indices of H+, OH-, Cl-, fictive aqueous O2(g), and
c     fictive aqueous e-. This will be repeated after data compression.
c
c     Calling sequence substitutions:
c       narn1a for narn1
c       narn2a for narn2
c       nchloa for nchlor
c       neleca for nelect
c       nhydra for nhydr
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
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Setup up the remaining array allocation variables required to
c     allocate variables needed to store the dta to be read from the
c     input file.
c
c     The following value of k_asv is special to EQ6.
c     A smaller value is adequate for EQ3NR.
c
      k_asv = 2*nbta_asv + 5
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Allocate arrays to store the data to be read from the input
c     file.
c
      ALLOCATE(jflgi(nbta_asv))
      ALLOCATE(nbaspi(nbta_asv))
c
      ALLOCATE(mtbaqi(nbta_asv))
      ALLOCATE(mtbi(nbta_asv))
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
      nmx_asv = 3*nmut_asv
      nsx_asv = 2*nslt_asv
c
      npet_asv = nbt_asv + 5
      nset_asv = 20*npet_asv
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
      npetmx = npet_asv
      nptmax = npt_asv
      nsetmx = nset_asv
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
c
      ALLOCATE(insgf(nat_asv))
c
      ALLOCATE(nbasp(nbt_asv))
      ALLOCATE(nbaspd(nbt_asv))
      ALLOCATE(nbaspx(nbt_asv))
      ALLOCATE(nbmap(nbt_asv))
      ALLOCATE(ndecsp(nbt_asv))
c
      ALLOCATE(ncmap(nct_asv))
c
      ALLOCATE(ndrs(ndrs_asv))
      ALLOCATE(ndrsd(ndrs_asv))
      ALLOCATE(ndrsx(ndrs_asv))
c
      ALLOCATE(ness(ness_asv))
c
      ALLOCATE(jpflag(npt_asv))
      ALLOCATE(ncmpr(2,npt_asv))
c
      ALLOCATE(jflag(nst_asv))
      ALLOCATE(jflagd(nst_asv))
      ALLOCATE(jsflag(nst_asv))
      ALLOCATE(jsitex(nst_asv))
      ALLOCATE(nphasx(nst_asv))
      ALLOCATE(nsmap(nst_asv))
c
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
      ALLOCATE(utitl1(ntit_asv))
      ALLOCATE(utitl2(ntit_asv))
c
      ALLOCATE(ubmtbi(nbt_asv))
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
      ALLOCATE(apx(iapx_asv,nxt_asv))
      ALLOCATE(bpx(ibpx_asv,nxt_asv))
      ALLOCATE(wfac(ikt_asv,nxt_asv))
c
      ALLOCATE(loph(npt_asv))
      ALLOCATE(moph(npt_asv))
c
      ALLOCATE(azero(nat_asv))
      ALLOCATE(atwt(nct_asv))
c
      ALLOCATE(zvclg1(k_asv))
      ALLOCATE(zvec1(k_asv))
c
      ALLOCATE(mtb(nbt_asv))
      ALLOCATE(mtbaq(nbt_asv))
c
      ALLOCATE(mte(nct_asv))
      ALLOCATE(mteaq(nct_asv))
c
      ALLOCATE(cdrs(ndrs_asv))
      ALLOCATE(cdrsd(ndrs_asv))
      ALLOCATE(cdrsx(ndrs_asv))
      ALLOCATE(losp(nst_asv))
      ALLOCATE(mosp(nst_asv))
      ALLOCATE(mwtsp(nst_asv))
      ALLOCATE(vosp0(nst_asv))
      ALLOCATE(zchar(nst_asv))
      ALLOCATE(zchcu6(nst_asv))
      ALLOCATE(zchsq2(nst_asv))
c
      ALLOCATE(cess(ness_asv))
      ALLOCATE(tf1(ntf1_asv))
      ALLOCATE(tf2(ntf2_asv))
      ALLOCATE(xlkmod(nxmd_asv))
c
      ALLOCATE(dadhh(ipch_asv))
      ALLOCATE(dadhv(ipcv_asv))
      ALLOCATE(dbdhh(ipch_asv))
      ALLOCATE(dbdhv(ipcv_asv))
      ALLOCATE(dbdth(ipch_asv))
      ALLOCATE(dbdtv(ipcv_asv))
c
      ALLOCATE(dhfe(ipch_asv))
      ALLOCATE(dvfe(ipcv_asv))
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Allocate standard-state data arrays for module mod6xf.
c     The coefficient arrays here correspond to the 'a' set
c     arrays read from the data file. The present coefficient
c     arrays (and the other arrays) correspond to a compressed
c     set of species.
c
      ALLOCATE(xhfs(nst_asv))
      ALLOCATE(xhfsd(nst_asv))
      ALLOCATE(xlks(nst_asv))
      ALLOCATE(xlksd(nst_asv))
      ALLOCATE(xvfs(nst_asv))
      ALLOCATE(xvfsd(nst_asv))
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
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Allocate Pitzer data arrays for module mod6pt.
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
      ALLOCATE(pslamn(0:ipbt_asv,nslt_asv))
c
      ALLOCATE(gpit(ipbt_asv,nap_asv))
      ALLOCATE(dgpit(2,ipbt_asv,nap_asv))
      ALLOCATE(palpha(ipbt_asv,nap_asv))
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
c     Allocate still more arrays.
c
      ALLOCATE(iffg(nffg_asv))
      ALLOCATE(jffg(nffg_asv))
      ALLOCATE(npchk(npt_asv))
c
      ALLOCATE(cscale(nst_asv))
c
      ALLOCATE(jcode(nrct_asv))
      ALLOCATE(jreac(nrct_asv))
      ALLOCATE(nsk(nrct_asv))
      ALLOCATE(nrndex(nrct_asv))
      ALLOCATE(nxridx(nrct_asv))
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
      ALLOCATE(ndac(ndct_asv,imch_asv,2,nrct_asv))
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
      ALLOCATE(cbsr(nbt1_asv,nsrt_asv))
      ALLOCATE(cbsri(nbt1_asv,nsrt_asv))
      ALLOCATE(cdac(ndct_asv,imch_asv,2,nrct_asv))
      ALLOCATE(cesr(nct_asv,nsrt_asv))
      ALLOCATE(cesri(nct_asv,nsrt_asv))
      ALLOCATE(csigma(imch_asv,2,nrct_asv))
      ALLOCATE(eact(imch_asv,2,nrct_asv))
      ALLOCATE(elecsr(nsrt_asv))
      ALLOCATE(hact(imch_asv,2,nrct_asv))
      ALLOCATE(rk(imch_asv,2,nrct_asv))
      ALLOCATE(rkb(imch_asv,2,nrct_asv))
      ALLOCATE(trkb(imch_asv,2,nrct_asv))
      ALLOCATE(fkrc(nrct_asv))
      ALLOCATE(modr(nrct_asv))
      ALLOCATE(morr(nrct_asv))
      ALLOCATE(mwtrc(nrct_asv))
      ALLOCATE(sfcar(nrct_asv))
      ALLOCATE(ssfcar(nrct_asv))
      ALLOCATE(vreac(nrct_asv))
      ALLOCATE(moffg(nffg_asv))
      ALLOCATE(xlkffg(nffg_asv))
      ALLOCATE(mprph(npt_asv))
      ALLOCATE(mprphi(nprp_asv))
      ALLOCATE(mprsp(nst_asv))
      ALLOCATE(mprspi(nprs_asv))
      ALLOCATE(ptk(nptk_asv))
      ALLOCATE(ttk(nttk_asv))
      ALLOCATE(rxbar(ikt_asv,nxrt_asv))
      ALLOCATE(rxbari(ikt_asv,nxrt_asv))
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
      call initiz(iopt,noptmx)
      call initiz(iopg,nopgmx)
      call initiz(iopr,noprmx)
      call initiz(iodb,nodbmx)
c
      call initcb(uzveci,kmax)
      call initcb(uzvec1,kmax)
      call initaz(zvclg1,kmax)
      call initaz(zvec1,kmax)
c
c     Note: cspi and uspi are not used in EQ6.
c
      call initiz(jflgi,nbtmax)
c
      call initcb(utitl1,ntitmx)
      call initcb(utitl2,ntitmx)
c
      call initcb(ubmtbi,nbtmax)
c
      nmax = 2*nbtmax
      call initcb(usbsw,nmax)
      call initcb(uobsw,nmax)
c
      nmax = nbt_asv
      call initiz(nbasp,nmax)
      call initiz(nbaspd,nmax)
      call initiz(nbaspx,nmax)
c
      nmax = nst_asv
      call initiz(jflag,nmax)
      call initiz(jflagd,nmax)
      call initiz(jsflag,nmax)
      call initcb(uspec,nmax)
      call initaz(vosp0,nmax)
      call initaz(xlks,nmax)
      call initaz(zchar,nmax)
      call initaz(zchsq2,nmax)
      call initaz(zchcu6,nmax)
c
cxxxxxxxxxxxx
      call initaz(cdrs,ndrsmx)
      call initaz(cdrsx,ndrsmx)
      call initiz(ndrs,ndrsmx)
      call initiz(ndrsx,ndrsmx)
c
cxxxxxxxxxxxx
      nmax = 2*nst_asv
      call initiz(ndrsr,nmax)
      call initiz(ndrsrx,nmax)
c
      nmax = narx_asv*ntpr_asv*nst_asv
      call initaz(axlks,nmax)
      call initaz(axhfs,nmax)
      call initaz(axvfs,nmax)
cxxxxxxxxxxxx
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
cxxxxxxxxxxxx
c
      call initiz(iffg,nffgmx)
      call initiz(jffg,nffgmx)
c
      call initcb(uffg,nffgmx)
      call initaz(moffg,nffgmx)
      call initaz(xlkffg,nffgmx)
c
      call initiz(jcode,nrctmx)
      call initiz(jreac,nrctmx)
      call initiz(nsk,nrctmx)
c
      call initcb(ureac,nrctmx)
      call initaz(fkrc,nrctmx)
      call initaz(modr,nrctmx)
      call initaz(morr,nrctmx)
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
      call initcb(ugexp,netmax)
c
      do ne = 1,netmax
        jgext(ne) = 0
        kern1(ne) = 0
        kern2(ne) = 0
      enddo
c
      nmax = jetmax*netmax
      call initiz(jern1,nmax)
      call initiz(jern2,nmax)
      call initiz(ngext,nmax)
      call initiz(ngexrt,nmax)
      call initaz(egexjf,nmax)
c
      nmax = ietmax*jetmax*netmax
      call initaz(cegexs,nmax)
      call initaz(cpgexs,nmax)
      call initaz(mrgexs,nmax)
      call initiz(ngexro,nmax)
      call initiz(ngexso,nmax)
      call initiz(ngexsa,nmax)
c
      nmax = ietmax*jetmax*nertmx
      call initaz(egersi,nmax)
      call initaz(xgersi,nmax)
c
      nmax = ketmax*netmax
      call initiz(kgexsa,nmax)
c
      nmax = iktmax*nxrtmx
      call initaz(rxbari,nmax)
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
      read (ninpts,1090,end=105,err=107) ux8
 1090 format(a8)
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
        call rd6inw(awmaxi,awmini,cbsri,cdac,cesri,cgexj,
     $  csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,
     $  dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,
     $  dlxplo,dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,
     $  iact,ibsrti,iesrti,ietmax,iktmax,imchmx,imech,iodb,iopg,
     $  iopr,iopt,igerti,itermx,ixrti,jcode,jgerti,jetmax,jflgi,
     $  jgext,jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,
     $  ksppmx,kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,
     $  mtbaqi,mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,
     $  nert,nertmx,net,netmax,nffg,nffgmx,ngexrt,ninpts,nobswt,
     $  nodbmx,nopgmx,noprmx,noptmx,nordmx,noutpt,nprob,nprpmx,
     $  nprpti,nprsmx,nprsti,nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,
     $  nsrt,nsrtmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,nttyo,nxmdmx,
     $  nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrt,nxrtmx,o2maxi,o2mini,
     $  phmaxi,phmini,pressb,pressi,ptk,qend,qgexsh,qrderr,rkb,
     $  rxbari,sfcar,ssfcar,tempcb,tempci,tgexp,timmxi,tistti,tolbt,
     $  toldl,tolsat,tolxsf,trkb,ttk,ubmtbi,ubsri,ucxri,udac,uesri,
     $  uffg,ugerji,ugermo,ugersi,ugexj,ugexmo,ugexp,ugexr,uhfgex,
     $  uobsw,uprphi,uprspi,ureac,usbsw,utitl1,utitl2,uvfgex,uxcat,
     $  uxkgex,uxmod,uxopex,uxopt,uzveci,vreac,xgersi,xhfgex,ximaxi,
     $  xistti,xlkffg,xlkgex,xlkmod,xvfgex,zgexj,zvclgi)
      else
c
c       Menu-style (D) format.
c
        call rd6ind(awmaxi,awmini,cbsri,cdac,cesri,cgexj,
     $  csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,
     $  dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,
     $  dlxplo,dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,
     $  iact,ibsrti,iesrti,ietmax,iktmax,imchmx,imech,iodb,iopg,
     $  iopr,iopt,igerti,itermx,ixrti,jcode,jgerti,jetmax,jflgi,
     $  jgext,jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,
     $  ksppmx,kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,
     $  mtbaqi,mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,
     $  nert,nertmx,net,netmax,nffg,nffgmx,ngexrt,ninpts,nobswt,
     $  nodbmx,nopgmx,noprmx,noptmx,nordmx,noutpt,nprob,nprpmx,
     $  nprpti,nprsmx,nprsti,nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,
     $  nsrt,nsrtmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,nttyo,nxmdmx,
     $  nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrt,nxrtmx,o2maxi,o2mini,
     $  phmaxi,phmini,pressb,pressi,ptk,qend,qgexsh,qrderr,rkb,
     $  rxbari,sfcar,ssfcar,tempcb,tempci,tgexp,timmxi,tistti,tolbt,
     $  toldl,tolsat,tolxsf,trkb,ttk,ubmtbi,ubsri,ucxri,udac,uesri,
     $  uffg,ugerji,ugermo,ugersi,ugexj,ugexmo,ugexp,ugexr,uhfgex,
     $  uobsw,uprphi,uprspi,ureac,usbsw,utitl1,utitl2,uvfgex,uxcat,
     $  uxkgex,uxmod,uxopex,uxopt,uzveci,vreac,xgersi,xhfgex,ximaxi,
     $  xistti,xlkffg,xlkgex,xlkmod,xvfgex,zgexj,zvclgi)
      endif
c
      go to 109
  105 qend = .true.
      go to 109
  107 qrderr = .true.
  109 continue
c
      if (qrderr) then
        write (noutpt,1100)
        write (nttyo,1100)
 1100   format(/' * Error - (EQ6/eq6) Encountered an error while',
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
        write (noutpt,1110)
        write (nttyo,1110)
 1110   format(/' No further input found.',/)
c
c       Make porting changes in the EQLIBU routines that are called in
c       this section. Do not make the porting changes here.
c
c       Get end time and date. Also get the run, user, and cpu times.
c
        call runtim(iexec0,jexec0,texec0,noutpt,nttyo,trun,
     $  tuser,tcpu,udate1,utime1)
c
        j2 = ilnobl(udate0)
        j3 = ilnobl(udate1)
        write (noutpt,1120) utime0,udate0(1:j2),utime1,udate1(1:j3)
        write (nttyo,1120) utime0,udate0(1:j2),utime1,udate1(1:j3)
 1120   format(10x,'Start time = ',a8,2x,a,/12x,
     $  'End time = ',a8,2x,a)
c
c       Print the run, user, and cpu times.
c
        write (noutpt,1140) trun
        write (nttyo,1140)  trun
 1140   format(/10x,' Run time = ',g10.3,' seconds')
        if (tuser .gt. 0.) then
          write (noutpt,1150) tuser
          write (nttyo,1150)  tuser
 1150     format(10x,'User time = ',g10.3,' seconds')
        endif
        if (tcpu .gt. 0.) then
          write (noutpt,1160) tcpu
          write (nttyo,1160)  tcpu
 1160     format(10x,' Cpu time = ',g10.3,' seconds')
        endif
c
        write (noutpt,1170)
        write (nttyo,1170)
 1170   format(/' Normal exit')
c
c       BEGIN_MACHINE_DEPENDENT_CODE
c
c         Clear the IEEE flag for floating-point underflow, if such a
c         flag is present, to avoid getting an unnecessary system
c         warning message. Underflow is a normal condition in EQ3/6.
c         Make porting changes in the EQLIBU subroutine that is called
c         in this section. Do not make the porting changes here.
c
          call cliefu()
c
c       END_MACHINE_DEPENDENT_CODE
c
c       Close and delete the stripped input file and the tabs file.
c
        close (ninpts,status='delete')
        if (ntabs .gt. 0) close (ntabs,status='delete')
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
        do n = 1,ntitl1
          ux80 = utitl1(n)
          call locase(ux80)
          i = index(ux80,'useoldpitzermu')
          if (i .gt. 0) then
            qhawep = .false.
            go to 106
          endif
        enddo
c
        do n = 1,ntitl2
          ux80 = utitl2(n)
          call locase(ux80)
          i = index(ux80,'useoldpitzermu')
          if (i .gt. 0) then
            qhawep = .false.
            if (ntitl1 .lt. ntitmx) then
              ntitl1 = ntitl1 + 1
              utitl1(ntitl1) = utitl2(n)
            else
              write (noutpt,1180)
              write (nttyo,1180)
 1180         format(/' * Note - (EQ6/eq6) Found the USEOLDPITZERMU',
     $        ' option string',/7x,'in the secondary input file title',
     $        ' but cannot add it',/7x,'to the main title. The option',
     $        ' will not carry over to any',/7x,'PICKUP file that may',
     $        ' be written at the end of the present run.')
            endif
            go to 106
          endif
        enddo
  106   continue
      endif
c
      if (.not.qhawep) then
        write (noutpt,1190)
        write (nttyo,1190)
 1190   format(/' * Note - (EQ6/eq6) Found the USEOLDPITZERMU',
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
        do n = 1,ntitl1
          ux80 = utitl1(n)
          call locase(ux80)
          i = index(ux80,'useoldpitzer75')
          if (i .gt. 0) then
            qpit75 = .true.
            go to 170
          endif
        enddo
c
        do n = 1,ntitl2
          ux80 = utitl2(n)
          call locase(ux80)
          i = index(ux80,'useoldpitzer75')
          if (i .gt. 0) then
            qpit75 = .true.
            if (ntitl1 .lt. ntitmx) then
              ntitl1 = ntitl1 + 1
              utitl1(ntitl1) = utitl2(n)
            else
              write (noutpt,1192)
              write (nttyo,1192)
 1192         format(/' * Warning - (EQ6/eq6) Found the USEOLDPITZER75',
     $        ' option string',/7x,'in the secondary input file title',
     $        ' but cannot add it',/7x,'to the main title. The option',
     $        ' will not carry over to any',/7x,'PICKUP file that may',
     $        ' be written at the end of the present run.')
            endif
            go to 170
          endif
        enddo
  170   continue
      endif
c
      if (qpit75) then
        write (noutpt,1193)
        write (nttyo,1193)
 1193   format(/' * Warning - (EQ6/eq6) Found the USEOLDPITZER75',
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
        do n = 1,ntitl1
          ux80 = utitl1(n)
          call locase(ux80)
          i = index(ux80,'writepitzerjtables')
          if (i .gt. 0) then
            qcwrpj = .true.
            go to 180
          endif
        enddo
c
        do n = 1,ntitl2
          ux80 = utitl2(n)
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
 1194   format(/' * Note - (EQ6/eq6) Found the WRITEPITZERJTABLES',
     $  ' option string',/7x,'in one of the input file titles.',
     $  ' Will calculate the Pitzer',/7x,"J(x) and J'(x) functions",
     $  ' for higher order electrostatic terms and',/7x,'write output',
     $  ' tables for both the Pitzer (1975) and Harvie (1981)',/7x,
     $  'approximations.')
      endif
c
c     Scan the input file title for the TURNOFFOPTIMIZER string.
c     If found, turn off the pre-Newton-Raphson optimizer.
c
      qoptmz = .true.
      do n = 1,ntitl1
        ux80 = utitl1(n)
        call locase(ux80)
        i = index(ux80,'turnoffoptimizer')
        if (i .gt. 0) then
          qoptmz = .false.
          go to 110
        endif
      enddo
c
      do n = 1,ntitl2
        ux80 = utitl2(n)
        call locase(ux80)
        i = index(ux80,'turnoffoptimizer')
        if (i .gt. 0) then
          qoptmz = .false.
          if (ntitl1 .lt. ntitmx) then
            ntitl1 = ntitl1 + 1
            utitl1(ntitl1) = utitl2(n)
          else
            write (noutpt,1195)
            write (nttyo,1195)
 1195       format(/' * Warning - (EQ6/eq6) Found the TURNOFFOPTIMIZER',
     $      ' option string',/7x,'in the secondary input file title',
     $      ' but cannot add it',/7x,'to the main title. The option',
     $      ' will not carry over to any',/7x,'PICKUP file that may',
     $      ' be written at the end of the present run.')
          endif
          go to 110
        endif
      enddo
  110 continue
c
      if (.not.qoptmz) then
        write (noutpt,1197)
        write (nttyo,1197)
 1197   format(/' * Note - (EQ6/eq6) Found the TURNOFFOPTIMIZER',
     $  ' option string',/7x,'in one of the input file titles.',
     $  ' Will turn off the',/7x,'pre-Newton-Raphson optimizer.',
     $  ' This may help if Newton-Raphson',/7x,'iteration fails',
     $  ' to converge due to insufficient optimizer',
     $  /7x,'performance. This is most likely to be a problem',
     $  ' when the',/7x,'phase assemblage changes.')
      endif
c
c     Scan the input file title for the TABFILEASTXT string. If found,
c     found, make the TAB file an ordinary text file (the old style
c     TAB file) instead of as a comma separated value (.csv) file.
c
      qtatxt = .false.
      do n = 1,ntitl1
        ux80 = utitl1(n)
        call locase(ux80)
        i = index(ux80,'tabfileastxt')
        if (i .gt. 0) then
          qtatxt = .true.
          go to 210
        endif
      enddo
c
      do n = 1,ntitl2
        ux80 = utitl2(n)
        call locase(ux80)
        i = index(ux80,'tabfileastxt')
        if (i .gt. 0) then
          qtatxt = .true.
          if (ntitl1 .lt. ntitmx) then
            ntitl1 = ntitl1 + 1
            utitl1(ntitl1) = utitl2(n)
          else
            write (noutpt,1196)
            write (nttyo,1196)
 1196       format(/' * Warning - (EQ6/eq6) Found the TABFILEASTXT',
     $      ' option string',/7x,'in the secondary input file title',
     $      ' but cannot add it',/7x,'to the main title. The option',
     $      ' will not carry over to any',/7x,'PICKUP file that may',
     $      ' be written at the end of the present run.')
          endif
          go to 210
        endif
      enddo
  210 continue
c
      if (qtatxt) then
        write (noutpt,1198)
        write (nttyo,1198)
 1198   format(/' * Note - (EQ6/eq6) Found the TABFILEASTXT',
     $  ' option string',/7x,'in one of the input file titles.',
     $  ' Will write the TAB',/7x,'file as an ordinary text file',
     $  ' instead of the default csv',/7x,'(comma-separated-variable)',
     $  ' file. The data content may not',/7x,'be the same for both',
     $  ' files.')
      endif
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
      if (iopt(16) .ge. 0) then
        inquire(file=bafile,opened=qop)
        if (.not.qop) call openou(noutpt,nttyo,bafile,'formatted',
     $  nrecl,nbkupa)
        if (iopt(16) .eq. 0) then
          inquire(file=bbfile,opened=qop)
          if (.not.qop) call openou(noutpt,nttyo,bbfile,'formatted',
     $    nrecl,nbkupb)
        endif
      endif
c
      if (iopt(17) .ge. 0) then
        inquire(file=pfile,opened=qop)
        if (.not.qop) call openou(noutpt,nttyo,pfile,'formatted',
     $  nrecl,newin)
      endif
c
      if (nprob.gt.1 .and. iopt(18).ge.1) iopt(18) = 0
      if (iopt(18) .ge. 1) then
        inquire(file=txfile,exist=qex)
        if (qex) then
          call openin(noutpt,nttyo,txfile,'formatted',ntabx)
          do i = 1,10000
            read (ntabx,1200,end=200) ux8
 1200       format(a8)
          enddo
  200     continue
        else
          write (noutpt,1210)
          write (nttyo,1210)
 1210     format(/' * Error - (EQ6/eq6) Have iopt(18)= 1, but there',
     $    /7x,'is no existing tabx file to which to append.')
          stop
        endif
      endif
c
      if (iopt(18) .ge. 0) then
c
c       Set up the TAB, TABS, and TABX files.
c
c       For the .csv form, the record length (nrecl) must match the
c       maximum character lengths of the table tag (8 characters)
c       and the ulinex variable (nllnmx characters), plus two
c       (for one comma trailing the table tag and another trailing
c       The ulinex variable).
c
        nrecl = nllnmx + 10
        if (qtatxt) nrecl = 128
c
        inquire(file=tfile,opened=qop)
        if (.not.qop) call openou(noutpt,nttyo,tfile,'formatted',
     $  nrecl,ntab)
c
        inquire(file=tsfile,opened=qop)
        if (.not.qop) then
          call openou(noutpt,nttyo,tsfile,'formatted',nrecl,ntabs)
        else
          rewind ntabs
        endif
c
        nrecl = nllnmx + 10
        if (qtatxt) nrecl = 129
c
        inquire(file=txfile,opened=qop)
        if (.not.qop) then
          call openou(noutpt,nttyo,txfile,'formatted',nrecl,ntabx)
        endif
        if (iopt(18) .eq. 0) rewind ntabx
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     qecon = .true. denotes economy mode.
c     qscon = .true. denotes super economy mode.
c
      qecon = .false.
      qscon = .false.
      if (iopt(1).lt.2 .and. iopt(2).le.0) then
        if (iopt(13) .ge. 1) qecon = .true.
        if (iopt(13) .ge. 2) qscon = .true.
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set default values of run control parameters read from the
c     data file. Put any out-of-range values into acceptable range.
c
      tolxst = 0.
      tolxsu = 0.
      call dfaltz(dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,
     $ dloplo,dloprn,dltplo,dltpll,dltprl,dltprn,dlxdmp,dlxmx0,
     $ dlxplo,dlxpll,dlxprl,dlxprn,iopt,itermx,ksplmx,ksppmx,kstpmx,
     $ net,noptmx,nordmx,nrct,noutpt,ntrymx,nttyo,prcinf,qecon,qscon,
     $ timmxi,tistti,tolbt,toldl,tolsat,tolxsf,tolxst,tolxsu,
     $ ximaxi,xistti)
      nrd1mx = nordmx + 1
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(2) .gt. 0) then
c
c       Disallow print and plot intervals in log time space if the
c       starting time is less than zero.
c
        if (tistti.lt.0. .or. timmxi.lt.0.) then
          dltprl = prcinf
          dltpll = prcinf
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set other run control parameters.
c
      call setrcp(aftarg,dlxmax,dlxmin,dlxmx0,npslmx,nsscmx,
     $ nsslmx,prcinf,sscrew,tolaft,tolsat,tolsst,zkfac,zklgmn,
     $ zklogl,zklogu)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set some parameters using the data read from the input file.
c
      nart = 0
      ngrt = 0
      nmrt = 0
      do nrc = 1,nrct
        if (jcode(nrc) .eq. 0) then
          nmrt = nmrt + 1
        elseif (jcode(nrc) .eq. 3) then
          nart = nart + 1
        elseif (jcode(nrc) .eq. 4) then
          ngrt = ngrt + 1
        endif
      enddo
c
      km1 = kbt + 1
      kx1 = kmt + 1
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set the starting and maximum values of reaction progress and
c     time to the input values. Note that the values of the input
c     variables will be changed prior to writing data on a backup
c     or pickup file.
c
      xistrt = xistti
      ximax = ximaxi
      tistrt = tistti
      timemx = timmxi
c
      xi1 = xistrt
      time1 = tistrt
c
c     Set the values of other run control parameters. These parameters
c     are also associated with mechanisms for terminating the current
c     reaction path; for example, if pH decreases to phmin, or increases
c     to phmax.
c
      phmin = phmini
      phmax = phmaxi
      ehmin = ehmini
      ehmax = ehmaxi
      o2min = o2mini
      o2max = o2maxi
      awmin = awmini
      awmax = awmaxi
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check the input reactant data for various kinds of errors
c     and inconsistencies.
c
      call chkinz(ier,imchmx,imech,iopt,jcode,kmax,kxt,nelect,
     $ noptmx,noutpt,no2gaq,nrct,nrctmx,nrk,nstmax,nttyo,rkb,ureac,
     $ uspeca,uzveci,zvclgi)
c
      if (ier .gt. 0) stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set some defaults.
c
      do nrc = 1,nrct
        if (fkrc(nrc) .le. 0.) fkrc(nrc) = 1.0
        if (sfcar(nrc) .lt. 0.) sfcar(nrc) = 0.
        if (ssfcar(nrc) .lt. 0.) ssfcar(nrc) = 0.
        do j = 1,2
          if (nrk(j,nrc) .eq. 2) then
            if (imech(j,nrc) .le. 0) imech(j,nrc) = 1
            do i = 1,imech(j,nrc)
              if (csigma(i,j,nrc) .le. 0.) csigma(i,j,nrc) = 1.
            enddo
          endif
        enddo
      enddo
c
      do n = 1,nsbswt
        if (usbsw(1,n)(1:24) .eq. ublk24(1:24))
     $  usbsw(1,n)(1:24) = uaqsln
        if (usbsw(2,n)(1:24) .eq. ublk24(1:24))
     $  usbsw(2,n)(1:24) = uaqsln
      enddo
c
      do n = 1,nobswt
        if (uobsw(1,n)(1:24) .eq. ublk24(1:24))
     $  uobsw(1,n)(1:24) = uaqsln
        if (uobsw(2,n)(1:24) .eq. ublk24(1:24))
     $  uobsw(2,n)(1:24) = uaqsln
      enddo
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
c     set will be redfined. This second stage form is equivalent to
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
c     Get the basis index of water (nbw).
c
c     Calling sequence substitutions:
c       nbaspd for nbasp
c       nbtd for nbt
c       narn1a for ns
c
      nbw = nbasis(nbaspd,nbtd,nbtmax,narn1a)
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
c     Initialize some arrays.
c
c     Note: csp and nsp are not used in EQ6.
c
      call initiz(ndecsp,nbtmax)
c
c     Note: mtb, mtbaq, mte, and mteaq are not used in EQ3NR.
c
      call initaz(mtb,nbtmax)
      call initaz(mtbaq,nbtmax)
      call initaz(mte,nctmax)
      call initaz(mteaq,nctmax)
c
c     Note: conc and conclg are initialized in EQ6/path.f.
c
      call initaz(moph,nptmax)
      call initaz(mosp,nstmax)
      av = -99999.
      call initav(loph,nptmax,av)
      call initav(losp,nstmax,av)
c
c     Note: mprsp is not used in EQ3NR.
c
      call initaz(mprsp,nstmax)
c
c     Note: acflg, acflgo, act, actlg, xbar, and xbarlg, are initialized
c     in EQ6/path.f.
c
c     Note: npchk and mprph are not used in EQ3NR.
c
      call initiz(npchk,nptmax)
      call initaz(mprph,nptmax)
c
c     Note: zvec0 and zvclg0 are initialized in EQ6/path.f.
c
      call initiz(iindx1,kmax)
      call initiz(ipndx1,kmax)
c
      call initaz(zvec1,kmax)
c
      av = -99999.
      call initav(zvclg1,kmax,av)
c
c     Note: the following are not used in EQ3NR.
c
      nmax = imchmx*2*nrctmx
      call initaz(rk,nmax)
c
      nmax = ndctmx*imchmx*2*nrctmx
      call initiz(ndac,nmax)
c
      nmax = iktmax*nxrtmx
      call initaz(rxbar,nmax)
c
      nmax = nctmax*nsrtmx
      call initaz(cesr,nmax)
c
      nmax = nbt1mx*nsrtmx
      call initaz(cbsr,nmax)
c
      nmax = ietmax*jetmax*nertmx
      call initaz(egers,nmax)
      call initaz(xgers,nmax)
      call initaz(mrgers,nmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Interpret the data file basis species listed on the input file.
c     "Data file" basis species to be created according to instructions
c     read from the input file (e.g., for generic ion exchangers)
c     will be ignored at this point, as creation occurs later in
c     this code.
c
      call intbs6(jflag,jflgi,kmax,narn1a,narn2a,nbaspd,nbtd,
     $ nbti,nbtmax,ndrsrd,ndecsp,noutpt,nsta,nstmax,nttyo,uspeca,
     $ ubmtbi)
c
c     Set jflag to -1 for species that can not appear in the system.
c
      call jflaux(jflag,nbaspd,nbtd,nbtmax,ndrsd,ndrsmx,
     $ ndrsrd,nstmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Test the temperature-tracking input.
c
      if (jtemp .eq. 0) then
        continue
      elseif (jtemp .eq. 1) then
        if (ttk(1) .eq. 0.) then
          write (noutpt,2110)
          write (nttyo,2110)
 2110     format(/' * Note (EQ6/eq6) The specified temperature',
     $    /7x,'derivative for computing the temperature as a linear',
     $    /7x,'function of Xi is zero. The temperature is therefore',
     $    /7x,'constant. Will reset the temperature-tracking flag',
     $    /7x,'jtemp from 1 to 0.')
          jtemp = 0
        endif
      elseif (jtemp .eq. 2) then
        if (ttk(1) .eq. 0.) then
          write (noutpt,2120)
          write (nttyo,2120)
 2120     format(/' * Note (EQ6/eq6) The specified temperature',
     $    /7x,'derivative for computing the temperature as a linear',
     $    /7x,'function of time is zero. The temperature is therefore',
     $    /7x,'constant. Will reset the temperature-tracking flag',
     $    /7x,'jtemp from 2 to 0.')
          jtemp = 0
        endif
      elseif (jtemp .eq. 3) then
        if (ttk(1) .eq. 0.) then
          write (noutpt,2130)
          write (nttyo,2130)
 2130     format(/' * Note (EQ6/eq6) The specified mass factor',
     $    /7x,'ratio (ttk(1) for computing the temperature from',
     $    /7x,'a fluid mixing ratio (jtemp = 3) was zero. This is',
     $    /7x,'not a valid value. Will change this to 1.0. This',
     $    /7x,'parameter affects the relation between the fluid',
     $    /7x,'mixing ratio and Xi. A value of 1.0 corresponds to',
     $    /7x,'a 1:1 mixing ratio.')
        endif
      else
        write (noutpt,2140) jtemp
        write (nttyo,2140) jtemp
 2140   format(/' * Error - (EQ6/eq6) The temperature tracking',
     $  /7x,'flag (jtemp) has an unknown value of ',i2,'.')
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Test the pressure-tracking input.
c
      if (jpress .eq. 0) then
        continue
      elseif (jpress .eq. 1) then
        continue
      elseif (jpress .eq. 2) then
        continue
      elseif (jpress .eq. 3) then
        if (ptk(1) .eq. 0.) then
          write (noutpt,2150)
          write (nttyo,2150)
 2150     format(/' * Note (EQ6/eq6) The specified pressure',
     $    /7x,'derivative for computing the pressure as a linear',
     $    /7x,'function of Xi is zero. The pressure is therefore',
     $    /7x,'constant. Will reset the pressure-tracking flag',
     $    /7x,'jpress from 3 to 2.')
          jpress = 0
        endif
      elseif (jpress .eq. 4) then
        if (ptk(1) .eq. 0.) then
          write (noutpt,2160)
          write (nttyo,2160)
 2160     format(/' * Note (EQ6/eq6) The specified pressure',
     $    /7x,'derivative for computing the pressure as a linear',
     $    /7x,'function of time is zero. The pressure is therefore',
     $    /7x,'constant. Will reset the pressure-tracking flag',
     $    /7x,'jpress from 4 to 2.')
          jpress = 0
        endif
      else
        write (noutpt,2170) jpress
        write (nttyo,2170) jpress
 2170   format(/' * Error - (EQ6/eq6) The pressure tracking',
     $  /7x,'flag (jpress) has an unknown value of ',i2,'.')
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the initial temperature.
c
      call gtemp(afcnst,al10,iopt,jtemp,noptmx,noutpt,nttkmx,
     $ nttyo,rconst,rtcnst,tempc,tempcb,tempk,time1,ttk,xi1)
c
c     Determine the corresponding temperature range flag.
c
      call gntpr(ntpr,ntprmx,ntprt,tempc,tempcu)
c
c     Check for constant temperature.
c
      qcntmp = jtemp .eq. 0
c
c     Check for a temperature jump.
c
      dt = tempc - tempci
      if (abs(dt) .gt. 1.e-4) then
        write (noutpt,1300) tempc,tempci
        write (nttyo,1300) tempc,tempci
 1300   format(/' * Note - (EQ6/eq6) The temperature is jumping from',
     $  /7x,'the previous run:',/,
     $  /'     Initial temperature      = ',f10.4,' C',
     $  /'     Previous run temperature = ',f10.4,' C')
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the data file reference grid pressure at the initial
c     temperature.
c
c     Calling sequence substitutions:
c       apresg for arr
c       presg for prop
c
      call evdat2(apresg,narxmx,narxt,ntpr,ntprmx,presg,tempc)
c
c     Compute the 1.013-bar/steam-saturation curve pressure at the
c     initial temperature.
c
      if (tempc .le. 100.) then
        ntprh = 1
      else
        ntprh = 2
      endif
c
c     Calling sequence substitutions:
c       apresh for arr
c       5 for narxmx
c       narxth for narxt
c       ntprh for ntpr
c       2 for ntprmx
c       presh for prop
c
      call evdat2(apresh,5,narxth,ntprh,2,presh,tempc)
c
c     Compute the initial pressure.
c
      call gpress(iopt,jpress,noptmx,noutpt,nptkmx,nttyo,presg,
     $ presh,press,pressb,time1,ptk,xi1)
c
c     Check for constant pressure.
c
      qcnpre = jpress.eq.2 .or. (jpress.eq.0 .and. qcntmp) .or.
     $ (jpress.eq.1 .and. qcntmp)
c
c     Check for a pressure jump.
c
      dp = press - pressi
      if (abs(dp) .gt. 1.e-4) then
        write (noutpt,1310) press,pressi
        write (nttyo,1310) press,pressi
 1310   format(/' * Note - (EQ6/eq6) The pressure is jumping from',
     $  /7x,'the previous run:',/,
     $  /'     Initial pressure      = ',1pg12.5,' bars',
     $  /'     Previous run pressure = ',1pg12.5,' bars')
      endif
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
      if (iopt(15) .ge. 1) then
c
c       Execute the option to suppress all redox reactions.
c
        call suprdx(jflag,jsflag,narn1a,narn2a,ndrsd,ndrsmx,
     $  ndrsrd,nrdxsp,nsta,nstmax,uspeca)
      endif
c
c     Execute the pure mineral subset selection suppression options.
c
      call mincsp(cdrsd,jpflag,jsflag,nbaspd,nbtd,nbtmax,ncmpra,
     $ ndrsd,ndrsmx,ndrsrd,nmrn1a,nmrn2a,noutpt,npta,nptmax,nstmax,
     $ nttyo,nxopex,nxopmx,nxopt,nxpemx,uspeca,uxcat,uxopex,uxopt)
c
c     Check the jpflag array to make sure it is consistent with the
c     jsflag array.
c
      call flgchk(jpflag,jsflag,ncmpra,npta,nptmax,nstmax,qclnsa)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Examine each active auxiliary basis species. Print a warning if
c     any other species in the corresponding dissociation reaction is
c     not present in the model.
c
      call bspchk(jsflag,nbaspd,nbtd,nbtmax,ndrsd,ndrsmx,ndrsrd,
     $ noutpt,nrdxsp,nstmax,nttyo,uspeca)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Do data array compression. Write working data arrays that
c     don't include phases and species that aren't necessary
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
c     Examine each active auxiliary basis species. Change the
c     corresponding log K polynomial coefficients so that log K is
c     fixed at a value of -99999 if any other species in the
c     corresponding dissociation reaction is not in the model.
c
      call bsplkp(axlks,narxmx,nbasp,nbt,nbtmax,ndrs,ndrsmx,
     $ ndrsr,nstmax,ntprmx)
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
c     Alter any log K values as directed by the input file
c     (actually, it is the set of interpolating polynomial
c     coefficients which is altered).
c
      if (nxmod .gt. 0) then
        call alters(afcnst,apresg,axlks,cdrs,kxmod,narxmx,narxt,
     $  ndrs,ndrsmx,ndrsr,noutpt,npt,nptmax,nst,nstmax,ntpr,ntprmx,
     $  nttyo,nxmdmx,nxmod,tempc,uphase,uspec,uxmod,xlkmod)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Create fictive minerals for the fixed fugacity option.
c
      nfrn1 = nst + 1
      nfrn2 = nst
      ifrn1 = npt + 1
      ifrn2 = npt
      if (nffg .gt. 0) then
        call nlkffg(axlks,cess,cdrs,iffg,ifrn1,ifrn2,jffg,jpflag,
     $  jsflag,mwtsp,narxmx,ncmpr,ndrs,ndrsmx,ndrsr,ness,nessmx,nessr,
     $  nffg,nffgmx,nfrn1,nfrn2,ngrn1,ngrn2,noutpt,nphasx,npt,nptmax,
     $  nst,nstmax,ntpr,ntprmx,nttyo,qcntmp,uffg,ufixf,uphase,uspec,
     $  vosp0,xlkffg)
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
c     Determine if the model to be calculated has a redox aspect.
c     Note: the iopt(15) option may be been excecuted above.
c
      call tstrdx(cdrs,iodb,iopt,jflag,jsflag,narn1,narn2,
     $ ndrs,ndrsmx,ndrsr,nodbmx,noptmx,noutpt,nrdxsp,nstmax,
     $ qredox,uspec)
c
c     If there is no redox aspect, make sure that all active redox
c     reactions are suppressed.
c
      if (.not.qredox .and. iopt(15).le.0) then
c
c       Calling sequence substitutions:
c         narn1 for narn1a
c         narn2 for narn2a
c         ndrs for ndrsd
c         ndrsr for ndrsrd
c         nst for nsta
c         uspec for uspeca
c
        call suprdx(jflag,jsflag,narn1,narn2,ndrs,ndrsmx,
     $  ndrsr,nrdxsp,nst,nstmax,uspec)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Interpret the mass balance totals. Construct the mtb and mtbaq
c     arrays.
c
      call intmtb(mtb,mtbaq,mtbaqi,mtbi,nbasp,nbt,nbti,nbtmax,
     $ noutpt,nstmax,nttyo,ubmtbi,uspec)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Interpret matrix variables. Construct the iindx1 and zvclg1
c     arrays. Note that new fictive fugacity-fixing phases will not
c     yet appear in the matrix after this has been done.
c
      call intmat(iaqsln,iindx1,ipndx1,kbt,kdim,kelect,khydr,
     $ khydx,kmax,km1,kmt,ko2gaq,kwater,kx1,kxt,narn1,narn2,nbasp,
     $ nbt,nbti,nbtmax,ncmpr,nelect,nern1,nern2,nhydr,nhydx,nobswt,
     $ noutpt,no2gaq,nphasx,npt,nptmax,nstmax,nttyo,qloffg,ubmtbi,
     $ ufixf,uobsw,uphase,uspec,uzveci,uzvec1,zvclgi,zvclg1,zvec1)
c
       krdxsp = ko2gaq
       nrdxsp = no2gaq
       if (krdxsp .eq. 0) then
         krdxsp = kelect
         nrdxsp = nelect
       endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nsrt .gt. 0) then
c
c       Construct reactions for special reactants in cases for
c       which a reaction was not read from the input file.
c
        call makrsr(cbsri,cesri,cess,eps100,ibsrti,iesrti,jcode,
     $  nbt1mx,nct,nctmax,ness,nessmx,nessr,noutpt,nrct,nrctmx,nsrtmx,
     $  nstmax,nttyo,ubsri,uelem,uesri,ureac,uspec,zchar)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Interpret reactant names and associated properties read from the
c     input file.
c
      call intrct(cbsr,cbsri,cesr,cesri,egers,egersi,ibsrti,
     $ iern1,iern2,iesrti,ietmax,igerti,iktmax,imrn1,imrn2,ixrn1,
     $ ixrn2,ixrti,jcode,jetmax,jgerti,jgext,narn1,narn2,nbaspd,
     $ nbt,nbtmax,nbt1mx,ncmpr,nct,nctmax,nertmx,netmax,ngexsa,
     $ ngext,ngrn1,ngrn2,noutpt,nptmax,nrct,nrctmx,nrndex,nsrtmx,
     $ nstmax,nttyo,nxridx,nxrtmx,rxbar,rxbari,ubsri,ucxri,uelem,
     $ uesri,ugerji,ugermo,ugersi,ugexj,ugexmo,uphase,ureac,uspec,
     $ xgers,xgersi)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nert .gt. 0) then
c
c       Calculate the mole ratios (mrgers) for species on exchange
c       sites of generic ion exchanger reactants. The mole ratio is the
c       moles of exchanger species (specific to a site) per mole of
c       exchanger substrate. These mole ratios are used to compute mass
c       balance increments. They are calculated from the corresponding
c       equivalent fractions (egers). Note that these mole ratios and
c       equivalent fractions are model-independent, unlike the
c       corresponding mole fractions.
c
        do nrc = 1,nrct
          if (jcode(nrc) .eq. 5) then
            ner = nxridx(nrc)
            np = nrndex(nrc)
            ne = np - iern1 + 1
            do je = 1,jgext(ne)
              ex = zgexj(je,ne)*cgexj(je,ne)
              do ie = 1,ngext(je,ne)
                nss = ngexsa(ie,je,ne)
                if (nss .gt. 0) then
                  mrgers(ie,je,ner) = -ex*egers(ie,je,ner)/zchar(nss)
                endif
              enddo
            enddo
          endif
        enddo
c
c       Calculate the mole fractions (xgers) for species on exchange
c       sites of generic ion exchanger reactants.
c
        do nrc = 1,nrct
          if (jcode(nrc) .eq. 5) then
            ner = nxridx(nrc)
            np = nrndex(nrc)
            ne = np - iern1 + 1
            do je = 1,jgext(ne)
              xx = 0.
              do ie = 1,ngext(je,ne)
                xx = xx  + mrgers(ie,je,ner)
              enddo
              do ie = 1,ngext(je,ne)
                xgers(ie,je,ner) = mrgers(ie,je,ner)/xx
              enddo
            enddo
          endif
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nsrt .gt. 0) then
c
c       Make sure that reactions for special reactants are charge
c       balanced.
c
        call chzrsr(cbsr,elecsr,eps100,jcode,nbasp,nbt,nbtmax,
     $  nbt1mx,noutpt,nrct,nrctmx,nrndex,nsrtmx,nstmax,nttyo,ureac,
     $  uspec,zchar)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Find any basis species in the matrix with jflag = 30 and
c     eliminate them from the matrix. Adjust as necessary the total
c     masses of those basis species which remain in the matrix.
c
      call combmb(cdrs,iindx1,ipndx1,jflag,kbt,kdim,km1,kmax,
     $ kmt,kx1,kxt,mtb,mtbaq,ndrsmx,nbasp,nbt,nbtmax,ndrs,ndrsr,
     $ noutpt,nstmax,nttyo,uspec,uzvec1,zvclg1,zvec1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize the error counter.
c
      nerr = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make sure that no mineral species in the equilibrium system (ES)
c     has been suppressed.
c
      do kcol = km1,kxt
        np = ipndx1(kcol)
        ns = iindx1(kcol)
        if (jsflag(ns) .gt. 0) then
c
c         Calling sequence substitutions:
c           uspec(ns) for unam48
c
          call fmspnm(jlen,uspec(ns),uspn56)
          write (noutpt,1350) uspn56(1:jlen)
          write (nttyo,1350) uspn56(1:jlen)
 1350     format(/' * Error - (EQ6/eq6) The species ',a," can't be",
     $    /7x,'suppressed because it is present in the equilibrium',
     $    ' system.')
          nerr = nerr + 1
        else
          jpflag(np) = -1
          jsflag(ns) = -1
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check to see that no species in the strict basis has a jflag
c     value of 30 at this point.
c
      do nb = 1,nbt
        ns = nbasp(nb)
        nt = ndrsr(2,ns) - ndrsr(1,ns) + 1
        if (nt .lt. 2) then
          if (jflag(ns) .eq. 30) then
            ux = uspec(ns)(1:24)
            j2 = ilnobl(ux)
            write (ux8,'(i5)') jflag(ns)
            call lejust(ux8)
            j3 = ilnobl(ux8)
            write (noutpt,1360) uspec(ns)(1:j2),ux8(1:j3)
            write (nttyo,1360) uspec(ns)(1:j2),ux8(1:j3)
 1360       format(/' * Error - (EQ6/eq6) The species ',a," can't have",
     $      /7x,'a jflag value of ',a,', because this species is in',
     $      ' the strict basis set.')
            nerr = nerr + 1
          endif
        endif
      enddo
c
      if (nerr .gt. 0) stop
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
c     Compute affinity scaling factors. Note that these are
c     defined in terms of the reactions before any eliminations
c     or ordinary basis switches are made.
c
      call gafscl(cdrsd,cscale,ndrsmx,ndrsrd,nst,nstmax)
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
c         Write tables concerning Pitzer interaction coefficients.
c
          call ptztab(iopr,narn1,narn2,natmax,nmutmx,nmux,nmxi,
     $    nmxmax,nmxx,noprmx,noutpt,nsltmx,nslx,nstmax,nsxi,nsxmax,
     $    nsxx,uspec)
        endif
c
c       Print warnings for species lacking Pitzer coefficients.
c
        call ptzchk(narn1,narn2,natmax,nmxi,noutpt,nstmax,nsxi,
     $  nttyo,uspec)
c
c       Transform conventional mu data to corresponding C, psi,
c       and zeta data (data originally defined in mu form is
c       not affected).
c
cxxxxxxxxxxx
        if (qhawep) then
          call rc3ocf(amu,jpfcmx,ifcphi1,ifcphi2,ifnnn,ifn2n,
     $    ifpsi1,ifpsi2,ifzeta,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,
     $    ilpsi2,ilzeta,iodb,nmux,nmut,nmutmx,nodbmx,noutpt,nstmax,
     $    nttyo,uspec,zchar)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the thermodynamic parameters that are functions of
c     temperature.
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
c     Compute the rate parameters that are functions of temperature.
c
      call evratc(eact,hact,iact,imchmx,imech,nrct,nrctmx,nrk,
     $ rk,rkb,rtcnst,tempk,trkb)
c
      tempcd = tempc
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
          write (noutpt,1370) press,presg,dp
          write (nttyo,1370) press,presg,dp
 1370     format(/' * Warning - (EQ6/eq6) The supporting data file',
     $    /7x,'contains no data to support making thermodynamic',
     $    /7x,'pressure corrections. No such corrections will be made.',
     $    /7x,'The current pressure is ',1pg12.5,' bars, the standard',
     $    /7x,'grid pressure is ',g12.5,' bars, and the pressure',
     $    /7x,'difference is ',g12.5,' bars.')
        else
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
      pressd = press
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1400)
 1400 format(//21x,'--- Inactive Species ---',/)
c
      k = 0
      klim = 20
      do ns = 1,nst
        if (jsflag(ns) .eq. 1) then
          if (ns .ne. no2gaq) then
            k = k + 1
            if (k .le. klim) then
c
c             Calling sequence substitutions:
c               uspec(ns) for unam48
c
              call fmspnm(jlen,uspec(ns),uspn56)
              write (noutpt,1410) uspn56(1:jlen)
 1410         format(4x,a)
            endif
          endif
        endif
      enddo
c
      if (k .gt. klim) then
        i = k - klim
        write (noutpt,1420) i
 1420   format(/6x,'plus ',i4,' others',/)
      elseif (k .le. 0) then
        write (noutpt,1430)
 1430   format(4x,'None',/)
      else
        write (noutpt,1440)
 1440   format(1x)
      endif
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
            mtb(nb) = 0.
            mtbaq(nb) = 0.
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
c
        do kcol = 1,kbt
          nb = iindx1(kcol)
          ns = nbasp(nb)
          if (ns.ge.nern1 .and. ns.le.nern2) uzvec1(kcol) = uspec(ns)
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nsrt .gt. 0) then
c
c       Check the reactions for the special reactants. Rewrite these
c       reactions to eliminate any basis species for which jflag = 30.
c
        call ckfrsr(cbsr,csts,jcode,jflag,nbaspd,nbtd,nbtmax,
     $  nbt1mx,noutpt,nrct,nrctmx,nrndex,nsrtmx,nstmax,nsts,nstsmx,
     $  nstsr,nttyo,ureac,uspec)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(5) .gt. 0) then
c
c       Clear the ES solids read from the input file. Fictive
c       fugacity-fixing minerals are not cleared.
c
        write (noutpt,1470)
        write (nttyo,1470)
 1470   format(/' * Note (EQ6/eq6) Clearing equilibrium system',
     $  ' (ES) solids',/7x,'read from the input file.')
c
        call clress(csts,iindx1,ipndx1,jpflag,jsflag,kdim,kmax,
     $  km1,kmt,kx1,kxt,loph,losp,moph,mosp,mtb,mtbaq,nbt,nbtmax,
     $  nptmax,nstmax,nsts,nstsmx,nstsr,ufixf,uzvec1,zvec1,zvclg1)
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
c     Put new fixed fugacity minerals into the matrix if the
c     corresponding added masses are positive. Keep old ones
c     in if the corresponding net masses are positive. If necessary,
c     recalculate the mass balance totals. Note: this call is
c     necessary even if nffg = 0, because it completes the
c     elimination of any fictive fugacity-fixing phases left over
c     from a previous run.
c
      call setffg(csts,iindx1,iffg,ipndx1,jpflag,jsflag,kbt,
     $ kdim,kmax,km1,kmt,kx1,kxt,losp,moffg,mtb,mtbaq,nbaspd,nbt,
     $ nbtmax,ncmpr,nffg,nffgmx,noutpt,nphasx,npt,nptmax,nstmax,
     $ nsts,nstsmx,nstsr,nttyo,qloffg,uffg,uspec,uzvec1,zvclg1,zvec1)
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
c
c     Calling sequence substitutions:
c       noutpt for nf
c
        ilevel = iopr(2)
        call echolk(axlks,cdrs,ilevel,jsflag,narxmx,ndrs,ndrsmx,
     $  ndrsr,noutpt,nst,ntprmx,nstmax,press,tempc,uspec,xlks)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nrct .le. 0) go to 350
c
c     Process the input reactant data.
c
c     Get local indices of reactants, check their names, and set up
c     their molecular weights and molar volumes.
c
      call rsetup(atwt,cbsr,cesr,iern1,ietmax,iindx1,iktmax,
     $ jcode,jern1,jern2,jetmax,jgext,kbt,kmax,mwtges,mwtrc,mwtsp,
     $ nbaspd,nbt,nbtmax,nbt1mx,ncmpr,nct,nctmax,nertmx,netmax,ngext,
     $ noutpt,nptmax,nrct,nrctmx,nrndex,nsrtmx,nstmax,nsts,nstsmx,
     $ nstsr,nttyo,nxridx,nxrtmx,rxbar,ureac,uspec,vosp0,vreac,xgers)
c
c     Calculate surface area parameters.
c
      do nrc = 1,nrct
        gx = morr(nrc)*mwtrc(nrc)
        if (nsk(nrc).eq.0 .or. nsk(nrc).eq.2) then
c
c         Calculate the surface area (cm2) from the specific surface
c         area (cm2/g).
c
          ssfcar(nrc) = 0.
          if (gx .gt. 0.) ssfcar(nrc) = sfcar(nrc)/gx
        elseif (nsk(nrc) .eq. 1) then
c
c         Calculate the specific surface area (cm2/g) from the
c         surface area (cm2).
c
          sfcar(nrc) = ssfcar(nrc)*gx
        else
          j2 = ilnobl(ureac(nrc))
          write (noutpt,3110) nsk(nrc),ureac(nrc)(1:j2)
          write (nttyo,3110) nsk(nrc),ureac(nrc)(1:j2)
 3110     format(/' * Error - (EQ6/eq6) The surface area code nsk',
     $    /7x,'has an unrecognized value of ',i3,' for reactant',
     $    /7x,a,'.')
          nerr = nerr + 1
        endif
      enddo
c
c     Set backward/precipitation kinetics flag array, npchk.
c
      do nrc = 1,nrct
        if (nrk(2,nrc) .ne. 0) then
          if (jcode(nrc).eq.0 .or. jcode(nrc).eq.1) then
            do np = 1,npt
            if (ureac(nrc) .eq. uphase(np)) then
              npchk(np) = 1
              go to 300
            endif
            enddo
          endif
        endif
  300   continue
      enddo
c
c     Check names for activity terms in rate equations.
c
      do nrc = 1,nrct
c
        if (nrk(1,nrc) .eq. 2) then
          do i = 1,imech(1,nrc)
            if (ndact(i,1,nrc) .gt. 0) then
              do n = 1,ndact(i,1,nrc)
                unamsp = udac(n,i,1,nrc)
                do ns = 1,nst
                  if (unamsp(1:24) .eq. uspec(ns)(1:24)) go to 310
                enddo
                j3 = ilnobl(unamsp)
                j2 = ilnobl(ureac(nrc))
                write (noutpt,1500) unamsp(1:j3),ureac(nrc)(1:j2)
                write (nttyo,1500) unamsp(1:j3),ureac(nrc)(1:j2)
 1500           format(/" * Error - (EQ6/eq6) Can't match the species",
     $          /7x,'name "',a,'" appearing in the rate law for',
     $          /7x,'with any species present in the model system.')
                nerr = nerr + 1
  310           ndac(n,i,1,nrc) = ns
              enddo
            endif
          enddo
        endif
c
        if (nrk(2,nrc) .eq. 2) then
          do i = 1,imech(2,nrc)
            if (ndact(i,2,nrc) .gt. 0) then
              do n = 1,ndact(i,2,nrc)
                unamsp = udac(n,i,2,nrc)
                do ns = 1,nst
                  if (unamsp(1:24) .eq. uspec(ns)(1:24)) go to 320
                enddo
                j3 = ilnobl(unamsp)
                j2 = ilnobl(ureac(nrc))
                write (noutpt,1500) unamsp(1:j3),ureac(nrc)(1:j2)
                write (noutpt,1500) unamsp(1:j3),ureac(nrc)(1:j2)
                nerr = nerr + 1
  320           ndac(n,i,2,nrc) = ns
              enddo
            endif
          enddo
        endif
c
      enddo
c
      if (nerr .gt. 0) stop
c
  350 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize non-aqueous phases in the physically removed system
c     (PRS).
c
      if (nprpti .gt. 0) then
        nerr = 0
        npi = 0
        unamph = ' '
        do nsi = 1,nprsti
          unamsp(1:24) = uprspi(nsi)(1:24)
          if (unamph(1:24) .ne. uprspi(nsi)(25:48)) then
c
c           Have found another phase.
c
            npi = npi + 1
            unamph(1:24) = uprphi(npi)(1:24)
            do np = 1,npt
              if (unamph(1:24) .eq. uphase(np)(1:24)) go to 370
            enddo
c
            j2 = ilnobl(unamph)
            write (noutpt,1530) unamph(1:j2)
            write (nttyo,1530) unamph(1:j2)
 1530       format(/' * Error - (EQ6/eq6) The phase ',a,' is listed',
     $      /7x,'as present in the physically removed system (PRS)',
     $      ' on the input file,',/7x,"but it isn't in the data file",
     $      ' set (after compression).',/7x,'You may be using the',
     $      ' wrong data file.')
            nerr = nerr + 1
            go to 380
c
  370       mprph(np) = mprphi(npi)
            nr1 = ncmpr(1,np)
            nr2 = ncmpr(2,np)
            do ns = nr1,nr2
              if (unamsp(1:24) .eq. uspec(ns)(1:24)) then
                mprsp(ns) = mprspi(nsi)
                go to 380
              endif
            enddo
c
            j2 = ilnobl(unamph)
            j3 = ilnobl(unamsp)
            write (noutpt,1540) unamsp(1:j3),unamph(1:j2)
            write (nttyo,1540) unamsp(1:j3),unamph(1:j2)
 1540       format(/' * Error - (EQ6/eq6) The species ',a,' of phase',
     $      /7x,a,' is listed as present in the physically removed',
     $      /7x,'system (PRS) on the input file, but this species',
     $      " isn't in the",/7x,'data file set (after compression).',
     $      ' You may be using the wrong',/7x,'data file.')
            nerr = nerr + 1
c
          else
c
c           Have found another species belonging to the same phase.
c
            do ns = nr1,nr2
              if (unamsp(1:24) .eq. uspec(ns)(1:24)) then
                mprsp(ns) = mprspi(nsi)
                go to 380
              endif
            enddo
c
            j2 = ilnobl(unamph)
            j3 = ilnobl(unamsp)
            write (noutpt,1540) unamsp(1:j3),unamph(1:j2)
            write (nttyo,1540) unamsp(1:j3),unamph(1:j2)
            nerr = nerr + 1
          endif
  380     continue
        enddo
c
        if (nerr .gt. 0) stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(9) .gt. 0) then
c
c       Clear the PRS solids read from the input file.
c
        write (noutpt,1570)
        write (nttyo,1570)
 1570   format(/' * Note (EQ6/eq6) Clearing physically removed',
     $  ' system (PRS) solids',/7x,'read from the input file.')
        call initaz(mprsp,nstmax)
        call initaz(mprph,nptmax)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write an echo of the input problem, including defaults, on
c     the output file.
c
      call echoz(axlks,awmaxi,awmini,azero,cbsr,cdac,cdrs,
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
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Roll over the main title to the secondary title.
c
      ntitl2 = ntitl1
      do n = 1,ntitl1
        utitl2(n) = utitl1(n)
      enddo
c
      if (utitl1(1)(7:37) .ne. 'This main title is a carry-over') then
        if (ntitl1 .lt. ntitmx) then
          do n = ntitl1,1,-1
            utitl1(n + 1) = utitl1(n)
          enddo
          utitl1(1) =
     $    'Note: This main title is a carry-over from a previous run.'
        endif
      endif
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
c     Trace the reaction path.
c
      call path(aadh,aadhh,aadhv,aaphi,abdh,abdhh,abdhv,
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
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Descramble the tabx file onto the scratch tab file, then copy
c     the scratch tab file onto the tab file.
c
      if (iopt(18) .ge. 0) then
        if (qtatxt) then
c
c         The TAB file is an ordinary text file.
c
c         Calling sequence substitutions:
c           ntabx for nf1
c           ntabs for nf2
c
          call dscrax(ntabx,ntabs,nllnmx,ulinex)
c
c         Calling sequence substitutions:
c           ntabs for nf1
c           ntab for nf2
c
          call fcopyx(ntabs,ntab,nllnmx,ulinex)
c
        else
c
c         The TAB file is a .csv file.
c
          call dscramc(ntabx,ntabs,nllnmx,ulinex)
c
c         Calling sequence substitutions:
c           ntabs for nf1
c           ntab for nf2
c
          call fcopyx(ntabs,ntab,nllnmx,ulinex)
        endif
      endif
c
      go to 20
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
