      program eqpt
c
c     This is the main program of the EQPT code. Configuration
c     identification, the copyright statement, legal disclaimers,
c     and similar statements are contained in EQPT/aaaeqt.f, the
c     lead-off subroutine in the EQPT source code. A short description
c     of this program is also contained in that subroutine.
c
c-----------------------------------------------------------------------
c
      implicit none


c-----------------------------------------------------------------------
c     File path parameters
c-----------------------------------------------------------------------
      integer :: numargs
      character(1024) :: temppath
      character(:), allocatable :: data0path
      integer :: pathindices(2)
      character(:), allocatable :: basename
      character(:), allocatable :: ofile
      character(:), allocatable :: d1file
      character(:), allocatable :: d1ffile
      character(:), allocatable :: sfile
      character(:), allocatable :: d0sfile
c
c-----------------------------------------------------------------------
c
c     Dimensioning parameters:
c
c       iapxpa = maximum number of interaction coefficients for
c                  computing activity coefficients in solid solutions
c       ibpxpa = maximum number of site-mixing parameters for
c                  computing activity coefficients in solid solutions
c
      integer iapxpa,ibpxpa,iktpar
c
      parameter(iapxpa = 20,ibpxpa = 20)
c
c-----------------------------------------------------------------------
c
c     Array allocation size variables used in EQPT.
c
c       ipch_asv = The order for pressure corrections to enthalpy
c                    functions
c       ipcv_asv = The order for pressure corrections to volume
c                    functions; the maximum order for pressure
c                    corrections to log K and other Gibbs-energy-based
c                    functions is one greater than this
c       nap_asv  = The maximum number of distinct sets of Pitzer alpha
c                    parameters
c       narx_asv = The maximum number of coefficients per temperature
c                    range
c       nat_asv  = The number of aqueous species on the data file
c       nazt_asv = the number of aqeuous species on the data file for
c                    which hard core diameters are specified
c       nbt_asv  = The number of basis species on the data file
c       nct_asv  = The number of chemical elements on the data file
c       ngt_asv  = The number of gas species on the data file
c       nlt_asv  = The number of pure liquids on the data file
c       nmt_asv  = The number of pure minerals on the data file
c       npx2_asv = The number of pairs of species not of the same
c                    charge sign for which Pitzer parameters are
c                    defined; typically, one of the pair is a cation
c                    and the other is an anion, but one or both
c                    species may also be electrically neutral
c       npx3_asv = The number of triplets of species corresponding to
c                    aqueous electrolyte mixtures for which Pitzer
c                    parameters are defined; generally, no more than
c                    two of these may have an electrical charge number
c                    that is postive, negative, or zero
c       ntid_asv = The number of lines in the data file title
c       ntpr_asv = The number of temperature ranges
c       ipbt_asv = The greatest index for a Pitzer beta(n) coefficient
c                     (the first such index is zero); also the number
c                     number of Pitzer alpha coefficients in a set
c       jpfc_asv = The number of coefficients in a function for
c                    representing a Pitzer coefficient
c
c       ndb_asv  = maximum number of distinct points on the temperature
c                    grid; ndb_asv = ntpr_asv*(narx_asv - 1) + 1
c
      integer ipbt_asv,ipch_asv,ipcv_asv,jpfc_asv,narx_asv,nat_asv,
     $ nazt_asv,nbt_asv,nct_asv,ngt_asv,nlt_asv,nmt_asv,npx2_asv,
     $ npx3_asv,ntid_asv,ntpr_asv
c
      integer ndb_asv
c
c-----------------------------------------------------------------------
c
c     Array allocation size variables which are determined in EQPT
c     but only written on the DATA1 file for use by EQ3NR and EQ6.
c
c       ikt_asv  = The maximum number of end-member component species
c                    of any solid solution on the data filen
c       nlat_asv = The number of members in the set of Pitzer lambda
c                    coefficients
c       nmut_asv = The number of members in the set of Pitzer mu
c                    coefficients
c       npt_asv  = The number of phases of all types on the data file
c       nst_asv  = The number of species of all types on the data file
c       nxt_asv  = The number of solid-solution phases on the data file
c
      integer ikt_asv,nap_asv,nlat_asv,nmut_asv,npt_asv,nst_asv,nxt_asv
c
c-----------------------------------------------------------------------
c
c     Global variable declarations.
c
c     Dimensioning variables.
c
      integer iapxmx,ibpxmx,iktmax,ipbtmx,ipchmx,ipcvmx,jpfcmx,narxmx,
     $ natmax,naztmx,nbtmax,nctmax,ngtmax,nltmax,nmtmax,npx2mx,npx3mx,
     $ ntidmx,ntprmx,nxtmax
      integer nbtmx1,nbtmx2,ndbmax
c
c     File unit numbers.
c
      integer ndata0,ndata1,ndat0s,ndat1f,noutpt,nslist,nttyo
c
c     Other variables.
c
      integer, dimension(:), allocatable :: insgf,ipivot,issot,nacdpr,
     $ narxt,nentei,nentri
c
      integer ier,ipch,ipcv,irang,nazt,nat,nbt,nct,ncvaz,ndbptg,ndbptl,
     $ nerr,ngt,nlt,nmt,nmodwr,nsb,nthdt,ntitld,ntprt,nwarn,nxt
c
      logical, dimension(:), allocatable :: qpdaz
c
      logical qelect
c
      character(len=80), dimension(:), allocatable :: utitld
      character(len=24), dimension(:), allocatable :: uaqsp,uazp,
     $ udrsi,ugassp,uliqsp,uminsp,uspec,ussoph
      character(len=24), dimension(:,:), allocatable :: ussosp
      character(len=16), dimension(:), allocatable :: udbval
      character(len=8), dimension(:), allocatable :: uelem,uessi
c
      character(len=16) udbfmt
      character(len=8) uplatc,uplatm,uveeqt,usteqt,uveelu,ustelu
      character(len=8) uakey,uethfl
c
      real(8), dimension(:), allocatable :: atwt,azero,cdrsi,cessi,
     $ mtotr,zaqsp,zchar
      real(8), dimension(:,:), allocatable :: cdrs,cess
c
      real(8) cco2(5)
c
      real(8), dimension(:,:,:), allocatable :: xhfs,xlks,xvfs
      real(8), dimension(:,:,:,:), allocatable :: dhfs,dvfs
c
      real(8), dimension(:,:), allocatable :: adh,adhh,adhv,aphi,bdh,
     $ bdhh,bdhv,bdot,bdoth,bdotv,prehw,presg,xhfe,xlke,xvfe
      real(8), dimension(:,:,:), allocatable :: dadhh,dadhv,dbdhh,
     $ dbdhv,dbdth,dbdtv,dhfe,dvfe
c
      real(8), dimension(:), allocatable :: tmpcmx,xdbval
      real(8), dimension(:,:), allocatable :: apr,avgrid,tempc,tempcs
c
      real(8), dimension(:), allocatable :: tvec,tvecs
      real(8), dimension(:), allocatable :: cof,xvec,yvec
      real(8), dimension(:,:), allocatable :: aamatr,gmmatr
c
      real(8) pcvaz
c
      real(8) eps100,smp100,tdamax,tdamin,tvecmx
c
c-----------------------------------------------------------------------
c
c     Pitzer interaction parameter arrays.
c
      integer naapr,ncapr,nccpr,nnapr,nncpr,nnnpr,nn2pr,
     $ naactr,na2ctr,nccatr,nc2atr,nncatr,nn2ntr,nn3tr
c
      integer, dimension(:), allocatable :: in2pr,in3tr
      integer, dimension(:,:), allocatable :: iaapr,icapr,iccpr,inapr,
     $ incpr,innpr
      integer, dimension(:,:), allocatable :: iaactr,ia2ctr,iccatr,
     $ ic2atr,incatr,in2ntr
c
      integer jassan,jassca,jassne,npx2t,npx3t
c
      integer ncvaa,ncvca,ncvcc,ncvna,ncvnc,ncvnn,ncvn2,
     $ ncvaac,ncvcca,ncvn2n,ncvnca
c
      integer npxca,npxth,npxni,npxnn,npxn2,npxpsi,npxzet,npxn2n
c
      logical, dimension(:), allocatable :: qpdaa,qpdca,qpdcc,qpdna,
     $ qpdnc,qpdnn,qpdn2
c
      logical, dimension(:), allocatable :: qpdaac,qpdcca,qpdnca,qpdn2n
c
      character(len=24), dimension(:,:), allocatable :: upair,uthdtr,
     $ utripl
c
      real(8), dimension(:,:), allocatable :: alpha,alphca,acphi,
     $ amua2c,amuaac,amucca,amuc2a,amunca,amun2n,amun3
      real(8), dimension(:,:,:), allocatable :: abeta,alamaa,alamca,
     $ alamcc,alamna,alamnc,alamnn,alamn2
c
      real(8), dimension(:,:), allocatable :: atheta,apsi
c
      real(8) pcvaa,pcvca,pcvcc,pcvna,pcvnc,pcvnn,pcvn2
c
      real(8) pcvaac,pcvcca,pcvn2n,pcvnca
c
c-----------------------------------------------------------------------
c
cXX
      real*8 apx(iapxpa),bpx(ibpxpa)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,iexec0,itgenf,j,jexec0,jpdblo,jptffl,j2,j3,j4,k,n,
     $ natm1,nch,nco,nmax,npx2,npx2r1,npx2r2,nrecl,nthd,ntpr
c
      integer ilnobl
c
      logical q500fl
c
      character(len=16) uacfst
      character(len=11) utime0,utime1
      character(len=8) ux8,ux8a,ux8b
      character(len=9) udate0,udate1
      character(len=80) uline
c
      real(8) tcpu,texec0,trun,tuser
c
c-----------------------------------------------------------------------
c
      data nrecl /0/
c
c-----------------------------------------------------------------------
c
c     BEGIN_MACHINE_DEPENDENT_CODE
c
c       On some systems, a BLOCK DATA subroutine must be declared in
c       an EXTERNAL statement to assure proper loading. On some other
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
        external bkdeqp
c
c     END_MACHINE_DEPENDENT_CODE
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
      ndata1 = 0
      ndat1f = 0
      nslist = 0
      ndata0 = 0
      ndat0s = 0
c
c     Open all output files.
c
      numargs = COMMAND_ARGUMENT_COUNT()
      if (numargs.eq.1) then
          call GET_COMMAND_ARGUMENT(1,temppath)
          data0path = TRIM(temppath)
      else
          write (0, *) 'usage: eqpt <data0>'
          stop
      end if

      call getbasename(data0path,pathindices)
      basename = data0path(pathindices(1):pathindices(2))

      call openin(noutpt,nttyo,data0path,'formatted',ndata0)

      ofile = basename // '.po'
      d1file = basename // '.d1'
      d1ffile = basename // '.d1f'
      sfile = basename // '.s'
      d0sfile = basename // '.d0s'

      call openou(noutpt,nttyo,ofile,'formatted',nrecl,noutpt)
      call openou(noutpt,nttyo,d1file,'unformatted',nrecl,ndata1)
      call openou(noutpt,nttyo,d1ffile,'formatted',nrecl,ndat1f)
      call openou(noutpt,nttyo,sfile,'formatted',nrecl,nslist)
      call openou(noutpt,nttyo,d0sfile,'formatted',nrecl,ndat0s)
c
c     Make a copy of the DATA0 file, stripped of comments.
c
      call stripl(ndata0,ndat0s)
      close(ndata0)
      ndata0 = 0
      rewind ndat0s
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get configuration identification data.
c
      call aaaeqt(usteqt,uveeqt)
      call aaaelu(ustelu,uveelu)
      call platfd(uplatc,uplatm)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write configuration identification data, the copyright statement,
c     and any remaining statements or disclaimers.
c
      i = index(uveeqt,' ') - 1
      j = index(uveeqt,'.') - 1
      k = index(uplatc,' ') - 1
      if (i .le. 0) i = 8
      if (j .le. 0) j = 8
      if (k .le. 0) k = 8
      write (nttyo,1000) uveeqt(1:i),uveeqt(1:j),uveeqt(1:i),
     $ uplatc(1:k)
      write (noutpt,1000) uveeqt(1:i),uveeqt(1:j),uveeqt(1:i),
     $ uplatc(1:k)
 1000 format(/' EQ3/6, Version ',a,' (EQ3/6-V',a,'-REL-V',a,'-',a,')')
c
      i = j
      j = index(usteqt,' ') - 1
      k = index(uplatm,' ') - 1
      if (j .le. 0) j = 8
      if (k .le. 0) k = 8
      write (nttyo,1010) uveeqt(1:i),usteqt(1:j),uplatm(1:k)
      write (noutpt,1010) uveeqt(1:i),usteqt(1:j),uplatm(1:k)
 1010 format(' EQPT Data File Preprocessor Code (EQ3/6-V',a,
     $ '-EQPT-EXE-',a,'-',a,')')
c
      write (noutpt,1020)
      write (nttyo,1020)
 1020 format(' Supported by the following EQ3/6 libraries:')
c
      i = index(uveelu,'.') - 1
      j = index(ustelu,' ') - 1
      if (i .le. 0) i = 8
      if (j .le. 0) j = 8
      write (nttyo,1030) uveelu(1:i),ustelu(1:j),uplatm(1:k)
      write (noutpt,1030) uveelu(1:i),ustelu(1:j),uplatm(1:k)
 1030 format('   EQLIBU (EQ3/6-V',a,'-EQLIBU-LIB-',a,'-',a,')',/)
c
      write (nttyo,1040)
      write (noutpt,1040)
 1040 format(' Copyright (c) 1987, 1990-1993, 1995, 1997, 2002 The',
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
      write(noutpt,1070) utime0,udate0(1:j2)
      write(nttyo,1070) utime0,udate0(1:j2)
 1070 format(' Run',2x,a8,2x,a,//)
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
      iapxmx = iapxpa
      ibpxmx = ibpxpa
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write header on the slist file.
c
      write (nslist,1090)
 1090 format(' EQPT Species List (SLIST) File:',//)
c
c     Initialize the cumulative error counter.
c
      nerr = 0
c
c     Initialize the cumulative warning counter.
c
      nwarn = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check the first line of the data file to ensure that the
c     required header is present.
c
      call hdrchk(ndat0s,noutpt,nttyo)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the aqueous species activity coefficient model type
c     (Extended Debye-Huckel, Pitzer, etc.) associated with this
c     data file. This subroutine scans the data file and counts
c     keywords to make the determination.
c
      call gakey(ndat0s,noutpt,nttyo,uakey)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Scan the data file to determine the necessary array dimensions
c     and any corresponding structural information that will be
c     required before "reading" the data file.
c
c     First scan the data file title for embedded data.
c
      call ggridp(ipch_asv,ipcv_asv,itgenf,jpdblo,jpfc_asv,
     $ jptffl,narx_asv,ndat0s,ndb_asv,noutpt,ntid_asv,ntpr_asv,
     $ nttyo,q500fl,uakey)
c
c     Set the number of parameters in a Pitzer alpha set.
c
      ipbt_asv = 1
      if (uakey(1:8) .eq. 'Pitzer  ') ipbt_asv = 2
c
c     Rescan the data file to get other array dimensions.
c     Also get the actual numbers of basis species and chemical
c     elements on the data file.
c
      call gnenb(ipbt_asv,ikt_asv,jpdblo,jpfc_asv,nap_asv,
     $ nat_asv,nazt_asv,nbt_asv,nct_asv,ndat0s,ngt_asv,nlt_asv,
     4 nmt_asv,noutpt,npt_asv,npx2_asv,npx3_asv,nsb,nst_asv,nttyo,
     $ nxt_asv,uakey)
c
      nbt = nbt_asv
      nct = nct_asv
c
      if (nsb .ne. (nct + 1)) then
        write (ux8a,'(i5)') nsb
        call lejust(ux8a)
        j2 = ilnobl(ux8a)
        write (ux8b,'(i5)') nct
        call lejust(ux8b)
        j3 = ilnobl(ux8b)
        write(noutpt,1100) ux8a(1:j2),ux8b(1:j3)
        write(nttyo,1100) ux8a(1:j2),ux8b(1:j3)
 1100   format(/' * Error - (EQPT/eqpt) The number of strict basis',
     $  ' species',/7x,'must be the number of chemical elements plus',
     $  ' one. The number of',/7x,'strict basis species is ',a,',',
     $  ' but the number of chemical elements',/7x,'is ',a,'.')
        stop
      endif
c
      ipbt_asv = 1
      if (uakey(1:8) .eq. 'Pitzer  ') then
c
        ipbt_asv = 2
c
c       The following array allocation variables are only used by
c       EQ3NR and EQ6.
c
        nlat_asv = npx2_asv + npx3_asv
        nmut_asv = npx2_asv + 2*npx3_asv
        nlat_asv = max(1,nlat_asv)
        nmut_asv = max(1,nmut_asv)
      endif
c
      iktmax = ikt_asv
      ipbtmx = ipbt_asv
      ipchmx = ipch_asv
      ipcvmx = ipcv_asv
      jpfcmx = jpfc_asv
      narxmx = narx_asv
      natmax = nat_asv
      naztmx = nazt_asv
      nbtmax = nbt_asv
      nctmax = nct_asv
      ndbmax = ndb_asv
      ngtmax = ngt_asv
      nltmax = nlt_asv
      nmtmax = nmt_asv
      npx2mx = npx2_asv
      npx3mx = npx3_asv
      ntidmx = ntid_asv
      ntprmx = ntpr_asv
      nxtmax = nxt_asv
c
      nbtmx1 = nbtmax + 1
      nbtmx2 = nbtmax + 2
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Allocate the necessary arrays.
c
      ALLOCATE(ipivot(narx_asv))
      ALLOCATE(issot(nxt_asv))
      ALLOCATE(nacdpr(ntpr_asv))
      ALLOCATE(narxt(ntpr_asv))
      ALLOCATE(nentei(nct_asv))
      ALLOCATE(nentri(nbt_asv + 1))
c
      ALLOCATE(qpdaz(nat_asv))
c
      ALLOCATE(uaqsp(nat_asv))
      ALLOCATE(zaqsp(nat_asv))
c
      ALLOCATE(ugassp(ngt_asv))
      ALLOCATE(uliqsp(nlt_asv))
      ALLOCATE(uminsp(nmt_asv))
      ALLOCATE(ussoph(nxt_asv))
c
      ALLOCATE(utitld(ntid_asv))
c
      ALLOCATE(uspec(nbt_asv + 1))
      ALLOCATE(udrsi(nbt_asv + 1))
c
      ALLOCATE(ussosp(ikt_asv,nxt_asv))
c
      ALLOCATE(uelem(nct_asv))
      ALLOCATE(uessi(nct_asv))
c
      ALLOCATE(atwt(nct_asv))
      ALLOCATE(cdrsi(nbt_asv + 1))
      ALLOCATE(cessi(nct_asv))
      ALLOCATE(zchar(nbt_asv + 1))
      ALLOCATE(mtotr(nct_asv))
c
      ALLOCATE(cdrs(nbt_asv + 2,nbt_asv + 1))
      ALLOCATE(cess(nct_asv,nbt_asv + 1))
c
      ALLOCATE(dhfs(narx_asv,ntpr_asv,ipch_asv,nbt_asv + 1))
      ALLOCATE(dvfs(narx_asv,ntpr_asv,ipcv_asv,nbt_asv + 1))
      ALLOCATE(xhfs(narx_asv,ntpr_asv,nbt_asv + 1))
      ALLOCATE(xlks(narx_asv,ntpr_asv,nbt_asv + 1))
      ALLOCATE(xvfs(narx_asv,ntpr_asv,nbt_asv + 1))
c
      ALLOCATE(adh(narx_asv,ntpr_asv))
      ALLOCATE(adhh(narx_asv,ntpr_asv))
      ALLOCATE(adhv(narx_asv,ntpr_asv))
      ALLOCATE(aphi(narx_asv,ntpr_asv))
      ALLOCATE(bdh(narx_asv,ntpr_asv))
      ALLOCATE(bdhh(narx_asv,ntpr_asv))
      ALLOCATE(bdhv(narx_asv,ntpr_asv))
      ALLOCATE(bdot(narx_asv,ntpr_asv))
      ALLOCATE(bdoth(narx_asv,ntpr_asv))
      ALLOCATE(bdotv(narx_asv,ntpr_asv))
      ALLOCATE(prehw(narx_asv,ntpr_asv))
      ALLOCATE(presg(narx_asv,ntpr_asv))
      ALLOCATE(xhfe(narx_asv,ntpr_asv))
      ALLOCATE(xlke(narx_asv,ntpr_asv))
      ALLOCATE(xvfe(narx_asv,ntpr_asv))
c
      ALLOCATE(dadhh(narx_asv,ntpr_asv,ipch_asv))
      ALLOCATE(dbdhh(narx_asv,ntpr_asv,ipch_asv))
      ALLOCATE(dbdth(narx_asv,ntpr_asv,ipch_asv))
      ALLOCATE(dhfe(narx_asv,ntpr_asv,ipch_asv))
c
      ALLOCATE(dadhv(narx_asv,ntpr_asv,ipcv_asv))
      ALLOCATE(dbdhv(narx_asv,ntpr_asv,ipcv_asv))
      ALLOCATE(dbdtv(narx_asv,ntpr_asv,ipcv_asv))
      ALLOCATE(dvfe(narx_asv,ntpr_asv,ipcv_asv))
c
      ALLOCATE(apr(narx_asv,ntpr_asv))
      ALLOCATE(avgrid(narx_asv,ntpr_asv))
      ALLOCATE(tempc(narx_asv,ntpr_asv))
      ALLOCATE(tempcs(narx_asv,ntpr_asv))
      ALLOCATE(tmpcmx(ntpr_asv))
      ALLOCATE(udbval(ndb_asv))
      ALLOCATE(xdbval(ndb_asv))
c
      ALLOCATE(tvec(narx_asv))
      ALLOCATE(tvecs(narx_asv))
      ALLOCATE(cof(narx_asv))
      ALLOCATE(xvec(narx_asv))
      ALLOCATE(yvec(narx_asv))
c
      ALLOCATE(aamatr(narx_asv,narx_asv))
      ALLOCATE(gmmatr(narx_asv,narx_asv))
c
      if (uakey(1:8) .eq. 'SEDH    ') then
        ALLOCATE(uazp(nazt_asv))
        ALLOCATE(azero(nazt_asv))
        ALLOCATE(insgf(nazt_asv))
      endif
c
      if (uakey(1:8) .eq. 'Pitzer  ') then
c
        ALLOCATE(alpha(ipbt_asv,npx2_asv))
        ALLOCATE(abeta(jpfc_asv,0:ipbt_asv,npx2_asv))
        ALLOCATE(acphi(jpfc_asv,npx2_asv))
c
        ALLOCATE(upair(2,npx2_asv))
        ALLOCATE(uthdtr(3,npx3_asv))
        ALLOCATE(utripl(3,npx3_asv))
c
        ALLOCATE(atheta(jpfc_asv,npx3_asv))
        ALLOCATE(apsi(jpfc_asv,npx3_asv))
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write some tallies on the OUTPUT file.
c
      write (noutpt,1120) nct,nbt
      write (nslist,1120) nct,nbt
 1120 format(11x,'Number of elements = ',i5,
     $ /11x,'Number of basis species = ',i5,/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the header on the DATA1 and DATA1f files. This begins with
c     the string 'data1' at the top of the files, followed by the
c     keystring for the type of aqueous species activity coefficient
c     model, and the array dimensions necessary to read the rest of
c     the DATA1 file.
c
      call wrhdr(ikt_asv,ipch_asv,ipcv_asv,jpfc_asv,nap_asv,
     $ narx_asv,nat_asv,nbt_asv,nct_asv,ndata1,ndat1f,ngt_asv,
     $ nlat_asv,nlt_asv,nmt_asv,nmut_asv,npt_asv,nst_asv,ntid_asv,
     $ ntpr_asv,nxt_asv,uakey)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the data file title on the various output files. Determine
c     certain embedded options apart from those that set dimensioning
c     parameters.
c
      ntitld = ntid_asv
      call rdwttl(ipch,ipcv,jpdblo,jpfcmx,jptffl,narxt,ndata1,
     $ ndat0s,ndat1f,noutpt,nslist,ntitld,ntidmx,ntprmx,ntprt,nttyo,
     $ uakey,utitld)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the miscellaneous parameters (write them after the chemical
c     elements block). Write the nominal temperature limits and the
c     upper limits of the temperature ranges here, however.
c
      call rdpar(adh,adhh,adhv,aphi,bdh,bdhh,bdhv,bdot,bdoth,
     $ bdotv,cco2,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,dhfe,dvfe,
     $ ipch,ipchmx,ipcv,ipcvmx,itgenf,nacdpr,narxmx,narxt,ndat0s,
     $ ndbmax,ndbptg,ndbptl,nerr,noutpt,ntprmx,ntprt,nttyo,nwarn,
     $ prehw,presg,q500fl,tdamax,tdamin,tempc,uakey,udbfmt,udbval,
     $ xdbval,xhfe,xlke,xvfe)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read and write the element data.
c
      call rdwele(atwt,nch,nco,nct,nctmax,ndata1,ndat0s,ndat1f,
     $ nerr,noutpt,nslist,nttyo,uelem)
c
c     Check the element names for uniqueness.
c
      call neleck(nct,nctmax,nerr,noutpt,nttyo,uelem)
c
      if (nerr .gt. 0) then
        write (ux8,'(i5)') nerr
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1160) ux8(1:j2)
        write (nttyo,1160) ux8(1:j2)
 1160   format(//' ',a,' errors(s) were encountered.',/)
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Scale the grid temperature values for subsequent interpolation.
c
      do ntpr = 1,ntprt
c
        nmax = narxt(ntpr)
        do n = 1,nmax
          tvec(n) = tempc(n,ntpr)
        enddo
c
c       Calling sequence substitutions:
c         tvec for avx
c         tvecmx for avxmax
c         tvecs for avxs
c
        call scalx1(tvec,tvecmx,tvecs,ier,nmax)
        do n = 1,nmax
          tempcs(n,ntpr) = tvecs(n)
        enddo
        tmpcmx(ntpr) = tvecmx
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set the interval parameter (nmodwr) for writing species and
c     phase names to nttyo. This is scaled to the number of chemical
c     elements.
c
      if (nct .gt. 50) then
        nmodwr = 20
      elseif (nct .gt. 30) then
        nmodwr = 10
      elseif (nct .gt. 20) then
        nmodwr = 5
      else
        nmodwr = 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Interpolate and write the miscellaneous parameters.
c
      call wrpar(aamatr,adh,adhh,adhv,aphi,apr,avgrid,bdh,bdhh,
     $ bdhv,bdot,bdoth,bdotv,cco2,cof,dadhh,dadhv,dbdhh,dbdhv,dbdth,
     $ dbdtv,dhfe,dvfe,eps100,gmmatr,ipch,ipchmx,ipcv,ipcvmx,ipivot,
     $ narxmx,narxt,ndata1,ndat1f,noutpt,ntprmx,ntprt,nttyo,presg,
     $ prehw,tdamax,tdamin,tempc,tempcs,tmpcmx,uakey,xhfe,xlke,xvfe,
     $ xvec,yvec)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read and write aqueous species.
c
      call pcraq(aamatr,apr,atwt,avgrid,cdrs,cdrsi,cess,cessi,
     $ cof,dhfe,dhfs,dvfe,dvfs,eps100,gmmatr,ipch,ipchmx,ipcv,ipcvmx,
     $ ipivot,itgenf,mtotr,nacdpr,narxmx,narxt,nat,natmax,nbt,nbtmx1,
     $ nbtmx2,nch,nco,nct,nctmax,ndata1,ndat0s,ndat1f,ndbmax,ndbptg,
     $ ndbptl,nentei,nentri,nerr,nmodwr,noutpt,nsb,nslist,ntprmx,ntprt,
     $ nttyo,nwarn,qelect,q500fl,tempc,tempcs,tmpcmx,uaqsp,udbfmt,
     $ udbval,udrsi,uelem,uessi,uspec,xdbval,xhfe,xhfs,xlke,xlks,xvfe,
     $ xvfs,xvec,yvec,zaqsp,zchar)
c
c     Read and write minerals, liquids, gases.
c
      call pcrsg(aamatr,apr,atwt,avgrid,cdrs,cdrsi,cess,cessi,
     $ cof,dhfe,dhfs,dvfe,dvfs,eps100,gmmatr,ipch,ipchmx,ipcv,ipcvmx,
     $ ipivot,itgenf,mtotr,nacdpr,narxmx,narxt,nat,natmax,nbt,nbtmx1,
     $ nbtmx2,nch,nco,nct,nctmax,ndata1,ndat0s,ndat1f,ndbmax,ndbptg,
     $ ndbptl,nentei,nentri,nerr,ngt,ngtmax,nlt,nltmax,nmodwr,nmt,
     $ nmtmax,noutpt,nsb,nslist,ntprmx,ntprt,nttyo,nwarn,qelect,
     $ q500fl,tempc,tempcs,tmpcmx,udbfmt,udbval,udrsi,uelem,uessi,
     $ ugassp,uliqsp,uminsp,uspec,xdbval,xhfe,xhfs,xlke,xlks,xvfe,
     $ xvfs,xvec,yvec,zchar)
c
c     Read and write solid solutions.
c
      call pcrss(apx,bpx,iapxmx,ibpxmx,iktmax,issot,nbtmx1,
     $ ndata1,ndat0s,ndat1f,nerr,nmodwr,nmt,nmtmax,noutpt,nslist,
     $ nttyo,nxt,nxtmax,uminsp,ussosp,ussoph)
c
c     Make sure that all aqueous species names are unique.
c
      call naqsck(nat,natmax,nerr,noutpt,nttyo,uaqsp)
c
c     Make sure that all pure mineral names are unique.
c
      call nminck(nerr,nmt,nmtmax,noutpt,nttyo,uminsp)
c
c     Make sure that all gas species names are unique.
c
      call ngasck(nerr,ngt,ngtmax,noutpt,nttyo,ugassp)
c
c     Make sure that all solid solution names are unique.
c
      call nssock(nerr,noutpt,nttyo,nxt,nxtmax,ussoph)
c
c     Make sure that all solid solution end-member names are unique
c     within each solid solution.
c
      call nxspck(iktmax,issot,nerr,noutpt,nttyo,nxt,nxtmax,
     $ ussoph,ussosp)
c
c     Validate the end-member names for the solid solutions. Each
c     end-member name must match that of a pure mineral.
c
      call vxspck(iktmax,issot,nerr,nmt,nmtmax,noutpt,nttyo,
     $ nxt,nxtmax,uminsp,ussoph,ussosp)
c
      if (nerr .gt. 0) then
        write (ux8,'(i5)') nerr
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1160) ux8(1:j2)
        write (nttyo,1160) ux8(1:j2)
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Rewind and position to read aqueous species activity coefficient
c     parameters.
c
      rewind(ndat0s)
c
      if (uakey(1:8) .eq. 'SEDH    ') then
        uacfst = 'bdot parameters '
      elseif (uakey(1:8) .eq. 'Pitzer  ') then
        if (jpdblo .eq. -1) then
          uacfst = 'single-salt para'
        else
          uacfst = 'ca combinations:'
        endif
      else
        j2 = ilnobl(uakey)
        write(noutpt,1200) uakey(1:j2)
        write(nttyo,1200) uakey(1:j2)
 1200   format(/' * Error - (EQPT/eqpt) Unrecognized data file',
     $  ' type = "',a,'".',/7x,'Allowed types are "SEDH" (Simple',
     $  ' Extended Debye-Huckel) and',/7x,'"Pitzer".')
        stop
      endif
c
  270 read(ndat0s,1210,end=280,err=280) uline
 1210 format(a)
      if (uline(1:16) .ne. uacfst(1:16)) go to 270
      go to 290
c
  280 j3 = ilnobl(uacfst)
      write(noutpt,1220) uacfst(1:j3)
      write(nttyo,1220) uacfst(1:j3)
 1220 format(/' * Error - (EQPT/eqpt) End-of-file hit or other read',
     $ /7x,'error occurred while searching for the block header',
     $ /7x,'beginning with "',a,'".')
      stop
c
  290 backspace(ndat0s)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Zero arrays for activity coefficient parameters.
c
      if (uakey(1:8) .eq. 'SEDH    ') then
c
c       Zero the azero array.
c
        call initaz(azero,naztmx)
c
c       Zero the insgf array.
c
        call initiz(insgf,naztmx)
      endif
c
      if (uakey(1:8) .eq. 'Pitzer  ') then
c
c       Zero the alpha array.
c
        nmax = ipbtmx*npx2mx
        call initaz(alpha,nmax)
c
c       Zero the abeta array.
c
        nmax = jpfcmx*(ipbtmx + 1)*npx2mx
        call initaz(abeta,nmax)
c
c       Zero the acphi array.
c
        nmax = jpfcmx*npx2mx
        call initaz(acphi,nmax)
c
c       Zero the atheta and apsi arrays.
c
        nmax = jpfcmx*npx3mx
        call initaz(atheta,nmax)
        call initaz(apsi,nmax)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the SEDH or Pitzer section.
c
      if (uakey(1:8) .eq. 'SEDH    ') then
c
c       Read hard core diameters and related parameters used
c       in the B-dot equation.
c
        call rdazp(azero,insgf,nazt,naztmx,ndat0s,nerr,noutpt,
     $  nttyo,uazp)
      endif
c
      if (uakey(1:8) .eq. 'Pitzer  ') then
c
c       Read Pitzer parameters.
c
        if (jpdblo .eq. -1) then
c
c         Read the data according to the old Pitzer data block
c         organization. There is one superblock for data
c         pertaining to pure aqueous electrolytes, another
c         for mixtures of such with a common ion. Parameters
c         involving neutral aqueous species are handled (awkwardly,
c         but with appropriate checks) within this structure as
c         noted below.
c
c         Read the beta(n) and Cphi parameters.
c           - normally for ca (cation-anion) combinations.
c           - lamda(0) values may be read in place of beta(0) values
c             for nn, nn', nc, and na combinations.
c           - mu values may be read in place of Cphi values for
c             nnn combinations.
c
          call rdpz2(abeta,alpha,acphi,ipbtmx,jpfcmx,nat,natmax,
     $    ndat0s,nerr,noutpt,npx2mx,npx2t,nttyo,nwarn,uaqsp,upair,
     $    zaqsp)
c
c         Read the S-theta and psi parameters.
c           - normally for cc'a and aa'c combinations.
c           - zeta values may be read in place of psi values for
c             nca combinations.
c           - mu values may be read in place of psi values for
c             nnn' and n'n'n combinations.
c
          call rdpz3(apsi,atheta,jpfcmx,nat,natmax,ndat0s,nerr,
     $    noutpt,npx3mx,npx3t,nthdt,nttyo,nwarn,uaqsp,uethfl,uthdtr,
     $    utripl,zaqsp)
        else
c
c         Read the data according to the new Pitzer data block
c         organization. There is one superblock for each allowed
c         species pair or triplet type. Note that for the "new"
c         organization, the E-theta flag is hardwired to 'on'.
c         That is because basically all modern implementations
c         of Pitzer's equations use the theoretically-based
c         higher-order electrostatic term formalism (which is
c         turned on by uethfl = 'on').
c
          uethfl = 'on'
          npx2t = 0
          npx3t = 0
c
c         Read ca (cation-anion) data.
c
          call rdpca(abeta,alpha,acphi,ipbtmx,jpfcmx,nat,natmax,
     $    ndat0s,nerr,noutpt,npxca,npx2mx,npx2t,nttyo,nwarn,uaqsp,
     $    upair,zaqsp)
c
c         Read cc' (cation-different cation) and aa' (anion-
c         different anion) data.
c
          call rdpth(abeta,ipbtmx,jpfcmx,nat,natmax,ndat0s,nerr,
     $    noutpt,npxth,npx2mx,npx2t,nttyo,nwarn,uaqsp,upair,zaqsp)
c
c         Read nc (neutral-cation) and na (neutral-anion) data.
c
          call rdpni(abeta,ipbtmx,jpfcmx,nat,natmax,ndat0s,nerr,
     $    noutpt,npxni,npx2mx,npx2t,nttyo,nwarn,uaqsp,upair,zaqsp)
c
c         Read nn (neutral-same neutral) data.
c
          call rdpn2(abeta,acphi,ipbtmx,jpfcmx,nat,natmax,
     $    ndat0s,nerr,noutpt,npxn2,npx2mx,npx2t,nttyo,nwarn,
     $    uaqsp,upair,zaqsp)
c
c         Read nn'(neutral-different neutral) data.
c
          call rdpnn(abeta,ipbtmx,jpfcmx,nat,natmax,ndat0s,nerr,
     $    noutpt,npxnn,npx2mx,npx2t,nttyo,nwarn,uaqsp,upair,zaqsp)
c
c         Read cc'a (cation-different cation-anion) and aa'c
c         (anion-different anion-cation) data.
c
          call rdppsi(apsi,jpfcmx,nat,natmax,ndat0s,nerr,noutpt,
     $    npxpsi,npx3mx,npx3t,nttyo,nwarn,uaqsp,utripl,zaqsp)
c
c         Read nca (neutral-cation--anion) data.
c
          call rdpzet(apsi,jpfcmx,nat,natmax,ndat0s,nerr,noutpt,
     $    npxzet,npx3mx,npx3t,nttyo,nwarn,uaqsp,utripl,zaqsp)
c
c         Read nnn' (neutral-neutral-different neutral) data.
c
          call rdpn2n(apsi,jpfcmx,nat,natmax,ndat0s,nerr,noutpt,
     $    npxn2n,npx3mx,npx3t,nttyo,nwarn,uaqsp,utripl,zaqsp)
c
c         Copy the theta data into arrays used with the old
c         Pitzer data block organization.
c
          nthdt = npxth
          nthd = 0
          npx2r1 = npxca + 1
          npx2r2 = npxca + npxth
          do npx2 = npx2r1,npx2r2
            nthd = nthd + 1
            uthdtr(1,nthd) = upair(1,npx2)
            uthdtr(2,nthd) = upair(2,npx2)
            uthdtr(3,nthd) = '<dummy>'
            do j = 1,jpfcmx
              atheta(j,nthd) = abeta(j,0,npx2)
            enddo
          enddo
c
        endif
      endif
c
      if (nerr .gt. 0) then
        write (ux8,'(i5)') nerr
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1160) ux8(1:j2)
        write (nttyo,1160) ux8(1:j2)
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Test and write the SEDH data section of DATA1.
c
      if (uakey(1:8) .eq. 'SEDH    ') then
c
c       Validate the names of aqueous species used to specify azero
c       ('bdot') data. Each such name should correspond to an aqueous
c       species for which there is a species block on the data file.
c       Write a note for any exceptions.
c
        call vazpck(nat,natmax,nazt,naztmx,noutpt,nttyo,
     $  uaqsp,uazp)
c
c       Test and process all azero and insgf ('bdot') data.
c
        call tpraz(nat,natmax,nazt,naztmx,ncvaz,nerr,noutpt,nttyo,
     $  pcvaz,qpdaz,uaqsp,uazp)
c
c       Write the azero and insgf data summary.
c
        ux8 = ' '
        write (ux8,'(f6.2)') pcvaz
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1230)
 1230   format(//2x,'Aqueous species hard core diameter coverage:')
        natm1 = nat - 1
        write (noutpt,1240) ncvaz,natm1,ux8(1:j2)
 1240   format(/4x,i5,' aqueous species have hard core diameters',
     $  ' specified on the data file',/4x,i5,' aqueous solute',
     $  ' species are present on this file',/4x,'Coverage is ',a,
     $  ' per cent')
c
        if (nerr .gt. 0) then
          write (ux8,'(i5)') nerr
          call lejust(ux8)
          j2 = ilnobl(ux8)
          write (noutpt,1160) ux8(1:j2)
          write (nttyo,1160) ux8(1:j2)
          stop
        endif
c
c       Write hard core diameters and related parameters used
c       for example in the B-dot equation.
c
        call wrazp(azero,insgf,nazt,naztmx,ndata1,ndat1f,
     $  noutpt,nttyo,uazp)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Test and write the Pitzer data section of DATA1.
c
      if (uakey(1:8) .eq. 'Pitzer  ') then
c
c       Count the numbers of aqueous solution cations, anions, and
c       neutrals, excluding any fictive redox species. These will
c       be used to compute the species pairs and triplets relevant
c       to Pitzer parameters.
c
        call coasst(jassan,jassca,jassne,nat,natmax,uaqsp,zaqsp)
c
c       Compute the number of each of the relevant types of these
c       species pairs and triplets (e.g., ca, cc', aa', nc, na, nn,
c       nn', cca, aac, nnn, nnn', and nca).
c
        call cpcomb(jassan,jassca,jassne,naapr,nat,ncapr,nccpr,
     $  nnapr,nncpr,nnnpr,nn2pr,naactr,na2ctr,nncatr,nccatr,nc2atr,
     $  nn2ntr,nn3tr)
c
c       Allocate the associated arrays.
c
        ALLOCATE(alphca(ipbt_asv,ncapr))
c
        ALLOCATE(alamaa(jpfc_asv,0:ipbt_asv,naapr))
        ALLOCATE(alamca(jpfc_asv,0:ipbt_asv,ncapr))
        ALLOCATE(alamcc(jpfc_asv,0:ipbt_asv,nccpr))
        ALLOCATE(alamna(jpfc_asv,0:ipbt_asv,nnapr))
        ALLOCATE(alamnc(jpfc_asv,0:ipbt_asv,nncpr))
        ALLOCATE(alamnn(jpfc_asv,0:ipbt_asv,nnnpr))
        ALLOCATE(alamn2(jpfc_asv,0:ipbt_asv,nn2pr))
c
        ALLOCATE(amua2c(jpfc_asv,na2ctr))
        ALLOCATE(amuc2a(jpfc_asv,nc2atr))
        ALLOCATE(amun3(jpfc_asv,nn3tr))
c
        ALLOCATE(amuaac(jpfc_asv,naactr))
        ALLOCATE(amucca(jpfc_asv,nccatr))
        ALLOCATE(amunca(jpfc_asv,nncatr))
        ALLOCATE(amun2n(jpfc_asv,nn2ntr))
c
        ALLOCATE(qpdaa(naapr))
        ALLOCATE(qpdca(ncapr))
        ALLOCATE(qpdcc(nccpr))
        ALLOCATE(qpdna(nnapr))
        ALLOCATE(qpdnc(nncpr))
        ALLOCATE(qpdnn(nnnpr))
        ALLOCATE(qpdn2(nn2pr))
c
        ALLOCATE(qpdaac(naactr))
        ALLOCATE(qpdcca(nccatr))
        ALLOCATE(qpdnca(nncatr))
        ALLOCATE(qpdn2n(nn2ntr))
c
        ALLOCATE(in2pr(nn2pr))
        ALLOCATE(in3tr(nn3tr))
c
        ALLOCATE(iaapr(2,naapr))
        ALLOCATE(icapr(2,ncapr))
        ALLOCATE(iccpr(2,nccpr))
        ALLOCATE(inapr(2,nnapr))
        ALLOCATE(incpr(2,nncpr))
        ALLOCATE(innpr(2,nnnpr))
c
        ALLOCATE(in2ntr(2,nn2ntr))
        ALLOCATE(ia2ctr(2,na2ctr))
        ALLOCATE(ic2atr(2,nc2atr))
        ALLOCATE(iaactr(3,naactr))
        ALLOCATE(iccatr(3,nccatr))
        ALLOCATE(incatr(3,nncatr))
c
c       Construct index arrays for those pairs and triplets.
c
        call bldspc(iaapr,icapr,iccpr,inapr,incpr,innpr,
     $  in2pr,iaactr,ia2ctr,iccatr,ic2atr,incatr,in2ntr,in3tr,
     $  jassan,jassca,jassne,nat,natmax,naapr,ncapr,nccpr,nnapr,
     $  nncpr,nnnpr,nn2pr,naactr,na2ctr,nncatr,nccatr,nc2atr,
     $  nn2ntr,nn3tr,uaqsp,zaqsp)
c
c       Note on subroutine and array naming for species pairs and
c       triplets associated with Pitzer coefficients:
c
c         ca   is represented by "ca"
c         cc'a is represented by "cca"
c         cca  is represented by "c2a"
c         aa'c is represented by "a2c"
c         aac  is represented by "aac"
c         cc' is represented by "cc"
c
c       Test and process all ca (cation-anion) parameters. This includes
c       the cca (repeated cation-anion) and aac (cation-repeated anion)
c       cases. The basic transformations are:
c
c         beta(n)(ca) -> lambda(n)(ca)   (n = 0,2)
c         Cphi(ca)    -> mu(cca) and mu(aac)
c
        call tprca(abeta,acphi,alamca,alpha,alphca,amua2c,amuc2a,
     $  icapr,ipbtmx,jpfcmx,natmax,na2ctr,ncapr,ncvca,nc2atr,nerr,
     $  noutpt,npx2mx,npx2t,nttyo,nwarn,pcvca,qpdca,uaqsp,upair,zaqsp)
c
c       Test and process the cc' (cation-different cation) parameters.
c       cases. The basic transformation is:
c
c         theta(cc') -> lambda(cc')
c
        call tprcc(alamcc,atheta,iccpr,ipbtmx,jpfcmx,natmax,
     $  nccpr,ncvcc,nerr,noutpt,npx3mx,nthdt,nttyo,nwarn,pcvcc,
     $  qpdcc,uaqsp,uthdtr)
c
c       Test and process the aa' (anion-different anion) parameters.
c       cases. The basic transformation is:
c
c         theta(aa') -> lambda(aa')
c
        call tpraa(alamaa,atheta,iaapr,ipbtmx,jpfcmx,naapr,
     $  natmax,ncvaa,nerr,noutpt,npx3mx,nthdt,nttyo,nwarn,pcvaa,
     $  qpdaa,uaqsp,uthdtr)
c
c       Test and process the nn (repeated-neutral) parameters. This
c       includes the nnn (doubly repeated neutral) cases. Note on
c       subroutine and array naming:
c
c         nn  is represented by "n2"
c         nn' is represented by "nn"
c
c         lambda(nn) -> lambda(nn)
c         mu(nnn)    -> mu(nnn)
c
        call tprn2(abeta,acphi,alamn2,amun3,in2pr,ipbtmx,
     $  jpfcmx,natmax,ncvn2,nerr,nn2pr,nn3tr,noutpt,npx2mx,npx2t,
     $  nttyo,nwarn,pcvn2,qpdn2,uaqsp,upair)
c
c       Test and process the nn' (neutral-different neutral) parameters.
c       See the above note on subroutine and array naming.
c
c         lambda(nn') -> lambda(nn')
c
        call tprnn(abeta,acphi,alamnn,innpr,in2pr,ipbtmx,
     $  jpfcmx,natmax,ncvnn,nerr,nnnpr,nn2pr,noutpt,npx2mx,npx2t,
     $  nttyo,nwarn,pcvnn,qpdnn,qpdn2,uaqsp,upair)
c
c       Test and process the nc (neutral-cation) parameters.
c
c         lambda(nc) -> lambda(nc)
c
        call tprnc(abeta,alamnc,incpr,ipbtmx,jpfcmx,natmax,
     $  ncvnc,nerr,nncpr,noutpt,npx2mx,npx2t,nttyo,nwarn,pcvnc,
     $  qpdnc,uaqsp,upair)
c
c       Test and process the na (neutral-anion) parameters.
c
c         lambda(na) -> lambda(na)
c
        call tprna(abeta,alamna,inapr,ipbtmx,jpfcmx,natmax,
     $  ncvna,nerr,nnapr,noutpt,npx2mx,npx2t,nttyo,nwarn,pcvna,
     $  qpdna,uaqsp,upair)
c
c       Test and process the cc'a (cation-different cation-anion)
c       parameters. Note that mu(cc'a) depends not only on psi(cc'a),
c       but also on mu(cca) and mu(c'c'a), hence indirectly on Cphi(ca)
c       and Cphi(c'a). Note on subroutine and array naming:
c
c         cc'a is represented by "cca"
c         cca  is represented by "c2a"
c
c         psi(cc'a) -> mu(cc'a)
c
        call tprcca(amucca,amuc2a,apsi,icapr,iccatr,ipbtmx,
     $  jpfcmx,natmax,ncapr,nccatr,ncvcca,nc2atr,nerr,noutpt,
     $  npx3mx,npx3t,nttyo,nwarn,pcvcca,qpdca,qpdcca,uaqsp,
     $  utripl,zaqsp)
c
c       Test and process the aa'c (anion-different anion-cation)
c       parameters. Note that mu(aa'c) depends not only on psi(aa'c),
c       but also on mu(aac) and mu(a'a'c), hence indirectly on Cphi(ca)
c       and Cphi(ca'). Note on subroutine and array naming:
c
c         aa'c is represented by "aac"
c         aac  is represented by "a2c"
c
c         psi(aa'c) -> mu(aa'c)
c
        call tpraac(amuaac,amua2c,apsi,iaactr,icapr,ipbtmx,
     $  jpfcmx,naactr,natmax,na2ctr,ncapr,ncvaac,nerr,noutpt,
     $  npx3mx,npx3t,nttyo,nwarn,pcvaac,qpdaac,qpdca,uaqsp,
     $  utripl,zaqsp)
c
c       Test and process the nnn' (repeated neutral-different neutral)
c       parameters. These are not processed with the nn' parameters
c       because there are two mu coefficients (nnn' and n'n'n) for
c       each lambda coefficient (nn'). Note on subroutine and array
c       naming:
c
c         nnn' is represented by "n2n"
c
c         mu(nnn') -> mu(nnn')
c
        call tprn2n(amun2n,apsi,innpr,in2pr,in2ntr,ipbtmx,
     $  jpfcmx,natmax,ncvn2n,nerr,nnnpr,nn2pr,nn2ntr,noutpt,npx3mx,
     $  npx3t,nttyo,nwarn,pcvn2n,qpdnn,qpdn2,qpdn2n,uaqsp,utripl)
c
c       Test and process the nca (neutral-cation-anion) parameters.
c
c         mu(nca) -> mu(nca)
c
        call tprnca(amunca,apsi,inapr,incatr,incpr,ipbtmx,
     $  jpfcmx,natmax,ncvnca,nerr,nnapr,nncatr,nncpr,noutpt,npx3mx,
     $  npx3t,nttyo,nwarn,pcvnca,qpdna,qpdnca,qpdnc,uaqsp,utripl)
c
c
c       What about the following combinations?
c
c         cc     (repeated cation)
c         aa     (repeated neutral)
c         ccc    (doubly repeated cation)
c         aaa    (doubly repeated cation)
c         ccc'   (repeated cation and a distinct cation)
c         aaa'   (repeated anion and a distinct anion)
c         cc'c'' (three distinct cations)
c         aa'a'' (three distinct anions)
c         nn'n'' (three distinct neutrals)
c         ncc'   (neutral and two distinct cations)
c         naa'   (neutral and two distinct anions)
c         nn'c   (two distinct neutrals and a cation)
c         nn'a   (two distinct neutrals and an anion)
c
c       Some (e.g., cc,  aa) are defined to be zero by convention.
c       Others (e.g., cc'c'' and nn'n'') correspond to quaternary
c       or higher mixtures. In any case, they are not used in the
c       normal Pitzer treatment of electrolyte solutions.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c       Write the cation-anion (ca) data summary.
c
        ux8 = ' '
        write (ux8,'(f6.2)') pcvca
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1300)
 1300   format(//2x,'Cation-anion (ca) pair coverage:')
        write (noutpt,1310) ncvca,ncapr,ux8(1:j2)
 1310   format(/4x,i5,' pairs have Pitzer parameters specified on',
     $  ' the data file',/4x,i5,' pairs can be constructed from',
     $  ' the species present on this file',/4x,'Coverage is ',a,
     $  ' per cent')
c
c       Write the cation-distinct cation (cc') data summary.
c
        ux8 = ' '
        write (ux8,'(f6.2)') pcvcc
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1320)
 1320   format(//2x,"Cation-distinct cation (cc') pair coverage:")
        write (noutpt,1310) ncvcc,nccpr,ux8(1:j2)
c
c       Write the anion-distinct anion (aa') data summary.
c
        ux8 = ' '
        write (ux8,'(f6.2)') pcvaa
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1330)
 1330   format(//2x,"Anion-distinct anion (aa') pair coverage:")
        write (noutpt,1310) ncvaa,naapr,ux8(1:j2)
c
c       Write the repeated-neutral (nn) data summary.
c
        ux8 = ' '
        write (ux8,'(f6.2)') pcvn2
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1340)
 1340   format(//2x,'Repeated-neutral (nn) pair coverage:')
        write (noutpt,1310) ncvn2,nn2pr,ux8(1:j2)
c
c       Write the neutral-distinct neutral (nn') data summary.
c
        ux8 = ' '
        write (ux8,'(f6.2)') pcvnn
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1350)
 1350   format(//2x,"Neutral-distinct neutral (nn') pair coverage:")
        write (noutpt,1310) ncvnn,nnnpr,ux8(1:j2)
c
c       Write the neutral-cation (nc) data summary.
c
        ux8 = ' '
        write (ux8,'(f6.2)') pcvnc
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1360)
 1360   format(//2x,'Neutral-cation (nc) pair coverage:')
        write (noutpt,1310) ncvnc,nncpr,ux8(1:j2)
c
c       Write the neutral-anion (na) data summary.
c
        ux8 = ' '
        write (ux8,'(f6.2)') pcvna
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1370)
 1370   format(//2x,'Neutral-anion (na) pair coverage:')
        write (noutpt,1310) ncvna,nnapr,ux8(1:j2)
c
c       Write the cation-distinct cation-anion (cc'a) data summary.
c
        ux8 = ' '
        write (ux8,'(f6.2)') pcvcca
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1380)
 1380   format(//2x,"Cation-distinct cation-anion (cc'a) triplet",
     $  ' coverage:')
        write (noutpt,1390) ncvcca,nccatr,ux8(1:j2)
 1390   format(/4x,i5,' triplets have Pitzer parameters specified on',
     $  ' the data file',/4x,i5,' triplets can be constructed from',
     $  ' the species present on this file',/4x,'Coverage is ',a,
     $  ' per cent')
c
c
c       Write the anion-distinct anion-cation (aa'c) data summary.
c
        ux8 = ' '
        write (ux8,'(f6.2)') pcvaac
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1400)
 1400   format(//2x,"Anion-distinct anion-cation (aa'c) triplet",
     $  ' coverage:')
        write (noutpt,1390) ncvaac,naactr,ux8(1:j2)
c
c       Write the repeated neutral-distinct neutral (nnn') data
c       summary.
c
        ux8 = ' '
        write (ux8,'(f6.2)') pcvn2n
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1410)
 1410   format(//2x,"Repeated neutral-distinct neutral (nnn')",
     $  ' triplet coverage:')
        write (noutpt,1390) ncvn2n,nn2ntr,ux8(1:j2)
c
c       Write the neutral-cation-anion data summary.
c
        ux8 = ' '
        write (ux8,'(f6.2)') pcvnca
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1420)
 1420   format(//2x,'Neutral-cation-anion (nca) triplet coverage:')
        write (noutpt,1390) ncvnca,nncatr,ux8(1:j2)
c
        if (nerr .gt. 0) then
          write (ux8,'(i5)') nerr
          call lejust(ux8)
          j2 = ilnobl(ux8)
          write (noutpt,1160) ux8(1:j2)
          write (nttyo,1160) ux8(1:j2)
          stop
        endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c       Write Pitzer parameters.
c
        call wrpz23(alphca,alamaa,alamca,alamcc,alamna,alamnc,
     $  alamnn,alamn2,amuaac,amua2c,amucca,amuc2a,amunca,amun2n,
     $  amun3,ipbtmx,iaapr,icapr,iccpr,inapr,incpr,innpr,in2pr,
     $  iaactr,ia2ctr,iccatr,ic2atr,incatr,in2ntr,in3tr,jpdblo,
     $  jpfcmx,natmax,naapr,ncapr,nccpr,nnapr,nncpr,nnnpr,nn2pr,
     $  naactr,na2ctr,nccatr,nc2atr,nncatr,nn2ntr,nn3tr,ndata1,
     $  ndat1f,noutpt,nttyo,uaqsp,uethfl)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      j2 = ilnobl(uakey)
      j3 = ilnobl(utitld(1))
      write (nttyo,1460) uakey(1:j2),utitld(1)(1:j3)
      write (noutpt,1460) uakey(1:j2),utitld(1)(1:j3)
 1460 format(//' Completed processing the ',a,' data file ',a,'.')
      write (noutpt,1470)
      write (nttyo,1470)
 1470 format(//' No errors were encountered.')
      if (nwarn .le. 0) then
        write (noutpt,1480)
        write (nttyo,1480)
 1480   format(/' No warnings were encountered.',//)
      else
        write (ux8,'(i5)') nwarn
        call lejust(ux8)
        j4 = ilnobl(ux8)
        write (noutpt,1490) ux8(1:j4)
        write (nttyo,1490) ux8(1:j4)
 1490   format(/' ',a,' warning(s) were encountered.',//)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get end time and date. Also get the run time. Optionally, on
c     a Unix platform, get the user and cpu times.
c
      call runtim(iexec0,jexec0,texec0,noutpt,nttyo,trun,
     $ tuser,tcpu,udate1,utime1)
c
      j2 = ilnobl(udate0)
      j3 = ilnobl(udate1)
      write (noutpt,1500) utime0,udate0(1:j2),utime1,udate1(1:j3)
      write (nttyo,1500) utime0,udate0(1:j2),utime1,udate1(1:j3)
 1500 format(10x,'Start time = ',a8,2x,a,/12x,
     $ 'End time = ',a8,2x,a)
c
c     Print the run, user, and cpu times.
c
      write (noutpt,1510) trun
      write (nttyo,1510)  trun
 1510 format(/10x,' Run time = ',g10.3,' seconds')
      if (tuser .gt. 0.) then
        write (noutpt,1520) tuser
        write (nttyo,1520)  tuser
 1520   format(10x,'User time = ',g10.3,' seconds')
      endif
      if (tcpu .gt. 0.) then
        write (noutpt,1530) tcpu
        write (nttyo,1530)  tcpu
 1530   format(10x,' Cpu time = ',g10.3,' seconds')
      endif
c
      write (noutpt,1540)
      write (nttyo,1540)
 1540 format(/' Normal exit')
c
c     Clear the IEEE flag for floating-point underflow, if such a
c     flag is present, to avoid getting an unnecessary system
c     warning message. Underflow is a normal condition in EQ3/6.
c     Make porting changes in the EQLIBU subroutine that is called
c     in this section. Do not make the porting changes here.
c
      call cliefu()
c
c     Close files.
c
      close(ndat0s,status='delete')
c
      close(noutpt)
      close(nslist)
      close(ndata1)
      close(ndat1f)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
