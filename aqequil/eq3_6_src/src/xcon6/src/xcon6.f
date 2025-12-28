      program xcon6
c
c     This is the main program of the XCON6 code. Configuration
c     identification, the copyright statement, legal disclaimers,
c     and similar statements are contained in XCON6/aaaxc6.f, the
c     lead-off subroutine in the XCON6 source code. A short description
c     of this program is also contained in that subroutine.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      include 'xcon6/x6op7.h'
      include 'eqlib/eqlo8.h'
c
c-----------------------------------------------------------------------
c
      integer ietpar,iktpar,imapar,imchpa,jetpar,kpar,nbtpar,nbt1pa,
     $ nctpar,ndctpa,nertpa,netpar,nffgpa,nodbpa,nopgpa,noprpa,noptpa,
     $ nordmx,nprppa,nprspa,nptkpa,nrctpa,nsscpa,nsrtpa,ntitpa,nttkpa,
     $ nxmdpa,nxoppa,nxpepa,nxrtpa
c
      parameter(ietpar = 10,iktpar = 20,imapar = 3,imchpa = 4)
      parameter(jetpar = 3,nertpa = 5,netpar = 12,ntitpa = 100)
      parameter(nctpar = 110,ndctpa = 4,nffgpa = 10)
      parameter(nprppa = 100,nrctpa = 25,nbtpar = 250)
      parameter(nptkpa = 3,nsrtpa = 5,nxmdpa = 40,nsscpa = 6)
      parameter(nxoppa = 20,nxpepa = 40,nxrtpa = 20,nttkpa  =  3)
      parameter(nodbpa = 20,nopgpa = 20,noprpa = 20,noptpa = 20)
      parameter(nbt1pa = nbtpar + 1,nprspa = 3*nprppa)
      parameter(kpar = 2*nbtpar)
c
      integer ietmax,iktmax,imamax,imchmx,jetmax,kmax,nbtmax,nbt1mx,
     $ nctmax,ndctmx,nertmx,netmax,nffgmx,nodbmx,nopgmx,noprmx,noptmx,
     $ nprpmx,nprsmx,nptkmx,nrctmx,nsscmx,nsrtmx,ntitmx,nttkmx,nxmdmx,
     $ nxopmx,nxpemx,nxrtmx
c
c-----------------------------------------------------------------------
c
      integer newin,ninpt,ninpts,noutpt,nttyo,nxcon
c
      integer iact(imchpa,2,nrctpa),ibsrti(nsrtpa),iesrti(nsrtpa),
     $ igerti(jetpar,nertpa),iktbt(nxrtpa),imech(2,nrctpa),iodb(nodbpa),
     $ iopg(nopgpa),iopr(noprpa),iopt(noptpa),ipivot(imapar),
     $ ixrti(nxrtpa),jcode(nrctpa),jflgi(nbtpar),jgerti(nertpa),
     $ jgext(netpar),jreac(nrctpa),jxmod(nxmdpa),kxmod(nxmdpa),
     $ ndact(imchpa,2,nrctpa),nelptr(nbtpar),nesrbt(nsrtpa),
     $ ngexrt(jetpar,netpar),nrk(2,nrctpa),nsk(nrctpa)
c
      integer ifile,ioscan,itermx,jpress,jtemp,kbt,kct,kdim,km1,kmt,
     $ kprs,ksq,ksplmx,ksppmx,kstpmx,kx1,kxt,nbti,nffg,nmodl1,nmodl2,
     $ nobswt,nordlm,npslmx,nprmn,nprmx,nrct,nsbswt,nsslmx,ntitl1,
     $ ntitl2,ntrymx,nxmod,nxopex,nxopt
c
      integer i,icount,ier,j,j2,j3,k,kl,ko2,krow,n,nbb,nn,nobsw,nbbh,
     $ nbbo,ncb,ncbh,ncbo,ncbt,nce,nceh,nceo,nelspt,nerr,nert,net,nmax,
     $ no2,nprob,nprpti,nprsti,nrc,nsr,nsrt,nxrt,ntest
c
      integer ilnobl
c
      logical qcntmp,qcnpre,qend,qex,qgexsh,qpr,qrderr,qxcon,q8bchk,
     $ q8beta
c
      character*80 utitl1(ntitpa),utitl2(ntitpa)
      character*56 ugexr(ietpar,jetpar,netpar)
      character*48 ubmtbi(nbtpar),uobsw(2,nbtpar),uprs(nprspa),
     $ uprspi(nprspa),usbsw(2,nbtpar),uxmod(nxmdpa),uzveci(kpar)
      character*24 ubsri(nbt1pa,nsrtpa),ucxri(iktpar,nxrtpa),
     $ udac(ndctpa,imchpa,2,nrctpa),uelspn(nctpar),
     $ uendb(iktpar,nxrtpa),uffg(nffgpa),ugermo(nertpa),
     $ ugersi(ietpar,jetpar,nertpa),ugexmo(netpar),ugexp(netpar),
     $ undms(kpar),unrms(kpar),uprphi(nprppa),ureac(nrctpa),
     $ uxcat(nxoppa),uxmd24(nxmdpa),uxopex(nxpepa)
      character*16 uxct16(nxoppa)
      character*8 uelemb(nctpar),uelnam(nctpar),uelnlc(nctpar),
     $ uesrb(nctpar,nsrtpa),uesri(nctpar,nsrtpa),ugerji(jetpar,nertpa),
     $ ugexj(jetpar,netpar),uhfgex(ietpar,jetpar,netpar),
     $ uvfgex(ietpar,jetpar,netpar),uxkgex(ietpar,jetpar,netpar),
     $ uxopt(nxoppa)
c
      character*80 uline
      character*24 uacion,uaqsln,unamph,ux24
      character*8 uplatc,uplatm,ustelu,ustxc6,uveelu,uvexc6
      character*8 ucode,urelno,ustage,ueqlrn,ueqlst,uv,ux
      character*3 uoldv,uoldvd,unewv
      character*1 uoldf,unewf,ux1
c
      real*8 apresh(5,2)
c
      real*8 aamatr(imapar,imapar),cdac(ndctpa,imchpa,2,nrctpa),
     $ cbsri(nbt1pa,nsrtpa),celspe(nctpar),celsph(nctpar),
     $ celspo(nctpar),cesrb(nctpar,nsrtpa),cesri(nctpar,nsrtpa),
     $ cgexj(jetpar,netpar),csigma(imchpa,2,nrctpa),
     $ eact(imchpa,2,nrctpa),egersi(ietpar,jetpar,nertpa),fk(nrctpa),
     $ fkrc(nrctpa),gmmatr(imapar,imapar),hact(imchpa,2,nrctpa),
     $ modr(nrctpa),moffg(nffgpa),morr(nrctpa),mprphi(nprppa),
     $ mprs(nprspa),mprspi(nprspa),mtbaqi(nbtpar),mtbi(nbtpar),
     $ mteaqb(nctpar),mteb(nctpar),mwtges(netpar)
c
      real*8 ptk(nptkpa),rkb(imchpa,2,nrctpa),rk0(imchpa,2,nrctpa),
     $ rxbarb(iktpar,nxrtpa),rxbari(iktpar,nxrtpa),sfcar(nrctpa),
     $ sscrew(nsscpa),sk(nrctpa),ssfcar(nrctpa),tgexp(netpar),
     $ ttk(nttkpa),trkb(imchpa,2,nrctpa),trk0(imchpa,2,nrctpa),
     $ vreac(nrctpa),xgersi(ietpar,jetpar,nertpa),
     $ xhfgex(ietpar,jetpar,netpar),xlkgex(ietpar,jetpar,netpar),
     $ xlkffg(nffgpa),xlkmod(nxmdpa),xvec(imapar),
     $ xvfgex(ietpar,jetpar,netpar),yvec(imapar),zvclgi(kpar),
     $ zelsp(nctpar),zgexj(jetpar,netpar)
c
      real*8 awmaxi,awmini,cplim,dlaplo,dlaprn,dleplo,dleprn,dlhplo,
     $ dlhprn,dloplo,dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,
     $ dlxpll,dlxplo,dlxprl,dlxprn,dlzmx1,dlzmx2,dlzidp,dzpllg,dzplot,
     $ dzprlg, dzprnt,ehmaxi,ehmini,electr,o2maxi,o2mini,phmaxi,phmini,
     $ tempci,tempc0,timemx,timmxi,presh,preshi,pressb,pressi,tempcb,
     $ tistti,tolbt,toldl,tolsat,tolsst,tolx,tolxsf,tstrt,ximaxi,
     $ xistti,zimax,zistrt,zkfac,zklogl,zklogu
c
      real*8 ctxh,ctxo,cx,electc,epstst,mtaqsm,mxc,verold,vernew,ztx
c
c-----------------------------------------------------------------------
c
c     BEGIN_MACHINE_DEPENDENT_CODE
c
c       On some systems, a BLOCK DATA routine must be declared in an
c       EXTERNAL statement to assure proper loading. On some other
c       systems, this is not necessary, but neither it is not harmful.
c       On yet some other systems, the EXTERNAL statement below may
c       cause a problem. If so, try commenting it out. If you still
c       have trouble, consult your local system documentation or
c       experiment to find out how to correctly handle a BLOCK DATA
c       routine on your system. The EXTERNAL statement below should not
c       cause a problem if you are using a compiler which is fully
c       compliant with the Fortran 90 standard. However, there is
c       no guarantee that it will be adequate to assure correct loading
c       of the BLOCK DATA routine.
c
        external bkdxc6
c
c     END_MACHINE_DEPENDENT_CODE
c
c-----------------------------------------------------------------------
c
      data noutpt /0/
c
      data uaqsln /'Aqueous solution        '/
c
      data apresh(1,1) / 1.013200000E+00/,
     $     apresh(2,1) / 0./,apresh(3,1) / 0./,
     $     apresh(4,1) / 0./,apresh(5,1) / 0./,
     $     apresh(1,2) /-4.345000000E-01/,
     $     apresh(2,2) / 7.632333333E-03/,
     $     apresh(3,2) / 5.514000000E-05/,
     $     apresh(4,2) /-1.263733333E-06/,
     $     apresh(5,2) / 1.396800000E-08/
c
      data epstst /1.e-12/
c
c-----------------------------------------------------------------------
c
c     BEGIN_MACHINE_DEPENDENT_CODE
c
c       Define the console device number.
c
c       The following works on UNIX, PC, and VAX machines.
c
        data nttyo  /6/
c
c       BEGIN_MAC_DEPENDENT_CODE
c
c         data nttyo  /9/
c
c       END_MAC_DEPENDENT_CODE
c
c     END_MACHINE_DEPENDENT_CODE
c
c-----------------------------------------------------------------------
c
      include 'xcon6/x6sbsp.h'
c
      data ninpt  /9/,ninpts /10/,newin /11/,nxcon /12/
c
c-----------------------------------------------------------------------
c
c     Get configuration identification data.
c
      call aaaxc6(ustxc6,uvexc6)
      call aaaelu(ustelu,uveelu)
      call platfd(uplatc,uplatm)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set dimensioning variables.
c
      ietmax = ietpar
      iktmax = iktpar
      imamax = imapar
      imchmx = imchpa
      jetmax = jetpar
      kmax = kpar
      nbtmax = nbtpar
      nbt1mx = nbt1pa
      nctmax = nctpar
      ndctmx = ndctpa
      nertmx = nertpa
      netmax = netpar
      nffgmx = nffgpa
      nodbmx = nodbpa
      nopgmx = nopgpa
      noprmx = noprpa
      noptmx = noptpa
      nprpmx = nprppa
      nprsmx = nprspa
      nptkmx = nptkpa
      nrctmx = nrctpa
      nsrtmx = nsrtpa
      nsscmx = nsscpa
      ntitmx = ntitpa
      nttkmx = nttkpa
      nxmdmx = nxmdpa
      nxopmx = nxoppa
      nxpemx = nxpepa
      nxrtmx = nxrtpa
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Does the old input file exist?
c
      inquire(file="input",exist=qex)
      if (.not.qex) then
        write (nttyo,1000)
 1000   format(' * Error- (XCON6/xcon6) The INPUT file (the old input',
     $  ' file).',/7x,"doesn't exist. Check the name that was",
     $  ' specified.')
        go to 999
      endif
c
c     Open the old input file.
c
      open(ninpt,file='input',status='old',err=700)
      go to 710
c
  700 write (nttyo,1010)
 1010 format(" * Error - (XCON6/xcon6) Can't open the INPUT file (the",
     $ ' old input file).')
      go to 999
c
  710 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Determine the old input file format.
c       W = compact
c       D = menu-style
c
      k = 0
      kl = 0
      ntest = 10
      uoldf = 'W'
      do 800 n = 1,100
        read (ninpt,1060) ux1
 1060   format(a1)
        if (ux1(1:1) .ne. 'c') kl = kl + 1
        if (ux1(1:1) .eq. '|') k = k + 1
        if (kl .ge. ntest) then
          if (k .ge. ntest) uoldf = 'D'
          go to 810
        endif
  800 continue
  810 continue
      rewind(ninpt)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Does the INPUTS file (the stripped input file) exist?
c     If so, kill it. This file will be used to contain a copy of the
c     old input file, stripped of comment lines.
c
      inquire(file="inputs",exist=qex)
      if (qex) then
        open(ninpts,file='inputs',status='old',err=720)
        close(ninpts,status='delete')
      endif
c
c     Open the INPUTS file.
c
      open(ninpts,file='inputs',status='new',err=720)
      go to 730
c
  720 write (nttyo,1020)
 1020 format(" * Error - (XCON6/xcon6) Can't open the INPUTS file",
     $ /7x,'(the stripped old input file).')
      close(ninpt)
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Strip comment lines from the input file. Also strip any blank
c     lines from a file in "D" format.
c
  730 if (uoldf(1:1) .eq. 'W') then
  740   read (ninpt,1030,end=750) uline
 1030   format(a80)
        if (uline(1:1) .ne. '*') write (ninpts,1030) uline
        go to 740
  750   continue
      elseif (uoldf(1:1) .eq. 'D') then
  755   read (ninpt,1030,end=757) uline
        j2 = ilnobl(uline)
        if (j2 .gt. 0) then
          if (uline(1:1).ne.'*' .and. uline(1:1).ne.'c')
     $    write (ninpts,1030) uline
        endif
        go to 755
  757   continue
      endif
c
c     Rewind the stripped input file.
c
      rewind(ninpts)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Does the NEWIN file (to be the new input file) exist?
c     If so, kill it.
c
      inquire(file="newin",exist=qex)
      if (qex) then
        open(newin,file='newin',status='old',err=760)
        close(newin,status='delete')
      endif
c
c     Open the new input file.
c
      open(newin,file='newin',status='new',err=760)
      go to 770
c
  760 write (nttyo,1050)
 1050 format(" * Error - (XCON6/xcon6) Can't open the NEWIN file",
     $ ' (to be the',/7x,'new input file.)')
      close(ninpt)
      close(ninpts,status='delete')
      go to 999
c
  770 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Input file formats:
c
c       W = compact
c       D = menu-style
c
c        uoldf  = format of the old file
c        unewf  = format of the new file
c
      unewf = ' '
c
c     Input file version levels.
c       '6.0' = version 6.0 (including version 6.1)
c       '7.0' = version 7.0 (including versions 7.0x and 7.1)
c       '7.2' = version 7.2
c       '8.0' = version 8.0
c
c         uoldv  = version level of the old file
c         uoldvd = default version level of the old file
c         unewv  = version level of the new file
c
      uoldv = ' '
      uoldvd = ' '
      unewv = ' '
c
c     Set ultimate default values.
c
      uoldvd = '7.2'
      unewv = '8.0'
      unewf = 'D'
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Does an options file exist?
c
      inquire(file="ixcon",exist=qxcon)
      if (qxcon) then
        open(nxcon,file='ixcon',status='old',err=780)
        go to 790
c
  780   write (nttyo,1055)
 1055   format(/" * Error - (XCON6/xcon6) Can't open the IXCON options",
     $  ' file.')
        close(ninpts,status='delete')
        close(ninpt)
        close(newin)
        go to 999
c
c       Read the IXCON options file.
c
  790  call rddixc(nxcon,uoldvd,unewf,unewv)
      else
        write (nttyo,1057)
 1057   format(/' * Error - (XCON6/xcon6) The IXCON options file',
     $  " doesn't exist.")
        close(ninpt)
        close(ninpts,status='delete')
        close(newin)
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Determine the old input file version level.
c       '6.0' = version 6.0 (including version 6.1)
c       '7.0' = version 7.0 (including versions 7.0x and 7.1)
c       '7.2' = version 7.2
c       '8.0' = version 8.0
c
      i = 0
      do 820 n = 1,ntitmx
        read (ninpts,1030) uline
        if (uoldf(1:1) .eq. 'W') then
          if (uline(1:8) .eq. 'endit.  ') go to 825
        elseif (uoldf(1:1) .eq. 'D') then
          if (uline(1:8) .eq. '|-------') i = i + 1
          if (i .ge. 3) go to 825
        endif
        j = index(uline,'Version level=')
        if (j .le. 0) j = index(uline,'version level=')
        if (j .gt. 0) then
          uv = uline(j + 14:j + 21)
          go to 840
        endif
  820 continue
  825 continue
c
      i = 0
      rewind(ninpts)
      do 830 n = 1,ntitmx
        read (ninpts,1030) uline
        if (uoldf(1:1) .eq. 'W') then
          if (uline(1:8) .eq. 'endit.  ') go to 835
        elseif (uoldf(1:1) .eq. 'D') then
          if (uline(1:8) .eq. '|-------') i = i + 1
          if (i .ge. 3) go to 835
        endif
        j = index(uline,'Version number=')
        if (j .le. 0) j = index(uline,'version number=')
        if (j .gt. 0) then
          uv = uline(j + 15:j + 22)
          go to 840
        endif
  830 continue
  835 continue
c
      write (nttyo,1070)
 1070 format(/' * Warning - (XCON6/xcon6) The title on the old input',
     $ /7x,"file doesn't contain a version level tag. If this code is",
     $ /7x,'unable to determine the correct version level by other',
     $ /7x,'means, add a version tag by putting "Version level= X.X",',
     $ /7x,'where "X.X" is the version level, in the title, preferably',
     $ /7x,'on line 3. If there is more than problem on a single input',
     $ /7x,'file, the tag need be put only in the title of the first',
     $ ' problem.')
      go to 880
c
  840 call lejust(uv)
c
      if (uv(1:3) .eq. '6.0') then
        uoldv = '6.0'
      elseif (uv(1:1) .eq. '6') then
        uoldv = '6.0'
      elseif (uv(1:3) .eq. '7.0') then
        uoldv = '7.0'
      elseif (uv(1:3) .eq. '7.1') then
        uoldv = '7.0'
      elseif (uv(1:3) .eq. '7.2') then
        uoldv = '7.2'
      elseif (uv(1:1) .eq. '7') then
        uoldv = '7.0'
      elseif (uv(1:3) .eq. '8.0') then
        uoldv = '8.0'
      elseif (uv(1:1) .eq. '8') then
        uoldv = '8.0'
      else
        j2 = ilnobl(uv)
        write (nttyo,1080) uv(1:j2)
 1080   format(/" * Warning - (XCON6/xcon6) Can't determine the",
     $  /7x,'version level of the old input file. The version level',
     $  /7x,'is specified by a tag in the title of the first problem',
     $  /7x,'on the file as "',a,'". This is not a valid version',
     $  /7x,'level descriptor. Valid descriptors include "6.0", "7.0",',
     $  /7x,'"7.2", and "8.0". If the code is unable to determine the',
     $  /7x,'correct version level, replace the invalid descriptor.')
      endif
c
  880 if (uoldv(1:3) .eq. '   ') then
        rewind(ninpts)
        do 837 n = 1,1000
          read (ninpts,1030,end=838) uline
          j = index(uline,'uacion= ')
          if (j .gt. 0) then
            uoldv = '6.0'
            go to 838
          endif
  837   continue
  838   continue
      endif
c
c     If necessary, use the default for the version level of the old
c     input file.
c
      if (uoldv(1:3) .eq. '   ') then
        uoldv = uoldvd
        write (nttyo,1090) uoldvd
 1090   format(/' * Warning - (XCON6/xcon6) Taking the default',
     $  /7x,'value of "',a3,'" for the version level of the old',
     $  /7x,'input file.')
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Count the number of chemical elements which can be mapped to
c     basis species. These are listed on the pseudo-data file
c     in x6sbsp.h.
c
      nelspt = 0
      do nce = 1,nctmax
        if (uelnam(nce)(1:8) .eq. 'endit.  ') go to 240
        nelspt = nelspt + 1
      enddo
  240 continue
c
c     Create a copy of the uelnam array in lower case.
c
      do nce = 1,nelspt
        ux = uelnam(nce)
        call locase(ux)
        uelnlc(nce) = ux
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Find the positions of O and H in the pseudo-data file in x6sbsp.h.
c
      nceo = 0
      do nce = 1,nelspt
        if (uelnam(nce)(1:2).eq.'O ' .or.
     $    uelnam(nce)(1:2).eq.'o ') then
          nceo = nce
          go to 300
        endif
      enddo
      ux = 'O'
      write (nttyo,1740) ux(1:1)
 1740 format(/" * Error - (XCON6/xcon6) Can't find the chemical",
     $ /7x,'element ',a,' listed in the pseudo-data file in x6sbsp.h.')
      stop
  300 continue
c
      nceh = 0
      do nce = 1,nelspt
        if (uelnam(nce)(1:2).eq.'H ' .or.
     $    uelnam(nce)(1:2).eq.'h ') then
          nceh = nce
          go to 310
        endif
      enddo
      ux = 'H'
      write (nttyo,1740) ux(1:1)
      stop
  310 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Process the old input file, working from the stripped copy.
c
      rewind(ninpts)
      q8bchk = .false.
      q8beta = .false.
      nprob = 0
  100 nprob = nprob + 1
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Zero or null various variables.
c
      call initcb(utitl1,ntitmx)
      call initcb(utitl2,ntitmx)
c
      call initiz(iopt,noptmx)
      call initiz(iopg,nopgmx)
      call initiz(iopr,noprmx)
      call initiz(iodb,nodbmx)
c
      call initaz(ptk,nptkmx)
      call initaz(ttk,nttkmx)
c
      call initcb(uxopt,nxopmx)
      call initcb(uxcat,nxopmx)
      call initcb(uxct16,nxopmx)
      call initcb(uxopex,nxpemx)
c
      call initcb(uxmod,nxmdmx)
      call initcb(uxmd24,nxmdmx)
      call initiz(jxmod,nxmdmx)
      call initiz(kxmod,nxmdmx)
      call initaz(xlkmod,nxmdmx)
c
      call initcb(uffg,nffgmx)
      call initaz(xlkffg,nffgmx)
      call initaz(moffg,nffgmx)
c
      call initaz(sscrew,nsscmx)
c
      call initcb(ureac,nrctmx)
      call initiz(jcode,nrctmx)
      call initiz(jreac,nrctmx)
      call initaz(morr,nrctmx)
      call initaz(modr,nrctmx)
      call initiz(nsk,nrctmx)
      call initaz(sk,nrctmx)
      call initaz(sfcar,nrctmx)
      call initaz(ssfcar,nrctmx)
      call initaz(fk,nrctmx)
      call initaz(fkrc,nrctmx)
      call initiz(imech,nrctmx)
c
      nmax = nrctmx*2
      call initiz(nrk,nmax)
c
      nmax = nrctmx*2*imchmx
      call initiz(ndact,nmax)
      call initiz(iact,nmax)
      call initaz(rkb,nmax)
      call initaz(rk0,nmax)
      call initaz(csigma,nmax)
      call initaz(trkb,nmax)
      call initaz(trk0,nmax)
      call initaz(eact,nmax)
      call initaz(hact,nmax)
c
      nmax = nrctmx*2*imchmx*ndctmx
      call initcb(udac,nmax)
      call initaz(cdac,nmax)
c
      call initiz(iktbt,nxrtmx)
c
      nmax = iktmax*nxrtmx
      call initcb(uendb,nmax)
      call initcb(ucxri,nmax)
      call initaz(rxbarb,nmax)
      call initaz(rxbari,nmax)
c
      call initiz(iesrti,nsrtmx)
      call initiz(nesrbt,nsrtmx)
c
      nmax = nctmax*nsrtmx
      call initcb(uesrb,nmax)
      call initcb(uesri,nmax)
      call initaz(cesrb,nmax)
      call initaz(cesri,nmax)
c
      nmax = nbt1mx*nsrtmx
      call initcb(ubsri,nmax)
      call initaz(cbsri,nmax)
c
      call initcb(uelemb,nctmax)
      call initaz(mteb,nctmax)
      call initaz(mteaqb,nctmax)
c
      call initcb(unrms,kmax)
      call initcb(undms,kmax)
      call initaz(zvclgi,kmax)
c
      call initcb(uzveci,kmax)
c
      nmax = 2*nbtmax
      call initcb(uobsw,nmax)
      call initcb(usbsw,nmax)
c
      call initcb(ubmtbi,nbtmax)
      call initaz(mtbi,nbtmax)
      call initaz(mtbaqi,nbtmax)
      call initiz(jflgi,nbtmax)
c
      call initcb(uprs,nprsmx)
      call initaz(mprs,nprsmx)
c
      call initcb(uprphi,nprpmx)
      call initaz(mprphi,nprpmx)
      call initcb(uprspi,nprsmx)
      call initaz(mprspi,nprsmx)
c
      dlzmx1 = 0.
      dlxmx0 = 0.
      tolx = 0.
      tolxsf = 0.
c
      dlzmx2 = 0.
      zkfac = 0.
      zklogl = 0.
      zklogu = 0.
      tolsst = 0.
      npslmx = 0
      nsslmx = 0
      nordlm = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The uacion variable only appears on version level '6.0' input
c     files. Provide a blank default value.
c
      uacion = ' '
c
c     The following parameters appear on 'W' format input files,
c     but not on 'D' format files. Provide the usual default values.
c
      dzplot = 1.e+38
      dzpllg = 1.e+38
      ksplmx = 10000
      ifile = 60
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The uacion variable only appears on version level '6.0' input
c     files. Provide a blank default value.
c
      uacion = ' '
c
c     The following parameters appear on 'W' format input files,
c     but not on 'D' format files. Provide the usual default values.
c
      dzplot = 1.e+38
      dzpllg = 1.e+38
      ksplmx = 10000
      ifile = 60
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the current problem on the stripped input file.
c
      if (uoldf(1:1) .eq. 'W') then
        if (uoldv(1:3) .eq. '6.0') then
c
          call rd6w6(cdac,cesrb,cplim,csigma,dlzmx1,dlzmx2,dlzidp,
     $    dzpllg,dzplot,dzprlg,dzprnt,electr,fk,ifile,iktbt,iktmax,
     $    imchmx,imech,iodb,iopg,iopr,iopt,ioscan,itermx,jcode,jreac,
     $    jtemp,jxmod,kct,kdim,kmax,kmt,kprs,ksq,ksplmx,ksppmx,kstpmx,
     $    kxmod,kxt,modr,moffg,morr,mprs,mteaqb,mteb,nctmax,ndact,
     $    ndctmx,nesrbt,nffg,nffgmx,ninpts,nmodl1,nmodl2,nodbmx,nopgmx,
     $    noprmx,noptmx,nordlm,npslmx,nprmn,nprmx,nprsmx,nrct,nrctmx,
     $    nrk,nsk,nsrtmx,nsscmx,nsslmx,ntitl1,ntitl2,ntitmx,ntrymx,
     $    nttkmx,nttyo,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrtmx,
     $    qend,qrderr,rk0,rxbarb,sscrew,sk,tempci,tempc0,timemx,tolbt,
     $    toldl,tolsat,tolsst,tolx,tstrt,ttk,uacion,udac,uelemb,uendb,
     $    uesrb,uffg,undms,unrms,uprs,ureac,utitl1,utitl2,uxct16,uxmd24,
     $    uxopex,uxopt,vreac,xlkffg,xlkmod,zimax,zistrt,zkfac,zklogl,
     $    zklogu,zvclgi)
          if (qrderr) go to 990
c
        elseif (uoldv(1:3) .eq. '7.0') then
c
          call rd6w7(cdac,cesrb,cplim,csigma,dlzmx1,dlzmx2,dlzidp,
     $    dzpllg,dzplot,dzprlg,dzprnt,eact,electr,fk,iact,ifile,iktbt,
     $    iktmax,imchmx,imech,iodb,iopg,iopr,iopt,ioscan,itermx,jcode,
     $    jreac,jtemp,jxmod,kct,kdim,kmax,kmt,kprs,ksq,ksplmx,ksppmx,
     $    kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprs,mteaqb,mteb,nctmax,
     $    ndact,ndctmx,nesrbt,nffg,nffgmx,ninpts,nmodl1,nmodl2,nodbmx,
     $    nopgmx,noprmx,noptmx,nordlm,npslmx,nprmn,nprmx,nprsmx,nrct,
     $    nrctmx,nrk,nsk,nsrtmx,nsscmx,nsslmx,ntitl1,ntitl2,ntitmx,
     $    ntrymx,nttkmx,nttyo,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,
     $    nxrtmx,qend,qrderr,rk0,rxbarb,sscrew,sk,tempci,tempc0,timemx,
     $    tolbt,toldl,tolsat,tolsst,tolx,trk0,tstrt,ttk,udac,uelemb,
     $    uendb,uesrb,uffg,undms,unrms,uprs,ureac,utitl1,utitl2,uxct16,
     $    uxmd24,uxopex,uxopt,vreac,xlkffg,xlkmod,zimax,zistrt,zkfac,
     $    zklogl,zklogu,zvclgi)
          if (qrderr) go to 990
c
        elseif (uoldv(1:3) .eq. '7.2') then
c
          call rd6w72(cdac,cesrb,cplim,csigma,dlzmx1,dlzmx2,dlzidp,
     $    dzpllg,dzplot,dzprlg,dzprnt,eact,electr,fk,iact,ifile,iktbt,
     $    iktmax,imchmx,imech,iodb,iopg,iopr,iopt,ioscan,itermx,jcode,
     $    jreac,jtemp,jxmod,kct,kdim,kmax,kmt,kprs,ksq,ksplmx,ksppmx,
     $    kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprs,mteaqb,mteb,nctmax,
     $    ndact,ndctmx,nesrbt,nffg,nffgmx,ninpts,nmodl1,nmodl2,nodbmx,
     $    nopgmx,noprmx,noptmx,nordlm,npslmx,nprmn,nprmx,nprsmx,nrct,
     $    nrctmx,nrk,nsk,nsrtmx,nsscmx,nsslmx,ntitl1,ntitl2,ntitmx,
     $    ntrymx,nttkmx,nttyo,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,
     $    nxrtmx,qend,qrderr,rk0,rxbarb,sscrew,sk,tempci,tempc0,timemx,
     $    tolbt,toldl,tolsat,tolsst,tolx,trk0,tstrt,ttk,udac,uelemb,
     $    uendb,uesrb,uffg,undms,unrms,uprs,ureac,utitl1,utitl2,uxct16,
     $    uxmd24,uxopex,uxopt,vreac,xlkffg,xlkmod,zimax,zistrt,zkfac,
     $    zklogl,zklogu,zvclgi)
          if (qrderr) go to 990
c
        elseif (uoldv(1:3) .eq. '8.0') then
c
          call rd6w8(awmaxi,awmini,cbsri,cdac,cesri,cgexj,
     $    csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,
     $    dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,
     $    dlxplo,dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,
     $    iact,ibsrti,iesrti,ietmax,iktmax,imchmx,imech,iodb,iopg,
     $    iopr,iopt,igerti,itermx,ixrti,jcode,jgerti,jetmax,jflgi,
     $    jgext,jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,
     $    ksppmx,kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,
     $    mtbaqi,mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,
     $    nert,nertmx,net,netmax,nffg,nffgmx,ngexrt,ninpts,nobswt,
     $    nodbmx,nopgmx,noprmx,noptmx,nordmx,noutpt,nprob,nprpmx,
     $    nprpti,nprsmx,nprsti,nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,
     $    nsrt,nsrtmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,nttyo,nxmdmx,
     $    nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrt,nxrtmx,o2maxi,o2mini,
     $    phmaxi,phmini,pressb,pressi,ptk,qend,qgexsh,qrderr,rkb,
     $    rxbari,sfcar,ssfcar,tempcb,tempci,tgexp,timmxi,tistti,tolbt,
     $    toldl,tolsat,tolxsf,trkb,ttk,ubmtbi,ubsri,ucxri,udac,uesri,
     $    uffg,ugerji,ugermo,ugersi,ugexj,ugexmo,ugexp,ugexr,uhfgex,
     $    uobsw,uprphi,uprspi,ureac,usbsw,utitl1,utitl2,uvfgex,uxcat,
     $    uxkgex,uxmod,uxopex,uxopt,uzveci,vreac,xgersi,xhfgex,ximaxi,
     $    xistti,xlkffg,xlkgex,xlkmod,xvfgex,zgexj,zvclgi)
          if (qrderr) go to 990
c
        else
c
          write (nttyo,1200) uoldv
 1200     format(/' * Error - (XCON6/xcon6) Coding to implement',
     $    /7x,'reading an input file in "W" format has not been',
     $    /7x,'implemented for version level "',a3,'."')
          go to 990
c
        endif
c
c       Try to recover the code and version data embedded in comments
c       on the old input file.
c
        call rd6wvz(ninpt,ucode,urelno,ustage,ueqlrn,ueqlst)
c
      elseif (uoldf(1:1) .eq. 'D') then
c
        if (uoldv(1:3).eq.'8.0' .and. .not.q8bchk) then
c
c         Distinguish 8.0 from 8.0 beta. 8.0 will include the
c         "Time print interval" and the "Log time print interval",
c         which 8.0 beta will not. Search for these only between
c         the "Temperature option" line and the "Steps print interval"
c         line.
c
          icount = 0
  120     read(ninpts,1030,end=140) uline
          j2 = index(uline(2:80),'|') - 1
          i = index(uline(2:j2),'Temperature option (jtemp)')
          if (i .le. 0) go to 120
c
  130     read(ninpts,1030,end=140) uline
          j2 = index(uline(2:80),'|') - 1
          i = index(uline(2:j2),'Time print interval')
          if (i .gt. 0) icount = icount + 1
          i = index(uline(2:j2),'Log time print interval')
          if (i .gt. 0) icount = icount + 1
          i = index(uline(2:j2),'Steps print interval')
          if (i .gt. 0) go to 140
          if (icount .lt. 2) go to 130
  140     q8beta = icount.lt.2
          q8bchk = .true.
          rewind(ninpts)
        endif
c
        if (uoldv(1:3) .eq. '6.0') then
c
          write (nttyo,1205)
 1205     format(/' * Error - (XCON6/xcon6) There is no "D" format',
     $    /7x,'for version level "6.0", hence',"can't read an input",
     $    /7x,'file in this format for this version level.')
          go to 990
c
        elseif (uoldv(1:3) .eq. '7.0') then
c
          call rd6d7(cdac,cesrb,csigma,dlzmx1,dlzmx2,dlzidp,
     $    dzprlg,dzprnt,eact,electr,fk,iact,iktbt,iktmax,imchmx,
     $    imech,iodb,iopg,iopr,iopt,ioscan,itermx,jcode,jreac,jtemp,
     $    jxmod,kct,kdim,kmax,kmt,kprs,ksq,ksppmx,kstpmx,kxmod,kxt,
     $    hact,modr,moffg,morr,mprs,mteaqb,mteb,nctmax,ndact,ndctmx,
     $    nesrbt,nffg,nffgmx,ninpts,nmodl1,nmodl2,nodbmx,nopgmx,
     $    noprmx,noptmx,nordlm,npslmx,nprmn,nprmx,nprsmx,nrct,nrctmx,
     $    nrk,nsk,nsrtmx,nsscmx,nsslmx,ntitl1,ntitl2,ntitmx,ntrymx,
     $    nttkmx,nttyo,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrtmx,
     $    qend,qrderr,rk0,rxbarb,sscrew,sk,tempci,tempc0,timemx,tolbt,
     $    toldl,tolsat,tolsst,tolx,trk0,tstrt,ttk,udac,uelemb,uendb,
     $    uesrb,uffg,undms,unrms,uprs,ureac,utitl1,utitl2,uxct16,uxmd24,
     $    uxopex,uxopt,vreac,xlkffg,xlkmod,zimax,zistrt,zkfac,zklogl,
     $    zklogu,zvclgi)
          if (qrderr) go to 990
c
        elseif (uoldv(1:3) .eq. '7.2') then
c
c         Note: version level '7.2" is not identical to version level
c         '7.0' for this format. However, the line parsing capability
c         used to handle this format allows this subroutine to read an
c         input file in this format for either version level.
c
          call rd6d7(cdac,cesrb,csigma,dlzmx1,dlzmx2,dlzidp,
     $    dzprlg,dzprnt,eact,electr,fk,iact,iktbt,iktmax,imchmx,
     $    imech,iodb,iopg,iopr,iopt,ioscan,itermx,jcode,jreac,jtemp,
     $    jxmod,kct,kdim,kmax,kmt,kprs,ksq,ksppmx,kstpmx,kxmod,kxt,
     $    hact,modr,moffg,morr,mprs,mteaqb,mteb,nctmax,ndact,ndctmx,
     $    nesrbt,nffg,nffgmx,ninpts,nmodl1,nmodl2,nodbmx,nopgmx,
     $    noprmx,noptmx,nordlm,npslmx,nprmn,nprmx,nprsmx,nrct,nrctmx,
     $    nrk,nsk,nsrtmx,nsscmx,nsslmx,ntitl1,ntitl2,ntitmx,ntrymx,
     $    nttkmx,nttyo,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrtmx,
     $    qend,qrderr,rk0,rxbarb,sscrew,sk,tempci,tempc0,timemx,tolbt,
     $    toldl,tolsat,tolsst,tolx,trk0,tstrt,ttk,udac,uelemb,uendb,
     $    uesrb,uffg,undms,unrms,uprs,ureac,utitl1,utitl2,uxct16,uxmd24,
     $    uxopex,uxopt,vreac,xlkffg,xlkmod,zimax,zistrt,zkfac,zklogl,
     $    zklogu,zvclgi)
          if (qrderr) go to 990
c
        elseif (uoldv(1:3).eq.'8.0' .and. q8beta) then
c
          call rd6d8b(cbsri,cdac,cesri,cgexj,csigma,dlxdmp,dlxmx0,
     $    dlxpll,dlxplo,dlxprl,dlxprn,eact,egersi,electr,fkrc,iact,
     $    ibsrti,iesrti,ietmax,iktmax,imchmx,imech,iodb,iopg,iopr,iopt,
     $    igerti,itermx,ixrti,jcode,jgerti,jetmax,jflgi,jgext,jpress,
     $    jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,ksppmx,kstpmx,
     $    kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,mtbaqi,mtbi,
     $    mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,nert,nertmx,net,
     $    netmax,nffg,nffgmx,ngexrt,ninpts,nobswt,nodbmx,nopgmx,noprmx,
     $    noptmx,noutpt,nprob,nprpmx,nprpti,nprsmx,nprsti,nptkmx,nrct,
     $    nrctmx,nrk,nsbswt,nsk,nsrt,nsrtmx,ntitl1,ntitl2,ntitmx,ntrymx,
     $    nttkmx,nttyo,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrt,
     $    nxrtmx,pressb,pressi,ptk,qend,qrderr,rkb,rxbari,sfcar,ssfcar,
     $    tempcb,tempci,tgexp,timmxi,tistti,tolbt,toldl,tolsat,tolxsf,
     $    trkb,ttk,ubmtbi,ubsri,ucxri,udac,uesri,uffg,ugerji,ugermo,
     $    ugersi,ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,uprphi,uprspi,
     $    ureac,usbsw,utitl1,utitl2,uvfgex,uxcat,uxkgex,uxmod,uxopex,
     $    uxopt,uzveci,vreac,xgersi,xhfgex,ximaxi,xistti,xlkffg,xlkgex,
     $    xlkmod,xvfgex,zgexj,zvclgi)
          if (qrderr) go to 990
c
        elseif (uoldv(1:3) .eq. '8.0') then
c
          call rd6d8(awmaxi,awmini,cbsri,cdac,cesri,cgexj,
     $    csigma,dlaplo,dlaprn,dleplo,dleprn,dlhplo,dlhprn,dloplo,
     $    dloprn,dltpll,dltplo,dltprl,dltprn,dlxdmp,dlxmx0,dlxpll,
     $    dlxplo,dlxprl,dlxprn,eact,egersi,ehmaxi,ehmini,electr,fkrc,
     $    iact,ibsrti,iesrti,ietmax,iktmax,imchmx,imech,iodb,iopg,
     $    iopr,iopt,igerti,itermx,ixrti,jcode,jgerti,jetmax,jflgi,
     $    jgext,jpress,jreac,jtemp,kbt,kct,kdim,kmax,kmt,kprs,ksplmx,
     $    ksppmx,kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprphi,mprspi,
     $    mtbaqi,mtbi,mwtges,nbti,nbtmax,nbt1mx,nctmax,ndact,ndctmx,
     $    nert,nertmx,net,netmax,nffg,nffgmx,ngexrt,ninpts,nobswt,
     $    nodbmx,nopgmx,noprmx,noptmx,nordmx,noutpt,nprob,nprpmx,
     $    nprpti,nprsmx,nprsti,nptkmx,nrct,nrctmx,nrk,nsbswt,nsk,
     $    nsrt,nsrtmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,nttyo,nxmdmx,
     $    nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrt,nxrtmx,o2maxi,o2mini,
     $    phmaxi,phmini,pressb,pressi,ptk,qend,qgexsh,qrderr,rkb,
     $    rxbari,sfcar,ssfcar,tempcb,tempci,tgexp,timmxi,tistti,tolbt,
     $    toldl,tolsat,tolxsf,trkb,ttk,ubmtbi,ubsri,ucxri,udac,uesri,
     $    uffg,ugerji,ugermo,ugersi,ugexj,ugexmo,ugexp,ugexr,uhfgex,
     $    uobsw,uprphi,uprspi,ureac,usbsw,utitl1,utitl2,uvfgex,uxcat,
     $    uxkgex,uxmod,uxopex,uxopt,uzveci,vreac,xgersi,xhfgex,ximaxi,
     $    xistti,xlkffg,xlkgex,xlkmod,xvfgex,zgexj,zvclgi)
          if (qrderr) go to 990
c
        else
c
          write (nttyo,1210) uoldv
 1210     format(/' * Error - (XCON6/xcon6) Coding to implement',
     $    /7x,'reading an input file in "D" format has not been',
     $    /7x,'implemented for version level "',a3,'."')
          go to 990
c
        endif
c
c       Try to recover the code and version data embedded in comments
c       on the old input file.
c
        call rd6dvz(ninpt,ucode,urelno,ustage,ueqlrn,ueqlst)
c
      else
c
        write (nttyo,1110) uoldf
 1110   format(/' * Error - (XCON6/xcon6) Have unknown format',
     $  /7x,'specifier "',a1,'" for the old input file.')
        go to 990
c
      endif
c
      if (qend) then
        close(ninpt)
        close(ninpts,status='delete')
        close(newin)
        close(nxcon)
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make changes in the content of the input file problem.
c
c     Look for the string "Version level= " or Version number= " in the
c     title. This identifies the version level of the input file.
c     If present, change it. If not, insert a new line containing
c     this at the beginning of the title.
c
      do 200 n = 1,ntitl1
        j = index(utitl1(n),'Version level=')
        if (j .eq. 0) j = index(utitl1(n),'version level=')
        if (j .gt. 0) then
          utitl1(n)(j + 14:j + 14) = ' '
          utitl1(n)(j + 15:j + 17) = unewv
          utitl1(n)(j + 18:80) = ' '
          go to 220
        endif
        j = index(utitl1(n),'Version number=')
        if (j .eq. 0) j = index(utitl1(n),'version number=')
        if (j .gt. 0) then
          utitl1(n)(j + 8:j + 14) = 'level= '
          utitl1(n)(j + 15:j + 17) = unewv
          utitl1(n)(j + 18:j + 18) = ' '
          go to 220
        endif
  200 continue
c
      if ((ntitl1 + 1) .gt. ntitmx) then
        write (nttyo,1115) nprob,ntitmx
 1115   format(/" * Error - (XCON6/xcon6) Can't add a version",
     $  /7x,'level marker to the first input file title of,'
     $  /7x,'problem ',i2,' because this title already has the',
     $  /7x,'maximum length of ',i3,' lines.')
        go to 990
      endif
c
      do 210 n = ntitl1,1,-1
        utitl1(n + 1) = utitl1(n)
  210 continue
      ntitl1 = ntitl1 + 1
c
      utitl1(1)(1:40) = '                                        '
      utitl1(1)(41:80) = '                                        '
      utitl1(1)(1:15)  = 'Version level= '
      utitl1(1)(16:18) = unewv
      utitl1(1)(21:27) = '(XCON6)'
c
  220 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get version levels in floating point format.
c
c     Calling sequence substitutions:
c       uoldv for ustr
c       verold for var
c
      call chreal(nttyo,qrderr,uoldv,verold)
c
c     Calling sequence substitutions:
c       unewv for ustr
c       vernew for var
c
      call chreal(nttyo,qrderr,unewv,vernew)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (verold.lt.8.0 .and. vernew.lt.8.0) then
c
        mtaqsm = 0.
        ncbt = kct
        do ncb = 1,ncbt
          mtaqsm = mtaqsm + mteaqb(ncb)
        enddo
c
        if (mtaqsm .le. 0.) then
          do ncb = 1,ncbt
            mteaqb(ncb) = mteb(ncb)
          enddo
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (verold.lt.8.0 .and. vernew.ge.8.0) then
c
c       Make translations from pre-Version 8.0 structures to Version 8.0
c       and up structures. All pre-Version 8.0 basis switching is mapped
c       to ordinary basis switching in Version 8.0 and up (none is
c       mapped to special basis switching).
c
        nerr = 0
c
        if (jtemp .le. 0) then
          if (ttk(1).eq.0. .and. ttk(2).eq.0. .and. ttk(3).eq.0.) then
c
c           Constant temperature.
c
            jtemp = 0
          elseif (ttk(2).eq.0. .and. ttk(3).eq.0.) then
            if (iopt(1) .le. 0) then
c
c             Linear tracking in Xi.
c
              jtemp = 1
            else
c
c             Linear tracking in time.
c
              jtemp = 2
            endif
          else
            write (nttyo,1505)
 1505       format(/" * Error - (XCON6/xcon6) Can't translate this",
     $      /7x,'pre-version level 8.0 input file to version 8.0 or',
     $      /7x,'higher because it uses a power series of order',
     $      /7x,'greater than one to track the temperature. At',
     $      /7x,'version level 8.0 or higher, only first order',
     $      /7x,'(linear) tracking is currently supported. Change',
     $      /7x,'the option on the old input file to a linear form',
     $      /7x,'by setting both ttk(2) and ttk(3) to zero.')
            nerr = nerr + 1
          endif
        else
c
c         Fluid mixing tracking.
c
          jtemp = 3
        endif
c
        jpress = 0
        pressb = 0.
        do n = 1,nptkmx
          ptk(n) = 0
        enddo
c
        if (tempci .le. 100.) then
          preshi = apresh(1,1)
        else
          preshi = 0.
          do nn = 1,5
            n = 6 - nn
            preshi = apresh(n,2) + tempci*preshi
          enddo
        endif
        pressi = preshi
c
        net = 0
        qgexsh = .false.
c
        iopt(17) = iopt(3)
        iopt(3) = iopt(2)
        iopt(2) = iopt(1)
        iopt(15) = iopt(11)
        iopt(11) = 0
        iopt(12) = iopt(7)
        iopt(7) = 0
        iopt(18) = iopt(13)
        iopt(9) = iopt(6)
        iopt(6) = iopt(5)
        if (iopt(6) .ge. 2) iopt(6) = 1
        iopt(5) = 0
c
        iopr(6) = iopr(11)
        iopr(11) = 0
        iopr(2) = iopr(3)
        iopr(3) = 0
        if (iopr(8) .ge. 1) then
          iopr(8) = 0
        else
          iopr(8) = -1
        endif
c
        iodb(7) = iodb(5)
        iodb(5) = iodb(3)
        iodb(3) = iodb(2)
        iodb(2) = iodb(9)
        iodb(9) = 0
c
        if (nmodl1 .ge. 3) then
          iopt(1) = 2
        elseif (nmodl1 .eq. 1) then
          iopt(1) = 1
        else
          iopt(1) = 0
        endif
c
        if (nmodl2 .le. 0) then
          iopt(13) = 0
        elseif (nmodl2 .eq. 1) then
          iopt(13) = 1
        else
          iopt(13) = 2
        endif
c
        nsbswt = 0
        nobswt = 0
c
        tempcb = tempc0
        tistti = tstrt
        timmxi = timemx
        xistti = zistrt
        ximaxi = zimax
        dlxprn = dzprnt
        dlxprl = dzprlg
        dlxplo = dzplot
        dlxpll = dzpllg
        dlxmx0 = dlzmx1
        dlxdmp = dlzidp
c
        dltprn = 1.e+38
        dltprl = 1.e+38
        dltplo = 1.e+38
        dltpll = 1.e+38
c
        dlaprn = 1.e+38
        dleprn = 1.e+38
        dlhprn = 1.e+38
        dloprn = 1.e+38
        dlaplo = 1.e+38
        dleplo = 1.e+38
        dlhplo = 1.e+38
        dloplo = 1.e+38
c
        awmini = -1.e+38
        awmaxi = 1.e+38
        ehmini = -1.e+38
        ehmaxi = 1.e+38
        o2mini = -1.e+38
        o2maxi = 1.e+38
        phmini = -1.e+38
        phmaxi = 1.e+38
c
        nordmx = nordlm
        if (nordmx .le. 0) nordmx = 6
c
        tolxsf = tolx
c
        do n = 1,nxopt
          if (uxopt(n)(1:1) .eq. ' ') uxopt(n) = 'None'
          if (uxopt(n)(1:5) .eq. 'none ') uxopt(n) = 'None'
          if (uxopt(n)(1:4) .eq. 'all ') uxopt(n) = 'All'
          if (uxopt(n)(1:7) .eq. 'alwith ') uxopt(n) = 'Allwith'
c
          if (uxct16(n)(1:5) .eq. 'none ') then
            uxcat(n) = 'None'
          else
c
c           Matching up other uxct16 strings is complicated
c           by the fact that the original inputs were converted
c           to lower case by the read subroutines.
c
            uxcat(n)(1:16) = uxct16(n)(1:16)
            do nn = 1,nctmax
              if (uxct16(n)(1:8) .eq. uelnlc(nn)(1:8)) then
                uxcat(n) = uelspn(nn)
                go to 225
              endif
            enddo
  225       continue
          endif
c
          if (unewf(1:1) .eq. 'W') then
            if (uxopt(n)(1:7) .ne. 'Allwith') then
              if (uxcat(n)(1:5) .eq. 'None ') uxcat(n) = ' '
            endif
          endif
        enddo
c
        do n = 1,nxopex
          if (uxopex(n)(1:1) .eq. ' ') uxopex(n) = 'None'
          if (uxopex(n)(1:5) .eq. 'none ') uxopex(n) = 'None'
        enddo
c
        do nrc = 1,nrct
          if (nsk(nrc) .eq. 0) then
            sfcar(nrc) = sk(nrc)
          elseif (nsk(nrc) .eq. 1) then
            j2 = ilnobl(ureac(nrc))
            write (nttyo,1510) ureac(nrc)(1:j2)
 1510       format(/" * Error - (XCON6/xcon6) Can't translate this",
     $      /7x,'pre-version level 8.0 input file to version 8.0 or',
     $      /7x,'higher because it has a reactant surface area flag',
     $      /7x,'(nsk) of 1 for the reactant ',a,'. This',
     $      /7x,'option fixes the specific surface area (cm2/g) of',
     $      /7x,'the reactant for both types of input files, but the',
     $      /7x,'corresponding input data for pre-version level 8.0',
     $      /7x,'is the total surface area (sk) in cm2, whereas that',
     $      /7x,'in version 8.0 or higher is the specific surface',
     $      /7x,"area in cm2/g. XCON6 can't calculate the latter from",
     $      /7x,"the former because it doesn't have the molecular",
     $      /7x,'weight of the reactant. Change the nsk option on',
     $      /7x,'the old input file to 0, convert the file, then',
     $      /7x,'change the surface area flag and data on the new',
     $      /7x,'file to the correct values.')
            nerr = nerr + 1
          else
            j2 = ilnobl(ureac(nrc))
            write (nttyo,1514) nsk(nrc),ureac(nrc)(1:j2)
 1514       format(/" * Error - (XCON6/xcon6) Can't translate this",
     $      /7x,'pre-version level 8.0 input file to version 8.0 or',
     $      /7x,'higher because it has an unknown reactant surface',
     $      /7x,'area flag (nsk) of ',i2,' for the reactant ',a,'.',
     $      /7x,'The only legal values are 0 and 1.')
            nerr = nerr + 1
          endif
        enddo
c
        do nrc = 1,nrct
          if (nrk(1,nrc) .ge. 4) then
            j2 = ilnobl(ureac(nrc))
            write (nttyo,1542) ureac(nrc)(1:j2),nrk(1,nrc)
 1542       format(/" * Error - (XCON6/xcon6) Can't translate this",
     $      /7x,'pre-version level 8.0 input file to version 8.0 or',
     $      /7x,'higher because the reactant ',a,' has a forward',
     $      /7x,'direction rate law code (nrk(1,nrc)) value of ',i3,'.',
     $      /7x,'This was defined at pre-version 8.0 levels as the',
     $      /7x,'activity-term rate law (e.g., the Plummer et al.',
     $      /7x,'(1978) rate law for carbonate mineral dissolution',
     $      /7x,"and precipitation). This isn't supported at version",
     $      /7x,'level 8.0 or higher because it is difficult to',
     $      /7x,'guarantee consistency between the rate constant',
     $      /7x,'data and the thermodynamic data. Use a TST-like',
     $      /7x,'formulation (nrk(1,nrc) = 2) instead.')
            nerr = nerr + 1
          endif
          if (nrk(2,nrc) .ge. 4) then
            j2 = ilnobl(ureac(nrc))
            write (nttyo,1544) ureac(nrc)(1:j2),nrk(2,nrc)
 1544       format(/" * Error - (XCON6/xcon6) Can't translate this",
     $      /7x,'pre-version level 8.0 input file to version 8.0 or',
     $      /7x,'higher because the reactant ',a,' has a backward',
     $      /7x,'direction rate law code (nrk(2,nrc)) value of ',i3,'.',
     $      /7x,'This was defined at pre-version 8.0 levels as the',
     $      /7x,'activity-term rate law (e.g., the Plummer et al.',
     $      /7x,'(1978) rate law for carbonate mineral dissolution',
     $      /7x,"and precipitation). This isn't supported at version",
     $      /7x,'level 8.0 or higher because it is difficult to',
     $      /7x,'guarantee consistency between the rate constant',
     $      /7x,'data and the thermodynamic data. Use a TST-like',
     $      /7x,'formulation (nrk(2,nrc) = 2) instead.')
            nerr = nerr + 1
          endif
        enddo
c
        call copyaa(fk,fkrc,nrctmx)
c
        nmax = imchmx*2*nrctmx
        call copyaa(rk0,rkb,nmax)
        call copyaa(trk0,trkb,nmax)
c
        call copyia(nesrbt,iesrti,nsrtmx)
c
        nmax = nctmax*nsrtmx
        call copyaa(cesrb,cesri,nmax)
        call copyca(uesrb,uesri,nmax)
c
        call copyia(iktbt,ixrti,nxrtmx)
c
        nmax = iktmax*nxrtmx
        call copyaa(rxbarb,rxbari,nmax)
        call copyca(uendb,ucxri,nmax)
c
c       Set the number of generic ion exchanger reactants.
c       There is no provision for these prior to Version 8.0,
c       so here this number must be zero.
c
        nert = 0
c
c       Count the number of special reactants.
c
        nsrt = 0
        do nrc = 1,nrct
          if (jcode(nrc) .eq. 2) then
            nsrt = nsrt + 1
          endif
        enddo
c
        if (nsrt .gt. 0) then
c
c         Create a reaction for each special reactant.
c
          nsr = 0
          do nrc = 1,nrct
            if (jcode(nrc) .eq. 2) then
              nsr = nsr + 1
              cbsri(1,nsr) = -1.0
              ubsri(1,nsr) = ureac(nrc)
              n = 1
c
              ncbo = 0
              do ncb = 1,nesrbt(nsr)
                if (uesrb(ncb,nsr)(1:2).eq.'O ' .or.
     $            uesrb(ncb,nsr)(1:2).eq.'o ') then
                  ncbo = ncb
                  go to 242
                endif
              enddo
  242         continue
c
              ncbh = 0
              do ncb = 1,nesrbt(nsr)
                if (uesrb(ncb,nsr)(1:2) .eq. 'H ' .or.
     $            uesrb(ncb,nsr)(1:2).eq.'h ') then
                  ncbh = ncb
                  go to 320
                endif
              enddo
  320         continue
c
              ctxo = 0.
              ctxh = 0.
              if (ncbo .gt. 0) ctxo = -cesri(ncbo,nsr)
              if (ncbh .gt. 0) ctxh = -cesri(ncbh,nsr)
              ztx = 0.
c
c             Map all chemical elements except H and O.
c
              do ncb = 1,nesrbt(nsr)
                if (ncb.ne.ncbo .and. ncb.ne.ncbh) then
                  do nce = 1,nelspt
                    if (uesrb(ncb,nsr)(1:8) .eq. uelnam(nce)(1:8)) then
                      n = n + 1
                      cx = cesri(ncb,nsr)/celspe(nce)
                      cbsri(n,nsr) = cx
                      ubsri(n,nsr) = uelspn(nce)
                      ctxo = ctxo + celspo(nce)*cx
                      ctxh = ctxh + celsph(nce)*cx
                      ztx = ztx + zelsp(nce)*cx
                      go to 243
                    endif
                  enddo
                  j2 = ilnobl(uesrb(ncb,nsr))
                  j3 = ilnobl(ureac(nrc))
                  write (nttyo,1126) uesrb(ncb,nsr)(1:j2),
     $            ureac(nrc)(1:j3)
 1126             format(/" * Error - (XCON6/xcon6) Can't map the",
     $            /7x,'chemical element ',a,', which appears in the'
     $            /7x,'composition for the special reactant ',a,',',
     $            /7x,'into a corresponding basis species to use in',
     $            /7x,'composing reaction for that reactant.')
                  nerr = nerr + 1
  243             continue
                endif
              enddo
c
c             Get H+ from charge balance.
c
              if (abs(ztx) .gt. 1.e-12) then
                n = n + 1
                cx = -ztx
                cbsri(n,nsr) = cx
                ubsri(n,nsr) = 'H+ '
                ztx = ztx + cx
                ctxh = ctxh + cx
              endif
c
c             Get H2O from H balance.
c
              if (abs(ctxh) .gt. 1.e-12) then
                n = n + 1
                cx = -0.5*ctxh
                cbsri(n,nsr) = cx
                ubsri(n,nsr) = 'H2O '
                ctxo = ctxo + cx
                ctxh = ctxh + 2.*cx
              endif
c
c             Get O2(g) from O balance.
c
              if (abs(ctxo) .gt. 1.e-12) then
                n = n + 1
                cx = -0.5*ctxo
                cbsri(n,nsr) = cx
                ubsri(n,nsr) = 'O2(g)'
                ctxo = ctxo + 2.*cx
              endif
c
c             Check the balances on H, O, and charge.
c
              if (abs(ctxh).gt.1.e-12 .or. abs(ctxo).gt.1.e-12 .or.
     $          abs(ztx).gt.1.e-12) then
                j2 = ilnobl(ureac(nrc))
                write (nttyo,1128) ureac(nrc)(1:j2)
 1128           format(/' * Error - (XCON6/xcon6) The reaction which',
     $          /7x,'was composed for the special reactant "',a,'"',
     $          /7x,'fails to satisfy all balance conditions.')
                nerr = nerr + 1
              endif
c
              ibsrti(nsr) = n
            endif
          enddo
        endif
c
        do n = 1,nxmod
          uxmod(n)(1:24) = uxmd24(n)
        enddo
c
        kbt = ksq
        km1 = kbt + 1
        kx1 = kmt + 1
        nbti = ksq
        nsbswt = 0
c
        do krow = 1,ksq
          nbb = krow
          ubmtbi(nbb)(1:24) = unrms(krow)
          ubmtbi(nbb)(25:48) = uaqsln(1:24)
        enddo
c
        nobsw = 0
        do krow = 1,ksq
          if (undms(krow)(1:24) .ne. unrms(krow)(1:24)) then
            nobsw = nobsw + 1
            uobsw(1,nobsw)(1:24) = undms(krow)(1:24)
            uobsw(1,nobsw)(25:48) = uaqsln(1:24)
            uobsw(2,nobsw)(1:24) = unrms(krow)(1:24)
            uobsw(2,nobsw)(25:48) = uaqsln(1:24)
          endif
        enddo
        nobswt = nobsw
c
        ncbt = kct
        ko2 = kct + 1
        no2 = ko2
c
        ncbo = 0
        do ncb = 1,ncbt
          if (uelemb(ncb)(1:2).eq.'O ' .or.
     $      uelemb(ncb)(1:2).eq.'o ') then
            ncbo = ncb
            go to 246
          endif
        enddo
  246   continue
c
        ncbh = 0
        do ncb = 1,ncbt
          if (uelemb(ncb)(1:2).eq.'H ' .or.
     $      uelemb(ncb)(1:2).eq.'h ') then
            ncbh = ncb
            go to 252
          endif
        enddo
  252   continue
c
        nbbo = ncbo
        nbbh = ncbh
c
        do ncb = 1,ncbt
          do nce = 1,nelspt
            if (uelemb(ncb)(1:8) .eq. uelnam(nce)(1:8)) then
              nelptr(ncb) = nce
              go to 247
            endif
          enddo
          j2 = ilnobl(uelemb(ncb))
          write (nttyo,1720) uelemb(ncb)(1:j2)
 1720     format(/" * Error - (XCON6/xcon6) Can't map the mass",
     $    /7x,'balance total for chemical element ',a,' to a',
     $    /7x,'mass balance total for the corresponding strict',
     $    /7x,"basis species because it doesn't match any element",
     $    /7x,'listed in the pseudo-data file in x6sbsp.h.')
          nerr = nerr + 1
  247     continue
        enddo
c
        do ncb = 1,ncbt
          if (ncb.ne.ncbo .and. ncb.ne.ncbh) then
            nbb = ncb
            nce = nelptr(ncb)
            mtbi(nbb) = mteb(ncb)/celspe(nce)
            mtbaqi(nbb) = mteaqb(ncb)/celspe(nce)
          endif
        enddo
c
        mtaqsm = 0.
        do ncb = 1,ncbt
          mtaqsm = mtaqsm + mteaqb(ncb)
        enddo
c
        if (mtaqsm .le. 0.) then
          do ncb = 1,ncbt
            if (ncb.ne.ncbo .and. ncb.ne.ncbh) then
              nbb = ncb
              mtbaqi(nbb) = mtbi(nbb)
            endif
          enddo
        endif
c
        yvec(1) = mteb(ncbo)
        yvec(2) = mteb(ncbh)
        yvec(3) = electr
        do ncb = 1,ncbt
          if (ncb.ne.ncbo .and. ncb.ne.ncbh) then
            nce = nelptr(ncb)
            mxc = mteb(ncb)/celspe(nce)
            yvec(1) = yvec(1) - celspo(nce)*mxc
            yvec(2) = yvec(2) - celsph(nce)*mxc
            yvec(3) = yvec(3) - zelsp(nce)*mxc
          endif
        enddo
c
        aamatr(1,1) = celspe(nceo)
        aamatr(1,2) = celspo(nceh)
        aamatr(1,3) = 2.0
        aamatr(2,1) = celsph(nceo)
        aamatr(2,2) = celspe(nceh)
        aamatr(2,3) = 0.0
        aamatr(3,1) = zelsp(nceo)
        aamatr(3,2) = zelsp(nceh)
        aamatr(3,3) = 0.0
        qpr = .false.
c
c       Calling sequence substitutions:
c         xvec for delvec
c         3 for kdim
c         imamax for kmax
c         0 for noutpt
c         yvec for rhsvec
c
        call msolvr(aamatr,xvec,gmmatr,ier,ipivot,3,imamax,
     $  0,nttyo,qpr,yvec)
c
        mtbi(nbbo) = xvec(1)
        mtbi(nbbh) = xvec(2)
        mtbi(no2) = xvec(3)
c
        if (mtaqsm .gt. 0.) then
          yvec(1) = mteaqb(ncbo)
          yvec(2) = mteaqb(ncbh)
          yvec(3) = electr
          do ncb = 1,ncbt
            if (ncb.ne.ncbo .and. ncb.ne.ncbh) then
              nce = nelptr(ncb)
              mxc = mteaqb(ncb)/celspe(nce)
              yvec(1) = yvec(1) - celspo(nce)*mxc
              yvec(2) = yvec(2) - celsph(nce)*mxc
              yvec(3) = yvec(3) - zelsp(nce)*mxc
            endif
          enddo
c
c         Calling sequence substitutions:
c           xvec for delvec
c           3 for kdim
c           imamax for kmax
c           0 for noutpt
c           yvec for rhsvec
c
          call msolvr(aamatr,xvec,gmmatr,ier,ipivot,3,imamax,
     $    0,nttyo,qpr,yvec)
c
          mtbaqi(nbbo) = xvec(1)
          mtbaqi(nbbh) = xvec(2)
          mtbaqi(no2) = xvec(3)
        else
          mtbaqi(nbbo) = mtbi(nbbo)
          mtbaqi(nbbh) = mtbi(nbbh)
          mtbaqi(no2) = mtbi(no2)
        endif
c
        do ncb = 1,ncbt
          nbb = ncb
          nce = nelptr(ncb)
          ubmtbi(nbb)(1:24) = uelspn(nce)
          ubmtbi(nbb)(25:48) = uaqsln
        enddo
        ubmtbi(no2)(1:24) = 'O2(g) '
        ubmtbi(no2)(25:48) = uaqsln
c
        do krow = 1,ksq
          uzveci(krow)(1:24) = unrms(krow)(1:24)
          uzveci(krow)(25:48) = uaqsln(1:24)
        enddo
c
        do krow = km1,kmt
          ux24 = undms(krow)
          if (undms(krow)(1:4) .eq. 'fix ') then
            ux24(1:5) = 'fix_f'
            ux24(6:24) = undms(krow)(9:24)
          endif
          uzveci(krow)(1:24) = ux24
          uzveci(krow)(25:48) = ux24
        enddo
c
        do krow = kx1,kxt
          uzveci(krow)(1:24) = unrms(krow)(1:24)
          uzveci(krow)(25:48) = undms(krow)(1:24)
        enddo
c
        nprpti = 0
        nprsti = 0
        if (kprs .gt. 0) then
          do n = 1,nprmn
            nprpti = nprpti + 1
            nprsti = nprsti + 1
            uprphi(nprpti) = uprs(n)(1:24)
            uprspi(nprsti)(1:24) = uprs(n)(1:24)
            uprspi(nprsti)(25:48) = uprs(n)(1:24)
            mprphi(nprpti) = mprs(n)
            mprspi(nprsti) = mprs(n)
          enddo
          unamph(1:24) = ' '
          do n = 1,nprmx
            if (uprs(n)(25:48) .ne. unamph(1:24)) then
              nprpti = nprpti + 1
              uprphi(nprpti) = uprs(n)(25:48)
              unamph(1:24) = uprs(n)(25:48)
            endif
            nprsti = nprsti + 1
            uprphi(nprpti) = uprs(n)(1:24)
            uprspi(nprsti)(1:24) = uprs(n)(1:24)
            uprspi(nprsti)(25:48) = uprs(n)(1:24)
            mprspi(nprsti) = mprs(n)
            mprphi(nprpti) = mprphi(nprpti) + mprs(n)
          enddo
        endif
c
        if (nerr .gt. 0) stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if ((abs(verold - 8.0).le.epstst .and. q8beta)
     $ .and. vernew.ge.8.0) then
c
c       Make additions if converting from Version 8.0 beta to
c       Version 8.0 and up.
c
        dltprn = 1.e+38
        dltprl = 1.e+38
        dltplo = 1.e+38
        dltpll = 1.e+38
c
        dlaprn = 1.e+38
        dleprn = 1.e+38
        dlhprn = 1.e+38
        dloprn = 1.e+38
        dlaplo = 1.e+38
        dleplo = 1.e+38
        dlhplo = 1.e+38
        dloplo = 1.e+38
c
        awmini = -1.e+38
        awmaxi = 1.e+38
        ehmini = -1.e+38
        ehmaxi = 1.e+38
        o2mini = -1.e+38
        o2maxi = 1.e+38
        phmini = -1.e+38
        phmaxi = 1.e+38
c
        nordmx = 6
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if ((abs(verold - 8.0).le.epstst .and. q8beta)
     $ .and. vernew.lt.8.0) then
c
c       Make additions if converting from Version 8.0 beta to
c       pre-Version 8.0.
c
        nordlm = 6
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (verold.ge.8.0 .and. vernew.lt.8.0) then
c
c       Make translations from Version 8.0 and up structures to
c       pre-Version 8.0 structures, and vice versa. Special basis
c       switching in Version 8.0 and up can't be mapped to any
c       pre-Version 8.0 structure.
c
        nerr = 0
c
        if (jtemp .eq. 0) then
c
c         Constant temperature.
c
           jtemp = 0
        elseif (jtemp .eq. 1) then
c
c         Linear tracking in Xi.
c
          jtemp = 0
          if (iopt(2) .gt. 0) then
            write (nttyo,1522)
 1522       format(/" * Error - (XCON6/xcon6) Can't translate this",
     $      /7x,'version level 8.0 or higher input file to a lower',
     $      /7x,'version level because the temperature is set to',
     $      /7x,'track as a linear function of Xi (jtemp = 1), while',
     $      /7x,'kinetic mode is selected (iopt(2) = 1). This',
     $      /7x,"combination can't be mapped to the lower version",
     $      /7x,'level, because temperature tracking in time instead',
     $      /7x,'of Xi is forced by selection of kinetic mode.')
            nerr = nerr + 1
          endif
        elseif (jtemp .eq. 2) then
c
c         Linear tracking in time.
c
          jtemp = 0
          if (iopt(2) .le. 0) then
            write (nttyo,1524)
 1524       format(/" * Error - (XCON6/xcon6) Can't translate this",
     $      /7x,'version level 8.0 or higher input file to a lower',
     $      /7x,'version level because the temperature is set to',
     $      /7x,'track as a linear function of time (jtemp = 2),',
     $      /7x,'while reaction progress mode is selected (iopt(2)',
     $      ' = 0).',/7x,"This combination can't be mapped to the",
     $      ' lower version',/7x,'level, because temperature tracking',
     $      ' in Xi instead',/7x,'of time is forced by selection of',
     $      ' reaction progress mode.')
            nerr = nerr + 1
          endif
        elseif (jtemp .eq. 3) then
c
c         Fluid mixing tracking.
c
          jtemp = 1
        else
          write (nttyo,1526) jtemp
 1526     format(/" * Error - (XCON6/xcon6) Can't translate this",
     $    /7x,'version level 8.0 or higher input file to a lower',
     $    /7x,'version level because the temperature tracing flag',
     $    /7x,'(jtemp) has an unkown value of ',i2,'.')
          nerr = nerr + 1
        endif
c
        if (jpress .ge. 1) then
          qcntmp = ttk(1).eq.0. .and. ttk(2).eq.0. .and. ttk(3).eq. 0.
          qcnpre = ptk(1).eq.0. .and. ptk(2).eq.0. .and. ptk(3).eq. 0.
          if (qcntmp .and. qcnpre) then
            if (tempcb .le. 100.) then
              presh = apresh(1,1)
            else
              presh = 0.
              do nn = 1,5
                n = 6 - nn
                presh = apresh(n,2) + tempcb*presh
              enddo
            endif
            if (abs(pressb - presh) .gt. 1.e-4) then
              write (nttyo,1700) pressb,presh
 1700         format(/" * Error - (XCON6/xcon6) Can't translate this",
     $        /7x,'version level 8.0 or higher input file to a lower',
     $        /7x,'version level because it specifies a pressure of',
     $        /7x,f9.4,' bars, which differs from the pre-Version 8',
     $        /7x,'standard grid pressure of ',f9.4,' bars. To convert',
     $        /7x,'this input file, set jpress = 0, which sets the',
     $        /7x,'pressure equal to the standard grid pressure.')
              nerr = nerr + 1
            endif
          else
              write (nttyo,1705) jpress
 1705         format(/" * Error - (XCON6/xcon6) Can't translate this",
     $        /7x,'version level 8.0 or higher input file to a lower',
     $        /7x,'version level because it specifies a pressure',
     $        /7x,'tracking option (jpress = ',i2,') which differs',
     $        /7x,'from the pre-Version 8 requirement that the',
     $        ' pressure',/7x,'be equal to the standard grid pressure',
     $        ' at any',/7x,'temperature. To convert this input file,',
     $        ' set jpress = 0.')
              nerr = nerr + 1
            endif
        endif
c
        if (net .gt. 0) then
          write (nttyo,1707) net
 1707     format(/" * Error - (XCON6/xcon6) Can't translate this",
     $    /7x,'version level 8.0 or higher input file to a lower',
     $    /7x,'version level because it defines ',i2,' generic',
     $    /7x,'ion exchangers. Pre-version 8 input files lack the',
     $    /7x,'generic ion exchange capability.')
          nerr = nerr + 1
        endif
c
        if (iopt(1) .eq. 2) then
          nmodl1 = 3
        elseif (iopt(1) .eq. 1) then
          nmodl1 = 1
        else
          nmodl1 = 2
        endif
c
        if (iopt(13) .le. 0) then
          nmodl2 = 0
        elseif (iopt(13) .eq. 1) then
          nmodl2 = 1
        else
          nmodl2 = 2
        endif
c
        iopt(1) = iopt(2)
        iopt(2) = iopt(3)
        iopt(3) = iopt(17)
        iopt(11) = iopt(15)
        iopt(15) = 0
        iopt(7) = iopt(12)
        iopt(12) = 0
        iopt(13) = iopt(18)
        iopt(18) = 0
        iopt(5) = iopt(6)
        if (iopt(5) .eq. 1) iopt(5) = 2
        iopt(6) = iopt(9)
        iopt(9) = 0
c
        if (iopt(20) .gt. 0) then
          write (nttyo,1520)
 1520     format(/" * Warning - (XCON6/xcon6) Can't translate the",
     $    /7x,'"Advanced EQ6 PICKUP File Options" (iopt(20)) option',
     $    /7x,'specified on this version level 8.0 or higher input',
     $    /7x,'file to a lower version level because this option does',
     $    /7x,'not exist at the lower version level.')
        endif
c
        iopr(11) = iopr(6)
        iopr(6) = 0
        iopr(3) = iopr(2)
        iopr(2) = 0
        if (iopr(8) .ge. 0) then
          iopr(8) = 1
        else
          iopr(8) = 0
        endif
c
        iodb(9) = iodb(2)
        iodb(2) = iodb(3)
        iodb(3) = iodb(5)
        iodb(5) = iodb(7)
        iodb(7) = 0
c
        tempc0 = tempcb
        tstrt = tistti
        timemx = timmxi
        zistrt = xistti
        zimax = ximaxi
        dzprnt = dlxprn
        dzprlg = dlxprl
        dzplot = dlxplo
        dzpllg = dlxpll
        dlzmx1 = dlxmx0
        cplim = 0.
        dlzidp = dlxdmp
        ioscan = 0
c
        tolx = tolxsf
c
        do n = 1,nxopt
          if (uxopt(n)(1:5) .eq. 'None ') uxopt(n) = ' '
          if (uxopt(n)(1:5) .eq. 'none ') uxopt(n) = ' '
          if (uxopt(n)(1:4) .eq. 'All ') uxopt(n) = 'all'
          if (uxopt(n)(1:7) .eq. 'Alwith ') uxopt(n) = 'alwith'
          if (uxopt(n)(1:8) .eq. 'Allwith ') uxopt(n) = 'alwith'
c
          if (uxcat(n)(1:5) .eq. 'none ') then
            uxct16(n) = ' '
          elseif (uxcat(n)(1:5) .eq. 'None ') then
            uxct16(n) = ' '
          else
            uxct16(n)(1:16) = uxcat(n)(1:16)
            do nn = 1,nctmax
              if (uxcat(n)(1:24) .eq. uelspn(nn)(1:24)) then
                if (vernew .ge. 7.0) then
                  uxct16(n) = uelnam(nn)
                else
                  uxct16(n) = uelnlc(nn)
                endif
                go to 257
              endif
            enddo
  257       continue
          endif
        enddo
c
        do n = 1,nxopex
          if (uxopex(n)(1:5) .eq. 'none ') uxopex(n) = ' '
          if (uxopex(n)(1:5) .eq. 'None ') uxopex(n) = ' '
        enddo
c
c       Count the number of generic ion exchanger reactants.
c       There is no provision for these prior to Version 8.0,
c       so it is not possible to translate a Version 8.0 or
c       greater input file to a pre-Version 8.0 format.
c
        nert = 0
        do nrc = 1,nrct
          if (jcode(nrc) .eq. 5) then
            nert = nert + 1
          endif
        enddo
c
        if (nert .gt. 0) then
          write (nttyo,1709) nert
 1709     format(/" * Error - (XCON6/xcon6) Can't translate this",
     $    /7x,'version level 8.0 or higher input file to a lower',
     $    /7x,'version level because it defines ',i2,' generic',
     $    /7x,'ion exchanger reactants. Pre-version 8 input files',
     $    /7x,'lack the generic ion exchange capability.')
          nerr = nerr + 1
        endif
c
        do nrc = 1,nrct
c
          if (nsk(nrc) .eq. 0) then
            sk(nrc) = sfcar(nrc)
          elseif (nsk(nrc) .eq. 1) then
            j2 = ilnobl(ureac(nrc))
            write (nttyo,1530) ureac(nrc)(1:j2)
 1530       format(/" * Error - (XCON6/xcon6) Can't translate this",
     $      /7x,'version level 8.0 or higher input file to pre-version',
     $      /7x,' 8.0 level because it has a reactant surface area',
     $      /7x,'flag (nsk) of 1 for the reactant ',a,'. This',
     $      /7x,'option fixes the specific surface area (cm2/g) of',
     $      /7x,'the reactant for both types of input files, but the',
     $      /7x,'corresponding input data for version level 8.0 or',
     $      /7x,'higher is the specific surface area in cm2/g,',
     $      /7x,'whereas that in pre-version 8.0 is the total surface'
     $      /7x,"area in cm2. XCON6 can't calculate the latter from",
     $      /7x,"the former because it doesn't have the molecular",
     $      /7x,'weight of the reactant. Change the nsk option on the',
     $      /7x,'old input file to 0, convert the file, then change',
     $      /7x,'the surface area flag and data on the new file to',
     $      /7x,'the correct values.')
            nerr = nerr + 1
          elseif (nsk(nrc) .eq. 2) then
            j2 = ilnobl(ureac(nrc))
            write (nttyo,1535) ureac(nrc)(1:j2)
 1535       format(/" * Error - (XCON6/xcon6) Can't translate this",
     $      /7x,'version level 8.0 or higher input file to pre-version',
     $      /7x,' 8.0 level because it has a reactant surface area',
     $      /7x,'flag (nsk) of 2 for the reactant ',a,'. This',
     $      /7x,'option provides for geometric growth or reduction',
     $      /7x,'in surface area for constant particle number. This',
     $      /7x,'option does not exist for versions earlier than',
     $      /7x,'version 8.0. Change the nsk option on the old input',
     $      /7x,'file to 0, convert the file, then change the surface',
     $      /7x,'area flag and data on the new file to the desired',
     $      ' values.')
            nerr = nerr + 1
          else
            j2 = ilnobl(ureac(nrc))
            write (nttyo,1540) nsk(nrc),ureac(nrc)(1:j2)
 1540       format(/" * Error - (XCON6/xcon6) Can't translate this",
     $      /7x,'version level 8.0 or higher input file to a lower',
     $      /7x,'version level because it has a reactant surface',
     $      /7x,'area flag (nsk) with a value of ',i2,' for the',
     $      /7x,'reactant ',a,'. This corresponds to an option that',
     $      /7x,"doesn't exist for pre-version level 8.0 files.",
     $      /7x,'The only legal values are such files are 0 and 1.')
            nerr = nerr + 1
          endif
        enddo
c
        call copyaa(fkrc,fk,nrctmx)
c
        nmax = imchmx*2*nrctmx
        call copyaa(rkb,rk0,nmax)
        call copyaa(trkb,trk0,nmax)
c
        call copyia(iesrti,nesrbt,nsrtmx)
c
        nmax = nctmax*nsrtmx
        call copyaa(cesri,cesrb,nmax)
        call copyca(uesri,uesrb,nmax)
c
        call copyia(ixrti,iktbt,nxrtmx)
c
        nmax = iktmax*nxrtmx
        call copyaa(rxbari,rxbarb,nmax)
        call copyca(ucxri,uendb,nmax)
c
        do 250 n = 1,nxmod
c
          uxmd24(n) = uxmod(n)(1:24)
c
c         Set the jxmod value. This is a flag for the type of species:
c
c           0 = aqueous species
c           1 = pure mineral
c           2 = gas
c           3 = solid solution
c
c         The jxmod array is not used in Version 8.0 and up. It is
c         not possible to design a perfect logic for determining the
c         correct value. Here there may be a problem in determining
c         whether a species is a pure mineral or a solid solution.
c         Also, in some cases, the Version 8.0 and up option may not
c         map back to the pre-Version 8.0 format.
c
          if (uxmod(n)(25:48) .eq. uaqsln(1:24)) then
            jxmod(n) = 0
            go to 250
          endif
c
          if (uxmod(n)(25:28) .eq. 'Gas ') then
            jxmod(n) = 2
            go to 250
          endif
c
          if (uxmod(n)(25:48).ne.uxmod(n)(1:24) .and.
     $      uxmod(n)(25:27).ne.'   ') then
            nerr = nerr + 1
            j2 = ilnobl(uxmod(n))
            write (nttyo,1120) uxmod(n)(1:j2)
 1120       format(/" * Error - (XCON6/xcon6) Can't translate a",
     $      ' suppress/alter option',/7x,'for "',a,'"',
     $      /7x,'to pre-Version 8.0 format.')
            go to 250
          endif
c
          i = index(uxmod(n),'(aq)')
          if (i .gt. 0) then
            jxmod(n) = 0
            go to 250
          endif
c
          ux24 = uxmod(n)(1:24)
          j2 = ilnobl(ux24)
          if (j2 .gt. 0) then
            if (ux24(j2:j2).eq.'+' .or. ux24(j2:j2).eq.'-') then
              jxmod(n) = 0
              go to 250
            endif
          endif
c
          i = index(uxmod(n),'(g)')
          if (i .le. 0) i = index(uxmod(n),'Hydrogen')
          if (i .le. 0) i = index(uxmod(n),'Oxygen')
          if (i .le. 0) i = index(uxmod(n),'Nitrogen')
          if (i .le. 0) i = index(uxmod(n),'Steam')
          if (i .gt. 0) then
            jxmod(n) = 2
            go to 250
          endif
c
          i = index(uxmod(n),'(ss)')
          if (i .le. 0) i = index(uxmod(n),'-ss')
          if (i .le. 0) i = index(uxmod(n),'Carbonate-Calcite')
          if (i .le. 0) i = index(uxmod(n),'Biotite')
          if (i .le. 0) i = index(uxmod(n),'Olivine')
          if (i .le. 0) i = index(uxmod(n),'Plagioclase')
          if (i .le. 0) i = index(uxmod(n),'Orthopyroxene')
          if (i .le. 0) i = index(uxmod(n),'Smectite-di')
          if (i .le. 0) i = index(uxmod(n),'Saponite-tri')
          if (i .gt. 0) then
            jxmod(n) = 3
            go to 250
          endif
c
          i = index(uxmod(n),'ite ')
          if (i .le. 0) i = index(uxmod(n),'ime ')
          if (i .le. 0) i = index(uxmod(n),'ine ')
          if (i .le. 0) i = index(uxmod(n),'Quartz')
          if (i .gt. 0) then
            jxmod(n) = 1
            go to 250
          endif
c
          j2 = ilnobl(uxmod(n))
          write (nttyo,1119) uxmod(n)(1:j2)
 1119     format(/" * Warning - (XCON6/xcon6) Can't unambiguously",
     $    /7x,'determine what kind of species in an alter/suppress',
     $    /7x,'option is "',a,'".',
     $    /7x,'Setting jxmod to 1. Check to see that this is correct',
     $    /7x,'(0 = aqueous species, 1 = pure mineral,',
     $    ' 3 = gas species',
     $    /7x,'3 = solid solution).')
          jxmod(n) = 1
  250  continue
c
        ksq = kbt
        km1 = ksq + 1
        kx1 = kmt + 1
c
        if (nsbswt .gt. 0) then
          write (nttyo,1710) nsbswt
 1710     format(/" * Error - (XCON6/xcon6) Can't translate this",
     $    /7x,'version level 8.0 or higher input file to a lower',
     $    /7x,'version level because it contains ',i3,' directives',
     $    /7x'for special basis switching.')
          nerr = nerr + 1
        endif
c
        ncbt = kct
        ko2 = kct + 1
        no2 = ko2
c
        nbbo = 0
        do nbb = 1,nbti
          if (ubmtbi(nbb)(1:24) .eq. uelspn(nceo)(1:24)) then
            nbbo = nbb
            go to 248
          endif
        enddo
  248   continue
c
        nbbh = 0
        do nbb = 1,nbti
          if (ubmtbi(nbb)(1:24) .eq. uelspn(nceh)(1:24)) then
            nbbh = nbb
            go to 249
          endif
        enddo
  249   continue
c
        ncbo = nbbo
        ncbh = nbbh
c
        do nbb = 1,nbti
          if (nbb .eq. no2) go to 253
          ncb = nbb
          do nce = 1,nelspt
            if (ubmtbi(nbb)(1:24) .eq. uelspn(nce)(1:24)) then
              nelptr(ncb) = nce
              go to 253
            endif
          enddo
          j2 = ilnobl(ubmtbi(nbb)(1:24))
          write (nttyo,1730) ubmtbi(nbb)(1:j2)
 1730     format(/" * Error - (XCON6/xcon6) Can't map the mass",
     $    /7x,'balance total for the strict basis species',
     $    /7x,a,' to a mass balance total for the corresponding',
     $    /7x,"chemical element because it doesn't match any strict",
     $    /7x,'basis species listed in the pseudo-data file in',
     $    ' x6sbsp.h.')
          nerr = nerr + 1
  253     continue
        enddo
c
        mteb(ncbo) = celspe(nceo)*mtbi(nbbo) + 2.*mtbi(no2)
        mteaqb(ncbo) = celspe(nceo)*mtbaqi(nbbo) + 2.*mtbaqi(no2)
        do nbb = 1,ncbt
          if (nbb.ne.nbbo .and. nbb.ne.nbbh) then
            nce = nelptr(nbb)
            mteb(ncbo) = mteb(ncbo) + celspo(nce)*mtbi(nbb)
            mteaqb(ncbo) = mteaqb(ncbo) + celspo(nce)*mtbaqi(nbb)
          endif
        enddo
c
        mteb(ncbh) = celspe(nceh)*mtbi(nbbh) + 2.*mtbi(nbbo)
        mteaqb(ncbh) = celspe(nceh)*mtbaqi(nbbh) + 2.*mtbaqi(nbbo)
        do nbb = 1,ncbt
          if (nbb.ne.nbbo .and. nbb.ne.nbbh) then
            nce = nelptr(nbb)
            mteb(ncbh) = mteb(ncbh) + celsph(nce)*mtbi(nbb)
            mteaqb(ncbh) = mteaqb(ncbh) + celsph(nce)*mtbaqi(nbb)
          endif
        enddo
c
        do ncb = 1,ncbt
          if (ncb.ne.ncbo .and. ncb.ne.ncbh) then
            nbb = ncb
            nce = nelptr(nbb)
            mteb(ncb) = mteb(ncb) + celspe(nce)*mtbi(nbb)
            mteaqb(ncb) = mteaqb(ncb) + celspe(nce)*mtbaqi(nbb)
          endif
        enddo
c
        electc = 0.
        do nbb = 1,ncbt
          nce = nelptr(nbb)
          electc = electc + zelsp(nce)*mtbaqi(nbb)
        enddo
c
        if (abs(electc - electr) .gt. 1.e-12) then
          write (nttyo,1735) electc,electr
 1735     format(/' * Warning - (XCON6/xcon6) The calculated',
     $    ' electrical balance',/7x,"doesn't match the input",
     $    ' electrical imbalance:',//9x,'Calculated= ',1pe22.15,
     $    /14x,'Input= ',e22.15)
        endif
c
        do ncb = 1,ncbt
          nbb = ncb
          nce = nelptr(ncb)
          uelemb(ncb) = uelnam(nce)
        enddo
c
        do krow = 1,ksq
          unrms(krow)(1:24) = uzveci(krow)(1:24)
          undms(krow)(1:24) = uzveci(krow)(1:24)
          do nobsw = 1,nobswt
            if (uobsw(1,nobsw)(1:24) .eq. undms(krow)(1:24)) then
              unrms(krow)(1:24) = uobsw(2,nobsw)(1:24)
              go to 255
            endif
          enddo
  255     continue
        enddo
c
        do krow = km1,kmt
          unrms(krow)(1:24) = ' '
          ux24 = uzveci(krow)(1:24)
          if (uzveci(krow)(1:5) .eq. 'fix_f') then
            ux24(1:8) = 'fix     '
            ux24(9:24) = uzveci(krow)(6:21)
          endif
          undms(krow)(1:24) = ux24
        enddo
c
        do krow = kx1,kxt
          unrms(krow)(1:24) = uzveci(krow)(25:48)
          undms(krow)(1:24) = uzveci(krow)(1:24)
        enddo
c
        nprmn = 0
        nprmx = 0
        if (kprs .gt. 0) then
          do n = 1,nprsti
            uprs(n)(1:48) = uprspi(n)(1:48)
            mprs(n) = mprspi(n)
            if (uprs(n)(1:24) .eq. uprs(n)(25:48)) then
              nprmn = nprmn + 1
            else
              nprmx = nprmx + 1
            endif
          enddo
        endif
c
        if (nerr .gt. 0) stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Trap incompatibilities between version level 6.0 and higher
c     version levels.
c
      if (abs(verold - 6.0).le.epstst  .and.
     $  vernew .gt. 6.0) then
c
c       Version level 6.0 contains a rate law code of 5 which does
c       not map to subsequent version levels.
c
        nerr = 0
        do nrc = 1,nrct
          if (nrk(1,nrc) .ge. 5) then
            j2 = ilnobl(ureac(nrc))
            write (nttyo,1117) ureac(nrc)(1:j2),nrk(1,nrc)
 1117       format(/' * Error - (XCON6/xcon6) Reactant ',a,
     $      /7x,'has a forward direction rate law code value',
     $      /7x,'of ',i3,'. This was defined at version level',
     $      /7x,"6.0. It can't be mapped to version level 7.0",
     $      /7x,'or higher.')
            nerr = nerr + 1
          endif
          if (nrk(2,nrc) .ge. 5) then
            j2 = ilnobl(ureac(nrc))
            write (nttyo,1118) ureac(nrc)(1:j2),nrk(2,nrc)
 1118       format(/' * Error - (XCON6/xcon6) Reactant ',a,
     $      /7x,'has a backward direction rate law code value',
     $      /7x,'of ',i3,'. This was defined at version level',
     $      /7x,"6.0. It can't be mapped to version level 7.0",
     $      /7x,'or higher.')
            nerr = nerr + 1
          endif
        enddo
        if (nerr .gt. 0) stop
      endif
c
      if (verold.gt.6.0 .and. abs(vernew - 6.0).le.epstst) then
c
c       Version level 6.0 has no treatment of the temperature dependence
c       of rate constants.
c
        nerr = 0
        do nrc = 1,nrct
          do i = 1,imech(1,nrc)
            if (iact(i,1,nrc).ne.0 .or. trk0(i,1,nrc).ne.0. .or.
     $        eact(i,1,nrc).ne.0. .or. hact(i,1,nrc).ne.0.) then
              j2 = ilnobl(ureac(nrc))
              write (nttyo,1122) ureac(nrc)(1:j2),i
 1122         format(/' * Error - (XCON6/xcon6) Reactant ',a,
     $        /7x,'has a temperature-dependent rate constant in',
     $        /7x,'term ',i2," of the forward rate law. This can't",
     $        /7x,"be mapped to version level 6.0.")
              nerr = nerr + 1
            endif
          enddo
          do i = 1,imech(2,nrc)
            if (iact(i,2,nrc).ne.0 .or. trk0(i,2,nrc).ne.0. .or.
     $        eact(i,2,nrc).ne.0. .or. hact(i,2,nrc).ne.0.) then
              j2 = ilnobl(ureac(nrc))
              write (nttyo,1124) ureac(nrc)(1:j2),i
 1124         format(/' * Error - (XCON6/xcon6) Reactant ',a,
     $        /7x,'has a temperature-dependent rate constant in',
     $        /7x,'term ',i2," of the backward rate law. This can't",
     $        /7x,'be mapped to version level 6.0.')
              nerr = nerr + 1
            endif
          enddo
        enddo
        if (nerr .gt. 0) stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the current problem on the new input file.
c
      if (unewf(1:1) .eq. 'W') then
        if (unewv(1:3) .eq. '6.0') then
c
          call wr6w6(cdac,cesrb,cplim,csigma,dlzmx1,dlzmx2,
     $    dlzidp,dzpllg,dzplot,dzprlg,dzprnt,electr,fk,ifile,iktbt,
     $    iktmax,imchmx,imech,iodb,iopg,iopr,iopt,ioscan,itermx,jcode,
     $    jreac,jtemp,jxmod,kct,kdim,kmax,kmt,kprs,ksq,ksplmx,ksppmx,
     $    kstpmx,kxmod,kxt,modr,moffg,morr,mprs,mteaqb,mteb,nctmax,
     $    ndact,ndctmx,nesrbt,newin,nffg,nffgmx,nmodl1,nmodl2,nodbmx,
     $    nopgmx,noprmx,noptmx,nordlm,npslmx,nprmn,nprmx,nprsmx,nrct,
     $    nrctmx,nrk,nsk,nsrtmx,nsscmx,nsslmx,ntitl1,ntitl2,ntitmx,
     $    ntrymx,nttkmx,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrtmx,
     $    rk0,rxbarb,sscrew,sk,tempci,tempc0,timemx,tolbt,toldl,tolsat,
     $    tolsst,tolx,tstrt,ttk,uacion,ucode,udac,uelemb,uendb,ueqlrn,
     $    ueqlst,uesrb,uffg,undms,unrms,uprs,ureac,urelno,ustage,
     $    utitl1,utitl2,uxct16,uxmd24,uxopex,uxopt,vreac,xlkffg,xlkmod,
     $    zimax,zistrt,zkfac,zklogl,zklogu,zvclgi)
c
        elseif (unewv(1:3) .eq. '7.0') then
c
          call wr6w7(cdac,cesrb,cplim,csigma,dlzmx1,dlzmx2,dlzidp,
     $    dzpllg,dzplot,dzprlg,dzprnt,eact,electr,fk,iact,ifile,iktbt,
     $    iktmax,imchmx,imech,iodb,iopg,iopr,iopt,ioscan,itermx,jcode,
     $    jreac,jtemp,jxmod,kct,kdim,kmax,kmt,kprs,ksq,ksplmx,ksppmx,
     $    kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprs,mteaqb,mteb,nctmax,
     $    ndact,ndctmx,nesrbt,newin,nffg,nffgmx,nmodl1,nmodl2,nodbmx,
     $    nopgmx,noprmx,noptmx,nordlm,npslmx,nprmn,nprmx,nprsmx,nrct,
     $    nrctmx,nrk,nsk,nsrtmx,nsscmx,nsslmx,ntitl1,ntitl2,ntitmx,
     $    ntrymx,nttkmx,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrtmx,
     $    rk0,rxbarb,sscrew,sk,tempci,tempc0,timemx,tolbt,toldl,tolsat,
     $    tolsst,tolx,trk0,tstrt,ttk,ucode,udac,uelemb,uendb,ueqlrn,
     $    ueqlst,uesrb,uffg,undms,unrms,uprs,ureac,urelno,ustage,utitl1,
     $    utitl2,uxct16,uxmd24,uxopex,uxopt,vreac,xlkffg,xlkmod,zimax,
     $    zistrt,zkfac,zklogl,zklogu,zvclgi)
c
        elseif (unewv(1:3) .eq. '7.2') then
c
          call wr6w72(cdac,cesrb,cplim,csigma,dlzmx1,dlzmx2,dlzidp,
     $    dzpllg,dzplot,dzprlg,dzprnt,eact,electr,fk,iact,ifile,iktbt,
     $    iktmax,imchmx,imech,iodb,iopg,iopr,iopt,ioscan,itermx,jcode,
     $    jreac,jtemp,jxmod,kct,kdim,kmax,kmt,kprs,ksq,ksplmx,ksppmx,
     $    kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprs,mteaqb,mteb,nctmax,
     $    ndact,ndctmx,nesrbt,newin,nffg,nffgmx,nmodl1,nmodl2,nodbmx,
     $    nopgmx,noprmx,noptmx,nordlm,npslmx,nprmn,nprmx,nprsmx,nrct,
     $    nrctmx,nrk,nsk,nsrtmx,nsscmx,nsslmx,ntitl1,ntitl2,ntitmx,
     $    ntrymx,nttkmx,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrtmx,
     $    rk0,rxbarb,sscrew,sk,tempci,tempc0,timemx,tolbt,toldl,tolsat,
     $    tolsst,tolx,trk0,tstrt,ttk,ucode,udac,uelemb,uendb,ueqlrn,
     $    ueqlst,uesrb,uffg,undms,unrms,uprs,ureac,urelno,ustage,utitl1,
     $    utitl2,uxct16,uxmd24,uxopex,uxopt,vreac,xlkffg,xlkmod,zimax,
     $    zistrt,zkfac,zklogl,zklogu,zvclgi)
c
        elseif (unewv(1:3) .eq. '8.0') then
c
          call wr6w8(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,
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
c
        else
          write (nttyo,1220) unewv
 1220     format(/' * Error - (XCON6/xcon6) Coding to implement',
     $    /7x,'writing an input file in "W" format has not been',
     $    /7x,'implemented for version level "',a3,'."')
          go to 990
c
        endif
      elseif (unewf(1:1) .eq. 'D') then
        if (unewv(1:3) .eq. '6.0') then
c
          write (nttyo,1230)
 1230     format(/' * Error - (XCON6/xcon6) There is no "D" format',
     $    /7x,'for version level "6.0", hence',"can't write an input",
     $    /7x,'file in this format for this version level.')
          go to 990
c
        elseif (unewv(1:3) .eq. '7.0') then
c
          call wr6d7(cdac,cesrb,csigma,dlzmx1,dlzmx2,dlzidp,
     $    dzprlg,dzprnt,eact,electr,fk,iact,iktbt,iktmax,imchmx,
     $    imech,iodb,iopg,iopr,iopt,ioscan,itermx,jcode,jreac,jtemp,
     $    jxmod,kct,kdim,kmax,kmt,kprs,ksq,ksppmx,kstpmx,kxmod,kxt,
     $    hact,modr,moffg,morr,mprs,mteaqb,mteb,nctmax,ndact,ndctmx,
     $    nesrbt,newin,nffg,nffgmx,nmodl1,nmodl2,nodbmx,nopgmx,noprmx,
     $    noptmx,nordlm,npslmx,nprmn,nprmx,nprsmx,nrct,nrctmx,nrk,
     $    nsk,nsrtmx,nsscmx,nsslmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,
     $    nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrtmx,rk0,rxbarb,
     $    sscrew,sk,tempci,tempc0,timemx,tolbt,toldl,tolsat,tolsst,
     $    tolx,trk0,tstrt,ttk,ucode,udac,uelemb,uendb,uesrb,ueqlrn,
     $    ueqlst,uffg,undms,unrms,uprs,ureac,urelno,ustage,utitl1,
     $    utitl2,uxct16,uxmd24,uxopex,uxopt,vreac,xlkffg,xlkmod,zimax,
     $    zistrt,zkfac,zklogl,zklogu,zvclgi)
c
        elseif (unewv(1:3) .eq. '7.2') then
c
          call wr6d72(cdac,cesrb,csigma,dlzmx1,dlzmx2,dlzidp,
     $    dzprlg,dzprnt,eact,electr,fk,iact,iktbt,iktmax,imchmx,
     $    imech,iodb,iopg,iopr,iopt,ioscan,itermx,jcode,jreac,jtemp,
     $    jxmod,kct,kdim,kmax,kmt,kprs,ksq,ksppmx,kstpmx,kxmod,kxt,
     $    hact,modr,moffg,morr,mprs,mteaqb,mteb,nctmax,ndact,ndctmx,
     $    nesrbt,newin,nffg,nffgmx,nmodl1,nmodl2,nodbmx,nopgmx,noprmx,
     $    noptmx,nordlm,npslmx,nprmn,nprmx,nprsmx,nrct,nrctmx,nrk,
     $    nsk,nsrtmx,nsscmx,nsslmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,
     $    nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrtmx,rk0,rxbarb,
     $    sscrew,sk,tempci,tempc0,timemx,tolbt,toldl,tolsat,tolsst,
     $    tolx,trk0,tstrt,ttk,ucode,udac,uelemb,uendb,uesrb,ueqlrn,
     $    ueqlst,uffg,undms,unrms,uprs,ureac,urelno,ustage,utitl1,
     $    utitl2,uxct16,uxmd24,uxopex,uxopt,vreac,xlkffg,xlkmod,zimax,
     $    zistrt,zkfac,zklogl,zklogu,zvclgi)
c
        elseif (unewv(1:3) .eq. '8.0') then
c
          call wr6d8(awmaxi,awmini,cbsri,cdac,cesri,cgexj,csigma,
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
c
        else
c
          write (nttyo,1240) unewv
 1240     format(/' * Error - (XCON6/xcon6) Coding to implement',
     $    /7x,'writing an input file in "D" format has not been',
     $    /7x,'implemented for version level "',a3,'."')
          go to 990
c
        endif
      else
c
        write (nttyo,1130) unewf
 1130   format(/' * Error - (XCON6/xcon6) Have unknown format',
     $  /7x,'specifier "',a1,'" for the new input file.')
        go to 990
c
      endif
c
c     Go back to do another problem on the same input file.
c
      go to 100
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  990 write (nttyo,1140) uoldv
 1140 format(/' * Error - (XCON6/xcon6) Have encountered a read',
     $ /7x,'error while reading the old input file. First check the',
     $ /7x,'version level. The version level was taken to be "',a3,'".',
     $ /7x,'If the version level is marked incorrectly in the primary',
     $ /7x,'title of the first problem, correct this, as the default',
     $ /7x,'value read from the IXCON file or set in the code itself',
     $ /7x,'does not override this. If the problem is not due to an',
     $ /7x,'error in the assumed version level, there is probably a',
     $ /7x,'format error on the old input file. Check the stripped',
     $ /7x,'copy of the old input file (INPUTS).')
c
c     Close all files. Don't delete the inputs file. Delete the NEWIN
c     file, which is not valid.
c
      close(ninpt)
      close(ninpts)
      close(newin,status='delete')
      close(nxcon)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
