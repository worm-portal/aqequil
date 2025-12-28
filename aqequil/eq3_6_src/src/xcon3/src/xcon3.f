      program xcon3
c
c     This is the main program of the XCON3 code. Configuration
c     identification, the copyright statement, legal disclaimers,
c     and similar statements are contained in XCON3/aaaxc3.f, the
c     lead-off subroutine in the XCON3 source code. A short description
c     of this program is also contained in that subroutine.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      include 'xcon3/x3op7.h'
c
      include 'eqlib/eqlj8.h'
      include 'eqlib/eqlo8.h'
c
c-----------------------------------------------------------------------
c
      integer ietpar,iktpar,jetpar,nbtpar,netpar,nodbpa,nopgpa,noprpa,
     $ noptpa,nsqpar,ntitpa,nxicpa,nxmdpa,nxtipa,nxtpar
c
      parameter(ietpar = 10,iktpar = 20,jetpar = 3,nbtpar = 250)
      parameter(netpar = 12,ntitpa = 100,nxmdpa = 40,nxtpar = 50)
      parameter(nodbpa = 20,nopgpa = 20,noprpa = 20,noptpa = 20)
c
      parameter(nsqpar = nbtpar,nxtipa = nxtpar,nxicpa = iktpar*nxtpar)
c
      integer ietmax,iktmax,jetmax,nbtmax,netmax,nodbmx,nopgmx,noprmx,
     $ noptmx,nsqmax,ntitmx,nxicmx,nxmdmx,nxtimx,nxtmax
c
c-----------------------------------------------------------------------
c
      integer newin,ninpt,ninpts,noutpt,nttyo
c
      integer iodb(nodbpa),iopg(nopgpa),iopr(noprpa),iopt(noptpa),
     $ jflagb(nsqpar),jflgi(nbtpar),jgext(netpar),jgexti(netpar),
     $ jxmod(nxmdpa),kxmod(nxmdpa),ncmpri(2,nxtipa),ncompb(nxtpar),
     $ ngexrt(jetpar,netpar),ngexti(jetpar,netpar)
c
      integer itermx,jpres3,net,neti,nprob,ntitl,nxcon,nxmod
      integer nsq,nxtb
      integer nbti,nobswt,nsbswt,nxti
c
      integer i,icount,iebal3,irdxc3,itdsf3,j,jfli,j2,ik,iktb,k,kl,n,
     $ nb,nerr,nmax,nn,nobsw,nr1,nr2,ns,ntest,nx,nxi,nxic
c
      integer ilnobl
c
      logical qend,qex,qgexsh,qrderr,qxcon,q8bchk,q8beta
c
      character*80 utitl(ntitpa)
      character*56 ugexr(ietpar,jetpar,netpar)
      character*48 uxmod(nxmdpa)
      character*48 ucospi(nbtpar),uobsw(2,nbtpar),usbsw(2,nbtpar),
     $ uspeci(nbtpar)
      character*24 ugexmo(netpar),ugexp(netpar),ugexpi(netpar),
     $ ugexsi(ietpar,jetpar,netpar),usoli(nxtipa),umemi(nxicpa)
      character*24 ubasis(nsqpar),umemb(iktpar,nxtpar),uphas1(nsqpar),
     $ uphas2(nsqpar),usolb(nxtpar),uspecb(nsqpar),uxmd24(nxmdpa)
      character*8 ugexj(jetpar,netpar),ugexji(jetpar,netpar),
     $ uhfgex(ietpar,jetpar,netpar),uvfgex(ietpar,jetpar,netpar),
     $ uxkgex(ietpar,jetpar,netpar)
c
      character*80 uline
      character*24 uacion,ublk24,uebal,uredox,ustr,ustr1,ux24
      character*8 uplatc,uplatm,ustelu,ustxc3,uveelu,uvexc3
      character*8 uv
      character*3 uoldv,uoldvd,unewv
      character*1 uoldf,unewf,ux1
c
      real*8 apresh(5,2)
      real*8 cspb(nsqpar),xbarb(iktpar,nxtpar)
      real*8 covali(nbtpar),xbari(nxicpa)
      real*8 cgexj(jetpar,netpar),cgexpi(netpar),
     $ egexsi(ietpar,jetpar,netpar),mwtges(netpar),
     $ tgexp(netpar),xgexsi(ietpar,jetpar,netpar),
     $ xhfgex(ietpar,jetpar,netpar),xlkgex(ietpar,jetpar,netpar),
     $ xvfgex(ietpar,jetpar,netpar),xlkmod(nxmdpa),zgexj(jetpar,netpar)
c
      real*8 ehi,epstst,fep,fo2lgi,pei,presh,press,rho,scamas,tempc,
     $ tdspkg,tdspl,tolbt,toldl,tolsat,tolspf
c
      real*8 verold,vernew
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
        external bkdxc3
c
c     END_MACHINE_DEPENDENT_CODE
c
c-----------------------------------------------------------------------
c
      data noutpt /0/
c
      data ublk24 /'                        '/
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
      data ninpt  /9/,ninpts /10/,newin /11/,nxcon /12/
c
c-----------------------------------------------------------------------
c
c     Get configuration identification data.
c
      call aaaxc3(ustxc3,uvexc3)
      call aaaelu(ustelu,uveelu)
      call platfd(uplatc,uplatm)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set dimensioning variables.
c
      ietmax = ietpar
      iktmax = iktpar
      jetmax = jetpar
      nbtmax = nbtpar
      netmax = netpar
      nodbmx = nodbpa
      nopgmx = nopgpa
      noprmx = noprpa
      noptmx = noptpa
      nsqmax = nsqpar
      ntitmx = ntitpa
      nxicmx = nxicpa
      nxmdmx = nxmdpa
      nxtimx = nxtipa
      nxtmax = nxtpar
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Does the INPUT file (the old input file) exist?
c
      inquire(file="input",exist=qex)
      if (.not.qex) then
        write (nttyo,1000)
 1000   format(' * Error- (XCON3/xcon3) The INPUT file (the old input',
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
 1010 format(" * Error - (XCON3/xcon3) Can't open the INPUT file (the",
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
      do n = 1,100
        read (ninpt,1060) ux1
 1060   format(a1)
        if (ux1(1:1) .ne. 'c') kl = kl + 1
        if (ux1(1:1) .eq. '|') k = k + 1
        if (kl .ge. ntest) then
          if (k .ge. ntest) uoldf = 'D'
          go to 810
        endif
      enddo
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
 1020 format(" * Error - (XCON3/xcon3) Can't open the INPUTS file",
     $ /7x,'(the stripped old input file).')
      close(ninpt)
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Strip comment lines from the input file. Also strip any blank
c     lines from a "D" format file.
c
  730 if (uoldf(1:1) .eq. 'W') then
        do n = 1,10000
          read (ninpt,1030,end=750) uline
 1030     format(a80)
          if (uline(1:1) .ne. '*') write (ninpts,1030) uline
        enddo
c
        write (nttyo,1040)
 1040   format(" * Error - (XCON3/xcon3) The old input file is too
     $  ' long.")
        close(ninpt)
        close(ninpts,status='delete')
        go to 999
c
  750   continue
      elseif (uoldf(1:1) .eq. 'D') then
        do n = 1,10000
          read (ninpt,1030,end=757) uline
          j2 = ilnobl(uline)
          if (j2 .gt. 0) then
            if (uline(1:1).ne.'*' .and. uline(1:1).ne.'c')
     $      write (ninpts,1030) uline
          endif
        enddo
c
        write (nttyo,1040)
        close(ninpt)
        close(ninpts,status='delete')
        go to 999
c
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
 1050 format(" * Error - (XCON3/xcon3) Can't open the NEWIN file",
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
c       '7.2' = version 7.2 (including versions 7.2a and 7.2b)
c       '8.0' = version 8.0
c
c         uoldv  = version level of the old file
c         uoldvd = default version level of the old file
c         unewv  = version level of the new file
c
      uoldv = '   '
      uoldvd = '   '
      unewv = '   '
c
c     Set ultimate default values.
c
      uoldvd = '7.0'
      unewv = '7.2'
      unewf = 'W'
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
 1055   format(/" * Error - (XCON3/xcon3) Can't open the IXCON options",
     $  ' file.')
        close(ninpt)
        close(ninpts,status='delete')
        close(newin)
        go to 999
c
c       Read the IXCON options file.
c
  790  call rddixc(nxcon,uoldvd,unewf,unewv)
      else
        write (nttyo,1057)
 1057   format(/' * Error - (XCON3/xcon3) The IXCON options file',
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
c       '7.2' = version 7.2 (including versions 7.2a and 7.2b)
c       '8.0' = version 8.0
c
      i = 0
      do n = 1,ntitmx + 4
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
      enddo
  825 continue
c
      i = 0
      rewind(ninpts)
      do n = 1,ntitmx
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
      enddo
  835 continue
c
      write (nttyo,1070)
 1070 format(/' * Warning - (XCON3/xcon3) The title on the old input',
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
 1080   format(/" * Warning - (XCON3/xcon3) Can't determine the",
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
        do n = 1,1000
          read (ninpts,1030,end=838) uline
          j = index(uline,'uacion= ')
          if (j .gt. 0) then
            uoldv = '6.0'
            go to 838
          endif
        enddo
  838   continue
      endif
c
c     If necessary, use the default for the version level of the old
c     input file.
c
      if (uoldv(1:3) .eq. '   ') then
        uoldv = uoldvd
        write (nttyo,1090) uoldvd
 1090   format(/' * Warning - (XCON3/xcon3) Taking the default',
     $  /7x,'value of "',a3,'" for the version level of the old',
     $  /7x,'input file.')
      endif
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
      call initcb(utitl,ntitmx)
c
      call initiz(iopt,noptmx)
      call initiz(iopg,nopgmx)
      call initiz(iopr,noprmx)
      call initiz(iodb,nodbmx)
c
      call initcb(uxmod,nxmdmx)
      call initcb(uxmd24,nxmdmx)
      call initiz(jxmod,nxmdmx)
      call initiz(kxmod,nxmdmx)
      call initaz(xlkmod,nxmdmx)
c
      call initcb(uspecb,nsqmax)
      call initcb(ubasis,nsqmax)
      call initiz(jflagb,nsqmax)
      call initaz(cspb,nsqmax)
      call initcb(uphas1,nsqmax)
      call initcb(uphas2,nsqmax)
c
      call initcb(usolb,nxtmax)
c
      nmax = iktmax*nxtmax
      call initcb(umemb,nmax)
      call initaz(xbarb,nmax)
c
      call initcb(uspeci,nbtmax)
      call initiz(jflgi,nbtmax)
      call initaz(covali,nbtmax)
      call initcb(ucospi,nbtmax)
c
      nmax = 2*nbtmax
      call initcb(usbsw,nmax)
      call initcb(uobsw,nmax)
c
      call initcb(usoli,nxtimx)
c
      call initcb(umemi,nxicmx)
      call initaz(xbari,nxicmx)
c
      nmax = 2*nxtimx
      call initiz(ncmpri,nmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The uacion variable only appears on version level '6.0' input
c     files. Provide a blank default value.
c
      uacion = ' '
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the current problem on the stripped input file.
c
      if (uoldf(1:1) .eq. 'W') then
        if (uoldv(1:3) .eq. '6.0') then
c
          call rd3w6(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,
     $    jflagb,jxmod,kxmod,ncompb,ninpts,nodbmx,nopgmx,noprmx,noptmx,
     $    nsq,nsqmax,ntitl,ntitmx,nttyo,nxmdmx,nxmod,nxtb,nxtmax,qend,
     $    qrderr,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,uacion,
     $    ubasis,uebal,umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,
     $    xbarb,uxmd24,xlkmod)
          if (qrderr) go to 990
c
        elseif (uoldv(1:3) .eq. '7.0') then
c
          call rd3w7(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,
     $    jflagb,jxmod,kxmod,ncompb,ninpts,nodbmx,nopgmx,noprmx,noptmx,
     $    nsq,nsqmax,ntitl,ntitmx,nttyo,nxmdmx,nxmod,nxtb,nxtmax,qend,
     $    qrderr,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,ubasis,
     $    uebal,umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,xbarb,
     $    uxmd24,xlkmod)
          if (qrderr) go to 990
c
        elseif (uoldv(1:3) .eq. '7.2') then
c
c         Note: version level 7.2 is identical to version level 7.0
c         for this format.
c
          call rd3w7(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,
     $    jflagb,jxmod,kxmod,ncompb,ninpts,nodbmx,nopgmx,noprmx,noptmx,
     $    nsq,nsqmax,ntitl,ntitmx,nttyo,nxmdmx,nxmod,nxtb,nxtmax,qend,
     $    qrderr,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,ubasis,
     $    uebal,umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,xbarb,
     $    uxmd24,xlkmod)
          if (qrderr) go to 990
c
        elseif (uoldv(1:3) .eq. '8.0') then
c
          call rd3w8(cgexj,cgexpi,covali,ehi,egexsi,fo2lgi,iebal3,
     $    ietmax,iodb,iopg,iopr,iopt,irdxc3,itdsf3,itermx,jgext,jetmax,
     $    jflgi,jgexti,jpres3,kxmod,mwtges,nbti,nbtmax,ncmpri,net,neti,
     $    netmax,ngexti,ninpts,ngexrt,nobswt,nodbmx,nopgmx,noprmx,
     $    noptmx,noutpt,nprob,nsbswt,ntitl,ntitmx,nttyo,nxicmx,nxmdmx,
     $    nxmod,nxti,nxtimx,pei,press,qend,qgexsh,qrderr,rho,scamas,
     $    tdspkg,tdspl,tempc,tgexp,tolbt,toldl,tolspf,ucospi,uebal,
     $    ugexj,ugexji,ugexmo,ugexp,ugexpi,ugexr,ugexsi,umemi,uobsw,
     $    uredox,usbsw,usoli,uspeci,utitl,uhfgex,uvfgex,uxkgex,uxmod,
     $    xbari,xgexsi,xhfgex,xlkgex,xvfgex,xlkmod,zgexj)
          if (qrderr) go to 990
c
        else
c
          write (nttyo,1200) uoldv
 1200     format(/' * Error - (XCON3/xcon3) Coding to implement',
     $    /7x,'reading an input file in "W" format has not been',
     $    /7x,'implemented for version level "',a3,'."')
          go to 990
c
        endif
      elseif (uoldf(1:1) .eq. 'D') then
c
        if (uoldv(1:3).eq.'8.0' .and. .not.q8bchk) then
c
c         Distinguish 8.0 from 8.0 beta. 8.0 will include the
c         "Advisory:" and the "Option: on further procesing" strings,
c         which 8.0 beta will not. Search for these only between
c         the "Default redox constraint (irdxc3):" line and the
c         "Alter/Suppress Options  | (nxmod)" line.
c
          icount = 0
  120     read(ninpts,1030,end=140) uline
          j2 = index(uline(2:80),'|') - 1
          i = index(uline(2:j2),'Default redox constraint (irdxc3):')                                            |
          if (i .le. 0) go to 120
c
  130     read(ninpts,1030,end=140) uline
          j2 = index(uline(2:80),'|') - 1
          i = index(uline(2:j2),'Advisory:')
          if (i .gt. 0) icount = icount + 1
          i = index(uline(2:j2),'Option: on further processing')
          if (i .gt. 0) icount = icount + 1
          i = index(uline(2:j2),'Alter/Suppress Options  | (nxmod)')
          if (i .gt. 0) go to 140
          if (icount .lt. 2) go to 130
  140     q8beta = icount.lt.2
          q8bchk = .true.
          rewind(ninpts)
        endif
c
        if (uoldv(1:3) .eq. '6.0') then
c
          write (nttyo,1205) uoldv
 1205     format(/' * Error - (XCON3/xcon3) There is no "D" format',
     $    /7x,'for version level "6.0", hence',"can't read an input",
     $    /7x,'file in this format for this version level.')
          go to 990
c
        elseif (uoldv(1:3) .eq. '7.0') then
c
          call rd3d7(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,
     $    jflagb,jxmod,kxmod,ncompb,ninpts,nodbmx,nopgmx,noprmx,noptmx,
     $    nsq,nsqmax,ntitl,ntitmx,nttyo,nxmdmx,nxmod,nxtb,nxtmax,qend,
     $    qrderr,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,ubasis,
     $    uebal,umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,xbarb,
     $    uxmd24,xlkmod)
          if (qrderr) go to 990
c
        elseif (uoldv(1:3) .eq. '7.2') then
c
c         Note: XCON3/rd3d7.f can read a "D" format input file at
c         version level '7.2' as well as at version level '7.0'.
c
          call rd3d7(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,
     $    jflagb,jxmod,kxmod,ncompb,ninpts,nodbmx,nopgmx,noprmx,noptmx,
     $    nsq,nsqmax,ntitl,ntitmx,nttyo,nxmdmx,nxmod,nxtb,nxtmax,qend,
     $    qrderr,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,ubasis,
     $    uebal,umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,xbarb,
     $    uxmd24,xlkmod)
          if (qrderr) go to 990
c
        elseif (uoldv(1:3).eq.'8.0' .and. q8beta) then
c
          call rd3d8b(cgexj,cgexpi,covali,ehi,egexsi,fo2lgi,iebal3,
     $    ietmax,iodb,iopg,iopr,iopt,irdxc3,itdsf3,itermx,jgext,jetmax,
     $    jflgi,jgexti,jpres3,kxmod,mwtges,nbti,nbtmax,ncmpri,net,neti,
     $    netmax,ngexti,ninpts,ngexrt,nobswt,nodbmx,nopgmx,noprmx,
     $    noptmx,noutpt,nprob,nsbswt,ntitl,ntitmx,nttyo,nxicmx,nxmdmx,
     $    nxmod,nxti,nxtimx,pei,press,qend,qrderr,rho,scamas,tdspkg,
     $    tdspl,tempc,tgexp,tolbt,toldl,tolspf,ucospi,uebal,ugexj,
     $    ugexji,ugexmo,ugexp,ugexpi,ugexr,ugexsi,umemi,uobsw,uredox,
     $    usbsw,usoli,uspeci,utitl,uhfgex,uvfgex,uxkgex,uxmod,xbari,
     $    xgexsi,xhfgex,xlkgex,xvfgex,xlkmod,zgexj)
          if (qrderr) go to 990
c
        elseif (uoldv(1:3) .eq. '8.0') then
c
          call rd3d8(cgexj,cgexpi,covali,ehi,egexsi,fo2lgi,iebal3,
     $    ietmax,iodb,iopg,iopr,iopt,irdxc3,itdsf3,itermx,jgext,jetmax,
     $    jflgi,jgexti,jpres3,kxmod,mwtges,nbti,nbtmax,ncmpri,net,neti,
     $    netmax,ngexti,ninpts,ngexrt,nobswt,nodbmx,nopgmx,noprmx,
     $    noptmx,noutpt,nprob,nsbswt,ntitl,ntitmx,nttyo,nxicmx,nxmdmx,
     $    nxmod,nxti,nxtimx,pei,press,qend,qgexsh,qrderr,rho,scamas,
     $    tdspkg,tdspl,tempc,tgexp,tolbt,toldl,tolspf,ucospi,uebal,
     $    ugexj,ugexji,ugexmo,ugexp,ugexpi,ugexr,ugexsi,umemi,uobsw,
     $    uredox,usbsw,usoli,uspeci,utitl,uhfgex,uvfgex,uxkgex,uxmod,
     $    xbari,xgexsi,xhfgex,xlkgex,xvfgex,xlkmod,zgexj)
          if (qrderr) go to 990
c
        else
c
          write (nttyo,1210) uoldv
 1210     format(/' * Error - (XCON3/xcon3) Coding to implement',
     $    /7x,'reading an input file in "D" format has not been',
     $    /7x,'implemented for version level "',a3,'."')
          go to 990
c
        endif
      else
c
        write (nttyo,1110) uoldf
 1110   format(/' * Error - (XCON3/xcon3) Have unknown format',
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
      do n = 1,ntitl
        j = index(utitl(n),'Version level=')
        if (j .eq. 0) j = index(utitl(n),'version level=')
        if (j .gt. 0) then
          utitl(n)(j + 14:j + 14) = ' '
          utitl(n)(j + 15:j + 17) = unewv
          utitl(n)(j + 18:80) = ' '
          go to 220
        endif
        j = index(utitl(n),'Version number=')
        if (j .eq. 0) j = index(utitl(n),'version number=')
        if (j .gt. 0) then
          utitl(n)(j + 8:j + 14) = 'level= '
          utitl(n)(j + 15:j + 17) = unewv
          utitl(n)(j + 18:80) = ' '
          go to 220
        endif
      enddo
c
      if ((ntitl + 1) .gt. ntitmx) then
        write (nttyo,1115) nprob,ntitmx
 1115   format(/" * Error - (XCON3/xcon3) Can't add a version",
     $  /7x,'level marker to the first input file title of,'
     $  /7x,'problem ',i2,' because this title already has the',
     $  /7x,'maximum length of ',i3,' lines.')
        go to 990
      endif
c
      do n = ntitl,1,-1
        utitl(n + 1) = utitl(n)
      enddo
      ntitl = ntitl + 1
c
      utitl(1)(1:40) = '                                        '
      utitl(1)(41:80) = '                                        '
      utitl(1)(1:15)  = 'Version level= '
      utitl(1)(16:18) = unewv
      utitl(1)(21:27) = '(XCON3)'
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
c     Patch up possible incompatibilities between "D" and "W" formats.
c     The relevant variables here do not appear on EQ6 input files.
c
      if (vernew .lt. 8.0) then
        if (unewf(1:1) .eq. 'D') then
          if (uebal(1:5) .eq. 'none ') uebal = ' '
          if (uebal(1:5) .eq. 'None ') uebal = ' '
          if (uredox(1:5) .eq. 'none ') uredox = ' '
          if (uredox(1:5) .eq. 'None ') uredox = ' '
        elseif (unewf(1:1) .eq. 'W') then
          if (uebal(1:5) .eq. 'None ') uebal = 'none'
          if (uebal(1:1) .eq. ' ') uebal = 'none'
          if (uredox(1:5) .eq. 'None ') uredox = 'none'
          if (uredox(1:1) .eq. ' ') uredox = 'none'
        endif
      else
        if (uebal(1:5) .eq. 'none ') uebal = 'None'
        if (uebal(1:1) .eq. ' ') uebal = 'None'
        if (uredox(1:5) .eq. 'none ') uredox = 'None'
        if (uredox(1:1) .eq. ' ') uredox = 'None'
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
        jpres3 = 0
c
        if (tdspkg .gt. 0.) then
          itdsf3 = 0
        elseif (tdspl .gt. 0.) then
          itdsf3 = 1
        else
          itdsf3 = 0
        endif
c
        j2 = ilnobl(uebal)
        if (j2 .eq. 0) then
          iebal3 = 0
        elseif (uebal(1:j2) .eq. 'None') then
          iebal3 = 0
        elseif (uebal(1:j2) .eq. 'none') then
          iebal3 = 0
        elseif (uebal(1:j2) .ne. 'pick1.') then
          iebal3 = 1
        else
          write (nttyo,1608)
 1608     format(/" * Error - (XCON3/xcon3) Can't translate this",
     $    /7x,'pre-version level 8.0 input file to version level',
     $    /7x,'8.0 or higher because it specifies the option that',
     $    /7x,'the code is to pick a species for electrical balancing',
     $    /7x,'(uebal= "pick1."). This option does not exist in',
     $    /7x,'version level 8.0 or higher. Change the uebal input',
     $    /7x,'on the old input file to blank or the name of a species',
     $    /7x,'to balance on.')
          nerr = nerr + 1
        endif
c
        scamas = 1.
        net = 0
c
        iopt(11) = iopt(2)
        iopt(2) = 0
        iopt(17) = iopt(3)
        iopt(3) = 0
        irdxc3 = iopt(1)
        iopt(1) = 0
        if (iopt(4) .ge. 2) iopt(4) = 1
c
        iopr(10) = iopr(9)
        iopr(9) = iopr(6)
        iopr(6) = iopr(5)
        iopr(5) = 0
        iopr(3) = iopr(8)
        iopr(8) = 0
c
        iodb(3) = iodb(2)
        iodb(2) = 0
c
        nbti = nsq
        nsbswt = 0
        nobswt = 0
c
        if (irdxc3 .eq. -2) then
          pei = fep
        elseif (irdxc3 .eq. -1) then
          ehi = fep
        elseif (irdxc3 .eq. 0) then
          fo2lgi = fep
        endif
c
        j2 = ilnobl(uredox)
        if (j2 .le. 0) uredox = 'None'
        j2 = ilnobl(uebal)
        if (j2 .le. 0) uebal = 'None'
c
        do n = 1,nxmod
          uxmod(n)(1:24) = uxmd24(n)
        enddo
c
        do ns = 1,nsq
c
c         Basis species and basis switching.
c
          nb = ns
          uspeci(nb) = uspecb(ns)
          jflgi(nb) = jflagb(ns)
          covali(nb) = cspb(ns)
          if (uphas2(ns)(1:3) .eq. '   ') then
            ucospi(nb)(1:24) = uphas1(ns)
            ucospi(nb)(25:48) = uphas2(ns)
          else
            ucospi(nb)(1:24) = uphas2(ns)
            ucospi(nb)(25:48) = uphas1(ns)
          endif
          if (ubasis(ns)(1:3) .ne. '   ') then
            nobswt = nobswt + 1
            uobsw(1,nobswt) = uspecb(ns)
            uobsw(2,nobswt) = ubasis(ns)
          endif
        enddo
c
        do nb = 1,nbti
          jfli = jflgi(nb)
          if (jfli.eq.19 .or. jfli.eq.20 .or. jfli.eq.21)
     $    jflgi(nb) = 25
        enddo
c
        do nb = 1,nbti
          jfli = jflgi(nb)
          if (jfli .eq. 31) jflgi(nb) = 20
          if (jfli .eq. 32) jflgi(nb) = 21
        enddo
c
        do nb = 1,nbti
          jfli = jflgi(nb)
          if (jfli.ge.4 .and. jfli.le.8) then
            j2 = ilnobl(uspeci(nb)(1:24))
            write (nttyo,1610) jfli,uspeci(nb)(1:j2)
 1610       format(/" * Error - (XCON3/xcon3) Can't translate this",
     $      /7x,'pre-version level 8.0 input file to version level',
     $      /7x,'8.0 or higher because it employs a jflag value of',
     $      /7x,i3,' for ',a,'. This jflag option for specifying',
     $      /7x,'a "free" concentration" is obsolete at the higher',
     $      /7x,'version level.')
            nerr = nerr + 1
          endif
        enddo
c
        qgexsh = .false.
c
        nxti = nxtb
        nxic = 0
        do nx = 1,nxtb
          nxi = nx
          usoli(nxi) = usolb(nx)
          iktb = ncompb(nx)
          ncmpri(1,nxi) = nxic + 1
          ncmpri(2,nxi) = nxic + iktb
          do ik = 1,iktb
            nxic = nxic + 1
            umemi(nxic) = umemb(ik,nx)
            xbari(nxic) = xbarb(ik,nx)
          enddo
        enddo
c
        tolspf = tolsat
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
        qgexsh = .false.
c
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
        if (tempc .le. 100.) then
          presh = apresh(1,1)
        else
          presh = 0.
          do nn = 1,5
            n = 6 - nn
            presh = apresh(n,2) + tempc*presh
          enddo
        endif
c
        if (jpres3 .le. 0) then
          press = 0.
        elseif (jpres3 .eq. 1) then
          write (nttyo,1692)
 1692     format(/' * Note - (XCON3/xcon3) This version level 8.0',
     $    /7x,'or higher input file specifies the option to compute',
     $    /7x,'the pressure from the temperature according to the',
     $    /7x,'1.013-bar/steam-saturation curve (jpres3 = 1). This',
     $    /7x,'option does not directly translate to pre-Version 8.0.',
     $    /7x,'level The pressure at that level implicitly lies on the',
     $    /7x,'data file reference pressure curve, which may or may',
     $    /7x,'not coincide with the 1.013-bar/steam-saturation curve.',
     $    /7x,'To avoid getting this message, change the option on the',
     $    /7x,'old input file to specify computing the pressure',
     $    /7x,'according to the data file reference pressure curve',
     $    /7x,'(jpres3 = 0).')
          press = 0.
        elseif (jpres3 .eq. 2) then
          if (press.gt.0 .and. abs(press - presh).gt.1.e-4) then
            write (nttyo,1700) press,presh
 1700       format(/" * Error - (XCON3/xcon3) Can't translate this",
     $      /7x,'version level 8.0 or higher input file to a lower',
     $      /7x,'version level because it specifies a pressure of',
     $      /7x,f9.4,' bars, which differs from the pre-Version 8',
     $      /7x,'standard grid pressure of ',f9.4,' bars. To convert',
     $      /7x,'this input file, set press = 0., which defaults to',
     $      /7x,'the standard grid pressure.')
            nerr = nerr + 1
          endif
        else
          write (nttyo,1702) jpres3
 1702     format(/" * Error - (XCON3/xcon3) Can't translate this",
     $    /7x,'version level 8.0 or higher input file to a lower',
     $    /7x,'version level because it specifies an unknown',
     $    /7x,'pressure option (jpres3= ',i2,'.')
          nerr = nerr + 1
        endif
c
        if (itdsf3 .le. 0) then
          continue
        elseif (itdsf3 .ge. 1) then
          continue
        endif
c
        if (iebal3 .le. 0) then
          uebal = ' '
        else
          continue
        endif
c
        if (abs(scamas - 1.0) .gt. 1.e-6) then
          write (nttyo,1710) scamas
 1710     format(/" * Error - (XCON3/xcon3) Can't translate this",
     $    /7x,'version level 8.0 or higher input file to a lower',
     $    /7x,'version level because it specifies a  scale factor',
     $    /7x,'of ',1pe11.4,' for the mass of aqueous solution to',
     $    /7x,'write on the PICKUP file. Pre-version 8 input files',
     $    /7x,'lack this parameter, effectively assuming a value of',
     $    ' 1.0.',/7x,'To convert this input file, set scamas = 1.0.')
          nerr = nerr + 1
        endif
c
        if (net .gt. 0) then
          write (nttyo,1720) net
 1720     format(/" * Error - (XCON3/xcon3) Can't translate this",
     $    /7x,'version level 8.0 or higher input file to a lower',
     $    /7x,'version level because it defines ',i2,' generic',
     $    /7x,'ion exchangers.Pre-version 8 input files lack the',
     $    /7x,'generic ion exchange capability.')
          nerr = nerr + 1
        endif
c
        iopt(2) = iopt(11)
        iopt(11) = 0
        iopt(3) = iopt(17)
        iopt(17) = 0
        iopt(1) = irdxc3
        if (iopt(4).ge.1 .and. nxti.gt.0) iopt(4) = 2
c
        if (iopt(19) .gt. 0) then
          write (nttyo,1730)
 1730     format(/" * Warning - (XCON3/xcon3) Can't translate the",
     $    /7x,'"Advanced EQ3NR PICKUP File Options" (iopt(19)) option',
     $    /7x,'specified on this version level 8.0 or higher input',
     $    /7x,'file to a lower version level because this option does',
     $    /7x,'not exist at the lower version level.')
        endif
c
        iopr(5) = iopr(6)
        iopr(6) = iopr(9)
        iopr(9) = iopr(10)
        iopr(10) = 0
        iopr(8) = iopr(3)
        iopr(3) = 0
c
        iodb(2) = iodb(3)
        iodb(3) = 0
c
        if (iopt(1) .eq. -2) then
          fep = pei
        elseif (iopt(1) .eq. -1) then
          fep = ehi
        elseif (iopt(1) .eq. 0) then
          fep = fo2lgi
        endif
c
cXX
        do nb = 1,nbti
          if (jflgi(nb) .eq. 22) then
c
c           Translate the pmH option to a pH option.
c           This is done at the version 8 level. Additional
c           translation may follow.
c
            jflgi(nb) = 20
            iopg(2) = 1
          endif
        enddo
c
        do nb = 1,nbti
          if (jflgi(nb) .eq. 23) then
c
c           The pmX option is not translatable.
c
            j2 = ilnobl(uspeci(nb))
            write (nttyo,1732) uspeci(nb)(1:j2)
 1732       format(/" * Error - (XCON3/xcon3) Can't translate this",
     $      /7x,'version level 8.0 or higher input file to a lower',
     $      /7x,'version level because it contains a pmX input for',
     $      /7x,a,'. The pmX option does not exist at a lower',
     $      /7x,'version level.')
            nerr = nerr + 1
          endif
        enddo
cXX
c
        nsq = nbti
c
        do nb = 1,nbti
c
c         Basis species.
c
          ns = nb
          uspecb(ns) = uspeci(nb)(1:24)
          jflagb(ns) = jflgi(nb)
          cspb(ns) = covali(nb)
          if (ucospi(nb)(25:27) .eq. '   ') then
            uphas1(ns) = ucospi(nb)(1:24)
            uphas2(ns) = ' '
          elseif (ucospi(nb)(25:33) .eq. 'Aqueous solution ') then
            uphas1(ns) = ucospi(nb)(1:24)
            uphas2(ns) = ' '
          elseif (ucospi(nb)(25:28) .eq. 'Gas ') then
            uphas1(ns) = ucospi(nb)(1:24)
            uphas2(ns) = ' '
          else
            uphas2(ns) = ucospi(nb)(1:24)
            uphas1(ns) = ucospi(nb)(25:48)
          endif
        enddo
c
        do ns = 1,nsq
          if (jflagb(ns) .eq. 19) then
            jflagb(ns) = 16
            cspb(ns) = -cspb(ns)
          endif
        enddo
c
        if (unewf(1:1) .eq. 'W') then
          do ns = 1,nsq
            if (jflagb(ns) .eq. 20) then
              jflagb(ns) = 16
              cspb(ns) = -cspb(ns)
            elseif (jflagb(ns) .eq. 21) then
              uphas1(ns) = 'Cl-'
              jflagb(ns) = 17
              cspb(ns) = -cspb(ns)
            endif
          enddo
        else
          do ns = 1,nsq
            if (jflagb(ns) .eq. 20) then
              jflagb(ns) = 31
            elseif (jflagb(ns) .eq. 21) then
              jflagb(ns) = 32
            endif
          enddo
        endif
c
        do ns = 1,nsq
          if (jflagb(ns) .eq. 25) then
            jflagb(ns) = 19
            if (uphas2(ns)(1:24).ne.ublk24(1:24) .and.
     $      uphas2(ns)(1:24).ne.uphas1(ns)(1:24)) jflagb(ns) = 20
            j = index(uphas1(ns),'(g)')
            if (j .gt. 0) jflagb(ns) = 21
          endif
        enddo
c
        if (nsbswt .gt. 0) then
          write (nttyo,1740) nsbswt
 1740     format(/" * Error - (XCON3/xcon3) Can't translate this",
     $    /7x,'version level 8.0 or higher input file to a lower',
     $    /7x,'version level because it contains ',i3,' directives',
     $    /7x'for special basis switching.')
          nerr = nerr + 1
        endif
c
        do nobsw = 1,nobswt
c
c         Basis switching.
c
          do nb = 1,nbti
            if (uobsw(1,nobsw)(1:48) .eq. uspeci(nb)(1:48)) then
              ns = nb
              ubasis(ns) = uobsw(2,nobsw)
              go to 244
            endif
          enddo
  244     continue
        enddo
c
        do n = 1,nxmod
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
          if (uxmod(n)(25:40) .eq. 'Aqueous solution') then
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
            write (nttyo,1117) uxmod(n)(1:j2)
 1117       format(/" * Error - (XCON3/xcon3) Can't translate a",
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
 1119     format(/" * Warning - (XCON3/xcon3) Can't unambiguously",
     $    /7x,'determine what kind of species in an alter/suppress',
     $    /7x,'option is "',a,'".',
     $    /7x,'Setting jxmod to 1. Check to see that this is correct',
     $    /7x,'(0 = aqueous species, 1 = pure mineral,',
     $    ' 3 = gas species',
     $    /7x,'3 = solid solution).')
          jxmod(n) = 1
  250    continue
       enddo
c
        nxtb = nxti
        do nxi = 1,nxti
          nx = nxi
          usolb(nx) = usoli(nxi)
          nr1 = ncmpri(1,nxi)
          nr2 = ncmpri(2,nxi)
          ncompb(nx) =  nr2 - nr1 + 1
          ik = 0
          do nxic = nr1,nr2
            ik = ik + 1
            umemb(ik,nx) = umemi(nxic)
            xbarb(ik,nx) = xbari(nxic)
          enddo
        enddo
c
        tolsat = tolspf
c
        if (nerr .gt. 0) stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (vernew.lt.8.0 .and. unewf(1:1).eq."D") then
c
c       Map any inputs equivalent to pH to pH, and any inputs equivalent
c       to pHCl to pHCl, but only if the new input file is to be in "D"
c       format. Use the hidden jflag values 31 (pH) and 32 (pHCl), which
c       do not exist in the case of "W" format.
c
        do ns = 1,nsq
          ustr = uspecb(ns)
          if (ustr(1:3).eq.'h+ ' .or. ustr(1:3).eq.'H+ ') then
            if (jflagb(ns) .eq. 16) then
              jflagb(ns) = 31
              cspb(ns) = -cspb(ns)
            elseif (jflagb(ns) .eq. 17) then
              ustr1 = uphas1(ns)
              if (ustr1(1:4).eq.'cl- ' .or. ustr1(1:4).eq.'Cl- ')
     $          then
                uphas1(ns) = ' '
                jflagb(ns) = 32
                cspb(ns) = -cspb(ns)
              endif
            endif
          endif
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (vernew .ge. 8.0) then
c
c       Map any inputs equivalent to pH to pH, and any inputs equivalent
c       to pHCl to pHCl.
c
        do nb = 1,nbti
          ustr = uspeci(nb)
          if (ustr(1:3) .eq. 'H+ ') then
            jfli = jflgi(nb)
            if (jfli .eq. 16) then
              jflgi(nb) = 20
              covali(nb) = -covali(nb)
            elseif (jfli .eq. 19) then
              jflgi(nb) = 20
            elseif (jfli .eq. 17) then
              ustr1 = ucospi(nb)
              if (ustr1(1:4).eq.'Cl- ') then
                ucospi(nb) = ' '
                jflgi(nb) = 21
                covali(nb) = -covali(nb)
              endif
            endif
          endif
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the current problem on the new input file.
c
      if (unewf(1:1) .eq. 'W') then
        if (unewv(1:3) .eq. '6.0') then
c
          call wr3w6(cspb,fep,iktmax,iodb,iopg,iopr,iopt,
     $    itermx,jflagb,jxmod,kxmod,ncompb,newin,nodbmx,nopgmx,
     $    noprmx,noptmx,nsq,nsqmax,ntitl,ntitmx,nxmdmx,nxmod,
     $    nxtb,nxtmax,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,
     $    uacion,ubasis,uebal,umemb,uphas1,uphas2,uredox,usolb,
     $    uspecb,utitl,xbarb,uxmd24,xlkmod)
c
        elseif (unewv(1:3) .eq. '7.0') then
c
          call wr3w7(cspb,fep,iktmax,iodb,iopg,iopr,iopt,
     $    itermx,jflagb,jxmod,kxmod,ncompb,newin,nodbmx,nopgmx,
     $    noprmx,noptmx,nsq,nsqmax,ntitl,ntitmx,nxmdmx,nxmod,
     $    nxtb,nxtmax,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,
     $    ubasis,uebal,umemb,uphas1,uphas2,uredox,usolb,uspecb,
     $    utitl,xbarb,uxmd24,xlkmod)
c
        elseif (unewv(1:3) .eq. '7.2') then
c
c         Note: version level 7.2 is identical to version level 7.0
c         for this format.
c
          call wr3w7(cspb,fep,iktmax,iodb,iopg,iopr,iopt,
     $    itermx,jflagb,jxmod,kxmod,ncompb,newin,nodbmx,nopgmx,
     $    noprmx,noptmx,nsq,nsqmax,ntitl,ntitmx,nxmdmx,nxmod,
     $    nxtb,nxtmax,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,
     $    ubasis,uebal,umemb,uphas1,uphas2,uredox,usolb,uspecb,
     $    utitl,xbarb,uxmd24,xlkmod)
c
        elseif (unewv(1:3) .eq. '8.0') then
c
          call wr3w8(cgexj,cgexpi,covali,ehi,egexsi,fo2lgi,iebal3,
     $    ietmax,iodb,iopg,iopr,iopt,irdxc3,itdsf3,itermx,jgext,jetmax,
     $    jflgi,jgexti,jpres3,kxmod,mwtges,nbti,nbtmax,ncmpri,net,neti,
     $    netmax,newin,ngexti,ngexrt,nobswt,nodbmx,nopgmx,noprmx,noptmx,
     $    nsbswt,ntitl,ntitmx,nttyo,nxicmx,nxmdmx,nxmod,nxti,nxtimx,pei,
     $    press,qgexsh,rho,scamas,tdspkg,tdspl,tempc,tgexp,tolbt,toldl,
     $    tolspf,ucospi,uebal,ugexj,ugexji,ugexmo,ugexp,ugexpi,ugexr,
     $    ugexsi,umemi,uobsw,uredox,usbsw,usoli,uspeci,utitl,uhfgex,
     $    uvfgex,uxkgex,uxmod,xbari,xgexsi,xhfgex,xlkgex,xvfgex,
     $    xlkmod,zgexj)
c
        else
c
          write (nttyo,1220) unewv
 1220     format(/' * Error - (XCON3/xcon3) Coding to implement',
     $    /7x,'writing an input file in "W" format has not been',
     $    /7x,'implemented for version level "',a3,'."')
          go to 990
c
        endif
      elseif (unewf(1:1) .eq. 'D') then
        if (unewv(1:3) .eq. '6.0') then
c
          write (nttyo,1230) unewv
 1230     format(/' * Error - (XCON3/xcon3) There is no "D" format',
     $    /7x,'for version level "6.0", hence',"can't write an input",
     $    /7x,'file in this format for this version level.')
          go to 990
c
        elseif (unewv(1:3) .eq. '7.0') then
c
          call wr3d7(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,
     $    jflagb,jxmod,kxmod,ncompb,newin,nodbmx,nopgmx,noprmx,noptmx,
     $    nsq,nsqmax,ntitl,ntitmx,nttyo,nxmdmx,nxmod,nxtb,nxtmax,
     $    rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,ubasis,uebal,
     $    umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,xbarb,
     $    uxmd24,xlkmod)
c
        elseif (unewv(1:3) .eq. '7.2') then
c
          call wr3d72(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,
     $    jflagb,jxmod,kxmod,ncompb,newin,nodbmx,nopgmx,noprmx,noptmx,
     $    nsq,nsqmax,ntitl,ntitmx,nttyo,nxmdmx,nxmod,nxtb,nxtmax,
     $    rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,ubasis,uebal,
     $    umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,xbarb,
     $    uxmd24,xlkmod)
c
        elseif (unewv(1:3) .eq. '8.0') then
c
          call wr3d8(cgexj,cgexpi,covali,ehi,egexsi,fo2lgi,iebal3,
     $    ietmax,iodb,iopg,iopr,iopt,irdxc3,itdsf3,itermx,jgext,jetmax,
     $    jflgi,jgexti,jpres3,kxmod,mwtges,nbti,nbtmax,ncmpri,net,neti,
     $    netmax,newin,ngexti,ngexrt,nobswt,nodbmx,nopgmx,noprmx,noptmx,
     $    nsbswt,ntitl,ntitmx,nttyo,nxicmx,nxmdmx,nxmod,nxti,nxtimx,pei,
     $    press,qgexsh,rho,scamas,tdspkg,tdspl,tempc,tgexp,tolbt,toldl,
     $    tolspf,ucospi,uebal,ugexj,ugexji,ugexmo,ugexp,ugexpi,ugexr,
     $    ugexsi,umemi,uobsw,uredox,usbsw,usoli,uspeci,utitl,uhfgex,
     $    uvfgex,uxkgex,uxmod,xbari,xgexsi,xhfgex,xlkgex,xvfgex,
     $    xlkmod,zgexj)
c
        else
c
          write (nttyo,1240) unewv
 1240     format(/' * Error - (XCON3/xcon3) Coding to implement',
     $    /7x,'writing an input file in "D" format has not been',
     $    /7x,'implemented for version level "',a3,'."')
          go to 990
c
        endif
      else
c
        write (nttyo,1130) unewf
 1130   format(/' * Error - (XCON3/xcon3) Have unknown format',
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
 1140 format(/' * Error - (XCON3/xcon3) Have encountered a read',
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
