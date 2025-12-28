      subroutine rd3w6(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,
     $ jflagb,jxmod,kxmod,ncompb,ninpts,nodbmx,nopgmx,noprmx,noptmx,
     $ nsq,nsqmax,ntitl,ntitmx,nttyo,nxmdmx,nxmod,nxtb,nxtmax,qend,
     $ qrderr,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,uacion,
     $ ubasis,uebal,umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,
     $ xbarb,uxmd24,xlkmod)
c
c     This subroutine reads the EQ3NR input file in compact ("W")
c     format for versions 6.0-6.1.
c
c     This subroutine is called by:
c
c       XCON3/xcon3.f
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer iktmax,nodbmx,nopgmx,noprmx,noptmx,nsqmax,ntitmx,
     $ nxmdmx,nxtmax
c
      integer iodb(nodbmx),iopg(nopgmx),iopr(noprmx),iopt(noptmx),
     $ jflagb(nsqmax),jxmod(nxmdmx),kxmod(nxmdmx),ncompb(nxtmax)
c
      integer itermx,ninpts,nsq,ntitl,nttyo,nxmod,nxtb
c
      logical qend,qrderr
c
      character*80 utitl(ntitmx)
      character*24 ubasis(nsqmax),umemb(iktmax,nxtmax),uphas1(nsqmax),
     $ uphas2(nsqmax),usolb(nxtmax),uspecb(nsqmax),uxmd24(nxmdmx)
      character*24 uacion,uebal,uredox
c
      real*8 cspb(nsqmax),xbarb(iktmax,nxtmax),xlkmod(nxmdmx)
c
      real*8 fep,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,iktb,j,j2,n
      integer ilnobl
c
      character*80 uline
      character*24 ux24
      character*8 uendit,ux8
c
      real*8 xx
c
c-----------------------------------------------------------------------
c
      data uendit /'endit.  '/
c
c-----------------------------------------------------------------------
c
      qrderr = .false.
c
c-----------------------------------------------------------------------
c
c     Main title.
c
c     Note: if there are exactly ntitpa lines in the title, the title
c     need not be terminated by an 'endit.'. The 'endit.' if present
c     is here not considered to be part of the title. The secondary
c     title is treated in the same manner.
c
      qend = .false.
      read (ninpts,1000,end=100,err=990) uline
 1000 format(a80)
      go to 105
c
  100 qend = .true.
      go to 999
c
  105 n = 1
      if (uline(1:8) .eq. uendit(1:8)) go to 120
      utitl(1) = uline
c
      do 110 n = 2,ntitmx
        read (ninpts,1000,err=990) uline
        if (uline(1:8) .eq. uendit(1:8)) go to 120
        utitl(n) = uline
  110 continue
      n = n + 1
c
      read (ninpts,1000,err=990) uline
      call locase(uline)
      j = index(uline,'tempc=')
      if (j .gt. 0) then
        backspace(ninpts)
      else
        write (nttyo,1030) ntitmx
 1030   format(/' * Error - (XCON3/rd3w6) Have too many lines in the',
     $  /7x,'secondary title. The code is only dimensioned for ',i4,
     $  /7x,'lines. Reduce the size of the title or increase the',
     $  /7x,'dimensioning parameter ntitpa.')
        go to 990
      endif
c
  120 ntitl = n - 1
c
c-----------------------------------------------------------------------
c
c     Temperature.
c
      read (ninpts,1040,err=990) tempc
 1040 format(12x,e12.5)
c
c     Density, total dissolved salts (per kg), and total dissolved
c     salts (per liter).
c
      read (ninpts,1050,err=990) rho,tdspkg,tdspl
 1050 format (3(12x,e12.5))
c
c     Redox parameter and name of redox species defining a controlling
c     couple, if any.
c
      read (ninpts,1060,err=990) fep,uredox
 1060 format (12x,e12.5,12x,a24)
c
c     Convergence tolerances (tolbt and toldl) and saturation
c     tolerance (tolsat).
c
      read (ninpts,1070,err=990) tolbt,toldl,tolsat
 1070 format (3(12x,e12.5))
c
c     Maximum number of iterations.
c
      read (ninpts,1080,err=990) itermx
 1080 format (12x,i2)
c
c-----------------------------------------------------------------------
c
c     Iopt option switches.
c     Note: iopt(1) = iopt1, etc.
c
      read (ninpts,1100,err=990) (iopt(i), i = 1,10)
 1100 format(12x,10i5)
c
c     Iopg option switches.
c     Note: iopg(1) = iopg1, etc.
c
      read (ninpts,1110,err=990) (iopg(i), i = 1,10)
 1110 format(12x,10i5)
c
c     Iopr option switches.
c     Note: iopr(1) = iopr1, etc.
c
      read (ninpts,1120,err=990) (iopr(i), i = 1,20)
 1120 format(12x,10i5)
c
c     Iodb option switches.
c     Note: iodb(1) = iodb1, etc.
c
      read (ninpts,1130,err=990) (iodb(i), i = 1,10)
 1130 format(12x,10i5)
c
c-----------------------------------------------------------------------
c
c     Species for electrical balancing.
c
      read (ninpts,1150,err=990) uebal
 1150 format(12x,a24)
c
c-----------------------------------------------------------------------
c
c     Species for defining equivalent stoichiometric ionic strength.
c
      read (ninpts,1160,err=990) uacion
 1160 format(12x,a24)
c
c-----------------------------------------------------------------------
c
c     Number of nxmod options.
c
      read (ninpts,1200,err=990) nxmod
 1200 format(12x,i2)
c
      if (nxmod .gt. nxmdmx) then
        write (nttyo,1210) nxmdmx
 1210   format(/' * Error - (XCON3/rd3w6) Have too many nxmod',
     $  /7x,'alter/suppress options. The code is only dimensioned',
     $  /7x,'for ',i3,' such options. Reduce the number of such',
     $  ' options',/7x,'or increase the dimensioning parameter',
     $  ' nxmdpa.')
        go to 990
      endif
c
c     Nxmod options.
c
      if (nxmod .gt. 0) then
        do 420 n = 1,nxmod
          read (ninpts,1220,err=990) uxmd24(n),jxmod(n),kxmod(n),
     $    xlkmod(n)
 1220     format(12x,a24,/12x,i2,22x,i2,22x,e12.5)
  420   continue
      endif
c
c-----------------------------------------------------------------------
c
c     Basis species and associated constraints.
c
      nsq = 0
c
  430 read (ninpts,1250,err=990) ux8,ux24
 1250 format(a8,18x,a24)
      if (ux8(1:8) .eq. uendit(1:8)) go to 440
c
      nsq = nsq + 1
c
      if (nsq .gt. nsqmax) then
        write (nttyo,1260) nsqmax
 1260   format(/' * Error - (XCON3/rd3w6) Have too many basis',
     $  /7x,'species present. The code is only dimensioned',
     $  /7x,'for ',i3,' basis species. Reduce the number of such',
     $  /7x,'species or increase the dimensioning parameter nsqpar.')
        go to 990
      endif
c
      uspecb(nsq) = ux24
c
      read (ninpts,1270,err=990) ubasis(nsq)
 1270 format(24x,a24)
c
      read (ninpts,1280,err=990) jflagb(nsq),cspb(nsq)
 1280 format(10x,i2,8x,e12.5)
c
      if (jflagb(nsq).ge.17 .and. jflagb(nsq).le.21) then
        read (ninpts,1290,err=990) uphas1(nsq),uphas2(nsq)
 1290   format(10x,a24,11x,a24)
      endif
      go to 430
c
  440 continue
c
c-----------------------------------------------------------------------
c
c     Mole fractions of solid solutions.
c
      nxtb = 0
      if (iopt(4) .ge. 2) then
c
  450   read (ninpts,1300,err=990) ux24
 1300   format(3x,a24)
        if (ux24(1:8) .eq. uendit(1:8)) go to 470
c
        nxtb = nxtb + 1
c
        if (nxtb .gt. nxtmax) then
          write (nttyo,1310) nxtmax
 1310     format(/' * Error - (XCON3/rd3w6) Have too many solid',
     $    /7x,'solutions present. The code is only dimensioned',
     $    /7x,'for ',i3,' solid solutions. Reduce the number of such',
     $    /7x,'phases or increase the dimensioning parameter nxtpar.')
          go to 990
        endif
c
        usolb(nxtb) = ux24
        iktb = 0
c
  460   read (ninpts,1320,err=990) ux24,xx
 1320   format(6x,a24,3x,f10.4)
        if (ux24(1:8) .eq. uendit(1:8)) then
          ncompb(nxtb) = iktb
          go to 450
        endif
c
        iktb = iktb + 1
c
        if (iktb .gt. iktmax) then
          j2 = ilnobl(usolb(nxtb))
          write (nttyo,1330) usolb(nxtb)(1:j2),iktmax
 1330     format(/' * Error - (XCON3/rd3w6) Solid solution',
     $    /7x,'"',a,'" has too many end-members present.',
     $    /7x,'This code is only dimensioned for ',i3,' end-members',
     $    /7x,'per solid solution. Reduce the number of end-members',
     $    /7x,'or increase the dimensioning parameter iktpar.')
          go to 990
        endif
c
        umemb(iktb,nxtb) = ux24
        xbarb(iktb,nxtb) = xx
        go to 460
c
      endif
c
  470 continue
c
c-----------------------------------------------------------------------
c
      go to 999
c
  990 qrderr = .true.
c
  999 continue
      end
