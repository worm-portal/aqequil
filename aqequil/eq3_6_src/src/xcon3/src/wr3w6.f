      subroutine wr3w6(cspb,fep,iktmax,iodb,iopg,iopr,iopt,
     $ itermx,jflagb,jxmod,kxmod,ncompb,newin,nodbmx,nopgmx,
     $ noprmx,noptmx,nsq,nsqmax,ntitl,ntitmx,nxmdmx,nxmod,
     $ nxtb,nxtmax,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,
     $ uacion,ubasis,uebal,umemb,uphas1,uphas2,uredox,usolb,
     $ uspecb,utitl,xbarb,uxmd24,xlkmod)
c
c     This subroutine writes the EQ3NR input file in compact ("W")
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
      integer iktmax,nodbmx,nopgmx,noprmx,noptmx,nsqmax,ntitmx,nxmdmx,
     $ nxtmax
c
      integer iodb(nodbmx),iopg(nopgmx),iopr(noprmx),iopt(noptmx),
     $ jflagb(nsqmax),jxmod(nxmdmx),kxmod(nxmdmx),ncompb(nxtmax)
c
      integer itermx,newin,nsq,ntitl,nxmod,nxtb
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
      integer i,ikb,iktb,j2,n,ns,nxb
      integer ilnobl
c
c-----------------------------------------------------------------------
c
c     Write the new input file.
c
c     Title.
c
      do 110 n = 1,ntitl
        j2 = ilnobl(utitl(n))
        write (newin,1000) utitl(n)(1:j2)
 1000 format(a)
  110 continue
      if (ntitl .lt. ntitmx) write (newin,1010)
 1010 format('endit.')
c
c-----------------------------------------------------------------------
c
c     Temperature.
c
      write (newin,1040) tempc
 1040 format(5x,'tempc= ',1pe12.5)
c
c     Density, total dissolved salts (per kg), and total dissolved
c     salts (per liter).
c
      write (newin,1050) rho,tdspkg,tdspl
 1050 format (7x,'rho= ',1pe12.5,4x,'tdspkg= ',1pe12.5,5x,
     $ 'tdspl= ',1pe12.5)
c
c     Redox parameter and name of redox species defining a controlling
c     couple, if any.
c
      j2 = ilnobl(uredox)
      write (newin,1060) fep,uredox(1:j2)
 1060 format (7x,'fep= ',1pe12.5,4x,'uredox= ',a)
c
c     Convergence tolerances (tolbt and toldl) and saturation
c     tolerance (tolsat).
c
      write (newin,1070) tolbt,toldl,tolsat
 1070 format (5x,'tolbt= ',1pe12.5,5x,'toldl= ',1pe12.5,4x,'tolsat= ',
     $ 1pe12.5)
c
c     Maximum number of iterations.
c
      write (newin,1080) itermx
 1080 format (4x,'itermx= ',i2)
c
c-----------------------------------------------------------------------
c
      write (newin,1090)
 1090 format('*',15x,'1    2    3    4    5    6    7    8    9   10')
c
c     Iopt option switches.
c     Note: iopt(1) = iopt1, etc.
c
      write (newin,1100) (iopt(i), i = 1,10)
 1100 format(2x,'iopt1-10= ',10i5)
c
c     Iopg option switches.
c     Note: iopg(1) = iopg1, etc.
c
      write (newin,1110) (iopg(i), i = 1,10)
 1110 format(2x,'iopg1-10= ',10i5)
c
c     Iopr option switches.
c     Note: iopr(1) = iopr1, etc.
c
      write (newin,1120) (iopr(i), i = 1,20)
 1120 format(2x,'iopr1-10= ',10i5,/1x,'iopr11-20= ',10i5)
c
c     Iodb option switches.
c     Note: iodb(1) = iodb1, etc.
c
      write (newin,1130) (iodb(i), i = 1,10)
 1130 format(2x,'iodb1-10= ',10i5)
c
c-----------------------------------------------------------------------
c
c     Species for electrical balancing.
c
      j2 = ilnobl(uebal)
      write (newin,1150) uebal(1:j2)
 1150 format(5x,'uebal= ',a)
c
c-----------------------------------------------------------------------
c
c     Species for defining equivalent stoichiometric ionic strength.
c
      j2 = ilnobl(uacion)
      write (newin,1160) uacion(1:j2)
 1160 format(4x,'uacion= ',a)
c
c-----------------------------------------------------------------------
c
c     Nxmod options.
c
      write (newin,1200) nxmod
 1200 format(5x,'nxmod= ',i2)
c
      if (nxmod .gt. 0) then
        do 420 n = 1,nxmod
          j2 = ilnobl(uxmd24(n))
          write (newin,1220) uxmd24(n)(1:j2),jxmod(n),kxmod(n),
     $    xlkmod(n)
 1220     format(3x,'species= ',a,/6x,'type= ',i2,14x,'option= ',
     $    i2,14x,'xlkmod= ',1pe12.5)
  420   continue
      endif
c
c-----------------------------------------------------------------------
c
c     Basis species and associated constraints.
c
      do 440 ns = 1,nsq
        j2 = ilnobl(uspecb(ns))
        write (newin,1240) uspecb(ns)(1:j2)
 1240   format('data file master species= ',a)
        j2 = ilnobl(ubasis(ns))
        write (newin,1250) ubasis(ns)(1:j2)
 1250   format(3x,'switch with species= ',a)
        write (newin,1260) jflagb(ns),cspb(ns)
 1260   format(3x,'jflag= ',i2,3x,'csp= ',1pg12.5)
        if (jflagb(ns).ge.17 .and. jflagb(ns).le.21) then
          j2 = ilnobl(uphas2(ns))
          write (newin,1270) uphas1(ns),uphas2(ns)(1:j2)
 1270     format(2x,'uphas1= ',a24,3x,'uphas2= ',a)
        endif
  440 continue
      write (newin,1280)
 1280 format('endit.')
c
c-----------------------------------------------------------------------
c
c     Mole fractions of solid solutions.
c
      if (iopt(4) .ge. 2) then
        write (newin,1300)
 1300   format('*   Solid solution compositions')
        do 460 nxb = 1,nxtb
          j2 = ilnobl(usolb(nxb))
          write (newin,1310) usolb(nxb)(1:j2)
 1310     format(3x,a)
          iktb = ncompb(nxb)
          do 450 ikb = 1,iktb
            write (newin,1320) umemb(ikb,nxb),xbarb(ikb,nxb)
 1320       format(6x,a24,3x,f10.4)
  450     continue
          write (newin,1330)
 1330     format(6x,'endit.')
  460   continue
        write (newin,1340)
 1340   format(3x,'endit.')
      endif
c
c-----------------------------------------------------------------------
c
  999 continue
      end
