      subroutine goback(acflg,acflg0,emop,emop0,emos,emos0,fje,
     $ fje0,fxi,fxi0,iemop,iemop0,iemos,iemos0,iindx0,iindx1,ipndx0,
     $ ipndx1,jpflag,jsflag,jreac,jreac0,kdim,kdim0,kmax,km1,km10,
     $ kmt,kmt0,kx1,kx10,kxt,kxt0,loph,losp,moph,moph0,mosp,mosp0,
     $ ncmpe,ncmpe0,npet,npetmx,npet0,npt,nptmax,nrct,nrctmx,nset,
     $ nsetmx,nset0,nst,nstmax,qreq,qriinf,sigmam,sigmm0,uzvec0,
     $ uzvec1,xi0,xi1)
c
c     This subroutine sets up to go back to the previous point of
c     reaction progress.
c
c     This subroutine is called by:
c
c       EQ6/path.f
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
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer kmax,npetmx,nptmax,nrctmx,nsetmx,nstmax
c
      integer iemop(npetmx),iemop0(npetmx),iemos(nsetmx),iemos0(nsetmx),
     $ iindx0(kmax),iindx1(kmax),ipndx0(kmax),ipndx1(kmax),
     $ jreac(nrctmx),jreac0(nrctmx),jpflag(nptmax),jsflag(nstmax),
     $ ncmpe(2,npetmx),ncmpe0(2,npetmx)
c
      integer kdim,kdim0,km1,km10,kmt,kmt0,kx1,kx10,kxt,kxt0,npet,npet0,
     $ npt,nrct,nset,nset0,nst
c
      logical qreq,qriinf
c
      character(len=48) uzvec0(kmax),uzvec1(kmax)
c
      real(8) acflg(nstmax),acflg0(nstmax),emop(npetmx),emop0(npetmx),
     $ emos(nsetmx),emos0(nsetmx),loph(nptmax),losp(nstmax),
     $ moph(nptmax),moph0(nptmax),mosp(nstmax),mosp0(nstmax)
c
      real(8) fje,fje0,fxi,fxi0,sigmam,sigmm0,xi0,xi1
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ileft,kcol,np,npe,nrc,ns,nse
c
c-----------------------------------------------------------------------
c
      xi1 = xi0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      qreq = .false.
      qriinf = .false.
c
      kdim = kdim0
      kxt = kxt0
      kmt = kmt0
      km1 = km10
      kx1 = kx10
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do kcol = 1,kdim
        uzvec1(kcol) = uzvec0(kcol)
        iindx1(kcol) = iindx0(kcol)
        ipndx1(kcol) = ipndx0(kcol)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do np = 1,npt
        moph(np) = moph0(np)
      enddo
c
c     Note that the loop is unrolled.
c
      ileft = (npt/8)*8
      do np = 1,ileft,8
        loph(np) = -99999.
        loph(np + 1) = -99999.
        loph(np + 2) = -99999.
        loph(np + 3) = -99999.
        loph(np + 4) = -99999.
        loph(np + 5) = -99999.
        loph(np + 6) = -99999.
        loph(np + 7) = -99999.
        jpflag(np) = max(jpflag(np),0)
        jpflag(np + 1) = max(jpflag(np + 1),0)
        jpflag(np + 2) = max(jpflag(np + 2),0)
        jpflag(np + 3) = max(jpflag(np + 3),0)
        jpflag(np + 4) = max(jpflag(np + 4),0)
        jpflag(np + 5) = max(jpflag(np + 5),0)
        jpflag(np + 6) = max(jpflag(np + 6),0)
        jpflag(np + 7) = max(jpflag(np + 7),0)
      enddo
c
      do np = ileft + 1,npt
        loph(np) = -99999.
        jpflag(np) = max(jpflag(np),0)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      fje = fje0
      fxi = fxi0
      sigmam = sigmm0
c
      do ns = 1,nst
        mosp(ns) = mosp0(ns)
        acflg(ns) = acflg0(ns)
      enddo
c
c     Note that the loop is unrolled.
c
      ileft = (nst/8)*8
      do ns = 1,ileft,8
        losp(ns) = -99999.
        losp(ns + 1) = -99999.
        losp(ns + 2) = -99999.
        losp(ns + 3) = -99999.
        losp(ns + 4) = -99999.
        losp(ns + 5) = -99999.
        losp(ns + 6) = -99999.
        losp(ns + 7) = -99999.
        jsflag(ns) = max(jsflag(ns),0)
        jsflag(ns + 1) = max(jsflag(ns + 1),0)
        jsflag(ns + 2) = max(jsflag(ns + 2),0)
        jsflag(ns + 3) = max(jsflag(ns + 3),0)
        jsflag(ns + 4) = max(jsflag(ns + 4),0)
        jsflag(ns + 5) = max(jsflag(ns + 5),0)
        jsflag(ns + 6) = max(jsflag(ns + 6),0)
        jsflag(ns + 7) = max(jsflag(ns + 7),0)
      enddo
c
      do ns = ileft + 1,nst
        losp(ns) = -99999.
        jsflag(ns) = max(jsflag(ns),0)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (kxt .ge. km1) then
        do kcol = km1,kxt
          ns = iindx1(kcol)
          np = ipndx1(kcol)
          jpflag(np) = -1
          jsflag(ns) = -1
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do npe = 1,npet0
        iemop(npe) = iemop0(npe)
        emop(npe) = emop0(npe)
        ncmpe(1,npe) = ncmpe0(1,npe)
        ncmpe(2,npe) = ncmpe0(2,npe)
      enddo
      npet = npet0
c
      do nse = 1,nset0
        iemos(nse) = iemos0(nse)
        emos(nse) = emos0(nse)
      enddo
      nset = nset0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do nrc = 1,nrct
        jreac(nrc) = jreac0(nrc)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
