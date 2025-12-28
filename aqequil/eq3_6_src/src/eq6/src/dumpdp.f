      subroutine dumpdp(csts,demop0,demos0,emop,emop0,emos,
     $ emos0,fdpe0,fdpem1,fdse0,fdsem1,iemop,iemos,iern1,iern2,
     $ ietmax,iindx0,iindx1,imrn1,imrn2,iodb,ipndx0,ipndx1,ixrn1,
     $ ixrn2,jcsort,jern1,jetmax,jgext,jpflag,jsflag,kbt,kdim,
     $ kdim0,kmax,km1,km10,kmt,kmt0,kord,kstep,kx1,kx10,kxt,
     $ kxt0,loph,losp,moph,mosp,mprph,mprsp,mrgexs,mtb,mtb0,
     $ nbasp,nbaspd,nbt,nbtmax,ncmpe,ncmpr,ndelay,netmax,
     $ ngext,nodbmx,nordmx,noutpt,npet,npetmx,npet0,npt,nptmax,
     $ npts,nset,nsetmx,nset0,nstmax,nsts,nstsmx,nstsr,nttyo,
     $ qbye,uaqsln,ufixf,uspec,uphase,uzvec0,uzvec1,xbar,
     $ xbarlg,zklgmn,zklogl,zvclg0,zvclg1,zvec0,zvec1)
c
c     This subroutine transfers the entire mass of eligible phases
c     in the Equilibrium System (ES) to the Physically Removed System
c     (PRS). The eligible phases exclude the aqueous solution phase
c     and any fictive fugacity-fixing phases. The ES phase membership
c     and the corresponding finite-difference representations are
c     changed accordingly.
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
      integer ietmax,jetmax,kmax,nbtmax,netmax,nodbmx,nordmx,npetmx,
     $ nptmax,nsetmx,nstmax,nstsmx
c
      integer noutpt,nttyo
c
      integer iemop(npetmx),iemos(nsetmx),iindx0(kmax),iindx1(kmax),
     $ iodb(nodbmx),ipndx0(kmax),ipndx1(kmax),jcsort(nstmax),
     $ jern1(jetmax,netmax),jgext(netmax),jpflag(nptmax),jsflag(nstmax),
     $ nbasp(nbtmax),nbaspd(nbtmax),ncmpe(2,npetmx),ncmpr(2,nptmax),
     $ ngext(jetmax,netmax),nsts(nstsmx),nstsr(2,nstmax)
c
      integer iern1,iern2,imrn1,imrn2,ixrn1,ixrn2,kbt,kdim,kdim0,
     $ km1,km10,kmt,kmt0,kord,kstep,kx1,kx10,kxt,kxt0,nbt,ndelay,
     $ npet,npet0,npt,npts,nset,nset0
c
      logical qbye
c
      character*48 uspec(nstmax),uzvec0(kmax),uzvec1(kmax)
      character*24 uphase(nptmax)
      character*24 uaqsln
      character*8 ufixf
c
      real*8 csts(nstsmx),demop0(nordmx,npetmx),demos0(nordmx,nsetmx),
     $ emop(npetmx),emop0(npetmx),emos(nsetmx),emos0(nsetmx),
     $ fdpe0(nordmx,npetmx),fdpem1(nordmx,npetmx),
     $ fdse0(nordmx,nsetmx),fdsem1(nordmx,nsetmx),loph(nptmax),
     $ losp(nstmax),moph(nptmax),mosp(nstmax),mprph(nptmax),
     $ mprsp(nstmax),mrgexs(ietmax,jetmax,netmax),mtb(nbtmax),
     $ mtb0(nbtmax),xbar(nstmax),xbarlg(nstmax),zvclg0(kmax),
     $ zvclg1(kmax),zvec0(kmax),zvec1(kmax)
c
      real*8 zklgmn,zklogl
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ier,j2,kcol,nmax,np,npe,nr1,nr2,ns,nse
c
      integer ilnobl
c
      logical qprflg,qshftd,qtotsh
c
c-----------------------------------------------------------------------
c
c     Do a total shift of any eligible ES phase whose mass has begun
c     to decrease.
c
      qtotsh = .true.
      qshftd = .false.
c
      do npe = 1,npet
        if (fdpe0(1,npe) .lt. 0.) then
          np = iemop(npe)
          if (uphase(np)(1:24) .ne. uaqsln(1:24)) then
            if (uphase(np)(1:5) .ne. ufixf(1:5)) then
              j2 = ilnobl(uphase(np))
              write (noutpt,1000) uphase(np)(1:j2)
 1000         format(/' --- ',a,' is no longer precipitating ---')
              if (iodb(1) .ge. 1) then
                write (noutpt,1010) emos(npe),fdpe0(1,npe)
 1010           format(/'   Mass= ',1pe12.5,' moles',
     $          /'   Two-point relative rate= ',1pe12.5,' mol/mol',/)
              endif
c
              call shftph(emop,emop0,emos,emos0,fdpe0,fdpem1,
     $        fdse0,fdsem1,iemop,iemos,iern1,iern2,ietmax,iindx1,
     $        imrn1,imrn2,ipndx1,ixrn1,ixrn2,jern1,jetmax,jgext,
     $        jpflag,jsflag,kbt,kmax,km1,kmt,kx1,kxt,loph,losp,moph,
     $        mosp,mprph,mprsp,mrgexs,nbtmax,ncmpe,ncmpr,netmax,
     $        ngext,nordmx,noutpt,np,npet,npetmx,nptmax,nsetmx,
     $        nstmax,nttyo,qshftd,qtotsh,uphase,xbar,xbarlg,zklgmn,
     $        zklogl,zvclg0,zvclg1,zvec0,zvec1)
c
              jpflag(np) = 0
              nr1 = ncmpr(1,np)
              nr2 = ncmpr(2,np)
              do ns = nr1,nr2
                if (jsflag(ns) .eq. -10) jsflag(ns) = 0
              enddo
            endif
          endif
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qshftd) then
        qbye = .true.
        kord = 0
        ndelay = 1
        npts = 1
        qprflg = iodb(10) .ge. 1
c
        call escalc(csts,iindx1,jcsort,kbt,kmax,moph,mosp,mtb,
     $  mtb0,nbaspd,nbt,nbtmax,ncmpr,noutpt,npt,nptmax,nstmax,nsts,
     $  nstsmx,nstsr,qprflg,uspec)
c
        call miidxz(ier,iindx1,ipndx1,jpflag,jsflag,kbt,kdim,
     $  kmax,km1,kmt,kx1,kxt,losp,ncmpr,noutpt,npt,nptmax,nstmax,
     $  nttyo,uspec,uzvec1,zvclg1,zvec1)
c
        if (ier .gt. 0) then
          write (noutpt,1020)
          write (nttyo,1020)
 1020     format(/' * Error - (EQ6/dumpdp) Programming error trap:',
     $    "Can't recover from",/7x,'having exceeded the maximum ',i4,
     $    ' elements of the iindx1 array. This should',/7x,'not be',
     $    ' possible in the present subroutine.')
          stop
        endif
c
        kdim0 = kdim
        kxt0 = kxt
        kmt0 = kmt
        km10 = km1
        kx10 = kx1
        do kcol = 1,kdim
          iindx0(kcol) = iindx1(kcol)
          ipndx0(kcol) = ipndx1(kcol)
          uzvec0(kcol) = uzvec1(kcol)
          zvclg0(kcol) = zvclg1(kcol)
          zvec0(kcol) = zvec1(kcol)
        enddo
c
c       The ES phase assemblage has changed. Re-set the index arrays
c       associated with finite-difference description of the numbers
c       of mole of phases and species in the Equilibrium System (ES).
c
        call iiemop(iemop,iemos,iindx1,ipndx1,jsflag,kdim,kmax,
     $  ncmpe,ncmpr,noutpt,npet,npetmx,npt,nptmax,nset,nsetmx,nstmax,
     $  nttyo,uaqsln,uspec,uphase)
c
        npet0 = npet
        nset0 = nset
c
        nmax = nordmx*npetmx
        call initaz(fdpe0,nmax)
        call initaz(fdpem1,nmax)
        call initaz(demop0,nmax)
c
        nmax = nordmx*nsetmx
        call initaz(fdse0,nmax)
        call initaz(fdsem1,nmax)
        call initaz(demos0,nmax)
c
        do npe = 1,npet
          np = iemop(npe)
          emop0(npe) = moph(np)
          emop(npe) = moph(np)
          nr1 = ncmpe(1,npe)
          nr2 = ncmpe(2,npe)
          do nse = nr1,nr2
            ns = iemos(nse)
            emos0(nse) = mosp(ns)
            emos(nse) = mosp(ns)
          enddo
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
c
      end
