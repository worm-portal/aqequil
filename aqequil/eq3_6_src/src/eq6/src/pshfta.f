      subroutine pshfta(csts,emop,emop0,emos,emos0,fdpe0,fdpem1,
     $ fdse0,fdsem1,iemop,iemos,iern1,iern2,ietmax,iindx1,imrn1,
     $ imrn2,iodb,ipndx1,ixrn1,ixrn2,jcsort,jern1,jetmax,jgext,
     $ jpflag,jsflag,kbt,km1,kmax,kmt,kx1,kxt,loph,losp,moph,mosp,
     $ mprph,mprsp,mrgexs,mtb,mtb0,nbasp,nbaspd,nbt,nbtmax,ncmpe,
     $ ncmpr,netmax,ngext,nodbmx,nordmx,noutpt,npet,npetmx,npt,
     $ nptmax,nsetmx,nstmax,nsts,nstsmx,nstsr,nttyo,uaqsln,ufixf,
     $ uphase,uspec,xbar,xbarlg,zklgmn,zklogl,zvclg0,zvclg1,
     $ zvec0,zvec1)
c
c     This subroutine oversees the partial transfer of eligible phases
c     (as a group) from the Equilibrium System (ES) to the Physically
c     Removed System (PRS). The eligible phases exclude the aqueous
c     solution phase and any fictive fugacity-fixing phases. Partial
c     transfer means that some mass of each phase transferred remains
c     in the ES.
c
c     This subroutine is called by:
c
c       EQ6/dumpdp.f
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
      integer iemop(npetmx),iemos(nsetmx),iindx1(kmax),iodb(nodbmx),
     $ ipndx1(kmax),jcsort(nstmax),jern1(jetmax,netmax),jgext(netmax),
     $ jpflag(nptmax),jsflag(nstmax),nbasp(nbtmax),nbaspd(nbtmax),
     $ ncmpe(2,npetmx),ncmpr(2,nptmax),ngext(jetmax,netmax),
     $ nsts(nstsmx),nstsr(2,nstmax)
c
      integer iern1,iern2,imrn1,imrn2,ixrn1,ixrn2,kbt,km1,kmt,
     $ kx1,kxt,nbt,npet,npt
c
      character*48 uspec(nstmax)
      character*24 uphase(nptmax)
      character*24 uaqsln
      character*8 ufixf
c
      real*8 csts(nstsmx),emop(npetmx),emop0(npetmx),emos(nsetmx),
     $ emos0(nsetmx),fdpe0(nordmx,npetmx),fdpem1(nordmx,npetmx),
     $ fdse0(nordmx,nsetmx),fdsem1(nordmx,nsetmx),loph(nptmax),
     $ losp(nstmax),moph(nptmax),mosp(nstmax),mprph(nptmax),
     $ mprsp(nstmax),mrgexs(ietmax,jetmax,netmax),mtb(nbtmax),
     $ mtb0(nbtmax),zvclg0(kmax),zvclg1(kmax),zvec0(kmax),
     $ zvec1(kmax),xbar(nstmax),xbarlg(nstmax)
c
      real*8 zklgmn,zklogl
c
c----------------------------------------------------------------------
c
c     Local variable declarations.
c
      logical qprflg,qshftd,qtotsh
c
      integer np,npe,nshftd
c
c----------------------------------------------------------------------
c
      nshftd = 0
      qtotsh = .false.
c
      write (noutpt,1000)
 1000 format(' --- Shifting (partial) ES solid(s) to the PRS ---',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Shift all ES phases except the aqueous solution phase and
c     any fictive fugacity-fixing phases.
c
      do npe = 1,npet
        np = iemop(npe)
        if (uphase(np)(1:24) .ne. uaqsln(1:24)) then
          if (uphase(np)(1:5) .ne. ufixf(1:5)) then
c
            call shftph(emop,emop0,emos,emos0,fdpe0,fdpem1,
     $      fdse0,fdsem1,iemop,iemos,iern1,iern2,ietmax,iindx1,
     $      imrn1,imrn2,ipndx1,ixrn1,ixrn2,jern1,jetmax,jgext,
     $      jpflag,jsflag,kbt,kmax,km1,kmt,kx1,kxt,loph,losp,moph,
     $      mosp,mprph,mprsp,mrgexs,nbtmax,ncmpe,ncmpr,netmax,
     $      ngext,nordmx,noutpt,np,npet,npetmx,nptmax,nsetmx,
     $      nstmax,nttyo,qshftd,qtotsh,uphase,xbar,xbarlg,zklgmn,
     $      zklogl,zvclg0,zvclg1,zvec0,zvec1)
c
            if (qshftd) nshftd = nshftd + 1
          endif
        endif
      enddo
c
      if (nshftd .gt. 0) then
c
c       Recompute the composition of the ES.
c
        qprflg = iodb(10) .ge. 1
        call escalc(csts,iindx1,jcsort,kbt,kmax,moph,mosp,mtb,
     $  mtb0,nbaspd,nbt,nbtmax,ncmpr,noutpt,npt,nptmax,nstmax,nsts,
     $  nstsmx,nstsr,qprflg,uspec)
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
