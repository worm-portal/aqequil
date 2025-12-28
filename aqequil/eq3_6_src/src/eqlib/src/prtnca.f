      subroutine prtnca(ctb,jflag,jsflag,mrmlra,mwtsp,narn1,narn2,
     $ nbasp,nbaspd,nbt,nbtmax,noutpt,nstmax,qrho,rho,uspec,wfh2o)
c
c     This subroutine computes and prints a table of the numerical
c     composition of the aqueous solution in terms of mass balance
c     totals for component (basis) species in the 'd' basis set.
c
c     This subroutine is called by:
c
c       EQ3NR/scripx.f
c       EQ6/scripz.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ctb    = array of molalities of the component (basis) species
c       rho    = the density of the aqueous solution
c       qrho   = .true. if a solution density value is available
c       uspec  = array of names of aqueous species
c
c     Principal output:
c
c       None
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nbtmax,nstmax
c
      integer noutpt
c
      integer jflag(nstmax),jsflag(nstmax),nbasp(nbtmax),nbaspd(nbtmax)
      integer narn1,narn2,nbt
c
      logical qrho
c
      character(len=48) uspec(nstmax)
c
      real(8) ctb(nbtmax),mwtsp(nstmax)
      real(8) mrmlra,rho,wfh2o
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer nb,ns1,ns2
c
      real(8) molrx,ppmv,ppmw
c
c-----------------------------------------------------------------------
c
      write (noutpt,1000)
 1000 format(/11x,'--- Numerical Composition of the',
     $ ' Aqueous Solution ---')
c
      if (qrho) then
c
c       Have density data, include mg/L and Molarity.
c
        write (noutpt,1010)
 1010   format(/3x,'Species',20x,'mg/L',7x,'mg/kg.sol',4x,
     $  'Molarity',5x,'Molality',/)
c
        do nb = 1,nbt
          ns1 = nbaspd(nb)
          ns2 = nbasp(nb)
          if (ns1.ge.narn1 .and. ns1.le.narn2) then
            if (jsflag(ns2).lt.2 .and. jflag(ns2).ne.30) then
              ppmw = 1000.*ctb(nb)*mwtsp(ns1)*wfh2o
              ppmv = ppmw*rho
              molrx = ctb(nb)*mrmlra
              write (noutpt,1020) uspec(ns1),ppmv,ppmw,molrx,ctb(nb)
 1020         format(1x,a24,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5,
     $        1x,1pe12.5)
            endif
          endif
        enddo
      else
c
c       Have no density data, so do not include mg/L and Molarity.
c
        write (noutpt,1030)
 1030   format(/3x,'Species',18x,'mg/kg.sol',4x,'Molality',/)
c
        do nb = 1,nbt
          ns1 = nbaspd(nb)
          ns2 = nbasp(nb)
          if (ns1.ge.narn1 .and. ns1.le.narn2) then
            if (jsflag(ns2).lt.2 .and. jflag(ns2).ne.30) then
              ppmw = 1000.*ctb(nb)*mwtsp(ns1)*wfh2o
              write (noutpt,1040) uspec(ns1),ppmw,ctb(nb)
 1040         format(1x,a24,1x,1pe12.5,1x,1pe12.5)
            endif
          endif
        enddo
      endif
c
      write (noutpt,1050)
 1050 format(/3x,'Some of the above data may not be physically',
     $ ' significant.',/)
c
      end
