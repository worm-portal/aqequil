      subroutine gssum(conc,dpslm,fxi,nslt,nsltmx,nslx,nstmax,pslm,
     $ ssumw,uspec)
c
c     This subroutine computes the following second order sum used
c     in Pitzer's equations:
c
c       SUM(ij) [ S-lambda(ij) + I*S-lambda'(ij) ]*m(i)*m(j)
c
c     This sum is used to calculate the activity coefficient
c     of water.
c
c     This subroutine is called by:
c
c       EQLIBG/gcoeff.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       conc   = array of species concentrations
c       fxi    = the ionic strength (the 2nd-order electrostatic
c                  moment function I)
c       nslt   = number of S-lambda functions
c       nslx   = array of S-lambda species index pairs
c       pslm   = array of S-lambda functions
c       dpslm  = array of ionic strength derivatives of S-lambda
c                functions
c
c     Principal output:
c
c       ssumw  = the sum:
c                  SUM(ij) [ S-lambda(ij) + I*S-lambda'(ij) ]*m(i)*m(j)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nsltmx,nstmax
c
      integer nslx(2,nsltmx)
      integer nslt
c
      character(len=48) uspec(nstmax)
c
      real(8) conc(nstmax),dpslm(2,nsltmx),pslm(nsltmx)
      real(8) fxi,ssumw
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer nsl,ns1,ns2
c
      real(8) cx
c
c-----------------------------------------------------------------------
c
      ssumw = 0.
      do nsl = 1,nslt
        ns1 = nslx(1,nsl)
        ns2 = nslx(2,nsl)
        cx = 2.
        if (ns1 .eq. ns2) cx = 1.
        ssumw = ssumw
     $  + cx*(pslm(nsl) + fxi*dpslm(1,nsl))*conc(ns1)*conc(ns2)
      enddo
c
      end
