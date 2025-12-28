      subroutine gsdsm(conc,dpslm,nslt,nsltmx,nslx,nstmax,
     $ spsum,spsump,uspec)
c
c     This subroutine computes the following second order sums used
c     in Pitzer's equations:
c
c       spsum:  SUM(jk) S-lambda'(jk)*m(j)*m(k)
c       spsump: SUM(jk) S-lambda''(jk)*m(j)*m(k)
c
c     Here jk always refers to a cation-anion pair. Note that:
c
c       spsum  = 2 * SUM(ca) B'(ca)m(a)m(c)
c       spsump = 2 * SUM(ca) B''(ca)m(a)m(c)
c
c     This subroutine is called by:
c
c       EQLIBG/gcoeff.f
c
c-----------------------------------------------------------------------
c
c      Principal input:
c
c       conc   = array of concentration values
c       nslt   = number of S-lambda functions
c       nslx   = array of S-lambda species index pairs
c       dpslm  = array of ionic strength derivatives of S-lambda
c                  functions
c
c     Principal output:
c
c       spsum  = the sum: SUM(jk) S-lambda'(jk)*m(j)*m(k)
c       spsump = the sum: SUM(jk) S-lambda''(jk)*m(j)*m(k)
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
      real(8) conc(nstmax),dpslm(2,nsltmx)
c
      real(8) spsum,spsump
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer nsl,ns2,ns3
c
      real(8) cp
c
c-----------------------------------------------------------------------
c
      spsum = 0.
      spsump = 0.
      do nsl = 1,nslt
        ns2 = nslx(1,nsl)
        ns3 = nslx(2,nsl)
        cp = conc(ns2)*conc(ns3)
        spsum = spsum + dpslm(1,nsl)*cp
        spsump = spsump + dpslm(2,nsl)*cp
      enddo
      spsum = 2.*spsum
      spsump = 2.*spsump
c
      end
