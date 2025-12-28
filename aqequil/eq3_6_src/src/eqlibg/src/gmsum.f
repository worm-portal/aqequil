      subroutine gmsum(conc,musumw,nmut,nmutmx,nmux,nstmax,pmu,uspec)
c
c     This subroutine computes the following triple sum used in Pitzer's
c     equations:
c
c       SUM(ijk) mu(ijk)*m(i)*m(j)*m(k)
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
c       conc   = array of concentration values
c       nmu    = number of mu values
c       nmux   = array of triples of species indices
c       pmu    = array of mu values
c
c     Principal output:
c
c       musumw =  the sum: SUM(ijk) mu(ijk)*m(i)*m(j)*m(k)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nmutmx,nstmax
c
      integer nmux(3,nmutmx)
      integer nmut
c
      character(len=48) uspec(nstmax)
c
      real(8) conc(nstmax),pmu(nmutmx)
      real(8) musumw
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer nmu,ns1,ns2,ns3
c
      real(8) cx
c
c-----------------------------------------------------------------------
c
      musumw = 0.
      do nmu = 1,nmut
        ns1 = nmux(1,nmu)
        ns2 = nmux(2,nmu)
        ns3 = nmux(3,nmu)
        cx = 6.
        if (ns1.eq.ns2 .and. ns1.ne.ns3) cx = 3.
        if (ns2.eq.ns3 .and. ns1.ne.ns3) cx = 3.
        if (ns1.eq.ns2 .and. ns1.eq.ns3) cx = 1.
        musumw = musumw + cx*pmu(nmu)*conc(ns1)*conc(ns2)*conc(ns3)
      enddo
c
      end
