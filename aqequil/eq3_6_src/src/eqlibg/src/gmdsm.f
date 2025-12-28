      subroutine gmdsm(conc,musum,na,natmax,nmutmx,nmxi,nmxmax,
     $ nmxx,ns,nstmax,pmu,uspec)
c
c     This subroutine computes the following second order sum used
c     in Pitzer's equations:
c
c       SUM(jk) mu(ijk)*m(j)*m(k)
c
c     where i is represented by the aqueous species index na and
c     overall species index ns.
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
c       na     = aqueous species index (i refers to ns, where ns is
c                  the corresponding species index)
c       nmxi   = range pointer array into the nmxx array:
c                  nmxi(1,na) and nmxi(2,na) are the first and last
c                  values of the second subscript (nmx) of the nmxx
c                  array for entries pertaining to the na-th
c                  aqueous species (ns-th species)
c       nmxx   = pointer array:
c                  nmxx(1,nmx) = the species index of the second
c                  species in the nmu-th triplet, nmxx(2,nmx) is the
c                  species index of the third species, and
c                  nmxx(3,nmx) = nmu
c       ns     = species index corresponding to the aqueous species
c                  index na
c       pmu    = array of mu values
c       uspec  = array of species names
c
c     Principal output:
c
c       musum  = the sum: SUM(jk) mu(ijk)*m(j)*m(k)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer natmax,nmutmx,nmxmax,nstmax
c
      integer nmxi(2,natmax),nmxx(3,nmxmax)
      integer na,ns
c
      character(len=48) uspec(nstmax)
c
      real(8) conc(nstmax),pmu(nmutmx)
      real(8) musum
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ifx,ilx,nmu,nmx,ns2,ns3
c
      real(8) cx
c
c-----------------------------------------------------------------------
c
      musum = 0.
      ifx = nmxi(1,na)
      ilx = nmxi(2,na)
      if (ilx .ge. ifx) then
        do nmx = ifx,ilx
          ns2 = nmxx(1,nmx)
          ns3 = nmxx(2,nmx)
          nmu = nmxx(3,nmx)
          cx = 2.
          if (ns2 .eq. ns3) cx = 1.
          musum = musum + cx*pmu(nmu)*conc(ns2)*conc(ns3)
        enddo
      endif
c
      end
