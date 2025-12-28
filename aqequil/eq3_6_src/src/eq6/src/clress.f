      subroutine clress(csts,iindx1,ipndx1,jpflag,jsflag,kdim,kmax,
     $ km1,kmt,kx1,kxt,loph,losp,moph,mosp,mtb,mtbaq,nbt,nbtmax,
     $ nptmax,nstmax,nsts,nstsmx,nstsr,ufixf,uzvec1,zvec1,zvclg1)
c
c     This subroutine clears equilibrium system (ES) solids. These
c     solids include pure minerals and solid solutions, but not fictive
c     fugacity-fixing minerals.
c
c     This subroutine is called by:
c
c       EQ6/eq6.f
c       EQ6/path.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       csts   = array of stoichiometric coefficients for mass
c                  balance expressions
c       losp   = array of log number of moles of species variables
c       mosp   = array of number of moles of species
c       mtb    = array of total numbers of moles of basis species
c       mtbaq  = array of total numbers of moles of basis species
c                  in the aqueous solution
c       uzvec1 = array of master variable names
c       zvec1  = array of master variables
c       zvclg1 = array of log master variables
c
c     Principal output:
c
c       losp   = array of log number of moles of species variables
c                  (modified)
c       mosp   = array of number of moles of species (modified)
c       mtb    = array of total numbers of moles of basis species
c                  (modified)
c       zvec1  = array of master variables (modified)
c       zvclg1 = array of log master variables (modified)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer kmax,nbtmax,nptmax,nstmax,nstsmx
c
      integer iindx1(kmax),ipndx1(kmax),jpflag(nptmax),jsflag(nstmax),
     $ nsts(nstsmx),nstsr(2,nstmax)
      integer kdim,km1,kmt,kx1,kxt,nbt
c
      character*48 uzvec1(kmax)
      character*8 ufixf
c
      real*8 csts(nstsmx),loph(nptmax),losp(nstmax),moph(nptmax),
     $ mosp(nstmax),mtb(nbtmax),mtbaq(nbtmax),zvec1(kmax),zvclg1(kmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer k,kcol,n,nb,nr1,nr2,np,ns
c
c-----------------------------------------------------------------------
c
c     Clear the phase and species status flags for the items to be
c     cleared from the ES.
c
      do kcol = km1,kmt
        if (uzvec1(kcol)(1:5) .ne. ufixf(1:5)) then
c
c         Do not have a fictive fugacity-fixing mineral.
c
          np = ipndx1(kcol)
          ns = iindx1(kcol)
          jpflag(np) = 0
          jsflag(ns) = 0
          moph(np) = 0.
          loph(np) = -99999.
          mosp(ns) = 0.
          losp(ns) = -99999.
        endif
      enddo
c
      do kcol = kx1,kxt
        np = ipndx1(kcol)
        ns = iindx1(kcol)
        jpflag(np) = 0
        jsflag(ns) = 0
        moph(np) = 0.
        loph(np) = -99999.
        mosp(ns) = 0.
        losp(ns) = -99999.
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Recompute the mass balance totals to account for the cleared
c     phases.
c
      do nb = 1,nbt
        mtb(nb) = mtbaq(nb)
      enddo
c
      do kcol = km1,kmt
        if (uzvec1(kcol)(1:5) .eq. ufixf(1:5)) then
c
c         Have a fictive fugacity-fixing mineral. Adjust the
c         mass balances to retain it in the system.
c         Retain it.
c
          ns = iindx1(kcol)
          nr1 = nstsr(1,ns)
          nr2 = nstsr(2,ns)
          do n = nr1,nr2
            nb = nsts(n)
            mtb(nb) = mtb(nb) + csts(n)*zvec1(kcol)
          enddo
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      k = km1 - 1
      do kcol = km1,kmt
        if (uzvec1(kcol)(1:5) .eq. ufixf(1:5)) then
c
c         Have a fictive fugacity-fixing mineral. Retain it in
c         the indexing scheme for the cleared system.
c
          k = k + 1
          zvec1(k) = zvec1(kcol)
          zvclg1(k) = zvclg1(kcol)
          ipndx1(k) = ipndx1(kcol)
          iindx1(k) = iindx1(kcol)
        endif
      enddo
c
      kmt = k
      kx1 = kmt + 1
c
      do kcol = kx1,kxt
c
c       Clear all other minerals, both pure minerals and solid
c       solutions.
c
        zvec1(kcol) = 0.
        zvclg1(kcol) = -99999.
        iindx1(kcol) = 0
        ipndx1(kcol) = 0
      enddo
c
      kxt = kmt
      kdim = kmt
c
      end
