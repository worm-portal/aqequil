      subroutine ncmpvh(acflg,act,actlg,cdrs,cgxj,jflag,
     $ jsflag,losp,mosp,mtxj,nbasp,nbt,nbtmax,ndrs,ndrsmx,
     $ ndrsr,nrr1,nrr2,nstmax,xbar,xbarlg,xlks)
c
c     This subroutine calculates part of the "expansion" of the basis
c     set description regarding a site of a generic ion exchange
c     phase. It assists EQLIB/ncmpve.h. The part of the expansion
c     performed here does not take into account the need for an
c     iterative process to complete the actual expansion; rather, it
c     calculates a tentative part of the expansion. The iterative
c     process is required (at least in the general case, for example
c     Na+ for Ca++ following the Vanselow exchange model) because
c     mtxj (the sum of the number of moles of the species on the site)
c     is known only approximately at the start of the expansion process.
c     This subroutine uses a tentative value of mtxj as an input.
c
c     This subroutine is called by:
c
c       EQLIB/ncmpvh.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       acflg  = array of logarithms of activity coefficients
c       actlg  = array of logarithms of species activities
c                  (input consists of valid results for aqueous
c                  basis species; output consists of the same
c                  for exchanger species belonging to Vanselow
c                  exchanger phases)
c       cdrs   = array of reaction coefficients
c       mosp   = array of numbers of moles of species
c                  (input contains valid results for basis species
c                  and carry-over values for non-basis species)
c       mtxj   = the sum of the number of moles of exchanger species
c                  in the current site of the current generic ion
c                  exchange phase
c       ndrs   = array parallel to cdrs giving the index of the
c                  corresponding species
c       ndrsr  = array giving the range in the cdrs/ndrs arrays
c                  corresonding to the reaction associated with a
c                  given species
c       nrr1   = the start of the species index range for the
c                  current site of the current ion exchange phase
c       nrr2   = the end of the species index range for the
c                  current site of the current ion exchange phase
c       xlks   = array of equilibrium constants
c
c     Principal output:
c
c       act    = array of species activities
c                  (output consists of the subset for exchanger
c                  species belonging to Vanselow exchanger phases)
c       actlg  = array of logarithms of species activities
c       mosp   = array of numbers of moles of species
c                  (output consists of valid results for non-basis
c                  species; valid results for basis species were
c                  input and are retained)
c       losp   = array of logarithms of numbers of moles of species
c       xbar   = array of mole fractions of species
c       xbarlg = array of logarithms of mole fractions of species
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nbtmax,ndrsmx,nstmax
c
      integer jflag(nstmax),jsflag(nstmax),nbasp(nbtmax),ndrs(ndrsmx),
     $ ndrsr(2,nstmax)
c
      integer nbt,nrr1,nrr2
c
      real*8 acflg(nstmax),act(nstmax),actlg(nstmax),cdrs(ndrsmx),
     $ losp(nstmax),mosp(nstmax),xbar(nstmax),xbarlg(nstmax),
     $ xlks(nstmax)
c
      real*8 cgxj,mtxj
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n,nb,nr1,nr2,ns,nss
c
      integer nbasis
c
      real*8 cxs,lax,lx,lxx,mx,xx
c
      real*8 texp,tlg
c
c-----------------------------------------------------------------------
c
c     Estimate the mole fractions and activities of the
c     basis species for the current site.
c
      do ns = nrr1,nrr2
        nb = nbasis(nbasp,nbt,nbtmax,ns)
        if (nb.gt.0 .and. jsflag(ns).le.0) then
          xx = mosp(ns)/mtxj
          xbar(ns) = xx
          lxx = tlg(xx)
          xbarlg(ns) = lxx
          lax = cgxj*(lxx + acflg(ns))
          actlg(ns) = lax
          act(ns) = texp(lax)
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Estimate the mole fractions and activities of the
c     non-basis species for the current site.
c
      do ns = nrr1,nrr2
        if (jflag(ns).eq.30 .and. jsflag(ns).le.0) then
c
c         Calculate the activities and mole fractions for
c         the exchanger species not in the basis set.
c
          nr1 = ndrsr(1,ns)
          nr2 = ndrsr(2,ns)
          cxs = cdrs(nr1)
          lax = -xlks(ns)
          do n = nr1 + 1,nr2
            nss = ndrs(n)
            lax = lax + cdrs(n)*actlg(nss)
          enddo
          lax = -lax/cxs
          actlg(ns) = lax
          act(ns) = texp(lax)
          lxx = (lax/cgxj) - acflg(ns)
          xbarlg(ns) = lxx
          xx = texp(lxx)
          xbar(ns) = xx
c
c         Also calculate the corresponding numbers of moles.
c
          mx = xx*mtxj
          lx = tlg(mx)
          mosp(ns) = mx
          losp(ns) = lx
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
