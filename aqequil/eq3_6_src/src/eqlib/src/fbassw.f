      subroutine fbassw(jcsort,jflag,mosp,narn1,narn2,nse,nsi,nsj,
     $ nstmax,weight,wsi)
c
c     This subroutine attempts to find a candidate species to switch
c     switch into the active basis set. The function here is similar
c     to that EQLIB/fdomsp.f, which finds the species that dominates
c     a mass balance.
c
c     As presently contructed, this routine considers switches only
c     among the set of aqueous species. This condition may be relaxed
c     in the future.
c
c     This subroutine is called by:
c
c       EQLIB/abswpk.f
c       EQ6/absswb.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       jcsort = array of aqueous species indices in order of
c                  increasing concentration
c       jflag  = jflag array defining constraint types imposed on
c                  aqueous species
c       mosp   = array of numbers of moles of species
c       narn1  = index of the first aqueous species
c       narn2  = index of the last aqueous species
c       nse    = the data file basis species defining the current
c                mass balance
c       nsj    = the active basis species for the same mass balance;
c                nsj may be nse
c       weight = stoichiometric weighting factor
c
c     Principal output:
c
c       nsi    = index of the aqeuous species that would represent the
c                  optimal basis switch
c       wsi    = weight(nsi)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nstmax
c
      integer jcsort(nstmax),jflag(nstmax)
c
      integer narn1,narn2,nse,nsi,nsj
c
      real(8) mosp(nstmax),weight(nstmax)
c
      real(8) wsi
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ns,nsc,nss
c
      real(8) ap0,ap1,m0,m1,p0,p1,rx,w0,w1
c
c-----------------------------------------------------------------------
c
c     Note that nsi = 0 is equivalent to nsi = nsj (the current
c     active basis species associated with the present mass balance).
c     Returning no candidate means keeping the existing basis species.
c
      nsi = 0
      wsi = 0.
c
c     Make sure that the nse-th species is an aqueous species.
c
      if (nse.lt.narn1 .or. nse.gt.narn2) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Start by testing against the current active basis species.
c
      m0 = mosp(nsj)
      w0 = weight(nsj)
      p0 = w0*m0
      ap0 = abs(p0)
c
c     Loop over all aqueous species, testing in order of decreasing
c     number of moles. The current active basis species is not tested
c     against itself in this loop because jflag(nsj) is not 30.
c     Therefore, the following code cannot return nsi equal to nsj. If
c     the corresponding data file basis species is not the current
c     active basis species for the present mass balance, then jflag(nse)
c     will be 30, and this species will be tested as a potential
c     candidate.
c
      do nsc = narn1,narn2
        nss = narn2 + narn1 - nsc
        ns = jcsort(nss)
        w1 = weight(ns)
        m1 = mosp(ns)
        if (m1 .le. 0.) go to 110
        if (w1 .ne. 0.) then
          if (jflag(ns) .eq. 30) then
            p1 = w1*m1
            ap1 = abs(p1)
            if (ap1 .gt. ap0) then
c
c             Have found a new leading candidate.
c
              m0 = m1
              w0 = w1
              p0 = p1
              ap0 = ap1
              nsi = ns
              go to 100
            endif
          endif
        endif
c
c       Stop the search if the mole number ratio is now so high that
c       it becomes unreasonable to expect the ratio of the weighting
c       factors to possibly overcome this so as to yield ap1 > ap0.
c       Note that the possibility of m1 being zero has been tested
c       above, so the ratio calculation below should be safe.
c
        rx = m0/m1
        if (rx .gt. 100.) go to 110
  100   continue
      enddo
c
  110 wsi = w0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
