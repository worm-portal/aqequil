      subroutine jflaux(jflag,nbaspd,nbtd,nbtmax,ndrsd,ndrsmx,
     $ ndrsrd,nstmax)
c
c     This subroutine sets jflag to -1 for auxiliary basis species that
c     can not appear in the model.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       jflag  = array of species control flags
c       nbaspd = array of indices of species in the data file basis
c                  set
c       nbtd   = number of species in the data file basis set
c
c
c     Principal output:
c
c       jflag  = array of species control flags
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
      integer jflag(nstmax),nbaspd(nbtmax),ndrsd(ndrsmx),
     $ ndrsrd(2,nstmax)
      integer nbtd
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n,nb,ncount,nr1,nr2,ns,nse
c
c-----------------------------------------------------------------------
c
c     Note that looping is required to accommodate multi-level
c     dependencies in the auxiliary basis set. If jflag is set
c     to -1 for any auxiliary basis species, another pass is required
c     to insure that jflag is also set to -1 for any higher-level
c     auxiliary basis species.
c
  100 ncount = 0
      do nb = 1,nbtd
        ns = nbaspd(nb)
        nr1 = ndrsrd(1,ns)
        nr2 = ndrsrd(2,ns)
        if (jflag(ns) .eq. 30) then
          do n = nr1 + 1,nr2
            nse = ndrsd(n)
            if (jflag(nse) .eq. -1) then
              jflag(ns) = -1
              ncount = ncount + 1
              go to 110
            endif
          enddo
        endif
  110   continue
      enddo
      if (ncount .gt. 0) go to 100
c
      end
