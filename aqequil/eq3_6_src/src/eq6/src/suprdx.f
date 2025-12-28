      subroutine suprdx(jflag,jsflag,narn1a,narn2a,ndrsd,ndrsmx,
     $ ndrsrd,nrdxsp,nsta,nstmax,uspeca)
c
c     This subroutine executes the option to suppress all redox
c     reactions. This is not the same as suppressing all redox
c     species. An auxiliary basis species (say Oxalate-) that is
c     in the active basis set (and hence has its own mass balance
c     relation) is detached from any other species in the active
c     basis set, including the redox species. It should not be
c     suppressed as part of a general suppression of redox. For
c     example, if Oxalate- is given its own mass balance relation,
c     then equilibirum of its associated reaction, which is of redox
c     type linking it to HCO3-, would be overridden. Tn effect,
c     Oxalate- would be treated as being composed of a pseudo-element.
c
c     This subroutine is called by:
c
c       EQ6/eq6.f
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
      integer ndrsmx,nstmax
c
      integer jflag(nstmax),jsflag(nstmax),ndrsd(ndrsmx),
     $ ndrsrd(2,nstmax)
c
      integer narn1a,narn2a,nrdxsp,nsta
c
      character(len=48) uspeca(nstmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer kcount,n,ncount,nr1,nr2,ns,nse
c
c-----------------------------------------------------------------------
c
      ncount = 0
c
c     First, suppress any species whose associated reaction involves
c     the redox species, unless equilibrium for that reaction has been
c     over-ridden. Currently, such an over-ride is possible only for
c     certain aqueous species. For them, this condition is marked by
c     jflag is not 30.
c
      do ns = 1,narn1a - 1
c
c       Species preceding the block of aqueous species. There
c       should be none.
c
        if (jsflag(ns) .lt. 2) then
          nr1 = ndrsrd(1,ns)
          nr2 = ndrsrd(2,ns)
          do n = nr1 + 1,nr2
            nse = ndrsd(n)
            if (nse .eq. nrdxsp) then
              jsflag(ns) = 2
              ncount = ncount + 1
              go to 100
            endif
          enddo
  100     continue
        endif
      enddo
c
      do ns = narn1a,narn2a
c
c       Aqueous species.
c
        if (jflag(ns) .eq. 30) then
          if (jsflag(ns) .lt. 2) then
            nr1 = ndrsrd(1,ns)
            nr2 = ndrsrd(2,ns)
            do n = nr1 + 1,nr2
              nse = ndrsd(n)
              if (nse .eq. nrdxsp) then
                jsflag(ns) = 2
                ncount = ncount + 1
                go to 110
              endif
            enddo
  110       continue
          endif
        endif
      enddo
c
      do ns = narn2a + 1,nsta
c
c       Species following the block of aqueous species.
c
        if (jsflag(ns) .lt. 2) then
          nr1 = ndrsrd(1,ns)
          nr2 = ndrsrd(2,ns)
          do n = nr1 + 1,nr2
            nse = ndrsd(n)
            if (nse .eq. nrdxsp) then
              jsflag(ns) = 2
              ncount = ncount + 1
              go to 120
            endif
          enddo
  120     continue
        endif
      enddo
c
      if (ncount .le. 0) go to 999
c
c     Now check to insure that any species that are actively linked
c     to species that are now suppressed are themselves suppressed.
c
  210 kcount = 0
c
      do ns = 1,narn1a - 1
c
c       Species preceding the block of aqueous species. There
c       should be none.
c
        if (jsflag(ns) .lt. 2) then
          nr1 = ndrsrd(1,ns)
          nr2 = ndrsrd(2,ns)
          do n = nr1 + 1,nr2
            nse = ndrsd(n)
            if (jsflag(nse) .ge. 2) then
              jsflag(ns) = 2
              kcount = kcount + 1
              ncount = ncount + 1
              go to 220
            endif
          enddo
  220     continue
        endif
      enddo
c
      do ns = narn1a,narn2a
c
c       Aqueous species.
c
        if (jflag(ns) .eq. 30) then
          if (jsflag(ns) .lt. 2) then
            nr1 = ndrsrd(1,ns)
            nr2 = ndrsrd(2,ns)
            do n = nr1 + 1,nr2
              nse = ndrsd(n)
              if (jsflag(nse) .ge. 2) then
                jsflag(ns) = 2
                kcount = kcount + 1
                ncount = ncount + 1
                go to 230
              endif
            enddo
  230       continue
          endif
        endif
      enddo
c
      do ns = narn2a + 1,nsta
c
c       Species following the block of aqueous species.
c
        if (jsflag(ns) .lt. 2) then
          nr1 = ndrsrd(1,ns)
          nr2 = ndrsrd(2,ns)
          do n = nr1 + 1,nr2
            nse = ndrsd(n)
            if (jsflag(nse) .ge. 2) then
              jsflag(ns) = 2
              kcount = kcount + 1
              ncount = ncount + 1
              go to 240
            endif
          enddo
  240     continue
        endif
      enddo
c
      if (kcount .gt. 0) go to 210
c
  999 continue
c
      end
