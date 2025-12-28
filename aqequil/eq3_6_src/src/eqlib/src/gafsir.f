      subroutine gafsir(actlg,afcnst,affsd,cdrsd,jflagd,ndrsd,ndrsmx,
     $ ndrsrd,nst,nstmax,sidrsp,uspec,xlksd)
c
c     This subroutine computes the saturation indices of the reactions
c     for the destruction of all species. The saturation index is
c     defined as SI = log Q/K, where Q is the activity product and K is
c     the equilibrium constant of the reaction.
c
c     This subroutine is very similar in function to EQLIB/afcalc.f.
c     They both compute affinties and saturation indices. However, the
c     present subroutine does this for all species; the other subroutine
c     does it for a single specified reaction. The former could be
c     written so that it calls the latter inside a loop; however, this
c     is not done for the sake of avoiding the overhead of repeatedly
c     making such a call.
c
c     This subroutine is called by:
c
c       EQLIB/gaffsd.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       actlg  = array of logarithms of thermodynamic activities of
c                  species
c       cdrsd  = coefficients of reactions in the 'd' set
c       ndrsd  = indices of species appearing in reactions in the
c                  'd' set
c       ndrsrd = range pointer array for the ndrsd array
c       xbar   = mole fraction of a species in the phase to which
c                  it belongs
c       xlksd  = logarithms of equilibrium constants for reactions
c                in the 'd' set
c
c     Principal output:
c
c       affsd  = affinities of reactions in the 'd' set for the
c                  dissociation or dissolution of species
c       sidrsp = saturation indices (log Q/K) of reactions in the 'd'
c                  set for the dissociation or dissolution of species
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
      integer jflagd(nstmax),ndrsd(ndrsmx),ndrsrd(2,nstmax)
      integer nst
c
      character(len=48) uspec(nstmax)
c
      real(8) actlg(nstmax),affsd(nstmax),cdrsd(ndrsmx),sidrsp(nstmax),
     $ xlksd(nstmax)
      real(8) afcnst
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n,nr1,nr2,ns,nse,nt
c
      real(8) sx
c
c-----------------------------------------------------------------------
c
c     Zero the sidrsp and affsd arrays. This sets the correct values
c     for species which have jflagd values of 30 and for species whose
c     reactions are identity reactions.
c
      do ns = 1,nst
        sidrsp(ns) = 0.
        affsd(ns) = 0.
      enddo
c
      do ns = 1,nst
        if (jflagd(ns).ne.30 .and. jflagd(ns).ne.27) then
c
          nr1 = ndrsrd(1,ns)
          nr2 = ndrsrd(2,ns)
          nt = nr2 - nr1 + 1
          if (nt .lt. 2) go to 120
c
          if (xlksd(ns) .le. -9999999.) then
            sidrsp(ns) = 9999999.
            affsd(ns) = 9999999.
          elseif (xlksd(ns) .ge. 9999999.) then
            sidrsp(ns) = -9999999.
            affsd(ns) = -9999999.
          else
            sx = -xlksd(ns)
            do n = nr1,nr2
              nse = ndrsd(n)
              if (actlg(nse) .gt. -99999.) then
                sx = sx + cdrsd(n)*actlg(nse)
              else
                sx = -9999999.
                go to 110
              endif
            enddo
  110       continue
            sidrsp(ns) = sx
            affsd(ns) = afcnst*sx
          endif
c
        endif
c
  120   continue
      enddo
c
      end
