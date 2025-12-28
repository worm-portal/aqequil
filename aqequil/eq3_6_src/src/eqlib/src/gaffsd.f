      subroutine gaffsd(actlg,afcnst,affpd,affsd,cdrsd,jflagd,jpflag,
     $ ncmpr,ndrsd,ndrsmx,ndrsrd,npt,nptmax,nst,nstmax,qxknph,sidrph,
     $ sidrsp,uphase,uspec,xbar,xlksd)
c
c     This subroutine computes affinities and saturation indices based
c     on reactions in the 'd' set (cdrsd/ndrsd/ndrsrd arrays).
c
c     This subroutine is called by:
c
c       EQ3NR/scripx.f
c       EQ6/cdappl.f
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
c                  in the 'd' set
c       qxknph = flag indicating if the composition of the phase is
c                  known; this is the composition which maximizes the
c                  affinity (the actual composition if the phase is
c                  in equilibrium with the aqueous solution)
c
c     Principal output:
c
c       affsd  = affinities of reactions in the 'd' set for the
c                  dissociation or dissolution of species
c       affpd  = affinities for the dissolution of phases, based on
c                  reactions in the 'd' set
c       sidrsp = saturation indices (log Q/K) of reactions in the 'd'
c                  set for the dissociation or dissolution of species
c       sidrph = saturation indices of phases, based on reactions in
c                  the 'd' set
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ndrsmx,nptmax,nstmax
c
      integer jflagd(nstmax),jpflag(nptmax),
     $ ncmpr(2,nptmax),ndrsd(ndrsmx),ndrsrd(2,nstmax)
      integer npt,nst
c
      logical qxknph(nptmax)
c
      character(len=48) uspec(nstmax)
      character(len=24) uphase(nptmax)
c
      real(8) actlg(nstmax),affpd(nptmax),affsd(nstmax),cdrsd(ndrsmx),
     $ sidrph(nptmax),sidrsp(nstmax),xbar(nstmax),xlksd(nstmax)
      real(8) afcnst
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer np,nr1,nr2,ns,nt
c
c-----------------------------------------------------------------------
c
c     Compute affinities and saturation indices for reactions in the
c     'd' set.
c
      call gafsir(actlg,afcnst,affsd,cdrsd,jflagd,ndrsd,ndrsmx,
     $ ndrsrd,nst,nstmax,sidrsp,uspec,xlksd)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute affinities and saturation indices for phases.
c
      do np = 1,npt
        affpd(np) = 0.
        sidrph(np) = 0.
        if (jpflag(np) .lt. 2) then
          nr1 = ncmpr(1,np)
          nr2 = ncmpr(2,np)
          nt = nr2 - nr1 + 1
          if (nt .eq. 1) then
            affpd(np) = affsd(nr1)
            sidrph(np) = sidrsp(nr1)
          else
            if (qxknph(np)) then
c
c             The composition of the phase is known.
c
              do ns = nr1,nr2
                sidrph(np) = sidrph(np) + xbar(ns)*sidrsp(ns)
              enddo
              affpd(np) = afcnst*sidrph(np)
            else
c
c             The composition of the phase is not known.
c
              affpd(np) = -9999999.
              sidrph(np) = -9999999.
            endif
          endif
        endif
      enddo
c
      end
