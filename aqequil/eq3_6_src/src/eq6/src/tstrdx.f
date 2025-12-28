      subroutine tstrdx(cdrs,iodb,iopt,jflag,jsflag,narn1,narn2,
     $ ndrs,ndrsmx,ndrsr,nodbmx,noptmx,noutpt,nrdxsp,nstmax,
     $ qredox,uspec)
c
c     This subroutine determines if the chemical model to be
c     computed has a redox aspect. This will be determined to be
c     so if a species in the model has an associated reaction
c     that is a redox reaction and the species is not in the
c     active basis set.
c
c     An auxiliary basis species (say Oxalate-) that is in the
c     model but is included in the active basis set is effectively.
c     treated as detached from a corresponding strict basis species
c     (e.g., the concentration/activity of Oxalate- is not determined
c     by the assumption of equilibrium for its associated reaction,
c     which would link it to HCO3-). In effect, Oxalate- is treated
c     as being composed of a pseudo-element, and its presence in the
c     model does not require a redox variable.
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
      integer ndrsmx,nodbmx,noptmx,nstmax
c
      integer noutpt
c
      integer iodb(nodbmx),iopt(noptmx),jflag(nstmax),jsflag(nstmax),
     $ ndrs(ndrsmx),ndrsr(2,nstmax)
c
      integer narn1,narn2,nrdxsp
c
      logical qredox
c
      character(len=48) uspec(nstmax)
c
      real(8) cdrs(ndrsmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jlen,ns
c
      character(len=56) uspn56
c
      real(8) cx
c
      real(8) coefdr
c
c-----------------------------------------------------------------------
c
      qredox = .false.
c
      if (iopt(15) .le. 0) then
        do ns = narn1,narn2
c
          if (jflag(ns) .eq. 30) then
c
c           The species is not in the active basis set. It is a
c           "dependent" species whose concentration/activity is
c           computed assuming its associated reaction is in a state
c           of equilibrium.
c
            if (jsflag(ns) .lt. 2) then
c
c             The species not hard suppressed.
c
c             Calling sequence substitutions:
c               nrdxsp for nse
c
              cx = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nrdxsp,ns,nstmax)
              if (cx .ne. 0.) then
                qredox = .true.
                if (iodb(1) .ge. 1) then
c
c                 Calling sequence substitutions:
c                   uspec(ns) for unam48
c
                  call fmspnm(jlen,uspec(ns),uspn56)
                  write (noutpt,1000) uspn56(1:jlen)
 1000             format(/' * Note - (EQ6/tstrdx) The reaction',
     $            ' associated with the species',/7x,a,' is a redox',
     $            ' reaction.')
                  go to 999
                endif
              endif
            endif
          endif
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
