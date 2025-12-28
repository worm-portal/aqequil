      subroutine gcsts(cdrs,csts,jflag,nbaspd,nbt,nbtmax,ndrs,ndrsmx,
     $ ndrsr,noutpt,nsts,nstsmx,nstsr,nst,nstmax,nttyo,uspec)
c
c     This subroutine computes the stoichiometric factors which relate
c     each aqueous species to the nb-th member of the data file
c     basis set. These stoichiometric factors permit calculation of
c     sums which correspond to physically meaningful 'total' masses
c     or concentrations of the basis species, except for three of
c     these species. It is not possible to define physically
c     meaningful mass balances for for water, hydrogen ion, and the
c     aqueous species oxygen gas. The reactions input to this subroutine
c     may be rewritten from the data file forms to reflect the
c     elmination of one or more auxiliary basis species from the
c     active basis set. They may not be rewritten to reflect basis
c     switching, except for a switch which exchanges a strict basis
c     species with an auxiliary basis species.
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
c       cdrs   = array of reaction coefficients
c       nbaspd = indices of the species in the 'd' basis set
c       ndrs   = array of indices of the species corresponding to the
c                  coefficients in the cdrs array
c       ndrsr  = array giving the range in the cdrs and ndrs arrays
c                  corresponding to the reaction for a given species
c       uspec  = array of species names
c
c     Principal output:
c
c       csts   = array of stoichiometric coefficients appearing in
c                  mass balance relations
c       nsts   = array of indices of the basis species corresponding
c                  to the coefficients in the csts array
c       nstsr  = array giving the range in the csts and nsts arrays
c                  corresponding to a given species
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nbtmax,ndrsmx,nstsmx,nstmax
c
      integer noutpt,nttyo
c
      integer jflag(nstmax),nbaspd(nbtmax),ndrs(ndrsmx),ndrsr(2,nstmax),
     $ nsts(nstsmx),nstsr(2,nstmax)
      integer nbt,nst
c
      character*48 uspec(nstmax)
c
      real*8 cdrs(ndrsmx),csts(nstsmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jlen,jlene,n,nb,nbb,nerr,nj,nn,nr1,nr2,ns,nse,nt,nts
c
      integer nbasis
c
      logical qbasis
c
      character*56 uspe56,uspn56
c
      real*8 cxs
c
c-----------------------------------------------------------------------
c
      nerr = 0
      n = 0
c
      do ns = 1,nst
c
c       Set first element of the range pointer array.
c
        nstsr(1,ns) = n + 1
c
        nr1 = ndrsr(1,ns)
        nr2 = ndrsr(2,ns)
        nt = nr2 - nr1 + 1
        nts = nt - 1
c
c       Check to see if the current species is a member of the active
c       basis set. If so, set qbasis to .true. and nts to 1.
c
c       Calling sequence substitutions:
c         nbaspd for nbasp
c
        nb = nbasis(nbaspd,nbt,nbtmax,ns)
c
        qbasis = .false.
        if (nb.gt.0 .and. jflag(ns).ne.30) then
          qbasis = .true.
          nts = 1
        endif
c
c       Check to see if the csts/nsts array size is sufficient.
c
        nj = n + nts
        if (nj .gt. nstsmx) then
c
c         Calling sequence substitutions:
c           uspec(ns) for unam48
c
          call fmspnx(jlen,uspec(ns),uspn56)
          write (noutpt,1000) nstsmx,uspn56(1:jlen)
          write (nttyo,1000) nstsmx,uspn56(1:jlen)
 1000     format(/' * Error - (EQLIB/gcsts) The maximum ',i7,
     $    ' entries have been exceeded',/7x,'computing the csts',
     $    ' array of stoichiometric coefficients. The last',
     $    /7x,'species for which the coefficients were being computed',
     $    /7x,'was ',a,'. Increase the dimensioning',/7x,
     $    'parameter nstspa')
          stop
        endif
c
        if (qbasis) then
c
c         Set mass balance coefficient for the ns-th species if it is
c         in the active basis set.
c
          n = n + 1
          nsts(n) = nb
          csts(n) = 1.
        else
c
c         Set mass balance coefficients for other species appearing in
c         the reaction for the ns-th species if the ns-th species is
c         not in the active basis set.
c
          if (nt .lt. 2) then
c           Calling sequence substitutions:
c             uspec(ns) for unam48
c
            call fmspnx(jlen,uspec(ns),uspn56)
            write (noutpt,1005) uspn56(1:jlen)
            write (nttyo,1005) uspn56(1:jlen)
 1005       format(/' * Error - (EQLIB/gcsts) The species ',a,
     $      /7x,'has no associated reaction on the data file, but it',
     $      /7x,'is not a strict basis species.')
            nerr = nerr + 1
c
c           Make a single null entry.
c
            n = n + 1
            nsts(n) = 0
            csts(n) = 0.
          else
            cxs = -cdrs(nr1)
            do nn = nr1 + 1,nr2
              n = n + 1
              nse = ndrs(nn)
c
c             Calling sequence substitutions:
c               nbaspd for nbasp
c               nse for ns
c
              nbb = nbasis(nbaspd,nbt,nbtmax,nse)
              if (nbb .le. 0) then
c               Calling sequence substitutions:
c                 uspec(ns) for unam48
c
                call fmspnx(jlen,uspec(ns),uspn56)
c
c               Calling sequence substitutions:
c                 jlene for jlen
c                 uspec(nse) for unam48
c                 uspe56 for uspn56
c
                call fmspnx(jlene,uspec(nse),uspe56)
                write (noutpt,1010) uspe56(1:jlene),uspn56(1:jlen)
                write (nttyo,1010) uspe56(1:jlene),uspn56(1:jlen)
 1010           format(/' * Error - (EQLIB/gcsts) The species ',a,
     $          /7x,'appears in the data file reaction for ',a,',',
     $          /7x,'but it is not in the data file basis set.')
                nerr = nerr + 1
              endif
              nsts(n) = nbb
              csts(n) = cdrs(nn)/cxs
            enddo
          endif
        endif
c
c       Set the second element of the range pointer array.
c
        nstsr(2,ns) = n
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nerr .gt. 0) then
        write (noutpt,1015)
        write (nttyo,1015)
 1015   format(/' * Error - (EQLIB/gcsts) One or more reactions',
     $  /7x,'on the data file are not consistent with the data file',
     $  /7x,'basis set.')
        stop
      endif
c
      end
