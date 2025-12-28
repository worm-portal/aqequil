      subroutine abswpk(beta,cdrs,csts,efac,ibswx,iebal,iindx1,
     $ jcsort,jflag,jssort,kbt,kmax,mosp,narn1,narn2,nbasp,nbaspd,
     $ nbt,nbtmax,ndrs,ndrsmx,ndrsr,nelect,nhydr,no2gaq,nstmax,
     $ nsts,nstsmx,nstsr,qbswx,q6mode,weight)
c
c     This subroutine determines the dominant species and the associated
c     factors required for continued fraction corrections.
c
c     This subroutine is called by:
c
c       EQLIB/absswa.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       beta   = array of normalized Newton-Raphson residual functions
c       q6mode = flag denoting usage for EQ3NR or EQ6:
c                  .false. = EQ3NR
c                  .true.  = EQ6
c
c     Principal output:
c
c       ibswx  = array of species indices of candidates for switching
c                  into the active basis set
c       qbswx  = flag denoting that candidates were found for
c                  basis switching
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer kmax,nbtmax,ndrsmx,nstmax,nstsmx
c
      integer ibswx(nbtmax),iindx1(kmax),jcsort(nstmax),jflag(nstmax),
     $ jssort(nstmax),nbasp(nbtmax),nbaspd(nbtmax),ndrs(ndrsmx),
     $ ndrsr(2,nstmax),nsts(nstsmx),nstsr(2,nstmax)
c
      integer iebal,kbt,narn1,narn2,nbt,nelect,nhydr,no2gaq
c
      logical qbswx,q6mode
c
      real(8) beta(kmax),cdrs(ndrsmx),csts(nstsmx),efac(nbtmax),
     $ mosp(nstmax),weight(nstmax)
c
      real(8) coefdr,coefst
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ix,iy,krow,nb,ns,nse,nsi,nsj
c
      logical qskip
c
      real(8) btest,cx,wsi
c
c-----------------------------------------------------------------------
c
c     Note: the following statements don't really do anything except
c     cause the compiler not to complain that nbt and jssort
c     are not used.
c
      ix = nbt
      iy = ix
      nbt = iy
c
      ix = jssort(1)
      iy = ix
      jssort(1) = iy
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set flag indicating that candidates for basis switching exist.
c
      qbswx = .false.
c
c     Clear the ibswx array.
c
      call initiz(ibswx,nbtmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Loop over all mass balance relations.
c
      do krow = 1,kbt
        nb = iindx1(krow)
        nse = nbaspd(nb)
        nsj = nbasp(nb)
c
        if (nse.ne.narn1 .and. nse.ne.nhydr .and. nb.ne.iebal .and.
     $    nse.ne.no2gaq .and.nse.ne.nelect) then
c
c         If the current species is not H2O, H+, O2(g,aq), or e-,
c         and is not being used for electrical balancing in EQ3NR,
c         consider a switch.
c
c         Pick a candidate for automatic basis switching if a
c         continued fraction correction would be large. Store
c         the index of the candidate, if any, in the array ibswx.
c
          btest = (beta(krow) + 1.)**efac(nb)
          if (btest .ge. 10.) then
c
            do ns = narn1,narn2
              weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
            enddo
c
c           Screen out species in the mass balance that are not
c           linked by current reaction with the current basis species.
c           Do this by setting the corresonding weights to zero.
c
            do ns = narn1,narn2
              if (jflag(ns) .eq. 30) then
c
c               Calling sequence substitutions:
c                 nsj for nse
c
                cx = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsj,ns,nstmax)
                if (cx .eq. 0.) weight(ns) = 0.
              endif
            enddo
c
            call fbassw(jcsort,jflag,mosp,narn1,narn2,nse,nsi,nsj,
     $      nstmax,weight,wsi)
c
            if (nsi .gt. 0) then
              ibswx(nb) = nsi
              qbswx = .true.
            endif
          endif
        endif
      enddo
c
      end
