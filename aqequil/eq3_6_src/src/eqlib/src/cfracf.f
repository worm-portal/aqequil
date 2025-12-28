      subroutine cfracf(cdrs,csts,efac,jcsort,jflag,jssort,kmax,mosp,
     $ narn1,narn2,nbasp,nbaspd,nbt,nbtmax,ndrs,ndrsmx,ndrsr,nern1,
     $ nern2,nfac,nst,nstmax,nsts,nstsmx,nstsr,q6mode,weight)
c
c     This subroutine determines the dominant species and the associated
c     exponents required for continued fraction corrections. The data
c     are returned without any filtering as to their applicability.
c     particularly when the calling prgram is EQ3NR. For example, data
c     may be returned for mass balances for H+ or O2(g,aq), when
c     no such balances are used.
c
c     This subroutine is called by:
c
c       EQ3NR/arrset.f
c       EQ6/optmzr.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       mosp   = array containing the number of moles variables for the
c                  species; the conc array may be substituted for this
c                  in some calls
c       jcsort = array of species indices in sorted order within phase
c                  ranges
c       jssort = array of species indices in sorted order, ignoring
c                  phase ranges
c       q6mode = flag denoting usage for EQ3NR or EQ6:
c                  .false. = EQ3NR
c                  .true.  = EQ6NR
c
c
c     Principal output:
c
c       efac   = array of exponents for continued fraction corrections
c       nfac   = array of indices of dominant species
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
      integer jcsort(nstmax),jflag(nstmax),jssort(nstmax),nbasp(nbtmax),
     $ nbaspd(nbtmax),ndrs(ndrsmx),ndrsr(2,nstmax),nfac(nbtmax),
     $ nsts(nstsmx),nstsr(2,nstmax)
c
      integer narn1,narn2,nbt,nern1,nern2,nst
c
      logical q6mode
c
      real*8 cdrs(ndrsmx),csts(nstsmx),efac(nbtmax),mosp(nstmax),
     $ weight(nstmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer nb,nr1,ns,nse,nsi,nss
c
      real*8 cx,stx,wsi
      real*8 coefdr,coefst
c
c-----------------------------------------------------------------------
c
c     Initialize the nfac and efac arrays.
c
      do nb = 1,nbt
        nfac(nb) = 0
        efac(nb) = 1.0
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Loop over all mass balances, used or not.
c
      do nb = 1,nbt
        nse = nbasp(nb)
c
        if (jflag(nse) .eq. 0) then
c
          do ns = 1,nst
            weight(ns) = 0.
          enddo
c
          do nss = narn1,narn2
            ns = jcsort(nss)
            weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
          enddo
c
          if (q6mode) then
            do nss = nern1,nern2
              ns = jcsort(nss)
              weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
            enddo
          endif
c
c         Find the species (nfac) that makes the largest contribution
c         to the mass balance. Get the exponent (efac) required for the
c         continued fraction correction.
c
          call fdomsp(jssort,mosp,nsi,nst,nstmax,weight,wsi)
c
c         Get the exponent (efac) required for the continued fraction
c         correction.
c
          nfac(nb) = nsi
          ns = nbaspd(nb)
          if (nse .eq. ns) then
            efac(nb) = 1./wsi
          else
            cx = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns,nstmax)
            nr1 = ndrsr(1,ns)
            stx = -cx/cdrs(nr1)
            efac(nb) = 1./(stx*wsi)
          endif
        endif
      enddo
c
      end
