      subroutine gegexw(cegexs,egexpc,egexpa,egexw,iern1,iern2,
     $ ietmax,jern1,jetmax,jgext,kern1,kern2,ketmax,kgexsa,moph,
     $ mosp,netmax,ngexsa,ngext,noutpt,nptmax,nstmax,nttyo,
     $ xgexw,zchar)
c
c     This subroutine computes the apparent "whole-phase" equivalent
c     fractions (egexw) and mole fractions (xgexw) of the exchange
c     ions present in generic ion exchanger phases. The exchange
c     fractions are calculated separately for cations and anions.
c     For example, the exchange fraction of a cation is the number
c     of equivalents of that cation divided by the sum of the
c     equivalents of all cations in the same exchange phase. The
c     mole fractions are calculated in similar fashion.
c
c     Subroutine EQLIB/gegexs.f computes the equivalent fractions
c     (egexs) and mole ratios (mrgexs) of exchanger species (moles
c     of exchanger species per mole of exchanger phase). Those
c     equivalent fractions are defined for individual sites.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ6/path.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       moph   = array of numbers of moles of phases
c       mosp   = array of numbers of moles of species
c
c     Principal output:
c
c       egexw  = array of apparent "whole-phase" equivalent fractions,
c                  for exchange ions
c       egexpc = array of apparent cationic exchange capacities
c       egexpa = array of apparent anionic exchange capacities
c       xgexw  = array of apparent "whole-phase" mole fractions,
c                  for exchange ions
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ietmax,jetmax,ketmax,netmax,nptmax,nstmax
c
      integer noutpt,nttyo
c
      integer jern1(jetmax,netmax),jgext(netmax),kern1(netmax),
     $ kern2(netmax),kgexsa(ketmax,netmax),ngext(jetmax,netmax),
     $ ngexsa(ietmax,jetmax,netmax)
c
      integer iern1,iern2
c
      real*8 cegexs(ietmax,jetmax,netmax),egexpa(netmax),
     $ egexpc(netmax),egexw(ketmax,netmax),moph(nptmax),mosp(nstmax),
     $ xgexw(ketmax,netmax),zchar(nstmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ie,je,ke,ne,np,ns,nss
c
      real*8 ex,expac,expcc,mx,mxpac,mxpcc,zx
c
c-----------------------------------------------------------------------
c
c     Compute the apparent "whole-phase" equivalent fractions (egexw)
c     and mole fractions (xgexw) for the exchange species in generic
c     ion exchangers. Note that for a cationic exchange species
c     (e.g., Na+), this equivalent fraction is defined as the total
c     number of equivalents of the species on the exchanger (in all
c     sites, regardless of charge) divided by the sum of the total
c     number of equivalents of all such cationic species (expcc).
c     For an anionic species, this equivalent fraction is defined
c     analogously (expac being analogous to expcc). The formally
c     declared exchange capacities for the sites and the phases are
c     not used here.
c
      do np = iern1,iern2
        ne = np - iern1 + 1
c
c       Loop on exchanging species (e.g., Na+). Get the total
c       number of equivalents of each such species on the phase
c       by summing over contributions from all sites.
c
        do ke = kern1(ne),kern2(ne)
          egexw(ke,ne) = 0.
          nss = kgexsa(ke,ne)
          if (nss .gt. 0) then
c
c           Loop on sites.
c
            do je = 1,jgext(ne)
              ns = jern1(je,ne) - 1
c
c             Loop on exchanger species (e.g., Na-Z).
c
              do ie = 1,ngext(je,ne)
                ns = ns + 1
                if (nss .eq. ngexsa(ie,je,ne)) then
c
c                 The current exchanger species (e.g., Na-Z) is
c                 comprised of the current exchanging species
c                 (e.g., Na+).
c
                  ex = mosp(ns)*cegexs(ie,je,ne)
                  egexw(ke,ne) = egexw(ke,ne) + ex
                endif
              enddo
            enddo
          endif
        enddo
c
c       Get the apparent cationic and anionic "whole-phase" exchange
c       capacities (egexpc and egexpa, respectively).
c
        expcc = 0.
        expac= 0.
c
        do ke = kern1(ne),kern2(ne)
          nss = kgexsa(ke,ne)
          if (nss .gt. 0) then
            ex = egexw(ke,ne)
            if (ex .gt. 0.) then
              expcc = expcc + ex
            elseif (ex .lt. 0.) then
              expac = expac + ex
            endif
          endif
        enddo
        egexpc(ne) = expcc/moph(np)
        egexpa(ne) = expac/moph(np)
c
c       Calculate the number of moles of exchange ions in the exchange
c       phase from the corresponding number of equivalents. Get the
c       total number of moles of exchange cations (mxpcc) and of
c       exchange anions (mxpac).
c
        mxpcc = 0.
        mxpac= 0.
c
        do ke = kern1(ne),kern2(ne)
          nss = kgexsa(ke,ne)
          if (nss .gt. 0) then
            zx = zchar(nss)
            if (zx .ne. 0.) then
              mx = egexw(ke,ne)/zx
            else
              mx = 0.
            endif
            xgexw(ke,ne) = mx
            if (ex .gt. 0.) then
              mxpcc = mxpcc + mx
            elseif (ex .lt. 0.) then
              mxpac = mxpac + mx
            endif
          endif
        enddo
c
c       Loop on exchanging species (e.g., Na+). This time calculate
c       the apparent "whole-phase" equivalent fractions.
c
        do ke = kern1(ne),kern2(ne)
          nss = kgexsa(ke,ne)
          if (nss .gt. 0) then
            ex = egexw(ke,ne)
            if (ex .gt. 0.) then
              egexw(ke,ne) = ex/expcc
            elseif (ex .lt. 0.) then
              egexw(ke,ne) = ex/expac
            endif
          endif
        enddo
c
c       Loop on exchanging species (e.g., Na+). This time calculate
c       the apparent "whole-phase" mole fractions.
c
        do ke = kern1(ne),kern2(ne)
          nss = kgexsa(ke,ne)
          if (nss .gt. 0) then
            zx = zchar(nss)
            mx = xgexw(ke,ne)
            if (zx .gt. 0.) then
              xgexw(ke,ne) = mx/mxpcc
            elseif (zx .lt. 0.) then
              xgexw(ke,ne) = mx/mxpac
            endif
          endif
        enddo
c
      enddo
c
      end
