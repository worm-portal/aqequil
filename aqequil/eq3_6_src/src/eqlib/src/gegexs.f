      subroutine gegexs(cegexs,cgexj,egexjc,egexjf,egexs,iern1,
     $ iern2,ietmax,jern1,jetmax,jgext,moph,mosp,mrgexs,netmax,
     $ ngexsa,ngext,noutpt,nptmax,nstmax,nttyo,zchar,zgexj)
c
c     This subroutine computes the equivalent fractions (egexs) and
c     mole ratios (mrgexs) of exchanger species belonging to generic
c     ion exchange phases.
c
c     The equivalent fractions for the exchanger species are the same
c     as those for the corresponding exchange ions. The mole ratios
c     are the number of moles of exchanger species per mole of
c     the exchange phase. The number of moles of an exchanger species
c     may or may not equal the number of moles of the corresponding
c     exchange ion. Thus, these mole ratios do not necessarily give
c     the number of moles of the exchange ions per mole of the exchange
c     phase.
c
c     The calculated equivalent fractions are based on formally
c     declared exchange capacities (egexjf). Thus, they need not sum
c     to unity, though they will do so for some models (e.g., Gapon,
c     Vanselow). If an anion should exchange onto a negatively charged
c     site, it will have a negative equivalent fraction.
c
c     Subroutine EQLIB/gegexw.f computes the apparent "whole-phase"
c     equivalent fractions (egexw) of cations and anions on generic ion
c     exchangers.
c
c     This subroutine is called by:
c
c       EQLIB/ncmpex.f
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
c       egexs  = array of equivalent fractions, by site, based on
c                  formally declared exchange capacities
c       mrgexs = array of mole ratios
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ietmax,jetmax,netmax,nptmax,nstmax
c
      integer noutpt,nttyo
c
      integer jern1(jetmax,netmax),jgext(netmax),
     $ ngexsa(ietmax,jetmax,netmax),ngext(jetmax,netmax)
c
      integer iern1,iern2
c
      real*8 cegexs(ietmax,jetmax,netmax),cgexj(jetmax,netmax),
     $ egexjc(jetmax,netmax),egexjf(jetmax,netmax),
     $ egexs(ietmax,jetmax,netmax),mrgexs(ietmax,jetmax,netmax),
     $ moph(nptmax),mosp(nstmax),zchar(nstmax),zgexj(jetmax,netmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ie,je,ne,np,ns,nss
c
      real*8 ex,exf,exsc
c
c-----------------------------------------------------------------------
c
c     Compute the equivalent fractions (egexs) for the sites of
c     the generic ion exchangers. Note that these fractions are
c     based on the formal declared exchange capacities (egexjf) of
c     these sites. These equivalent fractions will only sum to unity
c     for certain models, such as the Vanselow and Gapon models,
c     which require the actual capacity (egexjc) to equal the
c     corresponding formal capacity.
c
      do np = iern1,iern2
        ne = np - iern1 + 1
c
        do je = 1,jgext(ne)
c
c         First calculate the number of equivalents for each
c         component species from the corresponding number of moles.
c         Also calculate the total number of equivalents (egexjc)
c         actually on the current site of the current exchanger. This
c         may or may not equal the formal exchange capacity of the site
c         egexjf), depending on the exchange model. It will equal it in
c         the case of simple exchange models which allow exchange of
c         ions of only one charge sign on a given site. The Vanselow
c         and Gapon models satisfy this condition.
c
          exsc = 0.
          ns = jern1(je,ne) - 1
          do ie = 1,ngext(je,ne)
            ns = ns + 1
            ex = mosp(ns)*cegexs(ie,je,ne)
            egexs(ie,je,ne) = ex
            exsc = exsc + ex
          enddo
c
c         Note that the sign of the capacity is the opposite of that
c         of the charge on the site itself.
c
          egexjc(je,ne) = exsc/moph(np)
c
          if (egexjf(je,ne) .ne. 0.) then
c
c           Calculate the equivalent fractions.
c
            exf = egexjf(je,ne)*moph(np)
            do ie = 1,ngext(je,ne)
              egexs(ie,je,ne) = egexs(ie,je,ne)/exf
            enddo
          else
c
c           If the declared formal capacity is zero (not permitted
c           for the Vanselow and Gapon models), do not calculate
c           equivalent fractions. Set the corresponding elements of
c           the egexs array to zero.
c
            do ie = 1,ngext(je,ne)
              egexs(ie,je,ne) = 0.
            enddo
          endif
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the mole ratios (mrgexs) for exchanger species of
c     of generic ion exchanger reactants. The mole ratio is the
c     moles of exchanger species (specific to a site) per mole of
c     exchanger substrate.
c
      do np = iern1,iern2
        ne = np - iern1 + 1
        do je = 1,jgext(ne)
          ex = zgexj(je,ne)*cgexj(je,ne)
          do ie = 1,ngext(je,ne)
            nss = ngexsa(ie,je,ne)
            if (nss .gt. 0) then
              mrgexs(ie,je,ne) = -ex*egexs(ie,je,ne)/zchar(nss)
            endif
          enddo
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
