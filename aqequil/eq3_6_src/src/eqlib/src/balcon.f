      subroutine balcon(aamatr,aamgex,al10,cdrs,cjbasp,cnufac,
     $ conc,dmlge,eps100,ggmgex,iern1,ietmax,iimgex,iindx1,
     $ ixbasp,jcsort,jern1,jern2,jetmax,jjsort,jsitex,kbt,kmax,
     $ krow,narn1,narn2,nbasp,nbtmax,ndrs,ndrsmx,ndrsr,nern1,
     $ nern2,netmax,noutpt,nphasx,nstmax,nttyo,rhsgex,uspec,
     $ weight,xbar)
c
c     This subroutine computes a row of the EQ3NR Jacobian matrix
c     for one of the following:
c
c       Mass balance constraint
c       Charge balance (charge adjustment) constraint
c       Alkalinity balance constraint
c
c     This subroutine is not used to compute any part of the EQ6
c     Jacobian matrix.
c
c     This subroutine is called by:
c
c       EQLIB/matrix.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       conc   = array of species concentrations
c
c     Output:
c
c       aamatr = the Jacobian matrix
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer noutpt,nttyo
c
      integer ietmax,jetmax,kmax,nbtmax,ndrsmx,netmax,nstmax
c
      integer iimgex(ietmax),iindx1(kmax),ipvgex(ietmax),ixbasp(nbtmax),
     $ jcsort(nstmax),jern1(jetmax,netmax),jern2(jetmax,netmax),
     $ jjsort(nstmax),jsitex(nstmax),nbasp(nbtmax),ndrs(ndrsmx),
     $ ndrsr(2,nstmax),nphasx(nstmax)
c
      integer iern1,kbt,krow,narn1,narn2,nern1,nern2
c
      character*48 uspec(nstmax)
c
      real*8 aamatr(kmax,kmax),aamgex(ietmax,ietmax),cdrs(ndrsmx),
     $ cjbasp(nbtmax),cnufac(nstmax),conc(nstmax),dmlge(ietmax),
     $ ggmgex(ietmax,ietmax),rhsgex(ietmax),weight(nstmax),xbar(nstmax)
c
      real*8 al10,eps100
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer icol,idim,ier,irow,je,jlen1,jlen2,kcol,nb,nbb,
     $ ne,np,nrr1,nrr2,nr1,nsc,nsi,nsj,nss
c
      logical qej,qmj,qpr,qrhsgz,qwj
c
      character*56 usp156,usp256
c
      real*8 axfc,cgxj,cxcc,cxjc,dmlgej,sumx
c
      real*8 coefdr
c
c-----------------------------------------------------------------------
c
      data qpr    /.false./
c
c-----------------------------------------------------------------------
c
      nb = iindx1(krow)
      nsi = nbasp(nb)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nsi.gt.narn1 .and. nsi.le.narn2) then
c
c       Have an aqueous species mass balance.
c
c       Build columns 1 through kbt.
c
        do kcol = 1,kbt
          nbb = iindx1(kcol)
          nsj = nbasp(nbb)
          sumx = 0.
c
c         First loop over aqueous species.
c
          do nss = narn1,narn2
            nsc = jcsort(nss)
            if (weight(nsc) .ne. 0.) then
              if (cnufac(nsc) .ne. 0.) then
c
c               Calling sequence substitutions:
c                 nsj for nse
c                 nsc for ns
c
                cxjc = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsj,nsc,nstmax)
                sumx = sumx + weight(nsc)*cxjc*cnufac(nsc)
              endif
            endif
          enddo
c
c         Do not loop over generic ion exchanger species here.
c         They do not contribute to non-exchanger mass balances
c         in EQ3NR. This is by design. The situation is different
c         in EQ6.
c
c         Complete the calculation of the current matrix element.
c
          if (nsj .eq. narn1) then
            aamatr(krow,kcol) = al10*sumx
          else
            aamatr(krow,kcol) = al10*(weight(nsj)*conc(nsj) + sumx)
          endif
        enddo
c
      elseif (nsi.ge.nern1 .and. nsi.le.nern2) then
c
c       Have a generic ion exchanger species mass balance corresponding
c       to a given exchanger phase and site.
c
        cgxj = cjbasp(nb)
c
        je = jsitex(nsi)
        np = nphasx(nsi)
        ne = np - iern1 + 1
c
c       Build columns 1 through kbt.
c
        do kcol = 1,kbt
          nbb = iindx1(kcol)
          nsj = nbasp(nbb)
          qej = nsj.ge.nern1 .and. nsj.le.nern2
          qwj = nsj .eq. narn1
          qmj = nsj.gt.narn1 .and. nsj.le.narn2
c
c         Calculate the dmlge vector for the current site/exchanger
c         phase and the current basis species. In essence,
c
c           dmlge(i,j) is partial m(i)/partial log m(j),
c
c         where i denotes a non-basis species for the current site/
c         exchanger phase and j denotes the basis species for the
c         current column index kcol. Theoretically, dmlge(j,j) could
c         also be included, but it can be obtained from a
c         straightforward calculation (as dmlgej below), in contrast
c         to the others, which must be obtained simultaneously.
c
c         Do not loop over aqueous species as non-basis species.
c         They can not contribute to this kind of mass balance,
c         because they do not contain the exchanger ligand.
c
c         Loop over all relevant generic ion exchanger species.
c         Generate an entry only for an active non-basis species.
c         In the present context of EQ3NR, this includes only
c         those exchanger species belonging to the site and phase
c         associated with the current row. The do loop immediately
c         below is over all species in the current exchanger phase
c         (e.g., runs over all species on all sites). That is done
c         to permit the use of the jcsort array. The sorting in that
c         array is over the entire phase, not by individual site.
c         The weight array filters out contributions from species
c         not in the current site in all calculations except those
c         involving the idim variable and the iimgex array. There,
c         it is necessary to check the site to which a species
c         belongs.
c
          idim = 0
          nrr1 = jern1(je,ne)
          nrr2 = jern2(je,ne)
          do nss = nrr1,nrr2
            nsc = jjsort(nss)
            if (weight(nsc) .ne. 0.) then
              if (cnufac(nsc) .ne. 0.) then
c
c               Calling sequence substitutions:
c                 nsj for nse
c                 nsc for ns
c
                cxjc = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsj,nsc,nstmax)
                if (cxjc .ne. 0.) then
                  idim = idim + 1
                  iimgex(idim) = nsc
                endif
              endif
            endif
          enddo
c
          do irow = 1,idim
            nsc = iimgex(irow)
            nr1 = ndrsr(1,nsc)
            cxcc = cdrs(nr1)
c
c           Calling sequence substitutions:
c             nsj for nse
c             nsc for ns
c
            cxjc = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsj,nsc,nstmax)
c
            if (qej) then
c
c             The current column is for a generic ion exchanger
c             species.
c
              axfc = xbar(nsc)*(1.0 + (cxjc/cxcc))
              do icol = 1,idim
                aamgex(irow,icol) = -axfc
              enddo
              aamgex(irow,irow) = 1.0 + aamgex(irow,irow)
              rhsgex(irow) = al10*(-(cxjc/cxcc)*conc(nsc)
     $        + axfc*conc(nsj))
            elseif (qwj) then
c
c             The current column is for solvent water.
c
              axfc = xbar(nsc)
              do icol = 1,idim
                aamgex(irow,icol) = -axfc
              enddo
              aamgex(irow,irow) = 1.0 + aamgex(irow,irow)
              rhsgex(irow) = -al10*(cxjc/cxcc)*conc(nsc)
            elseif (qmj) then
c
c             The current column is for an aqueous solute species.
c
              axfc = xbar(nsc)*(1.0 + (cxjc/cxcc))
              do icol = 1,idim
                aamgex(irow,icol) = -axfc
              enddo
              aamgex(irow,irow) = 1.0 + aamgex(irow,irow)
              rhsgex(irow) = -al10*(cxjc/(cxcc*cgxj))*conc(nsc)
            else
c
c             Calling sequence substitutions:
c               jlen1 for jlen
c               uspec(nsi) for unam48
c               usp156 for uspn56
c
              call fmspnx(jlen1,uspec(nsi),usp156)
c
c             Calling sequence substitutions:
c               jlen2 for jlen
c               uspec(nsj) for unam48
c               usp256 for uspn56
c
              call fmspnx(jlen2,uspec(nsj),usp256)
c
              write (noutpt,1110) usp156(1:jlen1),usp256(1:jlen2)
              write (nttyo,1110) usp156(1:jlen1),usp256(1:jlen2)
 1110         format(/' * Error - (EQLIB/balcon) Programming',
     $        ' error trap: In calculating',/7x,'the Jacobian',
     $        ' matrix entry for the ',a,' row',/7x,' and the ',
     $        a,' column, the type of the column',/7x,' species',
     $        ' does not correspond to a programmed possibility.')
              stop
            endif
          enddo
c
c         Check for a zero rhsgex vector. Note that if idim is zero,
c         qrhsgz will remain set to .true.
c
          qrhsgz = .true.
          do irow = 1,idim
            if (rhsgex(irow) .ne. 0.) then
              qrhsgz = .false.
              go to 300
            endif
          enddo
  300     continue
c
          if (qrhsgz) then
c
c           Have a zero right-hand-side vector. Set all dmlge elements
c           to zero.
c
            do irow = 1,idim
              dmlge(irow) = 0.
            enddo
          else
c
c           Solve the matrix equation. If EQLIBU/msolvr.f can't solve
c           the matrix, it is because the matrix is either zero
c           (ier = 1) or non-zero, but computationally singular
c           (ier = 2).
c
c           Calling sequence substitutions:
c             aamgex for aamatr
c             dmlge for delvec
c             ggmgex for ggmatr
c             ipvgex for ipivot
c             idim for kdim
c             ietmax for kmax
c             rhsgex for rhsvec
c
            call msolvr(aamgex,dmlge,ggmgex,ier,ipvgex,idim,ietmax,
     $      noutpt,nttyo,qpr,rhsgex)
c
            if (ier .ne. 0) then
c
c             The matrix is zero or it is computationally singular.
c
c             Calling sequence substitutions:
c               jlen1 for jlen
c               uspec(nsi) for unam48
c               usp156 for uspn56
c
              call fmspnx(jlen1,uspec(nsi),usp156)
c
c             Calling sequence substitutions:
c               jlen2 for jlen
c               uspec(nsj) for unam48
c               usp256 for uspn56
c
              call fmspnx(jlen2,uspec(nsj),usp256)
c
              write (noutpt,1120) usp156(1:jlen1),usp256(1:jlen2)
              write (nttyo,1120) usp156(1:jlen1),usp256(1:jlen2)
 1120         format(/' * Error - (EQLIB/balcon) Have encountered a',
     $        ' zero or computationally',/7x,'singular matrix while',
     $        ' trying to calculate the Jacobian matrix entry',
     $        /7x,'for the ',a,' row and the ',a,' column.')
              stop
            endif
          endif
c
          if (qrhsgz) then
            sumx = 0.
          else
            sumx = 0.
            do irow = 1,idim
              nsc = iimgex(irow)
              sumx = sumx + weight(nsc)*dmlge(irow)
            enddo
          endif
c
          if (nsj .ne. narn1) then
            dmlgej = al10*conc(nsj)
            sumx = sumx + weight(nsj)*dmlgej
          endif
c
          aamatr(krow,kcol) = sumx
c
        enddo
      endif
c
      end
