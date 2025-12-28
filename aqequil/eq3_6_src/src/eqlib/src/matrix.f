      subroutine matrix(aamatr,al10,bpx,cdrs,cdrtw,cdrw,cjbasp,
     $ cnufac,conc,csts,dlogxw,eps100,ibpxmx,iebal,iern1,ietmax,
     $ iindx1,ipndx1,irdxc3,ixbasp,ixrn1,jcsort,jern1,jern2,jetmax,
     $ jflag,jjsort,jsitex,kbt,kction,kdim,kelect,khydr,kmax,kmt,
     $ km1,ko2gaq,kwater,kxt,kx1,mosp,narn1,narn2,nbasp,nbt,nbtmax,
     $ nbw,ncosp,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,netmax,
     $ noutpt,no2gaq,nphasx,nredox,nst,nstmax,nsts,nstsmx,nstsr,
     $ ntfx,ntfxmx,ntfxt,nttyo,nxtmax,omega,qredox,q6mode,tfx,
     $ ugexmo,uspec,weight,xbar,xbarw,xbarwc,zchar)
c
c     This subroutine computes the Jacobian matrix J[z] (aamatr).
c     This matrix is used by the algebraic equation solver.
c
c     This subroutine is called by:
c
c       EQLIB/nrstep.f
c       EQ6/optmzr.f
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
      integer ibpxmx,ietmax,jetmax,kmax,nbtmax,ndrsmx,netmax,nstmax,
     $ nstsmx,ntfxmx,nxtmax
c
      integer noutpt,nttyo
c
      integer iindx1(kmax),ipndx1(kmax),ixbasp(nbtmax),jcsort(nstmax),
     $ jern1(jetmax,netmax),jern2(jetmax,netmax),jflag(nstmax),
     $ jjsort(nstmax),jsitex(nstmax),kction(nbtmax),nbasp(nbtmax),
     $ ncosp(nbtmax),ndrs(ndrsmx),ndrsr(2,nstmax),nphasx(nstmax),
     $ nsts(nstsmx),nstsr(2,nstmax),ntfx(ntfxmx)
c
      integer iebal,iern1,irdxc3,ixrn1,kbt,kdim,kelect,khydr,kmt,km1,
     $ ko2gaq,kwater,kxt,kx1,narn1,narn2,nbt,nbw,nelect,nern1,nern2,
     $ no2gaq,nredox,nst,ntfxt
c
      logical qredox,qvansc,q6mode
c
      character*48 uspec(nstmax)
      character*24 ugexmo(netmax)
c
      real*8 aamatr(kmax,kmax),bpx(ibpxmx,nxtmax),cdrs(ndrsmx),
     $ cdrtw(nstmax),cdrw(nstmax),cjbasp(nbtmax),cnufac(nstmax),
     $ conc(nstmax),csts(nstsmx),dlogxw(nbtmax),mosp(nstmax),
     $ tfx(ntfxmx),weight(nstmax),xbar(nstmax),zchar(nstmax)
c
      real*8 al10,eps100,omega,xbarw,xbarwc
c
c-----------------------------------------------------------------------
c
c     Local pseudo-module.
c
c     This group contains scratch arrays for computing data (the
c     quantities 1/2.3026 d m(A)/d m(B), where A is a non-basis
c     generic ion exchanger species and B is a basis species. Here
c     x(B) replaces m(B) if B represents solvent water. These
c     quantities are required to compute the Jacobian matrix. The
c     quantities for all A sharing the same site of an exchanger
c     must be determined simultaneously. To do this, it is necessary
c     to solve a (usually small) matrix equation for each site.
c     If there are n non-basis exchanger species for a given site,
c     the matrix is n x n.
c
      integer isv_ietmax
c
      SAVE isv_ietmax
c
      integer, dimension(:), allocatable :: iimgex,ipvgex
c
      SAVE iimgex,ipvgex
c
      real(8), dimension(:), allocatable :: dmlge,rhsgex
c
      SAVE dmlge,rhsgex
c
      real(8), dimension(:,:), allocatable :: aamgex,ggmgex
c
      SAVE aamgex,ggmgex
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ix,iy,j,jlen,jfl,j2,kcol,krow,n,nb,nbb,nb2,ne,nmax,np,
     $ np1,nr1,ns,nsc,nsi,nsj,nss,ns1,ns2,nx
c
      integer ilnobl
c
      character*56 uspn56
c
      real*8 aax,aznsi,azns1,cvansf,cxcc,cxf,cxi1,cxjc,cxji,cxor,
     $ cx2r,cx21,stx,sumx,xx,zp
c
      real*8 coefdr,coefst
c
c-----------------------------------------------------------------------
c
c     Note: the following statements don't really do anything except
c     cause the compiler not to complain that kelect is not used.
c
      ix = kelect
      iy = ix
      kelect = iy
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (netmax .gt. 0) then
c
c       Have one or more generic ion exchange phases in the current
c       model. Allocate or reallocate scratch arrays as needed.
c
        if (.not.ALLOCATED(dmlge)) then
c
c         The arrays are not allocated. Zero the saved array size
c         variable. Note that only one array is tested to see if it
c         is allocated. It is assumed that all arrays in the set are
c         either allocated or not.
c
          isv_ietmax = 0
        else
c
c         The arrays are allocated. Check to see if the array size
c         variable hs changed. If so, deallocate the corresponding
c         arrays and zero the corresponding saved size variable.
c
          if (ietmax .ne. isv_ietmax) then
            DEALLOCATE (dmlge)
            DEALLOCATE (rhsgex)
            DEALLOCATE (aamgex,ggmgex)
            DEALLOCATE (iimgex,ipvgex)
            isv_ietmax = 0
          endif
        endif
c
c       At this point, the saved array size value is zero if the
c       corresponding arrays need to be allocated.
c
        if (isv_ietmax .eq. 0) then
          ALLOCATE (dmlge(ietmax))
          ALLOCATE (rhsgex(ietmax))
          ALLOCATE (aamgex(ietmax,ietmax),ggmgex(ietmax,ietmax))
          ALLOCATE (iimgex(ietmax),ipvgex(ietmax))
          isv_ietmax = ietmax
        endif
c
c       Zero these arrays.
c
        call initaz(dmlge,ietmax)
        call initaz(rhsgex,ietmax)
        nmax = ietmax*ietmax
        call initaz(aamgex,nmax)
        call initaz(ggmgex,nmax)
        call initiz(iimgex,ietmax)
        call initiz(ipvgex,ietmax)
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Zero the matrix aamatr.
c
      do i = 1,kdim
        do j = 1,kdim
          aamatr(i,j) = 0.
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Zero the cnufac array.
c
      do ns = 1,nst
        cnufac(ns) = 0.
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set up the cnufac array.
c
      do ns = narn1,narn2
        nr1 = ndrsr(1,ns)
        if (jflag(ns) .eq. 30) cnufac(ns) = -conc(ns)/cdrs(nr1)
      enddo
c
      do ns = nern1,nern2
        nr1 = ndrsr(1,ns)
        if (jflag(ns) .eq. 30) cnufac(ns) = -conc(ns)/cdrs(nr1)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (q6mode) go to 210
c
c     Write the matrix for EQ3NR.
c
c     Write rows 1 through kbt.
c
      do krow = 1,kbt
        nb = iindx1(krow)
        nsi = nbasp(nb)
        jfl = jflag(nsi)
        if (nsi.ge.narn1 .and. nsi.le.narn2) then
c
          if (krow.eq.kwater .and. jfl.eq.0) then
c
c           Water, mole fraction equation. Note that if this
c           block is executed, solvent water is a basis species,
c           log x(w) is a primary iteration variable, and that
c           cnufac(narn1) = 0.
c
            xx = xbarwc/omega
            do kcol = 1,kbt
              if (kcol .eq. kwater) then
c
c               Have the column for water itself.
c
                sumx = 0.
                do nss = narn1,narn2
                  nsc = jcsort(nss)
                  sumx = sumx + cdrw(nsc)*cnufac(nsc)
                enddo
                aamatr(krow,kcol) = -xx*sumx - 1.0
              else
                nbb = iindx1(kcol)
                nsj = nbasp(nbb)
                if (nsj.gt.narn1 .and. nsj.le.narn2) then
c
c                 Have a column for an aqueous solute basis species.
c
                  sumx = 0.
                  do nss = narn1,narn2
                    nsc = jcsort(nss)
c
c                   Calling sequence substitutions:
c                     nsc for ns
c                     nsj for nse
c
                    cxjc = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsj,nsc,
     $              nstmax)
                    sumx = sumx + cxjc*cnufac(nsc)
                  enddo
                  aamatr(krow,kcol) = -xx*(conc(nsj) + sumx)
                endif
              endif
            enddo
c
          elseif (nb .eq. iebal) then
c
c           Charge balance (charge adjustment).
c
            do ns = narn1,narn2
              weight(ns) = zchar(ns)
            enddo
c
            do ns = nern1,nern2
              weight(ns) = zchar(ns)
            enddo
c
            call balcon(aamatr,aamgex,al10,cdrs,cjbasp,cnufac,
     $      conc,dmlge,eps100,ggmgex,iern1,ietmax,iimgex,iindx1,
     $      ixbasp,jcsort,jern1,jern2,jetmax,jjsort,jsitex,kbt,kmax,
     $      krow,narn1,narn2,nbasp,nbtmax,ndrs,ndrsmx,ndrsr,nern1,
     $      nern2,netmax,noutpt,nphasx,nstmax,nttyo,rhsgex,uspec,
     $      weight,xbar)
c
          elseif (jfl .eq. 17) then
c
c           Log activity combination.
c
            ns1 = ncosp(nb)
            kcol = kction(nb)
            aamatr(krow,krow) = -1.
            aznsi = abs(zchar(nsi))
            azns1 = abs(zchar(ns1))
            zp = zchar(nsi)*zchar(ns1)
            if (zp .lt. 0.) then
              aamatr(krow,kcol) = -aznsi/azns1
            else
              aamatr(krow,kcol) = aznsi/azns1
            endif
c
          elseif (jfl .eq. 18) then
c
c           Log mean activity.
c
            ns1 = ncosp(nb)
            kcol = kction(nb)
            aamatr(krow,krow) = -1.
            aznsi = abs(zchar(nsi))
            azns1 = abs(zchar(ns1))
            aamatr(krow,kcol) = -aznsi/azns1
c
          elseif (jfl .eq. 21) then
c
c           pHCl.
c
            ns1 = ncosp(nb)
            kcol = kction(nb)
            aamatr(krow,krow) = -1.
            aznsi = abs(zchar(nsi))
            azns1 = abs(zchar(ns1))
            zp = zchar(nsi)*zchar(ns1)
            if (zp .lt. 0.) then
              aamatr(krow,kcol) = -aznsi/azns1
            else
              aamatr(krow,kcol) = aznsi/azns1
            endif
c
          elseif (jfl .eq. 25) then
c
c           Heterogeneous equilibrium.
c
            ns1 = ncosp(nb)
c
c           Calling sequence substitutions:
c             nsi for nse
c             ns1 for ns
c
            cxi1 = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsi,ns1,nstmax)
c
            do kcol = 1,kbt
              if (kcol.ne.krow .and. kcol.ne.kwater) then
                nb2 = iindx1(kcol)
                ns2 = nbasp(nb2)
c
c               Calling sequence substitutions:
c                 ns2 for nse
c
                cx21 = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns2,ns1,nstmax)
                aamatr(krow,kcol) = -cx21/cxi1
              endif
            enddo
            aamatr(krow,krow) = -1.0
c
          elseif (nsi .eq. no2gaq) then
c
c           Log fO2.
c
            if (irdxc3 .eq. 0) then
c
c             Identity.
c
              aamatr(krow,ko2gaq) = 1.
c
            elseif (irdxc3  .eq. -1) then
c
c             Eh residual (Note- if a pe- value was input, it has been
c             converted to an Eh value by EQ3NR/setup.f).
c
              aamatr(krow,kwater) = 2.0
              aamatr(krow,khydr) = aamatr(krow,khydr) - 4.0
              aamatr(krow,krow) = aamatr(krow,krow) - 1.0
c
            elseif (irdxc3  .eq. 1) then
c
c             Cross-linking (homogeneous aqueous redox) equilibrium.
c
              aamatr(krow,ko2gaq) = -1.
c
c             Calling sequence substitutions:
c               no2gaq for nse
c               nredox for ns
c
              cxor = coefdr(cdrs,ndrs,ndrsmx,ndrsr,no2gaq,nredox,nstmax)
              do kcol = 1,kbt
                if (kcol.ne.ko2gaq .and. kcol.ne.kwater) then
                  nb2 = iindx1(kcol)
                  ns2 = nbasp(nb2)
c
c                 Calling sequence substitutions:
c                   ns2 for nse
c                   nredox for ns
c
                  cx2r = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns2,nredox,
     $            nstmax)
                  aamatr(krow,kcol) = -cx2r/cxor
                endif
              enddo
            else
c
c             Error.
c
              write (noutpt,1000) irdxc3
              write (nttyo,1000) irdxc3
 1000         format(/' * Error - (EQLIB/matrix) Programming error',
     $        ' trap:',/7x,'Have encountered illegal irdxc3 value = ',
     $        i5,'.')
              stop
            endif
c
          elseif (jfl.ge.0 .and. jfl.le.3) then
c
c           Mass balance.
c
            do ns = narn1,narn2
              weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
            enddo
c
            do ns = nern1,nern2
              weight(ns) = 0.
            enddo
c
            call balcon(aamatr,aamgex,al10,cdrs,cjbasp,cnufac,
     $      conc,dmlge,eps100,ggmgex,iern1,ietmax,iimgex,iindx1,
     $      ixbasp,jcsort,jern1,jern2,jetmax,jjsort,jsitex,kbt,kmax,
     $      krow,narn1,narn2,nbasp,nbtmax,ndrs,ndrsmx,ndrsr,nern1,
     $      nern2,netmax,noutpt,nphasx,nstmax,nttyo,rhsgex,uspec,
     $      weight,xbar)
c
          elseif (jfl.ge.7 .and. jfl.le.11) then
c
c           Alkalinity balance.
c
c           Calling sequence substitutions.
c             weight for array
c             nstmax for nmax
c
            do ns = 1,nst
              weight(ns) = 0.
            enddo
c
            do n = 1,ntfxt
              ns = ntfx(n)
              weight(ns) = tfx(n)
            enddo
c
            call balcon(aamatr,aamgex,al10,cdrs,cjbasp,cnufac,
     $      conc,dmlge,eps100,ggmgex,iern1,ietmax,iimgex,iindx1,
     $      ixbasp,jcsort,jern1,jern2,jetmax,jjsort,jsitex,kbt,kmax,
     $      krow,narn1,narn2,nbasp,nbtmax,ndrs,ndrsmx,ndrsr,nern1,
     $      nern2,netmax,noutpt,nphasx,nstmax,nttyo,rhsgex,uspec,
     $      weight,xbar)
c
          elseif (jfl .eq. 16) then
c
c           Log activity.
c
            aamatr(krow,krow) = 1.0
c
          elseif (jfl.eq.19 .or. jfl.eq.20) then
c
c           pX (including pH).
c
            aamatr(krow,krow) = 1.0
c
          elseif (jfl.eq.22 .or. jfl.eq.23) then
c
c           pmX (including pmH).
c
            aamatr(krow,krow) = 1.0
c
          elseif (jfl .eq. 27) then
c
c           Aqueous homogeneous equilibrium.
c
            ns1 = nsi
c
c           Calling sequence substitutions:
c             nsi for nse
c             ns1 for ns
c
            cxi1 = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsi,ns1,nstmax)
            do kcol = 1,kbt
              if (kcol.ne.krow .and. kcol.ne.kwater) then
                nb2 = iindx1(kcol)
                ns2 = nbasp(nb2)
c
c               Calling sequence substitutions:
c                 ns2 for nse
c                 ns1 for ns
c
                cx21 = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns2,ns1,nstmax)
                aamatr(krow,kcol) = -cx21/cxi1
              endif
            enddo
            aamatr(krow,krow) = -1.0
c
          else
c
c           Have found a bad jflag value.
c
c           Calling sequence substitutions:
c             uspec(nsi) for unam48
c
            call fmspnx(jlen,uspec(nsi),uspn56)
            write (noutpt,1010) jfl,uspn56(1:jlen)
            write (nttyo,1010) jfl,uspn56(1:jlen)
 1010       format(/' * Error - (EQLIB/matrix) Programming error trap:',
     $      /7x,'Have encountered a bad jflag value of ',i4,' for',
     $      /7x,a,'.')
            stop
          endif
c
        elseif (nsi.ge.nern1 .and. nsi.le.nern2) then
c
c         Generic ion exchange species.
c
          if (jfl .eq. 0) then
c
c           Mass balance.
c
            do ns = narn1,narn2
              weight(ns) = 0.
            enddo
c
            do ns = nern1,nern2
              weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
            enddo
c
            call balcon(aamatr,aamgex,al10,cdrs,cjbasp,cnufac,
     $      conc,dmlge,eps100,ggmgex,iern1,ietmax,iimgex,iindx1,
     $      ixbasp,jcsort,jern1,jern2,jetmax,jjsort,jsitex,kbt,kmax,
     $      krow,narn1,narn2,nbasp,nbtmax,ndrs,ndrsmx,ndrsr,nern1,
     $      nern2,netmax,noutpt,nphasx,nstmax,nttyo,rhsgex,uspec,
     $      weight,xbar)
c
          else
c
c           Have found a bad jflag value.
c
c           Calling sequence substitutions:
c             uspec(nsi) for unam48
c
            call fmspnx(jlen,uspec(nsi),uspn56)
            write (noutpt,1010) jfl,uspn56(1:jlen)
            write (nttyo,1010) jfl,uspn56(1:jlen)
          endif
        endif
      enddo
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  210 continue
c
c     Write the matrix for EQ6.
c
c     Compute the dlogxw array (d log x(w)/d log m(s')).
c
      call gdlgxw(cdrs,cjbasp,cnufac,conc,dlogxw,eps100,
     $ ixbasp,jcsort,jflag,narn1,narn2,nbasp,nbt,nbtmax,nbw,
     $ ndrs,ndrsmx,ndrsr,nern1,nern2,noutpt,nstmax,nttyo,omega,
     $ xbar,xbarw)
c
c     Change the dlogxw array from derivatives with respect to
c     molalities to ones with respect to numbers of moles.
c     Currently, the definition is:
c
c       dlogxw = d log x(w)/d log m(s'), s' != w
c
c     The new definition will be:
c
c       dlogxw = d log x(w)/d log n(s'), s' includes w
c
c     There is now a derivative with respect to n(w). There was
c     none with respect to m(w).
c
c     Compute dlogxw(w) (d log x(w)/d log n(w)). For basis species
c     which are aqueous solute species, it happens that:
c
c       d log x(w)/d log n(s') = d log x(w)/d log m(s')
c
      sumx = 0.
      do nb = 1,nbt
        if (nb .ne. nbw) then
          ns = nbasp(nb)
          sumx = sumx + dlogxw(nb)
        endif
      enddo
      dlogxw(nbw) = -sumx
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Redefine the cnufac array in terms of numbers of moles in place of
c     concentrations.
c
      do ns = narn1,narn2
        nr1 = ndrsr(1,ns)
        if (jflag(ns) .eq. 30) cnufac(ns) = -mosp(ns)/cdrs(nr1)
      enddo
c
      do ns = nern1,narn2
        nr1 = ndrsr(1,ns)
        if (jflag(ns) .eq. 30) cnufac(ns) = -mosp(ns)/cdrs(nr1)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Build the rows corresponding to mass balance expressions for
c     the active basis species in the 'd' set.
c
      do krow = 1,kbt
        nb = iindx1(krow)
        nsi = nbasp(nb)
c
        if (nsi.eq.no2gaq .or. nsi.eq.nelect) then
          if (.not.qredox) then
            aamatr(krow,krow) = 1.0
            go to 250
          endif
        endif
c
        do ns = 1,nst
          weight(ns) = 0.
        enddo
c
c       Mass balance.
c
        do ns = narn1,narn2
          weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
        enddo
c
        do ns = nern1,nern2
          weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
        enddo
c
        do kcol = km1,kxt
          ns = iindx1(kcol)
          weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
        enddo
c
c       Build column kwater.
c
        aax = 0.
c
        do nss = narn1,narn2
          nsc = jcsort(nss)
          aax = aax - weight(nsc)*cnufac(nsc)
     $    *(cdrtw(nsc) + cdrw(nsc)*dlogxw(nbw))
        enddo
c
        aax = aax + weight(narn1)*mosp(narn1)
        aamatr(krow,kwater) = al10*aax
c
c       Build columns 1 to kbt, except kwater.
c
        do kcol = 1,kbt
          if (kcol .ne. kwater) then
            nbb = iindx1(kcol)
            nsj = nbasp(nbb)
c
c           Does the operational basis species associated with the
c           current column belong to a Vanselow ion exchange phase?
c
            qvansc = .false.
            if (nsj.ge.nern1 .and. nsj.le.nern2) then
              np = nphasx(nsj)
              ne = np - iern1 + 1
              j2 = ilnobl(ugexmo(ne))
              qvansc = ugexmo(ne)(1:j2).eq.'Vanselow' .or.
     $        ugexmo(ne)(1:9).eq.'Vanselow-'
            endif
c
            aax = 0.
c
            do nss = narn1,narn2
              nsc = jcsort(nss)
c
c             Calling sequence substitutions:
c               nsc for ns
c               nsj for nse
c
              cxjc = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsj,nsc,nstmax)
              aax = aax + weight(nsc)*cnufac(nsc)
     $        *(cxjc + cdrw(nsc)*dlogxw(nbb))
            enddo
c
cKK
            do nss = nern1,nern2
              nsc = jcsort(nss)
c
c             Calling sequence substitutions:
c               nsc for ns
c               nsj for nse
c
              cxjc = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsj,nsc,nstmax)
cKKKKKKKKK
c
              cxf = cxjc
              if (qvansc) then
c
c               The "Vanselow factor" cvansf calculated below is
c               a measure of the extent to which the mole fractions
c               of two exchanger species (the one, non-basis, the
c               other, basis) in a mass action equation for the
c               Vanselow ion exchange model can not be simply
c               replaced by corresponding molalities. This factor
c               is non-zero if the two exchanging ions have different
c               charges.
c
                nr1 = ndrsr(1,nsc)
                cxcc = cdrs(nr1)
                cvansf = -cxcc + (cxjc + cxcc)*xbar(nsj)
                cxf = cvansf
              endif
cKKKKKKKKK
              aax = aax + weight(nsc)*cnufac(nsc)
     $        *(cxf + cdrw(nsc)*dlogxw(nbb))
cOld          aax = aax + weight(nsc)*cnufac(nsc)
cOld $        *(cxjc + cdrw(nsc)*dlogxw(nbb))
            enddo
c
            aax = aax + weight(nsj)*mosp(nsj)
            aamatr(krow,kcol) = al10*aax
          endif
        enddo
c
c       Build the mineral columns, including those for solid solution
c       species.
c
        do kcol = km1,kxt
          nsj = iindx1(kcol)
          aamatr(krow,kcol) = al10*weight(nsj)*mosp(nsj)
        enddo
c
  250   continue
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Build the rows corresponding to mass action expressions
c     for pure minerals.
c
      do krow = km1,kmt
        nsi = iindx1(krow)
        do kcol = 1,kbt
          nbb = iindx1(kcol)
          nsj = nbasp(nbb)
          if (nsj .eq. narn1) then
            aamatr(krow,kcol) = -cdrtw(nsi) + cdrw(nsi)*dlogxw(nbb)
          else
c
c           Calling sequence substitutions:
c             nsi for ns
c             nsj for nse
c
            cxji = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsj,nsi,nstmax)
            aamatr(krow,kcol) = cxji + cdrw(nsi)*dlogxw(nbb)
          endif
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Build the rows corresponding to mass action expressions
c     for species belonging to solid solutions.
c
      if (kxt. ge. kx1) then
        do krow = kx1,kxt
          nsi = iindx1(krow)
          np = ipndx1(krow)
          nx = np - ixrn1 + 1
          stx = bpx(1,nx)
c
c         Build columns 1 through kbt.
c
          do kcol = 1,kbt
            nbb = iindx1(kcol)
            nsj = nbasp(nbb)
            if (nsj .eq. narn1) then
              aamatr(krow,kcol) = -cdrtw(nsi) + cdrw(nsi)*dlogxw(nbb)
            else
c
c           Calling sequence substitutions:
c             nsi for ns
c             nsj for nse
c
              cxji = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nsj,nsi,nstmax)
              aamatr(krow,kcol) = cxji + cdrw(nsi)*dlogxw(nbb)
            endif
          enddo
c
c         Build columns kx1 through kxt.
c
          nr1 = ndrsr(1,nsi)
          cxcc = cdrs(nr1)
          do kcol = kx1,kxt
            ns1 = iindx1(kcol)
            np1 = ipndx1(kcol)
            if (np1 .eq. np) then
              if (ns1 .eq. nsi) then
                aamatr(krow,kcol) = cxcc*stx*( 1. - xbar(nsi) )
              else
                aamatr(krow,kcol) = - cxcc*stx*xbar(ns1)
              endif
            endif
          enddo
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
