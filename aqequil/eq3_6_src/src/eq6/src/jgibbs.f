      subroutine jgibbs(aamatr,afcnst,affp,cdrs,cscale,csts,delvec,
     $ eps100,gmmatr,iindx1,iodb,ipivot,ipndx1,jpflag,kbt,kdim,kmax,
     $ km1,kmt,kx1,kxt,mtb,nbasp,nbtmax,ndrs,ndrsmx,ndrsr,nodbmx,
     $ noutpt,npadd,npdel,nptmax,nstmax,nsts,nstsmx,nstsr,nttyo,
     $ rhsvec,uphase,uspec,xlks)
c
c     This subroutine determines which one of a set of phases should be
c     removed from the equilibrium system when there is a mineralogic
c     phase rule violation. It has not yet been extended to handle
c     solid solutions.
c
c     It is assumed heres that a mineralogical phase rule violation
c     does in fact exist. This subroutine does not make an independent
c     confirmation.
c
c     This subroutine is called by:
c
c       EQ6/eqcalc.f.
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       npadd  = index of the last phase added to the equilibrium
c                  phase assemblage
c
c     Principal output:
c
c       npdel = index of the phase which should be deleted from the
c                 equilibrium phase assemblage in order to resolve
c                 the violation
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
      integer kmax,nbtmax,ndrsmx,nodbmx,nptmax,nstmax,nstsmx
c
      integer noutpt,nttyo
c
      integer iindx1(kmax),iodb(nodbmx),ipivot(kmax),ipndx1(kmax),
     $ jpflag(nptmax),nbasp(nbtmax),ndrs(ndrsmx),ndrsr(2,nstmax),
     $ nsts(nstsmx),nstsr(2,nstmax)
c
      integer kbt,kdim,km1,kmt,kx1,kxt,npadd,npdel
c
      character*48 uspec(nstmax)
      character*24 uphase(nptmax)
c
      real*8 aamatr(kmax,kmax),affp(nptmax),cdrs(ndrsmx),
     $ cscale(nstmax),csts(nstsmx),delvec(kmax),gmmatr(kmax,kmax),
     $ mtb(nbtmax),rhsvec(kmax),xlks(nstmax)
c
      real*8 afcnst,eps100
c
c-----------------------------------------------------------------------
c
c     Local variable declarations with global dimensions.
c
      integer isv_kmax
c
      SAVE isv_kmax
c
      integer, dimension(:), allocatable :: kmprv,ksprv
c
      SAVE kmprv,ksprv
c
      real(8), dimension(:), allocatable :: cscx
c
      SAVE cscx
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,idim,ier,irow1,irow2,j,jcol1,jcol2,jdim,j2,j3,j4,k,
     $ kcol,km,kmm,kp,kpadd,kpp,kpt,n,nb,np,npp,nr1,nr2,ns,nse,nss
c
      integer ilnobl
c
      logical qldep,qpr
c
      real*8 acx,acx2,affc,affmin,afscal,cx,mxx,mxxe,mxxm
c
      real*8 coefdr
c
c-----------------------------------------------------------------------
c
      data qpr    /.false./
c
c-----------------------------------------------------------------------
c
c     Allocate or reallocate local work arrays as needed.
c
      if (.not.ALLOCATED(kmprv)) then
c
c       Local work arrays are not allocated. Zero the saved
c       array size variables. Note that only one array is tested
c       to see if it is allocated. It is assumed that all local
c       work arrays are either allocated or not.
c
        isv_kmax = 0
      else
c
c       Local work arrays are allocated. Check to see if any of the
c       array size variables have changed. If so, deallocate
c       the corresponding local work arrays and zero the corresponding
c       saved size variables.
c
        if (kmax .ne. isv_kmax) then
          DEALLOCATE(kmprv)
          DEALLOCATE(ksprv)
          DEALLOCATE(cscx)
          isv_kmax = 0
        endif
      endif
c
c     At this point, the saved array size values are zero if the
c     corresponding arrays need to be allocated.
c
      if (isv_kmax .eq. 0) then
        ALLOCATE(kmprv(kmax),ksprv(kmax))
        ALLOCATE(cscx(kmax))
        isv_kmax = kmax
      endif
c
c     Zero the contents of the local work arrays.
c
      do k = 1,kmax
        kmprv(k) = 0
        ksprv(k) = 0
      enddo
c
      do k = 1,kmax
        cscx(k) = 0.
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(1) .gt. 0) then
        if (kxt .ge. kx1) then
          write (noutpt,1000)
          write (nttyo,1000)
 1000     format(/' * Note - (EQ6/jgibbs) This subroutine finds the',
     $    ' phases or phases',/7x,'involved in a violation of the',
     $    ' mineralogic phase rule, but it',/7x,'presently considers',
     $    ' only pure minerals in the analysis.')
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Clear the Jacobian matrix for use as workspace.
c
      do i = 1,kdim
        do j = 1,kdim
          aamatr(i,j) = 0.
        enddo
      enddo
c
      do i = 1,kbt + 1
        kmprv(i) = 0
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Find the phases which are directly involved in the mineralogic
c     phase rule violation. For each possible phase, test to see if
c     the violation disappears when the phase is removed from the
c     tested assemblage.
c
      kpt = 0
      do km = km1,kmt
        ns = iindx1(km)
c
c       For each of the other minerals involved in the current
c       assemblage, set up a corresponding row which contains the
c       reaction coefficients of the basis species.
c
        idim = 0
        do kmm = km1,kmt
          if (kmm .ne. km) then
            nss = iindx1(kmm)
            idim = idim + 1
            do j = 1,kbt
              nb = iindx1(j)
              nse = nbasp(nb)
c
c             Calling sequence substitutions:
c               nss for ns
c
              cx = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,nss,nstmax)
              aamatr(idim,j) = cx
            enddo
          endif
        enddo
c
c       Test for linear dependence.
c
        irow1 = 1
        irow2 = idim
        jcol1 = 1
        jcol2 = kbt
        call lindep(aamatr,eps100,irow1,irow2,jcol1,jcol2,kmax,qldep)
c
        if (.not.qldep) then
c
c         Have not obtained linear dependence. The phase being tested
c         is involved in the violation. Note its matrix index in the
c         kmprv array.
c
          kpt = kpt + 1
          kmprv(kpt) = km
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (kpt .ge. 1) then
        write (noutpt,1050)
        write (nttyo,1050)
 1050   format(/' * Note - (EQ6/jgibbs) The phases involved in the',
     $  ' mineralogic',/7x,'phase rule violation are:',/)
c
        do kp = 1,kpt
          km = kmprv(kp)
          ns = iindx1(km)
          np = ipndx1(km)
          j2 = ilnobl(uphase(np))
          write (noutpt,1060) np,uphase(np)(1:j2)
          write (nttyo,1060) np,uphase(np)(1:j2)
 1060     format(2x,i5,2x,a)
        enddo
c
        if (kpt .eq. 1) then
          write (noutpt,1070)
          write (nttyo,1070)
 1070     format(/' * Error (EQ6/jgibbs) Programming error trap:',
     $    ' Have a supposed',/7x,'mineralogical phase rule violation',
     $    ' involving only one phase.')
          stop
        endif
      else
        write (noutpt,1080)
        write (nttyo,1080)
 1080   format(/' * Error (EQ6/jgibbs) Programming error trap:',
     $  ' Have a supposed',/7x,'mineralogical phase rule violation',
     $  ' involving no phases.')
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Now determine the affinity of each involved phase under the
c     assumption that the solution is saturated with respect to
c     the others.
c
      do kp = 1,kpt
        km = kmprv(kp)
        ns = iindx1(km)
        np = ipndx1(km)
        if (iodb(1) .gt. 0) then
c
          j2 = ilnobl(uphase(np))
          write (noutpt,1090) uphase(np)(1:j2)
 1090     format(/1x,'Hypothetical affinity calculation for ',a,':')
c
          write (noutpt,1100)
 1100     format(/3x,'Cognate minerals:',/)
          do kpp = 1,kpt
            if (kpp .ne. kp) then
              kmm = kmprv(kpp)
              npp = ipndx1(kmm)
              j3 = ilnobl(uphase(npp))
              write (noutpt,1110) uphase(npp)(1:j3)
 1110         format(5x,a)
            endif
          enddo
        endif
c
c       Once again, clear the Jacobian matrix for use as workspace.
c
        do i = 1,kdim
          do j = 1,kdim
            aamatr(i,j) = 0.
          enddo
        enddo
c
c       For each of the other minerals involved in the violation,
c       set up a corresponding row which contains the reaction
c       coefficients of the basis species. The objective is to
c       find a set of such coefficients which does not exhibit linear
c       dependence. Such dependence is possible because the number of
c       basis species appearing in these reactions may exceed the
c       number of "other" minerals. For example, if three ferric iron
c       minerals are present in the violating assemblage, Fe+++ is
c       a dependent species, and Fe++ is the strict basis species for
c       iron, then the reactions will be written in terms of O2(g), H+,
c       H2O, and Fe++ instead of in terms of Fe+++. Depending on what
c       is exactly in a given set, one or more of O2(g), H+, and H2O
c       may be redundant in terms of relating the current mineral to
c       to the current set of other minerals involved in the violation.
c
        do i = 1,kbt
          ksprv(i) = 0
        enddo
c
        idim = 0
        do kpp = 1,kpt
          if (kpp .ne. kp) then
            kmm = kmprv(kpp)
            nss = iindx1(kmm)
            idim = idim + 1
            do j = 1,kbt
              nb = iindx1(j)
              nse = nbasp(nb)
c
c             Calling sequence substitutions:
c               nss for ns
c
              cx = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,nss,nstmax)
              aamatr(idim,j) = cx
            enddo
          endif
        enddo
c
c       Factor the rows.
c
        irow1 = 1
        irow2 = idim
        jcol1 = 1
        jcol2 = kbt
        call lindep(aamatr,eps100,irow1,irow2,jcol1,jcol2,kmax,qldep)
c
        if (qldep) then
          j2 = ilnobl(uphase(np))
          write (noutpt,1120) uphase(np)(1:j2)
          write (nttyo,1120) uphase(np)(1:j2)
 1120     format(/' * Error (EQ6/jgibbs) Programming error trap:',
     $    " Couldn't find a",/7x,'set of independent basis species',
     $    ' reaction coefficients needed',/7x,'to calculate the',
     $    ' hypothetical affinity of ',a,',',/7x,'which is involved',
     $    ' in a violation of the mineralogic phase rule.',/7x,
     $    'Previous testing should guarantee the existence of such',
     $    ' a set.',/7x,'Factoring shows linear dependence.')
          stop
        endif
c
        jdim = 0
        do i = 1,idim
          do kcol = jdim + 1,kbt
            if (abs(aamatr(i,kcol)) .ne. 0.) then
              jdim = jdim + 1
              ksprv(jdim) = kcol
              go to 100
            endif
          enddo
  100     continue
        enddo
c
        if (jdim .ne. idim) then
          j2 = ilnobl(uphase(np))
          write (noutpt,1130) uphase(np)(1:j2),jdim,idim
          write (nttyo,1130) uphase(np)(1:j2),jdim,idim
 1130     format(/' * Error (EQ6/jgibbs) Programming error trap:',
     $    " Couldn't find a",/7x,'set of independent basis species',
     $    ' reaction coefficients needed',/7x,'to calculate the',
     $    ' hypothetical affinity of ',a,',',/7x,'which is involved',
     $    ' in a violation of the mineralogic phase rule.',/7x,
     $    'Previous testing should guarantee the existence of such',
     $    ' a set.',/7x,'Found ',i3,' independent species for ',i3,
     $    ' "other" phases.')
          stop
        endif
c
        if (iodb(1) .gt. 0) then
          write (noutpt,1140)
 1140     format(/3x,'Cognate independent basis species:',/)
          do j = 1,jdim
            kcol = ksprv(j)
            nb = iindx1(kcol)
            nse = nbasp(nb)
            j4 = ilnobl(uspec(nse))
            write (noutpt,1110) uspec(nse)(1:j4)
          enddo
        endif
c
c       Write the matrix for calculating the hypothetical affinity.
c
        j = 0
        do kpp = 1,kpt
          if (kpp .ne. kp) then
            kmm = kmprv(kpp)
            nss = iindx1(kmm)
            j = j + 1
            do i = 1,idim
              kcol = ksprv(i)
              nb = iindx1(kcol)
              nse = nbasp(nb)
c
c             Calling sequence substitutions:
c               nss for ns
c
              cx = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,nss,nstmax)
              aamatr(i,j) = cx
            enddo
          endif
        enddo
c
c       Write the right-hand-side vector for calculating the
c       hypothetical affinity.
c
        do j = 1,jdim
          kcol = ksprv(j)
          nb = iindx1(kcol)
          nse = nbasp(nb)
          cx = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns,nstmax)
          rhsvec(j) = cx
        enddo
c
        if (iodb(1) .gt. 0) then
          write (noutpt,1150)
 1150     format(/3x,'Matrix:',/)
          do i = 1,idim
            write (noutpt,1160) (aamatr(i,j), j = 1,jdim)
 1160       format(5x,5(1pe10.3,2x),e10.3)
          enddo
c
          write (noutpt,1170)
 1170     format(/3x,'RHS vector:',/)
          do j = 1,jdim
            write (noutpt,1180) rhsvec(j)
 1180       format(5x,1pe10.3)
          enddo
        endif
c
c       Calling sequence substitutions:
c         idim for kdim
c
        call msolvr(aamatr,delvec,gmmatr,ier,ipivot,idim,kmax,
     $  noutpt,nttyo,qpr,rhsvec)
c
        if (ier .gt. 0) then
          j2 = ilnobl(uphase(np))
          write (noutpt,1220) uphase(np)(1:j2)
          write (nttyo,1220) uphase(np)(1:j2)
 1220     format(/" * Error (EQ6/jgibbs) Can't resolve a violation of",
     $    ' the mineralogical',/7x,'phase rule involving ',a,' because',
     $    ' of a failure',/7x,'in the calculation of the hypothetical',
     $    ' affinity of this phase',/7x,'in the presence of the',
     $    ' other phases involved in the violation.',/7x,"Couldn't",
     $    ' solve the requisite matrix equation.')
          stop
        endif
c
        if (iodb(1) .gt. 0) then
          write (noutpt,1230)
 1230     format(/3x,'The interphase reaction is:',/)
          cx = 1.0
          j2 = ilnobl(uphase(np))
          write (noutpt,1240) cx,uphase(np)(1:j2)
 1240     format(4x,f7.3,2x,a)
          j = 0
          do kpp = 1,kpt
            if (kpp .ne. kp) then
              kmm = kmprv(kpp)
              npp = ipndx1(kmm)
              j = j + 1
              cx = delvec(j)
              if (cx .lt. 0.) then
                acx = -cx
                j3 = ilnobl(uphase(npp))
                write (noutpt,1250) acx,uphase(npp)(1:j3)
 1250           format(2x,'+ ',f7.3,2x,a)
              endif
            endif
          enddo
          write (noutpt,1260)
 1260     format(10x,'==')
          j = 0
          do kpp = 1,kpt
            if (kpp .ne. kp) then
              kmm = kmprv(kpp)
              npp = ipndx1(kmm)
              j = j + 1
              cx = delvec(j)
              if (cx .gt. 0.) then
                j3 = ilnobl(uphase(npp))
                write (noutpt,1250) cx,uphase(npp)(1:j3)
              endif
            endif
          enddo
        endif
c
        affc = -xlks(ns)
        j = 0
        do kpp = 1,kpt
          if (kpp .ne. kp) then
            j = j + 1
            kmm = kmprv(kpp)
            nss = iindx1(kmm)
            affc = affc + xlks(nss)*delvec(j)
          endif
        enddo
        affp(np) = afcnst*affc
c
        if (iodb(1) .gt. 0) then
          j2 = ilnobl(uphase(np))
          write (noutpt,1270) uphase(np)(1:j2),affp(np)
 1270     format(/3x,'The hypothetical affinity for ',a,' is ',f10.5,/)
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Only phases with negative affinities are candidates to be
c     dropped. It is guaranteed that the last phase added will
c     not be one of these. The best approach to affinity scaling
c     is to find the phase which is limiting with respect to how
c     many moles can be created of the last phase added.
c
      kpadd = 0
      do kp = 1,kpt
        km = kmprv(kp)
        np = ipndx1(km)
        if (np .eq. npadd) then
          kpadd = kp
          go to 110
        endif
      enddo
c
      write (noutpt,1280)
      write (nttyo,1280)
 1280 format(/" * Note (EQ6/jgibbs) Don't know which phase was the",
     $ ' last one added'/7x,'to the ES phase assemblage. This',
     $ ' datum is needed to generate',/7x,'the optimum choice of',
     $ ' phase to delete from the assemblage',/7x,'to resolve the',
     $ ' present mineralogical phase rule violation.',/7x,'Will',
     $ ' generate a phase to delete using the available data.')
c
  110 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (kpadd .eq. 0) then
c
c       Use the regular affinity scaling factors.
c
        do kp = 1,kpt
          km = kmprv(kp)
          ns = iindx1(km)
          cscx(kp) = cscale(ns)
        enddo
      else
c
c       Normalize the inter-phase reaction so that the last phase
c       added has a reaction coefficient of -1.0. Note that the
c       reaction at this point is set up for the last phase involved
c       in the violation.
c
        delvec(kpt) = -1.0
        cx = delvec(kpadd)
        if (cx .lt. 0.) cx = -cx
        do kp = 1,kpt
          delvec(kp) = delvec(kp)/cx
        enddo
c
c       Calculate the maximum possible mass of each phase. Define
c       the scaling factor as equivalent to the maximum number
c       of mole equivalents of the last phase added.
c
        do kp = 1,kpt
          cscx(kp) = 1.
          km = kmprv(kp)
          ns = iindx1(km)
          np = ipndx1(km)
          mxxm = 1.e+38
          nr1 = nstsr(1,ns)
          nr2 = nstsr(2,ns)
          do n = nr1,nr2
            nb = nsts(n)
            cx = csts(n)
            mxx = mtb(nb)
            if (mxx.gt.0. .and. cx.gt.0.) then
              mxxe = mxx/cx
              if (mxxe .lt. mxxm) mxxm = mxxe
            endif
          enddo
          acx = abs(delvec(kp))
          acx2 = acx*acx
c
c         Note: dividing mxxm by acx once gives the maximum number
c         of mole equivalents of the last phase addded. Dividing
c         mxmm by acx a second time is the equivalent of multiplying
c         the hypothetical affinity by the same factor, which normalizes
c         the affinities with respect to the corresponding reaction
c         stoichiometries. Such normalization causes each of the
c         hypothetical affinities to take on the same magnitude.
c
          cscx(kp) = mxxm/acx2
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(1) .gt. 0) then
        if (kpadd .gt. 0) then
c
c         Write an example of the linking inter-phase reaction such
c         that the last phase added is on the right hand side. It is
c         presently on the left hand side, so reverse the reaction
c         coefficients.
c
          do kp = 1,kpt
            delvec(kp) = -delvec(kp)
          enddo
c
          write (noutpt,1340)
 1340     format(//3x,'The phases are linked by the reaction:',/)
c
          j = 0
          do kp = 1,kpt
            km = kmprv(kp)
            np = ipndx1(km)
            cx = delvec(kp)
            if (cx .lt. 0.) then
              acx = -cx
              j3 = ilnobl(uphase(np))
              if (j .eq. 0) then
                write (noutpt,1350) acx,uphase(np)(1:j3)
 1350           format(4x,f7.3,2x,a)
                j = j + 1
              else
                write (noutpt,1360) acx,uphase(np)(1:j3)
 1360           format(2x,'+ ',f7.3,2x,a)
              endif
            endif
          enddo
          write (noutpt,1370)
 1370     format(10x,'==')
          do kp = 1,kpt
            km = kmprv(kp)
            np = ipndx1(km)
            cx = delvec(kp)
            if (cx .gt. 0.) then
              j3 = ilnobl(uphase(np))
              write (noutpt,1360) cx,uphase(np)(1:j3)
            endif
          enddo
          write (noutpt,1380)
 1380     format(1x)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(1) .gt. 0) then
c
c       Write a table of the results.
c
        write (noutpt,1420)
 1420   format(/8x,'Phase',22x,'Affinity',2x,'Scale Factor',2x,
     $  'Scaled Affinity',/)
        do kp = 1,kpt
          km = kmprv(kp)
          np = ipndx1(km)
          afscal = affp(np)/cscx(kp)
          write (noutpt,1430) uphase(np),affp(np),cscx(kp),afscal
 1430     format(6x,a24,3x,f10.5,3x,1pe10.3,3x,e10.3)
        enddo
        if (kpadd .eq. 0) then
          write (noutpt,1440)
 1440     format(/3x,'Scaling employs the regular affinity scaling',
     $    ' factors.',/)
        else
          j2 = ilnobl(uphase(npadd))
          write (noutpt,1450) uphase(npadd)(1:j2)
 1450     format(/3x,'Scaling is based on stoichiometric normalization',
     $    ' of the reactions',/3x,'from which the hypothetical',
     $    ' affinities are calculated and on the',/3x,'maximum number',
     $    ' of mole equivalents of each phase in terms of',/3x,a,
     $    ', the last phase added to the ES phase assemblage.',/)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Determine the candidate to drop from the equilibrium phase
c     assemblage.
c
      npdel = 0
      affmin = 0.
      do kp = 1,kpt
        km = kmprv(kp)
        np = ipndx1(km)
        afscal = affp(np)/cscx(kp)
        if (afscal .lt. affmin) then
          if (affmin.ge.0. .or. npadd.ne.np) then
            affmin = afscal
            npdel = np
          endif
        endif
      enddo
c
      if (npdel .le. 0) then
        j2 = ilnobl(uphase(npdel))
        write (noutpt,1480) uphase(npdel)(1:j2)
        write (nttyo,1480) uphase(npdel)(1:j2)
 1480   format(/' * Error (EQ6/jgibbs) Programming error trap:',
     $  " Couldn't find",/7x,'a candidate phase to delete to',
     $  ' resolve the mineralogical',/7x,'phase rule violation.')
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
