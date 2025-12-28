      subroutine ncmpve(acflg,act,actlg,cdrs,cgexj,eps100,iern1,
     $ iern2,ietmax,jern1,jern2,jetmax,jflag,jgext,jsflag,losp,
     $ mgext,moph,mosp,nbasp,nbt,nbtmax,ndrs,ndrsmx,ndrsr,netmax,
     $ noutpt,nptmax,nstmax,nttyo,ugexj,uphase,uspec,xbar,
     $ xbarlg,xlks)
c
c     This subroutine computes part of the "expansion" of the
c     basis set variable data that yields the "total" system
c     description. This part is the expansion giving the mole
c     fractions and number of moles of components of generic ion
c     exchangers, for all exchange models. This part of the expansion
c     is necessarily iterative in the general case. The overall
c     expansion is conducted by EQLIB/ncmpex.f.
c
c     This subroutine is called by:
c
c       EQLIB/ncmpex.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       acflg  = array of logarithms of activity coefficients
c       actlg  = array of logarithms of species activities
c                  (input consists of valid results for aqueous
c                  basis species; output consists of the same
c                  for generic ion exchanger species)
c       cdrs   = array of reaction coefficients
c       mosp   = array of numbers of moles of species
c                  (input contains valid results for basis species
c                  and carry-over values for non-basis species)
c       moph   = array of numbers of moles of phases
c       ndrs   = array parallel to cdrs giving the index of the
c                  corresponding species
c       ndrsr  = array giving the range in the cdrs/ndrs arrays
c                  corresonding to the reaction associated with a
c                  given species
c       xlks   = array of equilibrium constants
c
c     Principal output:
c
c       act    = array of species activities
c                  (output consists of the subset for generic exchanger
c                  species)
c       actlg  = array of logarithms of species activities
c       mgext  = array of the sum of the number of moles of exchanger
c                  species in given sites of generic ion exchange
c                  phases
c       mosp   = array of numbers of moles of species
c                  (output consists of valid results for non-basis
c                  species; valid results for basis species were
c                  input and are retained)
c       losp   = array of logarithms of numbers of moles of species
c       xbar   = array of mole fractions of species
c       xbarlg = array of logarithms of mole fractions of species
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ietmax,jetmax,nbtmax,ndrsmx,netmax,nptmax,nstmax
c
      integer noutpt,nttyo
c
      integer jern1(jetmax,netmax),jern2(jetmax,netmax),jflag(nstmax),
     $ jgext(netmax),jsflag(nstmax),nbasp(nbtmax),ndrs(ndrsmx),
     $ ndrsr(2,nstmax)
c
      integer iern1,iern2,nbt
c
      character*48 uspec(nstmax)
      character*24 uphase(nptmax)
      character*8 ugexj(jetmax,netmax)
c
      real*8 acflg(nstmax),act(nstmax),actlg(nstmax),cdrs(ndrsmx),
     $ cgexj(jetmax,netmax),losp(nstmax),mgext(jetmax,netmax),
     $ moph(nptmax),mosp(nstmax),xbar(nstmax),xbarlg(nstmax),
     $ xlks(nstmax)
c
      real*8 eps100
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer icycle,icycmx,iterve,itvemx,je,j2,j3,j4,n,nb,ndlxfc,ne,
     $ np,nrr1,nrr2,nrsxfc,nr1,nr2,ns,nse,nsi,nss
c
      integer ilnobl,nbasis
c
      logical qprint
c
      character*8 ux8
c
      real*8 aajx,ardlx,ardlxo,aresx,aresxo,cgxj,delx,delxl,delxu,
     $ dlxfnc,dm,cxs,lax,lxx,mtxj,mtxjc,resx,rhsx,rsxfnc,resxo,xbarmx,
     $ xbarsm,xx
c
      real*8 texp,tlg
c
c-----------------------------------------------------------------------
c
c     Data statements.
c
      data icycmx /10/
      data itvemx /100/
c
c-----------------------------------------------------------------------
c
      qprint = .false.
c
c     Loop over all generic ion exchanger phases.
c
      do np = iern1,iern2
        ne = np - iern1 + 1
        if (moph(np) .gt. 0.) then
c
c         Have found an exchanger phase present in the system.
c         Loop over all sites.
c
          do je = 1,jgext(ne)
            nrr1 = jern1(je,ne)
            nrr2 = jern2(je,ne)
            cgxj = cgexj(je,ne)
c
c           Determine which species is the basis species for the
c           current site of the current exchanger phase.
c
            nse = 0.
            do ns = nrr1,nrr2
              nb = nbasis(nbasp,nbt,nbtmax,ns)
              if (nb.gt.0 .and. jsflag(ns).le.0) then
                nse = ns
              endif
            enddo
c
c           Here mtxj (saved as mgext(je,ne)) is the sum of the
c           numbers of moles of the exchanger species for the current
c           site. This is equivalent to the sum of the numbers of
c           moles of the corresponding exchange ions. Depending on
c           the ion exchange model in question, this sum may not be
c           simply proportional to the number of moles of the exchanger
c           phase. More generally, it may be a complex function of
c           the site composition. In this case, the value of mtxj
c           must be calculated using an iterative process.
c
c           Calculate an initial estimate from mostly inherited
c           values (use number of moles variables for all exchanger
c           species for the current site). Only the value for the
c           species in the basis set is actually current.
c
            mtxj = 0.
            do ns = nrr1,nrr2
              mtxj = mtxj + mosp(ns)
            enddo
c
c           Optional output.
c
            if (qprint) then
              j2 = ilnobl(ugexj(je,ne))
              j3 = ilnobl(uphase(np))
              write (noutpt,1040) ugexj(je,ne)(1:j2),
     $        uphase(np)(1:j3)
 1040         format(/' Pre-Newton-Raphson expansion for site ',a,
     $        /' of phase ',a,':',/)
            endif
c
            icycle = -1
            resx = 0.
c
c           Loop back here for another cycle of pre-Newton-Raphson
c           correction.
c
  110       icycle = icycle + 1
c
c           Make a tentative expansion to estimate the mole fractions,
c           activities, and numbers of moles of all species belonging
c           to the current site.
c
            call ncmpvh(acflg,act,actlg,cdrs,cgxj,jflag,
     $      jsflag,losp,mosp,mtxj,nbasp,nbt,nbtmax,ndrs,ndrsmx,
     $      ndrsr,nrr1,nrr2,nstmax,xbar,xbarlg,xlks)
c
c           Check the sum of the mole fractions.
c
            xbarsm = 0.
            do ns = nrr1,nrr2
              xbarsm = xbarsm + xbar(ns)
            enddo
c
c           Calculate the residual function.
c
            resxo = resx
            resx = xbarsm - 1.0
c
c           Calculate the improvement functions.
c
            aresxo = abs(resxo)
            aresx = abs(resx)
            if (aresxo .gt. 0.) then
              rsxfnc = (aresxo - aresx)/aresxo
            else
              rsxfnc = 0.
            endif
c
c           Optional output.
c
            if (qprint) then
              if (icycle .ge. 0) then
                write (noutpt,1080) icycle,mtxj,resx,rsxfnc
 1080           format(2x,'iter= ',i2,', mtxj= ',1pe10.3,
     $          ', resx= ',1e10.3,', rsxfnc= ',e10.3)
              endif
            endif
c
            if (resx.gt.0.2 .or. resx.lt.-0.2) then
              if (icycle .lt. icycmx) then
c
c               Make a new estimate of mtxj to use in a new
c               expansion calculation.
c
c               Normalize the mole fractions and recompute the
c               activities.
c
                do ns = nrr1,nrr2
                  xx = xbar(ns)/xbarsm
                  xbar(ns) = xx
                  xbarlg(ns) = tlg(xx)
                  actlg(ns) = cgxj*(xbarlg(ns) + acflg(ns))
                enddo
c
c               Find the species with the largest mole fraction.
c
                xbarmx = 0.
                nsi = 0
                do ns = nrr1,nrr2
                  if (xbar(ns) .gt. xbarmx) then
                    xbarmx = xbar(ns)
                    nsi = ns
                  endif
                enddo
c
c               Assume that mole fraction is correct. If the species
c               with the largest mole fraction is the basis species,
c               use its mole fraction to re-estimate mtxj. If it is
c               not the basis species, use the mass action relation
c               to calculate the mole fraction of the basis species,
c               and use that in turn to re-estimate mtxj.
c
                if (nsi .eq. nse) then
                  xx = xbar(nse)
                else
                  nr1 = ndrsr(1,nsi)
                  nr2 = ndrsr(2,nsi)
                  lax = xlks(nsi)
                  do n = nr1,nr2
                    nss = ndrs(n)
                    if (nss .eq. nse) then
                      cxs = cdrs(n)
                    else
                      lax = lax - cdrs(n)*actlg(nss)
                    endif
                  enddo
                  lax = lax/cxs
                  lxx = (lax/cgxj) - acflg(nse)
                  xx = texp(lxx)
                endif
                mtxjc = mosp(nse)/xx
                if (mtxjc .lt. mosp(nse)) then
c
c                 The value of mtxjc is less than the theoretical
c                 minimum of mosp(nse) (the calculated mole fraction
c                 xx is greater than 1.0). Move mtxj halfway to its
c                 theoretical lower limit.
c
                  dm = mosp(nse) - mtxj
                  mtxj = mtxj + 0.5*dm
                else
c
c                 The value of mtxjc is okay. Use it.
c
                  mtxj = mtxjc
                endif
                go to 110
              endif
            endif
c
c           Optional output.
c
            if (qprint) then
              write (noutpt,1090)
 1090         format(1x)
            endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c           Begin Newton-Raphson iteration. This presumes that a
c           starting value of mtxj is already provided. The initial
c           residual is calculated, even if it is already known.
c           This is intended to make the N-R block as independent
c           as possible from the above pre-N-R block.
c
c           Optional output.
c
            if (qprint) then
              j2 = ilnobl(ugexj(je,ne))
              j3 = ilnobl(uphase(np))
              write (noutpt,1140) ugexj(je,ne)(1:j2),
     $        uphase(np)(1:j3)
 1140         format(/' Newton-Raphson expansion for site ',a,
     $        /' of phase ',a,':',/)
            endif
c
            iterve = 0
            resx = 0.
            delx = 0.
            resxo = 0.
            rsxfnc = 0.
            ardlx = 0.
            ardlxo = 0.
            dlxfnc = 0.
            nrsxfc = 0
            ndlxfc = 0
c
c           Loop back here to continue more Newton-Raphson
c           iterations.
c
  120       continue
c
c           Make a tentative expansion to estimate the mole fractions,
c           activities, and numbers of moles of all species belonging
c           to the current site.
c
            call ncmpvh(acflg,act,actlg,cdrs,cgxj,jflag,
     $      jsflag,losp,mosp,mtxj,nbasp,nbt,nbtmax,ndrs,ndrsmx,
     $      ndrsr,nrr1,nrr2,nstmax,xbar,xbarlg,xlks)
c
c           Estimate the mole fractions and activities of the
c           basis species for the current site.
c
c           Check the sum of the mole fractions.
c
            xbarsm = 0.
            do ns = nrr1,nrr2
              xbarsm = xbarsm + xbar(ns)
            enddo
c
c           Calculate the residual function.
c
            resxo = resx
            resx = xbarsm - 1.0
c
c           Calculate the improvement functions.
c
            aresxo = abs(resxo)
            aresx = abs(resx)
            if (aresxo .gt. 0.) then
              rsxfnc = (aresxo - aresx)/aresxo
            else
              rsxfnc = 0.
            endif
c
            ardlxo = ardlx
            ardlx = abs(delx/mtxj)
            if (ardlxo .gt. 0.) then
              dlxfnc = (ardlxo - ardlx)/ardlxo
            else
              dlxfnc = 0.
            endif
c
            if (rsxfnc .ge. 0.) then
              nrsxfc = 0
            else
              nrsxfc = nrsxfc + 1
            endif
c
            if (dlxfnc .ge. 0.) then
              ndlxfc = 0
            else
              ndlxfc = ndlxfc + 1
            endif
c
c           Optional output.
c
            if (qprint) then
              if (iterve .ge. 0) then
                write (noutpt,1150) iterve,mtxj,resx,rsxfnc
 1150           format(2x,'iter= ',i2,', mtxj= ',1pe10.3,
     $          ', resx= ',1e10.3,', rsxfnc= ',e10.3)
              endif
            endif
c
c           Test for convergence.
c
            if (aresx .le. eps100) go to 170
c
c           Test for the maximum number of iterations.
c
            if (iterve .eq. itvemx) then
              j2 = ilnobl(ugexj(je,ne))
              j3 = ilnobl(uphase(np))
              write (ux8,'(i5)') itvemx
              call lejust(ux8)
              j4 = ilnobl(ux8)
              write (noutpt,1170) ugexj(je,ne)(1:j2),
     $        uphase(np)(1:j3),ux8(1:j4)
              write (nttyo,1170) ugexj(je,ne)(1:j2),
     $        uphase(np)(1:j3),ux8(1:j4)
 1170         format(/' * Warning - (EQLIB/ncmpve) Newton-Raphson',
     $        ' iteration failed to'/7x,'converge during an attempt',
     $        ' to expand the system description',/7x,'for site ',a,
     $        ' of phase ',a,'. It reached',/7x,'the limit of ',a,
     $        ' iterations.')
              write (nttyo,1172) aresx
              write (noutpt,1172) aresx
 1172         format(/9x,'abs. residual= ',1pe12.5)
              go to 170
            endif
c
c           Test for non-convergent behavior.
c
            if ((nrsxfc.ge.10 .and. ndlxfc.ge.10) .or.
     $        nrsxfc.ge.20 .or. ndlxfc.ge.20) then
              j2 = ilnobl(ugexj(je,ne))
              j3 = ilnobl(uphase(np))
              write (noutpt,1180) ugexj(je,ne)(1:j2),
     $        uphase(np)(1:j3)
              write (nttyo,1180) ugexj(je,ne)(1:j2),
     $        uphase(np)(1:j3)
 1180         format(/' * Warning - (EQLIB/ncmpve) Newton-Raphson',
     $        ' iteration failed to'/7x,'converge during an attempt',
     $        ' to expand the system description',/7x,'for site ',a,
     $        ' of phase ',a,'. Non-convergent',/7x,'behavior was',
     $        ' detected.')
              go to 170
            endif
c
c           Do another iteration.
c
            iterve = iterve + 1
c
c           Set the right-hand-side.
c
            rhsx = -resx
c
c           Calculate the derivative (1 x 1 Jacobian) of the
c           residual function.
c
            aajx = 0.
            do ns = nrr1,nrr2
              if (ns .eq. nse) then
                aajx = aajx + xbar(ns)
              else
                nr1 = ndrsr(1,ns)
                nr2 = ndrsr(2,ns)
                do n = nr1,nr2
                  nss = ndrs(n)
                  if (nss .eq. nse) then
                    aajx = aajx - cdrs(n)*xbar(ns)/cdrs(nr1)
                    go to 150
                  endif
                enddo
  150           continue
              endif
            enddo
            aajx = -aajx/mtxj
c
            if (aajx .eq. 0.) then
              j2 = ilnobl(ugexj(je,ne))
              j3 = ilnobl(uphase(np))
              write (noutpt,1190) ugexj(je,ne)(1:j2),
     $        uphase(np)(1:j3)
              write (nttyo,1190) ugexj(je,ne)(1:j2),
     $        uphase(np)(1:j3)
 1190         format(/' * Error - (EQLIB/ncmpve) Have encountered',
     $        ' a zero 1 x 1 Jacobian'/7x,'during Newton-Raphson',
     $        ' iteration to expand the description of',/7x,'site ',
     $        a,' of phase ',a,'.')
              stop
            endif
c
c           Solve for the Newton-Raphson correction term.
c
            delx = rhsx/aajx
c
c           Apply under-relaxation limits.
c
            delxu = 9999.*mtxj
            delxl = -0.9999*mtxj
            delx = min(delx,delxu)
            delx = max(delx,delxl)
c
c           Apply a special under-relaxation limit.
c
            mtxjc = mtxj + delx
            if (mtxjc .lt. mosp(nse)) then
c
c             The value of mtxjc is less than the theoretical
c             minimum of mosp(nse).  Set up to move mtxj halfway
c             to its theoretical lower limit.
c
              delx = 0.5*(mosp(nse) - mtxj)
            endif
c
c           Apply the correction term.
c
            mtxj = mtxj + delx
            go to 120
c
c           Save the sum of the number of moles of exchanger species
c           for the current site of the current exchanger phase.
c
  170       mgext(je,ne) = mtxj
c
c           Optional output.
c
            if (qprint) then
              write (noutpt,1200)
 1200         format(1x)
            endif
c
          enddo
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
