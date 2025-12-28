      subroutine betas(acflg,actlg,afcnst,alpha,amtb,bbig,beta,
     $ betamx,bneg,cdrs,conc,conclg,coval,csts,eh,ehfac,fo2lg,
     $ ibetmx,iebal,iindx1,irdxc3,jcsort,jflag,jsflag,jssort,kbt,
     $ kdim,kelect,khydr,kmax,km1,ko2gaq,kwater,kxt,mtb,mosp,
     $ narn1,narn2,nbasp,nbtmax,ncosp,ndrs,ndrsmx,ndrsr,nelect,
     $ nern1,nern2,nhydr,noutpt,no2gaq,nredox,nst,nstmax,nsts,
     $ nstsmx,nstsr,ntfx,ntfxmx,ntfxt,nttyo,omega,qredox,q6mode,
     $ tfx,ubbig,ubneg,ubetmx,uspec,uzvec1,weight,xbrwlg,xlke,
     $ xlks,zchar)
c
c     This subroutine computes the Newton-Raphson residual functions.
c
c     This subroutine is called by:
c
c       EQLIB/newton.f
c       EQLIB/nrstep.f
c       EQ3NR/arrset.f
c       EQ6/eqcalc.f
c       EQ6/optmzr.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       acflg  = array of log activity coefficients
c       q6mode = flag denoting usage for EQ3NR or EQ6:
c                  .false. = EQ3NR
c                  .true.  = EQ6NR
c
c     Principal output:
c
c       alpha  = array of raw Newton-Raphson residual functions
c       beta   = array of normalized Newton-Raphson residual functions
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer kmax,nbtmax,ndrsmx,nstmax,nstsmx,ntfxmx
c
      integer noutpt,nttyo
c
      integer iindx1(kmax),jcsort(nstmax),jflag(nstmax),
     $ jsflag(nstmax),jssort(nstmax),nbasp(nbtmax),ncosp(nbtmax),
     $ ndrs(ndrsmx),ndrsr(2,nstmax),nsts(nstsmx),nstsr(2,nstmax),
     $ ntfx(ntfxmx)
c
      integer ibetmx,iebal,irdxc3,kbt,kdim,kelect,khydr,km1,ko2gaq,
     $ kwater,kxt,narn1,narn2,nelect,nern1,nern2,nhydr,no2gaq,nredox,
     $ nst,ntfxt
c
      logical qredox,q6mode
c
      character*48 uspec(nstmax),uzvec1(kmax)
      character*48 ubbig,ubneg,ubetmx
c
      real*8 acflg(nstmax),actlg(nstmax),alpha(kmax),amtb(nbtmax),
     $ beta(kmax),cdrs(ndrsmx),conc(nstmax),conclg(nstmax),
     $ coval(nbtmax),csts(nstsmx),mtb(nbtmax),mosp(nstmax),
     $ tfx(ntfxmx),weight(nstmax),xlks(nstmax),zchar(nstmax)
c
      real*8 afcnst,bbig,betamx,bneg,eh,ehfac,fo2lg,omega,
     $ xbrwlg,xlke
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jfl,jlen,kcol,krow,n,nb,nr1,nr2,ns,nsc,nse,nss,ns1
c
      character*56 uspn56
      character*8 unone,uptgas
c
      real*8 abeta,af,alkc,atot,atx,ax,aznse,azns1,azsum,belect,bo2gaq,
     $ bx,bwater,cde,ctot,cx,c1,dx,fx,si,sigmmc,sigza,sigzc,sigzi,
     $ sigzm,sx,tx,xbarwc,xbrwlc,zp
c
      real*8 coefst,tlg
c
c-----------------------------------------------------------------------
c
      data unone  /'None    '/
      data uptgas /'Gas     '/
c
c-----------------------------------------------------------------------
c
      do krow = 1,kdim
        alpha(krow) = 0.
        beta(krow) = 0.
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (q6mode) go to 210
c
c     Compute residuals for EQ3NR.
c
      do krow = 1,kbt
        nb = iindx1(krow)
        nse = nbasp(nb)
        jfl = jflag(nse)
        if (nse.ge.narn1 .and. nse.le.narn2) then
c
c         Aqueous species.
c
          if (krow.eq.kwater .and. jfl.eq.0) then
c
c           The residual for water is based on the equation which
c           defines the mole fraction of water.
c
            call csigm(conc,jcsort,narn1,narn2,nstmax,sigmmc)
            xbarwc = omega/(omega + sigmmc)
            xbrwlc = tlg(xbarwc)
            ax = xbrwlc - xbrwlg
            alpha(krow) = ax
            beta(krow) = ax
c
          elseif (nb .eq. iebal) then
c
c           Charge balance.
c
            call gszm(conc,jcsort,narn1,narn2,nstmax,sigza,sigzc,
     $      sigzi,sigzm,zchar)
            alpha(krow) = sigzi
            if (sigzm .gt. 0.) beta(krow) = sigzi/sigzm
c
          elseif (jfl .eq. 17) then
c
c           Log activity combination.
c
            ns1 = ncosp(nb)
            aznse = abs(zchar(nse))
            azns1 = abs(zchar(ns1))
            ax = coval(nb)/azns1
            zp = zchar(nse)*zchar(ns1)
            if (zp .lt. 0.) then
              ax = ax - ( aznse/azns1 )*actlg(ns1)
              else
              ax = ax + ( aznse/azns1 )*actlg(ns1)
            endif
            dx = ax - actlg(nse)
            alpha(krow) = dx
            beta(krow) = dx
c
          elseif (jfl .eq. 18) then
c
c           Mean log activity.
c
            ns1 = ncosp(nb)
            aznse = abs(zchar(nse))
            azns1 = abs(zchar(ns1))
            azsum = aznse + azns1
            ax = ( azsum*coval(nb) - aznse*actlg(ns1) )/azns1
            dx = ax/azns1 - actlg(nse)
            alpha(krow) = dx
            beta(krow) = dx
c
          elseif (jfl .eq. 21) then
c
c           pHCl.
c
            ns1 = ncosp(nb)
            aznse = abs(zchar(nse))
            azns1 = abs(zchar(ns1))
            ax = -coval(nb)/azns1
            zp = zchar(nse)*zchar(ns1)
            if (zp .lt. 0.) then
              ax = ax - ( aznse/azns1 )*actlg(ns1)
              else
              ax = ax + ( aznse/azns1 )*actlg(ns1)
            endif
            dx = ax - actlg(nse)
            alpha(krow) = dx
            beta(krow) = dx
c
          elseif (jfl.eq.25 .or. jfl.eq.27) then
c
c           Heterogeneous or homogeneous equilibrium with a specified
c           species (computing fO2 from an aqueous redox couple is not
c           done here).
c
            if (jfl .eq. 25) then
              nsc = ncosp(nb)
            else
              nsc = nse
            endif
            cx = xlks(nsc)
            nr1 = ndrsr(1,nsc)
            nr2 = ndrsr(2,nsc)
c
            do n = nr1,nr2
              nss = ndrs(n)
              if (nss .eq. nse) then
                cde = cdrs(n)
                if (nse .ne. no2gaq) cx = cx - cde*acflg(nse)
              elseif (nss .eq. nsc) then
                if (uspec(nsc)(25:32) .eq. uptgas(1:8)) then
                  cx = cx - cdrs(n)*coval(nb)
                endif
              elseif (nss .eq. no2gaq) then
                cx = cx - cdrs(n)*fo2lg
              else
                cx = cx - cdrs(n)*actlg(nss)
              endif
            enddo
c
            cx = cx/cde
            c1 = conclg(nse)
            if (nse .eq. no2gaq) c1 = fo2lg
            dx = cx - c1
            alpha(krow) = dx
            beta(krow) = dx
c
          elseif (nse .eq. no2gaq) then
c
c           log fO2.
c
            if (irdxc3 .eq. -1) then
c
c             Eh residual (Note- if a pe- value was input, it has been
c             converted to an Eh value by EQLIB/setup.f). Here it is
c             assumed that water is the first species in the aqueous
c             phase (has the index narn1).
c
              fx = 4.*eh/ehfac  + xlke + 2.*actlg(narn1)
     $        - 4.*actlg(nhydr)
              dx = fx - fo2lg
              alpha(krow) = dx
              beta(krow) = dx
c
            elseif (irdxc3 .eq. 1) then
c
c             Cross-linking (homogeneous aqueous redox) equilibrium.
c
              cx = xlks(nredox)
              nr1 = ndrsr(1,nredox)
              nr2 = ndrsr(2,nredox)
              do n = nr1,nr2
                nss = ndrs(n)
                if (nss .eq. no2gaq) then
                  cde = cdrs(n)
                else
                  cx = cx - cdrs(n)*actlg(nss)
                endif
              enddo
              cx = cx/cde
              dx = cx - fo2lg
              alpha(krow) = dx
              beta(krow) = dx
c
            else
c
c             Log fO2 is directly specified (irdxc3 = 0).
c
              alpha(krow) = 0.
              beta(krow) = 0.
            endif
          elseif (jfl.ge.7 .and. jfl.le.11) then
c
c           Alkalinity balance.
c
            call calk(alkc,conc,nstmax,ntfx,ntfxmx,ntfxt,tfx)
            atot = coval(nb)
            dx = alkc - atot
            alpha(krow) = dx
            beta(krow) = dx/atot
c
          elseif (jfl .eq. 16) then
c
c           Log activity.
c
            dx = conclg(nse) + acflg(nse) - coval(nb)
            alpha(krow) = dx
            beta(krow) = dx
c
          elseif (jfl.eq.19 .or. jfl.eq.20) then
c
c           pX (including pH).
c
            dx = conclg(nse) + acflg(nse) + coval(nb)
            alpha(krow) = dx
            beta(krow) = dx
c
          elseif (jfl.eq.22 .or. jfl.eq.23) then
c
c           pmX (including pmH).
c
            dx = conclg(nse) + coval(nb)
            alpha(krow) = dx
            beta(krow) = dx
c
          elseif (jfl.ge.0 .and. jfl.le.3) then
c
c           Mass balance.
c
            sx = 0.
            do nss = narn1,narn2
              ns = jcsort(nss)
              weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
              sx = sx + weight(ns)*conc(ns)
            enddo
c           do nss = nern1,nern2
c             ns = jcsort(nss)
c             weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
c             sx = sx + weight(ns)*conc(ns)
c           enddo
            ctot = coval(nb)
            dx = sx - ctot
            alpha(krow) = dx
            beta(krow) = dx/ctot
c
          else
c
c           Have found a bad jflag value.
c
c           Calling sequence substitutions:
c             uzvec1(krow) for unam48
c
            call fmspnm(jlen,uzvec1(krow),uspn56)
            write (noutpt,1000) jfl,uspn56(1:jlen)
            write (nttyo,1000) jfl,uspn56(1:jlen)
 1000       format(/' * Error - (EQLIB/betas) Programming error trap:',
     $      /7x,'Have encountered a bad jflag value of ',i3,' for',
     $      /7x,'the species ',a,'.')
            stop
          endif
c
        elseif (nse.ge.nern1 .and. nse.le.nern2) then
c
c         Generic ion exchange species.
c
          if (jfl .eq. 0) then
c
c           Mass balance.
c
            sx = 0.
c           do nss = narn1,narn2
c             ns = jcsort(nss)
c             weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
c             sx = sx + weight(ns)*conc(ns)
c           enddo
            do nss = nern1,nern2
              ns = jcsort(nss)
              weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
              sx = sx + weight(ns)*conc(ns)
            enddo
            ctot = coval(nb)
            dx = sx - ctot
            alpha(krow) = dx
            beta(krow) = dx/ctot
c
          else
c
c           Have found a bad jflag value.
c
c           Calling sequence substitutions:
c             uzvec1(krow) for unam48
c
            call fmspnm(jlen,uzvec1(krow),uspn56)
            write (noutpt,1000) jfl,uspn56(1:jlen)
            write (nttyo,1000) jfl,uspn56(1:jlen)
          endif
        endif
      enddo
      go to 300
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  210 continue
c
c     Compute residuals for EQ6.
c
c     Compute mass balance elements.
c
      do krow = 1,kbt
        nb = iindx1(krow)
        nse = nbasp(nb)
c
        if ((nse.ne.no2gaq .and. nse.ne.nelect) .or. qredox) then
c
c         Note: this calculation is not made if the current basis
c         species is O2(g) or e-, and the current problem has no redox
c         aspect. In the absence of a redox aspect, the log fO2 or the
c         log a(e-) is treated as a "dead" iteration variable.
c
          sx = 0.
          do nss = 1,nst
            ns = jssort(nss)
            weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
            if (weight(ns) .ne. 0.) then
              if (mosp(ns) .gt. 0.) then
                tx = weight(ns)*mosp(ns)
                sx = sx + tx
              endif
            endif
          enddo
          sx = sx - mtb(nb)
          alpha(krow) = sx
          amtb(nb) = mtb(nb)
c
          if (krow.eq.kwater .or. krow.eq.khydr .or.
     $      krow.eq.ko2gaq .or. krow.eq.kelect) then
c
c           For the species H2O, H+, OH-, O2(g,aq), and e-, compute
c           "absolute" mass balance totals defined as:
c
c             sum of positive terms + abs(sum of negative terms)
c
c           Use this to calculate the corresponding beta residuals.
c
            sx = 0.
            do nss = 1,nst
              ns = jssort(nss)
              weight(ns) = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
              tx = weight(ns)*mosp(ns)
              atx = abs(tx)
              sx = sx + atx
            enddo
            amtb(nb) = sx
          endif
c
          beta(krow) = alpha(krow)/amtb(nb)
        endif
      enddo
c
c     Compute mass action elements.
c
      do krow = km1,kxt
        ns = iindx1(krow)
        call afcalc(actlg,af,afcnst,cdrs,jflag,jsflag,ndrs,ndrsmx,
     $  ndrsr,ns,nstmax,si,xlks)
        alpha(krow) = si
        beta(krow) = si
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get characteristic residual parameters.
c
c       betamx = max norm of beta
c       ubetmx = name of the corresponding species
c       bbig   = value of largest (positive) mass balance residual
c       ubbig  = name of the corresponding species
c       bneg   = value of largest (negative) mass balance residual
c       ubneg  = name of the corresponding species
c
  300 betamx = 0.
      ibetmx = 0
      ubetmx = unone
      do kcol = 1,kdim
        abeta = abs(beta(kcol))
        if (abeta .gt. betamx) then
          betamx = abeta
          ibetmx = kcol
          ubetmx = uzvec1(kcol)
        endif
      enddo
c
      bbig = 0.
      bneg = 0.
      ubbig = unone
      ubneg = unone
c
      if (.not.q6mode) then
        if (kwater .gt. 0) then
          bwater = beta(kwater)
          beta(kwater) = 0.
        endif
        if (kelect .gt. 0) then
          belect = beta(kelect)
          beta(kelect) = 0.
        endif
        if (ko2gaq .gt. 0) then
          bo2gaq = beta(ko2gaq)
          beta(ko2gaq) = 0.
        endif
      endif
c
      do kcol = 1,kbt
        nb = iindx1(kcol)
        nse = nbasp(nb)
        if (jflag(nse) .le. 15) then
          bx = beta(kcol)
          if (bx .lt. bneg) then
            bneg = bx
            ubneg = uzvec1(kcol)
          elseif (bx .gt. bbig) then
            bbig = bx
            ubbig = uzvec1(kcol)
          endif
        endif
      enddo
c
      if (.not.q6mode) then
        if (kwater .gt. 0) then
          beta(kwater) = bwater
        endif
        if (kelect .gt. 0) then
          beta(kelect) = belect
        endif
        if (ko2gaq .gt. 0) then
          beta(ko2gaq) = bo2gaq
        endif
      endif
c
      end
