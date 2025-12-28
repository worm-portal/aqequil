      subroutine ncmpex(acflg,act,actlg,cdrs,cegexs,cgexj,conc,
     $ conclg,cpgexs,egexjc,egexjf,egexs,eps100,fo2,fo2lg,fsort,
     $ fugac,fugalg,iern1,iern2,ietmax,ifrn1,ifrn2,igas,igstak,
     $ iindx1,ilrn1,ilrn2,imrn1,imrn2,istack,ixrn1,ixrn2,jcsort,
     $ jern1,jern2,jetmax,jflag,jgext,jgsort,jgstak,jjsort,jpflag,
     $ jsflag,jsitex,jssort,jstack,kbt,kdim,kelect,kmax,km1,ko2gaq,
     $ kwater,kxt,loph,losp,lsort,mgext,mrgexs,mtb,moph,mosp,narn1,
     $ narn2,nbasp,nbt,nbtmax,ncmpr,ndrs,ndrsmx,ndrsr,nelect,nern1,
     $ nern2,netmax,ngexsa,ngext,ngrn1,ngrn2,ngt,ngtmax,noutpt,
     $ no2gaq,nphasx,npt,nptmax,nst,nstmax,nttyo,omega,omeglg,
     $ press,qxbarw,q6mode,ugexj,ugexmo,uphase,uspec,xbar,xbarlg,
     $ xbarw,xbarwc,xbrwlc,xbrwlg,xlks,zchar,zgexj,zvclg1,zvec1)
c
c     This subroutine computes all parameters necessary to write the
c     Jacobian matrix from the zvclg1 array. It thus "expands" the
c     basis set variable data.
c
c     This subroutine is called by:
c
c       EQLIB/newton.f
c       EQLIB/nrstep.f
c       EQ3NR/arrset.f
c       EQ6/eqcalc.f
c       EQ6/exivar.f
c       EQ6/optmzr.f
c       EQ6/path.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       acflg  = array of logarithms of activity coefficients
c       cdrs   = array of reaction coefficients
c       ncmpr  = array giving the range in arrays corresponding to
c                  species of those species which belong to a given
c                  phase
c       ndrs   = array parallel to cdrs giving the index of the
c                  corresponding species
c       ndrsr  = array giving the range in the cdrs/ndrs arrays
c                  corresonding to the reaction associated with a
c                  given species
c       press  = pressure, bars
c       qxbarw = flag controlling iterative improvement of xbarw in
c                  the present subroutine:
c                  .false. = no iterative improvement
c                  .true.  = iterative improvement (relevant to EQ6 only)
c       q6mode = flag denoting usage for EQ3NR or EQ6:
c                  .false. = EQ3NR
c                  .true.  = EQ6NR
c       xbarlg = array of logarithms of mole fractions of species;
c                this is primarily an output of this subroutine, but
c                xbarlg(narn1), the mole fraction of water, is used
c                as an input if q6mode is .false.
c       xlks   = array of equilibrium constants
c       zvec1  = array of master variables
c       zvclg1 = array of logarithms of master variables
c
c     Principal output:
c
c       act    = array of species activities
c       actlg  = array of logarithms of species activities
c       conc   = array of species concentrations
c       conclg = array of logarithms of species concentrations
c       fo2    = the oxygen fugacity
c       fo2lg  = logarithm of the oxygen fugacity
c       fugac  = array of gas fugacities
c       fugalg = array of logarithms of gas fugacities
c       moph   = array of numbers of moles of phases
c       loph   = array of logarithms of numbers of moles of phases
c       mosp   = array of numbers of moles of species
c       losp   = array of logarithms of numbers of moles of species
c       xbar   = array of mole fractions of species
c       xbarlg = array of logarithms of mole fractions of species
c                  (but see above under 'principal input")
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ietmax,jetmax,kmax,nbtmax,ndrsmx,netmax,ngtmax,nptmax,
     $ nstmax
c
      integer noutpt,nttyo
c
      integer igstak(ngtmax),iindx1(kmax),istack(nstmax),
     $ jcsort(nstmax),jern1(jetmax,netmax),jern2(jetmax,netmax),
     $ jflag(nstmax),jgext(netmax),jgsort(ngtmax),jgstak(ngtmax),
     $ jjsort(nstmax),jpflag(nptmax),jsflag(nstmax),jsitex(nstmax),
     $ jssort(nstmax),jstack(nstmax),nbasp(nbtmax),ncmpr(2,nptmax),
     $ ndrs(ndrsmx),ndrsr(2,nstmax),ngexsa(ietmax,jetmax,netmax),
     $ ngext(jetmax,netmax),nphasx(nstmax)
c
      integer iern1,iern2,ifrn1,ifrn2,igas,ilrn1,ilrn2,imrn1,imrn2,
     $ ixrn1,ixrn2,kbt,kdim,kelect,km1,ko2gaq,kwater,kxt,narn1,narn2,
     $ nbt,nelect,nern1,nern2,ngrn1,ngrn2,ngt,no2gaq,npt,nst
c
      logical qxbarw,q6mode
c
      character*48 uspec(nstmax)
      character*24 ugexmo(netmax),uphase(nptmax)
      character*8 ugexj(jetmax,netmax)
c
      real*8 acflg(nstmax),act(nstmax),actlg(nstmax),cdrs(ndrsmx),
     $ cegexs(ietmax,jetmax,netmax),cgexj(jetmax,netmax),conc(nstmax),
     $ conclg(nstmax),cpgexs(ietmax,jetmax,netmax),
     $ egexjc(jetmax,netmax),egexjf(jetmax,netmax),
     $ egexs(ietmax,jetmax,netmax),fsort(ngtmax),fugac(ngtmax),
     $ fugalg(ngtmax),loph(nptmax),losp(nstmax),lsort(nstmax),
     $ mgext(jetmax,netmax),mrgexs(ietmax,jetmax,netmax),mtb(nbtmax),
     $ moph(nptmax),mosp(nstmax),xbar(nstmax),xbarlg(nstmax),
     $ xlks(nstmax),zchar(nstmax),zgexj(jetmax,netmax),zvclg1(kmax),
     $ zvec1(kmax)
c
      real*8 eps100,fo2,fo2lg,omega,omeglg,press,xbarw,xbarwc,
     $ xbrwlc,xbrwlg
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ix,iy,je,k,kcol,n,nb,ne,ng,np,nrr1,nrr2,nr1,nr2,ns,nss,
     $ nt,nxbarw
c
      character*48 ux48,uy48
      character*24 ux24,uy24
c
      real*8 ajx,ax,axw,axwfnc,axwmax,axwmxo,cx,cxs,cxw,cxx,dxw,fx,
     $ lcx,lx,mw,mx,sigmmc,sx,wconst,wcnstl,xx,xy,zx
c
      real*8 coefdr,texp,tlg
c
c-----------------------------------------------------------------------
c
c     Note: the following statements don't really do anything except
c     cause the compiler not to complain that igas, uspec, uphase, and
c     press are not used.
c
      ix = igas
      iy = ix
      igas = iy
c
      ux48 = uspec(1)
      uy48 = ux48
      uspec(1) = uy48
c
      ux24 = uphase(1)
      uy24 = ux24
      uphase(1) = uy24
c
      xx = press
      xy = xx
      press = xy
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize the master variables from their logarithms.
c
c     Clear the concentrations, activities, and numbers of moles
c     of the aqueous species.
c
      do ns = narn1,narn2
        conc(ns) = 0.
        act(ns) = 0.
        mosp(ns) = 0.
        xbar(ns) = 0.
      enddo
c
      do ns = narn1,narn2
        conclg(ns) = -99999.
        actlg(ns) = -99999.
        losp(ns) = -99999.
        xbarlg(ns) = -99999.
      enddo
c
c     Clear the concentrations, activities, etc., of the generic ion
c     exchange species. However, do not clear the number of moles
c     variables (mosp). These will be used to initiate the expansion.
c
      do ns = nern1,nern2
        conc(ns) = 0.
        act(ns) = 0.
        xbar(ns) = 0.
      enddo
c
      do ns = nern1,nern2
        conclg(ns) = -99999.
        actlg(ns) = -99999.
        losp(ns) = -99999.
        xbarlg(ns) = -99999.
      enddo
c
c     Clear the concentrations, activities, fugacities, and numbers
c     of moles of the gas species.
c
      do ns = ngrn1,ngrn2
        conc(ns) = 0.
        act(ns) = 0.
        mosp(ns) = 0.
        xbar(ns) = 0.
      enddo
c
      do ng = 1,ngt
        fugac(ng) = 0.
      enddo
c
      do ns = ngrn1,ngrn2
        conclg(ns) = -99999.
        actlg(ns) = -99999.
        losp(ns) = -99999.
        xbarlg(ns) = -99999.
      enddo
c
      do ng = 1,ngt
        fugalg(ng) = -99999.
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do k = 1,kdim
        zx = zvclg1(k)
        zvec1(k) = texp(zx)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (kwater .le. 0) then
        write (noutpt,1000)
        write (nttyo,1000)
 1000   format(/' * Error - (EQLIB/ncmpex) Programming error trap:',
     $  ' Have not yet',/7x,'implemented coding for the case of water',
     $  ' not in the basis set.')
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (q6mode) then
c
c       In EQ6, the mass of solvent water varies.
c
        wconst = omega/zvec1(kwater)
        wcnstl = tlg(wconst)
      else
c
c       In EQ3NR, the mass of solvent water is fixed at 1.0 kg.
c
        wconst = 1.0
        wcnstl = 0.
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the number of moles, concentrations, and activities of
c     the basis species. The following loop assumes that all such
c     species are normal aqueous solute species. Corrections must
c     be made after this loop for the species water, aqueous O2(g),
c     and aqueous e-, if they are in the basis set.
c
      do k = 1,kbt
        nb = iindx1(k)
        ns = nbasp(nb)
        lx = zvclg1(k)
        mx = zvec1(k)
        losp(ns) = lx
        mosp(ns) = mx
c
        if (ns.ge.narn1 .and. ns.le.narn2) then
          conclg(ns) = wcnstl + lx
          conc(ns) = wconst*mx
          ax = conclg(ns) + acflg(ns)
          actlg(ns) = ax
          act(ns) = texp(ax)
        elseif (ns.ge.nern1 .and. ns.le.nern2) then
c
          conclg(ns) = wcnstl + lx
          conc(ns) = wconst*mx
c
c         Mole fractions and activities of generic ion exchange species
c         will be calculated later. For some exchange models (e.g.,
c         Gapon), these quantities could be calculated here. However,
c         for others (e.g., Vanselow), that is not possible, at least
c         in the general case.
c
        else
          write (noutpt,1010)
          write (nttyo,1010)
 1010     format(/' * Error - (EQLIB/ncmpex) Programming error trap:',
     $    ' Have not yet',/7x,'implemented coding for the case of a',
     $    ' basis species which is neither',/7x,'an aqueous species',
     $    ' nor an ion-exchanger species.')
          stop
        endif
c
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The following is a return point for correcting the mole fraction
c     of water in EQ6. Here nxbarw is an iteration counter.
c
      nxbarw = 0
      xbarw = xbarwc
      xbrwlg= xbrwlc
      axwmxo = 0.
  200 continue
c
c     Make corrections for water if it is a basis species. The
c     following coding assumes that water is the first species in the
c     aqueous solution (i.e., it has species index narn1).
c
      if (kwater .gt. 0) then
        if (q6mode) then
          xbarlg(narn1) = xbrwlg
          xbar(narn1) = xbarw
        else
          xbrwlg = zvclg1(kwater)
          xbarlg(narn1) = xbrwlg
          xbarw = texp(xbrwlg)
        endif
c
        ax = xbarlg(narn1) + acflg(narn1)
        actlg(narn1) = ax
        act(narn1) = texp(ax)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make corrections for aqueous O2(g), if it is a basis species.
c
      if (ko2gaq .gt. 0) then
        fo2lg = zvclg1(ko2gaq)
        fo2 = zvec1(ko2gaq)
c
c       Note- act(no2gaq) and actlg(no2gaq) will contain the oxygen
c       fugacity and its logarithm to simplify looping over mass action
c       expressions. Technically, the aqueous O2(g) species has no
c       activity.
c
        actlg(no2gaq) = fo2lg
        act(no2gaq) = fo2
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make corrections for aqueous e-, if it is a basis species.
c
      if (kelect .gt. 0) then
        actlg(nelect) = zvclg1(kelect)
        act(nelect) = zvec1(kelect)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the mole fraction and activity of water, if this
c     species is not a basis species.
c
      if (kwater .le. 0) then
        nr1 = ndrsr(1,narn1)
        nr2 = ndrsr(2,narn1)
        cxs = cdrs(nr1)
        cxx = -xlks(narn1) + cxs*acflg(narn1)
        do n = nr1 + 1,nr2
          nss = ndrs(n)
          cxx = cxx + cdrs(n)*actlg(narn1)
        enddo
        cxx = -cxx/cxs
        xbarlg(narn1) = cxx
        xbar(narn1) = texp(cxx)
        ax = cxx + acflg(narn1)
        actlg(narn1) = ax
        act(narn1) = texp(ax)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute concentrations and activities of the non-basis
c     aqueous solute species.
c
      do ns = narn1 + 1,narn2
        if (jflag(ns) .eq. 30) then
          if (jsflag(ns) .le. 0) then
            nr1 = ndrsr(1,ns)
            nr2 = ndrsr(2,ns)
            cxs = cdrs(nr1)
            cxx = -xlks(ns) + cxs*acflg(ns)
            do n = nr1 + 1,nr2
              nss = ndrs(n)
              cxx = cxx + cdrs(n)*actlg(nss)
            enddo
            cxx = -cxx/cxs
            conclg(ns) = cxx
            conc(ns) = texp(cxx)
            ax = cxx + acflg(ns)
            actlg(ns) = ax
            act(ns) = texp(ax)
          endif
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (kwater.gt.0 .and. q6mode .and. qxbarw) then
c
c       In EQ6 mode, check to see that the mole fraction of water is
c       sufficiently well determined. If not, iteratively improve it.
c       This is done using a 1-variable Newton-Raphson method. This
c       is sometimes necessary because (1) to calculate the
c       concentrations of the dependent aqueous species, one may have
c       to have the mole fraction of water (xbarw) and (2) vice
c       versa. The process is initiated by using an "old" value
c       of the mole fraction of water (that of xbarwc on entering
c       the present subroutine).
c
c       Recalculate the mole fraction of water using the
c       updated concentrations of all dependent aqueous species.
c       The new value is stored as "xbarwc".
c
c       Note: here the jcsort array, which provides for a sorted
c       summation to calculate Sigma m (sigmmc), has not been updated.
c       This is extremely unlikely to cause a problem here. An unsorted
c       calculation would almost certainly be adequate.
c
        call csigm(conc,jcsort,narn1,narn2,nstmax,sigmmc)
        xbarwc = omega/(omega + sigmmc)
c
c       Calculate the alpha residual (axw) and the max norm (axwmax).
c
        axw = xbarwc - xbarw
        axwmax = abs(axw)
c
c       Calculate the improvement function (axwfnc).
c
        axwfnc = 0.
        if (axwmxo .gt. 0.) axwfnc = (axwmxo -axwmax)/axwmxo
c
        axwmxo = axwmax
c
        if (nxbarw .le. 0) write (noutpt,1200)
 1200   format(/1x)
        write (noutpt,1210) nxbarw,xbarw,axwmax,axwfnc
 1210   format(2x,'iter= ',i3,2x,'xbarw= ',f12.10,2x,
     $  'axwmax= ',1pe10.3,2x,'axwfnc= ',e10.3)
c
        if (axwmax .gt. eps100) then
          if (nxbarw .le. 20) then
c
c           Calculate an improved value of xbarw.
c
            nxbarw = nxbarw + 1
            sx = 0.
            do ns = narn1 + 1,narn2
              if (jflag(ns) .eq. 30) then
                if (jsflag(ns) .le. 0) then
                  nr1 = ndrsr(1,ns)
                  nr2 = ndrsr(2,ns)
                  cxs = cdrs(nr1)
c
c                 Calling sequence substitutions:
c                   narn1 for nse
c
                  cxw = coefdr(cdrs,ndrs,ndrsmx,ndrsr,narn1,ns,nstmax)
                  sx = sx + cxw*conc(ns)/cxs
                endif
              endif
            enddo
c
c           Here ajx is the 1 x 1 Jacobian, and dxw is (initially) the
c           raw correction term.
c
            ajx = (xbarwc/omega)*sx - 1.0
            dxw = -axw/ajx
c
c           Apply limits to the correction.
c
            xx = xbarw + dxw
            if (xx .gt. 1.0) dxw = 0.5*(1.0 - xbarw)
            if (xx .le. 0.0) dxw = -0.5*xbarw
            dxw = min(dxw,0.05)
            dxw = max(dxw,-0.05)
c
c           Make the correction and try again.
c
            xbarw = xbarw + dxw
            xbrwlg = tlg(xbarw)
            go to 200
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make corrections for water. The following coding assumes that
c     water is the first species in the aqueous solution (i.e., it
c     has species index narn1).
c
      conc(narn1) = omega
      conclg(narn1) = omeglg
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make corrections for aqueous O2(g).
c
      if (no2gaq .gt. 0) then
        conc(no2gaq) = 0.
        conclg(no2gaq) = -99999.
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make corrections for aqueous e-.
c
      if (nelect .gt. 0) then
        conc(nelect) = 0.
        conclg(nelect) = -99999.
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the numbers of moles of the non-basis aqueous solute
c     species.
c
      do ns = narn1 + 1,narn2
        if (jflag(ns).eq.30 .and. jsflag(ns).le.0) then
          losp(ns) = conclg(ns) - wcnstl
          mosp(ns) = conc(ns)/wconst
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make corrections for water. The following coding assumes that
c     water is the first species in the aqueous solution (i.e., it
c     has species index narn1).
c
      if (.not.q6mode) then
        mosp(narn1) = omega
        losp(narn1) = omeglg
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make corrections for aqueous O2(g).
c
      if (no2gaq .gt. 0) then
        mosp(no2gaq) = 0.
        losp(no2gaq) = -99999.
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make corrections for aqueous e-.
c
      if (nelect .gt. 0) then
        mosp(nelect) = 0.
        losp(nelect) = -99999.
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the mole fractions and activities of the basis and
c     non-basis ion-exchanger species, and the numbers of moles
c     of the non-basis species.
c
      call ncmpve(acflg,act,actlg,cdrs,cgexj,eps100,iern1,
     $ iern2,ietmax,jern1,jern2,jetmax,jflag,jgext,jsflag,losp,
     $ mgext,moph,mosp,nbasp,nbt,nbtmax,ndrs,ndrsmx,ndrsr,netmax,
     $ noutpt,nptmax,nstmax,nttyo,ugexj,uphase,uspec,xbar,
     $ xbarlg,xlks)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the concentrations (mol/kg.H2O) and numbers of moles of
c     the non-basis ion-exchanger species.
c
      do np = iern1,iern2
        if (moph(np) .gt. 0.) then
          ne = np - iern1 + 1
          do je = 1,jgext(ne)
            nrr1 = jern1(je,ne)
            nrr2 = jern2(je,ne)
            do ns = nrr1,nrr2
              if (jflag(ns).eq.30 .and. jsflag(ns).le.0) then
                lcx = wcnstl + losp(ns)
                cx = wconst*mosp(ns)
                conclg(ns) = lcx
                conc(ns) = cx
              endif
            enddo
          enddo
c
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the equivalent fractions (egexs) and mole ratios (mrgexs)
c     of exchanger species of generic ion exchanger phases.
c
      call gegexs(cegexs,cgexj,egexjc,egexjf,egexs,iern1,
     $ iern2,ietmax,jern1,jetmax,jgext,moph,mosp,mrgexs,netmax,
     $ ngexsa,ngext,noutpt,nptmax,nstmax,nttyo,zchar,zgexj)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (q6mode) then
c
c       Compute the numbers of moles and concentrations of the
c       species belonging to non-aqueous phases present in the
c       equilibrium sytem (ES).
c
        if (kxt .ge. km1) then
          do kcol = km1,kxt
            ns = iindx1(kcol)
            losp(ns) = zvclg1(kcol)
            mosp(ns) = zvec1(kcol)
            conclg(ns) = wcnstl + losp(ns)
            conc(ns) = wconst*mosp(ns)
          enddo
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Sort species according to log masses.
c
      call sortsp(iern1,iern2,istack,jcsort,jern1,jern2,
     $ jgext,jsitex,jetmax,jjsort,jssort,jstack,losp,lsort,ncmpr,
     $ nern1,nern2,netmax,noutpt,nphasx,npt,nptmax,nst,nstmax,nttyo)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the number of moles of water. The following coding
c     assumes that the aqueous solution is the first phase (i.e.,
c     has phase index 1).
c
      moph(1) = 0.
      loph(1) = -99999.
      mw = 0
      do nss = narn1,narn2
        ns = jcsort(nss)
        mw = mw + mosp(ns)
      enddo
      moph(1) = mw
      loph(1) = tlg(mw)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (q6mode) then
c
c       Compute the mole fractions of the aqueous species. The following
c       coding assumes that water is the first species in the aqueous
c       solution (i.e., has species index narn1).
c
        if (mw .gt. 0.) then
          do ns = narn1,narn2
            xx = mosp(ns)/mw
            xbar(ns) = xx
            xbarlg(ns) = tlg(xx)
          enddo
          xbarw = xbar(narn1)
          xbrwlg = xbarlg(narn1)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (q6mode) then
c
c       Compute the number of moles of non-aqueous phases. The following
c       coding assumes that the aqueous solution is the first phase
c       (i.e., has phase index 1).
c
c       Skip the exchanger phases. Start with the pure liquid phases.
c
        do np = ilrn1,ilrn2
          moph(np) = 0.
          loph(np) = -99999.
          if (jpflag(np) .le. 0) then
            ns = ncmpr(1,np)
            mx = mosp(ns)
            moph(np) = mx
            loph(np) = tlg(mx)
          endif
        enddo
c
c       Pure mineral phases.
c
        do np = imrn1,imrn2
          moph(np) = 0.
          loph(np) = -99999.
          if (jpflag(np) .le. 0) then
            ns = ncmpr(1,np)
            mx = mosp(ns)
            moph(np) = mx
            loph(np) = tlg(mx)
          endif
        enddo
c
c       Fixed fugacity phases.
c
        do np = ifrn1,ifrn2
          moph(np) = 0.
          loph(np) = -99999.
          if (jpflag(np) .le. 0) then
            ns = ncmpr(1,np)
            mx = mosp(ns)
            moph(np) = mx
            loph(np) = tlg(mx)
          endif
        enddo
c
c       Solid solutions.
c
        do np = ixrn1,ixrn2
          moph(np) = 0.
          loph(np) = -99999.
          if (jpflag(np) .le. 0) then
            mx = 0
            nr1 = ncmpr(1,np)
            nr2 = ncmpr(2,np)
            do nss = nr1,nr2
              ns = jcsort(nss)
              mx = mx + mosp(ns)
            enddo
            moph(np) = mx
            loph(np) = tlg(mx)
          endif
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (q6mode) then
c
c       Clear the activities of species belonging to non-aqueous
c       solution phases. The following coding assumes that the aqueous
c       solution is the first phase (i.e., has the phase index 1).
c
        do np = ixrn1,ixrn2
          nr1 = ncmpr(1,np)
          nr2 = ncmpr(2,np)
          do ns = nr1,nr2
            act(ns) = 0.
            actlg(ns) = -99999.
          enddo
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (q6mode) then
c
c       Compute the mole fractions and activities of the species which
c       belong to non-aqueous solution phases in the equilibrium system.
c
        do np = ixrn1,ixrn2
          mx = moph(np)
          if (mx .gt. 0.) then
            nr1 = ncmpr(1,np)
            nr2 = ncmpr(2,np)
            do ns = nr1,nr2
              xx = mosp(ns)/mx
              xbar(ns) = xx
              xbarlg(ns) = tlg(xx)
              ax = xbarlg(ns) + acflg(ns)
              actlg(ns) = ax
              act(ns) = texp(ax)
            enddo
          endif
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the fugacities of the gas species.
c
      do ns = ngrn1,ngrn2
        ng = ns - ngrn1 + 1
        if (jsflag(ns) .lt. 2) then
          nr1 = ndrsr(1,ns)
          nr2 = ndrsr(2,ns)
          nt = nr2 - nr1 + 1
          if (nt .ge. 2) then
            cxs = cdrs(nr1)
            fx = -xlks(ns)
            do n = nr1 + 1,nr2
              nss = ndrs(n)
              cxx = cdrs(n)
              fx = fx + cxx*actlg(nss)
            enddo
            fx = -fx/cxs
            fugalg(ng) = fx
            fugac(ng) = texp(fx)
          else
            fugalg(ng) = actlg(ns)
            fugac(ng) = act(ns)
          endif
        endif
      enddo
c
c     Sort all gas species according to equilibrium fugacities.
c     Put their indices in sorted order in the jgsort array.
c
c     Calling sequence substitutions:
c       fsort for asort
c       fugac for aval
c       jgsort for jsort
c       igstak for istack
c       jgstak for jstack
c       ngtmax for nmax
c       ngt for nval
c
c     Caution: the jgsort array from the last call is recycled as a
c     good starting point. Set jgsort(1) to 0 to make a sort starting
c     from scratch.
c
      call qsortw(fsort,fugac,igstak,jgsort,jgstak,ngtmax,noutpt,
     $ nttyo,ngt)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
