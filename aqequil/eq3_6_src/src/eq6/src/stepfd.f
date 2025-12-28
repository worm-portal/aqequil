      subroutine stepfd(acflg,acflg0,affp0,affp,afrc0,afrc1,aw0,aw1,
     $ delxi,dxsm00,eh0,eh1,emop,emop0,emos,emos0,fdafm1,fdaf0,
     $ fdarm1,fdar0,fdawm1,fdaw0,fdehm1,fdeh0,fdlim,fdo2m1,fdo20,
     $ fdpem1,fdpe0,fdphm1,fdph0,fdrem1,fdre0,fdrim1,fdri0,fdrrm1,
     $ fdrr0,fdsem1,fdse0,fdzvm1,fdzv0,fje,fje0,fo2lg0,fo2lg1,fxi,
     $ fxi0,iemop,iemop0,iemos,iemos0,iindx0,iindx1,iodb,iopt,ipndx0,
     $ ipndx1,jcode,jreac,jreac0,jpflag,kdim,kdim0,kmax,km1,km10,kmt,
     $ kmt0,kord,kx1,kx10,kxt,kxt0,modr,modr0,moph,moph0,morr,morr0,
     $ mosp,mosp0,mtb,mtb0,nbt,nbtmax,ncmpe,ncmpe0,nodbmx,noptmx,
     $ nordmx,noutpt,npet,npetmx,npet0,npt,nptmax,npts,nrct,nrctmx,
     $ nrd1mx,nrndex,nset,nsetmx,nset0,nstmax,ph0,ph1,qredox,rirec0,
     $ rirec1,rreac0,rreac1,rrelr0,rrelr1,sfcar,sfcar0,sigmam,sigmm0,
     $ tempc,tempc0,time0,time1,uzvec0,uzvec1,xi0,xi1,xim1,xirct,
     $ xirct0,zvclg0,zvclg1,zvec0,zvec1)
c
c     This subroutine saves information at the current point and
c     computes finite differences at the current point of reaction
c     progress.
c
c     This subroutine is called by:
c
c       EQ6/path.f
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
      integer kmax,nbtmax,nodbmx,noptmx,nordmx,npetmx,nptmax,nrctmx,
     $ nrd1mx,nsetmx,nstmax
c
      integer noutpt
c
      integer iemop(npetmx),iemop0(npetmx),iemos(nsetmx),iemos0(nsetmx),
     $ iindx0(kmax),iindx1(kmax),iodb(nodbmx),iopt(noptmx),ipndx0(kmax),
     $ ipndx1(kmax),jcode(nrctmx),jreac(nrctmx),jreac0(nrctmx),
     $ jpflag(nptmax),ncmpe(2,npetmx),ncmpe0(2,npetmx),nrndex(nrctmx)
c
      integer kdim,kdim0,km1,km10,kmt,kmt0,kord,kx1,kx10,kxt,kxt0,
     $ nbt,npet,npet0,npt,npts,nrct,nset,nset0
c
      logical qredox
c
      character*48 uzvec0(kmax),uzvec1(kmax)
c
      real*8 acflg(nstmax),acflg0(nstmax),affp(nptmax),affp0(nptmax),
     $ afrc0(nrctmx),afrc1(nrctmx),dxsm00(nrd1mx),emop(npetmx),
     $ emop0(npetmx),emos(nsetmx),emos0(nsetmx),fdafm1(nordmx,nptmax),
     $ fdaf0(nordmx,nptmax),fdarm1(nordmx,nrctmx),fdar0(nordmx,nrctmx),
     $ fdawm1(nordmx),fdaw0(nordmx),fdehm1(nordmx),fdeh0(nordmx),
     $ fdo2m1(nordmx),fdo20(nordmx),fdpem1(nordmx,npetmx),
     $ fdpe0(nordmx,npetmx),fdphm1(nordmx),fdph0(nordmx),
     $ fdrem1(nordmx,nrctmx),fdre0(nordmx,nrctmx),fdrim1(nrd1mx),
     $ fdri0(nrd1mx),fdrrm1(nrd1mx,nrctmx),fdrr0(nrd1mx,nrctmx),
     $ fdsem1(nordmx,nsetmx),fdse0(nordmx,nsetmx),fdzvm1(nrd1mx,kmax),
     $ fdzv0(nrd1mx,kmax)
c
      real*8 modr(nrctmx),modr0(nrctmx),moph(nptmax),moph0(nptmax),
     $ morr(nrctmx),morr0(nrctmx),mosp(nstmax),mosp0(nstmax),
     $ mtb(nbtmax),mtb0(nbtmax),rreac0(nrctmx),rreac1(nrctmx),
     $ rrelr0(nrctmx),rrelr1(nrctmx),sfcar(nrctmx),sfcar0(nrctmx),
     $ xirct(nrctmx),xirct0(nrctmx),zvclg0(kmax),zvclg1(kmax),
     $ zvec0(kmax),zvec1(kmax)
c
      real*8 aw0,aw1,delxi,eh0,eh1,fdlim,fje,fje0,fo2lg0,fo2lg1,fxi,
     $ fxi0,ph0,ph1,rirec0,rirec1,sigmam,sigmm0,tempc,tempc0,time0,
     $ time1,xi0,xi1,xim1
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer itrunc,j,k,kcol,n,nmax,np,npe,nrc,nse
c
      real*8 dfx,dfxl,dfxu
c
c-----------------------------------------------------------------------
c
      if (npts .eq. 1) then
c
c       Zero all finite differences.
c
        nmax = nrd1mx*kmax
        call initaz(fdzvm1,nmax)
        call initaz(fdzv0,nmax)
c
        nmax = nrd1mx*nrctmx
        call initaz(fdrrm1,nmax)
        call initaz(fdrr0,nmax)
c
        if (iopt(2) .gt. 0) then
          call initaz(fdrim1,nrd1mx)
          call initaz(fdri0,nrd1mx)
        endif
c
        call initaz(fdawm1,nordmx)
        call initaz(fdaw0,nordmx)
        call initaz(fdehm1,nordmx)
        call initaz(fdeh0,nordmx)
        call initaz(fdo2m1,nordmx)
        call initaz(fdo20,nordmx)
        call initaz(fdphm1,nordmx)
        call initaz(fdph0,nordmx)
c
        if (iopt(2) .gt. 0) then
          nmax = nordmx*nrctmx
          call initaz(fdrem1,nmax)
          call initaz(fdre0,nmax)
        endif
c
        nmax = nordmx*nptmax
        call initaz(fdafm1,nmax)
        call initaz(fdaf0,nmax)
c
        nmax = nordmx*npetmx
        call initaz(fdpem1,nmax)
        call initaz(fdpe0,nmax)
c
        nmax = nordmx*nsetmx
        call initaz(fdsem1,nmax)
        call initaz(fdse0,nmax)
c
        nmax = nordmx*nrctmx
        call initaz(fdarm1,nmax)
        call initaz(fdar0,nmax)
c
        go to 200
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the dxsm00 vector.
c
      do n = 1,kord
        j = kord - n + 1
        dxsm00(j + 1) = dxsm00(j) + delxi
      enddo
      dxsm00(1) = delxi
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Overflow protection is automatically activated in the code blocks
c     below in order to run on VAX machines with a small exponent range
c     (+/- 38) for real*8.
c
      itrunc = kord
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute finite differences for the elements of the z vector.
c     These are mostly the number of moles of basis species.
c
      do kcol = 1,kdim
        fdzv0(1,kcol) = (zvec1(kcol) - zvec0(kcol))/delxi
        do j = 2,kord + 1
          dfxu = fdlim*dxsm00(j)
          dfxl = -dfxu
          k = j - 1
          dfx = fdzv0(k,kcol) - fdzvm1(k,kcol)
          dfx = min(dfxu,dfx)
          dfx = max(dfxl,dfx)
          if (abs(dfx) .ge. dfxu) itrunc = min(itrunc,k)
          fdzv0(j,kcol) = dfx/dxsm00(j)
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute finite differences for relative rates.
c
      do nrc = 1,nrct
        if (jreac(nrc).eq.1 .or. jreac(nrc).eq.2) then
c
c         The relative rate must be a constant zero, because no
c         reactant mass remains:
c
c           jreac =  1: exhausted
c           jreac =  2: saturated; any remaining reactant mass is
c                         converted to the corresponding product
c                         phase, so the "reactant" is effectively
c                         exhausted.
c
c         Hence the corresponding finite differences must also be
c         zero. Avoid calculating finite differences from relative
c         rate values that may by non-zero owing to convergence
c         tolerances and the like.
c
          do j = 1,kord + 1
            fdrr0(j,nrc) = 0.
          enddo
        else
c
c         The relative rate is calculated for the reactant, which
c         is actively reacting according to a rate law:
c
c           jreac =  0: set to react
c           jreac = -1: saturated, but the remaining reactant mass
c                         continues to react irreversibly
c
c         Hence calculate the finite differences.
c
          fdrr0(1,nrc) = (rrelr1(nrc) - rrelr0(nrc))/delxi
          do j = 2,kord + 1
            dfxu = fdlim*dxsm00(j)
            dfxl = -dfxu
            k = j - 1
            dfx = fdrr0(k,nrc) - fdrrm1(k,nrc)
            dfx = min(dfxu,dfx)
            dfx = max(dfxl,dfx)
            if (abs(dfx) .ge. dfxu) itrunc = min(itrunc,k)
            fdrr0(j,nrc) = dfx/dxsm00(j)
          enddo
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(2) .gt. 0) then
c
c       Compute finite differences for the inverse rate.
c
        fdri0(1) = (rirec1 - rirec0)/delxi
        do j = 2,kord + 1
          dfxu = fdlim*dxsm00(j)
          dfxl = -dfxu
          k = j - 1
          dfx = fdri0(k) - fdrim1(k)
          dfx = min(dfxu,dfx)
          dfx = max(dfxl,dfx)
          if (abs(dfx) .ge. dfxu) itrunc = min(itrunc,k)
          fdri0(j) = dfx/dxsm00(j)
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute finite differences for the pH.
c
      fdph0(1) = (ph1 - ph0)/delxi
      do j = 2,kord
        dfxu = fdlim*dxsm00(j)
        dfxl = -dfxu
        k = j - 1
        dfx = fdph0(k) - fdphm1(k)
        dfx = min(dfxu,dfx)
        dfx = max(dfxl,dfx)
        if (abs(dfx) .ge. dfxu) itrunc = min(itrunc,k)
        fdph0(j) = dfx/dxsm00(j)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute finite differences for the Eh.
c
      if (qredox) then
        fdeh0(1) = (eh1 - eh0)/delxi
        do j = 2,kord
          dfxu = fdlim*dxsm00(j)
          dfxl = -dfxu
          k = j - 1
          dfx = fdeh0(k) - fdehm1(k)
          dfx = min(dfxu,dfx)
          dfx = max(dfxl,dfx)
          if (abs(dfx) .ge. dfxu) itrunc = min(itrunc,k)
          fdeh0(j) = dfx/dxsm00(j)
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute finite differences for the log fO2.
c
      if (qredox) then
        fdo20(1) = (fo2lg1 - fo2lg0)/delxi
        do j = 2,kord
          dfxu = fdlim*dxsm00(j)
          dfxl = -dfxu
          k = j - 1
          dfx = fdo20(k) - fdo2m1(k)
          dfx = min(dfxu,dfx)
          dfx = max(dfxl,dfx)
          if (abs(dfx) .ge. dfxu) itrunc = min(itrunc,k)
          fdo20(j) = dfx/dxsm00(j)
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute finite differences for the activity of water.
c
      fdaw0(1) = (aw1 - aw0)/delxi
      do j = 2,kord
        dfxu = fdlim*dxsm00(j)
        dfxl = -dfxu
        k = j - 1
        dfx = fdaw0(k) - fdawm1(k)
        dfx = min(dfxu,dfx)
        dfx = max(dfxl,dfx)
        if (abs(dfx) .ge. dfxu) itrunc = min(itrunc,k)
        fdaw0(j) = dfx/dxsm00(j)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute finite differences for the affinities of the
c     various phases.
c
      do np = 1,npt
        fdaf0(1,np) = (affp(np) - affp0(np))/delxi
        do j = 2,kord
          dfxu = fdlim*dxsm00(j)
          dfxl = -dfxu
          k = j - 1
          dfx = fdaf0(k,np) - fdafm1(k,np)
          dfx = min(dfxu,dfx)
          dfx = max(dfxl,dfx)
          if (abs(dfx) .ge. dfxu) itrunc = min(itrunc,k)
          fdaf0(j,np) = dfx/dxsm00(j)
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute finite differences for the number of moles of phases
c     present in the equilibrium system.
c
      do npe = 1,npet
        fdpe0(1,npe) = (emop(npe) - emop0(npe))/delxi
        do j = 2,kord
          dfxu = fdlim*dxsm00(j)
          dfxl = -dfxu
          k = j - 1
          dfx = fdpe0(k,npe) - fdpem1(k,npe)
          dfx = min(dfxu,dfx)
          dfx = max(dfxl,dfx)
          if (abs(dfx) .ge. dfxu) itrunc = min(itrunc,k)
          fdpe0(j,npe) = dfx/dxsm00(j)
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute finite differences for the number of moles of
c     selected species present in the equilibrium system. These
c     species include H2O(l) only for the aqueous phase and
c     all species of the non-aqueous phases.
c
      do nse = 1,nset
        fdse0(1,nse) = (emos(nse) - emos0(nse))/delxi
        do j = 2,kord
          dfxu = fdlim*dxsm00(j)
          dfxl = -dfxu
          k = j - 1
          dfx = fdse0(k,nse) - fdsem1(k,nse)
          dfx = min(dfxu,dfx)
          dfx = max(dfxl,dfx)
          if (abs(dfx) .ge. dfxu) itrunc = min(itrunc,k)
          fdse0(j,nse) = dfx/dxsm00(j)
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(2) .gt. 0) then
c
c       Compute finite differences for reaction rates.
c
        do nrc = 1,nrct
          if (jreac(nrc).eq.1 .or. jreac(nrc).eq.2) then
c
c           The reaction rate must be a constant zero, because no
c           reactant mass remains:
c
c             jreac =  1: exhausted
c             jreac =  2: saturated; any remaining reactant mass is
c                           converted to the corresponding product
c                           phase, so the "reactant" is effectively
c                           exhausted.
c
c           Hence the corresponding finite differences must also be
c           zero. Avoid calculating finite differences from reaction
c           rate values that may by non-zero owing to convergence
c           tolerances and the like.
c
            do j = 1,kord
              fdre0(j,nrc) = 0.
            enddo
          else
c
c           The reaction rate is calculated for the reactant, which
c           is actively reacting according to a rate law:
c
c             jreac =  0: set to react
c             jreac = -1: saturated, but the remaining reactant mass
c                           continues to react irreversibly
c
c           Hence calculate the finite differences.
c
            fdre0(1,nrc) = (rreac1(nrc) - rreac0(nrc))/delxi
            do j = 2,kord
              dfxu = fdlim*dxsm00(j)
              dfxl = -dfxu
              k = j - 1
              dfx = fdre0(k,nrc) - fdrem1(k,nrc)
              dfx = min(dfxu,dfx)
              dfx = max(dfxl,dfx)
              if (abs(dfx) .ge. dfxu) itrunc = min(itrunc,k)
              fdre0(j,nrc) = dfx/dxsm00(j)
            enddo
          endif
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute finite differences for the affinities of irreversible
c     reactions.
c
c     Zero the affinities of exhausted reactants which correspond
c     to product minerals present in the ES. This will prevent
c     non-zero values due to convergence tolerances, etc., from
c     adversely affecting the corresponding finite differences for
c     these reactants.
c
      do nrc = 1,nrct
        if (jreac(nrc) .eq. 1) then
          if (jcode(nrc).eq.0 .or. jcode(nrc) .eq. 1) then
c
c           Note that this is done only for reactants which are
c           pure minerals or solid solutions.
c
            np = nrndex(nrc)
            if (jpflag(np) .eq. -1) then
              afrc0(nrc) = 0.
              afrc1(nrc) = 0.
            endif
          endif
        endif
      enddo
c
c     Now compute the finite differences.
c
      do nrc = 1,nrct
        if (jreac(nrc).eq.-1 .or. jreac(nrc).eq.2) then
c
c         The affinity is zero by definition:
c
c           jreac = -1: saturated, but the remaining reactant mass
c                         continues to react irreversibly
c           jreac =  2: saturated; any remaining reactant mass is
c                         converted to the corresponding product phase,
c                         so the "reactant" is effectively exhausted.
c
c         Hence the corresponding finite differences must also be zero.
c         The corresponding computed affinities may be non-zero, owing
c         to finite convergence tolerances. Avoid calculating finite
c         differences from such computed affinity values.
c
          do j = 1,kord
            fdar0(j,nrc) = 0.
          enddo
        else
c
c         The affinity is generally a non-zero value calculated from
c         a rate law:
c
c           jreac =  0: set to react
c           jreac =  1: exhausted
c
c         Hence calculate the finite differences.
c
          fdar0(1,nrc) = (afrc1(nrc) - afrc0(nrc))/delxi
          do j = 2,kord
            dfxu = fdlim*dxsm00(j)
            dfxl = -dfxu
            k = j - 1
            dfx = fdar0(k,nrc) - fdarm1(k,nrc)
            dfx = min(dfxu,dfx)
            dfx = max(dfxl,dfx)
            if (abs(dfx) .ge. dfxu) itrunc = min(itrunc,k)
            fdar0(j,nrc) = dfx/dxsm00(j)
          enddo
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (itrunc .lt. kord) then
        if (iodb(1) .ge. 1) write (noutpt,1000) itrunc
 1000   format(/' * Note - (EQ6/stepfd) Cutting the order to ',i2,
     $  /7x,'to stay within the finite difference limit.')
        kord = itrunc
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make the current point the new base point.
c
  200 km10 = km1
      kmt0 = kmt
      kx10 = kx1
      kxt0 = kxt
      kdim0 = kdim
c
      xim1 = xi0
      xi0 = xi1
c
      time0 = time1
      tempc0 = tempc
c
      ph0 = ph1
      eh0 = eh1
      fo2lg0 = fo2lg1
      aw0 = aw1
c
      fje0 = fje
      fxi0 = fxi
      sigmm0 = sigmam
c
      call copyia(iindx1,iindx0,kdim)
      call copyia(ipndx1,ipndx0,kdim)
      call copyca(uzvec1,uzvec0,kdim)
c
      call copyaa(zvclg1,zvclg0,kdim)
      call copyaa(zvec1,zvec0,kdim)
c
      call copyaa(moph,moph0,nptmax)
      call copyaa(mosp,mosp0,nstmax)
      call copyaa(acflg,acflg0,nstmax)
c
      call copyaa(sfcar,sfcar0,nrct)
      call copyia(jreac,jreac0,nrct)
      call copyaa(xirct,xirct0,nrct)
      call copyaa(morr,morr0,nrct)
      call copyaa(modr,modr0,nrct)
c
      call copyaa(mtb,mtb0,nbt)
c
      call copyaa(affp,affp0,npt)
c
      npet0 = npet
      nset0 = nset
c
      call copyia(iemop,iemop0,npet)
      call copyaa(emop,emop0,npet)
c
      nmax = 2*npetmx
      call copyia(ncmpe,ncmpe0,nmax)
c
      call copyia(iemos,iemos0,nset)
      call copyaa(emos,emos0,nset)
c
      call copyaa(afrc1,afrc0,nrct)
c
      call copyaa(rrelr1,rrelr0,nrct)
c
      if (iopt(2) .gt. 0) then
        rirec0 = rirec1
        call copyaa(rreac1,rreac0,nrct)
      endif
c
      if (npts .gt. 2) then
c
c       Save the old finite differences.
c
        nmax = nrd1mx*kmax
        call copyaa(fdzv0,fdzvm1,nmax)
c
        nmax = nrd1mx*nrctmx
        call copyaa(fdrr0,fdrrm1,nmax)
c
        if (iopt(2) .gt. 0) then
          call copyaa(fdri0,fdrim1,nrd1mx)
        endif
c
        call copyaa(fdaw0,fdawm1,nordmx)
        call copyaa(fdeh0,fdehm1,nordmx)
        call copyaa(fdo20,fdo2m1,nordmx)
        call copyaa(fdph0,fdphm1,nordmx)
c
        nmax = nordmx*nptmax
        call copyaa(fdaf0,fdafm1,nmax)
c
        nmax = nordmx*npetmx
        call copyaa(fdpe0,fdpem1,nmax)
c
        nmax = nordmx*nsetmx
        call copyaa(fdse0,fdsem1,nmax)
c
        nmax = nordmx*nrctmx
        call copyaa(fdar0,fdarm1,nmax)
c
        if (iopt(2) .gt. 0) then
          nmax = nordmx*nrctmx
          call copyaa(fdre0,fdrem1,nmax)
        endif
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
