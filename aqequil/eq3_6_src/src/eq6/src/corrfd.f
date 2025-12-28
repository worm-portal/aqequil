      subroutine corrfd(delxi,dxsm11,fdlim,fdre0,fdre1,fdri0,
     $ fdri1,fdrr0,fdrr1,iodb,iopt,jreac,nodbmx,noptmx,nord,
     $ nordmx,noutpt,npts,nrct,nrctmx,nrd1mx,rirec0,rirec1,
     $ rreac0,rreac1,rrelr0,rrelr1)
c
c     This subroutine computes finite differences for use in ODE
c     corrector iteration. These finite differences are based at the
c     new point stepped to, as opposed to the point stepped from.
c     Finite differences based on the point stepped from are used to
c     generate predictor functions. Those finite differences are
c     computed by EQ6/stepfd.f and comprise a larger set. Here finite
c     differences are computed only for rate functions.
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
      integer nodbmx,noptmx,nordmx,nrctmx,nrd1mx
c
      integer noutpt
c
      integer iodb(nodbmx),iopt(noptmx),jreac(nrctmx)
c
      integer nord,npts,nrct
c
      real*8 dxsm11(nrd1mx),fdre0(nordmx,nrctmx),fdre1(nordmx,nrctmx),
     $ fdri0(nrd1mx),fdri1(nrd1mx),fdrr0(nrd1mx,nrctmx),
     $ fdrr1(nrd1mx,nrctmx)
c
      real*8 rreac0(nrctmx),rreac1(nrctmx),rrelr0(nrctmx),rrelr1(nrctmx)
c
      real*8 delxi,fdlim,rirec0,rirec1
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer itrunc,j,k,nmax,nrc
c
      real*8 dfx,dfxl,dfxu
c
c-----------------------------------------------------------------------
c
      if (npts .eq. 1) then
c
c       Zero all finite differences.
c
        nmax = nrd1mx*nrctmx
        call initaz(fdrr1,nmax)
c
        if (iopt(2) .gt. 0) then
          call initaz(fdri1,nrd1mx)
c
          nmax = nrd1mx*nrctmx
          call initaz(fdre1,nmax)
        endif
c
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Overflow protection is automatically activated in the code blocks
c     below in order to run on VAX machines with a small exponent range
c     (+/- 38) for real*8.
c
      itrunc = nord
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
          do j = 1,nord + 1
            fdrr1(j,nrc) = 0.
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
          fdrr1(1,nrc) = (rrelr1(nrc) - rrelr0(nrc))/delxi
          do j = 2,nord + 1
            dfxu = fdlim*dxsm11(j)
            dfxl = -dfxu
            k = j - 1
            dfx = fdrr1(k,nrc) - fdrr0(k,nrc)
            dfx = min(dfxu,dfx)
            dfx = max(dfxl,dfx)
            if (abs(dfx) .ge. dfxu) itrunc = min(itrunc,k)
            fdrr1(j,nrc) = dfx/dxsm11(j)
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
        fdri1(1) = (rirec1 - rirec0)/delxi
        do j = 2,nord + 1
          dfxu = fdlim*dxsm11(j)
          dfxl = -dfxu
          k = j - 1
          dfx = fdri1(k) - fdri0(k)
          dfx = min(dfxu,dfx)
          dfx = max(dfxl,dfx)
          if (abs(dfx) .ge. dfxu) itrunc = min(itrunc,k)
          fdri1(j) = dfx/dxsm11(j)
        enddo
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
            do j = 1,nord + 1
              fdre1(j,nrc) = 0.
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
            fdre1(1,nrc) = (rreac1(nrc) - rreac0(nrc))/delxi
            do j = 2,nord + 1
              dfxu = fdlim*dxsm11(j)
              dfxl = -dfxu
              k = j - 1
              dfx = fdre1(k,nrc) - fdre0(k,nrc)
              dfx = min(dfxu,dfx)
              dfx = max(dfxl,dfx)
              if (abs(dfx) .ge. dfxu) itrunc = min(itrunc,k)
              if (j .le. nord) fdre1(j,nrc) = dfx/dxsm11(j)
            enddo
          endif
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (itrunc .lt. nord) then
        if (iodb(1) .ge. 1) write (noutpt,1000) itrunc
 1000   format(/' * Note - (EQ6/corrfd) Cutting the order to ',i2,
     $  /7x,'to stay within the finite difference limit.')
        nord = itrunc
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
