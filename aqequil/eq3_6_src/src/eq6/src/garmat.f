      subroutine garmat(act,afrc1,aimatr,al10,armatr,cdac,cdacb,
     $ cdrs,csigma,csts,delvec,dlogxw,dvjdte,eact,eps100,fkrc,gmmatr,
     $ hact,iact,idirec,iindx1,iktmax,imchmx,imech,ipivot,jcode,jreac,
     $ jtemp,kbt,kdim,kmax,mmmatr,morr,mwtrc,nbasp,nbt,nbtmax,ncmpr,
     $ ndac,ndacb,ndact,ndctmx,ndrs,ndrsmx,ndrsr,noutpt,nptmax,nrct,
     $ nrctmx,nrct1,nrk,nrndex,nsk,nstmax,nsts,nstsmx,nstsr,nttkmx,
     $ nttyo,nxridx,nxrtmx,rirec1,rk,rkb,rreac1,rrelr1,rrxfi1,rtcnst,
     $ rxbar,sfcar,sgmatr,ssfcar,tempc,tempcb,tempk,ttk,ureac,whcfac,
     $ xi1,xlks,xxmatr,xymatr)
c
c     This subroutine calculates the Jacobian matrix J[r] (armatr)
c     used by the higher-order (stiff) ODE integrator.
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
c       armatr = the Jacobian matrix J[r]
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer iktmax,imchmx,kmax,nbtmax,ndctmx,ndrsmx,nrctmx,nptmax,
     $ nstmax,nstsmx,nttkmx,nxrtmx
c
      integer noutpt,nttyo
c
      integer jtemp,kbt,kdim,nbt,nrct,nrct1
c
      integer iact(imchmx,2,nrctmx),idirec(nrctmx),iindx1(kmax),
     $ imech(2,nrctmx),ipivot(kmax),jcode(nrctmx),jreac(nrctmx),
     $ nbasp(nbtmax),ndac(ndctmx,imchmx,2,nrctmx),ncmpr(2,nptmax),
     $ ndacb(nbt,imchmx,2,nrct),ndact(imchmx,2,nrctmx),
     $ ndrs(ndrsmx),ndrsr(2,nstmax),nrk(2,nrctmx),nrndex(nrctmx),
     $ nsk(nrctmx),nsts(nstsmx),nstsr(2,nstmax),nxridx(nrctmx)
c
      character*24 ureac(nrctmx)
c
      real(8) act(nstmax),afrc1(nrctmx),aimatr(kmax,kmax),
     $ armatr(nrct1,nrct1),cdac(ndctmx,imchmx,2,nrctmx),
     $ cdacb(nbt,imchmx,2,nrct),cdrs(ndrsmx),csigma(imchmx,2,nrctmx),
     $ csts(nstsmx),delvec(kmax),dlogxw(nbtmax),dvjdte(nrct),
     $ eact(imchmx,2,nrctmx),fkrc(nrctmx),gmmatr(kmax,kmax),
     $ hact(imchmx,2,nrctmx),mmmatr(nrct,kmax),morr(nrctmx),
     $ mwtrc(nrctmx),rk(imchmx,2,nrctmx),rkb(imchmx,2,nrctmx),
     $ rreac1(nrctmx),rrelr1(nrctmx),rrxfi1(imchmx,nrctmx),
     $ rxbar(iktmax,nxrtmx),sfcar(nrctmx),sgmatr(nrct,kmax),
     $ ssfcar(nrctmx),ttk(nttkmx),xlks(nstmax),xxmatr(kmax,nrct),
     $ xymatr(nrct,nrct1)
c
      real(8) al10,eps100,rirec1,rtcnst,tempc,tempcb,tempk,whcfac,xi1
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,id,idj,idn,ik,j,k,kcol,kk,krow,n,nb,np,nrc,nr1,nr2,
     $ ns,nxr
c
      real(8) ax,cx,cxs,dxy,rt2,sx,vx,xx
c
      real(8) coefst
c
c-----------------------------------------------------------------------
c
c     Zero the matrix J[r] (armatr).
c
      do j = 1,nrct1
        do k = 1,nrct1
          armatr(j,k) = 0.
        enddo
      enddo
c
c     Calculate the matrix M (mmmatr). This incorporates the
c     dependency of rate laws on the thermodynamic activities of
c     species (currently all such species must be of type aqueous).
c
      call gmmmat(act,afrc1,cdac,cdacb,cdrs,csigma,eps100,
     $ fkrc,idirec,iindx1,iktmax,imchmx,imech,jcode,jreac,kbt,kdim,
     $ kmax,mmmatr,nbasp,nbt,nbtmax,ncmpr,ndac,ndacb,ndact,ndctmx,
     $ ndrs,ndrsmx,ndrsr,noutpt,nptmax,nrct,nrctmx,nrk,nrndex,nstmax,
     $ nttyo,nxridx,nxrtmx,rk,rtcnst,rxbar,sfcar,ureac,xlks)
c
c     Now calculate the matrix SIGMA (sgmatr) from matrix M and the
c     array W-squiggle (dlogxw) array.
c
      do nrc = 1,nrct
        do k = 1,kmax
          sgmatr(nrc,k) = 0.
        enddo
      enddo
c
      do nrc = 1,nrct
        if (jreac(nrc) .le. 0) then
          sx = 0.
          do kk = 2,kbt
            sx = sx + mmmatr(nrc,kk)
          enddo
          sgmatr(nrc,1) = al10*(mmmatr(nrc,1)*dlogxw(1) - sx)
          do k = 2,kbt
            sgmatr(nrc,k) = al10*(mmmatr(nrc,1)*dlogxw(k)
     $      + mmmatr(nrc,k))
          enddo
        endif
      enddo
c
c     Calculate the matrix PHI (aimatr), which is the inverse of
c     the Jacobian matrix J[z] used by the algebraic equation solver.
c
      call ginvrt(aimatr,delvec,gmmatr,ipivot,kdim,kmax)
c
c     Emulate the matrix multiplication PHI*THETA. Put the result
c     in a scratch matrix (xxmatr). Note that xxmatr must have zero
c     rows for reactants that are not pure minerals or solid solutions.
c     Those are currently the only types that may have rate dependencies
c     on activities. That might change in the future. Also, in the row-
c     column multiplication, contributions are included only over the
c     index range (1, kbt), not (1, kdim), which would be the general
c     case. That is because the SIGMA array has no non-zero columns
c     beyond the kbt-th. That in turn results because all of the
c     currently allowed rate dependencies on thermodynamic activities
c     are restricted to those of aqueous species.
c
      do krow = 1,kmax
        do nrc = 1,nrct
          xxmatr(krow,nrc) = 0.
        enddo
      enddo
c
      do krow = 1,kdim
        do nrc = 1,nrct
          if (jreac(nrc) .le. 0) then
            if (jcode(nrc) .eq. 0) then
c
c             Pure mineral.
c
              np = nrndex(nrc)
              ns = ncmpr(1,np)
              xx = 0.
              do kcol = 1,kbt
                nb = iindx1(kcol)
                cx = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
                xx = xx + aimatr(krow,kcol)*cx
              enddo
              xxmatr(krow,nrc) = xx
            elseif (jcode(nrc) .eq. 1) then
c
c             Solid solution.
c
              nxr = nxridx(nrc)
              np = nrndex(nrc)
              nr1 = ncmpr(1,np)
              nr2 = ncmpr(2,np)
              xx = 0.
              do kcol = 1,kbt
                nb = iindx1(kcol)
                cx = 0.
                ik = 0
                do ns = nr1,nr2
                  ik = ik + 1
                  cxs = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
                  cx = cx + rxbar(ik,nxr)*cxs
                enddo
                xx = xx + aimatr(krow,kcol)*cx
              enddo
              xxmatr(krow,nrc) = xx
            else
              xxmatr(krow,nrc) = 0.
            endif
          endif
        enddo
      enddo
c
c     Calculate the product SIGMA*PHI*THETA, storing the result in
c     scratch matrix xymatr. This operation does not put anything
c     in the last (nrct1-th) column of xymatr. That column is
c     associated with explicit time dependency of the rate
c     equations. It is only used if the r vector of "rates
c     represented by finite differences" is composed of relative
c     rates plus the inverse rate.
c
      do nrc = 1,nrct
        do j = 1,nrct1
          xymatr(nrc,j) = 0.
        enddo
      enddo
c
      do nrc = 1,nrct
        if (jreac(nrc) .le. 0) then
          do j = 1,nrct
            xx = 0.
            do k = 1,kbt
              xx = xx + sgmatr(nrc,k)*xxmatr(k,j)
            enddo
            xymatr(nrc,j) = xx
          enddo
        endif
      enddo
c
c     Emulate the matrix addition SIGMA*PHI*THETA + ZETA.
c     The ZETA matrix addresses the dependencies of rates on
c     surface areas. This coordinates with EQ6\csfar.f, which
c     evaluates surface area models. Here again, the last (nrct1-th)
c     column of xymatr is not used.
c
      do nrc = 1,nrct
        if (jreac(nrc) .le. 0) then
c
          id = idirec(nrc)
          if (nrk(id,nrc) .eq. 1) then
c
c           The rate law is TST-type, therefore there is a surface
c           area dependency.
c
            if (jcode(nrc).eq.0 .or. jcode(nrc).eq.1) then
              np = nrndex(nrc)
              n = ndrsr(1,ns)
              cx = cdrs(n)
            else
              write (noutpt,"(' Illegal jcode value in garmat.f')")
              write (nttyo,"(' Illegal jcode value in garmat.f')")
              stop
            endif
c
            sx = 0.
            if (nsk(nrc) .eq. 0) then
c
c             Constant surface area.
c
              continue
            elseif (nsk(nrc) .eq. 1) then
c
c             Constant specific surface area.
c
              if (sfcar(nrc) .gt. 0.)
     $        sx = rreac1(nrc)*cx*ssfcar(nrc)*mwtrc(nrc)/sfcar(nrc)
            elseif (nsk(nrc) .eq. 2) then
c
c             Constant particle number surface area growth law.
c
              if (morr(nrc) .gt. 0.)
     $        sx = 2.*rreac1(nrc)*cx/(3.*morr(nrc))
            else
              write (noutpt,"(' Illegal nsk value in garmat.f')")
              write (nttyo,"(' Illegal nsk value in garmat.f')")
              stop
            endif
            xymatr(nrc,nrc) = xymatr(nrc,nrc) + sx
          endif
        endif
      enddo
c
c     Overlay for variable temperature. This coordinates with
c     EQ6\evratc.f, which calculates rate constants as a function
c     of temperature, and with EQ6\gtemp.f, which calculates the
c     temperature as a specified function of either Xi or t.
c
      do nrc = 1,nrct
        dvjdte(nrc) = 0.
      enddo
c
      if (jtemp .gt. 0) then
c
c       Calculate temperature derivatives of the rates, following
c       through the temperature dependency of the rate constants.
c       This is needed only if the temperature is a function of
c       Xi or t.
c
        rt2 = rtcnst*tempk
        do nrc = 1,nrct
          if (jreac(nrc) .le. 0) then
            id = idirec(nrc)
            if (nrk(id,nrc) .ge. 2) then
              do i = 1,imech(id,nrc)
                if (iact(i,id,nrc) .eq. 0) then
c
c                 No temperature dependence.
c
                  continue
                elseif (iact(i,id,nrc) .eq. 1) then
c
c                 Constant activation energy.
c
                  dvjdte(nrc) = dvjdte(nrc)
     $            + rrxfi1(i,nrc)*eact(i,id,nrc)/rt2
                elseif (iact(i,id,nrc) .eq. 2) then
c
c                 Constant activation enthalpy.
c
                  dvjdte(nrc) = dvjdte(nrc)
     $            + rrxfi1(i,nrc)*eact(i,id,nrc)/rt2
                endif
              enddo
            endif
          endif
        enddo
      endif
c
      if (jtemp .eq. 0) then
c
c       Constant temperature.
c
        continue
      elseif (jtemp .eq. 1) then
c
c       Temperature is a linear function of Xi.
c       tempc = tempcb + ttk(1)*Xi
c
        do nrc = 1,nrct
          id = idirec(nrc)
          if (jreac(nrc).le.0 .and. nrk(id,nrc).ge.2) then
            dxy = dvjdte(nrc)*ttk(1)
            do j = 1,nrct
              idj = idirec(j)
              if (jreac(j).le.0 .and. nrk(idj,j).ge.2) then
                xymatr(nrc,j) = xymatr(nrc,j) + dxy
              endif
            enddo
          endif
        enddo
      elseif (jtemp .eq. 2) then
c
c       Temperature is a linear function of time.
c       tempc = tempcb + ttk(1)*time1
c
        do nrc = 1,nrct
          id = idirec(nrc)
          if (jreac(nrc).le.0 .and. nrk(id,nrc).ge.2) then
            dxy = dvjdte(nrc)*ttk(1)
            xymatr(nrc,nrct1) = xymatr(nrc,nrct1) + dxy
          endif
        enddo
      elseif (jtemp .eq. 3) then
c
c       Fluid mixing tracking.
c       tempc = (tempcb*ttk(1) + Xi*ttk(2))/(Xi + ttk(1))
c
        do nrc = 1,nrct
          id = idirec(nrc)
          if (jreac(nrc).le.0 .and. nrk(id,nrc).ge.2) then
            dxy = dvjdte(nrc)*((ttk(2) - tempc)/(xi1 + ttk(1)))
            do j = 1,nrct
              idj = idirec(j)
              if (jreac(j).le.0 .and. nrk(idj,j).ge.2) then
                xymatr(nrc,j) = xymatr(nrc,j) + dxy
              endif
            enddo
          endif
        enddo
c
      endif
c
c     Multiply xymatr by the scalar w (whcfac).
c
      do nrc = 1,nrct
        do j = 1,nrct1
          xymatr(nrc,j) = whcfac*xymatr(nrc,j)
        enddo
      enddo
c
c     Emulate the matrix multiplication V*xymatr, which yields
c     the equivalent of J[r] + I. Note that V is nrct1 by nrct,
c     while xymatr is nrct by nrct1. Thus, J[r] and the I in question
c     are nrct1 by nrct1.
c
      do nrc = 1,nrct
        id = idirec(nrc)
        if (jreac(nrc) .le. 0) then
          if (nrk(id,nrc) .eq. 1) then
c
c           The current row of V corresponds to a reactant constrained
c           by a specified relative rate. V(nrc,nrc) = 1.0 and all
c           other elements of this row are zero.
c
            armatr(nrc,j) = xymatr(nrc,j)
          elseif (nrk(id,nrc).ge.2) then
c
c           The current row of V corresponds to a reactant constrained
c           by an absolute rate.
c
            do j = 1,nrct
              ax = 0.
              do n = 1,nrct
                idn = idirec(n)
                if (jreac(n) .le. 0) then
                  if (nrk(idn,n) .ge. 2) then
                    if (n .eq. nrc) then
                      vx = rirec1*(1.0 - rrelr1(nrc))
                    else
                      vx = -rirec1*rrelr1(nrc)
                    endif
                    ax = ax + vx*xymatr(n,j)
                  endif
                endif
                armatr(nrc,j) = ax
              enddo
            enddo
          endif
        endif
      enddo
c
c     Now get the nrct1-th row.
c
      vx = -rirec1**2
      do j = 1,nrct
        ax = 0.
        do n = 1,nrct
          idn = idirec(n)
          if (jreac(n) .le. 0) then
            if (nrk(idn,n) .ge. 2) then
              ax = ax + vx*xymatr(n,j)
            endif
          endif
          armatr(nrct1,j) = ax
        enddo
      enddo
c
c     Finish getting J[r] by substracting the unit matrix.
c
      do nrc = 1,nrct1
        armatr(nrc,nrc) = armatr(nrc,nrc) - 1.0
      enddo
c
      end
