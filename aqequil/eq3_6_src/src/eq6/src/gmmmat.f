      subroutine gmmmat(act,afrc1,cdac,cdacb,cdrs,csigma,eps100,
     $ fkrc,idirec,iindx1,iktmax,imchmx,imech,jcode,jreac,kbt,kdim,
     $ kmax,mmmatr,nbasp,nbt,nbtmax,ncmpr,ndac,ndacb,ndact,ndctmx,
     $ ndrs,ndrsmx,ndrsr,noutpt,nptmax,nrct,nrctmx,nrk,nrndex,nstmax,
     $ nttyo,nxridx,nxrtmx,rk,rtcnst,rxbar,sfcar,ureac,xlks)
c
c     This subroutine calculates the matrix M (mmmatr) used in turn to
c     calculate the Jacobian matrix J[r] (armatr) used by the higher-
c     order (stiff) ODE integrator. The dependence of J[r] on exact
c     kinetic rate law expressions (coded in EQ6\crrate.f) is handled
c     through the M matrix. This matrix expresses the dependency of the
c     rate laws on the thermodynamic activities of species, which
c     presently must all be of type aqueous.
c
c       The present subroutine must be tightly coordinated with
c     EQ6\ccrate.f in order to maintain proper consistency and thus
c     the correctness of the calculated M matrix. If such coordination
c     and consistency is compromised, the higher-order (stiff) ODE
c     solver may fail to converge.
c
c     This subroutine is called by:
c
c       EQ6/garmat.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c
c     Principal output:
c
c       mmmatr = the matrix M, which is used to calculate the Jacobian
c                  matrix J[r] used by the higher-order (stiff) ODE
c                  integrator
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
     $ nstmax,nxrtmx
c
      integer noutpt,nttyo
c
      integer kbt,kdim,nbt,nrct
c
      integer idirec(nrctmx),iindx1(kmax),imech(2,nrctmx),jcode(nrctmx),
     $ jreac(nrctmx),nbasp(nbtmax),ndac(ndctmx,imchmx,2,nrctmx),
     $ ncmpr(2,nptmax),ndacb(nbt,imchmx,2,nrct),ndact(imchmx,2,nrctmx),
     $ ndrs(ndrsmx),ndrsr(2,nstmax),nrk(2,nrctmx),nrndex(nrctmx),
     $ nxridx(nrctmx)
c
      character*24 ureac(nrctmx)
c
      real(8) act(nstmax),afrc1(nrctmx),cdac(ndctmx,imchmx,2,nrctmx),
     $ cdacb(nbt,imchmx,2,nrct),cdrs(ndrsmx),csigma(imchmx,2,nrctmx),
     $ fkrc(nrctmx),mmmatr(nrct,kmax),rk(imchmx,2,nrctmx),sfcar(nrctmx),
     $ rxbar(iktmax,nxrtmx),xlks(nstmax)
c
      real(8) eps100,rtcnst
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ik,j2,j3,k,n,nb,np,nr1,nr2,nrc,ns,nse,nxr
c
      integer ilnobl
c
      logical qovstp
c
      character(len=8) ux8
c
      real(8) affr,aprod,cx,cxs,efx,fs,mijk,skfq,thetij,vij,xk,xlk
c
      real(8) coefdr,texp
c
c-----------------------------------------------------------------------
c
c     First zero the matrix M.
c
      do nrc = 1,nrct
        do k = 1,kmax
          mmmatr(nrc,k) = 0.
        enddo
      enddo
c
c     Loop on reactants (rows of M).
c
      do nrc = 1,nrct
c
        affr = afrc1(nrc)
        qovstp = affr.lt.0. .and. nrk(2,nrc).eq.0 .and. jreac(nrc).le.0
        if (jreac(nrc).le.0 .and.
     $    (abs(affr).gt.eps100 .or. qovstp))then
c
c         Have an active reactant. It is not exhausted or formally
c         saturated as indicated by its jreac flag. Also, it either
c         has a finite affinity or the overstep condition pertains.
c         Compute the corresponding row of M. For a reactant not
c         satisfying these conditions, the rate is taken to be zero,
c         and the corresponding row of M should be left filled with
c         zeroes.
c
          fs = fkrc(nrc)*sfcar(nrc)
c
          if (idirec(nrc) .eq. 1) then
c
c           The rate law expression used follows the net forward rate
c           format.
c
            if (nrk(1,nrc) .eq. 1) then
c
c             Have a specified relative rate (arbitrary kinetics).
c             Leave the current row of M filled with zeroes.
c
              continue
            elseif (nrk(1,nrc) .eq. 2) then
c
c             Have the TST-like rate equation (multi-term). Compute
c             the contribution of each successive term to the current
c             row of M and accumulate all contributions.
c
c             First get the equilibrium constant K (xk).
c
              np = nrndex(nrc)
              xlk = 0.
              if (jcode(nrc) .eq. 0) then
                ns = ncmpr(1,np)
                xlk = xlks(ns)
              elseif (jcode(nrc) .eq. 1) then
                nxr = nxridx(nrc)
                ik = 0
                nr1 = ncmpr(1,np)
                nr2 = ncmpr(2,np)
                do ns = nr1,nr2
                  ik = ik + 1
                  xlk = xlk + rxbar(ik,nxr)*xlks(ns)
                enddo
              else
                write (noutpt,"(' Illegal jcode value in gmmmat.f')")
                write (nttyo,"(' Illegal jcode value in gmmmat.f')")
                stop
              endif
              xk = texp(xlk)
c
c             Loop over terms in the rate law.
c
              do i = 1,imech(1,nrc)
c
c               Calculate the kinetic activity product (aprod).
c
                aprod = 1.
                do n = 1,ndact(i,1,nrc)
                  ns = ndac(n,i,1,nrc)
                  aprod = aprod*act(ns)**cdac(n,i,1,nrc)
                enddo
c
                skfq = rk(i,1,nrc)*fs*aprod
                efx = affr/(csigma(i,1,nrc)*rtcnst)
                vij = skfq*(1.0 - exp(-efx))
c
                do k = 1,kbt
                  nb = iindx1(k)
                  nse = nbasp(nb)
                  if (jcode(nrc) .eq. 0) then
                    ns = ncmpr(1,np)
                    cx = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns,nstmax)
                  elseif (jcode(nrc) .eq. 1) then
                    nxr = nxridx(nrc)
                    cx = 0.
                    ik = 0
                    nr1 = ncmpr(1,np)
                    nr2 = ncmpr(2,np)
                    do ns = nr1,nr2
                      ik = ik + 1
                      cxs = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns,nstmax)
                      cx = cx + rxbar(ik,nxr)*cxs
                    enddo
                  else
                    write (noutpt,
     $              "(' Illegal jcode value in gmmmat.f')")
                    write (nttyo,"(' Illegal jcode value in gmmmat.f')")
                    stop
                  endif
                  thetij = -skfq*exp(-efx)/csigma(i,1,nrc)
                  mijk = -vij*ndacb(nse,i,1,nrc) + cx*thetij
                  mmmatr(nrc,k) = mmmatr(nrc,k) + mijk
                enddo
              enddo
c
            elseif (nrk(1,nrc) .eq. 3) then
c
c             Have a specified rate (no affinity dependence).
c             Leave the current row of M filled with zeroes.
c
            else
c
c             Unrecognized rate law code.
c
              j2 = ilnobl(ureac(nrc))
              write (ux8,'(i5)') nrk(1,nrc)
              call lejust(ux8)
              j3 = ilnobl(ux8)
              write (noutpt,1000) ureac(nrc)(1:j2),ux8(1:j3)
              write (nttyo,1000) ureac(nrc)(1:j2),ux8(1:j3)
 1000         format(/' * Error - (EQ6/gmmmat) The reactant ',a,
     $        ' has an',/7x,'unrecognized forward rate law code'
     $        ' of ',a,'.')
              stop
            endif
          else
c
c           The rate law expression used follows the net backward rate
c           format.
c
            if (nrk(2,nrc) .eq. 1) then
c
c             Have specified the case of instantaneous equilibrium.
c             For now, leave the current row of M filled with zeroes.
c
              continue
            elseif (nrk(2,nrc) .eq. 2) then
c
c             Have the TST-like rate equation (multi-term). Compute
c             the contribution of each successive term to the current
c             row of M and accumulate all contributions.
c
c             First get the equilibrium constant K (xk).
c
              np = nrndex(nrc)
              xlk = 0.
              if (jcode(nrc) .eq. 0) then
                ns = ncmpr(1,np)
                xlk = xlks(ns)
              elseif (jcode(nrc) .eq. 1) then
                nxr = nxridx(nrc)
                ik = 0
                nr1 = ncmpr(1,np)
                nr2 = ncmpr(2,np)
                do ns = nr1,nr2
                  ik = ik + 1
                  xlk = xlk + rxbar(ik,nxr)*xlks(ns)
                enddo
              else
                write (noutpt,"(' Illegal jcode value in gmmmat.f')")
                write (nttyo,"(' Illegal jcode value in gmmmat.f')")
                stop
              endif
              xk = texp(xlk)
c
c             Loop over terms in the rate law.
c
              do i = 1,imech(2,nrc)
c
c               Calculate the kinetic activity product (aprod).
c
                aprod = 1.
                do n = 1,ndact(i,2,nrc)
                  ns = ndac(n,i,2,nrc)
                  aprod = aprod*act(ns)**cdac(n,i,2,nrc)
                enddo
c
                skfq = rk(i,2,nrc)*fs*aprod
                efx = affr/(csigma(i,2,nrc)*rtcnst)
                vij = skfq*(1.0 - exp(efx))
c
                do k = 1,kbt
                  nb = iindx1(k)
                  nse = nbasp(nb)
                  if (jcode(nrc) .eq. 0) then
                    ns = ncmpr(1,np)
                    cx = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns,nstmax)
                  elseif (jcode(nrc) .eq. 1) then
                    nxr = nxridx(nrc)
                    cx = 0.
                    ik = 0
                    nr1 = ncmpr(1,np)
                    nr2 = ncmpr(2,np)
                    do ns = nr1,nr2
                      ik = ik + 1
                      cxs = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,ns,nstmax)
                      cx = cx + rxbar(ik,nxr)*cxs
                    enddo
                  else
                    write (noutpt,
     $              "(' Illegal jcode value in gmmmat.f')")
                    write (nttyo,"(' Illegal jcode value in gmmmat.f')")
                    stop
                  endif
                  thetij = -skfq*exp(efx)/(csigma(i,2,nrc))
                  mijk = -vij*ndacb(nse,i,2,nrc) + cx*thetij
                  mmmatr(nrc,k) = mmmatr(nrc,k) + mijk
                enddo
              enddo
c
            elseif (nrk(2,nrc) .eq. 3) then
c
c             Have a specified rate (no affinity dependence).
c             Leave the current row of M filled with zeroes.
c
            else
c
c             Unrecognized rate law code.
c
              j2 = ilnobl(ureac(nrc))
              write (ux8,'(i5)') nrk(2,nrc)
              call lejust(ux8)
              j3 = ilnobl(ux8)
              write (noutpt,1020) ureac(nrc)(1:j2),ux8(1:j3)
              write (nttyo,1020) ureac(nrc)(1:j2),ux8(1:j3)
 1020         format(/' * Error - (EQ6/gmmmat) The reactant ',a,
     $        ' has an',/7x,'unrecognized backward rate law code'
     $        ' of ',a,'.')
              stop
            endif
          endif
        endif
      enddo
c
      end
