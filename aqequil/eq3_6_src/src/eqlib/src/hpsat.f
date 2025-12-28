      subroutine hpsat(acflg,act,actlg,afcnst,affp,affs,apx,bpx,
     $ cdrs,eps100,iapxmx,ibpxmx,ier,iktmax,ixrn1,jflag,jpflag,
     $ jsflag,jsol,ncmpr,ndrs,ndrsmx,ndrsr,noutpt,np,nptmax,
     $ nstmax,nttyo,nxrn1,nxrn2,nxtmax,sidrsp,sidrph,uphase,
     $ uspec,wfac,xbar,xbarlg,xlks)
c
c     This subroutine calculates the most stable (least soluble)
c     composition of a given solid solution, given the composition
c     of the aqueous phase it is in equilibrium with.
c
c     This subroutine is called by:
c
c       EQ3NR/scripx.f
c       EQ6/satchk.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       apx    = array of temperature-independent coefficients for
c                  solid solution activity coefficients
c       bpx    = array of site-mixing parameters for solid solution
c                  activity coefficients
c       ixrn1  = start of solid solution phases in the list of phases
c       ixrn2  = end of solid solution phases in the list of phases
c       jsol   = array of activity coefficient model indices for solid
c                  solutions
c       ncmpr  = species range pointer array for phases
c       uphase = array of phase names
c       uspec  = array of species names
c       wterm  = array of temperature-dependent coefficients for
c                  solid solution activity coefficients
c
c     Principal output:
c
c       acflg  = array of activity coefficients
c       act    = array of activities
c       actlg  = array of log activity values
c       affp   = array of affinities for phases
c       affs   = array of affinities for species
c       sidrph = array of saturation indices for phases
c       sidrsp = array of saturation indices for species
c       xbar   = array of mole fractions
c       xbarlg = array of log mole fraction values
c
c     Local:
c
c       sisppu = saturation index of a pure end-member, normalized
c                  so that one mole of end-member is destroyed in the
c                  corresponding reaction
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
      integer iapxmx,ibpxmx,iktmax,ndrsmx,nptmax,nstmax,nxtmax
c
      integer noutpt,nttyo
c
      integer jflag(nstmax),jpflag(nptmax),jsflag(nstmax),jsol(nxtmax),
     $ ncmpr(2,nptmax),ndrs(ndrsmx),ndrsr(2,nstmax)
c
      integer ier,ixrn1,np,nxrn1,nxrn2
c
      character*48 uspec(nstmax)
      character*24 uphase(nptmax)
c
      real*8 acflg(nstmax),actlg(nstmax),act(nstmax),affp(nptmax),
     $ affs(nstmax),apx(iapxmx,nxtmax),bpx(ibpxmx,nxtmax),
     $ cdrs(ndrsmx),sidrsp(nstmax),sidrph(nptmax),wfac(iktmax,nxtmax),
     $ xbar(nstmax),xbarlg(nstmax),xlks(nstmax)
c
      real*8 afcnst,eps100
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer nactc,nrn1,nr1,nr2,ns,nx
c
      real*8 af,cx,si,siph,sisp,sisppu,stx,stxi,stxm1,xl,xqk,xqks,
     $ xqksum,xx
c
      real*8 texp,tlg
c
c-----------------------------------------------------------------------
c
      ier = 0
      nx = np - ixrn1 + 1
c
      nr1 = ncmpr(1,np)
      nr2 = ncmpr(2,np)
c
      affp(np) = -9999999.
      sidrph(np) = -9999999.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the affinities and saturation indices of the pure
c     end-members.
c
      do ns = nr1,nr2
        xbar(ns) = 1.0
        xbarlg(ns) = 0.
        acflg(ns) = 0.
        act(ns) = 1.0
        actlg(ns) = 0.
c
        call afcalc(actlg,af,afcnst,cdrs,jflag,jsflag,ndrs,ndrsmx,
     $  ndrsr,ns,nstmax,si,xlks)
        affs(ns) = af
        sidrsp(ns) = si
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Count the number of active components.
c
      nactc = 0
      do ns = nr1,nr2
        if (jsflag(ns) .le. 0) then
          nactc = nactc + 1
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check the value of the site-mixing parameter.
c
      if (bpx(1,nx) .le. 0.) bpx(1,nx) = 1.0
      stx = bpx(1,nx)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (abs(stx - 1.0) .gt. eps100) go to 200
c
c     Calculate the mole fractions for the case of an ideal molecular
c     mixing solution.
c
      xqksum = 0.
      do ns = nr1,nr2
        if (jsflag(ns) .le. 0) then
          nrn1 = ndrsr(1,ns)
          cx = (-1.0/cdrs(nrn1))
          sisppu = cx*sidrsp(ns)
          xqk = texp(sisppu)
          xqksum = xqksum + xqk
        endif
      enddo
c
      do ns = nr1,nr2
        if (jsflag(ns) .le. 0) then
          nrn1 = ndrsr(1,ns)
          cx = (-1.0/cdrs(nrn1))
          sisppu = cx*sidrsp(ns)
          xqk = texp(sisppu)
          xx = xqk/xqksum
          xbar(ns) = xx
          xbarlg(ns) = tlg(xx)
        endif
      enddo
c
c     Calculate the activities, activity coefficients, affinities,
c     and saturation indices for the case of an ideal molecular
c     mixing solution.
c
      siph = 0.
      do ns = nr1,nr2
        if (jsflag(ns) .le. 0) then
          nrn1 = ndrsr(1,ns)
          cx = (-1.0/cdrs(nrn1))
          sisppu = cx*sidrsp(ns)
          xx = xbar(ns)
          xl = xbarlg(ns)
          sisp = sisppu - xl
          acflg(ns) = 0.
          act(ns) = xx
          actlg(ns) = xl
          sisp = sisp/cx
          sidrsp(ns) = sisp
          affs(ns) = afcnst*sisp
          siph = siph + xx*sisp
        endif
      enddo
      sidrph(np) = siph
      affp(np) = afcnst*siph
c
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  200 continue
c
c     Calculate the mole fractions for the case of an ideal site-mxing
c     solution.
c
      stxi = 1./bpx(1,nx)
      stxm1 = stx - 1.0
c
      xqksum = 0.
      do ns = nr1,nr2
        if (jsflag(ns) .le. 0) then
          nrn1 = ndrsr(1,ns)
          cx = (-1.0/cdrs(nrn1))
          sisppu = cx*sidrsp(ns)
          xqk = texp(sisppu)
          xqks = xqk**stxi
          xqksum = xqksum + xqks
        endif
      enddo
c
      do ns = nr1,nr2
        if (jsflag(ns) .le. 0) then
          nrn1 = ndrsr(1,ns)
          cx = (-1.0/cdrs(nrn1))
          sisppu = cx*sidrsp(ns)
          xqk = texp(sisppu)
          xqks = xqk**stxi
          xx = xqks/xqksum
          xbar(ns) = xx
          xbarlg(ns) = tlg(xx)
        endif
      enddo
c
c     Calculate the activities, activity coefficients, affinities,
c     and saturation indices for the case of an ideal site-mixing
c     solution.
c
      siph = 0.
      do ns = nr1,nr2
        if (jsflag(ns) .le. 0) then
          nrn1 = ndrsr(1,ns)
          cx = (-1.0/cdrs(nrn1))
          sisppu = cx*sidrsp(ns)
          xx = xbar(ns)
          xl = xbarlg(ns)
          sisp = sisppu - stx*xl
          acflg(ns) = stxm1*xl
          act(ns) = xx**stx
          actlg(ns) = stx*xl
          sisp = sisp/cx
          sidrsp(ns) = sisp
          affs(ns) = afcnst*sisp
          siph = siph + xx*sisp
        endif
      enddo
      sidrph(np) = siph
      affp(np) = afcnst*siph
c
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
